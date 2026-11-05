# Codex実装・検証指示：声門流量の再アンカー処理を修正し、全有効区間Jitterと目視検証を確立する

> 作成：2026-11-05 / 対象ローカルリポジトリ：`simulation-main`  
> **今回の最優先目的は「周期イベント検出とJitter算出の信頼性の確立」**。剛性比0.30で観察された後半の振動安定化の**物理的原因分析は今回行わない**。  
> **既存の`tools/analyze_airflow_metrics.py`を必要最小限変更する。Pythonソースを別ファイルへ分割しない。C++ソルバ・元の`airflow_vt.dat`・既存解析出力・ユーザーの未コミット差分は非改変。**  
> **再アンカー改善・自動検証・目視用データ生成・手動注釈の評価機能まで一括実装**する。ユーザーが注釈をまだ提供していない時点で「人手検証済み」とは宣言せず、手動入力が必要な部分だけ明示する。

## 0. 今回の入力根拠と検証用ベースライン

ユーザー提供の`investigation_report.md`に基づく。900 Pa、解析窓`[0.50,1.00) s`、`dt=1e−5 s`、出力`stride=5`（20 kHz）。条件は**剛性比`r=0.25`と`r=0.30`**。対応する`airflow_vt.dat`はローカルの既存runを探索し、実パスと`SHA256`を特定する。**不明なファイルパスを推測して誤入力を使わない**。

既知の厳格版（改修前）の値：

| 指標（主法A） | r=0.25 | r=0.30 |
|---|---:|---:|
| 既存の最長区間Jitter local [%] | 9.828 | 1.741 |
| **分断をまたがない全有効区間の集計Jitter [%]** | **10.873** | **6.486** |
| 全有効区間RAP [%] | 6.679 | 3.885 |
| 全有効区間PPQ5 [%] | 5.851 | 3.458 |
| 全追跡周期数 | 36 | 54 |
| 最長区間周期数 | 8 | 26 |
| covered_time_fraction | 54.52% | 86.13% |

- r=0.30の最長区間は`0.766345–0.989736 s`の26周期。前半の有効区間には10–13%のJitterも存在。**最長区間だけを全窓の代表として扱うのは不適切**。
- AとBの同一窓Jitterはr0.30のA最長窓で**1.741% / 2.727%**、r0.25のA最長窓で**9.828% / 4.369%**。異なる位相のイベントからJitterを計算するため、無理に同値にしない。
- 中断原因：r0.25主法Aは相関不足15判定・「候補なし」15判定（**候補なしのうち13件は生ピークが実在するが再アンカー除外済み**）。r0.30主法Aは相関不足4・候補なし4（**3件は再アンカー除外済み、残り1件は端の照合窓不足**）。
- **候補間スコア差を0～0.15で変更しても両条件のA検出結果は不変**。主ボトルネックは再アンカー除外と形状相関の失敗。
- 相関閾値を`0.65→0.55`に下げるとr0.30で被覆率94.11%、最長53周期のJitter6.76%。ただしr0.25では別の局所極大へ位相移動する例があるため**一律に閾値を緩めない**。

**改修前の結果を`baseline`として固定して保存**。まずローカルの`investigation_report.md`、`summary.json`、`break_diagnostics.csv`、`break_candidates.csv`、`reanchor_exclusions.csv`、`new_connections_vs_strict.csv`、`same_window_metrics.csv`、`preservation_check.json`の有無を調べ、あれば現物を参照する。GitHub mainに診断スクリプトの最新差分がない可能性があるので、**ローカルの`git status`と実際のコードを正とし、mainで機械的に上書きしない**。

## 1. 目的と成功判定

**目的**：
1. すでに見つかっている生ピークを、再アンカー除外の履歴だけが理由で使えなくなる問題を修正する。
2. **波形の同じ位相に対応するイベントの追跡を維持**し、誤ピーク採用・1周期飛ばしを増やさず、有効周期の時間被覆を改善する。
3. r0.25/r0.30の**全有効区間Jitter local、RAP、PPQ5**を、区間境界をまたがない共通ルールで求める。
4. **人が実際の流量波形を見てイベント時刻を指定し、自動検出を検証できる一式**を生成する。手動注釈が未提出でも、注釈用の波形図・CSVテンプレート・照合CLIまで完成させる。
5. 閾値感度・サブサンプルON/OFF・A/B同一時間窓・既存のCPPを保持し、**数値だけでなくイベント対応の信頼性**を報告する。

**「全区間からデータ」の意味**：0.50–1.00秒のすべてのサンプルを入力として探索し、追跡できたすべての**独立した有効連続区間**からJitter項を集計すること。**被覆率100%や全周期の必須接続を要求しない**。曖昧な区間の無理な接続、実在しないイベントの挿入、既知周期への時刻スナップは明確に禁止。

## 2. Phase A：再アンカー処理の原因箇所の特定と修正

### A1. 現行データ構造を調査

現行`tools/analyze_airflow_metrics.py`の次の各関数/状態を追う：

- `find_peaks`で検出した全`raw candidates`。
- 追跡で`accepted`となったイベントID、`rejected`候補、`reanchor excluded`候補、追跡区間`segment`、anchor選択規則。
- 一方向の失敗後、反対方向/別anchorから同じピークを再検討できない条件。
- すでに採用済みイベントと単に「失敗地点の近傍にある」だけの候補が混同されていないか。
- 同一ピークの二重採用、異なるトラックの時間的重複、前後両方向の競合、二峰の位相切替の防止条件。

修正前のr0.25/r0.30について**exactな採用ID・区間・結果・CPP一致を回帰テストとして保存**。既存の診断結果は上書きしない。

### A2. **再アンカーの除外規則を変更**

以下を設計要件として実装する。関数名・内部表現は既存コードを優先。

- **`raw candidate`は解析中削除しない**。`used_by_accepted_track`、`rejected_for_context`、`conflicted_with_accepted_event`、`unattempted`など**理由別の状態**を持つ。
- 直前イベント`i`から候補`j`への辺（追跡接続）で相関`ρ`やスコア差が不合格だった場合、**その失敗は原則`(i,j,方向,参照テンプレート)`という文脈依存の試行履歴**とする。**候補`j`そのものを以後の全anchor候補から永久除外しない**。
- 従来の「失敗地点から±1.45`T_ref`内の全候補を将来の再アンカーから除外」は廃止/限定する。**広い時間近傍にいるだけで未採用ピークを排除してはならない**。
- 明確に排他するのは、既存トラックで採用済みのイベントそのもの、または同じ物理周期の別の局所峰である根拠がある候補。単純な時刻距離だけで別周期のピークまで除外しない。
- 失敗後は近傍の**未採用候補を別のanchorとして再評価できる**。新しいanchorから前向き/後ろ向きに評価し直し、テンプレート/形状相関が別文脈で妥当なら採用可能にする。
- 同一anchor状態・同一試行を繰り返して無限ループしないよう、**contextをキーにしたvisited/attempt cap**を設ける。任意に上限到達で探索を中断する場合は`SEARCH_LIMIT`をQCに記録。
- 異なる追跡区間を統合する際は**時間順序・イベントID一意性・位相/形状の連続性・ACF基本周期候補との整合**を評価する。**中間の1周期を欠落させて2Tを1周期にする接続は禁止**。統合が曖昧なら元の`segments`を分離のまま残す。
- A（極大）とB（極小）は同一周期の別位相であり、片方の採用イベントをもう片方の`ground truth`として強制しない。

### A3. 周期の位相飛び・倍周期・形状変化への対策

- `prominence=10%`、最低形状相関`ρ=0.65`、スコア差`0.10`、探索窓`[0.65,1.45]T_ref`、ACF基本周期、Hann/CPP設定は**まず既定値を変更しない**。
- r0.25で報告された`0.7975s`付近の局所峰への位相移動（相関`ρ=0.50`時）、r0.30のprominence3%での`0.5567s`の後方局所峰、r0.25 Bの`0.5857→0.600428s`（約1.95倍`T_ACF`）は**必須の誤接続回帰例**として図とCSVで検証。
- イベントの採否には波形形状相関だけでなく、同位相であることの根拠（前後の谷/山の位置関係、周期内相対位相、ピーク形状の特徴）も**診断値**として持たせる。どの形状ゲートが未実装かを明示。ユーザーの許可なく恣意的な新スコア重みを本番に導入しない。
- 形状が激変して本当に対応が曖昧なら区間を分断する。**不規則振動のJitterを人工的に小さくする規則禁止**。
- `MISSING_EVENTS_POSSIBLE`を、少なくとも`NO_RAW_PEAK` / `BELOW_PROMINENCE` / `REANCHOR_EXCLUSION_CONFLICT` / `LOW_SHAPE_CORR` / `INSUFFICIENT_TEMPLATE_EDGE` / `TRACK_CONFLICT`のように分解し、単なる候補管理の失敗を「物理的な流量パルス消失」と表現しない。元のフラグと対応可能にする。

## 3. Phase B：**全有効連続区間からのJitter/RAP/PPQ5集計**

**既存の最長区間`longest_segment`は後方互換のため残し、**研究用の主候補を`pooled_within_segments`として新設/明確化する。既存フィールドを説明なしに意味変更しない。

有効な連続区間`s`で周期列`T_{s,1},…,T_{s,n_s}`を得る。**区間境界を跨ぐ周期差は計算しない**。

**Jitter local [%]**：
\[
J_{\rm pooled}[\%]=100
\frac{\displaystyle\frac{\sum_{s:n_s\ge2}\sum_{i=1}^{n_s-1}|T_{s,i+1}-T_{s,i}|}{\sum_{s:n_s\ge2}(n_s-1)}}
{\displaystyle\frac{\sum_{s:n_s\ge2}\sum_{i=1}^{n_s}T_{s,i}}{\sum_{s:n_s\ge2}n_s}}.
\]

**RAP**：各区間に3周期以上ある場合のみ`i=2…n_s−1`の局所3周期平均との差の絶対値を計算し、すべての有効項の平均を、**RAPに使用した3周期以上の区間に含まれる周期の平均`T`**で正規化する。

**PPQ5**：各区間に5周期以上ある場合のみ`i=3…n_s−2`の5周期平均との差の絶対値を計算し、すべての有効項の平均を、**PPQ5に使用した5周期以上の区間に含まれる周期の平均`T`**で正規化する。

**差分上の重要点**：従来の調査レポートでは**1周期しかない区間をJitterの分母平均周期に含める**集計だった。今回の研究用`pooled_within_segments_v2`では、**各指標で有効な項を生成できる区間だけ**を分母平均へ含める。従来の`10.873/6.486%`等の**`pooled_v1`は過去数値の回帰確認用に保持**し、v2と混同しない。各`method/version`、周期数、差分項数、mean period、coverageを出力して差の原因を明示。
- ユーザーに提示する結果表は`longest_legacy`、`pooled_v1_baseline`、**`pooled_v2_current`**を列で区別。
- `pooled`は**追跡できた部分を集計した値**であり、未追跡区間を含む「0.50–1.00秒全体のJitter真値」と呼ばない。
- 複数segmentで平均周期が大きく違う場合、全体集計は過渡/非定常の影響を受ける。**period mean・各segmentのJitter・seg duration・時間的位置**を必ず出す。今回の目的は**アルゴリズムの信頼性確認**であり、r0.30の後半安定化の物理的機構には踏み込まない。
- 評価時間窓外の不完全な周期を集計しない。各区間`[t_start,t_end]`を`periods.csv`のイベントIDで検証可能にする。

## 4. Phase C：**目視検証を一度の実装で可能にする**

目視で周期位置が分かることを利用し、**ユーザーが画面上で確認・修正した周期イベントを真値の暫定基準`manual reference`として、追跡器のaccuracyを直接比較できる機能**を同じ`tools/analyze_airflow_metrics.py`に実装する。

### C1. レビュー資料を自動作成（注釈がなくても実行）

- r0.25/r0.30の評価窓`[0.50,1.00) s`を、**50～70 ms程度の連続パネルに区切って全体をカバー**。重複は例として1周期分（イベント見逃し防止）。**端で評価窓外に出るものをJitterに含めない**。
- 各パネルに**(a)生流量だけのブラインド注釈用**と**(b)候補、採用A/B、中断点、segment境界、予測`T0`、ρ/scoreを重ねた診断用**の2図を生成。図の名称/時間範囲を一致させ、凡例を分かりやすくする。
- 必須の拡大診断地点：
  - r0.30 `0.5497–0.5644s`（`ρ≈0.199`）、`0.6815–0.6965s`（`ρ≈0.5625`）、`0.7491–0.7664s`（`ρ≈0.631`）。
  - r0.25 `0.6141s`付近（`ρ≈−0.389`）、`0.7975–0.8000s`付近（位相飛びの危険）。
  - 比較的安定して追跡できた**対照区間**も各条件1～2箇所。
- どの図も**波形に勝手に平滑化をかけない**。元の20kHz流量`Ug`、時間軸[s]または[ms]、極大/極小と採用イベントを識別する記号を使う。
- `manual_review_index.csv`：`case_id,panel_id,t_start_s,t_end_s,figure_blind,figure_diagnostic,break_ids,notes`。
- 目視時に必要なら`--manual-review-start/end`で一部区間だけ生成可能。**全窓版を作るのが既定**。

### C2. **手動イベントの注釈形式と評価CLI**

最低限、追加依存なしの`CSV入力`を必ずサポートする。GUIでのマウスクリック入力が容易なら**同じ1ファイルに`--annotate-interactive`（Matplotlib）を追加してよいが、GUIがない環境でもCSVで完結すること**。

**`manual_events.csv`の形式**：

~~~csv
case_id,phase,event_time_s,confidence,reviewer,notes
r0p30,A,0.512344,high,manual,large-flow-peak
r0p30,A,0.519820,high,manual,next-corresponding-peak
r0p30,B,0.516210,high,manual,flow-trough
~~~

- 上の数値は**形式の例示であって正解時刻ではない**。ユーザー用テンプレートには架空のイベント時刻をあらかじめ真値として埋めない。
- `phase=A`は**同種の流量極大**、`phase=B`は**同種の流量極小**。同じ条件内で同じ位相基準を保つこと。中断点の前後をまたぐ手動イベントも、目視で追えるなら続けて記録する。
- `confidence=high/uncertain`を選べる。曖昧な箇所は`uncertain`として残し、**highのみの精度と全注釈の精度を分けて報告**。
- 最初に**r0.25/r0.30それぞれ少なくとも30周期分の高信頼イベント**を確認できるようレビュー図・テンプレートを作る。望ましくは0.50–1.00 sの全イベント。**中断地点と安定区間を両方含める**。この人手作業をCodexが勝手に実施済みと見なさない。
- `--manual-events <path>`でCSVを読み、単位・増加時刻・重複・解析窓・流量イベント型の整合を検証。
- 注釈対応は**同一phaseで時刻差最小の一対一対応**（近傍イベントが複数ある場合は全体で一意になる割当）。**許容誤差±100μs（20kHzで2サンプル）**を基準に、50/100/150μs感度も併記。別phase間の対応はprecision/recallに使わない。
- 人手高信頼イベントを`TP/FP/FN`に分類して`precision/recall/F1`、時刻誤差中央値/95パーセンタイル、周期の欠落・2重検出・位相飛び、**人手周期列由来のJitterと自動検出Jitterの差**を報告。比較は**同一時刻範囲・同じphase・完全に含まれる周期だけ**。区間境界を跨ぐ差を作らない。
- ユーザーが手動注釈をまだ入力していない段階でも**目視レビュー図・CSVテンプレートの作成と`--self-test`までを完了**。手動検証の統計欄は`PENDING_MANUAL_ANNOTATION`と表示。人間のレビューが未実施でも達成済みと報告しない。

## 5. Phase D：回帰テスト、ストレステスト、比較

### D1. 必須の合成テスト（既存`--self-test`に追加）

単峰の既知イベント列、複数極大/周期、振幅変調、時間とともに波形形状が変わる、基準相関0.65前後の接続、見かけの2周期飛び、局所峰位相ジャンプ、欠落候補、窓端で照合窓不足、前向き・後ろ向きの競合、テンプレート変更による再評価、ループ/再試行上限をテスト。

**期待**：
- 中断候補が存在するときは履歴の除外だけで`NO_CANDIDATE`にしない。
- 正解の基本周期を不適切に2分割または2倍接続しない。
- 再アンカー改善後も`events`は時間的に一意、トラック境界は交差しない。
- Jitter/RAP/PPQ5の数式が既知周期列に一致、**区間を跨ぐ差0件**。
- `pooled_v1`で旧診断値を再現、`pooled_v2`では短区間の分母ルールを明示。
- 本体のCPP実装と設定（10kHz、60ms、5ms hop、2048FFT、Theil–Sen）は**不変**。
- `manual_events`の入力検証、許容時刻差と一対一対応、precision/recall計算の合成テスト。

### D2. 実データr0.25/r0.30に必ず実行する比較

| 比較 | 出すもの |
|---|---|
| baseline strict vs new strict | A/B採用イベントID・候補数・追跡被覆・break理由・段階別スコア |
| Jitter longest vs pooled_v1 vs pooled_v2 | 4指標と平均周期、周期数、差分項数、time coverage |
| 前後半`[0.5,0.75)`/`[0.75,1.0)` | 全有効区間Jitterと追跡coverage。**安定化の物理的原因は考察しない** |
| CPP before vs after | **完全同一**（浮動小数点丸め差以外）。`CPP`を変更するなら止めて報告 |
| Jitter A/B same absolute windows | 同じ完全部分に含まれるイベント周期のみ比較。異phaseの差を誤検出とは即断しない |
| threshold sensitivity | `ρ`、margin、prominence、補間ON/OFFによる**一致率とJitter変動** |
| 増えた接続のチェック | 新規イベントの前後波形・`T/T_ACF`・候補順位・A/B対応・phase維持を図にまとめる |

**成功の条件は「被覆率が上がった」だけではない**。**全新規接続で位相飛びや2周期飛ばしがないかを検証し、曖昧なものはacceptedにしない**。基準を保ったまま改善できない条件は`QC_WARN/PENDING_MANUAL_VALIDATION`で保存し、改善したと偽らない。例えばr0.25の被覆が目標90%に達しないこと自体は**失敗扱いにしない**。検出信頼度を維持する方が優先。

## 6. 出力先とレポート、CLI

新しい出力先は既存の`analysis_runs/`内の**専用新規ディレクトリ**（例`analysis_runs/airflow_reanchor_validation_20261105/`）とし、既存の`r0p25`/`r0p30`出力へ上書きしない。必要であれば各条件を別サブディレクトリへ。既存の`analysis_runs`が`.gitignore`で無視されているか確認。

推奨する追加CLIオプション（既存オプションと矛盾する場合は既存形に合わせて実装）：

~~~bash
# 既存の入力パスを実際に確認してから置換。新規出力ディレクトリに診断込みで保存
python3 tools/analyze_airflow_metrics.py \
  --input "/absolute/path/to/r0p30/airflow_vt.dat" \
  --dt 1e-5 --start 0.50 --end 1.00 \
  --reanchor-policy recover-rejected \
  --diagnose-tracking \
  --export-manual-review \
  --out-dir "analysis_runs/airflow_reanchor_validation_20261105/r0p30"

# 同じ設定でr0p25を解析し、値/被覆/誤採用を比較
python3 tools/analyze_airflow_metrics.py \
  --input "/absolute/path/to/r0p25/airflow_vt.dat" \
  --dt 1e-5 --start 0.50 --end 1.00 \
  --reanchor-policy recover-rejected \
  --diagnose-tracking \
  --export-manual-review \
  --out-dir "analysis_runs/airflow_reanchor_validation_20261105/r0p25"

# 人手イベント記入後に実行。人手注釈が無いときはテンプレートを作って終了
python3 tools/analyze_airflow_metrics.py \
  --input "/absolute/path/to/r0p30/airflow_vt.dat" \
  --dt 1e-5 --start 0.50 --end 1.00 \
  --reanchor-policy recover-rejected \
  --manual-events "/absolute/path/to/manual_events_r0p30.csv" \
  --out-dir "analysis_runs/airflow_reanchor_validation_20261105/manual_check_r0p30"
~~~

`--reanchor-policy legacy`で元の挙動を再現できるようにする。**今回の変更で旧`summary.json`のキーや意味を黙って変更しない**。追跡器version、policy、解析窓、ソルバdt、stride、閾値、データSHAを`summary.json`に保存する。

**最低出力**：
- `reanchor_comparison.csv`：条件・旧/新イベント・新規接続・coverageと破綻理由
- `tracking_breaks_detailed.csv` / `reanchor_candidate_states.csv`：候補の実在/不採用/失敗文脈・anchor再利用
- `periods.csv` / `segments.csv`：イベントID、各周期、各区間、採否・phase・平均周期
- `jitter_comparison.csv`：longest、pooled_v1、pooled_v2のJitter/RAP/PPQ5、採用項数・covered_fraction
- `manual_review_index.csv`、`manual_events_template.csv`、blind/diagnosticの**全窓パネル画像**、中断点拡大図
- 注釈後のみ`manual_validation.csv` / `manual_validation_report.md`（未注釈なら明示的`PENDING`）
- `report.md`：0.25/0.30で何が改善し、何が依然曖昧で、何が人手待ちか。**最長区間と全有効区間集計を区別したまま**報告。
- `summary.json`：保存されたCPP、設定hash、解析器version、baseline比較、QC

## 7. 作業環境・コミット・報告ルール

- まず`pwd`、`git status --short`、`git branch --show-current`、`AGENTS.md`、ローカルで本当に使われている`tools/analyze_airflow_metrics.py`の最新版と既存テスト・診断を確認。**事前のユーザー未コミット差分や手元の解析結果を勝手に削除しない**。
- **新規Pythonファイルを増やさず**、既存`tools/analyze_airflow_metrics.py`の変更が中心。必要なテストは既存の`--self-test`へ統合。`numpy/scipy/matplotlib`等の既存依存を使用。
- シミュレーションは**再実行しない**、`input/`・`output`・既存`analysis_runs`は非改変。`CPP`と他の物理指標は同一。新規`results`と図は指定`out-dir`に隔離。
- 正常時・異常時とも一時ファイルを残さず、ユーザーのGit作業状態を保存する。`git reset/clean/checkout`や無断`push`は禁止。コミットは通常のユーザー/Codex運用に従う。
- 最終報告には**(a) baseline vs new strictの数値表、(b) 新規接続の誤位相検査、(c) pooled_v2のJitter/RAP/PPQ5とcoverage、(d) 目視確認用パネルの場所、(e) manual precision/recallが未完ならPENDINGの明示、(f) 残る失敗例とその理由**を必ず含める。
- **本実装の完成と「人手目視による正解照合が完了した」は別概念**。今ある情報だけでできるのは`auto diagnosis`、`review sheet生成`、`annotation comparisonコードのテスト`まで。人手で正解イベントを確認・注釈した後、実データで`precision/recall`等を確定する。

## 8. Codexへの最終依頼

> 今回は剛性比0.30の後半安定化の物理的原因は調べず、**流量の周期イベント追跡器の品質を完成させることに専念**してください。既存`tools/analyze_airflow_metrics.py`の再アンカー処理が、実在ピークを失敗履歴の近傍除外によって候補から失わせている問題を修正してください。ただし相関の既定閾値0.65を一律に緩めたり、周期をT0へ強制補正したりしないでください。剛性比0.25・0.30の既存0.5–1.0秒流量に対し、改修前後の全イベント/区間/被覆率と、最長区間および**分断を跨がない全有効区間Jitter/RAP/PPQ5**を比較してください。あわせて、**全窓のブラインド/診断用拡大図、手動周期イベント入力CSVテンプレート、その後の一対一時刻照合とprecision/recall/Jitter誤差評価機能**を同じPythonファイルにまとめて実装してください。人手イベント未入力なら`PENDING_MANUAL_ANNOTATION`と正確に報告してください。新規採用の周期に誤った位相移動/周期飛ばしがないことを個別に検証し、既存のCPP・ソルバ・元流量・既存結果・ユーザー差分は保護してください。
