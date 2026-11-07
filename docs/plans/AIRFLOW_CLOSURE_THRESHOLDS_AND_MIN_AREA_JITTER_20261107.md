# Codex実装・検証指示：流量下降時の相対閾値イベント（0.5/2/5%）と最小声門断面積波形によるJitter

作成：2026-11-07。対象：`simulation-main`、既存`tools/analyze_airflow_metrics.py`。

## 目的と基本方針

剛性比`0.25`、`0.30`、声門下圧900 Paの既存の同一シミュレーション結果を用い、**2系統の周期イベント定義を追加してJitter local・RAP・PPQ5を評価**する。

1. **声門流量`Ug(t)`：流量が閉鎖に向けて低下するとき、直前の局所ピーク流量の`0.5%`、`2%`、`5%`を下向きに通過する時刻を、周期イベントとする。** これは`Ug=0`の平坦部の任意点ではなく、**下降側の閾値交差**の時刻である。3閾値は必ず同じデータで比較する。
2. **最小声門断面積`Amin(t)`：既存の`global-cycle`法と同様の考え方で、**極大値および極小値**を別々に追跡する。データがないときは発明せず、入力を検出・報告し、CLIで指定可能にする。
3. **解析は`0.50–1.00 s`の全区間を対象**とする。旧式の「難しい箇所でtrackingを打ち切り、最長の連続区間だけでJitterを出す」という方法は禁止。局所ACFを使って全区間で探索を継続し、**全区間の時刻被覆を出力**する。
4. ただし本当に閾値通過/物理的周期がない区間に**イベントを捏造したり、時刻を自己相関の予測周期`T0`にスナップしてはいけない**。未解決箇所を可視化・数値化し、報告した上で、可能な限り全窓の観測イベントを回復する。
5. **今回の第一目的は評価手法の比較**。0.30で時間とともに振動が安定化する物理的理由の説明や、全36条件への適用は行わない。CPPの算出法、C++ソルバ、元データは変更しない。

既存の`report(10).md`と前回指示書`docs/plans/AIRFLOW_FULL_WINDOW_EVENT_INFERENCE_THRESHOLD_REVIEW_CLEANUP_20261106.md`の新実装を前提にする。**現行ローカルコードがGitHub mainより進んでいる可能性があるため、`git status`・`AGENTS.md`・最新版`tools/analyze_airflow_metrics.py`を必ず先に確認**。未コミット変更は保護。

## 1. データとbaseline

### 1.1 既知の対象と時間窓

- 両条件：`Ps=900 Pa`、`solver dt=1e−5 s`、流量`step stride=5`（20 kHz）、解析窓`[0.50,1.00) s`。
- 流量r0p25：`analysis_runs/manual_run_archive/r0p25_P900_20261104T064537_9fba0102/snapshot/dat/airflow_vt.dat`。既知SHA256`936e26e7d8284e11b76093f19480a46cdb123e33b6457677bd44c8139383cdca`。
- 流量r0p30：`analysis_runs/manual_run_archive/r0p30_P900_20261104T071153_3b5765f4/snapshot/dat/airflow_vt.dat`。既知SHA256`061df933beb293e803d2f50b4c924516803b17dc503f618669fe56e742b7b840`。
- **最小断面積の入力`area.dat`（または現行出力名）は実ファイル・列構造を現物で確認**：`snapshot/dat/`以下、既存の`calcArea` / `output`出力、ヘッダ/単位/step番号/複数面積列/欠測値を調査。**`Amin`を「最小断面積」以外の総面積/平均面積で代用しない**。column index/nameが不明なら明示的CLIで選べるようにし、検証なしに推測で数値を作らない。
- 時刻変換は各データの**`t=step×dt`**を使用し、ファイルごとのstrideを実測。流量と面積が同じサンプリング間隔とは限らない。片方をもう片方の時間点に無理に再サンプルせず、**イベントは元の時刻系で算出**してから比較する。タイムオフセットがあれば検証して明記。
- `r0p25`、`r0p30`の`Ug`と`Amin`が**同一run**由来かmanifest/hash/step分布で確認。もし`Amin`が見つからなければ、流量3閾値を先に完成させ、面積系は`MISSING_INPUT`とし、後で`--area-input`を渡して再実行可能にする。**ソルバ再実行はしない**。

### 1.2 比較baseline（改変禁止）

既存`global-cycle balanced-v1`の流量極大A/極小B、`recover-rejected`、`pooled_v2`、CPPを再現・保全。報告値：

| 条件・入力イベント | global-cycle Jitter local [%] | RAP [%] | PPQ5 [%] | Jitter対象時間 |
|---|---:|---:|---:|---:|
| r0p25 Ug極大A | 14.728 | 8.127 | 8.746 | 96.61% |
| r0p25 Ug極小B | 10.783 | 6.259 | 6.598 | 98.48% |
| r0p30 Ug極大A | 7.028 | 4.403 | 4.625 | 98.78% |
| r0p30 Ug極小B | 7.014 | 4.050 | 4.851 | 95.01% |

旧CPPは`r0p25=15.6578928242dB`、`r0p30=16.9117537400dB`。**今回CPPの計算関数と設定を変えない**。

## 2. 新方式A：下降側の局所ピーク流量比に基づく周期イベント

### 2.1 正確な定義（重要）

各流量パルス`i`における、**主たる上昇・下降を伴う局所最大流量`Upeak,i`**を特定する。

\[
U_{\theta,i}=(\theta/100)\,U_{\mathrm{peak},i},\quad
\theta\in\{0.5,2,5\}\%.
\]

そして`Upeak,i`の**後の下降枝**で、`Ug(t)`が`Uθ,i`を**上から下へ初めて安定的に横切る**時刻`t_close,i(θ)`を周期イベントとする：

\[
Ug(t_{\mathrm{close},i}(\theta))=U_{\theta,i},\quad
dUg/dt<0\quad\text{(descending crossing)}.
\]

- **流量が0になっている区間の中心、最小流量のサンプル時刻、流量の極大時刻はこの指標のイベントではない**。
- 「局所ピーク」は**このパルスに対応する主ピーク`Upeak,i`**。長い区間の最大、全解析区間の最大、局所パルス内の小さな二次峰を無条件に採用するわけではない。**0.25の多峰流量**ではどのピークを選んだかと代替候補を出力。
- 流量のベースラインはユーザー指定の**0 (m³/s)**。`(Upeak−Umin)`の0.5/2/5%ではない。局所流量が常に正で、`Uθ,i`まで下がらない場合、**その閾値に対する実際の下降交差は存在しない**ため、フラグ`NO_THRESHOLD_CROSSING`を付ける。勝手に最小点を`Ug=0`の交差時刻扱いしない。
- `Ug`の数値クリッピングで0に張り付く場合も、**下降枝で最初に閾値を通過する瞬間**を拾う。これを**声帯の物理接触時刻と断定しない**（音源流量のclosure-like marker）。
- 生波形のサンプリング点間の上下交差から**線形補間**で`t_cross`を求める。極大位置の放物線補間とは別。各閾値の`t_cross`と隣接サンプルID/流量をCSVに残す。

### 2.2 波ごとの候補作成・二次峰対策

- 局所ACF`T0(t)`（既存100ms窓/20ms hop）を使用して、全`0.50–1.00s`を予測周期に沿って走査する。ただし`T0`は**探索範囲のガイドのみ**、算出時刻に使わない。
- 先行する開放相の主流量ピークから閉鎖に向けた**主下降枝**を探索。`0.25`の多峰では極大候補を複数保持し、前後の谷、ピークprominence、局所上下振幅、主パルスの面積/幅、同じ位相の前後パルスとの整合で主ピークを選択。**既存global-cycleのAピークは候補には利用するが、全てを正解seedとして固定しない**。
- **反発・小さな再上昇**で同じ閾値を何度も通過する場合は、波を1周期とした**主閉鎖下降枝の安定交差**を選ぶ。初期案：`crossing`後に少なくとも`0.03T0`（ただし生標本2点以上）の間、`Ug≤Uθ`が維持されること、または直後に真の零流量`plateau`へ至ることを確認。局所に小さな`bounce`があり無理な条件になる場合は**ALTERNATIVE_CROSSINGS**として複数保持、候補を捨てて強制選択しない。
- 閾値ごとに周期数が変わって見えないよう、**共通の候補パルスID/パルス区間をベースに3閾値を評価**する。閾値が低いほど交差なしになることはあり得るので、それぞれの`event_found`・理由を記録し、足りないイベントを補間生成しない。
- `0.5%`は急峻な閉鎖部・数値ノイズ/0クリップに近く不安定な可能性、`5%`はより早い下降時刻を取る可能性がある。**どれが正しいと最初から決めず、同一波形で比較**。
- 各パルスで`Upeak_i, Uθ,i, Umin_i, peak/prominence, closing-slope, multiple_crossings, sample indices, cross_time, no_crossing`など全情報を`flow_threshold_events.csv`に出力。

### 2.3 3閾値で共通の周期列・結果

各`θ`について、採用下降交差時刻を昇順に並べ、

\[
T_i^{\theta} = t^{\theta}_{i+1}-t^{\theta}_{i}
\]

から従来と同じJitter local、RAP、PPQ5の式を使う。**連続している同種パルス間のみ**を隣接周期とし、`NO_THRESHOLD_CROSSING`を挟んで2Tを1周期として数えることは禁止。失敗区間はそこで分割・再開し、`pooled_v2`のように**全有効segmentにある差項**を集約（最長だけに限定しない）。

**全窓の欠測を小さくする追加規則**：初めに極大ピークベースのglobal-cycleで周期候補を定義しても、既存方式が途切れたところでは`Ug(t)`の**閾値下降交差の候補を直接探索**し、別アンカー/後ろ向きから回復。`NO_THRESHOLD_CROSSING`が残れば主値と別に`coverage`、`missing pulses`、未交差理由を表示し、その区間が不規則であっても結果から消さない。**「評価全域」=全時刻をチェックしたことであって、存在しない交差を追加したことではない**。

可能なら**3閾値それぞれで閉鎖イベント数・同一pulseIDへの対応・イベント時刻差の分布**を評価し、閾値変更が実際のイベント時刻だけでなくpulse選択まで変えていないか分離する。

## 3. 新方式B：最小声門断面積`Amin(t)`の極大・極小追跡

### 3.1 入力解析

- 実際の`area.dat`のフォーマットを先に調査。ステップ整数+面積値の複数列か、ヘッダ・単位`mm²`/`m²`か、最小面積の列名を確認する。**既存`area.dat`の`minA`以外の値を誤選択しない**。
- 単位を`mm²`等へ統一し、負値、0クリップ、欠測、step重複、時刻の非単調、可変stride、異常な開口値をQC記録。**最小断面積が接触penalty等で0に丸められている場合は測定下限を識別**。
- `--area-input`、必要に応じ`--area-column`、`--area-time-column`、`--area-unit`をサポート。自動検出できる場合も選んだ列を報告。

### 3.2 イベントの抽出

- `Amin(t)`の局所**極大（open-area peak）**と**極小（closed-area trough）**を独立に`find_peaks`と`find_peaks(-Amin)`で抽出する。
- **現在のglobal-cycleの局所ACF＋局所prominence/上下振幅＋前後の波形形状・時間整合＋全窓候補列最適化**を、信号`Amin`に適用する（可能なら信号別設定の汎用化）。**面積波形専用の都合の良い閾値調整をせず、既存balanced-v1を基準に感度比較**。
- `Amin=0`の平坦な谷は**min-plateau**。`find_peaks(-Amin)`が単一の点を返す場合でも、**plateau開始・終了時刻・幅・代表時刻（初期は中央）**を出力。**plateau最小値のJitterは代表時刻の選択に依存する**ため、開始/中央/終了の感度を必ず併記。極大はピーク補間ON/OFFを比較。
- **Aminの極大Jitter**と**Aminの極小Jitter**を別々に出し、流量閾値3方式とも区別したラベルで保存。
- 流量下降側`tclose(θ)`との絶対遅延`Δt=t_Aminmin−tclose`の分布、同一pulseID/phaseの対応・欠落率を**診断として**示す。**流量閾値到達時刻と最小面積が一致することを強制しない**。時間遅延の正負を物理的な因果と断定しない。

## 4. 必須：解析範囲を途切れさせない設計

- **どの方式でも全窓`[0.50,1.00) s`に対してイベント候補の抽出を必ず行う**。曖昧な位置があってもその後の探索を止めず、再アンカー・局所ACFの更新・グローバル列選択で全時刻へ復帰。
- **解析対象時間の100%を**、少なくとも`VALID_EVENT_PERIOD` / `PROVISIONAL_EVENT_PERIOD` / `UNRESOLVED` / `WINDOW_EDGE`として分類する。`classification_coverage=100%`は`valid_jitter_coverage=100%`とは別。
- Jitterは**全ての有効連続区間を使い、区間境界をまたぐ周期差を絶対に作らない**。最長区間のみの値を「主値」と表示することは禁止。指標別`used_time_fraction`、`N_events`、`N_periods`、`N_Jitter_terms`、`mean_period`、`N_unresolved`、`N_provisional`を出す。
- **低信頼イベントを採用して全窓を連結する場合**はイベント候補と不確かさを`PROVISIONAL`で記録し、例えば`confirmed-only`と`including-provisional`の**両方**のJitterを出す。**架空のイベント挿入・予測周期の強制補正・異なる位相間の乗換えは禁止**。
- 閾値未達でも面積極値や流量ピークが明瞭な場合、その情報から区間をまたいだイベントの存在可能性は検討できる。ただし**「他の指標から補完した流量0.5%の閾値通過」として同じJitterへ混入させてはいけない**。
- 本当に閾値交差のない波や無周期状態では、**閾値式Jitterはその区間について定義不能**。未解決率とパルス振幅、局所ACF信頼度、CPP等を合わせて報告。**難しい区間を「綺麗な区間だけ」へ切り捨てていないかを検証**。

## 5. 解析条件、閾値感度、結果比較

### 5.1 固定する設定

- 解析窓`0.50–1.00s`、局所ACF主`100ms/hop20ms`、周期周波数候補`60–400Hz`、サンプル時刻`step×dt`。
- 流量イベント定義`θ=0.005,0.02,0.05`（ピーク比）、共通`Upeak_i`選択規則、下降交差。
- 面積イベント定義`maxima/minima`、主balanced-v1相当の候補選択、plateau-minは代表点の感度比較。
- 同一設定を`r0p25`、`r0p30`で使う。値が良い条件に合わせて片側だけ変えない。

### 5.2 比較対象と数値表（**一つの統合テーブル**）

各rにつき最低**5種類**のJitter：

1. `Ug_downcross_0p5pct`（閉鎖下降側0.5%）
2. `Ug_downcross_2pct`
3. `Ug_downcross_5pct`
4. `Amin_maxima`
5. `Amin_minima`

比較の基準として旧`Ug_global_maxima_A` / `Ug_global_minima_B`を併記する。各々に`Jitter local[%]、RAP[%]、PPQ5[%]、mean T[ms]、N periods、event count、event confidence、time coverage、ambiguous/unresolved fraction、使用項数、QC`を出す。CPPは既存の**流量CPPを参照比較のみ**。面積CPPを勝手に流量CPPと同一扱いしない。

- **流量3閾値同士のイベント時刻差とJitter差**を可視化：特に0.5%が0/ノイズ床の影響を受けやすいか、5%が早期交差に偏らないか。
- **Amin極大/極小のイベント対応差**とJitterの違い、plateau代表時刻感度。
- **UgイベントとAminイベントの対応**：完全同時性を求めず、同じ周期として対応する率、時刻差の分布、対応不能ケースの説明。
- 全窓**前半`0.50–0.75`、後半`0.75–1.00`**の同条件での集計も出す（物理的安定化原因について考察不要）。
- 面積データが足りない場合も、流量3方式の解析だけは完了し、面積3指標に架空の値を入れない。

## 6. 波形図・目視用成果物

両条件で以下を出す。**全窓図だけで終わらせず、難しい場所の拡大を必ず出力**。

- `flow_threshold_crossings_full.png`：生`Ug`、各パルスの`Upeak`、3つの下降側閾値交差を色/記号で重ね、UNRESOLVED/LOW_CONFIDENCEを影付き表示。
- `flow_threshold_crossings_zoom/`：`0.574–0.600`、`0.614`、`0.7975–0.8000`（r0p25）、`0.549–0.565`、`0.681–0.697`、`0.749–0.767`（r0p30）。波の前後約5周期も含めた拡大図。**0.5%/2%/5%の交差が異なる箇所を矢印表示**。
- `min_area_extrema_full.png`：生`Amin`、極大・極小、plateauの始端/中央/終端、UNRESOLVEDを表示。
- `min_area_extrema_zoom/`：流量側と同じ時刻帯の波形を拡大。**実際にその時刻帯の`Amin`が存在する場合のみ**。
- `flow_area_aligned_panels/`：上下段で同じ実時間軸の`Ug`と`Amin`、イベント対応線（無理に同期させない）。
- `jitter_comparison.png`：r0p25/r0p30×各イベント定義のJitter/RAP/PPQ5。値だけでなく**各方式のJitter評価時間割合を併記**。
- `threshold_event_timeline.png`：全窓のVALID/PROVISIONAL/UNRESOLVED/WINDOW_EDGE、3閾値でどの周期が存在/欠落したか。
- `closure_threshold_sensitivity.png`：0.5,2,5%のイベント数、時刻差、Jitter、閾値依存性。
- `manual_review/`：画像(a)イベントを表示しないblind版と(b)候補・選択理由・採用時刻のdiagnostic版を50–70ms程度の重複パネルに分割し、全0.5秒を覆う。ユーザーが**どの下降交差・面積の山谷を採用したかを直接確認**できるようにする。
- 閾値交差時刻を±数サンプル拡大した図で、**実データの隣接2点と補間時刻**を見せ、計算処理の正当性を確認。
- 注釈CSV（`case_id,signal,event_type,threshold_pct,event_time_s,confidence,reviewer,notes`）と既存`--manual-events`機能を拡張または互換維持。**人手が未注釈なら`PENDING_MANUAL_ANNOTATION`、真値のprecision/recall等を捏造しない**。

## 7. 合成テスト・実データ品質チェック

**実装の新機能を使う前に、次の正解既知の小さなテストを追加（`--self-test`に統合可）**：

1. 単峰＋完全閉鎖の周期一定波`Ug`で、`θ=0.5/2/5%`のイベント数が正しく、`Jitter≈0`。
2. **既知の周期ゆらぎ**を加えた流量波形で、各閾値のJitterが既知のイベント時刻差から得た値と一致。`T0`に吸着してJitterを過小評価しない。
3. **多峰/反跳**：主パルスの下降交差が正しく選択され、同一周期に複数`t_close`を作らない。
4. **流量0のplateau**：0で張り付く区間の中心でなく、**下降側の交差点**を採用。
5. **流量が閾値まで下がらない不完全閉鎖**：`NO_THRESHOLD_CROSSING`となり、極小を偽の通過時刻に置かない。すぐ先の周期では追跡再開。
6. **全体振幅が変化**：局所`Upeak`比の閾値で比較可能なことを確認。微小振幅や分母0の例はQC。
7. **Aminが単峰・多峰・平坦極小・不完全閉鎖**：極大極小の独立追跡、plateau代表時刻の感度。
8. **流量と面積で異なるstride**：時刻・イベント対応に誤差がないことを確認。
9. **評価窓全域の分類**、連続区間が途切れたときも**再開**、`periods`の重複なし、**境界跨ぎJitter項0**、旧`global-cycle`とCPPの再現性。
10. 閾値0.5/2/5%に対し、手動で注釈した交差時刻と±1～2標本精度で一致するようロジック/補間を検証。**この自動テストは実データ人手正解とは区別**。

出力する全`periods`/`events`/`metric_terms`の時刻・イベントIDを別のチェック処理で再計算して確認する。選択と集計の数値の追跡可能性を優先。

## 8. CLI・出力先・ファイル保護

現行`tools/analyze_airflow_metrics.py`を改修の中心とし、入出力部分だけ必要に応じて汎用化する。既存のCLIや`summary.json`の意味を説明なく破壊しない。新規ソースを増やすより同じファイルの関数追加を優先。

想定CLI（実装時は既存の引数命名に合わせてよい）：

~~~bash
python3 tools/analyze_airflow_metrics.py \
  --input "analysis_runs/manual_run_archive/r0p25_P900_20261104T064537_9fba0102/snapshot/dat/airflow_vt.dat" \
  --area-input "/path/to/real/area.dat" \
  --dt 1e-5 --start 0.50 --end 1.00 --case-id r0p25 \
  --tracking-method closure-threshold-and-area \
  --closure-thresholds "0.5,2,5" \
  --area-event-method global-cycle-extrema \
  --out-dir "analysis_runs/airflow_closure_area_jitter_20261107/r0p25"
~~~

`--area-input`は**実在パスに限る**。areaがない場合でも`--input`と`--closure-thresholds`だけで**流量3閾値の評価を完走**できるようにする。r0p30も同様。

今回の結果は`analysis_runs/airflow_closure_area_jitter_20261107/{r0p25,r0p30,combined}`等の**新規out-dirに隔離**。既存の`analysis_runs/airflow_global_cycle_validation_20261106`、旧`airflow_vt.dat`/`area.dat`、baseline、図、manual注釈、ソルバ、ユーザーの未コミット差分を削除・上書き・再計算しない。**今回は追加の`analysis_runs`整理/削除は依頼していない**（前回整理済み）。必要に応じmanifest/元ファイルSHAを出す。

最低限の保存物：

- `flow_threshold_events.csv`（3閾値×全pulse候補、time interpolation、採否・QC）
- `min_area_candidates.csv`・`min_area_events_A_B.csv`・`min_area_plateaus.csv`（areaありの場合）
- `jitter_metrics_by_event_definition.csv`（5方法+旧2方法の統合比較、各指標の有効周期数・項数・分母・使用時間）
- `periods_by_definition.csv`・`metric_terms_by_definition.csv`（全使用時刻と区間、Jitter/RAP/PPQ5項の監査）
- `event_matching_flow_vs_area.csv`・`threshold_crossing_differences.csv`
- `qc_intervals.csv`（全窓の時刻分類・ギャップ・信頼度・理由）
- 図群と`manual_review_index.csv`・空の注釈CSV
- `summary.json` / `analysis_metadata.json`（input SHA、実際のarea列・単位、dt、stride、閾値/方法/コード版）
- `report.md`：**r0p25/r0p30の3流量閾値＋Amin極大/極小の一覧、旧A/Bとの差、評価範囲の途切れ、未評価時間、時刻選択根拠、手動検証が必要な点**を簡潔にまとめる。

## 9. 完了条件とCodexへの最終依頼

- [ ] `r0p25`・`r0p30`の`0.50–1.00s`全窓について、**局所Upeakの0.5/2/5%を下降通過する実時刻**を可能な限り抽出し、Jitter local/RAP/PPQ5を同じ定義で算出する。
- [ ] 最小声門断面積`Amin(t)`が実データにあれば、**極大値・極小値でそれぞれ独立に**周期推定と同じ指標を算出する。**入力列/単位/strideを確定**。
- [ ] 旧方式のように追跡しやすい短い区間だけを主値としない。**全窓を走査し、全有効区間の差項を集計する**。真の欠落/閉鎖なしの区間は定義不能と明示し、**無理に接続しない**。
- [ ] すべてのイベントID・対応する局所ピーク・3閾値交差・Amin極値・plateauを**生波形画像で目視検証できる**。人手未実施なら`PENDING_MANUAL_ANNOTATION`。
- [ ] 同じ共通設定を2条件に適用。新旧Jitterの差だけでなく**イベント時刻/位相と評価時間割合**がどこで変わったかを報告する。
- [ ] 既存の`global-cycle`・CPP・ソルバ・一次データ・前回の解析結果・ユーザーGit作業状態を保護し、合成テストと回帰チェックを通す。

**最終Codex依頼：**

> 今回は剛性比0.25と0.30の既存900Pa結果で、流量`Ug`が**各局所流量ピークの0.5%/2%/5%を下降側に通過する実時刻**を新しい周期イベントとして抽出・追跡し、Jitter local/RAP/PPQ5を比較してください。並行して**最小声門断面積`Amin(t)`の極大・極小イベント**からも同じ指標を計算し、3流量閾値＋2面積イベント＋既存流量A/Bを一つの表にまとめてください。前のように難しい部分で最長区間だけを使うことはせず、**0.50–1.00秒の全時刻を探索**して全ての有効周期を評価し、欠測や不確実なイベントは明示してください。各時刻の局所Upeak・閾値・下降交差・採否をCSVと拡大波形図に示し、面積は極大・極小・plateauの時刻定義を明確にしてください。旧方式とCPP、元データ、ソルバ、ユーザー作業状態は保護。目視図・手動注釈テンプレート・検証レポートまで**一括で作成・実行**してください。**実在しないイベントを補完して全窓を埋めることは禁止**します。
