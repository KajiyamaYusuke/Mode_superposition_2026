# Codex実装指示書：対称声帯の流量から Jitter / RAP / PPQ5 / CPP を計算する単一Pythonツール

> 作成：2026-11-05  
> 対象：ローカル作業ディレクトリ `simulation-main` / リポジトリ `KajiyamaYusuke/Mode_superposition_2026`  
> **今回の作業は既存runの後処理だけ**。COMSOL、固有モード、C++ソルバ、物理条件の変更、シミュレーション再実行はしない。  
> **原則として新規実装は `tools/analyze_airflow_metrics.py` の1ファイルで完結**させる。出力結果（CSV/JSON/PNG）以外に散らばった新規ソースや設定ファイルを増やさない。  
> 対象ファイルはまず手元にある**対称条件の実測（シミュレーション）声門流量**。最終的には任意の `*.dat` をコマンドラインで指定でき、非対称runにも同じコードを適用できること。

## 1. 目的と範囲（必読）

今回の目的は、声帯シミュレーションが出力する声門流量 `Ug(t)` の周期性を、**Jitter local / RAP / PPQ5 / CPP** の4つで評価し、各値の**信頼性（検出率・閾値感度・サブサンプル推定・サンプリング間隔）**を可視化・定量化することである。

- **最初に実施**：既存の**対称条件**（左右同じBody/Cover）の`airflow_vt.dat`をユーザーの指定パスから読み、1条件の4指標・QC・図を作る。左/右それぞれの流量ではなく、**1つの声門流量Ug**を解析。
- **今回やらない**：未完成の剛性比36条件runの起動、仮のJitterを使った科学的結論、任意の流量ファイルの自動探索、入力データの書換え。
- **主出力**は解析窓ごとに4指標とQC。**Jitter/RAP/PPQ5の周期イベントが怪しい場合も値を生成可能な範囲で算出したうえで暫定/不合格と明記**。完全な周期不明は無理に0や補間で埋めない。**CPPは周期イベント時刻が不明でも独立に算出できる**。
- ソースにある`tools/jitter.py`は実際には声門面積の振幅点列への5点移動平均差を処理しており、**正しい流量周期Jitterではない**。`tools/compare_stiffness_metrics.py`は変位波形由来。流用しない。`tools/dynamics_common.py`、`tools/analyze_irregularity.py`、`tools/output_layout.py`にピーク/時刻/出力の参考実装があるが、**依存モジュールを新増設せず1ファイル内で必要最小限の処理を実装**してよい。

## 2. 入力仕様・ファイル読込

ユーザーの既存`.dat`は次のようなASCII 2列形式。`#`行はコメントであり、**1列目は時間[s]ではなく計算step番号**。

~~~text
# step  airflow[m^3/s]
0 0
5 0
10 0
15 0
...
99875 5.34497e-05
99880 5.83582e-05
...
99995 0.000305447
~~~

- `time_s = step * solver_dt_s`。例：`dt=1.0e-5 s`、step差5なら**有効標本間隔`Δt=5e−5s`、fs=20000 Hz**。先頭の`step=0`と初期流量0は正常（ただし振動前の静止区間はJitterの評価窓から除く）。
- `--input <path-to-airflow_vt.dat>`は必須。**パス/ファイル名に固定文字列を組み込まない**。絶対・相対パス・空白を含むパスに対応。`--dt 1e-5`指定可能。`--manifest path`を与えたときは既存manifestの`dt_s`から推定可能とする。dtが特定できなければ**推測せず入力エラー**。両方指定されたら一致確認し、異なる場合は警告のみで続行せず**不一致として停止**する。
- `numpy.loadtxt`または同等の堅牢な読取を使い、`#`、空行、数値科学表記、2列、有限値、stepの単調増加、重複、等間隔性（末尾だけの例外を含め厳格記録）を検証。任意の抜け行を**勝手に等間隔へ補間しない**。
- 読込後に`step_stride`、`dt_s`、`sampling_interval_s`、`sample_rate_hz`、総長、window内サンプル数、min/max、静止区間有無を`analysis_metadata.json`へ記録。ファイルSHA256も記録。
- **時間窓はCLI指定**：`--start 0.30 --end 0.50`や`--start 0.50 --end 1.00`。1s/0.5sの既存データの両方を想定し、**固定の0.5–1.0秒を勝手に適用しない**。境界`[start,end)`（左閉右開）。省略時は全区間解析可だが、**初期圧力ランプが含まれJitterに不適切となる可能性を警告**し、主評価としてはCLIで定常窓明示を推奨。
- 評価窓だけのピーク検出による端点欠落を抑えるため、前後に数周期分のpaddingを持って候補抽出し、**完成した周期区間が窓内に含まれるものだけを主計算に採用**。window外イベントはCSVに区別して保存。
- `--out-dir`で出力先を指定。**入力ファイルと同じパスへ上書きしない**。デフォルト出力先は明示された入力名に基づく安全な`analysis_results/<stem>_...`等でよいが、既存同名出力があるなら上書き拒否（`--overwrite`が明示された場合のみ生成結果を上書き、**元`.dat`は絶対に保護**）。
- 各図の流量表示は`mL/s`（`1 m³/s=10^6 mL/s`）でもよいが、計算内部はm³/sで保持。時間単位もs/msの混同禁止。

## 3. **イベント検出：自己相関は周期候補、最終Jitterは実際の時刻差**

**最重要**：自己相関で求めた`T0`だけでは**Jitterは算出できない**。`Jitter/RAP/PPQ5`は`Ug`波形上の**同種の周期イベント`t_i`**を追跡し、`T_i=t_{i+1}-t_i`から計算する。**単純な全極大値間隔は多峰波形で誤る**。Praat等の手法を参考にするが、実装をPraatそのものと偽らない。

### 3.1 周期候補推定（ACFとスペクトル）

1. 窓内の流量波形からDCを除く。過度な移動平均や強い平滑化で本来の不規則性を消さない。窓ごとにゼロ流量やフラットな信号は適切に`NO_OSCILLATION`。
2. 正規化自己相関`R(τ)`を求める。**ラグによって重なるサンプルが減るため、重なり部分の分散やエネルギーで正規化**。τ=0以外の局所極大候補を取得。
3. `--fmin/--fmax`（初期例50–500Hz、固定の人声正常値に拘束しない）から`τ∈[1/fmax,1/fmin]`を探索。`T0/2`、`T0`、`2T0`、サブハーモニクスに対応する候補を、波形類似性とスペクトルの根拠とともに列挙。**最大ACFの1候補を絶対真値としない**。`pitch_ambiguity`フラグを用意。
4. F0候補とその根拠（ACF lag/強度・スペクトルの倍音/分数倍の位置）をCSV/JSONへ。周期倍化/非定常性の疑いがある場合は主周期候補を選ぶ根拠を可視化。
5. 局所時間窓のACF（例0.08–0.12 s、設定可能）で時間変動を検出するのは補助。**局所窓ごとに`T0`を強制補正してJitterを平滑化しない**。

### 3.2 主イベント追跡法

- 主法A：`scipy.signal.find_peaks(Ug)`でprominenceベースに候補抽出（後述の感度分析）。`distance`は`T0`の候補と関連づけるが、**固定距離で実際の周期変動や倍周期を無理に消さない**。
- 大/小の多峰が混在する場合、ACFの`T0`近傍の次の同種パルスを探索し、**局所形状の相互相関（windowed waveform matching）やピークの前後の立上り・下降相**を使って系列をつなぐ。対応候補を複数保持し、各周期の`time`と`match_score`を返す。
- 波形全体の形状が変化しても一貫した位相に対応する時刻を選ぶ。**単にt_i+T0に近いピークを常に採用しない**。周期長の急変をJitter計算から恣意的に削除しない。
- 対照法B：流量最小値、急峻な立上り/下降相など、**異なる位相基準のイベント**で周期列を別計算する。主法と異なる場合は原因を`event_alignment.csv`に残す。検出法Bの時刻をAの欠落点に無断で補填しない。
- ±1サンプル程度の丸め誤差を減らすため`subsample=quadratic`（ピーク近傍3点放物線）を実装。凹曲率が足りない・平坦・境界外なら補間を無効として記録。**単なる線形補間/ゼロ埋めアップサンプリングを「高分解能計測」と呼ばない**。
- デフォルトは、**明らかな周期外れ値でも消さずに生のイベント列と周期をCSV保存**。period floor/ceiling、隣接周期比（Praatでの`maximum period factor`）による除外・トリミングをするなら、**その条件と除外周期数、除外前後の各指標を明示**。1.3倍等を無言のデフォルトフィルタにしない。
- **同じ`t_i`をJitter local/RAP/PPQ5で共有**。分析者に都合よい指標ごとの再検出や窓切替は禁止。
- 推定`T0`と追跡イベント間隔が整合しなくても、無理やりイベントを増減しない。`PERIOD_AMBIGUOUS`/`DOUBLE_PEAKS`/`MISSING_EVENTS`/`PERIOD_MULTIPLIED_POSSIBLE`等で明示し、**Jitterの暫定値と信頼性を切り離す**。

## 4. 指標の計算式・単位と欠測規則

有効な同種イベント時刻`t_1,\ldots,t_M`と連続する`N=M−1`周期`T_i=t_{i+1}−t_i`、`\overline T=mean(T_i)`を用いる。

### 4.1 Jitter local [%]（主指標）

\[
J_{\mathrm{local}}[\%]
=100\frac{\frac{1}{N-1}\sum_{i=1}^{N-1}|T_{i+1}-T_i|}{\overline T}.
\]

### 4.2 RAP [%] (Relative Average Perturbation)

\[
\mathrm{RAP}[\%]
=100\frac{\frac{1}{N-2}\sum_{i=2}^{N-1}\left|T_i-\frac{T_{i-1}+T_i+T_{i+1}}{3}\right|}{\overline T}.
\]

### 4.3 PPQ5 [%] (five-point Period Perturbation Quotient)

\[
\mathrm{PPQ5}[\%]
=100\frac{\frac{1}{N-4}\sum_{i=3}^{N-2}\left|T_i-\frac{T_{i-2}+T_{i-1}+T_i+T_{i+1}+T_{i+2}}{5}\right|}{\overline T}.
\]

- 指標に必要な最低周期数はそれぞれ`J_local N≥2, RAP N≥3, PPQ5 N≥5`。ただし**研究上の主評価QC**はたとえば`N≥10`を別条件として設ける（設定可能）。少数でも計算式上可能な値と信頼できる値を区別。
- 分母`\overline T`、採用周期数、周期CV、波形からのmedian period/median f0もメタに保存。
- 3指標は**同じ流量イベントからの統計量**であり、RAP/PPQ5はJitterのピーク誤判定を修正しない。標準病的閾値を流量に適用しない。
- フィルタあり/なし、補間あり/なし、閾値違いは**別解析条件**として保存。欠測は`null/NaN+status+reason`、0にしない。

### 4.4 CPP [dB] (Cepstral Peak Prominence)

**CPPは声門流量の「周期構造の明瞭さ」を見る別指標**。個々のイベント`t_i`の検出は不要。ただし**録音音声の臨床CPPと数値互換であると主張しない**。今回要求は`CPP`であり、時間/ケフレンシー平滑化を伴う`CPPS`とは区別。実装ラベルは`flow_CPP_dB`。

**必須のアルゴリズム仕様をコード・出力metadataに固定**：
1. `Ug`の評価窓を重複フレームへ分割（初期例`frame_length=0.060 s`、`hop=0.005 s`、Hann窓、window単位のDC除去、必要なら線形トレンド除去）。主設定をCLIで上書き可。
2. フレームの**対数パワースペクトル**を求め、逆フーリエ変換により実ケプストラム`c(q)`を取得。**power cepstrum（例：`|c(q)|²`）をdB表記**する場合、`10log10`する。`rfft`と`irfft`の向き、FFT長、周波数帯域（例上限5kHz、処理方法）、窓、dB基準、`epsilon`処理を固定して実装し、方法を`metadata`に全記録。
3. `quefrency q`の探索域は`[1/fmax,1/fmin]`（例50–500 Hz→2–20 ms）。複数ピーク、Tと2Tの候補を保持し、**CPPピーク位置が主F0に対応していなければflag**。CPPはイベント追跡が不合格でも計算可能だが、ピーク位置の信頼性は別に評価。
4. ケプストラムの背景傾向（`q≥1ms`～上限など、固定した範囲）を**回帰線**で近似。必要なら`scipy.stats.theilslopes`等のロバスト回帰を検討し、**ピークの影響を受けすぎない背景**にする。実ケプストラムと回帰線の縦軸は共にdBにする。
5. 主要ピークの高さ[dB]と**同じquefrencyでの背景回帰線[dB]の差**を`CPP_frame_dB`とし、フレーム値の平均（中央値・標準偏差も補助）を`flow_CPP_dB`として出す。**ピークの絶対dB値そのものをCPPと誤称しない**。
6. **`cepstrum_example.png`で実際のケプストラム曲線、探索域、選択されたピーク、背景回帰線、CPPの垂直差、推定`1/q`を図示**。ピッチ分数倍の誤認、窓長による差を確認。
7. `frame_length`が極端に短く基本周期を数回含まない、音量ゼロ、低帯域優勢、数値劣化、frame不足、FFTナイキスト制限などをQCで検出。CPPも`NO_SIGNAL/INSUFFICIENT_FRAMES/CEPSTRAL_PEAK_AMBIGUOUS`等を持つ。
8. **コードで実際に採用した数学的定義・PowerCepstrum/単位/dBの根拠**をmetadataへ書く。PraatのCPPSに近い数値を名乗らない。

参考：
- Praat CPP/CPPSマニュアル：`https://praat.org/manual/PowerCepstrogram__Get_CPPS___.html`
- Praat smoothing：`https://praat.org/manual/PowerCepstrogram__Smooth___.html`
- Praat Jitter local：`https://fon.hum.uva.nl/praat/manual/PointProcess__Get_jitter__local____.html`
- Praat PPQ5：`https://www.fon.hum.uva.nl/praat/manual/PointProcess__Get_jitter__ppq5____.html`

**注意**：`CPP`は算出手順の違いで絶対値が変わる。上の定義を固定し、対称条件で得られた値の相対比較から始める。臨床基準値やPraat/ADSV値と直接比較しない。

## 5. **信頼性評価（ユーザー要求の4項目）**

### 5.1 基本周期検出率／追跡率

**何を分母にしたか必ず記載すること。実際の「正解周期数」が未知なのに真の検出率と呼ばない**。

- ACF主周期候補から`N_expected≈analysis_window_length/T_ACF`を推定。`N_detected`＝追跡した有効完全周期数。暫定的な`coverage_proxy=N_detected/N_expected`を記録する。ただし周期倍化・非定常・ACF候補間違いでは比率が100%超や過小値になり得るため、**これは真の検出率ではない**。
- 検出の厳密なQCとして、**連続追跡に成功した時間長/全窓長**の`covered_time_fraction`、平均波形相関/マッチスコア、短周期/長周期の異常・重複候補・未採用ピーク数を記録。
- 人間が正解時刻を与えた場合だけ`--manual-events path.csv`で真のprecision/recall（許容時間誤差明示、誤差分布）を算出できる**任意オプション**。必須ではない。**正解ラベルなしの「検出率」は推定proxy**であることを図にも記載。
- ACF`T0`と中央値`median(T_i)`の比、`T/2`と`2T`の曖昧さを可視化。
- 規則的な対称条件ではイベント対応が連続しているかを主に見る。**対称だからJitter=0と決めつけない**。

### 5.2 閾値によるJitterのばらつき

- 候補`--prominences 3 5 8 10 12 15 20` （単位：`(P95−P5)`で正規化した%）。**JitterのみでなくRAP/PPQ5とイベント数・coverage・ACF整合・QCも各閾値で記録**。
- **主設定は事前固定（初期値10%）**。最良Jitterや最小Jitterを与える閾値を結果を見て主値に選ばない。
- `threshold_sensitivity.csv`に閾値別`jit_local_percent,rap_percent,ppq5_percent,n_cycles,coverage_proxy,covered_time_fraction,acf_consistency,flags`等。
- 信頼性要約：同一窓内の閾値別Jitterの`median,min,max,IQR,CV`（無効閾値は統計から除外し、その数を必ず報告）。**これは方法選択に対する感度であり、統計学的CIではない**。
- 閾値間イベントの時刻対応率（許容ずれのサンプル数またはT0比率を明示）も計算。**Jitter値が似ていてもイベントが全然違う場合に見逃さない**。
- 図`threshold_sensitivity.png`はx=prominence[%]、y=Jitter/RAP/PPQ5[%]、系列は主法/対照法など。上下にサンプル数・QCを併記。

### 5.3 サブサンプル補間 ON/OFF

- 同じ候補・同じ採用イベントIDに対して、**整数サンプル時刻 vs 補間時刻**を比較し、`jit_local_percent`、`RAP`、`PPQ5`、平均周期、時刻オフセット（µs）の変化を記録。
- **補間のON/OFFでイベントの選択自体が変わったら、その差を補間による精度差と混同しない**。まず同一イベント対応で比較する。候補が変わる場合は別途`event_identity_changed`。
- 20kHz（50 µs）程度のデータで200Hz周期5msなら1サンプルは周期の1%。補間が数字をなめらかにしても**正しい真値に近づいた保証はない**。Jitterの評価時には`timestamp_resolution_us`と`median_samples_per_period`を示す。
- `subsample_comparison.csv`、`subsample_comparison.png`を保存。

### 5.4 1ステップ vs 5ステップ出力間隔（**データが存在する場合のみ**）

- まず入力`.dat`の`step`差を自動判定。**5ステップ間隔しかないなら1ステップ波形は復元不能**。ゼロ詰め・補間upsamplingを「1stepの測定値」として比較しない。
- `--reference-fullrate path/to/airflow_stride1.dat`（または`--compare-strides`）の**独立した実際の1step流量データ**がある場合のみ厳密な比較を行う。該当データがなければ`NOT_AVAILABLE`をCSV/JSONに記載し、解析全体は成功させる。
- 1step参照がある場合、`dt`・時刻原点・入力の同一条件/同一状態が一致するか確認。**条件の異なる2本のsimulation出力をサンプリング比較として扱わない**。可能であれば**1step出力という同じ波形から`phase=0,1,2,3,4`で5通りに5step間引き**、時間分解能のみによる差を評価。独立に取得した5step版と数値が合うかも確認。
- 各rateで主法のイベント対応、`J_local/RAP/PPQ5`、CPP（フレーム長は秒で同一、Nyquistと上限周波数を揃える）、補間ON/OFF、F0、検出率proxyを比較。**phase依存のばらつきを数値と図で提示**。
- `sampling_comparison.csv`と`sampling_comparison.png`を保存。1stepデータがない場合は**省略を明示**し、将来の出力頻度変更は別計画として提案のみ。

## 6. 必須の図・成果物（すべて1回のCLIから生成）

出力先`--out-dir`以下に作成。**本体コードは1ファイル**、CSV/JSON/PNGは実行結果であり追加ソースではない。

| ファイル名 | 内容 |
|---|---|
| `summary.json` | 4指標、時間窓・単位・検出法/閾値・ACF周期候補・QCフラグ、入力SHA、計算式とCPP実装定義 |
| `summary.csv` | `jitter_local_percent,rap_percent,ppq5_percent,flow_CPP_dB`、status、選択窓、period/フレーム数 |
| `events.csv` | 各イベント`t_i`、sample ID、subsample時刻、prominence、`match_score`、検出法、採否・理由 |
| `periods.csv` | 周期`T_i`、時刻端点、イベントID、`ΔT_i`、外れ値候補・欠落疑い |
| `threshold_sensitivity.csv` | 7閾値ごとのJitter/RAP/PPQ5・イベント数とQC |
| `subsample_comparison.csv` | ON/OFFの3指標・時刻差・同一イベントかどうか |
| `sampling_comparison.csv` | 1stepがある場合のみ（なければ`NOT_AVAILABLE`） |
| `flow_overview.png` | `Ug`全窓とイベント位置・評価区間。静止～ランプ部分もあれば区別 |
| **`flow_events_zoom.png`** | **人間が視認できる数周期（例50ms）を拡大し、採用ピーク番号・周期`T_i`・棄却候補を描く**。これが最も重要なQC図 |
| `autocorrelation_period_candidates.png` | ACFのピークと基本周期候補 `T0/2`, `T0`, `2T0` を、対応周波数・強度とともに示す |
| `cepstrum_example.png` | **ケプストラム・背景回帰線・CPPの高さ・ピーク`q`/`f`** |
| `cpp_frames.png` | 時間フレームごとのCPPとピーク周波数。定常性 |
| `threshold_sensitivity.png` | 閾値変化による周期数/Jitter/RAP/PPQ5・QC |
| `subsample_comparison.png` | サブサンプル補間ON/OFFの差 |
| `sampling_comparison.png` | 1step参照があれば5step間引き比較 |
| `report.md` | 結果要約、計算条件、図の読み方、無発振・曖昧波形の扱い、信頼性と限界 |

注：`autocorrelation_period_candidates.png`では自己相関の主周期候補だけでなく半分・2倍候補が実際に探索範囲外ならその旨を記す。必要な図の枚数を減らしたい場合は**意味の異なるグラフを一枚に詰め込まず**、いくつかの図を省いてもよいが**`flow_events_zoom.png`と`cepstrum_example.png`は必須**。

**例示可視化の要求**：
- ピーク採用点（丸）と棄却点（×）、検出法A/Bのマーカー、自己相関推定周期`T0`の目安区間を波形に重ねる。
- PPQ5に使った**連続5周期がどこに存在するか**をイベント番号・範囲で示す例を作る。Jitter/RAP/PPQ5の計算に使う周期が肉眼で確かめられること。
- CPP図は実際のデータから作り、架空の概念図で代用しない。
- 背景やグラフ装飾より、**どのデータ点を使って何を計算したかの検証可能性**を優先。

## 7. 実装インターフェースと単一ファイル要件

**新規コード**：`tools/analyze_airflow_metrics.py`。最低限`argparse, pathlib, csv, json, hashlib, numpy, scipy, matplotlib`を使用可能。既存依存を優先し、不要な`pip install`や新規`requirements`ファイルを増やさない。**ソルバのC++/入力ファイルを編集しない**。

`--input`、`--dt`（または`--manifest`）、`--start`、`--end`、`--out-dir`、`--prominences`、`--fmin`、`--fmax`、`--primary-prominence`、`--min-cycles`、`--zoom-start`、`--zoom-duration`、`--cpp-frame-length`、`--cpp-hop`、`--reference-fullrate`、`--overwrite`、`--self-test`を用意する。すべて初期値と単位をヘルプに明記。必要以上の多数のCLIオプションが複雑になるなら、低頻度に使う項目は1ファイル内の`DEFAULT_CONFIG`にまとめてよいが**`--input`と時間窓、output先は常にCLIで指定可能**。

**CLI例（実行済みではなく、Codexが実装する仕様）**：

~~~bash
# 既存対称run（0.5秒）の解析例。パスは実際の対称ケースに差し替える。
python tools/analyze_airflow_metrics.py \
  --input "/absolute/path/to/symmetric_case/output/dat/airflow_vt.dat" \
  --dt 1e-5 \
  --start 0.30 --end 0.50 \
  --out-dir "/absolute/path/to/airflow_metrics_symmetric"

# 同じ解析を1秒runの後半へ
python tools/analyze_airflow_metrics.py \
  --input "/absolute/path/to/another_symmetric_case/dat/airflow_vt.dat" \
  --dt 1e-5 \
  --start 0.50 --end 1.00 \
  --prominences 3 5 8 10 12 15 20 \
  --out-dir "/absolute/path/to/another_symmetric_metrics"

# もし実測の1step出力が別に存在するときだけ
python tools/analyze_airflow_metrics.py \
  --input "/path/to/airflow_stride5.dat" --dt 1e-5 \
  --reference-fullrate "/path/to/airflow_stride1.dat" \
  --start 0.30 --end 0.50 \
  --out-dir "/path/to/compare"

# 合成データテスト（入力ファイルを必要とせず外部の既存runを変えない）
python tools/analyze_airflow_metrics.py --self-test
~~~

**処理フロー**：read/validate → time window → ACF/F0 candidates → event extraction A/B → subsample ON/OFF → Jitter/RAP/PPQ5 → CPP frames → threshold QC → sampling comparison if available → CSV/JSON → PNG → report。

複数`--input`を同時に受け付けることは今回の必須ではない（単一波形で完成させてから拡張）。余裕があれば対称runを別入力として逐次実行できる小さな`--batch`を追加してよいが**別Pythonファイルに分割しない**。

## 8. **テスト仕様：`--self-test`を同じPython内に入れる**

少なくとも以下の合成・例外条件を**再現可能な固定乱数seed**で生成し、`assert`や`unittest`で合格しなければ完成としない。**合成波形のJitter真値は「ピーク時刻が既知」の場合、既知イベント列の式でも別計算して照合**。

1. 一定周期・単峰の合成パルス：Jitter/RAP/PPQ5≈0、ACFが既知`T`を推定、各イベントが1周期1件。サブサンプル推定の誤差を記録。
2. 意図的な既知周期変動：`T_i`を直接構成してJitter/RAP/PPQ5の式値に一致。**真のJitterを閾値処理で人工的に低下させていないか**。
3. **各周期に大小2峰があるが全体は周期的**：検出周期数を2倍にせず、正しい基本周期を採用し、ピークの取り違えを`DOUBLE_PEAKS`として発見する。
4. 周期倍化`T,2T`/交互振幅：ACFの倍周期曖昧性を報告し、無根拠に一定周期化しない。
5. 連続流量トレンド、最初の静止・ランプ、有限の振幅変調、波形形状の緩やかな変化、実波形の多峰ノイズ。
6. 欠落ピーク・追加偽ピーク・ゼロ波形・NaN・重複step・不等間隔step・単なる短い窓。**誤周期を正常として報告しない**。
7. 時間間隔1step/5stepの合成参照と5つの間引きphase。元`5step`しかないときは`NOT_AVAILABLE`で成功。
8. ピーク時刻の補間ON/OFFで同一イベントのまま統計を比較。
9. CPP：単一周波数または高調波を含む周期波形でケプストラムのピーク`q≈T`に出ること、ホワイトノイズやゼロ信号との違い、CPPのdB単位・窓長/epsilon/線形スケーリングに対する頑健性。**CPPの閾値を医学的正常値としてテストしない**。
10. 2列以外の入力、コメント、指数表記、step stride5、`--dt`未指定/manifest不在、`--dt`不一致、既存出力ディレクトリの上書き拒否、標準出力の指標報告のテスト。

可能なときは**既存の対称波形1つに対して可視化検証を行う**。`--self-test`がPASSしても、**実波形の周期イベントが目視と一致しているか**は独立に評価する。目視正解の入力は任意。

## 9. 作業環境の清潔さ・保護ルール

**最初に**`pwd`、`git status --short`、ブランチ、作業中差分、`AGENTS.md`、既存Python環境/依存の有無を確認。**ユーザーの未コミット変更を削除しない**。既存のC++、`tools/jitter.py`、`input/param.txt`、既存`analysis_runs`や`result/output`は不変にする。

- 追加する新規**ソースファイルは基本1つ**：`tools/analyze_airflow_metrics.py`。テストはその末尾の`--self-test`に内蔵。新規ドキュメントの重複や新しいモジュール・依存パッケージは極力禁止。
- 入力`.dat`やrun内部にレポートを作らず、全結果を`--out-dir`以下にまとめる。暫定ファイルや生成png/csvは`.gitignore`対象の既存`analysis_runs`/任意の外部ディレクトリなどへ配置する。入力を上書きしない。
- Pythonのキャッシュ、テストの一時生成物は`tempfile.TemporaryDirectory`等の一時領域を使い、正常時・例外時とも後片付け。汚れていた既存worktreeを勝手に`git clean/reset`で消さない。
- 物理式の変更・再ビルド・高コストのsimulation実行は必要なし。解析だけで結果を出す。
- 実装後に`git status --short`と`git diff --stat`を確認し、新規単一スクリプトと意図した出力以外に変化がないことを報告。**既存差分があれば維持**。Codexの実装commit/pushはユーザーの指示に従う（この指示書そのものをリポジトリへ保存する操作とは別）。
- 既存の1秒スイープ用計画書`docs/plans/JITTER_VS_BODY_STIFFNESS_RATIO_M5TEST_261103_20261103.md`を**書き換えず**、将来この独立ツールを36条件へ適用できるようにする。

## 10. 成果報告の書式と完了判定

**Codexには、まずローカルの対称`airflow_vt.dat`を実際に1つユーザーが指定しているなら解析**してほしい。ファイルパス未指定なら、**ツールを完成させCLI使い方を提示し、架空の結果を出さない**。

`report.md`の必須構成：
- ファイル（basename・hash・step stride・dt・選択時間窓・流量レンジ・ACF周期候補）。
- **Jitter local / RAP / PPQ5 [%]、CPP [dB]、周期数、CPPフレーム数**。各数値のQC状態と算出方法。
- **基本周期追跡率**：推定expected、追跡actual、coverage、ACF整合、period doubling判定、注意点。
- **閾値感度**：閾値の全値・検出イベント数・採用率、Jitter中央値/範囲/IQR、主設定の値。
- **サブサンプルON/OFF**：同じイベント集合での3指標の差、検出時刻差。
- **1step/5step**：比較できた場合は全phaseの値。できなければ`NOT_AVAILABLE`とその理由。
- **図**：流量の採用イベント／周期に実際に使った領域、自己相関、CPP背景線/ピーク。
- **研究上の解釈**：対称モデルでも完全周期とは限らない。自励振動のJitter、検出器の揺れ、サンプリング誤差を区別し、流量CPPの臨床的意味を過大主張しない。

**最低完了条件**：
- [ ] 追加ソースが原則`tools/analyze_airflow_metrics.py`の1ファイルのみ。
- [ ] 任意パスの`# step airflow[m^3/s]`付き2列`.dat`を正確に読める。
- [ ] `dt`・strideが正しく時間軸に反映される。定常評価窓をユーザーが指定可能。
- [ ] 4指標を算出：`Jitter local/RAP/PPQ5`は**イベント間隔**から、`CPP`は**対数スペクトル→ケプストラム→背景回帰線との差**から。
- [ ] ACF・波形照合で多峰波形の基本周期誤認を抑え、あいまいな周期は`QC`で警告。Rawイベントを必ず残す。
- [ ] **基本周期追跡率の定義と分母**が明示され、正解未知時はproxyと呼ぶ。
- [ ] 閾値7設定の指標値・周期数・QCとJitter感度を保存。
- [ ] サブサンプル補間OFF/ONを**同一イベント集合**で比較。
- [ ] 1step真値データがあれば5phase間引き比較。なければ復元したふりをせず`NOT_AVAILABLE`。
- [ ] `flow_events_zoom.png`（周期番号の可視化）、`cepstrum_example.png`（CPPの差の可視化）などから算出根拠を確認できる。
- [ ] `--self-test`がPASSし、既存ソース・入力/出力・未コミット変更を保全。新しい本体ソース/依存/作業残骸を散らさない。

**Codexへの最終指示**：上記を実装し、対称の`.dat`1件で実際の4指標と信頼性を検証する。現時点で具体的なローカル対称流量パスが指定されていなければ、それを決め打ちせず**スクリプト・合成テスト・実行例まで仕上げる**。特に検出イベントを可視化し、**Jitterの値だけでなく「どの周期から求めたか」「その値をどこまで信用できるか」を必ず数値と図で残す**。
