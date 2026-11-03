# Body剛性比 0.10–0.95 と流量Jitterの関係：M5test/261103固有データの18条件スイープ（Codex実装・解析仕様）

> 作成：2026-11-03 / 対象ローカルリポジトリ：`simulation-main`（対応するGitHubリポジトリ：`KajiyamaYusuke/Mode_superposition_2026`）  
> この文書は**Codexに実装・計算・評価を依頼するための指示**であり、計算を実施済みと表現しない。  
> **調査の主図：横軸＝右Body/左Bodyのヤング率比、縦軸＝声門流量`Ug`から求めるJitter [%]。**  
> **最優先：ローカル`input/m5test/261103`（または実際の大小文字表記）にユーザーが追加した`.txt`/`.vtk`を最初に調査し、ファイル形式と対応を確定する。新データは現時点のGitHub mainでは確認できないため、ファイルが実在することや中身を推測しない。**

## 0. 目的・今回の実験設計の固定事項

これまでBody 12/12, 12/6, 12/4, 12/2 kPaの計算（各Ps=900/1400/1900 Pa）で、左右剛性比と振幅・周期性の関係が単調でないこと、特に12/6・900 Paでは柔らかい右が小振幅になったことを確認した。今回の研究目的は**流量波形の周期のゆらぎが、Body層の剛性比によってどう変化するか**を、系統的な条件群で調べること。

- **左Body層**：`E_B,L = 12 kPa`で固定。
- **右Body層**：`r=E_B,R/E_B,L = 0.10,0.15,...,0.95`の18条件。
- **右Body[kPa]**：`1.2,1.8,2.4,3.0,3.6,4.2,4.8,5.4,6.0,6.6,7.2,7.8,8.4,9.0,9.6,10.2,10.8,11.4`。浮動小数の加算ループで作らず、整数index`k=2..19`、`r=k/20`、`E=12*k/20`（Decimal/整数で厳密に計画）。
- **Cover層**：左右とも`3 kPa`固定を**既存の900 Pa参照計算に合わせる仮設定**とする。**新しいモードファイルが別Cover値を示す場合は勝手に使わない**。
- **左右の形状・密度・モード質量規格化・接触・減衰・流体回路など**：Bodyの右固有データ以外すべて同一。**左Body=12 kPa用の周波数・モードは右条件によらず同一ファイル**。
- **圧力の初期案**：**Ps=900 Pa、静止初期状態から各条件を独立に開始**。まず18条件の**単一圧力系列**を完成させる。以前の900/1400/1900 Paの圧力別3系列に拡張できる設計にするが、**1400/1900 Paを無条件に追加実行して計54本としない**。追加圧力はユーザー判断の`--pressures`で選択。
- **数値基準（既存12/6・900Pa参照に合わせる初期案）**：`dt=1.0e−5 s`、`nstep=50000`（0.5 s）、`nmode=20`/側、`zetaL=zetaR=0.06`、`kc1=5, kc2=6000, kc3=5000`、`ncont=100`、初期間隙0 mm、圧力立上げ0.10 s、`iforce=0`。ただし必ず元の**完成済み**`L12_R6_P900`の`params_used.txt/manifest.txt`と照合し、**`input/param.txt`の現在値を無検証で引き継がない**。
- **Jitter評価主窓**：`[0.30,0.50) s`。立上げ後の`[0.10,0.20)`、`[0.20,0.30)`、`[0.30,0.40)`、`[0.40,0.50)`の窓別結果を補助とする。主窓は直近の振幅比較と同じであり、**本当に定常かどうかは別途検査**。
- **主観測波形**：**`output/dat/airflow_vt.dat`の`Ug [m³/s]`**。左/右の変位・声門面積では**Jitterの主値を求めない**。変位/面積はピークの誤検出判定と機構との関連を見る補助に限る。

### この依頼における「剛性比」の意味

`r`は**設定したBodyヤング率の右/左比**であって、左右の構造固有周波数比、声帯全体の有効剛性比、左右の自励振動周波数比ではない。この3つを図や本文で取り違えない。x軸の表記例：`Body Young's modulus ratio E_R/E_L`。`r=1`の対称参照を追加する場合は、**新規18条件の外側にある参考点**として区別し、同一条件の既存runに限り採用する。

## 1. 最初に行う安全確認・ローカル実態の確認

Codexは以下を、実装や実行の前に必ず確認・報告：

1. `pwd`、`git status --short`、`git branch --show-current`、`git rev-parse HEAD`、現在の`AGENTS.md`、既存のrunや未コミット変更。**git reset/clean/rebase、強制チェックアウト、元runや入力の上書きを禁止**。未コミットの比較器、接触診断、直近の後処理コードを保護する。
2. `find input -iname '*261103*'`等で**新データの正確なディレクトリパスと大小文字**を確定（Linux/WSLでは`m5test`と`M5_test`は異なる）。新しい入力を勝手に`input/M5_test`へ移動しない。
3. 新フォルダの`.txt`、`.vtk`、必要なら`.vtu`をすべて一覧し、**完全なファイル名**、サイズ、SHA256、最初の数行/ヘッダ、個数、同名衝突を記録。**実データ未確認の段階で命名規則・VTK種別を確定したかのように書かない**。
4. 過去の正しい12/6・900Pa run（候補`analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900`）から`txt/params_used.txt`、`txt/manifest.txt`、実行環境・左右入力SHA256を取得する。
5. 元の`leftFrequencyFile`、`leftModeFile`、`leftSurfaceNasFile/rightSurfaceNasFile`が現在も存在するか。**新データは右Body=1.2–11.4 kPa用で、左Body=12 kPaの対応モデルは別の既存ファイルが必要**。見つからなければ右入力だけで実行しない。
6. 本シミュレーションの入力は、`src/simulation/readParams.cpp`の`leftFrequencyFile/rightFrequencyFile/leftModeFile/rightModeFile/leftSurfaceNasFile/rightSurfaceNasFile`による。生成`param.txt`の置き場所が変わると相対パス解釈も変わるので、**必要な実入力パスは絶対パス化**し、元`input/param.txt`を変更しない。
7. 計算本体は`src/simulation/simulation.cpp`、VTU読取は`src/structure/modedata.cpp`。`tools/run_right_stiffness_sweep.py`は**整数kPa・古い`M5_test`・`.vtu`前提**でそのまま使用不能の可能性あり。`tools/compare_stiffness_metrics.py`は**左右変位からJitter**を計算しており、今回は主値を流用しない。`tools/jitter.py`は**声門面積を平滑平均との差で評価するコード**であり、**流量の周期Jitterではない**。共通関数のみ検証して再利用する。

### 1.1 対応表をまず完成（実行禁止ゲート）

ファイル名中の`b`直後がBody層のヤング率、というユーザー指定を基本とする。例`...b6c3...`ならBody=6 kPaと解釈できるか確認する。ただし**`.vtk`と`.txt`の名付け方が完全に同じとは決め打ちしない**。`b`の後に小数点、アンダースコア、`p`等が使われる可能性があるので、**実ファイル一覧から明示的なマッピング規則を確定**し、`b1.2`と`b12`などの誤対応を防止する。

**必須`mode_inventory.csv`列：**
`ratio, body_left_kpa, body_right_kpa, frequency_file, mode_file, actual_extension, detected_format, cover_kpa, point_count, available_modes, first_mode_hz, last_selected_mode_hz, freq_hash, mode_hash, pair_status, notes`。

**18行すべて**に一意の周波数txt＋モードファイルを対応づけ、各`r`について`E_R=12r`を一致確認。欠けた比率・重複・Cover不明・ファイル名に含まれるモデル条件の矛盾は`INPUT_BLOCKED`。推定の最近傍剛性ファイルを代用したり、足りない中間剛性で補間して実行しない。

## 2. **重要：`.vtk`は現行のVTU読取器にそのまま渡せるとは限らない**

現行`ModeData::loadFromVTU`は、**テキストのXML DataArray名**（例えば`Displacement_field,_X-component_@_...`など）を正規表現で読んでいる。**拡張子を`.vtu`へ変えるだけでは、legacy VTKをXML VTUに変換したことにならない**。

前処理ゲート：
1. 各`.vtk`の**実フォーマット**をヘッダで確認。`# vtk DataFile Version`はlegacy、`<VTKFile type="UnstructuredGrid"...>`はXML系。バイナリ/ASCII、データ配列名、浮動小数精度、PointData数、geometry節点数を検査。
2. Legacy VTK等で既存`loadFromVTU`が読み込めない場合は、**元ファイルを変更せず**、`staging/converted_modes/`へ対応する**XML ASCII VTU**を作る読み取り変換器を実装（必要なら`vtk`/`meshio`を利用し、環境に存在する依存関係を先に調査）。**バイナリXML VTUのままでは現行の行テキストパーサが読めない可能性**があることに注意。
3. **単純なファイル名置換・DataArray名の正規化でモードや方向を変えない**。`modedata.cpp`は読み込んだCOMSOL列の`(z,x,y)`への軸順入替を内部で行う。変換時にさらに座標変換しない。PointDataのX/Y/Z成分・周波数・順序・節点IDを**損失なく保持**し、同じ`Geometry`と整合することを検査。
4. 既存`loadFromVTU`はモードの名前内のHzで並べ替え、`loadFreqDamping`は周波数txtの**行順**を使用する。**両者のモード順が一致することが最重要ゲート**。各モードのHzを比較し、ファイル名の周波数とtxtの差が丸め精度に収まること、20モード以上あること、3方向のPointDataがすべて同じ点数・有限値を持つことを確認。**二重固有値やmode重複を雑に並べ替えない**。
5. Txtは現行ローダで`周波数Hz, ω, damping`の3列を要求している（`%`始まり行は読み飛ばす）。COMSOL出力が1列、別列や見出し付きなら元を保持して`staging/normalized_freq/`へ変換し、各変換後行を検算する。`omega≈2πf`も確認。**列値や固有周波数を推測補充しない**。
6. VTKモード数、節点座標・インデックス、固定節点/表面NAS、左右対称の形状対応を点検。モデル剛性だけが異なると期待するが、実際に形状/節点順が変わっているなら**中止**。可視化座標のmm/mとCOMSOL座標は既存正常ケースと比較。
7. `convert_vtk_modes.py`等は**必要な場合のみ**実装し、出力メタ`conversion_manifest.json`に元/変換後hash・配列名/単位・順序検査・最大数値差を記録。変換済みファイルをGitへコミットしない。

**不可欠な事前動作確認**：代表3条件`r=0.10,0.50,0.95`で、読込可能な20モードの有限値・周波数一致・計算初期化が通るかを確認。ソルバでのフルラン前に周波数とmode形状の誤対応を解消する。

## 3. 固定条件を厳密に揃える

### 3.1 参照条件を`12/6・900Pa`からコピー（単一変更原則）

`run_config.json`（新規）を一度作り、`base_case`と参照`params_used.txt`/`manifest.txt`への完全なパスとhashを記載する。**原則、条件間の変更は右の周波数txtと右モードVTK→VTUファイルだけ**。

固定する項目：`leftBody=12, coverL=coverR=3, Ps=900, rho, mu, mass, nu, nsurfz, surface NAS, modeSelection, nmode, dt, nstep, nwrite, zetaL/R, kc1/2/3, ncont, contactReferenceFrequencyHz, initialGapMm, pressureRampTimeSec, vocal tract/subglottal/inlet acoustic profile, separation, viscous pressure, flow options, glottal method profile, AP contact exclusion, writeDeformationVtu`、プログラムrevisionとCMakeオプション。

**絶対に固定しないもの：**右の固有周波数、右のモード形状（新データに依存）。**ファイル名中のBody剛性だけをパラメータ`E`へ入力する機能ではない**ため、右の固有データの入替が正しく反映されたか記録する。

`input/param.txt`にある`nsurfz=18`等のサンプル値や左b4/右b20の入力は、**今回の基準とは異なる可能性が高い**。新runは旧12/6完成runを基準に、実際のgeometryに合う`nsurfz`とモードファイルを明示する。旧runのヘッダ由来の設定は元ファイルを読むまで**確認前**と区別する。

### 3.2 対称点と既存runの再利用

- **`r=1.00`**は18条件に含まれない。過去`L12_R12_P900`が**同一幾何・Cover・接触・減衰・音響・圧力ランプ・計算方法・Jitter定義・評価窓**であれば、**参考の白抜き点**として再解析し追加可能。**新しい対称計算を自動実行しない**。
- 過去`r=0.50`（右Body 6）など既存runがあっても、**新しいCOMSOL固有データのSHAが異なれば別データセット**。その場合は旧runを新系列の代用にしない。
- 完了済み条件をスキップするには、`ratio, Ps, config_hash, input_hashes, executable_hash, complete`の全一致・必要な`airflow_vt.dat`存在・解析QC合格を要求。`status.json`だけで`completed`と判定しない。
- **既存`analysis_runs/`・`output/`・`result/`を絶対に上書きしない**。新`analysis_runs/jitter_body_ratio_261103_.../`下に独立run。容量・計算時間・並列スレッド上限を確認し、同時複数runで出力ディレクトリ競合を起こさない。
- 以前の診断の`FAILED_QC`（フル旧新bit不一致）はこのスイープの全18点へ機械的に伝播させないが、今回の採用ビルドが既存12/6基準を近似再現するかの検査は記録する。

## 4. **流量`Ug`からJitterを測定する方法（最重要）**

### 4.1 入力の事実と単位を固定

現行`src/simulation/simulation.cpp`は`dat/airflow_vt.dat`へ`n%5==0`のとき`step Ug`を出力する（**`nwrite`とは独立**）。`dt=1e−5s`なら実サンプル間隔は`5e−5 s`、サンプル周波数20 kHz。**一列目は時間[s]ではなくstep番号**：
\[
t_i = \mathrm{step}_i\,dt,\qquad U_i=U_g(t_i)\ [\mathrm{m^3/s}].
\]

`numpy.loadtxt(..., comments='#')`等で読み、`step`単調増加、stride、有限値、重複なし、期間、終了時刻を検証。**従来の`tools/jitter.py`の声門最小面積を使わず**、必ず`airflow_vt.dat`の2列目からピーク時刻`t_k`を検出する。流量がすべてほぼ0、非発振、サンプルが欠落した場合は適切にステータス付け。

**注意**：ピーク検出の時間量子化だけで擬似的なJitterが出ることがある。100–200 Hz級の周期が100–200サンプル程度なので、**整数サンプルpeakの結果と3点放物線などのサブサンプル補間結果を両方計算**。補間は明瞭な局所ピークで有限・妥当な補間シフトのときのみ使用し、フラットピーク・段差・クリッピングでは無効フラグを付ける。元データの有効分解能を超える精度を主張しない。

### 4.2 主Jitterの定義を明記

検出した**連続する同種の声門流量周期イベント**の時刻を`t_1,\ldots,t_M`とすると`T_i=t_{i+1}-t_i`、`N=M−1`。主指標は

\[
\mathrm{Jitter}_{local}[\%]
=100\frac{\frac{1}{N-1}\sum_{i=1}^{N-1}|T_{i+1}-T_i|}
{\frac{1}{N}\sum_{i=1}^{N}T_i}.
\]

必要な連続周期数（暫定`N≥10`）を満たさなければ数値を出さず`INSUFFICIENT_CYCLES`。主図は**`Ug`由来の`Jitter_local[%]`**に統一する。

**補助指標**：
- `Jitter_sample_grid[%]`：補間前の整数サンプル時刻のみで同じ計算（サンプリング依存性の見積り）。
- `Jitter_5period[%]`（**参考値であり主指標と混ぜない**）：5周期の局所平均と各`T_i`との差を取り、定義と端点処理を明示する。ユーザーが過去に使った「5周期移動平均」定義に対応するなら、主`Jitter_local`と別系列/別列にする。
- `CV_T=std(T)/mean(T)`、周期列の自己相関、ピーク数・有効周期数、流量振幅とDC、主な周波数ピーク、前後半の比較。

**推奨**：SciPy `find_peaks`等の既存依存を使い、ピークprominence（例えば各窓の`P95−P5`に対して`5%,10%,15%`）とmin peak distance（参考帯域80–400Hzに対して過剰ピークを防ぐ値）を**固定ルールで探索**する。**主結果に最も都合のよい閾値を条件ごとに選ばない**。まず対称12/12と右6の既存`Ug`でピークの妥当性を検証し、主設定を固定して全18条件へ適用。

### 4.3 最大の誤解：高Jitter ≠ 非周期・カオス

- 流量`Ug`は、1振動周期中に複数の局所ピークを持つ場合がある。**生の極大をすべて周期境界にしない**。目視のピーク印付き`Ug`波形、ピーク間隔・基本周期候補（自己相関とスペクトル）、トラフ/閉鎖相を合わせて**同じ種類のパルスを抽出できたか**を検査する。
- 1:2や2:2などの**規則的な周期倍化やサブハーモニクス**では、交互の流量ピーク間隔の変動が大きくても、長い反復周期自体は規則的でありうる。その場合、`Jitter_local`を出すなら**パルス列依存の値**と明示し、`PERIOD_MULTIPLIED`等の状態ラベルも付ける。
- 波形が非定常、発振開始中、周波数ドリフト、ピーク欠落、多峰性、振幅変調、ゼロ流量・未発振、計算停止、幾何不正、異常ピークの場合は、それぞれ`NONSTATIONARY`、`AMBIGUOUS_PEAKS`、`NO_OSCILLATION`、`RUN_FAILED`等の**理由をCSVへ保存**。特に`NO_OSCILLATION`をJitter=0%としない。
- ピーク検出が有効でない条件を**単なる巨大Jitterとして描画しない**。数値的なJitterだけで「非周期振動／カオス」と判定しない。
- 0.30–0.40sと0.40–0.50sの平均周期と`Jitter`、振幅包絡、イベント数を比較して**定常性／時間窓依存性**を判定。最終窓で変化が続くなら`WINDOW_DEPENDENT`として主図で区別し、必要ならより長いrunを**別途提案**（今回無断で延長しない）。
- **感度分析**：prominence設定（3段階）、ピーク補間ON/OFF、平滑化条件（無し／ノイズ対策に必要な場合のみ）、時刻窓を変えて、Jitterの傾向が頑健かを確認。**主設定は前もって固定し、変化が大きい条件にQCフラグ**。測定の推定誤差と計算条件間の変動を混同しない。

### 4.4 メイン図の見せ方

1. **主図**`jitter_vs_body_ratio.png` / `.svg`：x=`E_B,R/E_B,L`（0.10–0.95、0.05刻みの数値軸）、y=`Ug`由来`Jitter_local[%]`。各条件を点表示、**有効・定常な点だけ薄い線で接続**。値がない・QC不合格な点は`×`や下部小帯に**理由付き**で示し、0%へ押し込まない。対称`r=1.0`の既存基準が検証できたら参照の白抜き点として付けてもよい。
2. 別図`jitter_qc_and_cycles.png`：ratio vs peak count、周期CV、状態フラグ。Jitterの極大がピーク誤検出や期間不足でないか確認。
3. 代表波形`flow_examples.png`：最小Jitter・最大Jitter・中間/不良QCの`Ug`（必ず流量[mL/s]への変換`1 m³/s = 10^6 mL/s`を記録）、ピーク位置と周期列。ピークと有効周期を重ね、図だけで算出根拠が追えるようにする。
4. 感度`jitter_sensitivity.png`：各ratioで閾値・補間・前半/後半の値の散らばり。**周期bootstrapを独立なrun間誤差・信頼区間と呼ばない**。
5. 将来圧力を増やしたら、**同一ピーク定義**で`Ps`ごとの線/マーカーを追加し、pressureごとに完走数・QCを別記。900 Paしか実行しない段階では架空の圧力系列を描かない。

## 5. Codexに作ってほしいプログラム構成

`tools/body_ratio_jitter_261103/`を新設するか、既存`tools/run_right_stiffness_sweep.py`を壊さない独立スクリプトとして追加：

- `inventory_modes.py`：フォルダ列挙、Body値パーサ、freq/modeペア・Cover・ハッシュ・20モード・NAS/節点数/周波数検査、`mode_inventory.csv`。
- `convert_modes.py`：入力がLegacy VTK/非対応XMLの場合のみlosslessな変換、DataArray名・座標順・点数・周波数順の検証。元ファイルは読み取り専用。
- `prepare_sweep.py`：参照`params_used.txt`から18条件の専用`param.txt`、`run_config.json`、`cases_manifest.csv`を生成。`--dry-run`をデフォルトにし、18件の変更差分をユーザー/ログへ一覧表示。
- `run_sweep.py`：`--execute`を明示した時だけ開始。`--ratios`/`--pressures`/`--resume`/`--fail-fast`/`--threads`を持ち、独立runごとのログ/終了コード/manifest/hash、**既存データの上書き拒否**。途中終了からのresumeは条件・入力・実行ファイルが完全一致する場合だけ。**ソルバ物理式を変更しない**。
- `analyze_flow_jitter.py`：`dat/airflow_vt.dat`から時間へ変換、`Ug`ピーク・周期・`Jitter_local`、QC、窓・感度解析、`flow_peaks.csv`、`flow_cycle_metrics.csv`、各runの`jitter_metrics.json`。
- `plot_jitter_sweep.py`：全runの`summary.csv`を集約して指定した主図/補助図を作成。未完了runを黙って除外せず凡例・QC表へ残す。
- `tests/`または既存`tools/tests/`へ、ファイル名パーサ、VTK/VTU配列一致、周波数モード順、param diff、Ug時間軸、Jitter精度、ピーク多重/脱落、QC、resume、シミュレーション実行防止のテストを追加。

**流用すべき既存部品**：`tools/dynamics_common.py`の時刻/ピーク/period helper、`tools/left_right_ratios.py`のイベント処理、`tools/analyze_irregularity.py`のピーク可視化とスペクトル、`tools/output_layout.py`。ただし主信号は**必ず`Ug`**にし、旧`tools/jitter.py`や`compare_stiffness_metrics.py`の出力列を「Ug jitter」と名前だけ変更しない。

### 5.1 CLIイメージ（パスの実体を調べた後、Codexが確定）

~~~bash
# まず新しい固有データだけを検証。まだ計算しない。
python tools/body_ratio_jitter_261103/inventory_modes.py \
  --source-dir input/m5test/261103 \
  --left-body-kpa 12 --ratios 0.10:0.95:0.05 \
  --output-dir analysis_runs/jitter_body_ratio_261103_900Pa/inventory

# 固定900Pa・18条件の入力を新規生成（dry-run）
python tools/body_ratio_jitter_261103/prepare_sweep.py \
  --inventory analysis_runs/jitter_body_ratio_261103_900Pa/inventory/mode_inventory.csv \
  --base-run analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900 \
  --pressures 900 --dry-run \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900Pa

# 実行は上記の入力整合とテストをPASS後、まず代表1条件でUg Jitter QC。
# 問題なければ同一固定ルールで残り18条件を順に処理。
python tools/body_ratio_jitter_261103/run_sweep.py \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900Pa \
  --pressures 900 --execute --resume --threads 4

# 条件横断集計／プロット
python tools/body_ratio_jitter_261103/plot_jitter_sweep.py \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900Pa
~~~

**注意**：コマンドはインターフェース仕様例であり、**現時点ではプログラムは存在しない**。`--source-dir`は実際の大小文字に合わせる。**原則として先にdry-run/データ照合を終えてから重い計算**。別圧力実行は`--pressures 900 1400 1900`等を**ユーザーが明示したときのみ**行う。

## 6. 実行順序・停止条件・コスト管理

**Phase 0: Preflight（実行0）**
- Git/入力inventory・本体`modedata.cpp`の読取互換・周波数／Mode順／Geometry・設定比較・旧runの状態を`preflight_report.md`に記録。
- 新モード変換が必要ならそのコードと検証を完了。
- 計画18条件を一覧化し、ratio→`b`値→周波数txt/モードvtk→変換後VTU→生成paramのマッピングが**完全1対1**であることを確認。
- **1つでも`INPUT_BLOCKED`なら全件一括実行に進まない**。ユーザーが追加したデータに欠損がある場合、その場で欠損一覧を報告する。

**Phase 1: 最小動作確認**
- テストで合成`Ug`波形（一定周期Jitter≈0、意図的周期交替、周波数ドリフト、多峰パルス、欠落ピーク、無振動、丸め誤差）を検証。
- 必要なら`--dry-run`のパラメータを使った**短いsmoke**で初期化・出力を検査するが、本体の接触計算/流体の式を変更しない。代表`r=0.10,0.50,0.95`の入力読込が検査できればよい。
- **1条件だけフル`0.50s`計算**（例えば`r=0.50`）。`Ug`ピーク検出を可視化し、主窓の周期抽出・過去の12/6の基本的な波形傾向との一致を確認。ただし新旧モードSHAが異なる場合はbit一致を期待しない。
- スモーク短縮runのJitterを本番Jitterとして混ぜない。

**Phase 2: 900Pa・18条件**
- Phase1が`PASS`なら、すべて**独立初期条件・同一圧力立上げ**で順番に本番計算。失敗時にはログ・停止step・原因を保存し、**隣の剛性値のデータを代入しない**。
- 先に全18本を走らせるのではなく、実行済みケースから順次`Ug`と図・QCを検査。重大なピーク解釈・入力形式の不備があれば**残りを中止**して修正。
- 質量や剛性に依存して出現する大きな動揺、接触異常、幾何のAP順序反転、`NaN`等は記録し、計算が最後まで走っただけで`VALID_JITTER`にしない。
- 最後まで`r`の18点をすべて一覧に残し、`VALID_JITTER / INVALID_JITTER / NO_OSCILLATION / RUN_FAILED`を統一管理。

**Phase 3: 解析・解釈**
- 主図・QC図・代表`Ug`波形・周期系列・感度図・集計CSV・レポートを作成。
- **中心的な問い**：剛性比に対し`Jitter`は単調か／極大や急変があるか／低剛性比で多峰・周期倍化・非定常が多いか。**現象が見つかった条件だけ**次段の精細化やPs拡張を提案する（今回自動実行しない）。
- 既存12/6の局所振幅差と新スイープのJitterを関連づける場合、**異なる固有データ／実行条件／波形期間の混同を避ける**。

## 7. テスト・QCの明示要件

- **ファイル名**：小数`1.2`,`11.4`,整数`3`等を正確にparseし、`b12`との混同なし。重複/欠損にfail closed。
- **VTK/モード**：各ケース3成分、少なくとも20モード、入力周波数列との同順一致、節点数/点順・軸順・幾何整合、有限値、変換前後一致。`.vtk`を単に`.vtu`へrenameするテストは**失敗する**こと。
- **設定**：左入力と音響/接触/時間/初期化条件は18条件で不変。差は右固有データとメタのみ。`Ps`だけは別系列として変化可能、同圧力系列内で固定。
- **波形読取**：`#`コメント/step→s/stride=5/単位`m³/s`、部分run、重複step、NaN、時間軸切れ、ゼロ流量を検出。
- **Jitter数値**：既知のpeak時刻`t_k`から公式が厳密に合う、一定周期≈0（サンプリング補間誤差と区別）、周期交替に高Jitter、欠落ピーク/多峰条件にQCフラグ、最低周期数未満はNaN/invalid、閾値感度・窓依存・サブサンプル比較。
- **再現性**：同一`Ug`を繰返し解析して全CSV数値・図が再現可能、設定JSONと入力SHA固定、スクリプトseed固定（使用する場合）。未測定条件のJitterを0や内挿で埋めない。
- **中止ゲート**：参照run基準が未確定、Cover不一致、左12のファイル不足、入力VTKが正しく読めない、モード周波数不一致、run出力の上書き危険、Jitterピークイベントが物理的に解釈不能、ディスク不足はそれぞれ失敗ステータスとして計算を止める。
- 現在のC++の既知`input_file_params`テストの陳腐化が残っている場合、**本件と無関係な改修で解消しようとしない**。新規のPythonファイル名/解析テストを優先。

## 8. 成果物（新`analysis_runs/jitter_body_ratio_261103_.../`内）

| 成果物 | 中身 |
|---|---|
| `mode_inventory.csv` | 18組のBody剛性、ratio、txt/vtk実ファイル、フォーマット・shape・hash・対応QC |
| `conversion_manifest.json` | Legacy VTK等を変換した場合のlossless検査・入出力SHA・配列/点順 |
| `run_config.json` / `cases_manifest.csv` | 参照12/6の条件、左右入力、各ケースのparams、全固定条件・実行ファイルSHA |
| `cases/rXX_P900/{input,dat,csv,txt,...}` | 1条件1フォルダ、生成param、流量ほか元の通常出力、`stdout.log`・`stderr.log`・`status.json` |
| `jitter_summary.csv` | 全18条件のratio/Body/Ps、`Jitter_local_percent`、補助`Jitter_5period`、周期数・周期CV・定常性・peak QC・状態・失敗理由 |
| `flow_peaks.csv` / `flow_cycle_metrics.csv` | 条件・流量極大時刻・ピーク属性・個別周期・隣接周期変動・採否 |
| `jitter_sensitivity.csv` | prominence・補間有無・評価窓を変えたときの指標の変化、QC |
| `jitter_vs_body_ratio.png/svg` | **主成果：x剛性比、y流量Jitter[%]**。18条件とQC、任意の一致確認済み対称参照 |
| `jitter_qc_and_cycles.png` | 周期数、周期CV、QCステータス |
| `flow_examples.png` | 代表3〜4条件の`Ug`流量波形・検出ピーク・周期 |
| `jitter_sensitivity.png` | 前後半・補間・閾値変化に対する頑健性 |
| `preflight_report.md` / `report.md` | 実装上の対応と入力ゲート／主結果、失敗理由、数値を伴う解釈・限界、次に調べる剛性比領域 |

**Markdownレポート`report.md`の必須構成**：
1. **結論（3～5行）**：何圧で何条件を計算し、どの剛性比で`Ug` Jitterが上昇/低下したか（**実測した値のみ**）。単調性がない場合はその事実。
2. **方法**：独立定圧、左12kPa固定、右18種、新モードデータの読取/変換、`Ug`ピーク検出基準、定義、評価窓、各種QC。
3. **結果**：全18条件表（Jitter[%]、N cycles、期間、period CV、peak QC、fail状況）、図の読み方と代表波形の特徴。
4. **頑健性**：閾値/補間/時間窓での変化と不安定条件。期間不足・非定常条件を数値で明示。
5. **解釈と限界**：「高Jitter」は周期ゆらぎの指標。周期倍化や多峰性・ピーク誤認・圧力立上げ・時間刻み・COMSOLモード固有データ違いを区別し、**カオスの確定とは書かない**。
6. **次に絞る領域**：急増条件の近くの比率/圧力の再計算、長時間評価、別周期検出法、初期値依存の提案を**今回実行せず**示す。

## 9. 実装・実行時の作法

- `AGENTS.md`のC++17,`icpx`+OpenMP、OpenMP配置固定、`-ffast-math`禁止を守る。**本調査ではC++の物理式は原則無変更**。
- `input/m5test/261103`のファイルは読み取り専用。VTK→VTUや周波数正規化は**別staging**。大容量固有モードや生run、生成VTUをgit add/commitしない。
- 18条件の定常圧力系列だけでも計算量は大きい。**並列18本起動を避け、原則1本ずつ**計算。`OMP_NUM_THREADS,OMP_PROC_BIND,OMP_PLACES`とCPU・ディスク・メモリの状態を保存。
- `--dry-run`では本当にソルバを起動しない。`--execute`は明示時だけ。`--resume`は入力/実行バイナリhashの完全照合が条件。中途半端なrunを`completed`にしない。ソルバが`output/`を使う可能性があるので、`SIMULATION_RUN_DIR`等の既存分離機構を必ず使用。
- **解析中にピーク検出ルールを都合よく比率ごとに変えない**。探索用設定と最終報告の主設定を分け、主設定の事前定義を`manifest`に記録する。
- コードの変更ファイルと検証コマンドをまとめ、`git diff --stat`を報告。**Codexはユーザーの許可なくcommit/pushしない**。
- 0.10–0.95全18条件で入力が見つからないなどの問題が生じたら、**計算をしない状態でも`preflight_report.md`と欠損一覧・解決案を残す**。品質不明の点を補って「完成図」とは呼ばない。

## 10. Codexへの最終依頼・完了判定

**依頼**：ローカル`simulation-main`の`input/m5test/261103`内の新固有周波数txt・固有モードvtkを実際に確認し、左Body12kPa固定、右Body比0.10～0.95を0.05刻みとする18条件で、既存3D自励声帯振動ソルバの流量`Ug`を計算してください。まず**Ps=900 Pa**を従来の完成runと同じ条件で固定し、独立runを安全に構築し、ピーク検出QCに合格した`Ug`波形から`Jitter_local[%]`を算出して、剛性比 vs Jitterの図、全条件CSV、代表波形と周期テーブル、感度/QC図、研究的考察を作ってください。実データがLegacy VTKなら検証付きで変換してください。既存run/入力やC++モデル式は変更しないでください。1400/1900 Paなど追加圧力の全条件実行は、実行前にユーザーへ確認してください。

**完了チェック**：
- [ ] `input/m5test/261103`の**実物**で18個のratioをtxt+vtk対応・Body/Cover/節点/20モード/周波数整合まで検証した。
- [ ] 左Body12kPa・Cover3/3・その他固定条件を完成済み900Paケースと照合した。`param.txt`原本を編集していない。
- [ ] `.vtk`形式と現行`loadFromVTU`の互換性を確認し、必要な場合だけlossless変換した。
- [ ] 少なくとも1条件で`Ug`流量のピーク時刻・周期・`Jitter_local`計算を検証した。
- [ ] 全18条件は個別manifest/hash/status付きでrunし、終了した各runの`dat/airflow_vt.dat`を読んだ。未完了/失敗も集計に残した。
- [ ] 主図のx軸は`E_B,R/E_B,L`、y軸は**`Ug`から計算した`Jitter_local [%]`**。無発振や不確実な周期を0にしていない。
- [ ] 時間窓依存・prominence依存・サンプリング補間依存を検証し、周期倍化と非周期性を混同していない。
- [ ] `jitter_summary.csv`、検出周期表、PNG/SVG、`preflight_report.md`、`report.md`を作り、数値をもとにどの剛性比で周期ゆらぎが増えるかを報告した。
- [ ] 元`output/`、`result/`、既存`analysis_runs/`、`input/m5test/261103`は非改変。ソルバ式変更0、無断commit/push 0。

**未確定として明記するもの**：実際の`261103`内の名前とVTKのファイル形式、Cover値、左Body12入力の場所、今回の対象声門下圧の拡張範囲。これらの資料はローカル環境で確認し、**実物が示さない事実を作らない**。
