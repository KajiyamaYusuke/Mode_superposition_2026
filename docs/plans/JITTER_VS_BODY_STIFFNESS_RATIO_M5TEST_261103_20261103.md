# Body剛性比0.10–0.95と流量Jitter：261103固有モードによる900/1400 Pa・各1秒の36条件スイープ

> 2026-11-03 更新版 / Codexへの実装・計算・解析指示書  
> 対象ローカルリポジトリ：`simulation-main`。本書は連携GitHub `KajiyamaYusuke/Mode_superposition_2026` の`main`に保存する。  
> **確定事項：剛性比18点 × 声門下圧2点（900、1400 Pa）＝36条件。各1.0 s。Jitterの主データは声門流量`Ug`。**  
> **ユーザー確認済み命名：`b1p2c3.txt`はBody=1.2 kPa、Cover=3 kPa（`p`は小数点）。`b11p4c3`はBody=11.4 kPaと解釈。**  
> 実行前にローカル`input/m5test/261103`（実在する大文字小文字を優先）の中身を確認する。GitHubの`main`にはこの新入力群が載っていないため、未知の完全なファイル名や実データ形式は捏造しない。  
> **対象はJitter–剛性比の2系列グラフと、その測定の信頼性評価。追加圧力や時間延長を勝手に行わない。**

## 0. 研究の目的・主図・解釈の範囲

左Body層のヤング率`E_B,L=12 kPa`を固定し、右Body層`E_B,R`を`1.2–11.4 kPa`に変えたとき、**声門流量の周期変動**がどう変化するかを調査する。比較圧力は**900 Paと1400 Pa**、それぞれ同じ18剛性比を試す。各runは同じ初期状態から独立に計算し、圧力上昇後の流量波形からJitterを評価する。

\[
r=E_{B,R}/E_{B,L},\qquad r=0.10,0.15,\ldots,0.95,\qquad E_{B,R}=12r\ {\rm kPa}.
\]

対象18点の対応表：

| ratio `r` | 右Body [kPa] | filename内Bodyトークン（例） |
|---:|---:|---|
| 0.10 | 1.2 | `b1p2c3` |
| 0.15 | 1.8 | `b1p8c3` |
| 0.20 | 2.4 | `b2p4c3` |
| 0.25 | 3.0 | `b3c3` または `b3p0c3`（実在を優先） |
| 0.30 | 3.6 | `b3p6c3` |
| 0.35 | 4.2 | `b4p2c3` |
| 0.40 | 4.8 | `b4p8c3` |
| 0.45 | 5.4 | `b5p4c3` |
| 0.50 | 6.0 | `b6c3` または `b6p0c3` |
| 0.55 | 6.6 | `b6p6c3` |
| 0.60 | 7.2 | `b7p2c3` |
| 0.65 | 7.8 | `b7p8c3` |
| 0.70 | 8.4 | `b8p4c3` |
| 0.75 | 9.0 | `b9c3` または `b9p0c3` |
| 0.80 | 9.6 | `b9p6c3` |
| 0.85 | 10.2 | `b10p2c3` |
| 0.90 | 10.8 | `b10p8c3` |
| 0.95 | 11.4 | `b11p4c3` |

**注意**：表のトークンは解釈規則を示す例。**`b3c3`/`b3p0c3`のどちらが実在するかはローカルで確認**。ファイル名を近似検索で選ばない。

- 主図x軸：`Body Young's modulus ratio E_B,R/E_B,L`。これはBodyヤング率比であり、全声帯の等価剛性比・モード周波数比ではない。
- 主図y軸：**声門流量`Ug`の`Jitter_local [%]`**。**900 Paと1400 Paを別の色/線として重ねる**。各圧力18点の完成・未完成・無発振・周期抽出不合格を区別。
- Cover：左右とも3 kPa（`c3`）で統一。左Body12 kPaの対応freq/modeは既存の正常なファイルを特定し、**全36runで完全同一**とする。
- 自励振動は流体/構造/接触の連成系。**Jitterの変化をただちに右剛性のみの直接因果、乱流、カオス、声質の絶対評価と結びつけない**。
- 対称比`r=1.00`は今回の18点に含まれない。既存12/12 runが**各圧力・時間・固有データ・数値条件の点で互換**な場合に限って「別データの参考点」として使用。今回新規に対称runを追加しない。

## 1. **確定実行条件（以前の0.5秒／900Paだけという仕様を置換）**

| 項目 | 指定 |
|---|---|
| Bodyヤング率 | 左12 kPa、右`12r` kPa（18点） |
| Coverヤング率 | 左右3 kPa（両側のfreq/modeも`c3`） |
| 声門下圧 | **900 Paと1400 Pa、2系列を両方計算** |
| 条件数 | **18 × 2 = 36個の1秒フルrun** |
| 1条件の計算長 | **1.000 s**、`dt=1.0e−5 s`、`nstep=100000` |
| 圧力ランプ | 既存12/6モデルと揃え、当面`pressureRampTimeSec=0.10`（元の`params_used`で再確認） |
| モード数 | 各側20（元run設定から照合し、新VTKに20以上含まれること） |
| 減衰率・接触等 | 元`L12_R6_P900`を基準：`zetaL=zetaR=0.06`、`kc1=5,kc2=6000,kc3=5000`、`ncont=100`、初期間隙0 mm等（**元manifestと照合が先**） |
| 形状・密度・音響条件 | 元runの左右NAS、質量、空気物性、声道・声門下管、流量モデル等を固定 |
| 初期条件 | 各比・各圧力で初期状態から独立計算。前ケースの最終状態を引き継がない |
| **主Jitter窓** | **`[0.50,1.00) s`**（500 ms） |
| **定常性比較窓** | `[0.50,0.75)`および`[0.75,1.00) s` |
| 副参考窓 | `[0.30,0.50) s`（以前の研究との比較用、**1秒主結果には混入しない**） |

**重要な時間窓ルール**：先行の0.5 s runに対する主窓`0.30–0.50s`を、今回の1s run主結果としてそのまま流用してはならない。今回の主結果は`0.50–1.00s`。その間にも包絡や平均周期がドリフトしていれば`NONSTATIONARY`などのQC警告を付ける。1秒計算しただけで定常とは見なさない。

- `Ps`以外の設定を圧力系列間で変えない。`900`と`1400`の違いは肺圧のみ。
- 同じ圧力内の18条件間で変えるのは**右Bodyのfreq/modeペアのみ**。
- `dt,nstep`の変更は**全run共通**（従来との比較では1秒版として明記）。
- 高負荷の`writeDeformationVtu`やコンパイル診断出力を省略する場合は**全36runで同じ出力設定**にし、物理ソルバの状態が変わらないことを既存テスト／短い検証で確認。既存の`writeDeformationVtu=1`は1sの末尾250VTU/runを作り大容量になり得る。今回のJitterでは**`writeDeformationVtu=0`を省容量の選択肢**とし、採用の有無をrunメタに明記。**物理式・接触・流量は変更しない**。

## 2. 作業開始時のローカル確認・データ棚卸し

1. `pwd`、`git status --short`、`git branch --show-current`、`git rev-parse HEAD`、`AGENTS.md`、`df -h`、利用可能メモリ、ビルド`icpx`/OpenMP環境を確認。**未コミット差分と既存runは保存**。`git reset/clean/rebase/checkout`や既存出力の上書き禁止。無断commit/push禁止。
2. `find input -iname '*261103*'`等で**実際のパス・大小文字**を確定。`input/m5test/261103`か`input/M5_test/261103`か等を確認。入力は読み取り専用。
3. `.txt`と`.vtk`および別形式があれば全列挙し、**実際のファイル名・サイズ・SHA256・先頭ヘッダ・周波数列・モード配列数・節点数**を記録。新しいデータをGitHubから取れると仮定せずローカルを探索。
4. **命名パーサ**：`b(?P<int>\d+)(?:p(?P<frac>\d+))?c(?P<cover>\d+)` をファイル名の適切なトークン境界で全一致。`b1p2c3` → `Decimal("1.2")`、`b10p8c3` → `10.8`、`b11p4c3` → `11.4`、`b12c3` → `12`。**`p`を圧力、指数、区切りとは解釈しない**。整数小数部`p0`可。ただし`b1p20c3`を`1.20`とする等桁数はDecimalとして解釈し、**同じ値に複数ファイルが対応したら重複検出**。
5. 新入力のBody18種×Cover3を**それぞれのfreq txt＋mode vtkへ1対1対応**。拡張子を含めると名前の共通prefixは異なるかもしれないので、ファイルの`b...c3`とデータの構造を併用。**同じbトークンを含む別メッシュ/厚さ/δ等のモデルを混ぜない**。
6. 左Body12の元freq+modeおよび左右`M5_surface...nas`の実在・点順・座標・質量規格化を確認。新データに左12がなければ、既存の12kPaモデルを使ってよいが**Cover、幾何、座標、COMSOLの計算法とメッシュ同一性を確認してから**。未知なら`BLOCKED`で中止。
7. 既存`analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900`等から`txt/params_used.txt`、`txt/manifest.txt`、完走状態・ランプ・音響モデル・接触・減衰・ソルバrevisionを回収。**現行`input/param.txt`にある左b4/右b20、nsurfz=18等の値を今回のベースとして勝手に採用しない**。
8. **最初の成果物`mode_inventory.csv`は18行（固有データの種類）。`cases_manifest.csv`は36行（固有データ×圧力）。** 身元が不確かなファイル、欠損、重複、必要20モード不足があればフル計算開始禁止。

**mode_inventory.csvの列**：`ratio,body_left_kpa,body_right_kpa,body_token,cover_kpa,frequency_txt,mode_vtk,raw_format,npoints,nmodes,frequency_1_hz,frequency_20_hz,raw_hashes,prepared_hashes,geometry_check,mode_order_check,status,notes`。

**cases_manifest.csvの列**：`ratio,body_right_kpa,pressure_pa,case_id,param_path,run_dir,left_inputs_hash,right_inputs_hash,mesh_hash,solver_hash,physics_config_hash,nstep,dt_s,analysis_window,status,reason`。

## 3. VTK・周波数の正規化（**ファイル名を直すだけは禁止**）

現行`src/structure/modedata.cpp`の`ModeData::loadFromVTU`は`Displacement_field,_X-component_@_...`等の**XMLテキストDataArray**を解析し、読み込んだXYZモードを内部で`{z,x,y}`へ並べ替える。また`loadFreqDamping`は周波数txtの`Hz omega damping`3列を入力し、TXTは**行順**、VTUは**DataArray名のHz順**を用いる。

1. `.vtk`の先頭が`# vtk DataFile Version`ならlegacy、`<VTKFile...>`ならXML。ASCII/binary/各DataArray名・PointDataかFieldDataか・座標/点数を確認。**`.vtk`を`.vtu`にrenameするだけでは互換にならない**。
2. 現行のロード形式と互換でない場合のみ、原本を一切改変せずに`staging/converted_modes/`へ`XML ASCII VTU`を生成。`vtk`や`meshio`等の既存環境にあるライブラリを利用する。必要な`PointData`のX/Y/Z成分を**周波数名を含めて損失なく保持**し、ローダが実際に読める形式にする。
3. **変換器で勝手にXYZの軸を入れ替えない**。軸変換は現行ローダが行うので、二重変換しない。点の並べ替えも禁止。nPoints、各モードのXYZの値、節点順、モード順、欠損/NaN/Infを厳格検証。
4. 周波数`.txt`が現行3列形式でないときだけ`staging/normalized_freq/`へ変換。`f,omega,damping`に無い物理値を推測して埋めない。原本の`omega`が既知なら`2πf`と差を検査。
5. **最大の品質ゲート**：VTUの周波数DataArray順とTXTの固有周波数行順が**20モードすべて同じモード**を指していることを確かめる。COMSOLの縮退/近接モードがある場合に黙って対応付けしない。頻度名の丸め精度とmode形状を記録し、危険なケースでは`MODE_ORDER_AMBIGUOUS`として中止。
6. 節点順・幾何が左12既存モデル／右モデルで一致しない場合は止める。`ModeData::normalizeModes`は質量正規化をソルバで行うため、変換器がさらに規格化しない。
7. `conversion_manifest.json`にraw/prepared hash、変換前後の最大数値差、PointData/周波数名一覧、対応関係、ソルバでの読込確認を残す。

**代表3ケース（例`r=0.10,0.50,0.95`）の読込・初期化ゲート**を先に通す。`nmode=20`、3方向のモード形状が有効、正常な幾何であることを確認する。**ここで失敗したら36runに進まない**。

## 4. 流量`Ug`由来のJitterの定義と比較方法

### 4.1 流量の正しい時間軸

現行`src/simulation/simulation.cpp`の`dat/airflow_vt.dat`は`n%5==0`の`step Ug`（先頭`# step airflow[m³/s]`）。`dt=1e−5s`なので出力間隔`5e−5s`、流量データは20 kHz相当。**第1列は秒でなくstep**。`t=step*dt`、`Ug`はm³/s。読み込み後に`nstep`、step連続性（5刻み）、時刻範囲、有限性、0に近い流量、ヘッダ、サンプル欠落を検査する。`nwrite`をサンプル間隔と混同しない。

Jitter主信号は**声門流量`Ug`だけ**。声門面積や左/右表面変位で計算した値は、検出の妥当性をみる**独立な補助指標**として区別する。現行`tools/jitter.py`は声門面積の5点移動平均との差であり**今回の流量周期Jitterではない**。`tools/compare_stiffness_metrics.py`は左右変位のJitterなので、`Ug`とラベルを入れ替えて使わない。

### 4.2 流量パルス時刻`t_i`の取得

- パルス検出候補を**至少2通り**持つ。**A）流量`Ug`の局所最大値（`scipy.signal.find_peaks`＋prominence/最小距離）**、**B）同じ流量の下降相の急峻な変化、最小値または包絡に基づく頑健な周期境界**。Bは「閉鎖イベント」と断言せず、選択根拠を図示。2方法で取り出したイベントは同じ物理周期かどうかを波形と自己相関/スペクトルで確認する。
- Aはパルスが1周期に複数峰を持つ時に過大なピーク数になる。Bも部分閉鎖や連続開口時には存在しないことがある。どちらも**検出不能なら無理に周期を作らない**。
- 最初に`Ug`の生波形とピーク位置を代表条件で可視化し、peak prominenceの固定基本値（例：各主窓の`P95−P5`の10%）、最小距離、基礎周波数帯域を**事前に固定**。追加の感度検証は`5%,10%,15%`で行い、条件ごとに都合のよい閾値を選んで主結果にしない。
- 流量が基礎周波数100〜300 Hz程度なら20kHz相当のサンプルで各周期が約67〜200点。整数sample単位の周期差が擬似Jitterを生むため、**サブサンプル補間**（局所3点放物線等）ON/OFFで両方算出。補間が局所形状に適さなければ除外フラグ。
- `t_i`が厳密に増加し、欠落ピーク・重複・明らかな偽ピークがないか、`Δt`・prominence・包絡・自己相関・FFTから検査。**単一のFFT主ピークは周期列検出の代わりにならない**。
- 主窓`[0.50,1.00)`の境界では**窓を少し越える波形を読んでから境界内に有効な全周期を切り出す**。窓端をまたぐ中途半端な周期を誤って含めない。周期に使用したピークの数と区間を記録。

### 4.3 周期列と**複数Jitter指標**（すべて`Ug`から）

`M`個の有効な流量パルス時刻から`N=M−1`個の周期`T_i=t_{i+1}-t_i`を得る。**主図は一般的な`Jitter_local`で統一**し、他の指標を補助比較として別図/表へ。

**(1) Jitter local（主指標）：**
\[
J_{\rm local}[\%]
=100\,\frac{\frac{1}{N-1}\sum_{i=1}^{N-1}|T_{i+1}-T_i|}{\bar T}.
\]

**(2) RAP（3周期平均との差）：**
\[
J_{\rm RAP}[\%]
=100\,\frac{\frac{1}{N-2}\sum_{i=2}^{N-1}\left|T_i-\frac{T_{i-1}+T_i+T_{i+1}}{3}\right|}{\bar T}.
\]

**(3) PPQ5（5周期平均との差）：**
\[
J_{\rm PPQ5}[\%]
=100\,\frac{\frac{1}{N-4}\sum_{i=3}^{N-2}\left|T_i-\frac{T_{i-2}+T_{i-1}+T_i+T_{i+1}+T_{i+2}}{5}\right|}{\bar T}.
\]

**(4) CV of periods（参考）：**
\[
CV_T[\%]=100\,\frac{\mathrm{std}(T_i)}{\bar T}.
\]

**(5) 整数sampleのJitter local／補間後のJitter local**を両方記録し、離散サンプリングによる見かけの不規則性を定量化。

**(6) 流量の別パルス検出法BからのJitter local**：同じ期間・周期定義で検出法差を計算し、A/Bの有効周期数・イベント対応と一緒に示す。Bのイベントが不明瞭なら`METHOD_B_NOT_APPLICABLE`であり、検出法Aに強制的に揃えない。

**(7) 任意の独立波形による検算**：声門最小面積・音響圧波形（保存され実際に意味のある信号に限る）から周期を検出し、Ugの代表的周期長や不規則なイベントの時刻と照合してよい。ただし**違う物理量から得たJitterを同じ数値として平均しない**。表には`signal`を必ず含める。

**注意**：RAPとPPQ5は期間ゆらぎの異なる平滑窓による指標であり、互いに独立な物理測定ではない。変動の少ない純周期合成波、交互周期、緩やかな周波数ドリフト、ダブルピーク、失われた周期、振幅変調、強い雑音、無発振の合成テストで正しく区別できること。過去の「5点移動平均」とは**対象が周期列か信号振幅かを厳密に区別**。

### 4.4 周期性と品質判定

**Jitter値の算出可能性**と**不規則な振動状態**は同じではない。高Jitterをすぐに「カオス」や「周期破綻」と呼ばない。以下の状態を全条件に付与：

- `VALID_STATIONARY`：一定の手順でパルスが抽出され、周期数と半窓安定性を満たす。
- `VALID_WINDOW_DEPENDENT`：周期抽出可能だが前後半でJitterまたは平均周期/包絡が大きく変わる。
- `PERIOD_MULTIPLIED`：規則的な交互振動・サブハーモニクス等を検出。Jitterは定義とイベント選択に依存、注記して掲載。
- `AMBIGUOUS_PEAKS`：多峰性・検出法A/Bの食い違い・prominence感度などにより同一種類の周期境界が不明。
- `NO_OSCILLATION`、`INSUFFICIENT_CYCLES`、`RUN_FAILED`、`INPUT_BLOCKED`、`NUMERIC_INVALID`等。

主指標は暫定`N≥10`個の**連続有効周期**と、十分な流量振幅/SNRを必要とする。少なければJitterは`NaN + reason`で記録。**無発振を0%にしない**。1:2の規則的周期倍化は、交互のパルス間隔を使ったlocal Jitterが高くても**それだけで「ランダムゆらぎ」ではない**。周期倍化の検出は自己相関・FFTの半周波数や2周期繰返しで慎重に判断する。

**定常性**：主窓`[0.50,1.00)`を`[0.50,0.75)`と`[0.75,1.00)`に分け、平均周期、`J_local/RAP/PPQ5`、流量RMS/ピーク幅、包絡のドリフトを比較。閾値は事前に`analysis_config.json`で固定し、実際の判定率とflagを報告。主窓が非定常なら主図に斜線／マーカー変更などで区別し、勝手に0.75–1.0に窓を切り替えて本値を差し替えない。

**感度検査**：Aのprominence 5/10/15%、peak時刻補間ON/OFF、ピーク候補A/B、主窓と半窓、平滑化なし／必要なときの対照を比較し、相対変化・周期数・ランキングがどの程度安定か記録。既存1sデータで終える。**解析結果が都合のよいような条件別手法変更は不可**。

## 5. 36runの生成・実行管理

**本体を改造しない。** 既存`src/simulation/readParams.cpp`は`leftFrequencyFile/rightFrequencyFile/leftModeFile/rightModeFile/leftSurfaceNasFile/rightSurfaceNasFile`を受け付ける。過去の`tools/run_right_stiffness_sweep.py`は**整数Body値、従来`M5_test`と`.vtu`前提**で今回にはそのまま使えない。旧ツールは破壊せず新しい実装で`p`表記とVTK変換後の入力を扱う。

推奨新設ディレクトリ：`tools/body_ratio_jitter_261103/`。
- `inventory_modes.py`：小数BodyトークンをDecimalでparse、右18対のtxt/vtk照合・hash・nModes・点順・freq順・Cover/幾何チェック。
- `convert_modes.py`：必要な場合だけ`.vtk`→`.vtu`準備と`.txt`の正規化、変換前後のモード値・周波数の検算。
- `prepare_sweep.py`：完成済み12/6条件のparamを正とし、`Ps=900/1400`、`nstep=100000`、右Body入力だけを変更した**36独立`param.txt`**を新規生成。`--dry-run`を標準動作にする。相対パスは生成param移動により壊れるので必要なfreq/mode/NASを絶対パス化。`input/param.txt`原本に書き込まない。
- `run_sweep.py`：`--execute`をユーザー指定でのみ受け取り実行、`--resume`、`--threads`、`--fail-fast`、`--pressures 900 1400`。条件別runディレクトリの`SIMULATION_RUN_DIR`と保存先独立性を実際に確認してから実行。
- `analyze_flow_jitter.py`：`Ug`のA/Bイベント検出、subsample、`local/RAP/PPQ5/CV`、ピークQC、主窓・半窓・感度、図・生周期CSV。
- `plot_jitter_sweep.py`：全36条件集計、900と1400の同一x軸への表示、QC理由の表示、手法別指標の比較図、代表`Ug`波形、圧力差`ΔJ=J_1400−J_900`の比較（**共通有効剛性比だけ**）。

**実行前チェック**：完成・再利用条件は`ratio,Ps,nstep,dt,full physics config hash,left/right raw input hashes,converted hashes,solver hash,analysis method version,expected output`が一致する場合のみ。既存`0.5s`のrunや別COMSOL入力の`r=0.5`を**1s新系列として流用しない**。

**ディレクトリ例**：
~~~text
analysis_runs/jitter_body_ratio_261103_900_1400_1s/
  mode_inventory.csv
  conversion_manifest.json
  analysis_config.json
  cases_manifest.csv
  cases/
    r0p10_P900/
      input/param.txt
      dat/airflow_vt.dat
      txt/params_used.txt
      txt/manifest.txt
      jitter_metrics.json
      flow_peaks.csv
      flow_cycles.csv
      status.json
      stdout.log
      stderr.log
    r0p10_P1400/
      ...
  jitter_summary.csv
  jitter_method_comparison.csv
  jitter_sensitivity.csv
  figures/
    jitter_vs_ratio_900_1400.png
    jitter_vs_ratio_900_1400.svg
    jitter_method_comparison.png
    jitter_qc.png
    flow_peak_examples.png
    pressure_difference.png
  preflight_report.md
  report.md
~~~

`status.json`は`not_started/running/completed/analysis_failed/input_blocked/run_failed/jitter_invalid`等の状態とコード、理由、入力hash、計算時間を保持。**ソルバ完走`completed`と有効なJitter`VALID_STATIONARY`は別の判定**。進捗が中断しても**既存runを消さず再開**。生VTKや大容量VTU、波形の全コピーをGitへ加えない。

## 6. メイングラフと結果の解釈

### 図1：必須の主図

`figures/jitter_vs_ratio_900_1400.png`/`.svg`：
- x`=E_B,R/E_B,L`＝`0.10..0.95`、数値軸、tickは意味のある位置`0.1,0.2,...,1.0`等。
- y`=Flow Jitter_local [%]`。**900 Paと1400 Paを2系列**で同じ図へ（色だけでなくマーカー/線種も変える）。値のない点は0として線で結ばない。曲線は**隣接するQC合格点に限り**連結し、欠測で区間を橋渡ししない。
- `VALID_WINDOW_DEPENDENT`や`PERIOD_MULTIPLIED`は通常と違うマーカーにして、定常Jitterと混同しない。`AMBIGUOUS_PEAKS/NO_OSCILLATION`等は別のQC帯・テーブルに配置し、意味のないy=0表示を避ける。
- **主図の主系列はUg–local(A,補間ON/主設定)で固定**。複数のJitter定義を同じ主図の900/1400各曲線に混ぜない。解析手法を変えた比較は別図。
- `r=1.0`参考点は各Psで同条件1sの既存runがある場合のみ注記付きで追加（原則36点の外）。

### 図2：周期ゆらぎ指標の比較

`jitter_method_comparison.png`：
- 900/1400ごとに`J_local`、`RAP`、`PPQ5`を同じ`Ug`の**同じ有効パルス列**で算出して比較。**値の尺度が指標間で異なるので絶対的同値は期待しない**。剛性比に対する傾向、ピーク位置、順位・相関、検出法差、`N cycles`を評価。
- 補助で`CV_T`、subsample ON/OFF、流量検出法A/B、参考の面積/音圧からの周期検証を提示。**他信号の値を流量Jitterの平均に混入しない**。

### 図3～5：品質・代表波形・圧力差

- `jitter_qc.png`：両圧力、valid cycle count、各種flag、窓安定性、prominence感度、整数sample vs補間で差が大きい箇所。
- `flow_peak_examples.png`：900/1400の低/高Jitter、周期倍化/多峰/不合格条件があれば代表波形（mL/sへの表示換算：`1 m³/s=10^6 mL/s`）、ピーク時刻・`T_i`が追跡できるように。
- `pressure_difference.png`：同一`r`で有効な`J_1400−J_900`。**一方が無効な点は空欄**。平均差・分布・どこで圧力依存が強いかを分析。両曲線の交差をデータのない区間へ外挿しない。

### レポートの主題

- 剛性比によってJitterが単調増減するか、局所的な増大域があるか。**急変点の前後は波形/ピーク抽出が妥当か**。
- 900 Paと1400 PaでJitterの変化が違うか。同じ比で圧力上昇が周期ゆらぎを増すか減らすか。**圧力を変えたことで発振状態自体が変わる可能性**。
- 複数方法で傾向が共通しているか。`local`だけ異常に大きいならダブルピーク/周期倍化/イベント抽出の影響を検討。
- これまでの12/6・900Paの**振幅左右比0.302**と、今回の流量Jitterは異なる指標。局所振幅非対称≠Jitter大とする。音響学的な声質劣化の定量的保証とも言わない。
- 36ケースの数値が得られるまで**どの剛性比が最大Jitterになるかを予言したり、架空のスムーズな曲線を描いたりしない**。

## 7. 実行順序（無駄な36本を避ける）

**Phase 0：実行なし・Preflight**
1. 作業ツリー/既存run保全、入力36条件のリストと新固有データ18種の厳格検査。
2. モード変換の必要性判定、20モード・周波数順・幾何の照合。
3. 参照`L12_R6_P900`から物理設定を確定し、`Ps`と`nstep`以外の差を右freq/modeに限定。全36個のparamを生成し、差分を表にする。
4. Jitter解析器の合成波形テストを実行。`preflight_report.md`で全ゲートを`PASS`としてから本番へ。

**Phase 1：代表条件で波形QC**
1. `r=0.50`の900 Paと1400 Paを代表として**各1秒**のrun、`Ug`ピーク/周期/主窓/前後半比較を可視化。**この2本は36本の内数**。
2. もしpressure 1400 Paで別物理式・モード出力の不適合・数値異常・発散が生じたら、残りを自動実行せず原因と入力差分を報告。
3. `r=0.10`と`0.95`の入力の**読み込み動作確認**。初期化だけでよい場合はフルrunに数えない。確認済み条件以外のフルrunを無条件に走らせない。
4. 主ピーク手順を固定し`analysis_config.json`に保存。主指標は`Ug` Local jitter、その他手法の感度比較は副次。

**Phase 2：残り34本**
- Phase1がPASSなら全36点を**同じ物理設定・同じ方法**で順次実行。大量並列の競合を避ける。`OMP_NUM_THREADS`／`OMP_PROC_BIND`／`OMP_PLACES`を固定、CPU/メモリ/ディスク上限を確認。途中resume可。
- `SIMULATION_RUN_DIR`で出力をケースごとに隔離し、`output/`、`result/`、`input/`、過去`analysis_runs/`を改変しない。`stdout/stderr`、`params_used.txt`、`manifest`、実行結果`status`を残す。
- 完走したrunから順次`Ug`を解析し、明白なピーク誤認や数値異常が出たら**追加runを止める**。見かけの`Jitter=0`で穴を埋めない。

**Phase 3：図と論文向けの評価**
- 36条件の`jitter_summary.csv`、local/RAP/PPQ5/CV/A-B差の`jitter_method_comparison.csv`、感度、周期表を作る。
- 900/1400の2曲線主図、Jitter手法比較、代表Ug波形、QC、圧力差を図にまとめる。
- `report.md`に、定義・計算条件・有効ケース数、各系列の最大/最小、ピーク検出の妥当性、比率と圧力による変化、未確定な機構を記録する。
- 急変域などがあれば追加精細化は**提案のみ**。今回の36本以外を勝手に計算しない。

## 8. テストと品質ゲート（必須）

**ファイルパーサ**
- `b1p2c3→1.2`、`b11p4c3→11.4`、`b10p8c3→10.8`、`b12c3→12`。不正`p`、異なる`c`、重複`b3c3/b3p0c3`、別モデル名の衝突を検出。
- 18種のfreqとmodeを一意対応。計36ケース。**同名の頻度txtを複数モデルで使い回さない**。両圧力で右モードファイルhashが同一、左12が全条件で同一。

**VTK/固有モード**
- legacy/XML判定、20以上の周波数とXYZモード、座標・節点順・幾何、元/変換後の一致。
- DataArray名の`Hz`と周波数txt20行のmode対応。近接固有値・重複キーは不明なら`MODE_ORDER_AMBIGUOUS`で止める。

**設定**
- `dt=1e−5`、`nstep=100000`、各独立run`1s`、`Ps∈{900,1400}`、`r`は18種のexact decimal。
- 圧力系列間の差は`Ps`のみ。剛性比条件間の差は**右周波数+右モード**のみ。既存12/6の参照設定（Cover3/3、ランプ、減衰、接触、管路等）と照合。出力専用オプションを切った場合は別途記録。
- 同じ条件とソルバhashならresume可、他条件を「完了」と誤認しない。

**UgとJitter**
- 1列目はstep、`n%5==0`を検証。m³/s→mL/sの換算と時間軸。
- 合成一定周期`J_local/RAP/PPQ5≈0`、period交互例、単調drift、多峰波形、ピーク欠落、無振動、スケール違い、NaN、短い窓。
- RAP/PPQ5の端点計算、必要周期数、窓左閉右開、subsample補間の有効性、全36点同一prominence/判定基準。
- 検出法A/Bは別の周期候補を**無理に合流しない**。A/Bが異なるなら波形・周期数を見せ、`AMBIGUOUS_PEAKS`を出す。
- メイングラフは**2圧力系列、流量Jitter localのみ**。全18比率について成功・失敗が表に残る。

## 9. 想定コマンド（ローカルの実ファイルを確認して最終確定する）

~~~bash
# まず新COMSOL入力の完全な実在確認
python tools/body_ratio_jitter_261103/inventory_modes.py \
  --source-dir input/m5test/261103 \
  --left-body-kpa 12 --ratios 0.10:0.95:0.05 \
  --output-dir analysis_runs/jitter_body_ratio_261103_900_1400_1s/inventory

# 36条件を準備（dry-runでソルバは起動させない）
python tools/body_ratio_jitter_261103/prepare_sweep.py \
  --inventory analysis_runs/jitter_body_ratio_261103_900_1400_1s/inventory/mode_inventory.csv \
  --base-run analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900 \
  --pressures 900 1400 --duration 1.0 --dry-run \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900_1400_1s

# 読込と合成流量Jitterテストが通ってから、代表2条件→残り34条件
python tools/body_ratio_jitter_261103/run_sweep.py \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900_1400_1s \
  --pressures 900 1400 --execute --resume --threads 4

# 流量Jitter方法比較、QC、900/1400の主図
python tools/body_ratio_jitter_261103/plot_jitter_sweep.py \
  --sweep-dir analysis_runs/jitter_body_ratio_261103_900_1400_1s
~~~

これらは**完成後のCLI仕様例**。フォルダの実際の大小文字を調査し、依存関係・実行バイナリを確認してCodexが具体的コマンドを確定する。**各runを実行するかどうかは入力ゲートで制御**し、dry-runで書かれた未検証パラメータを実計算に回さない。

## 10. Codexへの最終依頼・完了条件

**実際の`simulation-main/input/m5test/261103`にある新しい固有データ（`b1p2c3`が1.2 kPaである`p`小数表記）を取り込み、左Body12/右1.2–11.4kPaの18点について、`Ps=900Pa,1400Pa`で各`1秒`、計`36条件`の自己振動計算を行ってください。前に完成した12/6モデルの物理条件を基準とし、COMSOLモードと周波数を厳密に対応づけてください。主結果は`Ug`ピーク周期からの`Jitter_local [%]`を横軸`E_B,R/E_B,L`で2圧力系列に描き、`Ug`由来のRAP・PPQ5・周期CV・別イベント検出法・前後半比較で信頼性を検証してください。欠損入力やQC不合格を補わず、既存runや元入力を上書きせず、全条件の数値表・品質表・図・研究的レポートを完成させてください。**

**完了チェック**：
- [ ] `b1p2c3`など小数命名の実ファイル18ペアに正しいBody値を割当、Cover=3、20モード、NAS/幾何・Hz順・VK/VTU互換を検証。
- [ ] 36条件は`18 ratio × 900/1400Pa`、各`dt=1e−5`・`nstep=100000`・1.00 s、ランプ/接触/減衰/流体系は同一設定。
- [ ] 実行前に参照12/6の正しい条件を確認、param原本・新VTK原本・旧runは非改変。
- [ ] `Ug`から有効パルスを検出してlocal/RAP/PPQ5/CV、subsample ON/OFF、検出法A/B、定常性・ピークQCを全条件で統一。
- [ ] **主窓`[0.50,1.00)`、半窓`[0.50,0.75)`・`[0.75,1.00)`を厳守**。旧0.5sの主窓を流用していない。
- [ ] x=Body剛性比、y=流量Jitter_local[%]の主図に**900/1400Paの2系列**。失敗点・無発振は0%で埋めない。
- [ ] `mode_inventory.csv`18行、`cases_manifest.csv`36行、`jitter_summary.csv`36行、`jitter_method_comparison.csv`、`flow_cycles.csv`、図・`preflight_report.md`・`report.md`を保存。
- [ ] 追加物理条件・圧力・計算時間の無断拡張なし、既存入力/出力の上書きなし、C++物理式変更なし、巨大runのGitコミットなし。
