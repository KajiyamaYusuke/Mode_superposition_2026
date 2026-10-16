# Phase 1 機能・入出力棚卸し

調査日: 2026-10-14（JST）。指示書: `CODEX_cleanup_phase1_2.md`。本書のソース行番号は、特記しない限り**整理前の作業ツリー**に対応する。整理前コピーは `/tmp/simulation-cleanup-1cGHFz/work`、開始時の差分・SHA256・既存出力のハッシュは同ディレクトリの親にある `start.json` に保存した。Phase 2 の結果は `cleanup_change_report.md` に記録する。

## 1. 作業基準と保護

- HEAD: `3ecd141786d2ee3384276853f48a84a7c14ad38b`。HEADではなく、ユーザーが指定した未コミット変更込みの作業ツリーが比較基準。
- 既存差分: `input/param.txt` の nstep=50000、zetaL=zetaR=0.02、および `src/simulation.cpp:1624` の有効なVTU呼び出し。維持必須。
- 停止ゲート第2項の例外、および既存CTest失敗1件を維持したまま続行することは、ユーザーの後続指示で明示的に許可された。
- ネストした `AGENTS.md` は検索で検出されず、ルートの `AGENTS.md` を適用。CMake・指示書に編集上の競合なし。ただしAGENTSの「自動テストなし」は現物と異なり、CMakeには3件のCTestが存在する。
- `input/` は55ファイル、約7.6GiB。モードVTU、NAS、周波数txt、param、Zone.Identifierが存在する。削除・移動・変更しない。
- 既存 `output/` 約1.9GiB、`result/` 約5.5GiBの計351ファイルをSHA256・サイズ・更新日時で記録した。元リポジトリではビルド・計算を実行しない。
- ソース・ヘッダー・CMake・tests・tools・必要入力を `/tmp` にコピー。実行cwdもコピー内の `work/build` とし、`geometry.cpp` の `../output` 直接出力を隔離する。`SIMULATION_RUN_DIR` はテスト環境から除外。
- 新規リポジトリ成果物は本書と `cleanup_change_report.md` のみ。既存同名ファイルは開始時に存在しなかった。ディレクトリ再構成、移動、リネーム、commit、push、依存追加は行わない。

## 2. 入力・ビルド・動作モード

`main.cpp:5` → `Simulation::initialize` → パラメータ検証 → 左右のGeometry/ModeData → 右側反転・初期間隙設定 → State/ForceCalculator初期化 → `Simulation::run`。

時間ループは、流路断面 → 流体または正弦波の表面荷重 → 流体モード力保存・時系列出力 → 接触反復（モード射影、Newmark、表面予測変位、接触、収束判定）→ 最終荷重で再射影・Newmark → 全節点変位 → 診断専用最終接触評価 → 状態確定 → VTU・音響サンプル → 終了後WAV。順序を変更しない。

| 分岐 | 現物の条件・参照 | 現在値・検証方針 |
|---|---|---|
| 通常加振 | `forceCalculator.cpp:267` iforce=0 | paramの実効スロット=0。流路・音響回路・圧力・表面力を計算 |
| 正弦波加振 | `forceCalculator.cpp:244` iforce=1 | forcef=100Hz、famp=-1e-5。forceDirection=0は左右ともx方向、1は左右逆向きy方向。試験は実効iforceだけ1に変更 |
| 旧iforceスロット | `readParams.cpp:89` legacyIforce、同:96 実効iforce | 両スロットを読み込む。互換処理を維持 |
| 接触反復 | `simulation.cpp:1230` icont<=ncont | ncont=100、残差で早期終了。ncont=0でも最終Newmarkと診断専用接触評価が存在 |
| 接触収束 | `simulation.cpp:1345` | contactFlag、力残差<1e-4、めり込み残差<1e-7。変更禁止 |
| 声道 | `flowModel.cpp:33` L_vt>1e-6 かつ断面>1e-12 | `SimulationParams.h:95` L_vt=0。今回の実行は声道なし。声道ありの経路は維持するが実行未検証 |
| VTU | `simulation.cpp:1623` n%20==0 && nstep-n<=5000 | **呼び出し有効**。指示書の「コメントアウト」は現在の実体には当てはまらない。初期専用呼び出し:992はコメントのまま |
| プロファイル | `CMakeLists.txt:8,49` PROFILE_CALCDIS | 既定OFF。既存の計測・コンパイル条件を残す |
| 擾乱・荷重凍結 | `simulation.cpp:950,998`、`readParams.cpp:314` | perturbationEnabled=false、diagnosticLoadMode=live。freeze_contact/fluid/allの実処理は維持 |
| AP端除外 | `forceCalculator.cpp:733` contactIsExcluded | apContactDiagnosticEnabled=false。有効時は荷重にも影響するため診断ログ扱いで削除しない |

現行入力は `simulation.cpp:44–68` がパラメータファイル基準で解決する。指示書の「モード入力パスがハードコード」は現物では該当しない。左右のVTU/NAS/周波数の設定は独立だが、本番paramでは同一ファイルを指定している。

| 入力 | 実在パス | 内容・参照 |
|---|---|---|
| パラメータ | `input/param.txt` | nmode=100、nsurfz旧スロット=18、nstep=50000、nwrite=20、dt=1e-5、zeta左右0.02、kc=(5,6000,5000)、ncont=100、Ps=2900Pa、mass=1.07155257e-3kg |
| 左右モード・幾何 | `input/M5_test/M5_mode_kawahara_mesh7_soft.vtu` | 約287MiB。片側45333点、248268セル。Geometry/ModeDataがそれぞれ読み込む |
| 左右表面 | `input/M5_test/M5_surface_kawahara_mesh7_soft.nas` | NAS→VTU点対応、接続性とzレベルから表面格子を構築 |
| 左右周波数 | `input/M5_test/M5_freq_kawahara_mesh7_soft.txt` | 選択モード番号に対応する周波数等を読み込む |
| 他の研究入力 | `input/mode.vtu`、`input/M5_test/M5_mode_T*.vtu`、`M5_surface_T*.nas`、`M5_freq*.txt`等 | 現行paramが使わないものも研究入力として全保持。今回実行は上記3モデルファイルのみ |

Release: IntelLLVM/icpx 2025.1.0、C++17、`-O3 -std=gnu++17 -fiopenmp`、OpenMP必須、CMake 3.26.4。`-ffast-math`なし。`BUILD_TESTING`既定ON。Python依存は `pyproject.toml` と `requirements-analysis.txt` に存在し、今回追加しない。

その他の設定・文書は、Python 3.10を指定する `.python-version`、依存固定の `uv.lock`、出力等を除外する `.gitignore`、圧力掃引の `tools/pressure_sweep_config.json`、解析説明の `tools/README_dynamics_analysis.md`、既存研究の `MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md`。全て保持する。AGENTSが挙げる `nonperiodic_vibration_analysis_plan.md` は現物のファイル検索では見つからなかった。`analysis_runs/` と `diagnostic_runs/` に既存の研究設定・成果物があり、参照調査のみで変更しない。

## 3. 機能・削除候補分類

分類は指示書の5種類のみ。`REMOVE_SAFE`は許可リスト内で、リポジトリのソース・ヘッダー・tools・テスト・文書の参照を調べた範囲で読み手がなく、出力以外の副作用がないもの。リポジトリ外の利用は調査対象外。実施は両モードの基準実行・維持出力比較が成立する場合に限る。

| ID | パス:関数/変数/出力ファイル | 機能と呼出し条件 | 依存元/参照先 | 判定 | 理由（参照箇所） | 計算影響リスク | Phase 2 実施可否 |
|---|---|---|---|---|---|---|---|
| C01 | `src/main.cpp:5` main | 初期化・実行・例外の終了コード | Simulation | KEEP_CORE | ERROR表示とreturn1を保持 | 高 | 不可 |
| C02 | `src/readParams.cpp:50,221` loadFromFile/validate、`include/SimulationParams.h` | 旧/実効iforce、位置・名前付き設定、検証 | 全構成要素、input_file_params | KEEP_CORE | 互換と安全判定。既定値不整合も今回修正禁止 | 高 | 不可 |
| C03 | `src/geometry.cpp:323,584,698` loadFromVTK/surfArea/surfExtractFromNAS、`include/Geometry.h` | 幾何・面積・節点対応 | 初期化→流路・接触・射影・VTU | KEEP_CORE | NASの探索・格子検査を含む | 高 | 不可 |
| C04 | `src/modedata.cpp:10,31,188,229`、`include/ModeData.h` | モード選択・VTU読込・周波数・質量正規化 | initialize、State、ModalProjector | KEEP_CORE | 座標置換と正規化も物理計算 | 高 | 不可 |
| C05 | `src/channelSections.cpp:11,54,106,171`、`include/ChannelSections.h` | 共通x平面切断、正間隙と三乗積分 | updateChannelSections→FlowModel/壁圧 | KEEP_CORE | 表面と断面の補間・閾値を維持 | 高 | 不可 |
| C06 | `src/flowModel.cpp:8,42,65`、`include/FlowModel.h` | 音響回路、Newton流量更新、声道分岐 | applyFluidLoads | KEEP_CORE | 圧力ランプ含め維持 | 高 | 不可 |
| C07 | `src/separationSelection.cpp:6`、`include/SeparationSelection.h` | 下流面積比による分離、閉鎖・fallback | applyFluidLoads、CTest | KEEP_CORE | 分離条件を維持 | 高 | 不可 |
| C08 | `src/forceCalculator.cpp:220` applyFluidLoads iforce=0 | 圧力回復クランプ・表面への荷重写像 | run、FlowModel、SurfaceLoad | KEEP_CORE | :351のクランプを含む | 高 | 不可 |
| C09 | `src/forceCalculator.cpp:244` iforce=1 | 正弦波・左右配分・方向 | sp.forcef/famp/forceDirection | KEEP_TEST | 試験モードを全保持 | 高 | 不可 |
| C10 | `src/forceCalculator.cpp:620` applyContactLoads | 候補・cache・fallback、法線・接触力、AP診断と除外 | run→残差・表面力・モード力 | KEEP_CORE | :1168以降の検出・配分は編集しない | 高 | 不可（専用ログ行のみ別項） |
| C11 | `include/ForceCalculator.h:115–124` 残差/接触力、`forceCalculator.cpp:1620–1675,1845` | 更新差、正規化力残差、めり込み差 | simulation.cpp:1345収束 | KEEP_CORE | 名前が診断風でも分岐に直結 | 高 | 不可 |
| C12 | `src/modalProjector.cpp:3`、`include/ModalProjector.h`、`src/surfaceLoad.cpp:4`、`include/SurfaceLoad.h` | 左右の表面力保存・モード射影 | projectLoadsToModes、接触加算 | KEEP_CORE | モード内加算順・OpenMP維持 | 高 | 不可 |
| C13 | `src/rungeStep.cpp:25` newmarkStep、`include/TimeIntegrator.h` | β=0.275625、γ=0.55、予測修正 | runの反復・最終再解法 | KEEP_CORE | 積分係数・順序維持 | 高 | 不可 |
| C14 | `src/state.cpp:5,51,91,133`、`include/State.h` | 初期化、全点/表面モード重畳、uf2u | initialize/run/VTU | KEEP_CORE | q/qdot/qddot、vel、surfacePointIds等維持 | 高 | 不可 |
| C15 | `include/Displacement.h`、`include/FluidConstants.h` | 変位構造・閉鎖閾値 | State/ModeData/flow/separation | KEEP_CORE | 単位換算・共通閾値 | 高 | 不可 |
| C16 | `src/simulation.cpp:1663` writeVTKCombined、`include/Simulation.h:33` | 左右Points、connectivity、offsets、types | :1624有効呼び出し | KEEP_TEST | 既存変更を含め維持必須 | 高 | 不可 |
| C17 | `include/wavwrite.h:36` WavWriter::save | 終了後の16bit mono PCM | simulation.cpp:1636 | KEEP_OUTPUT | 音響信号の正規化・ヘッダー保持 | 中 | 不可 |
| C18 | `src/simulation.cpp:310–329,1638–1657` timing/反復集計 | 各計算と終了時cout | run | KEEP_OUTPUT | 指示書で明示維持。時間は数値一致対象外 | 低 | 不可 |
| D01 | `src/simulation.cpp:420` fstepdbg/debug_step_summary.csv | 初期化時ヘッダーのみ | `tools/kcExplore.py:12,118` | REVIEW_ONLY | ヘッダーのみだが解析読み手が実在。許可条件「参照なし」を満たさない | 入出力互換 | 保留 |
| D02 | `src/simulation.cpp:431,1252` fmodedbg/debug_mode_summary.csv | 接触反復の n%10==0 またはt>0.12 | ログ専用maxAbsIndex | REMOVE_SAFE | 全参照は生成・ヘッダー・after_f2mode行。読み手なし、fiL/Rを読むだけ | 低、回帰比較必須 | 可 |
| D03 | `src/simulation.cpp:890` maxAbsIndex、maxFiL/R・imaxFiL/R | D02の引数計算 | :1253–1260のみ | REMOVE_SAFE | D02除去で全参照消滅。fi・状態を変更しない | 低 | 可 |
| D04 | `src/simulation.cpp:853,869,909` maxAbsNodeDisp/minMaxIndex/maxAbsPsurf | 宣言のみ、呼び出しなし | 全文・識別子参照確認 | REMOVE_SAFE | 許可リストに明記された未使用ローカルラムダ | 低（実行処理なし） | 可 |
| D05 | `src/simulation.cpp:814` maxAbsVector | fiL/R最大値 | :1336 contact_coupling_debug.csv | KEEP_OUTPUT | 未使用ではない。維持出力に依存 | 中 | 不可 |
| D06 | `src/forceCalculator.cpp:112,324` debugForceFile/debug_force.txt、`include/ForceCalculator.h:131` | 初期化・流体の毎step minA/currentPgログ | 実行系への読み戻しなし | REMOVE_SAFE | open/header/行のみ。setw等は専用streamへの操作 | 低 | 可 |
| D07 | `src/forceCalculator.cpp:184,367,559` flowDebugFile/debug_flow_detail.csv | 流体のn%10==0またはt>0.12、異常圧力追記 | `tools/analyze_self_oscillation.py:191` | REVIEW_ONLY | legacy fallbackとして実際に読む。異常値検出とも混在 | 互換・異常解析 | 保留 |
| D08 | `src/forceCalculator.cpp:117,672,1356,1583` contactDebugFile/contact_debug.csv | 0.120–0.135sの区間統計 | 専用debug変数・合計値 | REVIEW_ONLY | 0.02s試験では区間を通過せず、安全性を実証できない削除は指示書どおり保留 | 未検証経路 | 保留 |
| D09 | `src/forceCalculator.cpp:129,1677` contactMonitorFile/contact_monitor.csv、`include/ForceCalculator.h:133` | 非diagnosticOnlyの接触呼び出しごと | monitorI/J、gap、fdis、残差を読み出すのみ | REMOVE_SAFE | CSVの読み手なし。open/header/専用出力ブロックだけを対象とする | 低、接触比較必須 | 可 |
| D10 | `src/forceCalculator.cpp:139,1700` contactIterationFile/contact_iteration_debug.csv | 接触反復ごと | `tools/run_pressure_sweep.py:32`、README:41 | REVIEW_ONLY | スイープ保存対象、研究上の参照あり | 入出力互換 | 保留 |
| D11 | `src/forceCalculator.cpp:153,1717` contactSearchFile/contact_search_debug.csv、`include/ForceCalculator.h:135` | 非diagnosticOnly接触呼び出しごと | 探索統計を読むだけ | REMOVE_SAFE | CSVの読み手なし。open/header/出力ブロックだけを対象 | 低 | 可 |
| D12 | `src/forceCalculator.cpp:693–696,1158,1459–1492` 候補カウンタ | 探索ラムダの参照引数にも使用 | runCandidate/cache/two-pointer/fallback | REVIEW_ONLY | D11削除後も探索の構造と参照引数を維持。カウンタ・加算を削らない | 探索に隣接 | 保留 |
| D13 | `src/forceCalculator.cpp:213,691,1068,1160` setContactMonitor、monitor統計、contactForceL/R_ij | 接触力配分と監視 | addSegmentForceXY、他のcontact CSV | REVIEW_ONLY | 配列・API共用。D09削除後も保持 | 高 | 不可 |
| D14 | `src/forceCalculator.cpp:629等` PROFILE_CALCDIS | コンパイル時ONの計測 | CMake/AGENTS | REVIEW_ONLY | 明示的な保持対象 | 低〜中 | 不可 |
| D15 | `src/geometry.cpp:546,688,1233,1265` CSV/DAT/表示 | 初期化・幾何検査 | surfp.py、初期化のreturn経路 | REVIEW_ONLY | 許可リスト外。幾何と表示を分離しない | 中 | 不可 |
| D16 | `src/rungeStep.cpp:5` rungeStep、`geometry.cpp:457` surfExtract、`forceCalculator.cpp:1828` findNsep、`ModeData.h:32` loadAll | 定義/宣言あり、現在の呼び出しは検出されず | 公開ヘッダー、CMakeに実装を含む | REVIEW_ONLY | 定義・宣言・全文検索確認。今回の削除許可外 | 未調査の外部利用 | 不可 |
| D17 | `src/forceCalculator.cpp:351,590` pressureRecoveryClampCount等 | クランプ統計・非有限フラグ | simulation.cpp:1145 flow_diagnostics_v2.csv | KEEP_OUTPUT | 他の維持出力・異常判定へ参照あり | 高 | 不可 |
| D18 | `include/FlowModel.h:29` findSeparation、同:17–18 pressure/area、`src/surfaceLoad.cpp:13` clear | findSeparationは宣言のみ、他は定義あり。現行呼び出し検出なし | 公開/非公開API、保存メンバー | REVIEW_ONLY | 許可リスト外。宣言や未使用に見える保存域を独断で削除しない | API・状態保持 | 不可 |

## 4. 書き込みファイルの全目録（C++実行系）

通常の保存先はparam親の親の `output/`、VTUは `result/`。環境変数 `SIMULATION_RUN_DIR` 設定時は `runDir` と `runDir/result`（幾何の直接相対出力は例外）。単位「generalized」は質量正規化した一般化座標・力・仕事であり、N/J等と断定しない。CSVはヘッダーに個別単位がある。

| ID | ファイル／保存先 | 書込元（整理前） | 条件・頻度 | 内容・単位 | 解析参照・判定 |
|---|---|---|---|---|---|
| O01 | `txt/params_used.txt` | simulation.cpp:80 initialize/copy_file | 毎初期化1回 | 使用param原文 | dynamics_common、left_right_ratios / KEEP_OUTPUT |
| O02 | `txt/manifest.txt` | simulation.cpp:82 initialize | 毎初期化1回 | 条件、パス、単位・設定 | dynamics_common、analyze_modal_exit_test、スイープ / KEEP_OUTPUT |
| O03 | `csv/fixed_boundary_mode_audit.csv` | simulation.cpp:224 initialize | fixedNodeIdsFile指定時、左右各モード | 固定点と全点の最大モード値・比 | AP診断参照文書 / KEEP_OUTPUT |
| O04 | `dat/area.dat` | simulation.cpp:351,1170 run | n%5==0 | step、50断面面積mm² | show/freq/jitter/dynamics/AP/sweep等 / KEEP_OUTPUT |
| O05 | `dat/displace.dat` | simulation.cpp:352,1192 run | n%5==0 | 秒、左右uy変位mm | show2/show_displace/ratios/sweep/self_oscillation / KEEP_OUTPUT |
| O06 | `dat/pressure.dat` | simulation.cpp:353,1176 run | n%5==0 | step、50断面圧力Pa | pressure、AP解析 / KEEP_OUTPUT |
| O07 | `dat/pressure_vt.dat` | simulation.cpp:354,1200 run | n%5==0 | step、出口圧Pa | pressure_vt/fft/fft_vt/dynamics / KEEP_OUTPUT |
| O08 | `dat/airflow_vt.dat` | simulation.cpp:355,1201 run | n%5==0 | step、流量m³/s | airflow/close/spectol/sweep/dynamics / KEEP_OUTPUT |
| O09 | `dat/section_x.dat` | simulation.cpp:356,1177 run | n%5==0 | step、断面位置mm | 形状説明用、許可リスト外 / KEEP_OUTPUT |
| O10 | `dat/gap_cubed.dat` | simulation.cpp:357,1178 run | n%5==0 | step、∫g₊³dz mm⁴ | 粘性・間隙の説明用 / KEEP_OUTPUT |
| O11 | `dat/separation.dat` | simulation.cpp:358,1184 run | n%5==0 | step、index、x/blend mm、Pa、面積mm²、fallback | 分離点解析 / KEEP_OUTPUT |
| O12 | `dat/displace_xy.dat` | simulation.cpp:359,1195 run | n%5==0 | 秒、uxL/uyL/uxR/uyR mm | plot_displacement_xy / KEEP_OUTPUT |
| O13 | `csv/modal_contribution.csv` | simulation.cpp:360,552,1212 run | n%nwrite==0、左右全モード | Hz、q/qdot、probeとsurface RMS mm | analyze_dynamics/plot_modal_energy_share/dynamics_common / REVIEW_ONLY（保持） |
| O14 | `csv/modal_dominant.csv` | simulation.cpp:361,582 run | n%nwrite==0、左右 | probe/surface優勢モード・mm | run_pressure_sweep保存 / REVIEW_ONLY（保持） |
| O15 | `csv/modal_top10.csv` | simulation.cpp:362,601,1214 run | n%nwrite==0、左右×3scope×最大10 | モード寄与mm、比、相殺指標 | plot_modal_top10 / REVIEW_ONLY（保持） |
| O16 | `csv/irregularity_timeseries.csv` | simulation.cpp:363,1202 run | n%5==0 | 秒、左右変位mm、面積mm²、流量m³/s、圧Pa、接触 | analyze_irregularity/self_oscillation / KEEP_OUTPUT |
| O17 | `csv/flow_diagnostics_v2.csv` | simulation.cpp:364,1107 run | n%diagnosticOutputIntervalSteps==0 | 秒、Pa、m³/s、m³/s²、mm²、分離・異常・クランプ | analyze_self_oscillation/AP / KEEP_OUTPUT |
| O18 | `csv/modal_energy.csv` | simulation.cpp:365,1560 run | 診断間隔、左右全モード | q/qdot、一般化力・パワー・エネルギー | 研究用エネルギー出力 / KEEP_OUTPUT |
| O19 | `csv/energy_summary.csv` | simulation.cpp:366,1594 run | 診断間隔 | 左右・総エネルギー、仕事、減衰、接触 | analyze_self_oscillation / KEEP_OUTPUT |
| O20 | `csv/gap_field.csv` | simulation.cpp:367,767 run | initial、基準窓終端、擾乱直後 | 共通x/z、左右y、符号付き間隙mm | AP解析照合 / KEEP_OUTPUT |
| O21 | `csv/gap_field_summary.csv` | simulation.cpp:370,804 run | O20と同時、i毎 | 正間隙面積mm²、最小間隙mm、負間隙率 | AP診断研究 / KEEP_OUTPUT |
| O22 | `csv/modal_force_split.csv` | simulation.cpp:374,1517 run | 診断間隔、左右各モード | q/qdot/qddot、流体・接触・復元・減衰・残差 | summarize_modal_force_split/AP/modal_exit / KEEP_OUTPUT |
| O23 | `csv/diagnostic_events.csv` | simulation.cpp:381,1029,1056 run | 擾乱/荷重切替イベント | step、秒、説明 | 再現性 / KEEP_OUTPUT |
| O24 | `csv/perturbation_energy.csv` | simulation.cpp:383,1412,1479 run | 参照標本あり、切替時刻以後、診断間隔 | dq/dv、一般化仕事、残差、平衡残差 | analyze_ap_contact/modal_exit / KEEP_OUTPUT |
| O25 | `csv/perturbation_energy_summary.csv` | simulation.cpp:391,1505 run | O24条件、左右集計 | 秒、一般化エネルギー・仕事 | 研究用 / KEEP_OUTPUT |
| O26 | `csv/contact_coupling_debug.csv` | simulation.cpp:413,1334 run | 接触反復毎 | モード力差、予測変位差mm、めり込みm、収束残差 | run_pressure_sweep保存 / KEEP_OUTPUT |
| O27 | `csv/debug_step_summary.csv` | simulation.cpp:420 | ヘッダーのみ | 各種診断列（行なし） | D01 / REVIEW_ONLY |
| O28 | `csv/debug_mode_summary.csv` | simulation.cpp:431,1256 | 反復中n%10==0またはt>0.12 | fiL/R最大値・index、接触差 | D02 / REMOVE_SAFE |
| O29 | `txt/perturbation_manifest.txt` | simulation.cpp:1032 run | 擾乱を実際に適用時1回 | 時刻・モード・delta_q・probe mm | 再現性 / KEEP_OUTPUT |
| O30 | `wav/test_sound.wav` | simulation.cpp:1636、wavwrite.h:69 | 終了時1回、毎stepサンプル | 16bit mono、Fs=round(1/dt)Hz | wav_power_spectrum/compare（パス指定） / KEEP_OUTPUT |
| O31 | `result/deform_combinedNNNN.vtu` | simulation.cpp:1624,1663 | n%20==0 && 残り<=5000 | 左右変形座標mmとセル情報 | 可視化用 / KEEP_TEST |
| O32 | `../output/csv/surfp_output.csv` | geometry.cpp:1233 surfExtractFromNAS、:546 surfExtract | 左右の表面抽出時（同じパスで右が上書き） | i/j、0-based VTU点ID、原形座標mm | surfp.py / REVIEW_ONLY |
| O33 | `../output/dat/SurfArea.dat` | geometry.cpp:688 surfArea | 左右面積構築時（右が上書き） | 表面格子面積mm² | 名前直接のtools参照なし、保持指定 / REVIEW_ONLY |
| O34 | `txt/debug_force.txt` | forceCalculator.cpp:112,324 | 初期化・流体毎step | step、面積mm²、圧Pa | D06 / REMOVE_SAFE |
| O35 | `csv/contact_debug.csv` | forceCalculator.cpp:117,1583 | 接触反復かつ0.120<=t<=0.135 | 件数、めり込みm、圧Pa、力N、法線等 | D08 / REVIEW_ONLY |
| O36 | `csv/contact_monitor.csv` | forceCalculator.cpp:129,1677 | 有効monitorかつ非diagnosticOnly呼出し | 間隙mm、接触力N、残差 | D09 / REMOVE_SAFE |
| O37 | `csv/contact_iteration_debug.csv` | forceCalculator.cpp:139,1700 | 非diagnosticOnlyの各反復 | 接触・圧力・力・最大めり込み等 | D10 / REVIEW_ONLY |
| O38 | `csv/contact_search_debug.csv` | forceCalculator.cpp:153,1717 | 非diagnosticOnlyの各反復 | 候補・cache・fallback件数 | D11 / REMOVE_SAFE |
| O39 | `csv/contact_final_state.csv` | forceCalculator.cpp:160,1554 | diagnosticOnly（最終状態）の毎step | 接触件数、最大めり込みm、圧Pa、力N | run_pressure_sweep保存 / KEEP_OUTPUT |
| O40 | `csv/contact_region_summary.csv` | forceCalculator.cpp:168,1536 | AP診断ON、diagnosticOnly、診断間隔、3領域 | 候補/接触、mm、Pa、N、間隙率 | analyze_ap_contact / KEEP_OUTPUT |
| O41 | `csv/contact_pairs_detail.csv` | forceCalculator.cpp:175,1330 | AP診断ON、diagnosticOnly、指定detail時間内 | 節点ID、位置mm、間隙、法線、圧Pa、力N、除外 | analyze_ap_contact / KEEP_OUTPUT |
| O42 | `csv/debug_flow_detail.csv` | forceCalculator.cpp:184,379,574 | 流体n%10==0またはt>0.12、異常値も追記 | mm²/m²、流量、圧Pa、クランプ等 | D07 / REVIEW_ONLY |

コメントのみの `debug_fluid.dat`、旧contact_iteration_debugストリーム、初期VTU呼び出しは実行されない。コメントを一括削除しない。C++の全 `ofstream`/`.open`/`copy_file`/WavWriter/VTU経路を調査した。

## 5. tools の全スクリプト依存表

`tools/`直下は43本のPythonスクリプト、`tools/tests/`は8本。下表の互換性は読み込みコードと現行出力仕様の静的照合であり、全スクリプトの実行成功を意味しない。全て保持。類似用途は統合許可ではない。

| スクリプト（参照行） | 主な読込ファイル・依存 | 目的 | 現行出力との両立・類似用途 |
|---|---|---|---|
| airflow.py:8 | dat/airflow_vt.dat | 流量波形 | 形式対応、dt固定。close/spectolと類似 |
| analyze_ap_contact.py:81,145–147,186,231,276 | modal_force_split、flow_diagnostics_v2、contact_region_summary/pairs_detail、perturbation_energy、area、pressure | AP端接触・仕事・仮説比較 | 対応。AP診断を有効にした出力が必要。modal_exitと関連 |
| analyze_dynamics.py:14,239 | dynamics_common経由area/displace/airflow/pressure/manifest/params、modal_contribution | 周期・変調・Poincaré、左右比 | 対応。analyze_irregularity、ratiosと一部類似 |
| analyze_irregularity.py:42 | irregularity_timeseries.csv | 非周期指標、位相・スペクトル | 対応。analyze_dynamicsと関連 |
| analyze_modal_exit_test.py:100,103,701,718 | manifest、run_manifest.json、modal_force_split、perturbation_energy | live/freeze_fluidの一般化仕事比較 | 対応。長い基準/擾乱窓が必要 |
| analyze_self_oscillation.py:162–194 | displace、area、irregularity、flow_diagnostics_v2、energy_summary、debug_flow_detail fallback | 自励応答・成長率 | 対応。D07の読み手を確認 |
| close.py:6 | dat/airflow_vt.dat | 流量FFT・窓解析 | 形式対応、時間窓固定。spectolと類似 |
| compare_stiffness_metrics.py:14–23,128,166 | displace、manifest、ratios・steinecke_classification | 剛性別jitter/shimmer/周波数・振幅 | 対応。複数runを比較 |
| compare_wav_power_spectrum.py:26,33,153 | 指定した2WAV（既定simulation.wav/subject2.wav） | Welchスペクトル比較 | test_sound.wavは明示パス指定が必要。wav_power_spectrum2と類似 |
| dynamics_common.py:34,80,234–251 | manifest/params、area、airflow、pressure_vt、displace、modal_contribution、output_layout | 解析共通読込・信号処理 | 対応。analyze_dynamics/run_pressure_sweep/compare_stiffnessがimport |
| fft.py:11 | dat/pressure_vt.dat | 出口圧FFT | 形式対応、dt/窓固定。fft_vtと類似 |
| fft_vt.py:8 | dat/pressure_vt.dat | 出口圧FFT・比較表示 | 形式対応、dt/窓固定。fftと類似 |
| freq.py:7 | dat/area.dat | 面積ピーク周波数・FFT | 形式対応、時間窓固定。analyze_dynamicsと類似 |
| jitter.py:8 | dat/area.dat | 最小面積の局所変動指標 | 形式対応、dt/窓固定。周期jitterと定義は同一ではない |
| kcExplore.py:9–12,118 | 旧絶対パスparam/displace/area/debug_step_summary | 接触係数探索・診断値参照 | 旧パス・現行debug_stepはヘッダーのみ。参照があるのでD01削除保留。paramを書き換えるため今回実行しない |
| left_right_ratios.py:96,280,486 | displace、params、スイープのleft_right_ratio_metrics | 左右の周期・振幅比 | 対応。compare_stiffness/spectol_sweepと関連 |
| monitor_point.py:12,36 | output_pair_offset/monitor_point.csv | x/z変位表示 | 現行C++に生成元なし。contact_monitor.csvとは別物。保持 |
| output_layout.py:11,20 | 拡張子別canonical/旧flatパス | 入出力パス互換 | 多数の解析がimport。legacy比較ソルバではない |
| plot_bifurcation.py:51–56 | bifurcation_points.csv、sweep_summary.csv | 分岐図 | run_pressure_sweepの派生出力と対応 |
| plot_displacement_xy.py:24,46 | dat/displace_xy.dat | 左右xy波形 | 対応。show2/show3と関連 |
| plot_modal_energy_share.py:126 | modal_contribution.csv、output_layout | モードエネルギー比 | 対応。modal_top10とは指標が異なる |
| plot_modal_top10.py:32,76 | csv/modal_top10.csv | 平均変位寄与の上位3モード | 対応。エネルギー比とは別 |
| plot_pair_offset_summary.py:14,47 | output_pair_offset/pair_offset_summary.csv | 対応位置差の最大値 | 現行C++に生成元なし。保持 |
| plot_poincare.py:35 | poincare_points.csv、output_layout | Poincaré図 | analyze_dynamics派生出力と対応 |
| pressure.py:4 | dat/pressure.dat | 圧力の第1断面表示 | 形式対応、dt固定 |
| pressure_vt.py:4 | dat/pressure_vt.dat | 出口圧波形 | 形式対応、dt固定。fft系と関連 |
| run_ap_contact_diagnostics.py:24,46,136 | paramとそのモデル入力、analyze_ap_contactを起動 | AP除外・荷重凍結の隔離A/B | 現行オプション対応。今回実行しない |
| run_diagnostic_sweep.py:31,146,195 | param、モデル入力、analysis_summary.json | 圧力/initialGapの2段階探索 | 現行オプション対応。analyze_self_oscillationを起動 |
| run_pressure_sweep.py:27–33,109,135,219 | param、config、基本出力、modal、contact_iteration/final/coupling、metrics等 | 圧力掃引・保存・解析 | 対応。D10は保存参照があるため保持 |
| run_right_stiffness_sweep.py:18,69–79,147 | param、左右モデル、M5_freq_T3_d2_b*c*.txt/M5_mode_T3_b*c*.vtu | 右剛性モデル掃引 | 入力が存在する条件で利用。compare_stiffnessを起動 |
| show.py:9 | dat/area.dat | 最小声門面積波形 | 形式対応、dt固定。show_area_lastと類似 |
| show2.py:5 | dat/displace.dat | 左右uy波形 | 対応。変数dispXも右uy列を読む（修正しない） |
| show3.py:6–7 | 旧displace_c5.dat、displace2_c5.dat | 2観測点の変位軌跡 | 現行はこの2ファイルを生成しない。保持 |
| show_area_last.py:12,46 | 指定area.dat | 終盤の最小面積 | 対応。showと類似 |
| show_displace.py:14–16,51 | param、dat/airflow_vt.dat、dat/displace.dat | 変位波形 | 対応。パス/cwd前提あり |
| spectol.py:12,21 | dat/airflow_vt.dat | 流量スペクトル | 形式対応、窓固定。closeと類似 |
| spectol_sweep.py:18–21 | 現リポジトリ絶対パスparam、airflow/displace/area、steinecke_classification | 圧力等の掃引・波形・locking | 現行出力形式対応。paramを書き換えるため今回実行しない |
| steinecke_classification.py:1 | 配列引数（直接ファイル読込なし） | 左右周波数locking判定 | spectol_sweep/compare_stiffnessがimport |
| summarize_modal_force_split.py:148,162,216 | 指定modal_force_split.csv | chunk集計・選択モード抽出 | 対応。modal_exitと関連 |
| surfp.py:14–15,26,59 | csv/surfp_output.csv | 表面格子と監視点表示 | 形式対応。ただし既定targetX=9.6、現行本体は9.2。変更せず提案 |
| wav_power_spectrum.py:142 | 指定WAV | 片側パワースペクトル | test_sound.wavに対応 |
| wav_power_spectrum2.py:177 | 指定WAV | Welch平均スペクトル | 対応。前者と推定法が異なる |
| zeta_sweep.py:39–43,116,184 | 旧Mode_superposition_2026絶対パスparam/displace | 減衰係数掃引・左右比 | パスが現リポジトリと異なる。paramを書き換えるため今回実行しない |

toolsはCSV/DATを読み、派生CSV/JSON/PNG、スイープparamとログを書き出す。C++出力とは別の派生成果物であり、削除対象にしない。依存は全スクリプトのAST（ファイル文字列・ローカルimport）と読込箇所を照合した。

## 6. テスト・既存不整合・保留

- `tests/test_separation_selection.cpp`: 面積比、exact閾値、fallback、不正/NaN断面、閉鎖。CTest名 `separation_selection`。
- `tests/test_mode_selection.cpp`: 既定・単一・複数モード選択、sourceModeIndex。CTest名 `mode_selection`。Releaseフラグは `-O3` で、CMakeのRelease既定を `-DNDEBUG` に変更する処理はない。
- `tests/test_input_file_params.cpp`: 左右入力パスとAP診断既定値。CTest名 `input_file_params`。**整理前から失敗**。:34–35は0.15/0.19秒を期待するが、`include/SimulationParams.h:79–80`は0.25/0.30秒。paramの対応行はコメント。両者とも変更禁止。ユーザー許可に従い、整理後も成功2件＋同原因の失敗1件だけを許容する。
- `tools/tests/` は8本: `test_analyze_irregularity.py`、`test_analyze_modal_exit_test.py`、`test_compare_stiffness_metrics.py`、`test_dynamics_common.py`、`test_left_right_ratios.py`、`test_modal_energy_share.py`、`test_right_stiffness_sweep.py`、`test_spectol_sweep.py`。C++診断streamの削除とは独立した解析関数のテスト。今回のコード整理でtoolsは変更しない。
- 専用のlegacy/current二重ソルバ比較実装は検出しなかった。実在する比較はWAV、剛性条件、AP A/B、live/freeze、入力・出力の旧形式互換である。名前だけで除去しない。
- 基準実行は隔離環境で `./simulation ../input/fluid.param` と `./simulation ../input/sine.param`。nstepだけ2000へ短縮、sineのみ実効iforce=1。dt、全100モード、zeta=0.02、ncont=100、kc、Ps、famp/forcef、形状を保持。OpenMPは4 threads、close、cores。
- 基準・整理後の比較と完了状態は変更報告に記載する。短時間試験は0.02秒であり、0.120–0.135秒のcontact_debug、長期自励/非周期応答、擾乱、荷重凍結、声道あり、PROFILE_CALCDIS=ONは実行カバレッジ外。関連削除は保留。
- 基準計算・全VTU構造検査は両モードとも完了。DAT各400行、VTU各100ファイル、WAV各2000サンプル/100000Hz/mono/16bit。最終状態の接触ありstepは通常802、正弦波1964。主な物理時系列・全変形座標は有限。診断欄の既存NaN（初回の変位差、無接触の最大位置等、ゼロ変位時の比、正弦波の未更新分離状態）は変更報告に明示する。
- Phase 3の構造変更は行わない。将来候補は診断出力の責務分離、旧解析のパス/列仕様整合、幾何直接相対出力の保存先整理。いずれも今回実装しない。

Phase 2ではD02/D03/D04/D06/D09/D11を実施し、両加振モードで最終維持ファイル各133件がバイト一致した。別途判断が必要な候補: D01/D07/D08/D10/D12/D13/D14/D15/D16/D18。現物の読み手や共用計算があるものは保持した。詳細・比較条件・未検証事項は `cleanup_change_report.md` を参照。
