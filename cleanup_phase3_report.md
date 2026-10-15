# Phase 3 CSV・擾乱診断・指定旧関数の整理報告

作業日: 2026-10-15（JST）。指示書は `CODEX_cleanup_phase3_remove_csv_and_legacy.md`。指定のソース整理後、正弦波のVTU座標に数値不一致を検出したため安全停止した。**整理は未合格・未完了であり、同一結果を保証できない。** 原因未確定のため、追加のコード変更・再計算・許容誤差緩和は行わず指示を待つ。過去の `cleanup_inventory.md` / `cleanup_change_report.md` は履歴として変更しない。

## 基準状態・隔離・保護

- 基準HEADは `1eb8e23c0a27e56ced75a39c1fd8ff038636885f`。開始時のstatus・通常diff・cached diffは全て空。GitHubの参照コミットではなく、この実作業ツリーを基準にした。
- 本番 `input/param.txt` は nstep=100000、zetaL=zetaR=0.02、nmode=100、dt=1e-5、Ps=2900、kc=(5,6000,5000)、ncont=100、forcef=100、famp=-1e-5、実効iforce=0。全入力を変更しない。
- 削除対象の名前付きオプションは、本番paramでは全てコメントアウトされており、有効な指定はない。既存コメントは入力保護のため保持する。これらを将来有効化すると、削除後は `unknown named parameter` のエラーになる。黙って読み飛ばす処理は追加しない。
- 隔離ディレクトリは `/tmp/simulation-phase3-Wmk6qH`。開始時のdiff、cached diff、trackedファイル/全inputのSHA256、元output/resultのSHA256・サイズ・mtime、研究アーカイブのサイズ・mtime、コンパイラ・OpenMP条件を `start.json` に保存。
- ソース・ヘッダー・tests・CMakeの編集前コピーは `original/`、検証ビルドは `work/build`。テストごとに `cases/fluid/build` / `cases/sine/build` を固定cwdとし、それぞれ独立したinput/output/resultを使う。幾何の `../output` 直接書込も隔離される。`SIMULATION_RUN_DIR` はテスト環境から除外。
- 元のinput55ファイル、output/result351ファイル、analysis_runs/diagnostic_runsの研究データは削除しない。元のbuildも検証には使わない。ディレクトリ移動・rename・依存追加・commit/push・git reset/restore/checkout/stash/cleanは行わない。
- 停止後に再照合し、全input55ファイルと元output/result351ファイルのSHA256・サイズ・mtimeは開始時と全て一致した。研究アーカイブのサイズ・mtimeも一致（内容ハッシュは対象外）。`protection.json` に結果を保存。本番paramのSHA256は `0e4a085c96f028f8935b9701d00af831c82ed61c96f4e942f27fae4e50ae7778`。変更ソース14件は隔離ビルドのソースと一致し、許可対象以外の既存trackedファイルの変更はない。

## C++ CSVの削除一覧

下記の場所は編集前コピーの行番号。単に書込を無効化するのでなく、CSV専用stream、ヘッダー、行、集計・中間状態・設定を除去する。

| CSV | 編集前の書込元 | 対応 |
|---|---|---|
| fixed_boundary_mode_audit.csv | simulation.cpp:224 initialize | 固定点監査と専用設定を削除済み |
| modal_contribution.csv | simulation.cpp:360 / writeModalDiagnostics | probe/surface RMSと専用処理を削除済み |
| modal_dominant.csv | simulation.cpp:361 / writeModalDiagnostics | 優勢モード集計を削除済み |
| modal_top10.csv | simulation.cpp:362 / writeTop10ModalContributions | 上位寄与・射影比・相殺比と専用構造体を削除済み |
| irregularity_timeseries.csv | simulation.cpp:363 / run | CSV専用の面積極値・観測変位の集計を削除済み。DAT観測点は保持 |
| flow_diagnostics_v2.csv | simulation.cpp:364 / run | CSV専用の診断走査・getter・統計を削除済み |
| modal_energy.csv | simulation.cpp:365 / run | CSV専用のモードエネルギーを削除済み |
| energy_summary.csv | simulation.cpp:366 / run | CSV専用の総仕事・減衰集計を削除済み |
| gap_field.csv | simulation.cpp:367 / writeGapSnapshot | 専用snapshot・初期/基準窓/擾乱イベント処理を削除済み |
| gap_field_summary.csv | simulation.cpp:370 / writeGapSnapshot | 専用間隙面積・負間隙率を削除済み |
| modal_force_split.csv | simulation.cpp:374 / run | 専用モード力分解と残差の出力を削除済み。接触収束残差は別に保持 |
| diagnostic_events.csv | simulation.cpp:381 / run | 擾乱・荷重切替イベントを削除済み |
| perturbation_energy.csv | simulation.cpp:383 / run | 擾乱仕事・エネルギー評価を削除済み |
| perturbation_energy_summary.csv | simulation.cpp:391 / run | 擾乱集計を削除済み |
| contact_coupling_debug.csv | simulation.cpp:413 / run | fi差・予測変位差snapshotと専用ラムダを削除済み |
| debug_step_summary.csv | simulation.cpp:420 / run | ヘッダーのみのstreamを削除済み |
| contact_debug.csv | forceCalculator.cpp:111 initialize / applyContactLoads | 0.120–0.135秒の診断統計・最悪値を削除済み |
| contact_iteration_debug.csv | forceCalculator.cpp:123 / applyContactLoads | 専用接触位置・圧力・力の統計を削除済み |
| contact_final_state.csv | forceCalculator.cpp:137 / applyContactLoads | CSV生成を削除済み。diagnosticOnlyの計算・副作用は保持 |
| contact_region_summary.csv | forceCalculator.cpp:145 / applyContactLoads | AP分析領域・専用統計を削除済み。接触除外は保持 |
| contact_pairs_detail.csv | forceCalculator.cpp:152 / applyContactLoads | 詳細ペア出力・時間窓を削除済み |
| debug_flow_detail.csv | forceCalculator.cpp:161 / applyFluidLoads | CSV専用診断値を削除済み。圧力回復クランプ本体は保持 |
| surfp_output.csv | geometry.cpp:1233 surfExtractFromNAS | **維持**。パス・ヘッダー・内容・左右での生成順を変更しない |

旧 `Geometry::surfExtract` にも同名surfpライター（geometry.cpp:546）があるが、指定旧関数の削除に含める。NAS側とSurfArea.datは全保持。旧phaseで既に削除したdebug_mode/debug_force/contact_monitor/contact_searchは再導入しない。

## 削除機能と保持する共用計算

削除対象は、擾乱のmodal q加算、擾乱manifest、freeze_fluid/contact/all、基準窓の集計・平均・荷重適用、擾乱エネルギー/仕事/平衡残差、CSV専用のprojectFluidLoadsToModes・fiFluid/fiContactの分解バッファ・projectPair、CSV保存先と監視setterフック。

専用パラメータは `diagnosticOutputIntervalSteps`、`perturbationEnabled/TimeSec/ModeIndex/AmplitudeAtProbeMm/Pattern`、`diagnosticLoadMode`、`diagnosticReferenceWindowStartSec/EndSec`、`fixedNodeIdsFile`、`apContactAnalysisRowsPerEnd/DistanceMm`、`contactDetailStartSec/EndSec`。宣言・読込・検証・表示・manifestを揃えて除去する。

保持する設定は旧/実効iforce、forceDirection、contactReferenceFrequencyHz、flowBlendLengthMm、flowSeparationAreaRatio、initialGapMm、pressureRampTimeSecと省略互換、左右モデル設定、`apContactDiagnosticEnabled` と `apContactExcludedRowsPerEnd/DistanceMm`。最後の3項目は実際の接触力除外条件なので保持する。

| 指定旧関数 | 編集前の場所 | 参照調査・処理 |
|---|---|---|
| TimeIntegrator::rungeStep | src/rungeStep.cpp:5、include/TimeIntegrator.h:8 | 定義/宣言のみ。削除済み。同CPPとnewmarkStepを保持 |
| Geometry::surfExtract | src/geometry.cpp:457、include/Geometry.h:44 | 定義/宣言のみ。削除済み。surfExtractFromNAS/surfAreaを保持 |
| ForceCalculator::findNsep | src/forceCalculator.cpp:1768、include/ForceCalculator.h:158 | 定義/宣言のみ。削除済み。現行面積比分離選択を保持 |
| ModeData::loadAll | include/ModeData.h:32 | inline定義のみ。削除済み。個別読込を保持 |
| FlowModel::findSeparation | include/FlowModel.h:29 | 宣言のみ。削除済み。FlowModelの計算と現行分離処理は保持 |
| SurfaceLoad::clear | src/surfaceLoad.cpp:13、include/SurfaceLoad.h | 定義/宣言のみ。削除済み。resizeと保存配列を保持 |

`rg`でsrc/include/tests/toolsの呼び出し・参照を調べ、6件に呼び出しは検出されなかった。他の未使用に見える関数は追加削除しない。CSV専用のラムダ・setter/getter・分解用メソッドの除去は、指示書3/4節の対象である。

`applyContactLoads(..., true)` は残す。内部はfdis/contactForce配列をゼロ化・再計算し、6個の接触フラグ/残差/最大めり込み/前回めり込みスカラーを保存・復元してreturnするが、全配列を元に戻す処理ではない。従って純粋なファイル出力と同一視せず、元と同じ時点で呼び出す。diagnosticOnly分岐、cache更新条件、配列への加算、スカラー復元を維持し、CSV専用部分だけ削除する。

保持する物理計算は、流路/面積/流量/圧力/分離/表面力、接触候補の探索順とcache・fallback、法線・ペナルティ・減衰・力分配、maxPenetration更新、前回力の更新と収束残差、モード射影、Newmarkと状態確定、VTU/WAV/DAT。モード射影からNewmarkへのfiL/R受け渡し、β=0.275625・γ=0.55、全OpenMP pragmaを変更しない。CSV候補件数カウンタは物理分岐に影響しない加算・ログ引数だけを除去し、processPair呼び出しやcache条件は保持する。

## ビルドと数値検証

IntelLLVM/icpx 2025.1.0、CMake 3.26.4、Release、`-O3 -std=gnu++17 -fiopenmp`、PROFILE_CALCDIS=OFF、BUILD_TESTING=ON。既存CMakeを変更しない。

```bash
# cwd: /tmp/simulation-phase3-Wmk6qH/work/build
cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 4
ctest --output-on-failure

# 各caseの固定cwd: /tmp/simulation-phase3-Wmk6qH/cases/{fluid,sine}/build
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
  /tmp/simulation-phase3-Wmk6qH/work/build/simulation ../input/test.param
```

| 項目 | 整理前 | 整理後 |
|---|---|---|
| Releaseビルド | 成功 | 成功 |
| separation_selection | 成功 | 成功 |
| mode_selection | 成功 | 成功 |
| input_file_params | 既存失敗（診断窓の期待値不整合） | 成功（削除設定の期待値のみ除去） |

既存失敗は0.15–0.19秒の期待値と0.25–0.30秒の既定値の不整合。今回は指示書6節に従い、削除するAP分析行・loadMode・基準窓の期待値だけをテストから除去する。左右入力パス、接触除外設定のテストは保持し、skip/無効化は行わない。

整理前CTestの終了コードは8（成功2件・既存失敗1件）、整理後は0（3件全成功）。PROFILE_CALCDIS=ONも別の隔離ビルド `profile-build` でコンパイルとCTest3件の成功を確認した。ONでのシミュレーション数値比較は実施していない。削除オプション14種類を一つずつ有効にした一時paramをパーサーテストへ渡し、全て終了コード1と `unknown named parameter` の明示エラーを確認した（`retired-options.json`）。

同一caseの整理前後で、入力ファイル、全100モード、dt/Ps/zeta/kc/初期間隙/力の条件、コンパイラ・フラグ、cwd・OpenMPを固定する。試験paramの変更はnstep短縮と正弦波側の実効iforce=1のみ。本番paramは変更しない。

| 条件 | ステップ／時間 | 基準計算 | 整理後計算・比較 |
|---|---|---|---|
| 通常流体加振 iforce=0 | 13620／0.1362秒 | 終了コード0、VTU250 | 正弦波の不一致を受け途中中断（実行セッション終了コード130）。VTU149の部分出力が残る。比較未完了・合否未判定 |
| 正弦波加振 iforce=1、direction=0 | 2000／0.02秒 | 終了コード0、VTU100 | 計算終了コード0、VTU100。ただし97ファイルの座標が不一致のため回帰失敗 |

基準の通常加振では、contact_debug.csvにstep12000–13500の3866行があり、0.120–0.135秒の実行を確認した。接触反復31422回、非接触1198step、接触収束12422step、最大反復3、上限到達0。正弦波は反復5672回、非接触36step、接触収束1964step、最大反復3、上限到達0。両ケースで接触と非接触を通過した。

## ソースの差分・静的監査

| ファイル | 追加行 | 削除行 |
|---|---:|---:|
| include/FlowModel.h | 0 | 1 |
| include/ForceCalculator.h | 0 | 33 |
| include/Geometry.h | 0 | 1 |
| include/ModeData.h | 0 | 5 |
| include/SimulationParams.h | 0 | 14 |
| include/SurfaceLoad.h | 0 | 1 |
| include/TimeIntegrator.h | 0 | 3 |
| src/forceCalculator.cpp | 8 | 567 |
| src/geometry.cpp | 0 | 123 |
| src/readParams.cpp | 2 | 97 |
| src/rungeStep.cpp | 0 | 20 |
| src/simulation.cpp | 0 | 1004 |
| src/surfaceLoad.cpp | 0 | 6 |
| tests/test_input_file_params.cpp | 1 | 5 |
| 合計（14ファイル、レポートを除く） | 11 | 1880 |

追加行はログ専用の戻り値保存/カウンタ引数の除去に伴う構文、残すパラメータ検証条件とテストの括弧、不要になった診断へのコメント言及の修正であり、物理式の新規追加はない。ソース監査（`source-audit.json`）では、全OpenMP pragma、VTU関数本体と呼び出し条件、DAT書込ブロック、左右Newmark/最終Newmark、正弦波分岐、contactIsExcluded、接触力と収束残差の更新ブロックを編集前と照合した。NAS/面積は旧surfExtractより後の実装全体、newmarkStepの本体も同一である。CMake、flowModel.cpp、channelSections.cpp、separationSelection.cpp、modalProjector.cpp、modedata.cpp、state.cpp、wavwrite.hは未編集。

比較対象は全DAT（SurfArea含む）、WAV、全VTU、surfp_output.csv、params_used.txtと保持manifest項目。ファイルのバイナリ比較は `filecmp.cmp(..., shallow=False)`、許容誤差は絶対/相対とも0。manifestだけは許可された削除設定17項目の除去を期待差分とし、他の行の変更を許容しない。正弦波のファイル集合では、基準で生成されていた削除CSV19本だけが消え、予期しない追加/欠落はない。未有効の診断CSV3本を含め、ソースでは22種類の生成箇所を除去した。

基準VTUは全ファイルをXMLとして読み、90666点・496536セル、全座標の有限性、connectivity範囲、offsets末尾・types数、トポロジーを検査した。整理後の正弦波は全100ファイルを比較したが全バイト一致は達成していない。CSVを撤去するため、整理後に全接触力配列を別ファイルへ新規出力する実装は加えていない。

正弦波では流体分離状態を更新しないため、維持するseparation.datの一部に元からNaNがある。これも位置・表記を含め同一であることを要求し、新しいNaN/Infや差分を許容しない。

## 回帰不一致と安全停止

最初の比較処理は `result/deform_combined0014.vtu` のバイト不一致で停止した。原因を推測して修正せず、実行中だった通常加振も中断。その後は保存済みの正弦波出力の読取照合と保護確認だけを行った。再計算はしていない。証拠は隔離先の `sine-difference.json`、`measure_difference.py`、`final-sine.log/json`、`final-fluid.log`、`protection.json` に残す。

| 保持出力／検証（正弦波） | 結果 | 最大絶対差／相対差 |
|---|---|---|
| DAT10本（SurfArea、airflow_vt、area、displace、displace_xy、gap_cubed、pressure、pressure_vt、section_x、separation） | 全てバイト一致。各時系列400行の基準列数・時刻軸・既存NaN位置も維持 | 保存された有限数値は0／0 |
| test_sound.wav | バイト一致。mono/16bit/100000Hz、2000サンプル | サンプル差0 |
| surfp_output.csv、params_used.txt | バイト一致。再帰探索でCSVはsurfpの1本のみ | 数値差0／0 |
| manifest.txt | 許可された削除設定17行だけ除去、保持行は一致 | 保持内容の差なし |
| VTU100本 | 0000、0001、0003だけバイト一致。残り97本にPoints座標差。Points以外の全文（Cells・点/セル数等）は全ファイル一致 | **1.0658141036401503e-14／1.1100638479121887e-15** |
| 接触集計 | 全9項目が一致（反復5672、平均2.836、最大3、非接触36、収束1964、上限0、f2mode7672、mode2uf15344、calcDis5672） | 集計差0。接触力の内部配列一致は未証明 |
| 通常加振の保持出力 | 中断したため比較未完了 | 未測定 |

VTUの差分はファイル番号順の最初が `deform_combined0002.vtu`、0-based節点39927のx：`3.321402542165793` → `3.321402542165792`。座標成分の差は全100ファイル累計327329件。最大絶対差は `deform_combined0007.vtu` の節点24817・z：`14.50091115681474` → `14.50091115681475`。最大相対差は `deform_combined0079.vtu` の節点25752・z：`1.000143394195989` → `1.00014339419599`。相対差は `abs(after-before)/abs(before)`、基準0の成分に非0となった例は0件。新たな非有限座標はない。座標の単位は出力値の単位に従う。

差は小さいが、絶対/相対許容誤差0の条件を満たさないため許容しない。DATの出力精度で一致しても内部状態が完全一致する証明にはならない。計算式・主な計算ブロックの静的保持だけでは、実行結果の同一性を証明できなかった。コンパイラの丸め・並列性・削除対象の副作用等のどれが原因かは未特定であり、断定しない。許容誤差の緩和、計算式/ビルド条件の変更、独断の復旧や追加切り分けは行わない。

試験paramのSHA256は整理前後で一致：通常 `c9c1fa0aef73391620bbe92f7f819d08af46d32e1126e97456c80b0662e94f84`、正弦波 `abf95ba6ad9fd331cc2a1a1fe57b7b3a7ba3c0d3d51cd8b265c27cf3581d28ec`。同一のビルド設定・実行環境を使用したが、原因を確定する追加実験は未実施。

## Python解析への影響

静的な読込/設定生成の照合結果であり、各解析の実行成功を表すものではない。toolsのコード・過去CSV・派生CSVを削除/変更しない。過去の保存済みCSVに対する解析はその入力を維持する。

| Pythonスクリプト | 依存していたC++ CSV | 代替DAT等の有無 | 変更後の動作可否 | 次回作業候補 |
|---|---|---|---|---|
| analyze_irregularity.py:42 | irregularity_timeseries | displace/area/airflow/pressure_vtに部分情報あり、接触フラグなし | 新runは必須入力欠落で不可 | 明示依頼後にDATベースの指標を設計 |
| analyze_ap_contact.py:145–147,186 | modal_force_split、flow_diagnostics_v2、contact_region_summary/pairs_detail、perturbation_energy | 一般化力・仕事の代替DATなし | 新runのAP/仕事解析は不可 | 必要な研究指標と保存方法を再検討 |
| run_ap_contact_diagnostics.py:78以降 | AP/擾乱/freeze系の設定と上記CSV | 計算機能自体を削除 | 削除設定を渡すため新実行不可 | 残すAP荷重除外と新試験方式を分離して設計 |
| analyze_modal_exit_test.py:701,718 | modal_force_split、perturbation_energy | displaceからモード一般化力/仕事は復元不可 | 新runは不可 | 解析継続の要否を別途判断 |
| analyze_self_oscillation.py:162–194 | irregularity、flow_diagnostics_v2、energy_summary、debug_flow_detail fallback | displace/areaは保持。圧・流量DATはあるが自動fallback未実装 | 変位の部分解析は可能、接触率/流れ/仕事評価は欠落 | 必要ならDAT読込と欠落指標の扱いを明示設計 |
| run_diagnostic_sweep.py | self_oscillationの結果、生成する削除診断設定 | initialGap/Psは保持、擾乱/freeze設定は削除 | 現状の生成paramは削除オプションを含むため不可 | 通常計算に限定したsweepを別途設計 |
| run_pressure_sweep.py:27–33 | modal_contribution/dominant、contact_iteration/final/coupling | 基本DATは保持、CSV保存部は存在チェックあり | 基本の圧力掃引・DAT解析は可能、モード/接触の保存は欠落 | 保存リストと必要な指標を後続依頼で整理 |
| kcExplore.py:12,118 | debug_step_summary | displace/areaあり、診断列の代替なし | 元から旧絶対パス・ヘッダーのみ。新runの診断読込は不可 | パスと探索の要件を別途見直す |
| plot_modal_top10.py:32 | modal_top10 | 代替なし | 新runでは不可 | モード解析の要否を判断 |
| plot_modal_energy_share.py:126 | modal_contribution | q/qdotのDATなし | 新runでは不可 | モード状態の必要性を判断 |
| summarize_modal_force_split.py:148 | modal_force_split | 代替なし | 新runでは不可 | 仕事/力解析の継続可否を判断 |
| analyze_dynamics.py:239、dynamics_common.py:249 | modal_contribution（任意） | area/displace/airflow/pressure_vt/manifestあり | 基本周期・左右比は可能、モード関連だけ省略 | 省略表示が必要なら別途対応 |
| surfp.py:26,59 | surfp_output | CSV自体を保持 | 従来どおり。既存targetX/nxsup既定の不整合は今回変更しない | 必要なら別途合わせる |
| monitor_point.py、plot_pair_offset_summary.py | monitor_point、pair_offset_summary | 今回のC++生成対象ではない | 今回の削除による追加変更なし。元から別入力が必要 | 生成元・用途を別途確認 |
| plot_bifurcation.py、plot_poincare.py | Python派生のbifurcation_points/sweep_summary/poincare_points | 派生CSVは削除対象外 | 元のDAT解析経路が成立する範囲で維持 | 今回変更なし |
| その他DAT/WAV解析、left_right_ratios、compare_stiffness_metrics、spectol_sweep等 | 今回削除するC++ CSVへの必須依存なし | 維持DAT/WAV/params/manifest | 今回の削除による基本入力の欠落なし。既存のパス・時間窓前提は残る | 今回変更なし |

## 未検証・完了状態

**安全停止・未完了。** ソース14ファイルと本レポートの未コミット変更を残す。ビルド/CTestは成功したが、数値回帰は正弦波で失敗、通常は途中中断で未判定。利用者の「全く同じ結果」という完了条件は未達である。追加の原因切り分け・復旧は利用者の指示を待つ。

0.02秒の正弦波を長時間検証とは呼ばない。通常側の基準のみ0.1362秒で診断区間を通過し、整理後の全区間は未検証。本番1秒の長時間自励/非周期挙動、0.25–0.30秒の旧基準窓、声道あり、forceDirection=1、AP除外ON、ncont=0、PROFILE_CALCDIS=ONのシミュレーションは未検証。接触力の内部配列・全一般化状態の厳密一致も未証明。元の本番build/input/output/resultおよび研究データは保持。`git diff --check` は指摘なし、src/includeの `.csv` 参照はNASのsurfpだけである。
