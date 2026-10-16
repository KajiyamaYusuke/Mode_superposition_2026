# Phase 4 最終ディレクトリ整理報告

作業日：2026-10-16 JST。指示書：ルートの `CODEX_phase4_final_directory_reorganization.md`。

## 結果

指定された配置整理を完了した。CPP14個・ヘッダー15個・文書9個の計38ファイルを `git mv` で移動し、`src/main.cpp` は据え置いた。**全15CPP・全15ヘッダーおよび移動文書9個は移動前とバイト完全一致し、移動先の存在を確認した。** C++ファイルの内容、API、物理式、演算、関数呼出し、Phase 4段階Cには手を加えていない。

Releaseで `simulation` と既存3テストターゲットのコンパイル・リンクに成功した。**依頼により数値計算・数値回帰・CTest・テスト実行ファイルの実行は未実施**であり、今回数値的同等性を再判定したものではない。

## 開始状態

- HEAD：`0e701d8721a336c1ac75af5fb7d23f5529a41696`（`docs: add Phase 4 directory reorganization instructions`）。
- 指定基準 `573614721258d1c31ac18d492e11c94093196a77`（`clean_AB`）との差は今回の指示書追加だけだった。
- `git status --short`：`?? tools/testFluid.py` のみ。作業ツリー・インデックスの既存tracked差分は空。
- `tools/testFluid.py` は移動・編集・削除・実行・ステージングをしていない。開始時SHA256は `adf0ca7cfd47f166665825e980873c47cb0eb12810c65ffc43c495837f8c86dd`。
- 開始HEAD・status・diff/cached diff・trackedファイルSHA256・移動対応表・保護対象の記録は `/tmp/vocal-layout-audit-T2LfNp/start.json` に保存した。

## CPPの配置

| 旧パス | 配置後 |
|---|---|
| `src/main.cpp` | `src/main.cpp`（据え置き） |
| `src/simulation.cpp` | `src/simulation/simulation.cpp` |
| `src/readParams.cpp` | `src/simulation/readParams.cpp` |
| `src/geometry.cpp` | `src/structure/geometry.cpp` |
| `src/modedata.cpp` | `src/structure/modedata.cpp` |
| `src/state.cpp` | `src/structure/state.cpp` |
| `src/rungeStep.cpp` | `src/structure/rungeStep.cpp` |
| `src/channelSections.cpp` | `src/fluid/channelSections.cpp` |
| `src/flowModel.cpp` | `src/fluid/flowModel.cpp` |
| `src/separationSelection.cpp` | `src/fluid/separationSelection.cpp` |
| `src/forceCalculator.cpp` | `src/forces/forceCalculator.cpp` |
| `src/forceCalculatorFluid.cpp` | `src/forces/forceCalculatorFluid.cpp` |
| `src/forceCalculatorContact.cpp` | `src/forces/forceCalculatorContact.cpp` |
| `src/modalProjector.cpp` | `src/forces/modalProjector.cpp` |
| `src/surfaceLoad.cpp` | `src/forces/surfaceLoad.cpp` |

`forceCalculatorFluid.cpp` はForceCalculatorの荷重実装としてforcesに配置し、クラスを分けていない。`rungeStep.cpp` の名前も変更していない。

## ヘッダーの配置

| 旧パス | 配置後 |
|---|---|
| `include/Simulation.h` | `include/simulation/Simulation.h` |
| `include/SimulationParams.h` | `include/simulation/SimulationParams.h` |
| `include/wavwrite.h` | `include/simulation/wavwrite.h` |
| `include/Geometry.h` | `include/structure/Geometry.h` |
| `include/ModeData.h` | `include/structure/ModeData.h` |
| `include/State.h` | `include/structure/State.h` |
| `include/Displacement.h` | `include/structure/Displacement.h` |
| `include/TimeIntegrator.h` | `include/structure/TimeIntegrator.h` |
| `include/ChannelSections.h` | `include/fluid/ChannelSections.h` |
| `include/FlowModel.h` | `include/fluid/FlowModel.h` |
| `include/FluidConstants.h` | `include/fluid/FluidConstants.h` |
| `include/SeparationSelection.h` | `include/fluid/SeparationSelection.h` |
| `include/ForceCalculator.h` | `include/forces/ForceCalculator.h` |
| `include/ModalProjector.h` | `include/forces/ModalProjector.h` |
| `include/SurfaceLoad.h` | `include/forces/SurfaceLoad.h` |

ヘッダー名15個に重複はない。内容・メンバー順・クラス宣言・可視性を変更せず、既存の非接頭辞include形式を維持した。

## 文書の配置

| 旧パス | 配置後 |
|---|---|
| `CODEX_cleanup_phase1_2.md` | `docs/cleanup/CODEX_cleanup_phase1_2.md` |
| `CODEX_cleanup_phase3_remove_csv_and_legacy.md` | `docs/cleanup/CODEX_cleanup_phase3_remove_csv_and_legacy.md` |
| `CODEX_cleanup_phase4_forcecalculator.md` | `docs/cleanup/CODEX_cleanup_phase4_forcecalculator.md` |
| `cleanup_inventory.md` | `docs/cleanup/cleanup_inventory.md` |
| `cleanup_change_report.md` | `docs/cleanup/cleanup_change_report.md` |
| `cleanup_phase3_report.md` | `docs/cleanup/cleanup_phase3_report.md` |
| `cleanup_phase4_forcecalculator_report.md` | `docs/cleanup/cleanup_phase4_forcecalculator_report.md` |
| `MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md` | `docs/plans/MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md` |
| `MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md:Zone.Identifier` | `docs/plans/MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md:Zone.Identifier` |

特殊名の `:Zone.Identifier` はtrackedファイルで、Linuxの現在環境では同名のまま問題なく `git mv` できた。内容を推測して変更したり削除したりしていない。過去の文書の行番号・ソースパス・検証判定は当時の記録としてそのまま維持した。

本報告書のみ新規作成した。今回の指示書自体は移動対象に指定されていないためルートに残した。それ以外の指定外ファイルも移動・削除していない。

## CMake・構成説明・参照確認

`CMakeLists.txt` の変更は次に限定した。

- simulationターゲットの15CPPを指示書の新パスと列挙順へ更新。全CPPがちょうど1回登録され、存在することを機械的に照合した。GLOBは導入していない。
- `SIMULATION_INCLUDE_DIRS` にinclude直下と4分類の検索パスを定義し、simulationと既存3テストターゲットで共用した。
- test_separation_selectionのseparationSelection.cpp、test_mode_selectionのmodedata.cpp、test_input_file_paramsのreadParams.cppの登録パスを更新した。
- icpx、C++17、Releaseの実効 `-O3`、OpenMP、PROFILE_CALCDISの定義・条件、include(CTest)、BUILD_TESTING条件、add_testと引数は維持した。input_file_paramsの `${CMAKE_SOURCE_DIR}/input/param.txt` も変更していない。

AGENTS.mdはディレクトリ構成説明だけを更新した。入口、4分類、docs、testsの説明を反映し、ビルド/実行コマンド例・パラメータ仕様・その他のガイドラインはそのまま残した。

移動前後にrgで旧CPP/Hパスと移動文書名の参照を検索した。`.gitignore` の `simulation` 規則が同名の新ディレクトリにも適用されるため、移動後の検索は `--no-ignore` を用い、src/simulationとinclude/simulationも含めて確認した。CMake以外に変更が必要な有効な固定パス参照は見つからず、Pythonツール・テストソース・include文は変更していない。過去のcleanup文書は検索対象から除外して履歴を保持し、今回の指示書の旧→新対応表も書き換えていない。

AGENTS.mdの `nonperiodic_vibration_analysis_plan.md` は実在しない既存参照であるため、その旨だけを構成説明に明記した。類似の計画書への読み替え、新規作成はしていない。未確認資料の既存記述として報告し、現行のビルド・実行コードに解決不能な旧パス参照は確認していない。移動保留ファイルはない。

## 内容一致・配置確認

開始HEADを固定し、全30CPP/Hについて、旧パスの `git show 0e701d8721a336c1ac75af5fb7d23f5529a41696:<旧パス>` の生バイトと配置後ファイルの `read_bytes()` を機械的に比較した。全件一致した。移動文書9個も同じ方法で全件一致した。

すべての移動先が通常ファイルとして存在し、移動元には残っていない。src直下のCPPはmain.cppだけ、include直下のヘッダーは0個であり、再帰的な総数はCPP15・ヘッダー15のまま。報告書も指定先に配置した。確認スクリプトと対応別SHA256・結果は `/tmp/vocal-layout-audit-T2LfNp/audit.py`、同先 `verification.json` に保存した。

## Releaseビルド・リンク

一時ビルド先は `/tmp/vocal-layout-build-IRIFGr`。既存build/は使用・変更していない。

```bash
cmake -S . -B /tmp/vocal-layout-build-IRIFGr -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build /tmp/vocal-layout-build-IRIFGr --target simulation test_separation_selection test_mode_selection test_input_file_params -j 4
```

構成生成・ビルドとも終了0。IntelLLVM/icpx 2025.1.0とOpenMPを使用した。simulationの生成されたflags.makeでも、実効フラグは `-O3 -std=gnu++17 -fiopenmp`、インクルード検索パスは指定5ディレクトリであることを確認した。PROFILE_CALCDISは既定OFFのまま。

| ターゲット | コンパイル・リンク | 一時ビルド先での実行ファイル存在 |
|---|---|---|
| simulation | 成功 | 確認済み |
| test_separation_selection | 成功 | 確認済み |
| test_mode_selection | 成功 | 確認済み |
| test_input_file_params | 成功 | 確認済み |

実行権限の存在も確認したが、いずれの実行ファイルも起動していない。CTest、通常/正弦波計算、数値比較、PROFILE_CALCDIS=ONの追加ビルド・プロファイル調査は依頼により未実施。

## データ保護・Git状態

開始時と配置後・ビルド後の保護対象を機械的に照合した。本番input/はSHA256・サイズ・mtime・ファイル種別、output/・result/・analysis_runs/・diagnostic_runs/・既存build/は全エントリのサイズ・mtime・ファイル種別・存在一覧が一致した。tools/testFluid.pyもSHA256・サイズ・mtimeが一致した。既存データ、ビルド成果物、計算出力は消去・移動・上書きしていない。

Gitの変更は38件のR100リネーム、CMakeLists.txtとAGENTS.mdの更新、本報告書の新規作成だけで、別途利用者のtools/testFluid.pyが未追跡のまま残る。git mvに伴う38件のリネームはステージ済みだが、CMake/AGENTSの編集と本報告書は追加ステージングしていない。tools/testFluid.pyもステージングしていない。

`git diff --check` と `git diff --cached --check` は指摘なし。HEADを変更せず、commit/push、reset/restore/stash/clean、広範囲削除は行っていない。指定の配置整理・存在確認・内容一致・Releaseビルド/リンク確認までで作業を終了し、次の機能変更や段階Cには進まない。
