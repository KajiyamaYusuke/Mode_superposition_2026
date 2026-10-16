# Codex作業指示書：Phase 4 最終作業 — ディレクトリ構成の整理

> **作業目的：ファイル配置を整理してPhase 4を一区切りにする。**  
> **許可された変更：ファイルの移動、CMakeのパス更新、必要最小限のインクルードパス調整、AGENTS.mdの構成説明更新、文書の配置整理。**  
> **許可されていない変更：物理式・計算処理・データ構造・メンバー・関数呼び出し・アルゴリズム・数値出力の変更。**  
> **数値計算の実行・数値回帰テストは不要。Releaseビルド／リンク確認のみ行う。**

## 0. 前提と作業方針

この指示書は、声帯振動計算コード `Mode_superposition_2026` のディレクトリ構成を整理するためのもの。作業の開始基準は、Phase 4の段階A・Bの完了コミット：

- リポジトリ：`KajiyamaYusuke/Mode_superposition_2026`
- 基準コミット：`573614721258d1c31ac18d492e11c94093196a77`
- コミットメッセージ：`clean_AB`
- この指示書の作成時点：2026-10-16

**実際に編集するのはCodexが開いているローカル作業ツリーである。** まず `AGENTS.md`、`git status`、現在のHEADを確認し、ローカルの状態を優先する。基準より新しい変更が存在する場合、その変更を上書き・破棄せず内容を確認してから作業する。

目的は**コードの見通しを改善することだけ**であり、Phase 4の段階C（`applyContactLoads()` の内部関数分割）は実施しない。`ForceCalculator` のメンバー、API、コンストラクタ、ヘッダー内容、関数の呼び出し順、ファイル中の演算は一切変更しない。

今回の作業は独立したフォルダー分類であり、**数値的同等性の再判定やコンパイラ最適化の調査を始めないこと**。全CPP/Hの内容を変えず、ビルドが成功すれば作業を終了してよい。

### 作業開始前の確認

```bash
git status --short
git rev-parse HEAD
git diff --name-status
git diff --cached --name-status
```

- 既存の未コミット差分、未追跡ファイルを一覧化する。
- 特に、ユーザーが保持している未追跡の `tools/testFluid.py` を**絶対に削除・移動・編集・ステージングしない**。
- 既存の `input/`、`output/`、`result/`、`analysis_runs/`、`diagnostic_runs/`、`build/` には手を触れない。
- `git reset --hard`、`git clean`、`git restore`、`git stash`、広範囲の一括削除を使用しない。
- 本書の移動対象と競合する既存変更がある場合、勝手に上書きせず該当ファイルを保留して理由を報告する。

---

## 1. 完成後のディレクトリ構成

`src/` と `include/` をそれぞれ同じ4領域に分ける。`main.cpp` のみ `src/` 直下に置く。

```text
Mode_superposition_2026/
├── CMakeLists.txt
├── AGENTS.md
├── .gitignore
├── .python-version
├── pyproject.toml
├── uv.lock
├── requirements-analysis.txt
│
├── src/
│   ├── main.cpp
│   ├── simulation/
│   │   ├── simulation.cpp
│   │   └── readParams.cpp
│   ├── structure/
│   │   ├── geometry.cpp
│   │   ├── modedata.cpp
│   │   ├── state.cpp
│   │   └── rungeStep.cpp
│   ├── fluid/
│   │   ├── channelSections.cpp
│   │   ├── flowModel.cpp
│   │   └── separationSelection.cpp
│   └── forces/
│       ├── forceCalculator.cpp
│       ├── forceCalculatorFluid.cpp
│       ├── forceCalculatorContact.cpp
│       ├── modalProjector.cpp
│       └── surfaceLoad.cpp
│
├── include/
│   ├── simulation/
│   │   ├── Simulation.h
│   │   ├── SimulationParams.h
│   │   └── wavwrite.h
│   ├── structure/
│   │   ├── Geometry.h
│   │   ├── ModeData.h
│   │   ├── State.h
│   │   ├── Displacement.h
│   │   └── TimeIntegrator.h
│   ├── fluid/
│   │   ├── ChannelSections.h
│   │   ├── FlowModel.h
│   │   ├── FluidConstants.h
│   │   └── SeparationSelection.h
│   └── forces/
│       ├── ForceCalculator.h
│       ├── ModalProjector.h
│       └── SurfaceLoad.h
│
├── docs/
│   ├── cleanup/
│   │   ├── CODEX_cleanup_phase1_2.md
│   │   ├── CODEX_cleanup_phase3_remove_csv_and_legacy.md
│   │   ├── CODEX_cleanup_phase4_forcecalculator.md
│   │   ├── cleanup_inventory.md
│   │   ├── cleanup_change_report.md
│   │   ├── cleanup_phase3_report.md
│   │   ├── cleanup_phase4_forcecalculator_report.md
│   │   └── directory_reorganization_report.md  # 今回作成
│   └── plans/
│       └── MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md
│
├── tests/
│   ├── test_separation_selection.cpp
│   ├── test_mode_selection.cpp
│   └── test_input_file_params.cpp
├── tools/                 # Python解析ツール：今回一括移動しない
└── input/                 # 研究入力：変更しない
```

上記に記載した以外の既存ファイルは**勝手に削除・移動しない**。`output/`、`result/`、`build/` 等は環境によって存在するが、生成物としてそのまま残す。

### 分類ルール

| フォルダー | 役割 |
|---|---|
| `simulation/` | 計算全体の実行制御と入力パラメータ |
| `structure/` | 声帯の形状、モードデータ、構造状態、時間積分 |
| `fluid/` | 流路断面、1D流れ、剥離点の選択 |
| `forces/` | 流体力、接触力、表面荷重、モード力への射影 |
| `docs/cleanup/` | 過去のコード整理の指示書・検証報告と今回の作業報告 |
| `docs/plans/` | 研究・解析の設計資料、実装計画 |

`forceCalculatorFluid.cpp` は名前にFluidが入るが、`ForceCalculator` クラスの**荷重生成実装**であるため `forces/` に配置する。これを `fluid/` へ分離してクラスやAPIを変えることはしない。

`rungeStep.cpp` は現在Newmark法を含む古い名前だが、**今回はリネームしない**。`timeIntegrator.cpp` への改名は必要なら将来の独立作業とする。

---

## 2. 移動するC++ソースファイル（15ファイル）

**`git mv` で移動し、ファイル内容を編集しないこと。** 以下のうち `src/main.cpp` は移動しない。

| 旧パス | 新パス |
|---|---|
| `src/main.cpp` | **変更なし** |
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

## 3. 移動するヘッダーファイル（15ファイル）

| 旧パス | 新パス |
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

**ヘッダーの中身・名前・クラス宣言・可視性・メンバー順序は変更しない。** 既存の非接頭辞形式 `#include "ForceCalculator.h"` や `#include "Geometry.h"` は、後述するCMakeの検索パス設定で解決し、原則書き換えない。

---

## 4. CMakeLists.txt の修正（今回の中心作業）

現行の実行ファイル名 `simulation`、`icpx`、C++17、Release `-O3`、OpenMP、`PROFILE_CALCDIS` オプション、およびテストターゲットの名前と有効/無効条件を保持する。

### 4.1 ソースパスを新しいパスへ変更

`add_executable(simulation ...)` のソース列挙を以下に対応するパスへ更新する。

```cmake
add_executable(simulation
    src/main.cpp
    src/simulation/simulation.cpp
    src/simulation/readParams.cpp
    src/structure/geometry.cpp
    src/structure/modedata.cpp
    src/structure/state.cpp
    src/structure/rungeStep.cpp
    src/fluid/channelSections.cpp
    src/fluid/flowModel.cpp
    src/fluid/separationSelection.cpp
    src/forces/forceCalculator.cpp
    src/forces/forceCalculatorFluid.cpp
    src/forces/forceCalculatorContact.cpp
    src/forces/modalProjector.cpp
    src/forces/surfaceLoad.cpp
)
```

すべての既存CPPが**ちょうど1回**含まれるようにする。存在しないCPPの登録、ファイルの重複、`GLOB`での無条件自動探索への変更は行わない。

### 4.2 ヘッダーの検索パスを設定

今回は **C++全体の `#include` を書き換えるのではなく、CMakeのインクルードディレクトリで既存コードを維持**する。以下をCMakeに定義し、実行ファイル・テストターゲット双方に使う。

```cmake
set(SIMULATION_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/include/structure
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fluid
    ${CMAKE_CURRENT_SOURCE_DIR}/include/forces
)

target_include_directories(simulation PRIVATE ${SIMULATION_INCLUDE_DIRS})
```

既存ヘッダーはファイル名が重複していないことを確認すること。`#include "Geometry.h"` 等が従来のまま解決するようにする。全ソースのincludeを `#include "structure/Geometry.h"` に一括変換する必要はない。

### 4.3 テスト用ビルドターゲットのソースパスを修正

実行している既存のCTestは3ターゲット。**今回はテストを実行しないが、ビルド対象から消さないこと。**

```cmake
add_executable(test_separation_selection
    tests/test_separation_selection.cpp
    src/fluid/separationSelection.cpp
)
target_include_directories(test_separation_selection PRIVATE ${SIMULATION_INCLUDE_DIRS})

add_executable(test_mode_selection
    tests/test_mode_selection.cpp
    src/structure/modedata.cpp
)
target_include_directories(test_mode_selection PRIVATE ${SIMULATION_INCLUDE_DIRS})

add_executable(test_input_file_params
    tests/test_input_file_params.cpp
    src/simulation/readParams.cpp
)
target_include_directories(test_input_file_params PRIVATE ${SIMULATION_INCLUDE_DIRS})
```

※ これはCMakeの該当箇所の参考抜粋。実装時には既存の `include(CTest)`、`if(BUILD_TESTING)`、`add_test(...)` とその引数を**そのまま保持し、重複して追加しない**。 `input_file_params` の `${CMAKE_SOURCE_DIR}/input/param.txt` もそのまま保持する。

### 4.4 コンパイラ設定は触らない

- `set(CMAKE_CXX_COMPILER icpx)` を変更しない。
- Release最適化 `-O3`、OpenMP、コンパイル定義、リンクフラグを変更しない。
- 新しいライブラリ／依存パッケージ／実行形式は導入しない。
- `CMakeLists.txt` 全体の不要な整形や大規模な書換えは行わず、パス変更に集中する。

---

## 5. docs/ の整理

### 5.1 過去の整理文書を `docs/cleanup/` に移動

次のファイルを内容を変更せずに `git mv` する。

```text
CODEX_cleanup_phase1_2.md
CODEX_cleanup_phase3_remove_csv_and_legacy.md
CODEX_cleanup_phase4_forcecalculator.md
cleanup_inventory.md
cleanup_change_report.md
cleanup_phase3_report.md
cleanup_phase4_forcecalculator_report.md
```

これらは各Phaseの**過去時点の検証記録**であるため、記載されている古い行番号やソースパスを現在の構成に合わせて機械的に修正しない。歴史的記録の意味を守る。

### 5.2 研究・実装計画を `docs/plans/` に移動

```text
MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md
```

さらに、ルートには次の特殊名のファイルがGitHub上で確認されている。

```text
MINIMAL_MODAL_WORK_EXIT_TEST_SPEC.md:Zone.Identifier
```

これは通常の設計資料とは異なるメタデータ風の名前を持つ。**勝手に削除したり、内容を推測して変更しない。** 現在の作業環境でtrackedであり、問題なく `git mv` できるなら `docs/plans/` に同名で移してよい。Windows側のファイル名制約などで支障がある場合は手を触れず、その理由を報告する。削除が必要なら別途許可を得る。

存在しない `nonperiodic_vibration_analysis_plan.md` などの文書を、過去の記述を根拠に新規作成しない。

### 5.3 新規報告書

今回の変更報告書を `docs/cleanup/directory_reorganization_report.md` として新規作成。内容は第9節に従う。

---

## 6. AGENTS.md とパス参照の更新

`AGENTS.md` はルートに残し、**現在のディレクトリ構成の説明だけ**更新する。最低限、以下を反映する。

- `src/main.cpp` は入口、`src/simulation/` は実行制御・設定。
- `src/structure/` は形状・モード・構造状態・積分。
- `src/fluid/` は流路断面、流れ、分離点。
- `src/forces/` は流体力、接触力、モード力射影。
- `include/` は同名4分類でヘッダーを管理。
- `docs/cleanup/` は整理報告、`docs/plans/` は仕様・計画。
- `tools/`、`tests/`、`input/` の機能と場所は変わらない。
- `rungeStep.cpp` の名前に反して現在Newmark法が含まれることは必要なら補足できる。ただしコードの改名は禁止。

コンパイラコマンド例・実行例・研究パラメータの仕様は今回変更しない。現在実在しない資料への参照があれば**勝手に類似文書へ読み替えず**、「未確認／既存記述」として報告する。

### 6.1 全文検索による影響範囲確認

最低限、次を検索する。

```bash
rg -n 'src/(simulation|readParams|geometry|modedata|state|rungeStep|channelSections|flowModel|separationSelection|forceCalculator|forceCalculatorFluid|forceCalculatorContact|modalProjector|surfaceLoad)\.cpp' --glob '!docs/cleanup/**' .
rg -n 'include/(Simulation|SimulationParams|wavwrite|Geometry|ModeData|State|Displacement|TimeIntegrator|ChannelSections|FlowModel|FluidConstants|SeparationSelection|ForceCalculator|ModalProjector|SurfaceLoad)\.h' --glob '!docs/cleanup/**' .
rg -n 'CODEX_cleanup_|cleanup_(inventory|change_report|phase3_report|phase4_forcecalculator_report)|MINIMAL_MODAL_WORK_EXIT_TEST_SPEC' --glob '!docs/cleanup/**' .
```

`rg` が無い場合は `grep -R` 等で代用する。`.git/`、`build/`、計算結果ディレクトリ、巨大な入力VTUは検索対象から除外してよい。

- CMake、AGENTS、**現在の作業に使用するリンク／スクリプト**が旧パスを参照する場合は、必要な文字列パスのみ修正してよい。
- Python `tools/` にC++のソースファイルを直接開くコードがあれば、その固定パスだけ最小修正してよい。ただし解析アルゴリズム、依存ライブラリ、出力内容、入力ファイル形式は変更しない。
- 過去の調査・検証レポート中の古いソースパス、行番号は当時の記録なので修正しない。
- 歴史的記録以外に解決不能な参照が残った場合、隠さず報告する。

---

## 7. 実施手順（推奨順序）

1. 作業ツリー状態とHEADを記録し、対象ファイルの有無を確認する。
2. `mkdir -p` で `src/{simulation,structure,fluid,forces}`、`include/{simulation,structure,fluid,forces}`、`docs/{cleanup,plans}` を作る。
3. 第2〜3節の表に従って **`git mv` でCPP/Hを移動**する。コピーしてから元を削除するより、履歴の追跡を容易にする。
4. 第5節の文書を移動する（特殊 `:Zone.Identifier` は上記条件に従う）。
5. `CMakeLists.txt` を第4節のとおり変更する。
6. 既存の `#include` のままヘッダーが解決することを優先する。C++の内容を変えることが必要になった場合、まずCMakeの検索パス不足や同名衝突を確認する。
7. `AGENTS.md` と**現在有効な**参照先だけを最小修正する。
8. 第8節の確認を行う。成功したら第9節の報告書を作成する。
9. **自動コミット・pushせず**、変更点と結果をユーザーへ報告して停止する。

---

## 8. 確認方法：数値テストなし、ビルドのみ

### 8.1 移動が完全か検査

以下を確認する。

- もとの `src/` 直下には `main.cpp` 以外のCPPが残っていない。
- もとの `include/` 直下に移動忘れのヘッダーが残っていない。
- `src/` のCPP総数は元と同じ15（`main.cpp` 含む）。
- `include/` のヘッダー総数は元と同じ15。
- ファイル名の重複・ソースのCMake登録漏れ／重複がない。
- **C++ソースとヘッダーの内容が移動前とバイト単位で一致**する。

移動前と現行の比較には、例えば次のような手順を用いる。

```bash
# 例：移動前に追跡されていたコードの内容と移動後の内容を比較
git show HEAD:src/forceCalculatorContact.cpp | cmp - src/forces/forceCalculatorContact.cpp
git show HEAD:include/ForceCalculator.h | cmp - include/forces/ForceCalculator.h

# 変更内容の監査
git diff --check
git diff --cached --check
git diff --name-status -M
git diff --cached --name-status -M
```

上記の `cmp` は例示。**30個すべてのCPP/Hについて旧パス→新パスの対応を使って機械的に比較**し、一件でも不一致なら確認して停止。`git mv` したファイルは通常インデックスにステージされるため、`git diff --cached` も確認する。

### 8.2 Releaseビルド／リンク検査

既存の `build/` に影響しないよう、別の一時ビルドディレクトリを使用する。ビルドコマンド例：

```bash
build_dir=$(mktemp -d /tmp/vocal-layout-build-XXXXXX)
cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "$build_dir" --target simulation test_separation_selection test_mode_selection test_input_file_params -j
```

- Intel `icpx` およびOpenMPは元のCMake設定のまま使用する。
- `simulation` と既存の3テスト実行ファイルが**コンパイル・リンクできること**を確認する。
- **`ctest`、`./simulation`、その他の数値計算の実行はしない。** ユーザーが数値回帰テストを不要と明示しているため。
- `PROFILE_CALCDIS=ON` の別ビルドやプロファイリング実行も今回の必須条件にしない。コンパイルフラグを変更する必要はない。
- ビルドに失敗した場合は、パス・include検索の問題に限って修正する。物理計算コードを書き換えてコンパイルエラーを隠してはならない。
- `icpx` 自体が使えない環境なら、理由を報告してビルド未検証とする。ユーザーの許可なしにCMakeのコンパイラ指定を変えない。

### 8.3 データ保護検査

- 本番 `input/` のパラメータファイル・VTU/NAS/周波数ファイルは変更しない。
- `output/`、`result/`、`analysis_runs/`、`diagnostic_runs/`、既存のビルド成果物は消去・移動・上書きしない。
- `tools/testFluid.py` は未追跡で保持し、`git add -A` 等で意図せずステージングしない。
- `git status --short` に現れるものは意図したファイル移動・CMake・AGENTS・必要最小限の参照修正・新規報告書だけか確認する。

---

## 9. 完了報告書の作成

`docs/cleanup/directory_reorganization_report.md` を新規作成し、最低限次を報告する。

1. 作業開始時の `HEAD` と `git status`、既存未追跡ファイルの保護状況。
2. 移動したCPP/Hの件数、具体的な旧パス→新パス一覧。`main.cpp` は据え置き。
3. 移動した文書、`:Zone.Identifier` の扱い。
4. CMakeで変更したソース登録・インクルード検索パス・テストターゲットのソースパス。
5. AGENTSやPython等で、実際に変更した参照パス。
6. **全CPP/Hが元とバイト一致したか**、確認方法と結果。もし不一致なら理由と作業停止状況。
7. Releaseビルド・リンクの成否（`simulation` と既存3テストターゲット）。数値テストは**依頼により未実施**と明記。
8. 本番入力・生成出力・未追跡 `tools/testFluid.py` の保護確認。
9. 未解決の参照、移動保留ファイル、任意の今後の改善案（必要な場合のみ）。

## 10. 完了判定と禁止事項の再確認

**完成条件：**

- 第1節のディレクトリ構造へ移動できている（例外は報告済み特殊ファイルのみ）。
- 全15CPP・全15ヘッダーの内容は移動前と完全一致する。
- CMakeのパスがすべて有効で、Releaseのビルド・リンクが成功する。
- 既存の研究データや計算結果は変更されていない。
- 変更内容を `docs/cleanup/directory_reorganization_report.md` に記録した。

**以下は今回行わない：**

- `applyContactLoads()` の内部分割、`ContactSolver` 等の新クラス導入。
- クラスメンバー・可視性・状態の所有者・関数引数・呼出順の変更。
- ソースファイル名の変更（`rungeStep.cpp` → `timeIntegrator.cpp` 等）。
- 物理式・数値モデル・接触判定・減衰・Newmark積分・出力形式の変更。
- Python解析コードの全面移動や再設計。
- 数値計算実行、CTest実行、VTU/DAT/WAVの再生成や数値比較。
- 自動commit/push、ブランチ履歴の書換え、無関係な削除。

**このディレクトリ整理が終了したら、Phase 4の作業はここで一区切りとし、次の機能変更へ勝手に進まないこと。**
