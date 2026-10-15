# Phase 1・2 変更・検証報告

作業日: 2026-10-14（JST）。棚卸し・参照根拠は `cleanup_inventory.md`。Phase 1と許可条件を満たすPhase 2の限定整理を完了した。判断が分かれる対象は保持した。

## 作業基準

Git HEAD: `3ecd141786d2ee3384276853f48a84a7c14ad38b`。比較基準はHEADではなく、開始時の未コミット変更を含む現在コードのコピー。

既存変更は以下の通りで、整理の差分とは分けて保護する。

```diff
 input/param.txt
-20000
+50000
-0.03  (zetaL)
+0.02
-0.03  (zetaR)
+0.02
 src/simulation.cpp （n%20==0 && nstep-n<=5000 のブロック内）
-//writeVTKCombined(num, geomL, stateL, geomR, stateR, vtuResultDir.string(), 20);
+writeVTKCombined(num, geomL, stateL, geomR, stateR, vtuResultDir.string(), 20);
```

開始時の完全な `git diff --binary`、status、ソース/入力SHA256、コンパイラ情報、既存出力のハッシュは `/tmp/simulation-cleanup-1cGHFz/start.json`。編集前のソース・ヘッダー・CMakeは同ディレクトリの `original/` に保存。既存変更の退避・破棄のためのgit操作は行わない。

## 対象ごとの判断

この表の「編集前の場所」は開始時コピーの行番号。「バイト一致」は削除対象の非生成を確認したうえで、全維持ファイルを比較した結果。

| 項目 | 編集前の場所 | 対応（削除/維持/保留） | 根拠 | ビルド | iforce=0 比較 | iforce=1 比較 | 留意点 |
|---|---|---|---|---|---|---|---|
| debug_step_summary.csv | simulation.cpp:420 | 保留 | kcExplore.py:12,118に読込参照 | 成功 | 維持ファイル一致 | 維持ファイル一致 | ヘッダーのみだが読み手が存在 |
| debug_mode_summary.csv・maxAbsIndex・maxFi等 | simulation.cpp:431,890,1252 | 削除 | 許可リスト内、専用出力への参照のみ | 成功 | 中間・最終バイト一致 | 中間・最終バイト一致 | 計算のfiL/Rを読み出すだけ |
| maxAbsNodeDisp/minMaxIndex/maxAbsPsurf | simulation.cpp:853,869,909 | 削除 | 許可リスト内の未使用ラムダ、呼び出しなし | 成功 | 中間・最終バイト一致 | 中間・最終バイト一致 | 実行処理なし |
| maxAbsVector | simulation.cpp:814,1336 | 維持 | contact_coupling_debug.csvで使用 | 成功 | 維持出力一致 | 維持出力一致 | 未使用ではない |
| debug_force.txt/debugForceFile | forceCalculator.cpp:112,324、ForceCalculator.h:131 | 削除 | 許可リスト内、ログ生成だけ、読み手なし | 成功 | 最終バイト一致 | 最終バイト一致 | 流体ログのブロックのみ除去 |
| debug_flow_detail.csv/flowDebugFile | forceCalculator.cpp:184,379,559 | 保留 | analyze_self_oscillation.py:191のfallback読込 | 成功 | 維持ファイル一致 | 維持ファイル一致 | 異常圧力の記録・クランプ統計も維持 |
| contact_debug.csv/debugContact区間統計 | forceCalculator.cpp:117,672,1356,1583 | 保留 | 0.120–0.135秒の経路は短時間試験の範囲外 | 成功 | ヘッダー一致、区間未検証 | ヘッダー一致、区間未検証 | 安全性を実証できない削除を保留する指示に従う |
| contact_monitor.csv/contactMonitorFile | forceCalculator.cpp:129,1677、ForceCalculator.h:133 | 削除 | 許可リスト内、監視値の読出しのみ、読み手なし | 成功 | 最終バイト一致 | 最終バイト一致 | 配列・監視設定API・統計カウンタは保持 |
| contact_iteration_debug.csv/contactIterationFile | forceCalculator.cpp:139,1700 | 保留 | run_pressure_sweep.py:32が保存、READMEでも参照 | 成功 | 維持ファイル一致 | 維持ファイル一致 | 接触力の比較出力としても利用 |
| contact_search_debug.csv/contactSearchFile | forceCalculator.cpp:153,1717、ForceCalculator.h:135 | 削除 | 許可リスト内、統計のCSV書込だけ、読み手なし | 成功 | 最終バイト一致 | 最終バイト一致 | 探索・cache・fallback・カウンタと参照引数を保持 |
| contactForceL/R_ij、setContactMonitor、monitorI/J | forceCalculator.cpp:108,213,1068、ForceCalculator.h | 維持 | 接触力配分・他の保持出力と共用、削除禁止 | 成功 | 接触・維持出力一致 | 接触・維持出力一致 | 力の加算インターフェースを変更しない |
| モード/エネルギー/非周期/流路/最終接触出力 | simulation.cpp:351以降、forceCalculator.cpp:160以降 | 維持 | 許可リスト外または研究上の保持指定 | 成功 | 生成された全ファイル一致 | 生成された全ファイル一致 | modal_*、irregularity、gap、separation等。条件付き経路は後述 |
| VTU・正弦波・WAV・収束・物理計算 | simulation.cpp:1623,1663、forceCalculator.cpp:244等 | 維持 | 絶対に残す機能 | 成功 | 維持出力一致 | 維持出力一致 | 有効なVTU呼び出しの既存差分を保持 |
| PROFILE_CALCDIS、Geometry診断、rungeStep等 | CMake/forceCalculator/geometry/rungeStep | 保留 | 今回の削除許可なし | 未編集 | 対象外 | 対象外 | 構造変更や整理を実施しない |
| CTest期待値・SimulationParams既定値 | tests/test_input_file_params.cpp:34–35、SimulationParams.h:79–80 | 維持 | 既存不整合の維持をユーザーが指示 | 基準ビルド成功 | CTestで別途比較 | CTestで別途比較 | テスト無効化・期待値修正・実装値変更なし |

## ビルドとCTest

隔離作業ディレクトリ `/tmp/simulation-cleanup-1cGHFz/work/build` で実行。

```bash
cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 4
ctest --output-on-failure
```

Intel icpx/IntelLLVM 2025.1.0、CMake 3.26.4。フラグは `-O3 -std=gnu++17 -fiopenmp`。OpenMP検出成功、PROFILE_CALCDIS=OFF、BUILD_TESTING=ON。既存CMakeは変更しない。

| 段階 | Releaseビルド | separation_selection | mode_selection | input_file_params | CTest終了コード |
|---|---|---|---|---|---|
| 整理前baseline | 成功 | 成功 | 成功 | 既存失敗 | 8 |
| simulation.cpp整理後 | 成功 | 成功 | 成功 | 同原因の既存失敗 | 8 |
| debug_force除去後 | 成功 | 成功 | 成功 | 同原因の既存失敗 | 8 |
| contact_monitor除去後 | 成功 | 成功 | 成功 | 同原因の既存失敗 | 8 |
| contact_search除去後（最終） | 成功 | 成功 | 成功 | 同原因の既存失敗 | 8 |

既存失敗の原因: `tests/test_input_file_params.cpp:34–35` は診断基準窓0.15–0.19秒を期待する一方、`include/SimulationParams.h:79–80` は0.25–0.30秒を設定する。失敗メッセージは `AP contact diagnostic defaults were not parsed correctly`。ユーザーの許可に従い、この1件が同じ原因で失敗する場合だけ許容する。**CTest全体は成功していない。**

## 基準計算と比較条件

全て同じ隔離ソース/ビルドパス・入力パス・実行cwdを使う。段階間で試験結果だけを `/tmp/simulation-cleanup-1cGHFz/{baseline,stage-simulation,final}/{fluid,sine}/` に保存し、次の実行には空の出力領域を使う。元のoutput/resultは使用しない。

```bash
# cwd: /tmp/simulation-cleanup-1cGHFz/work/build
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores ./simulation ../input/fluid.param
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores ./simulation ../input/sine.param
```

試験paramは本番のコピーで、nstepのみ2000へ短縮。sine側は後段の実効iforceだけ1にする。旧iforceスロットは0のまま、読込互換性を保持する。dt=1e-5、nmode=100、nwrite=20、zeta左右0.02、ncont=100、kc=(5,6000,5000)、Ps=2900Pa、初期間隙0mm、forcef=100Hz、famp=-1e-5、forceDirection=0、声道なし、live荷重、擾乱無効。シミュレーション時間0.02秒。

| 入力 | SHA256 |
|---|---|
| 本番input/param.txt（保護対象） | `8b2dab251ecbb0316a2137f03673de1f045ff5fc8aa1746be7552c58b35548f7` |
| 試験fluid.param | `9316655edebdec158482022caf2aa87c1cc5a353faab0826b8fd0445431b11ab` |
| 試験sine.param | `abf95ba6ad9fd331cc2a1a1fe57b7b3a7ba3c0d3d51cd8b265c27cf3581d28ec` |
| 左右モードVTU | `c91a8ec36c367e56a8909881634a01fdd8ec55ab84d8e915b5bb866852ede16f` |
| 左右表面NAS | `481801e62372983bb4aabba710a22739b7c0db93e7b66b8301d0ff46645b9f92` |
| 左右周波数txt | `17252caffd3dd493f6b37eb1d0ea3c4280b351c3de443e8a8a701a36a546dbd0` |

| 条件 | 基準終了コード | 接触反復数 | 接触収束step | 非接触step | VTU生成数 |
|---|---|---|---|---|---|
| iforce=0 | 0 | 2802 | 802 | 1198 | 100 |
| iforce=1 | 0 | 5672 | 1964 | 36 | 100 |

基準は両モードで `FIRST NORMAL CONTACT DETECTED` を出力。最大反復は流体2、正弦波3、反復上限に到達したstepは両方0。計算所要時間は基準流体121.02秒・正弦波149.01秒だが、VTU書出し・初期化も含む単発計測であり、性能改善の結論には使わない。

最終版も両モードの終了コード0、VTU各100ファイル、上表と同じ反復・接触集計。最終実行時間は流体127.38秒、正弦波169.63秒。通常加振では最大絶対流量4.66706e-4 m³/s、最大絶対断面圧791.097 Pa、観測点の最大絶対uyは左右とも0.550662 mmで、非ゼロの流量・圧力・変位と接触を通過した（出力の保存精度で記載）。

## 数値・入出力検証

比較方針は全維持ファイルのバイナリ比較（Python標準 `filecmp.cmp(..., shallow=False)`）。同時に各ファイルのSHA256も保存する。削除指定ファイルだけがなくなり、他のファイルの追加・欠損・内容差がないことを確認する。これにより数値だけでなくヘッダー、列、行数、時刻、出力間隔、WAVバイト列、全VTU座標とセル情報も比較する。stdoutは所要時間の行を除いて一致を要求する。

- DAT9種類（area/displace/displace_xy/pressure/pressure_vt/airflow_vt/section_x/gap_cubed/separation）の行数・列数・時刻/step列を検査。実際の基本出力間隔はnwriteではなくn%5==0。試験では各400行。
- 声門圧は `pressure.dat` と `flow_diagnostics_v2.csv`、接触力は保持する `contact_iteration_debug.csv`、`contact_final_state.csv`、`modal_force_split.csv` で比較する。
- VTUは全100ファイルのXMLを基準側で読み、片側45333点を合わせた90666点と496536セル、全座標の有限性、connectivityの範囲、offsets末尾とconnectivity長、types数、全ファイルでのトポロジー不変を検査。整理後は全ファイルのバイト一致で両側全節点・全セルの同一性を確認する。
- WAVは2000サンプル、100000Hz、mono、16bit PCMのヘッダーと全サンプルバイトを比較。声道なし・正弦波試験の音響波形を通常加振と同一にする要求は置かない。各モードの整理前後を比較する。
- 基準から存在する診断未定義値を隠さない。contact_couplingの初回予測変位差、接触なし時のcontact_iteration/final最大接触位置・圧力等、ゼロ変位時のmodal_top10の寄与/相殺比はNaNを出す。正弦波モードでは流体分離状態を更新せず、separation.datのx_sep/x_blend_end/p_sepおよびflow_diagnosticsの対応欄がNaN（has_nonfinite=1）。主な物理時系列・変形座標の非有限値とは区別し、整理前後の位置・表記まで一致を要求する。

| 比較段階 | iforce=0 | iforce=1 | 許容誤差 | 最大絶対差／最大相対差（有限数値） |
|---|---|---|---|---|
| simulation.cpp整理後 vs baseline | 136ファイル一致 | 136ファイル一致 | バイト一致（abs=0、rel=0） | 0／0 |
| 最終整理後 vs baseline | 133ファイル一致 | 133ファイル一致 | バイト一致（abs=0、rel=0） | 0／0 |

baselineは各137ファイル。中間段階はdebug_mode_summary.csvだけ非生成、最終段階はdebug_mode_summary.csv/debug_force.txt/contact_monitor.csv/contact_search_debug.csvの4ファイルだけ非生成であり、他の欠落・追加はない。各モード133ファイルには、VTU100、DAT10（SurfAreaを含む）、CSV20、TXT2、WAV1を含む。NaNの位置・文字列、整数・時刻列、WAVのヘッダー/全サンプル、VTUの全Points/Cellsもバイト単位で一致した。非時間計測のstdoutも両段階・両モードで一致した。新たなコンパイル警告は各ビルドログで検出されなかった。

検証証拠は `/tmp/simulation-cleanup-1cGHFz/` の `baseline-validation.json`、`stage-simulation-comparison.json`、`final-comparison.json`、`final-{fluid,sine}-output-sha256.json`、各configure/build/ctest/実行のlog/json、`source-audit.json`、`protection.json`。再検査用の一時スクリプトは同所の `check_cleanup.py`（リポジトリへの追加なし）。一時結果は同所に保持し、研究用の恒久出力ディレクトリには追加していない。

## 未検証・保留と完了状態

0.02秒を超える長期挙動、0.120–0.135秒のcontact_debug区間、声道あり、forceDirection=1、ncont=0、PROFILE_CALCDIS=ON、AP除外・detail、固定境界ID監査、擾乱・荷重凍結の実行は今回の比較条件に含めない。関連処理は保持する。全Python解析スクリプトの実行検証も行わない。

Phase 1目録と、許可リスト内のソース整理・両モードの中間/最終回帰比較を完了した。0.02秒、全100モード、上記入力・OpenMP条件において、維持出力の最大有限数値差は0、全ファイルはバイト一致。これを長期挙動や未実行のオプションへの実行検証と拡張して解釈しない。

既存不整合の修正、参照先のあるログの削除、未検証区間のログ削除、構造変更は今回実施せず、別途指示が必要な保留事項として残す。ソースは3ファイルのみ編集、新規Markdownは指定の2件のみ。本番param・モデル入力・既存の意図的変更・AGENTS・指示書・CMake・tests・tools・.gitignoreは開始時から不変。元のoutput/result全351ファイルもSHA256・サイズ・更新日時が一致した。`git diff --check`は合格。commit/pushは行っていない。

## ソース差分の監査

開始時コピーとの比較では、`src/simulation.cpp` の `Simulation::run` から86行、`src/forceCalculator.cpp` の initialize/applyFluidLoads/applyContactLoadsから60行、`include/ForceCalculator.h`からログ専用ofstreamメンバー3行を除去した。合計149行の削除だけで、ソースの追加・置換はない。

`source-audit.json`で、VTU関数本体と有効呼び出しブロック、正弦波分岐、全OpenMP pragmaの不変を照合した。接触探索・力の加算・クランプ・収束・積分・更新順はdiffの目視でも変更なし。既存パラメータ・全テスト・SimulationParams・CMake・tools等のハッシュは開始時と一致し、検証したコピーのソース/ヘッダー/tests/CMakeも元リポジトリの最終状態と一致した。整理対象の専用メンバー/ラムダは実装参照から消えている。旧コメント中のdebug_forceパス等は動作しないためそのまま残した。

元リポジトリの `build/` は更新していない。検証ビルドは全て隔離コピー内で行ったため、手元の既存実行ファイルに今回の整理を反映するには通常のRelease再ビルドが必要。
