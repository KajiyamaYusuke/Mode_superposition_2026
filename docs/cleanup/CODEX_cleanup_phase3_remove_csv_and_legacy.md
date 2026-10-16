# Codex作業指示：CSV全廃・擾乱診断・旧関数の一括削除

## 0. 権限と基準
このファイルは**コードを実際に変更する指示**であり、調査・提案だけで終わらせない。最新の実作業ツリーとルート`AGENTS.md`が優先。参照時点のGitHub mainは`7b31b648e84bc574cc4b1b0b8055d5bbce1a2568`。以前の`CODEX_cleanup_phase1_2.md`と衝突する場合、今回明示的に許可した削除だけ本書を優先する。

目標：C++ランタイムが出すCSVを`surfp_output.csv`以外**全部削除**し、専用の出力・集計・メンバー・不要な計算も撤去。さらに擾乱・荷重凍結・基準窓・擾乱エネルギー評価および指定旧関数6件を撤去する。関連する宣言・パラメータ解析・設定表示・テストも整合させる。

範囲はC++シミュレーターのCSV。Python解析ツールが独立生成する派生CSVを全部消す依頼ではない。既存の`output/`、`result/`、`analysis_runs/`、`diagnostic_runs/`の実データを絶対に削除しない。

## 1. 必ず維持するもの
- `Geometry::surfExtractFromNAS`が書く唯一のCSV `surfp_output.csv`。出力内容、パス、生成タイミングを変更しない。
- `Simulation::writeVTKCombined`と有効な時間ループ呼出し。出力条件（前回は`n%20==0 && nstep-n<=5000`）、変形座標、セル情報、VTU名を維持。
- `iforce=0`の通常流体・構造連成、`iforce=1`の正弦波加振（`famp`、`forcef`、`forceDirection`、旧/実効iforce入力互換）。
- `*.dat`（`area.dat`、`displace.dat`、`displace_xy.dat`、`pressure.dat`、`pressure_vt.dat`、`airflow_vt.dat`、`section_x.dat`、`gap_cubed.dat`、`separation.dat`、`SurfArea.dat`等）、`*.wav`、`txt/manifest.txt`、`txt/params_used.txt`。
- 幾何/NAS/モード読込、モード射影、流体計算、分離点、圧力/表面力、接触候補探索・力加算・キャッシュ・**収束残差**、Newmark積分、声道/声門下モデル。
- `input/param.txt`および全研究入力の内容。本書作成時のGitHub値は`nstep=100000, zetaL=zetaR=0.02`だが実作業ツリーの値を優先。数値式、`beta=0.275625`、`gamma=0.55`、dt、ループ/演算順、OpenMP配置も変更禁止。

コードのディレクトリ移動・rename・最適化・依存追加・全体のPython書換え・git reset/clean/restore/stash・自動commit/pushは禁止。

## 2. 事前検査、基準実行
1. `git status --short`、`git rev-parse HEAD`、`git diff --binary`、`git diff --cached --binary`を保存。未コミット変更を勝手に書き換えない。今回の変更箇所と競合する予期しない差分があれば停止。
2. `cleanup_inventory.md`と現行実装の入出力一覧・参照検索を照合。コードが違えば現行を採用。
3. **編集前**にReleaseビルドとCTest、`iforce=0`/`iforce=1`の基準計算を行う。nstepを短縮した試験専用paramコピーを使い、nmode、dt、Ps、zeta、kc、力の条件、OpenMP条件は統一。
4. **隔離した作業コピー／実行cwd**で検証。`geometry.cpp`は`../output`へ直接書くため`SIMULATION_RUN_DIR`の指定だけでは隔離不十分。本番`output/`・`result/`既存ファイルのSHA256、サイズ、mtimeを前後で確認する。
5. 比較対象は「**残す**DAT/WAV/VTU/TXT/`surfp_output.csv`」。削除CSVとの差分は正しい変更として扱う。

## 3. CSV出力の完全撤去
C++の`std::ofstream`、`.open()`、`.csv`、ファイル名定数、関係するコードとヘッダーを全文検索し、**surfpを除くC++ランタイムCSVが一本も生成されない**状態まで削除する。単に書込行をコメントアウトせず、関連する**専用の**ストリーム・初期化・ヘッダー・集計ラムダ・構造体・変数・演算も撤去。

### 3a. `src/simulation.cpp`
少なくとも次のCSVの生成と専用処理を削除する。
- `fixed_boundary_mode_audit.csv`
- `modal_contribution.csv`、`modal_dominant.csv`、`modal_top10.csv`
- `irregularity_timeseries.csv`、`flow_diagnostics_v2.csv`
- `modal_energy.csv`、`energy_summary.csv`
- `gap_field.csv`、`gap_field_summary.csv`
- `modal_force_split.csv`、`diagnostic_events.csv`
- `perturbation_energy.csv`、`perturbation_energy_summary.csv`
- `contact_coupling_debug.csv`、`debug_step_summary.csv`
- ほかに実際のソースで発見したC++生成CSVのすべて。

`writeModalDiagnostics`、`writeTop10ModalContributions`、`precomputeSurfaceModeRms`、`writeGapSnapshot`、これら専用の`SurfaceModeRms`、`ModalContributionEntry`等は、他の参照がなければ削除。**変位DATを生成する`nearestIdxL/R`や、`area.dat`を生成する面積処理は残す**。

### 3b. `src/forceCalculator.cpp` / `include/ForceCalculator.h`
少なくとも`contact_debug.csv`、`contact_iteration_debug.csv`、`contact_final_state.csv`、`contact_region_summary.csv`、`contact_pairs_detail.csv`、`debug_flow_detail.csv`を全廃。`debugContact`の固定時間窓（0.120–0.135秒）と、その専用最悪値/回数統計を含む。すでに削除された`contact_monitor.csv`等の残留フックも、他用途がない場合だけ除去。

**保持必須の共用値**：`contactFlag`、`max_force_diff`、`contact_force_residual`、`contact_penetration_residual`、`max_contact_penetration`、`previous_contact_penetration_`、`contactPairCache`、接触力・前回力配列、接触候補探索・接触力配分・キャッシュ・収束判定に使う変数。CSVに書いていたからという理由で物理計算を削除しない。

`setContactMonitor`、`contactMonitorI/J`は、ログ撤去後に完全未参照と証明できれば削除可能。

**重要：`Simulation::run`は前回`fCalc.applyContactLoads(n, icont_used_this_step+1, true)`を毎ステップ呼ぶ。** `diagnosticOnly`という名称だけで即削除しない。内部で接触荷重配列・残差・cache/前回値を変更しうるので、全副作用と整理前後の数値一致を監査する。物理的に同等と確認できた場合のみ呼出し/分岐を削除。確認できなければCSV書込だけ消して該当計算を保持し、レポートで理由を示す。

### 3c. `src/geometry.cpp`
`Geometry::surfExtractFromNAS`の`surfp_output.csv`生成は絶対に維持。旧`Geometry::surfExtract`を削除する際、同名CSVライターのうちNAS側の必要な実装を残す。`SurfArea.dat`も維持。

### 3d. 設定・説明
削除CSV**専用**の設定（例：`diagnosticOutputIntervalSteps`、`contactDetailStartSec`/`contactDetailEndSec`、`fixedNodeIdsFile`）は、完全未参照になった場合に`SimulationParams.h`、`readParams.cpp`の読み込み・検証・表示と`manifest.txt`から整合して削除。入力中に**有効な**削除済みオプションが残る場合、無言で読み飛ばさず、検証前に明示的に報告する。

一方、`apContactDiagnosticEnabled`や`apContactExcludedRowsPerEnd`/`apContactExcludedDistanceMm`は接触荷重の**実際の除外条件**に関係する可能性が高い。`contactIsExcluded`と荷重除外動作は維持。CSV専用のAP分析行/詳細窓設定だけ不要なら削除する。

## 4. 擾乱・freeze・基準窓・擾乱エネルギーを機能ごと削除
`src/simulation.cpp`および関連ヘッダー/パラメータから、次を根から削除する。
1. `perturbationEnabled`、`perturbationTimeSec`、`perturbationModeIndex`、`perturbationAmplitudeAtProbeMm`、`perturbationPattern`、`stateL/R.q`への擾乱加算、イベント、専用`perturbation_manifest.txt`。
2. `diagnosticLoadMode`（`live`以外の分岐すべて）、`freeze_fluid`/`freeze_contact`/`freeze_all`、`applyDiagnosticLoadMode`、`diagnosticReferenceWindowStartSec/EndSec`、`referenceLoadSamples`、`referenceFluid*`、`referenceContact*`、`referenceQ*`等の窓集計・平均・凍結荷重適用。
3. `PerturbationTotals`、`initialPertEnergy*`、`previousPert*`、`cumulativePert*`、擾乱仕事・エネルギー・残差と出力。
4. CSV・凍結・擾乱の**ためだけ**の`projectFluidLoadsToModes`/`fiFluidL/R`、モード力分解や仕事計算等の中間状態は、他用途がないことを確認したものだけ削除。メソッド宣言・定義を揃える。
5. `SimulationParams.h`、`readParams.cpp`、`manifest.txt`の削除機能専用項目/検証も削除。

**通常時（旧`diagnosticLoadMode=live`、擾乱無効）の物理計算をそのまま保持**。接触後の`fiL/R`からNewmarkへの受け渡しを、ログ削除ついでに書き換えない。

## 5. 指定旧関数6件を実際に削除
- `TimeIntegrator::rungeStep()`：`src/rungeStep.cpp`＋`include/TimeIntegrator.h`。**同じCPPの`newmarkStep()`とファイルは保持。**
- `Geometry::surfExtract()`：`src/geometry.cpp`＋`include/Geometry.h`。**`surfExtractFromNAS`/`surfArea`保持。**
- `ForceCalculator::findNsep()`：`src/forceCalculator.cpp`＋`include/ForceCalculator.h`。現行分離点処理保持。
- `ModeData::loadAll()`：`include/ModeData.h`のinline定義。`loadFromVTU`、`loadFreqDamping`保持。
- `FlowModel::findSeparation()`：`include/FlowModel.h`に残る宣言。物理的な分離選択は保持。
- `SurfaceLoad::clear()`：`src/surfaceLoad.cpp`＋`include/SurfaceLoad.h`。`resize`は保持。

参照検索で本当に使用されているものが見つかれば、該当関数だけ保留・報告。他の「未使用そうな」関数を追加削除しない。

## 6. 変更後の検証（作業終了条件）
- ReleaseビルドとCTestを実施。既存`test_input_file_params.cpp`には削除対象`diagnosticReferenceWindowStart/EndSec`を期待する不整合（`0.15–0.19`と`0.25–0.30`）がある。**削除した機能の期待値だけ**最小修正してよい。無関係なテストの削除・skip・テスト偽装禁止。
- `iforce=0`と`iforce=1`で同一基準条件を実行。**保持DAT、WAV、VTU、唯一の`surfp_output.csv`、保持TXTは原則バイト一致**（差があれば理由と最大絶対/相対差を示し、不明なら中断）。接触回数/収束結果と主要時系列も照合。
- クリーンな隔離出力先を再帰探索し、`*.csv`は`surfp_output.csv`のみであることを証明。削除CSVを出していないことを確認。本番の昔のCSVを消して条件達成に見せかけない。
- `VTU`の全ファイル数とPoints/Cells、生成条件、WAV、DAT列数・時刻軸を確認。`git diff --check`・全文参照検索を実施。
- `0.02秒`の短時間検証だけで長時間自励振動や0.12秒以降の挙動を検証済みと書かない。可能なら追加計算、無理なら未検証として明記。
- 新しい数値不一致、VTU欠損、原因不明のコンパイル/テスト失敗、既存データ損傷は安全停止。物理式変更やテスト無効化で隠さない。

## 7. Python解析の互換性を必ず報告
Python `tools/`自体やその派生CSVは今回**削除/変更しない**。ただしC++ CSV撤去により、`analyze_irregularity.py`、`analyze_ap_contact.py`、`run_ap_contact_diagnostics.py`、`analyze_modal_exit_test.py`、`analyze_self_oscillation.py`、`run_pressure_sweep.py`、`kcExplore.py`、モード分析等で入力が失われる可能性がある。現行ツールを全文検索して以下の表を報告：

`Pythonスクリプト | 依存していたC++ CSV | 代替DAT等の有無 | 変更後の動作可否 | 次回作業候補`

ユーザーが明示したのは**C++のCSV出力と列挙した旧計算機能の削除**であり、Pythonの解析方式・jitter定義・波形評価を推測で書き換えたり、動かなくなるスクリプトを黙って削除したりしない。

## 8. 提出物
`cleanup_phase3_report.md`を作成：
1. 削除したCSV全件とC++書込元、唯一残る`surfp_output.csv`。
2. 削除した擾乱/freeze/基準窓/擾乱エネルギー/6旧関数、残した共用物とその理由。
3. ソース差分のファイル別行数、ビルド/CTest結果、通常加振・正弦波加振の保持出力一致表、最大差、未検証項目。
4. Python解析への破壊的変更と今後の対応候補。
5. `input/param.txt`と元の`output/`・`result/`が保護されている証拠。

**実装と基準/回帰計算まで一気に実施**。個々の削除ごとに承認を求めない。数値安全性が確認できない箇所は、その部分のみ保持・明記し、勝手に物理挙動を変えない。**commit/push、ディレクトリ再構成はしない。**
