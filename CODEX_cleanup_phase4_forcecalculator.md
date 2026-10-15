# Phase 4：ForceCalculatorの整理・分割（Codex実装指示）

基準はGitHub `main` の `3e410ad2b901f3425e223380a0cb955d2cf24007`（2026-10-15）。実作業ツリーと`AGENTS.md`、`cleanup_phase3_report.md`を必ず先に確認する。これは**実装修正・ビルド・回帰比較まで行う指示**。調査だけで終了しない。Phase 3の許容差は今回に自動適用しない。

## 目的と優先順位

現行`src/forceCalculator.cpp`は約1245行。`applyFluidLoads`が約270行、`applyContactLoads`が約810行を占める。物理式を変更せず、(A)明白な残骸や未使用状態を削り、(B)関数の定義位置を複数CPPへ移す。その上で(C)巨大接触関数の補助処理を**1件ずつ**抽出する。Bまでで十分に見やすくなりCの数値的安全性が担保できない場合、Cは無理に実施しない。

## 絶対維持

- `iforce=0`の流体・構造連成、`iforce=1`の正弦波加振、左右の力の向き、`forcef/famp/forceDirection`。
- `FlowModel`、`ChannelSections`、分離点選択とfallback、圧力回復クランプ、smoothstep、変形面法線による流体圧荷重。
- 接触点探索・候補キャッシュ・two-pointer・非単調形状fallback、負間隙評価、接触代表点、原形状に基づく法線方向、`kc1/kc2/kc3`と`contactReferenceFrequencyHz`、力の分配。
- `contactFlag`、`max_force_diff`、`contact_force_residual`、`contact_penetration_residual`、`max_contact_penetration`、`previous_contact_penetration_`、`fdis*`、`prevFdis*`、`contactPairCache`と全接触収束ロジック。
- `diagnosticOnly=true` の呼び出しと副作用。名称は診断でも、Phase 3の調査で`fdis`・`contactForce`の再計算が確認された。今回は**この呼出し・分岐を削除禁止**。
- `apContactDiagnosticEnabled`とAP端の**物理的荷重除外**分岐。単なる診断用と誤認して削除しない。
- `Simulation::run`の呼出順、Newmark法、モード射影、OpenMP`schedule(static)`、積分定数、数値演算の評価・加算順。
- `*.dat`、`*.wav`、`*.vtu`、`surfp_output.csv`、`params_used.txt`、`manifest.txt`の内容と生成条件。再び他のCSVを作らない。本番`input/`、`output/`、`result/`、`analysis_runs/`、`diagnostic_runs/`は変更しない。
- `PROFILE_CALCDIS`は現段階では残す。将来の削除は別途判断。コメントアウトされた旧デバッグ表示は削除してよい。

## 0. 事前記録・ベースライン

1. `git status --short`、`git rev-parse HEAD`、`git diff --binary`、`git diff --cached --binary`を記録。未コミット変更を誤って破棄しない。編集箇所に競合する未知の既存変更があれば停止。
2. `src/include/tests/tools`で参照検索し、`ForceCalculator`の公開メンバーがどこで読まれ、どこで書かれるか一覧化。`FlowModel`や`Simulation`との依存も確認。
3. 編集前にReleaseビルド・CTest・通常/正弦波の基準計算を実施。以前の比較条件を参考に、通常は少なくとも0.1362秒、正弦波は0.04秒まで、可能なら長時間条件も隔離して比較。`nmode`、`dt`、`Ps`、`kc`、減衰、OpenMP数/配置、コンパイラを同一に。実行は`../output`の幾何直接出力も含め、完全に隔離した作業コピー/実行cwdを使う。
4. 本番`input`と既存`output/result`のSHA256、サイズ、mtimeを前後で監査。旧出力を消さない。本番paramを短縮するのではなく、検証専用コピーを作る。

## A. 明確な不要コードの削除（最初に完了・回帰比較）

- `src/forceCalculator.cpp`の`debug_force.txt`に関する大量のコメントアウト済コード（旧320〜365行付近）、不必要なコメント済`std::cout`、崩れた空行・古い手順だけのコメントを削除。ただし法線の向き、単位、分離点・数値保護など、物理的根拠のあるコメントは保持。
- `contactForceL_ij`/`contactForceR_ij`：現行では確保・ゼロ化・接触力の大きさの格納のみで、内部読出しが見つからない。**`src/include/tests/tools`および公開API依存を確認して完全に不要と証明できた場合だけ**削除。`addSegmentForceXY`から`contactForce`引数と`sqrt(Fx*Fx+Fy*Fy)`の格納のみ外し、`fdisXL/YL/XR/YR`への重み付き加算は不変。
- `minHarea`の全ステップ履歴：現行では書込みはあるが読込が見つからない。外部参照まで確認して不要ならメンバー、サイズ確保と履歴代入を削除。ただし`findMinHarea()`、`minA`の物理計算は残す。
- `lg`（初期化時のz範囲計算）など、一度書いて参照されない変数は、全参照検索で立証できた場合だけ撤去。`nxsup`などログ出力にも使われるメンバーを機械的に削除しない。
- `[DEBUG] FIRST NORMAL CONTACT DETECTED`の初回メッセージと専用`is_first_contact`は、制御・接触判定と無関係と確認した場合に撤去可。標準出力の意図的差分を報告し、実際の接触開始ステップは不変であることを検証。
- 真の未使用ローカル変数・`#include`のみ削除。`psurf`、`Ug`、`currentPg`、`previousUg`、`rho`、`mu`、`harea`、`fdis*`等を見た目で不要扱いしない。

この段階でビルド/CTestと両加振モードの全保持出力を基準と比較。数値差が生じたら原因を特定するまでBに進まない。

## B. 物理式を触らない関数単位のCPP分割

段階A合格後、**関数の中身・ラムダ・ループの演算を変えずに**定義を移す。

- `src/forceCalculator.cpp`：コンストラクタ、`initialize`、`projectLoadsToModes`、`updateChannelSections`、必要な小関数。
- `src/forceCalculatorFluid.cpp`：`applyFluidLoads`。`pressureAtSurfaceX`と`applyPressureForce`のラムダは**当初そのまま**移す。
- `src/forceCalculatorContact.cpp`：`applyContactLoads`、`resetPreviousContactForce`。`processPair`、`computeXYNormal`、`computePenetrationCentroid`、`buildSegmentsForJ`等のラムダは**当初そのまま**移す。
- `findMinHarea`は所属するデータ依存に応じて適切なCPPに一度だけ定義。`CMakeLists.txt`へ新CPPを追加し、多重定義や未解決参照を防ぐ。
- `#include`と`PROFILE_CALCDIS`のビルド条件を各CPPへ移し、ON/OFF両構成を維持する。

関数を新しい`ContactSolver`クラスに書き直したり、構造体のレイアウト・値渡し・ループ順・浮動小数点式・OpenMPを変更してはいけない。**ファイルを分けるだけでもコンパイラ最適化に伴う丸め差が生じ得るため**、B単独で再度回帰比較する。物理式に差を入れて出力を合わせることは禁止。

## C. 800行の接触関数の内部整理（安全が確認できた場合のみ）

B合格後、以下から小さな1機能ずつ抽出し、その都度テスト：
1. `ContactSeg`と`buildSegmentsForJ`（線分生成）。
2. `evalPredAtX`、`computePenetrationCentroid`（補間、負間隙区間）。
3. `computeXYNormal`（左右の**原形状方向**に基づく法線）。
4. `rangesAreOrdered`と候補列挙（cache、two-pointer、fallbackの**列挙順を維持**）。
5. `processPair`と`addSegmentForceXY`（接触圧、面積、荷重分配、**加算順を維持**）。
6. `fdis`から`fx/fy`への反復差加算と残差更新。

一度にすべてを抽出しない。外部関数化による引数コピー・一時オブジェクト・演算再結合・コンパイル最適化差を検査。性能改善や接触アルゴリズム改良は今回の対象外。Cで説明できない差があればCを止め、A/Bまでを成果として保持。

## 必須の検証と停止条件

A/B/CそれぞれでReleaseビルド・CTest（全3件）と`iforce=0/1`を同一条件で比較。保持DAT/WAV/`surfp_output.csv`/TXTは原則バイト一致。VTUは全ファイル・全節点xyzの最大差/RMS、セル構造、名前・生成step・非有限値を比較。接触はステップごとの有無、反復回数、終了理由、最大めり込み、force/penetration residualを検証。Phase 3でCSVを削除したので、必要なら**隔離コピーにだけ一時計測**を作り、本体にCSVを復活させない。

**Phase 3で例外許容した`1e-12 mm`のVTU閾値は、今回の自動合格基準ではない。** まずバイト一致を試み、差が生じたら最初の差の場所・接触力・コンパイラ挙動を調べる。原因不明・差分拡大・接触判定変更・CTestの新規失敗があれば安全停止。意図的なデバッグ標準出力削除だけは文章で説明すればよい。`git diff --check`、`PROFILE_CALCDIS=ON`のビルド確認も行う。

## 報告・完了

`cleanup_phase4_forcecalculator_report.md`に、A/B/Cの実施範囲、削除したフィールド・コメント・関数、実際の参照検索根拠、新CPPへの移動、基準・整理後のビルド/CTest、通常/正弦波の数値・接触・VTU比較、未検証条件、本番データの保護証拠を記載。

**各小作業のたびに削除許可を取り直さず、合格した段階は続行する。** 数値非一致の原因が未確定の場合や保持機能に問題が出た場合は止めて報告。`git commit`・`push`・`reset`・`clean`・`stash`・ディレクトリの広範な再構成はしない。
