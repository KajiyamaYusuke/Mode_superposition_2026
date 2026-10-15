# Phase 4 ForceCalculator整理・段階検証報告

作業日：2026-10-15 JST。指示書：`CODEX_cleanup_phase4_forcecalculator.md`。**段階Aの隔離実装・両加振の検証で数値差とRMS差の増加を確認し、安全停止した。Aは未合格・本体未反映、B/Cは未着手、Phase 4は未完了。** Phase 3の許容差・RMS微増の例外承認は今回へ自動適用しない。

## 基準とデータ保護

`git fetch origin main` を実行し、HEAD・origin/main・FETCH_HEADが全て `5fe929b4fe76a8e4a28f193343fef88f7d8409ff` と一致した。開始時はmainブランチで、status・diff・cached diffは空。共有worktreeのGit管理領域が通常サンドボックスでは読取専用だったため、許可された拡張権限でfetchのみ実行した。更新するコミット差がなく、merge/rebase等は不要だった。commit/push/reset/restore/checkout/stash/cleanは行わない。

隔離検証先は `/tmp/simulation-phase4-S0Na8t`。`start.json` に開始Git状態、trackedファイルのSHA256、本番input55ファイルと既存output/result351ファイルのSHA256・サイズ・mtime、analysis_runs/diagnostic_runsのサイズ・mtime、コンパイラ・環境を保存した。既存データ・本番buildは検証に使わず、上書きしない。

編集前ソースは `original/`、候補は `work/`、段階ごとのソースと実行ファイルは各段階ディレクトリへ保存する。通常と正弦波それぞれの実行cwdは固定した `cases/{fluid,sine}/build`。幾何の `../output` 直接出力もこの中で隔離される。実行終了ごとにその段階の出力を別ディレクトリへ保存し、基準・前段階を上書きしない。候補は検証に合格してから本体作業ツリーへ反映する。

## 公開メンバー・依存関係の参照調査

`src/include/tests/tools` の全C++/ヘッダー/PythonとCMakeを検索。全43公開データメンバーについて出現ファイル・行・内容とForceCalculator外の参照を `reference-inventory.json` に記録した。以下は意味を読んで分類した一覧であり、単に外部で名前が出ないことだけで共用物を削除しない。

| メンバー群 | 書込み／初期化 | 読出し・依存 | 方針 |
|---|---|---|---|
| fxL/fyL/fzL、fxR/fyR/fzR | SurfaceLoadへの参照、流体荷重・接触反復差 | ModalProjectorが同じSurfaceLoad保存配列を読む | 保持 |
| fiL/fiR | モード射影 | Simulationの左右Newmark | 保持 |
| psurf | 流体圧再構成、クランプ・ブレンド | 表面圧荷重とpressure.dat | 保持 |
| Ug、previousUg | 流量履歴、前step取得 | 次のFlowModel.advance入力 | 保持 |
| currentUg/currentPg | FlowModelの流量・声門圧 | 圧再構成、履歴、流量DAT | 保持 |
| harea | ChannelSectionsの断面積への参照 | findMinHarea・FlowModel・分離選択・area.dat | 保持 |
| fdisXL/YL/XR/YR | 接触の重み付き加算、診断時再計算 | 反復差・残差・表面力更新 | 保持 |
| prevFdisXL/YL/XR/YR | reset、反復力更新 | 次反復との差・収束 | 保持 |
| contactFlag | 接触判定 | Simulationの収束・非接触終了、診断保存復元 | 保持 |
| max_force_diff、contact_force_residual | 反復差の最大値・正規化 | 接触収束、診断保存復元 | 保持 |
| contact_penetration_residual、max_contact_penetration、previous_contact_penetration_ | めり込み最大値・前回値更新 | 接触収束、診断保存復元 | 保持 |
| contactPairCache | 第1反復で候補を登録 | 後続反復の候補順、diagnosticOnly時更新禁止 | 保持 |
| rho/mu | spから初期化 | Bernoulli・粘性圧力式 | 保持 |
| geomL/R、modeDataL/R、stateL/R、sp | コンストラクタの参照 | 全計算・モデル・境界設定 | 保持 |
| nxsup | initialize | 初期化標準出力にも使用 | 保持 |
| contactForceL_ij/R_ij | 各5箇所：宣言・確保・ゼロ化・lambdaへの参照渡しと大きさ格納 | ForceCalculator外0件。lambda内+=は累積用自己読出しだけで、その結果を物理計算へ渡す消費者なし | Aの削除候補 |
| minHarea | 4箇所：宣言・確保・両加振で履歴代入 | ForceCalculator外0件、履歴読出し0件。findMinHareaと物理的minAは別 | Aの削除候補 |
| lg | 2箇所：宣言・z範囲初期化 | ForceCalculator外0件、値の読出し0件 | Aの削除候補 |

ForceCalculatorを所有するのはSimulation。Geometry/ModeData/State/SimulationParamsへの参照を保持し、ChannelSectionBuilder、SurfaceLoad、ModalProjector、FlowModel、SeparationSelectionを内部で使用する。CMakeは実行ファイルターゲットのみ、include設定はPRIVATEで、公開ライブラリ/export/install/bindingはない。削除候補を利用する公開API消費コードはリポジトリ内に見つからなかった。管理外の外部プログラムの不存在までは保証しない。

## 検証方法

Intel icpx 2025.1.0、C++17、Release `-O3 -std=gnu++17 -fiopenmp`。コンパイラ既定の浮動小数点最適化を変更せず、fast-mathやstrict等の追加を主比較に使用しない。PROFILE_CALCDIS=OFF/ONをそれぞれビルドしてCTest3件を確認する。

全比較でOpenMPは `OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores OMP_WAIT_POLICY=PASSIVE`、SIMULATION_RUN_DIR未設定。待機方針も基準から全段階で固定し、性能評価は行わない。

試験paramは本番ファイルのコピーで、通常nstep=13620（0.1362秒）、正弦波nstep=4000（0.04秒）と実効iforce=1だけを変更する。100モード、dt=1e-5、Ps=2900、zetaL/R=0.02、kc=(5,6000,5000)、ncont=100、forcef=100、famp=-1e-5、forceDirection=0は統一。本番paramは変更しない。3種類の利用モデルは隔離先へコピーしてハッシュを照合した。

保持DAT/WAV/surfp/TXT/VTUはまずバイト一致で比較し、現時点の絶対・相対許容誤差は0。差があれば全VTU座標の最大絶対差/RMS・最初の差・時間推移、接触記録を測定して中断判断する。Phase 3の1e-12 mmを自動合格基準にはしない。

接触の全反復計測は隔離先だけのリンク時ラッパー `contact_probe.cpp` を使用。通常実行ファイルと同じコンパイル済み本体オブジェクトを使い、applyContactLoadsの本体実行後、step・反復番号・diagnosticOnly・contactFlag・5残差/めり込み値と配列の生バイトハッシュを96byteレコードへ保存する。ラッパーは物理配列を書き換えず、浮動小数点演算も行わず、記録にヒープを確保しない。C++本体ソース・Simulationの呼出順・計算式は計測のために編集しない。配列ハッシュはfdis・prevFdis・合成表面力・一般化力・候補cacheの読取照合であり、数学的な無衝突証明ではない。CSVを復活させない。

基準VTU全ファイルの座標有限性、90666点・496536セル、connectivity範囲、offsets/types数とトポロジーを検査。出力条件・番号は通常n=8620–13600で250本、正弦波n=0–3980で200本。接触は全反復の有無・回数・終了理由・最大めり込み・force/penetration residualを照合し、diagnosticOnlyが保存復元する値とcache・一般化力・合成表面力もチェックする。

## 段階別の実施・結果

編集前の基準計算と出力検査は完了した。基準ReleaseとPROFILE_CALCDIS=ONのビルド・CTestはどちらも成功3件。通常は終了0（634.37秒）・VTU250本、正弦波は終了0（295.06秒）・VTU200本。時間は並行実行・一時計測込みの参考値であり、性能比較に使わない。

基準の接触計測は通常反復31422・診断呼出し13620・接触開始n=1198（0.01198秒）・収束12422step・非接触1198step・上限0・最大反復3。正弦波は通常反復10947・診断呼出し4000・接触開始n=5（0.00005秒）・収束3837step・非接触163step・上限0・最大反復3。全残差・めり込み値は有限で、全VTU・DAT・WAVの構造検査を通過した。

一時計測による出力変化がないことの追加観察として、同一のsrc/include/tests/CMakeと同じ試験paramで既に保存されていた計測なし実行と比較した。通常の保持263ファイル、正弦波の保持213ファイル（VTU含む、cwd依存manifestは除外）が全てバイト一致した。これは今回の基準計算を省略・再利用したものではない。過去実行とのOMP_WAIT_POLICYの違いも `probe-check-{fluid,sine}.json` に明記した。

段階Aを隔離コピーへ実装し、Release/PROFILE_CALCDIS=ONともビルド・CTest3件に成功。通常/正弦波の回帰計算は両方終了0、通常624.47秒・VTU250本、正弦波293.82秒・VTU200本。しかし数値比較が不合格であり、A/B/Cは本体へ反映していない。

段階A候補は、旧debug_force/std::cout表示コメント、初回接触メッセージと専用bool、不要な履歴・大きさ格納配列・lg、未使用ローカル5件とinclude2件、古い手順コメント・空行/字下げの整理。接触荷重の4つの重み付き加算、最小面積の物理計算、法線・ペナルティ・減衰・残差・cache・AP物理除外・diagnosticOnly・PROFILE計測は保持する。候補差分生成時に重要演算ブロックの空白を除いた同一性を検査した。

削除候補の未使用ローカル5件はinitializeのnsurfl_L/nsurfz_L/nsurfl_R/nsurfz_RとapplyFluidLoadsのnsurfl。別途lg専用のzmin/zmaxと範囲走査も除去。include2件はDisplacement.h（直接型参照なし、State等から定義を取得）とsstream（使用なし）。候補の行数はforceCalculator.cppが1244→1109（追加30/削除165）、ForceCalculator.hが126→123（追加1/削除4）。追加行は引数削除・字下げ等に伴う差分で、物理式の新規追加はない。PROFILE_CALCDISの計測・集計・条件は残し、コメントアウト済みの表示だけ除去した。contactForce配列のゼロ化にしか使わないOpenMP領域2件は配列と共に除去し、保持する物理処理のOpenMP配置とschedule(static)は維持した。

検証ヘルパーの初期実行では、通常計算終了前の出力を開こうとしたため未生成WAVで停止した。これは基準実行の異常ではなく、終了後に再検査した。また、最終通常接触呼出しとdiagnosticOnlyの間には元から最終モード射影があるため、その間のfiハッシュまで不変とする当初のヘルパー仮定を修正した。fiは同じ呼出し位置の基準/候補間で比較し、診断の保存復元値・cache・表面力・prevFdisは区間内でも検査する。計算コードや出力は変更していない。

作業中に未追跡 `tools/testFluid.py` の追加を検出した。内容は保持DATによる独立した流体抵抗比較で、削除候補の公開メンバーへの参照はなく、整理対象と競合しない。変更・削除・実行せず保護し、参照検索を再実行して削除候補の外部参照が依然0件であることを確認した。

### 段階別判定

| 段階 | 実施範囲 | ビルド・CTest | 通常/正弦波比較 | 本体への反映 |
|---|---|---|---|---|
| A | 上記の不要配列・履歴・lg・局所変数・旧デバッグ残骸を隔離コピーで削除 | Release/PROFILE_CALCDIS=ONとも成功、各CTest3件成功 | 両方でVTUと接触残差等に微差、RMS増加。STOP_DIFFERENCE | 未反映。候補を隔離先に保持 |
| B | A不合格のため未着手 | 未実施 | 未実施 | 新CPP追加・CMake変更なし |
| C | A不合格のため未着手 | 未実施 | 未実施 | 補助処理抽出なし |

### Aの保持出力・VTU・接触比較

通常/正弦波とも保持DAT10本、WAV、surfp_output.csv、params_used.txt、manifest.txtは全てバイト一致（14ファイル）。通常はVTU250本全て、正弦波は198本に座標差があり、正弦波の0000/0001だけVTUもバイト一致した。全ファイルの名前・数・生成step・節点対応は一致し、Points以外のXML全文が各対応ファイルで完全一致した。従って点数・セル数・connectivity・offsets・typesは一致。新たな非有限座標はなく、基準で存在する正弦波separation.datのNaN位置・表記も一致した。

| 加振 | VTU最大絶対差 mm | 全座標RMS差 mm | x/y/z最大絶対差 mm |
|---|---:|---:|---|
| 通常0.1362秒 | 2.1316282072803006e-14 | 1.874230910203388e-15 | 2.1316282072803006e-14 / 2.1316282072803006e-14 / 1.0658141036401503e-14 |
| 正弦波0.04秒 | 1.0658141036401503e-14 | 3.2832062849008464e-16 | 全成分1.0658141036401503e-14 |

通常のxyz別RMSは `2.3775751413510844e-15 / 1.7100723367863343e-15 / 1.40036194052572e-15 mm`。正弦波は `3.5490271456906004e-16 / 3.3259496105938595e-16 / 2.9463190653048824e-16 mm`。

最初の保存VTU差は通常 `deform_combined0000.vtu`、n=8620（0.0862秒）、0-based節点4・x：`6.874258468253350` → `6.874258468253349`。通常の保存開始はn=8620なので、それ以前の最初の座標差はVTUからは分からない。正弦波は `deform_combined0002.vtu`、n=40（0.0004秒）、節点9203・z：`0.7500458158412939` → `0.7500458158412938`。VTUは状態確定後の出力で、公称状態時刻は(n+1)dt。

| 実ループn範囲 | 最大座標差 mm | 全成分RMS差 mm |
|---|---:|---:|
| 通常8620–9860 | 1.0658141036401503e-14 | 1.5946288559295667e-15 |
| 通常9880–11120 | 2.1316282072803006e-14 | 1.7267290107410787e-15 |
| 通常11140–12360 | 2.1316282072803006e-14 | 1.9012431893315657e-15 |
| 通常12380–13600 | 2.1316282072803006e-14 | 2.2217147775976408e-15 |
| 正弦波0–980 | 1.0658141036401503e-14 | 1.7202010355557245e-16 |
| 正弦波1000–1980 | 1.0658141036401503e-14 | 3.1307529061109353e-16 |
| 正弦波2000–2980 | 1.0658141036401503e-14 | 3.594700181595679e-16 |
| 正弦波3000–3980 | 1.0658141036401503e-14 | 4.175547781960663e-16 |

通常のRMSは最初→最後の区間で約39%増加し、最大差の包絡も最初の区間の約2倍となった。正弦波は最大差の包絡が一定だがRMSは増加した。値が小さいことだけで合格にせず、指示書の差分拡大の停止条件にも該当すると判断する。

全反復のstep・反復番号・diagnosticOnly・接触フラグ・候補cacheのハッシュは一致し、**全stepの最終反復数と終了理由（収束/非接触/上限）も一致**。上記基準の接触開始時刻・反復回数・収束/非接触/上限集計はAでも同じ。標準出力のContact Iteration Summary全項目も一致した。ただし接触力配列のハッシュと以下の残差・めり込み値は完全一致ではない。

| 接触記録の最大絶対差 | 通常 | 正弦波 |
|---|---:|---:|
| max_force_diff N | 1.0977546996415732e-18 | 1.9809136020962695e-17 |
| contact_force_residual（無次元） | 1.079302005706133e-13 | 4.826490088534951e-12 |
| contact_penetration_residual m | 4.797594613248357e-18 | 6.776263578034403e-18 |
| max_contact_penetration m | 3.903127820947816e-18 | 4.336808689942018e-18 |
| previous_contact_penetration_ m | 3.903127820947816e-18 | 4.336808689942018e-18 |

この表は基準/候補とも生double記録による比較で、十進CSVの表示丸めは含まない。表のmax_force_diff差は残差計算に用いる反復力差の最大値の差であり、全接触力配列の要素別最大差と同じ意味ではない。fdis全体はハッシュで非一致を検出したが、全要素を個別出力していないため、配列の要素別最大差は未測定。DATの表示精度で同じであることを内部状態の完全一致の証明にはしない。

最初に差のある接触計測箇所は、通常n=1198・反復1・diagnosticOnly=false（接触開始と同じ）。max_force_diffは `1.8104076455191964e-09` → `1.8104076455192104e-09 N`、force residualは `0.5507263330674618` → `0.5507263330674633`、最大めり込みは `1.2315877366830818e-08` → `1.2315877367431034e-08 m`。この時点で接触フラグ・合成表面力/一般化力/cacheのハッシュは一致し、fdis・prevFdisと残差に差があった。

正弦波の最初の差はn=5・反復1・diagnosticOnly=false。max_force_diffの値は一致したが、force residualは `0.26223631466036473` → `0.2622363146603718`、最大めり込みは `2.1249385895316663e-10` → `2.1249385895317935e-10 m`。fdis・prevFdis・合成表面力のハッシュも異なり、一般化力/cacheのハッシュはこの時点では一致した。

両計測は接触呼出しの**直後**であり、予測節点座標やqの全値を呼出し直前に記録したわけではない。従って差の生成箇所を接触関数内の単一演算と断定しない。記録した範囲では、これより前の接触記録は一致している。

### 原因の調査範囲と停止理由

保持する物理式・重み付き接触加算・法線・負間隙代表点・残差更新・cache順・AP荷重除外・diagnosticOnlyの保存復元を手作業で改変した差は見つからず、重要ブロックは空白を除いて元と一致する。Aでは未使用の大きさ格納に使うsqrt・2つのクリア領域・初回接触メッセージを撤去し、公開フィールド/割当て数も変わった。これらがコンパイラのインライン化・演算再結合・共通式最適化や配置に影響した可能性がある。

同じ主ビルドフラグの `icpx -###` には既定の `-funsafe-math-optimizations`、`-mreassociate`、`-ffp-contract=fast` が含まれた（`compiler-driver.json`）。この設定を今回追加/変更したのではない。本体オブジェクトのSHA256が基準とAで変わったのはforceCalculator.cpp、ForceCalculator.h依存のsimulation.cppとmain.cppのみで、State/FlowModel/ModalProjector/Newmark等のオブジェクトは同一だった。

観測した微差は最適化の変化に伴う丸め差と整合的だが、**どの削除項目が起点か、コンパイル結果と実行時の配置のどちらが寄与したかは未特定**。Phase 3の原因切り分けを今回の証明として流用しない。今回は候補を複数項目まとめてAとして作成したため、単一項目への追加分離やstrict等の対照実行は未実施。差と増加を検出した段階でB/Cおよび本体反映を止め、その後は保存済みファイル・記録・ソース・コンパイラ設定の読み取り整理だけを行った。新しい許容差の設定、出力桁数の変更、物理式による合わせ込みは行っていない。

意図的な標準出力変更は `[DEBUG] FIRST NORMAL CONTACT DETECTED` と専用初回boolの撤去。実際の接触開始stepは両モードで不変だった。タイミング値は実行ごとに変化するが、計算結果の差を隠すために標準出力を変更していない。

### 保存した成果と保護確認

- A候補は `/tmp/simulation-phase4-S0Na8t/A/src/forceCalculator.cpp` と `A/include/ForceCalculator.h` に保存。work内の候補はこの実際に検証したスナップショットと一致する。本体の同2ファイルはoriginalと完全一致し、候補を未合格のまま本体へ移していない。
- `comparison-A-fluid.json` / `comparison-A-sine.json` に全ファイル一致/不一致一覧、全VTUフレームの最大差/RMS、接触記録比較を保存。`A-analysis.json` に時間区間集計、最初の接触記録値、変更されたオブジェクト一覧を保存。`final-source-and-reasons.json` に全step終了理由とソース照合結果。
- `baseline-validation-{fluid,sine}.json`、各段階のbuild/ctest/profile-build/profile-ctestログ、object-sha256.json、`reference-inventory.json`、`params.json`、`probe-check-{fluid,sine}.json`、`compiler-driver.json` を隔離先に保存。計測方法と比較処理は同先 `contact_probe.cpp` / `compare.py` / `verify.py`。
- 停止後の `protection.json` でも、本番input55・既存output/result351のSHA256・サイズ・mtimeは開始時と全て一致。研究アーカイブのサイズ・mtimeも一致。既存trackedファイルの変更は0件。本番buildは使わず、rootの追加ファイルは本報告書と作業中に検出した利用者側のtools/testFluid.pyのみ。後者は未変更・未実行で維持した。
- `git diff --check` は指摘なし。新しい本体CPP/CMake変更はなく、追加CSV生成もない。commit/push、既存差分の退避・破棄・復旧コマンドは実行していない。

## 未検証・完了状態

**段階Aで安全停止・Phase 4未完了。** 原因を厳密に特定せず数値的同等性を合格にせず、A候補は隔離先へ保持、本体のForceCalculator実装は整理前のまま。段階B/Cの分割・抽出は未実施。追加の原因切り分けや許容基準について利用者の指示を待つ。

比較を実施した範囲は通常0.1362秒、正弦波0.04秒。これを超える長時間・声道あり・forceDirection=1・AP除外ON・ncont=0・PROFILE_CALCDIS=ONでの時間計算は未検証。今回の候補に対して新たな非0の許容差や例外承認は設定していない。物理精度や長時間の数値同等性を保証しない。自動コミット・pushは行っていない。
