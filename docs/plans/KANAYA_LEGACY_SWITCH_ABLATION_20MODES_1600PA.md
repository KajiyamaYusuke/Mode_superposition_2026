# Codex実装指示：左右独立C++にYoshinaga旧方式を並設し、自励振動消失の原因をアブレーション比較する

作成日：2026-10-21 / 対象：`Mode_superposition_2026`、左右対称の基準ラン。  
**目的：旧方式で自励振動を再現できたという既存の実績を活用し、左右独立C++での変更のどこが振動特性を変えたか特定する。** モデルを旧版へ全面的に戻す作業ではない。将来の左右非対称計算に必要な現行コードは維持する。

## 0. ユーザーの確定条件・重要な前提

- **左右対称**、同一形状・同一材料・同一固有モード／質量／初期ギャップ／流路の条件で比較。
- **nmode = 20、声門下圧 Ps = 1600 Pa**。20モードは計算短縮、1600 Paは過去に良好だった条件として固定。
- **主比較の構造減衰は左右とも ζ=0.06 で固定**。以前の発振成功計算はζ=0.1だったので、ζ=0.1の確認ランは**別の明示的比較**にし、流体処理の変更と同時に変えない。
- `dt=1e-5 s`、まず `nstep=20000` (0.2 s) を暫定評価窓とする。圧力ランプ、接触係数、接触参照周波数、声道・声門下管路・質量・固有モードファイルは**変更せず、実際の入力をmanifestに記録**。旧成功ランとの条件差も列挙する。
- ユーザーはすでにYoshinagaらのコードの自励振動を再現済み。**再現できた旧実装／入力／出力を比較の真の基準とする**。公開Fortran `yossie333/modeFold_v1.3` は実装参照であって、ユーザーが発振を再現した版と一致するとは限らない。公開版の現時点`input/param.txt`は例としてPs=1800Pa、ζ=0.06、20モード、声道17.5cm等であり、そのまま今回の1600Pa条件とはしない。
- **高速化方針：ローカルPCで使用可能なCPU資源を積極利用し、実測で最速のOpenMPスレッド数とジョブ並列数を選ぶ。** 全論理CPUを常時占有させれば最速とは限らない。物理コア／論理CPU・スレッド割当と計算時間を測定する。
- 本タスクは**実装＋ユニットテスト＋有効性確認＋必要な比較ランの実行と報告**まで行う。数十ケースを盲目的に投下しない。

## 1. 開始前の保護・参照確認（必須）

1. `pwd; git branch --show-current; git status --short; git rev-parse HEAD; git worktree list`。`AGENTS.md`の指示を読む。既存未コミット差分を保護。前回の `/tmp/simulation-main-kanaya-diagnostic-20261021` 等に作業成果がある場合、そこへ勝手に上書きせず**新規専用worktree**（例 `/tmp/simulation-main-kanaya-ablation-20261021`）を推奨。現行コードの状態を取り違えない。
2. 対象ソースの現物とCMake実体を確認する。主要箇所（存在する実際のパスを優先）：  
   `src/fluid/channelSections.cpp`、`src/fluid/flowModel.cpp`、`src/fluid/separationSelection.cpp`、`src/forces/forceCalculatorFluid.cpp`、`src/forces/forceCalculatorContact.cpp`、`src/forces/modalProjector.cpp`、`src/structure/rungeStep.cpp`、`src/structure/modedata.cpp`、`src/simulation/simulation.cpp`、`src/simulation/readParams.cpp`、対応する`include/`。
3. 旧ソース候補（ユーザーの既存成功版があれば優先）と論文系Fortran実装を照合。公開版の参照箇所：  
   [calcArea.f90](https://github.com/yossie333/modeFold_v1.3/blob/main/src/calcArea.f90)、[calcForce.f90](https://github.com/yossie333/modeFold_v1.3/blob/main/src/calcForce.f90)、[calcFlow.f90](https://github.com/yossie333/modeFold_v1.3/blob/main/src/calcFlow.f90)、[step.f90](https://github.com/yossie333/modeFold_v1.3/blob/main/src/step.f90)、[contactForce.f90](https://github.com/yossie333/modeFold_v1.3/blob/main/src/contactForce.f90)。**成功実績のあるソースのcommit/hash、入力値、runディレクトリ、対応出力**を調べ、入手できなければ「未確認」と記載し公開ソースとの違いを隠さない。
4. 旧成功条件と今回条件の照合表を作る：COMSOLモデルと周波数・節点質量正規化・20モードの順番／鏡映、声門下圧、密度、質量、粘性、音速、入口/声門下管/声道の長さと面積、初期gap、ζ、接触係数、接触参照ω、ncont、dt、圧力立上がり、初期流体状態。**「基本同じ」と仮定せず、記録から確認できた/できないを区別**。
5. 既存`output/`と`result/`は絶対に上書きしない。独立runディレクトリを使う。ユーザーの既存実行環境がない場合はパスと不足入力を報告し、虚偽の実行結果を作らない。

## 2. 推奨アーキテクチャ：1本の計算ループ＋個別の実装戦略

ユーザー案のとおり**既存src内に旧実装専用ファイルを追加**し、`SimulationParams`と計算本体から切替える。**巨大な旧ソルバ丸ごとコピー、if(legacy)の全ソースへの散乱、個別のmain複製は避ける。** モダンがデフォルトのまま、旧版相当の計算をアダプタで実装する。

実際のクラス・ファイル構造を確認して柔軟に決めてよいが、推奨の概念構成：

```text
include/fluid/
  GlottalMethodOptions.h          # enumと組合せ妥当性チェック
  LegacyGlottalMethods.h           # 旧方式の純粋計算・read-only geometry API
src/fluid/
  legacyChannelArea.cpp            # 旧断面相当 (左右対称条件で厳密一致テスト)
  legacyFlowCoefficients.cpp       # 旧Lg/Rvの近似
  legacyWallPressure.cpp           # 旧壁面圧力分布・剥離以降
src/forces/
  legacySurfacePressureLoad.cpp    # 旧material-grid(i,j)→力。左右符号に注意
src/simulation/
  ...                              # 既存runループは1本、最小の切替のみ
tools/kanaya_ablation/
  compare_runs.py                  # 全ランから同一定義の指標と図を生成
  run_matrix.py または run_matrix.sh # ケースとCPU割当を管理
tests/
  test_legacy_glottal_methods.*    # 固定形状の単体テスト
```

実装上必要なら `LegacyFlowPipeline` クラス1個とサブメソッドに集約してもよい。**各機構の実装は別関数**とし、main計算の切替可能性、C++単体テスト、独立形状比較を維持。

### 2.1 設定インターフェース

既存`param.txt`の順序付きパラメータを壊さず、末尾の**名前付きオプション**で追加（`readParams.cpp`に既存のoptional/namedパーサあり）。名称の例（名前は既存命名規則に合わせてもよい）：

```ini
# defaults: all modern, no behavior change
glottalMethodProfile = modern
legacyArea = 0
legacyFlowCoefficients = 0
legacyWallViscosity = 0
legacySeparation = 0
legacyPostSeparation = 0
legacySurfacePressureLoad = 0
# optional future controls; separate stages
legacyContact = 0
legacyIntegrator = 0
```

- プロファイル`modern`: 全0。**未指定入力でも今のmainと数値結果が変わらないことが絶対条件**。
- `legacy_fluid_all`: 面積、Lg/Rv、壁面の粘性式、剥離・下流圧力、表面圧力荷重を**一貫した旧方式の組で適用**。構造積分Newmark、左右独立構造、現行接触、ζは変えない。「旧コード完全同一」ではなく**旧流体機構相当**と明示。
- `custom`: 各フラグで一要素だけ差し替える。**式の依存・添字・次元・初期条件が不整合になる組合せはエラーで拒否**し、必要なら比較単位をペアに統合。実際に5つの切替を全て完全独立にできなくてもよいが、理由をmanifestに記載。
- `legacyContact`／`legacyIntegrator`は**第2段階**。最初の主比較には入れず、必要なときに別コミット/別ケースで検証。無効なフラグを黙って無視しない（未実装なら入力を拒否）。
- すべてのランのmanifestへ有効なフラグと旧ソース出典を完全記録する。

## 3. 旧方式で再現すべき「物理式」と境界（旧成功版優先）

### 3.1 声門断面 A(x)

旧`calcArea.f90`では左右対称を仮定し、x方向の**物質格子行i**に対し1枚の声帯と中央平面との距離を2倍する：
\[
A_i^{old}=\sum_j 2\left[y_{mid}-\frac{y_{i,j}+y_{i,j+1}}{2}\right]\,|z_{i,j+1}-z_{i,j}| .
\]
元コードは正負gapを原則そのまま積分し、現行は左右実際の切断面上で`max(h,0)`を積分する。両者の閉鎖/負面積条件の違いに注意。**左右独立を保った実装**として、対称時に上式と一致する左右ペア幅`y_R−y_L`によるmaterial-row面積を作る。ただし非対称時は左右rowが物理xで異なり得るので**旧方式は比較用・対称に限定**するか、物理的な対応条件を明示して安全に使用する。現行x=const切断方式は保存。
- 50断面に旧面積を無理やり補間して「旧方式一致」とは言わない。`nsurfl/nxsup`と添字変換を明示。
- 初期、固定小振幅変形、正弦の既知位相、接触の各例で旧式とC++が一致するかテスト。

### 3.2 流量を決める Lg と Rv

公開`calcFlow.f90`では、最小面積`A_min`、声門の幾何長さ`lg`、流路代表長さ`lis`から
\[
L_{g1}^{old}=\rho\,0.5\,lis/A_{\min},\quad
R_{v1}^{old}=12\mu\,lg^2\,lis/A_{\min}^3 .
\]
`R_{k1}=\beta\rho/A_{\min}^2`は同一係数/パラメータなら不変と考えるが、実際の旧成功版・単位と確認する。現行`flowModel.cpp`は`0.5ρ∫dx/A(x)`と`12μ∫dx/G_3(x)`を採用。旧`lis`は固定形状等での定義をsourceから確認し、`xSep−xIn`と勝手に同一視しない。`lg`はAP長、単位m/mm換算を検証。
- 旧`Lg/Rv`と新`Lg/Rv`を**同じA_min,同じ凍結形状**で計算し、値・比・差・流量NR残差をCSVに残す。
- **壁面圧力勾配の粘性式と流量回路のRvは別の切替**。片方だけ旧方式にしたときの意味・数値整合性を明示。
- `Rk`、入口・声門下管・声道の回路`L/C/R`、β(Ps)、NR収束許容値、初期`U_g`も検証する。旧成功版が声道ありで新ランが声道なしなら**別の発振条件であり「手法差」ではない**。

### 3.3 剥離判定、psurf、圧力下流処理

旧`calcForce.f90`には、圧力更新中に`psurf(i)>psurf(i-1)`となれば`exit`する判定がある。添字iは1始まり。旧実装の「break以降に残る初期化ゼロ」を、現行の`argmin(A)`剥離や`areaRatio=1`、下流`Pd[0]`への遷移、0.5 mm smoothstepと区別する。
- `legacySeparation`: 旧ソースのbreak/セクション検索をその**検証済みの実際の定義**で再現。旧実装の`nsep`/出口挙動に曖昧さがあれば観測して文書化し、不定値のコピーは禁止。
- `legacyPostSeparation`: 旧版の壁面圧力割当（典型例：未計算行ゼロ）を別に再現。下流音響状態`Pd[0]`や`Pout`と壁面荷重の圧力を混同しない。**壁圧を0とするのか音響Pdとするかは出典に合わせ、恣意的に選ばない。**
- 現行`psurf`は0.5mmのsmoothstepで`Psep→Pd0`補間。旧`psurf`出力も統一の物理位置x・material-rowを記録し、補間効果を比較可能にする。
- 声門流量のNRで使う流体係数の積分終端と壁面の剥離終端が一致するか、旧／新各プロファイルで明示。旧では定義が異なればそのまま独立に比較。

### 3.4 圧力→表面荷重写像（流体仕事に直結）

旧`calcForce.f90`は物質格子`(i,j)`への`psurf(i)`直接割当＋`degree(1/2)`に依存する表面力：
\[
\boldsymbol f_{ij}=p_i\,\boldsymbol n_{ij}\,dS_{ij}
\]
現行は**変形後の物理x**から`p(x)`を補間し、局所接線の外積による面積ベクトルで`f=p\,d\boldsymbol A`を構成。
- 旧式の`cos/sin(degree)`、`ds*dz`の式と既存sourceのi/j境界を**明示的にC++化**。例：旧式`fx=+p ds dz cosθ_z sinθ_x`、`fy=−p ds dz cosθ_z cosθ_x`、`fz=+p ds dz sinθ_z`（元片側の向き）。左右への符号反転・正負圧の作用向きは旧単側と左右鏡映に整合させる。
- 旧Fortranの`calcForce.f90`には、圧力勾配ループ中の`j`参照などバージョン特有の危険箇所がある。**未定義変数・添字外・NaN・競合を「忠実再現」しない**。真に発振成功した版での意図を確認し、曖昧な箇所はコメント・テストで根拠を残す。
- 荷重の総面積、開大向き合力、重心位置、モード一般化力`Q_m=Σ \boldsymbol f\cdot\phi_m`、仕事`ΣQ_mΔq_m`を出力し比較。左右対称時`A_L=A_R`、`U_g`共有、y力左右逆向き、一般化力の鏡映が一致することを数値検査。

### 3.5 そのほかの差（主アブレーションからは固定）

- 旧RK4 vs 新Newmark(β=0.275625,γ=0.55) は別の軸。現行`rungeStep.cpp`の振動応答と旧`runge.f90`を**単独モード・無負荷・正弦外力・既知減衰**で比較。必要なら第2段階で`legacyIntegrator`追加。**RK4と流体写像を同時に変えない。**
- 旧接触は単側ペナルティ／中心面、新方式は左右独立面の接触。`contactReferenceFrequencyHz`の現行40 Hzと旧版の`omg(1)`の関係（旧ソースでは係数`kc1*omg²`）を照合。元成功時の接触値と初期gap/初期接触数を確認。第2段階まで新接触を維持。
- **圧力立上がり**：現行`flowModel.cpp`は省略時0.05秒のcosine ramp、公開`calcFlow.f90`は`Ps`を直接回路に投入。初期音響状態も一致するようにし、必要なら`pressureRampTimeSec`や別の明示的ramp modeを導入するが、**初期立上がりと流体機構を同時に変更しない**。
- **声道なし特例**：現行`flowModel.cpp`は`L_vt=0`で`La_=0.1,Ca_=1e30`の代替、`Pd,Ud`に特別な初期化がある。旧成功版の声道条件と一致させる。これが異なるならまず境界条件を揃えない限り差分の帰属は禁止。
- **構造モード**：現在のモード形状mass normalizationは各モード単独で、Gramの非対角がゼロとは限らない。旧成功実行のモード・周波数・質量と**対応順序**を確認。20モードで100モードと変形が変わることは比較条件の意図した差であり、旧方式効果と混同しない。
- **計算時間**：`nstep=20000`の最終1ステップで`State::mode2uf(step>=nSteps)`がスキップされる既知境界問題を把握。最終区間を無自覚に仕事解析へ入れず、できれば**このタスクの物理比較と独立に**小さな回帰テストで境界条件を確認。修正が必要なら既存のmodern再現性を再検証する。

## 4. 先に「固定形状・固定流量」で各モジュールの正しさを検証

本番の計算の前に、ソルバを時間発展させず**固定形状を使う比較ハーネス**を作る。既存の旧成功ランの時刻スナップショットが利用可能ならそれを優先。無ければ初期、均一gap、小変形、上下で異なるgap、x変形あり、部分接触の人工形状を使う。

段階と比較:
1. `A_old` vs 旧成功出力／旧式解析解。完全左右対称では同じ材料行のギャップ2倍の関係が成立するか。`A_modern`は厳密に同値とは限らない（切断方法が異なる）。
2. 同一`A_min`と`U_g,Pu,Pd`を凍結し、`Lg/Rv/Rk`旧/新、NR残差を比較。**完全連成で出た異なる`U_g`を入力して比較しない**。
3. 同一の気道状態と形状から旧/新`psurf`。圧力回復・壁圧のゼロ域・剥離x・負圧負荷の差を数値化。
4. **同一`p(x)`を与えるテスト**と「各方式固有の`p`」のテストを分離し、圧力値の差と空間写像差を区別。
5. 旧／新`Q_m`と`∑Q_mΔq_m`を比較。**時刻整合：Q_n×(q_{n+1}−q_n)**、右鏡映、接触と流体を別記録。
6. 旧実装で保証される量には強いassertを設定。旧成功バイナリの内部精度の違い、印字桁丸めによる誤差、異なる幾何定義の差は区別。

このハーネスで**旧手法と数式上一致すると分かってから**、初めてフルの連成ランを投入。比較CSV案：
`snapshot_comparison.csv` (source_hash,step,time,mechanism,side,x_i,Aold,Anew,Lgold,Lgnew,Rvold,Rvnew,Pg,Psep,xSep,Fx,Fy,Fz,Qm,diff,relative_diff);
量が異質なら正規化されたlong formの別CSVへ分割可。

## 5. 連成比較ケース：最初にR0/R1、次に識別力のある少数ケース

**共通条件**：
- 左右同一の「実際に振動した」材質・幾何・周波数ファイル、20モード、Ps1600Pa、ζL=ζR=0.06、初期gap、音響回路、圧力ランプ、接触、時間刻み、`nstep=20000`固定。
- ※ユーザーが「旧成功条件と他は同じ」と言っているので、**現行値・成功版値の照合を必ず実施**。違いが判明したら実行前にレポートに記録し、勝手に成功版条件へ書き換えない。
- 1ケースごとに独立runディレクトリ＆完全`manifest.json`。

| case | 使用方式 | 一つの問い |
|---|---|---|
| R0 | modernすべて（既定） | 現行20モード・1600 Paの振動状態を確定する |
| R1 | legacy_fluid_all、構造／接触／減衰はR0と同じ | 旧**流体処理一式**に戻すと発振するか |
| R2 | R1に新`A(x)`（互換性ある組で） | 断面評価方法に感度があるか |
| R3 | R1に新`Lg/Rv`のみ | 流量の慣性／粘性係数に感度があるか |
| R4 | R1に新`wall viscous pressure`／剥離規則（依存性あれば組で） | 壁面圧力形成に感度があるか |
| R5 | R1に新`post-separation pressure`のみ | 下流の壁面圧力の取り扱いに感度があるか |
| R6 | R1に新`surface pressure mapping`のみ | material-rowとx補間が励振を変えるか |

**順序**：
1. R0とR1を実行し結果と比較。R1が発振しなくても実装不良と即断しない。旧成功版との違い（構造積分、接触、ランプ、音響、モード）を固定形状＋小テストで先に検証する。
2. **R1で発振・R0で不十分なら** R2–R6を実行、各1因子差。ただし互換性のない組は棄却し、明示的に「組み合わせでしか変更できない」扱い。
3. R1/R0が共に発振する場合は、**振幅・平均面積・CQ・モード仕事などの差**で判断し、単純なyes/noだけで満足しない。
4. 単独差分で主要候補が見つかったら、上位2要因の組合せを**最大2ケース**実施。全2^Nの網羅実験はしない。
5. R0・R1がともに不十分なら、追加の大量R2–R6は保留し、**旧成功版との差、ζ、接触・ランプ・音響・積分**を切り分ける少数ケースを優先。
6. **ζ0.1／0.06の比較は単独因子**。同一流体方式・同一接触・同一境界で行い、過去の成功版ζ0.1という事实の影響を独立確認。時間積分や旧接触も同様に後段の別要因。

### 発振の判定

「最後に少し波打つ」「スペクトルに133Hzがある」だけでは発振成立としない。圧力ランプ終了後の同じ長さの窓、例 0.10–0.15s と 0.15–0.20s（±端点ずれを統一）を比較する：

- `Amin`平均、AC RMS、最小、閉鎖割合（例`Amin ≤ 0.01 mm²`での離散割合、閾値の妥当性明記）、CQ（閉鎖が周期的に成立する場合だけ）、固定材料点左右y変位・x変位の平均・AC RMS・成長率。
- `U_g(t)`、`P_g(t)`、`p_{\rm sep}(t)`、`x_{\rm sep}(t)`、`Lg/Rv/Rk`、正負圧の作用面積と左右各合力。
- モードごとの`Q_n`、`q_n`、`Δq_n`、仕事`∑Q_nΔq_n`、平均力／変動力の分解、減衰。**合計流体仕事の負だけで自励なしと決めない**。フェーズが不安定なら仕事の時窓依存性を明記。
- 同じ最終窓で振動が持続・成長しているか、過渡減衰か、一定開大に滞在しているかを分類。「閉鎖しない」と「励振がない」は別々に報告。
- 左右対称ランの鏡映誤差を強制検査（面積・流量共有、変位・圧力力の符号、一般化力、contact含む）。差が有意なら結果の物理解釈前に修正。

**判断基準**：旧流体一式へ戻しても発振しないなら「流体変更だけが原因」という仮説は弱まり、接触・積分器・圧力ランプ・音響回路などの差を検証。R1は発振してR0が発振しなければ流体差が強い要因だが、**その中のどの項目かはR2–R6で特定**。各変更加算は非線形系なので相互作用を認める。

## 6. CPUを全力で使い、**実際の壁時計時間を最小化**

単に`OMP_NUM_THREADS=$(nproc)`を固定するのではなく、利用可能な全CPUを活用する戦略を実測で選ぶ。

1. ローカルWindows+WSL2の`nproc`、`lscpu -e=CPU,CORE,SOCKET,ONLINE`、`lscpu | grep -E 'CPU\(s\)|Core\(s\)|Thread\(s\)'`、`grep Cpus_allowed_list /proc/self/status` を取得。WSL2の割当制限・他ジョブを確認。
2. まず小規模な**各同条件（20モード, 2000ステップなど）**で1/4/物理コア数/利用可能論理CPU数のスレッド候補をテスト。runごとの結果SHA、計算時間、CPU時間、peak RSS、CPU使用率を保存。診断版の多数ファイル出力がボトルネックなら軽量比較ログへ切替。
3. OpenMPの起動例：
   ```bash
   export OMP_NUM_THREADS="$(nproc)"
   export OMP_PROC_BIND=close
   export OMP_PLACES=cores
   export OPENBLAS_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export VECLIB_MAXIMUM_THREADS=1
   export NUMEXPR_NUM_THREADS=1
   /usr/bin/time -v ./simulation /absolute/path/to/case_param.txt
   ```
   **ただしこのスレッド数が最速と検証されなければ最速候補を採用**。設定値・ピークCPU使用率をrun manifestへ保存。WSL2で`OMP_PLACES=cores`が有効かログで確認。
4. R0/R1や複数アブレーションを**並行投入**する場合、ジョブごとのOpenMP数を分割し、合計使用スレッドとCPU affinityが許容CPU資源を超えないようにする。例：利用可能16コアなら、16-thread×1または8-thread×2の実測で短い方。**16-thread×全ジョブの過剰並列は禁止**。
5. `icpx`＋Release＋OpenMPを優先。`-ffast-math`は禁止、浮動小数点削減順序を勝手に変えない。`schedule(static)`、数値再現性の維持を優先。ビルド設定・CPU機種・スレッド数・実行時間を記録。
6. 高負荷長時間実行はユーザーPC上を前提とし、他者共有サーバ・授業計算資源では占有許可を確認。**全コア使用の目的は短時間終了であり、数字としてCPU 100%を追うことではない**。

## 7. テストと既存modernの後方互換性

**停止ゲート（テスト失敗時は本番大量計算禁止）**
- オプション未指定と`glottalMethodProfile=modern`の**bitwise一致**を基本目標。開発差により不可避なら同一コンパイラ同一スレッドで、複数時刻の全観測量が高精度に一致し、差の理由を提示。旧フラグを入れただけでmodernが変わるなら原則不可。
- `old_area`: 同一の左右鏡映・material-row平面・既知gapでFortran式に一致。左右反転で結果不変、異なるAP幅・部分接触・端点処理は期待動作を明示。
- `old_coeff`: 単位をm統一し、一定`A`・`G3`の人工例で解析解、各種Aで`A^-1,A^-3`のスケーリングを確認。
- `old_pressure`: 同じ形状・`Ug,Pg`に対し旧成功版の`psurf`・break位置と一致。下流wall pressureと下流音響`Pd`を混同しない。
- `old_force`: 直線、傾斜、鏡映面で表面積とFx/Fy/Fz、モード一般化力が一致。左右に同方向の開大符号間違いがない。
- 切替オプションのparser、無効値、不整合モード、既存入力形式の互換性。ディレクトリ別出力と別runの混入防止。
- 連成プロファイル`modern`のすべてのforce/coefficient/work/areaの不変性テスト。20モード20k-stepのR0は既存modernの正解として保存する。
- コンパイル`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`、`cmake --build build -j`。`ctest`または対象テストを実行。ビルドはインテル`icpx`構成を既存環境に合わせ、勝手に別の物理結果へ差し替えない。

## 8. 成果物と報告

追加コード、CMake登録、設定parser、独立した固定形状比較、run matrix、比較解析を含める。原則として`src`へ旧方式の本体を隔離する。

```text
docs/plans/KANAYA_LEGACY_SWITCH_ABLATION_20MODES_1600PA.md  # 本書
src/fluid/legacy*.cpp           # 旧断面・流量係数・圧力
src/forces/legacy*.cpp          # 旧壁面荷重
include/...                     # オプション・宣言
tools/kanaya_ablation/...       # 実行制御と比較可視化
tests/...                       # 合成テストと既存modern回帰
output/kanaya_ablation/         # 各case独立・既存outputを変更しない
  R0_modern/
  R1_legacy_fluid_all/
  R2_.../
  ...
  comparison_summary.csv
  comparison_report.md
  comparisons.png
  run_matrix.json
  provenance.json
```

`comparison_report.md`に必ず以下を入れる：
1. 比較に用いた「成功実績ある旧コード」のバージョン・実際の入力と、公開Fortranの関係。参照不能なら明記。
2. どのサブ機構の旧式を忠実に再現したか（対応する元ソースと式）、互換でない組合せと制限。
3. 固定形状比較での絶対・相対差、旧再現テスト合否、modernデフォルトの回帰誤差。
4. ケース別に`平均Amin, AC RMS, 最小Amin, CQ, f0, 位相/周期性, Ug, 圧力/荷重, モード仕事, 速度RMS`等を表と図で比較。
5. **「最も発振を弱くした候補」・「根拠」・「分離できなかった差」・「次に変更すべき条件（最大2）」**。旧方式が発振したことを「旧式が物理的に正しい証拠」と誤認しない。
6. OpenMPスレッド数、CPU利用可能コア数、1/4/物理/論理の測定結果、ジョブ並列数、計算所要時間。
7. 原本の保護（元`/home/.../simulation-main`の変更なし）、既存成果物の上書きなし、Git status、作業branch、出力絶対パス。
8. 原因が未確定なら未確定と記す。途中までなら**何が実装済み／検証済み／未実行かを正確に区別**。

### Codexに要求する完了条件

- **旧流体計算機構を別`src`で実装し、主ループから設定で選択できる。**
- **旧／新の固定形状テストが通る。** 旧成功版がアクセス可能ならそこに対する数値同等性を検証。
- **現行modern入力は無変更・挙動不変**、左右対称の数値対称性が維持される。
- **20モード・1600PaのR0/R1比較を実行**、条件とCPUを記録、必要に応じR2–R6を追加。未完了なら理由を報告。
- **最速の実測CPU構成で実行し、出力を分離、結果・原因候補を定量的に整理**。
- **このタスクで旧モデルを主ソルバのデフォルトにしない。物理モデル改善は原因が判明した後の別タスク。**
