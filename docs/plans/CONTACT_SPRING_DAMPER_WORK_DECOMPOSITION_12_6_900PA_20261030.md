# 12/6 kPa非対称声帯：接触ばね・速度依存項の仕事分解診断（Codex実装仕様）

> 2026-10-30 / 対象ブランチ：`main`（文書の配置先）  
> **実装・実行は未着手。本書はCodexがローカル作業ツリーで実施するための指示書。**  
> **今回の追加本番計算は Body 12/6 kPa・Cover 3/3 kPa・Ps=900 Pa の1条件を最大1回。** 他の剛性・圧力・接触係数の探索、12/12の再計算は禁止。

## 0. 目的と研究上の位置づけ

これまでの3Dモード重ね合わせ計算では、左右のBody剛性差によって、柔らかい右声帯のほうが小振幅になる条件を確認した。12/12対称条件と12/6非対称条件の比較から、右のモード応答・観測点でのモード相殺・接触参考仕事の変化が明らかになった。

**今回の問いは「非対称12/6条件で観測された右側の大きく負の接触参考仕事が、接触ばね反力と速度依存力のどちらから生じているか」である。**

次の量を、**ソルバが実際に採用した接触荷重**について分ける：

1. ばね項（`kc1`, `kc2`）による接触仕事 `W_contact,spring,ref`
2. 閉鎖速度に比例する項（`kc3`）による接触仕事 `W_contact,damp,ref`
3. 両者の和 `W_contact,total,ref` と、前回診断の `W_contact_reference` の整合
4. 左右別・両側合計・主な構造モード別・AP端/中央別・接触時刻別の分布

**解釈上の線引き**：これはまず「採用された接触荷重による離散参照仕事」の分解であり、`W_damp<0` が出てもその絶対値を厳密な不可逆散逸と断定しない。接触ばね項は位置依存でも、接触面積・法線・対応点が時間変化し、さらに荷重が予測接触反復由来なので、単純に保存ポテンシャルの変化と等置できない。1条件の受動観測では「接触が右の振幅低下を引き起こした」という因果までは立証しない。

## 1. 完了済み診断の固定基準（新しい結果と混ぜない）

**報告された事実**（2026-10-29の`report(3).md`、`report(4).md`）：

- 非対称：左Body 12 kPa／右Body 6 kPa、左右Cover 3/3 kPa、Ps=900 Pa、ランプ0→900 Pa/0.10 s、`dt=1e-5 s`、`nstep=50,000`、左右20モード、`zetaL=zetaR=0.06`、`kc1=5, kc2=6000, kc3=5000`、最大接触反復100、初期間隙0 mm。
- 評価主窓：`[0.30,0.50) s`。代表点のhalf-ptp R/L=**0.302355**、表面half-ptp空間RMS R/L=**0.4906**。
- 12/6接触参考仕事：**L=−0.152632 mJ、R=−1.20558 mJ、左右合計=−1.358212 mJ**（報告数値の桁数での和）。**今回も同一窓・同一離散定義でこの値を回帰基準にする**。
- 同じ窓のraw流体仕事：L=2.93198 mJ、R=1.83817 mJ。**raw流体仕事と動的流体仕事は別であり、今回は流体仕事の再同定をしない。**
- 対称12/12の接触参考仕事：L=R=−0.95918 mJ、両側合計=−1.91836 mJ。ただし**対称runにはばね／速度依存項の別々の採用時系列がない**。既存の対称値を比較の背景とするにとどめ、対称条件のばね/速度依存内訳を捏造しない。
- 元run／診断runは全9個の`dat`がSHA256一致、診断ON/OFFはbit一致、全99,998 side-interval/各runの仕事整合性は合格。各runの保存250 VTUで四面体反転0。
- `12/6`の前回診断では「ばね/速度依存項は未保存のため分離できない」。**このため今回のみ新しい1回の12/6実行が必要になり得る**。

探索すべきローカルの元データ（**名前は既報、存在や配置を開始時に再確認**）：

- `analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900`
- `amplitude_asymmetry/production` と `amplitude_asymmetry/` 以下の既存`report.md`、`manifest.json`、`checks.json`、`modal_work.csv`、`contact_geometry.csv`
- 既存12/6診断`production/kanaya_diagnostic/modal_native.bin`、`run_metadata.json`、`contact_regions.csv`、`consistency.csv`。
- 対称比較の`compare_symmetric_asymmetric.py`／`report.md`があれば参照。**GitHub mainにこれらの最新未コミットコードや大容量runファイルが必ず存在するとは仮定しない**。

実行条件の正は元runの`txt/params_used.txt`、`txt/manifest.txt`、入力ファイルハッシュ、既存計算ログ、ローカル実装で確認する。`input/param.txt`を勝手に使わない。

## 2. 安全・再現性の絶対条件

1. 開始前に`pwd`、`git status --short`、`git branch --show-current`、`git rev-parse HEAD`、`git remote -v`を報告し、現在の`AGENTS.md`を読む。**未コミットの差分・既存runを絶対に消さない／上書きしない**。`git reset`／`checkout`／`rebase`／`clean`／無断commit・pushは禁止。必要なら独立したworktree/作業ブランチ/ビルド先を使用。
2. ソルバの**接触圧の算式、流体、Newmark、モード投影、時間刻み、係数、接触反復の判定と順番、OpenMP順序、既定出力を変更しない**。追加するのは診断observerと後処理のみ。
3. 診断は**既定OFF**。既存`ENABLE_KANAYA_DIAGNOSTICS`(CMake)＋`SIMULATION_KANAYA_DIAGNOSTICS=1`(実行時)を利用し、必要なら追加の`SIMULATION_CONTACT_COMPONENTS=1`等の**明示opt-in**を設ける。OFFでは新しい出力も計測コストも発生させない。
4. **本番フルランは0.50 sの12/6を最大1回**。最初に既存データの棚卸し→小さな単体/スモーク検査→基準差分ゼロ確認→必要なときだけ本番を実行。対称12/12の再計算、係数`kc3=0`などの物理アブレーションは別タスク。
5. 新規runは`SIMULATION_RUN_DIR`で固有の空ディレクトリへ。元`output/`、`result/`、`analysis_runs/`、前回の`production/`には書き込まない。VTU/VTU入力は読み取り専用として利用。出力保存・ビルドに十分なディスク空きを先に確認。
6. 元runのhashと新runのhash、診断コンパイルオプション・入力ファイルSHA256・ビルド情報・スレッド環境・最終採用反復、評価窓の定義を`manifest.json`で記録。
7. 作業はどこまでできたかを報告。数値不一致・入力欠損でゲートが通らなければ**重いフルランを止め、失敗原因を明記**。成功したことにしない。

## 3. 現行ソースで先に確認すべき実装上の事実

2026-10-30にGitHub`main`で確認できる主要部位（**ローカルの現行内容と差分確認が必須**）：

- `src/forces/forceCalculatorContact.cpp`：`ForceCalculator::applyContactLoads`、`computeXYNormal`、`processPair`、`addSegmentForceXY`、`prevFdisX/Y(L/R)`。
- `src/simulation/simulation.cpp`：`resetPreviousContactForce()`→接触反復（モード積分・表面予測変位→`applyContactLoads`）→**最後に適用した接触荷重で再投影・Newmark最終積分**→`kanaya.adoptedContact()`→`applyContactLoads(...,diagnosticOnly=true)`→`kanaya.endStep()`→状態確定`uf2u()`。
- `include/diagnostics/KanayaDiagnostics.h`、`src/diagnostics/KanayaDiagnostics.cpp`：`beginStep`、`afterFluid`、`beginContact`、`contact`、`adoptedContact`、`endStep`、`finish`。現行の`contact(j,x,z,pen,FxL,FyL,FxR,FyR,excluded)`は**総接触力しか受け取らず、成分別の節点割当を持たない**。
- `tools/amplitude_asymmetry/analyze_production.py`、比較器`tools/amplitude_asymmetry/compare_symmetric_asymmetric.py`（ローカルに存在すると報告済み）、既存テスト`tests/test_kanaya_diagnostics.cpp`。
- 現行`modal_native.bin`は`[step,side,mode,7]`、フィールドは`Q_pos,Q_neg,Q_fluid,Q_contact_adopted,q_n,qdot_n,delta_q`。**既存の7列schemaを壊さず、新成分は別ファイルにする**。

### 3.1 接触圧はソルバの**既存値をその場で採取**

現在の接触ペアで、単位mに換算した法線方向貫入量`δ≥0`、閉鎖側の相対法線速度から`δdot_+=max(0,-v_rel,n)`を計算する。実装の形は：

\[
p_{\rm spring}
 =kc1\,\omega_c^2\,\delta\,\bigl(1+kc2\,\omega_c^2\delta^2\bigr),
\qquad
p_{\rm damp}=kc3\,\dot\delta_+,
\]
\[
p_{\rm contact}=\max(0,p_{\rm spring}+p_{\rm damp}),
\qquad \omega_c=2\pi\,f_{\rm contact,ref}.
\]

- ここで`p_spring`/`p_damp`は現行コードの局所変数`f_spring`/`f_damp`に相当する**圧力[Pa]**であり、まだ節点力[N]ではない。`kc`の単位・実係数・速度の評価時刻をコードから確認して記録。
- `f_damp`は、単純な現在/次時刻変位差ではなく、接触ルーチンがその反復で参照した`stateL.vel/stateR.vel`を補間した値による。**ソルバの採用したその値を採取**し、後処理で時刻の違う`qdot`から作り直さない。
- `normalGap`と相対速度は同じ接触法線の定義だが、接触法線・面積・対応点・閉鎖枝は時変。`δdot_+`を`dδ/dt`と常に同一視しない。
- 合力：`areaL`、`areaR`を重なり長さに応じて計算し、`area=0.5*(areaL+areaR)`、`F_contact=p_contact*area`、`FxL/FyL=−F*(nx,ny)`、`FxR/FyR=+F*(nx,ny)`。各側`addSegmentForceXY`で隣接2節点へ重み`(1−sL,sL)`・`(1−sR,sR)`で分配している。
- 「総圧の非負clamp」「非有限値除外」「`pen<=0`」「`area<=0`」「`contactIsExcluded`」の各枝を通った**実際の適用力**だけを成分分解。`p_spring+p_damp=p_contact`が成立しない枝があればその理由を記録し、強制的に等式を捏造しない。今回の正の`kc`では通常どちらも非負のはずだがテストする。

### 3.2 **採用接触反復**と診断専用呼び出しを混同しない

現行`KanayaDiagnostics::beginContact(diagnosticOnly,iteration)`は、`diagnosticOnly=false`の各反復で地域集計をリセットし、最後に実行された反復を保持する。`adoptedContact()`は最終的に`prevFdisXL/YL/XR/YR`へ入った**適用済み接触荷重**をモード投影する。

**実装の優先策：**

1. 各`applyContactLoads(...,diagnosticOnly=false)`開始時に、observer専用の左右の**springとdamp節点力バッファ**をゼロにする（小さな`nsurfl×nsurfz`のX/Y成分を各側×2種）。
2. `processPair`内で、既存の`f_spring`、`f_damp`、`area`、法線、`iL/iR`、`j`、`sL/sR`をobserverへ**読み取り専用で渡す**。observer側で**同じ節点への分配規則**を使い、ばね力と速度依存力を**別々のバッファ**に積む。**既存の`contact_pressure`、`F`、`fdis*`への演算は変えない。**
3. 追加情報には`j,xContact,z,pen,penDot,pressureSpring,pressureDamp,normalX,normalY,area,excluded,contactIter`、左右のsegment節点indexと補間重みを含める。**必要最小限の引数で読み取る**。新callbackは`#ifdef ENABLE_KANAYA_DIAGNOSTICS`でのみビルドし、runtime無効時はno-opにする。
4. `excluded=true`や`diagnosticOnly=true`では、**適用されていない荷重を成分バッファに加算しない**。接触ペア統計（候補/除外）は任意の別カウンタに分ける。
5. 各反復開始でバッファをリセットするため、次の反復の候補に加算された仕事を累積せず、最後の反復のバッファだけが残る。**診断専用の最終`applyContactLoads(...,true)`でリセットしない・上書きしない。**
6. `kanaya.adoptedContact()`でspring/dampバッファの和が、実際の`prevFdis*`（最終採用接触）と節点別・成分別に一致することをまず検算する。**これに合格した場合だけ**各モードへ個別投影する。
7. `kanaya.endStep()`で、最終採用の総`Q_contact`と`q_{n+1}-q_n`に合わせて各成分の仕事を記録する。既存`Q_contact_adopted`/`W_contact_reference`との**同じ段階・同じ`Δq`**を使う。
8. `ncont`で打ち切り、`no_contact`で終了、反復で接触が消える場合、連続接触、複数接触ペア、非単調形状のfallback、`diagnosticOnly`を全て扱う。接触がなければ両成分をゼロとし、欠測`NaN`とは区別する。

**注意**：既存`applyContactLoads`は、接触荷重更新時に`fx/fy += (fdis−prevFdis)`と差分加算する。ここで記録すべきは**反復間の差分`dFx/dFy`ではなく、最終採用の全接触荷重`prevFdis`**。差分を積算して重複計数しないこと。

## 4. 接触成分の仕事と代数的な整合性

左右`s∈{L,R}`、20モード`m`、採用接触力の成分`a∈{spring,damp}`について：

\[
Q^{a}_{s,m,n}
 =\sum_{k\in\mathrm{surface}}
   \mathbf f^{a}_{s,k,n}\cdot\boldsymbol\phi_{s,m,k},
\]
\[
\Delta W^{a}_{s,m,n}
 =Q^{a}_{s,m,n}\Delta q_{s,m,n},
\qquad
\Delta q=q_{n+1}-q_n.
\]

**流体荷重の`Q_fluid`・総接触`Q_contact`と同じ投影方法**を使用し、`state.q`（確定前`t_n`）と最終`state.qf`（`t_{n+1}`）を組み合わせる。`qdot_n`との単純乗算へ置き換えない。

各`n`で必須の検算：

1. **圧力**：適用した各接触ペアで`p_spring+p_damp=p_contact`（clampで例外があれば明記）。両方の圧力は有限・非負。
2. **節点力**：各左右節点・X/Y成分で`f_spring+f_damp=prevFdis`。Z成分の接触力は現行実装では使用していない。`FxL+FxR=0`、`FyL+FyR=0`は**ペア総和レベル**の作用反作用として検算（左右の節点分配先は異なってよい）。
3. **モード力**：各`m`で`Q_spring+Q_damp=Q_contact_adopted`。また`Q_fluid+Q_contact_adopted=Q_total_applied`（既存`fi`）を検算。
4. **仕事**：各側・各モード・各区間で`W_spring+W_damp=W_contact_reference`。全モード和・時間窓和・左右和についても成立。
5. **節点とモードの仕事**：
\[
\sum_m Q^a_{s,m,n}\Delta q_{s,m,n}
 \stackrel{?}{=}
\sum_k\mathbf f^a_{s,k,n}\cdot\Delta\mathbf u_{s,k,n}^{SI},
\]
実際の`State::mode2uf`の再構成・`mm→m`換算と同じ節点集合で検算。**単位未確認なら物理Jとして断言しない**。総接触仕事の既存検算に新たに加える。
6. **既存診断との互換**：`modal_native.bin`の7列、`modal_energy_summary.csv`、既存`consistency.csv`、本体`dat`の意味や既定フォーマットを変えない。

既存observerの`scaledResidual`（相対`1e-8`＋絶対床`1e-12`等）を参照し、**N・一般化力・Jに応じて絶対許容値を明示**。ゼロ近傍で相対誤差だけを使わない。**整合が通らない場合は物理的な正負の議論を止め、`FAILED_QC`とする**。

### 重要：速度依存項の負仕事を自動的に仮定しない

`p_damp≥0`であっても、
\[
W^{damp,ref}
 =\sum_{n,k}\mathbf f^{damp}_{n,k}\cdot
 \bigl(\mathbf u_{k,n+1}-\mathbf u_{k,n}\bigr)
\]
の符号は、**力を評価した予測反復と、最終Newmark解の時間的なずれ**、接触位置・法線・節点分配が変わることなどに依存する可能性がある。`W_damp>0`が現れても自動的に実装バグとしない。個々の採用荷重の評価時刻、相対法線速度の符号、適用反復と仕事の符号を確認する。

物理的な「散逸」の追加評価を**条件付き・参考**で行うなら、接触ペア`F_damp`とそのペア位置での採用反復の**相対法線速度`v_rel,n`**から瞬間候補`P_pair,damp=F_damp v_rel,n`を計算し、閉鎖枝で非正になるか検査。ただし、採用反復の予測速度に関する局所的な参考量であり、**`Σ P_pair Δt` を本体の有限増分仕事と無条件に同一視しない**。接触点が移動するため厳密な散逸率とは別扱いにする。

また`f_spring`には非線形項`kc2`が含まれる。さらに原因を絞る必要がある場合だけ`linear_spring`/`nonlinear_spring`の補助小計を採取してもよいが、今回の**必須成果物はspring合計とdampの2成分**である。

## 5. 解析窓・記録戦略（大量出力を避ける）

**最優先の解析窓**：`[0.30,0.50) s`。副窓：`[0.10,0.20)`、`[0.20,0.30)`、`[0.30,0.40)`、`[0.40,0.50)`。`dt=1e-5`に合わせ、実際に有効な整数step集合と積分区間を`manifest`へ書く。既存observerは`nstep-1`区間のみを仕事計算対象とする既知の終端除外がある。**同一ルールを維持**し、`0.50s`端を勝手に補完しない。

**小さい追加バイナリを推奨：**

- `contact_component_native.bin`：各有効ステップ、side L/R、active mode順に**2つのlittle-endian float64**（`Q_contact_spring_adopted`、`Q_contact_damp_adopted`）。`delta_q`は既存`modal_native.bin`を参照する。計測対象20モード/side、49,999区間なら目安約32 MB。**レコードschemaと正確なバイト数はメタ情報で検証**。
- 既存`modal_native.bin`を複製・再生成することは避けられないが、**既存バイナリの7フィールドには追記せず**、別ファイルにして過去の比較器と互換性を保つ。
- `contact_component_summary.csv`：終了時に解析窓×side×mode別（`W_spring,W_damp,W_total_ref`、各正負や事件/非事件の加算、使った区間数、品質フラグ）を記録。重複する既存集計があるなら後処理だけで生成してよい。
- `contact_component_regions.csv`（必要最小）：**AP端/中央の成分別節点合力**、接触pair数、`p_spring/p_damp`、penetration、`penDot`、`x,z`の区間平均・ピーク、採用反復番号。毎stepの重い全節点ダンプは不要。**AP端・中央の仕事**まで報告するなら、`Q Δq`に対応する各地域の節点力を成分別に保持し、`Σ_k f_{region,k}·Δu_k`を内部dtで累積する。地域ごとの合力×代表点変位を「地域別仕事」と呼ばない。
- `run_metadata.json`は旧フィールドの意味を保ったまま、新しい独立ファイル名、成分順、dtype/endian、shape、step除外規約、単位、接触評価段階、`spring+damp=total`の最大残差を**追加メタデータとして記録**。既存解析器が未知キーを許すか確認する。安全なら別`contact_component_metadata.json`でも可。

**効率**：observer内で毎stepの成分別`Q`と`W`を計算し、ペア別の全履歴・全節点×全stepのCSVは作らない。後処理で`numpy.memmap`またはchunked readを使う。不要な追加依存や`input`への書込みは禁止。

## 6. 既存12/6結果との同条件回帰・1回だけの本番

### A. 先に実行条件を完全照合

- 元12/6の`manifest`/`params_used`を読み、Body/Coverと対応する左右COMSOL固有周波数・VTU入力のハッシュ、geometry/NAS、その他パラメータ、音響回路、圧力履歴、モード選択、`contactReferenceFrequencyHz`、`AP除外フラグ`、コンパイルマクロを記録。
- 現行ソルバが元runと同じ数値軌道を再現できるか、Git SHA/dirty/コンパイラ/OpenMP/アルゴリズム差を評価。**`main`とローカルの違いを勝手に解消しない**。
- 前回12/6診断の`W_contact_reference`の窓集計値に`L≈−0.152632, R≈−1.20558 mJ`があることを検証。**元runの生バイナリのあるものを優先**。

### B. 追加コードの単体テスト（本番前）

- 接触0のケース：`spring=damp=total=0`。候補ペアがあっても`excluded`なら適用ゼロ。
- `kc3=0`を**単体の合成接触計測関数/小規模スモーク入力**でテスト：速度成分0、総力=ばね。**本番の接触係数は変えない**。
- `kc1=0`等の合成検査：ばね成分0、総力=速度依存項。**本番条件は不変**。
- 閉鎖速度<=0（closing）のときの`p_damp`、opening時の`p_damp=0`、penetration/areaのゼロ枝、非線形`kc2`、左右作用反作用、分配重み`s=0/0.5/1`、AP端/中央、複数ペア・重複節点、clamp/除外枝。
- **2回以上の接触反復**：1回目の力を上書きし最終採用のみを残す。`ncont`終了・`no_contact`終了・診断専用`diagnosticOnly=true`はobserverバッファを汚染しない。
- synthetic`Δq`で`Q_springΔq+Q_dampΔq=Q_totalΔq`、力分解/節点仕事一致/左右総和/時間窓集合が成立。
- 旧`modal_native.bin`解析器と`compare_symmetric_asymmetric.py`が、新しいobserverのファイル追加でも動く（旧schemaの互換性を守る）。
- 既存C++テストの`test_input_file_params`等に古いサンプル名の期待値がある場合は、既知の陳腐化として**本タスクと無関係な修正で範囲を広げない**。結果に記載。

### C. 80～200 step程度のスモーク（同一条件）

①診断をコンパイル時OFF、②コンパイルONだが新診断実行時OFF、③新診断ON、を**独立した空のrunディレクトリ**で実行。元・診断用の同一入力・OMP配置・コンパイラで比較する。

**ゲート**：`regression_state.bin`または既存比較器でstate・`dat`・接触反復回数・代表変位・流量/圧力・VTUがbit一致するか確認する。機械精度差があるなら`max_abs/max_rel`と原因を報告し、**bit一致を勝手に近似判定へ置き換えない**。新observerの各stepのspring+damp=accepted totalの全検査が合格してから本番。

### D. 本番は12/6だけ・最大1回

- フルラン`nstep=50000`、`dt=1e-5`、`nmode=20`、`Ps=900`、`zeta=0.06`、`kc=(5,6000,5000)`、ランプ0.1秒、初期間隙0、その他元条件に厳密一致。
- `ENABLE_KANAYA_DIAGNOSTICS=ON`＋新成分診断実行時ONを設定し、**本番run用に新しく確保した`SIMULATION_RUN_DIR`**を使用。保存先は`output`の既存最新結果を上書きしない。実行のために元`input/param.txt`を変更しない。
- 完走後、**今回runの全9`dat`と、元12/6 runの全9`dat`のSHA256を比較**。同一結果なら`trajectory_identical=true`。違えばsource/binary/configと波形を調査し、旧負仕事と同一軌道の分解として断言しない。**再実行で合わせ込まない**。
- 新runの`W_contact,total,ref`が前回L/Rと一致すること、`W_spring+W_damp=W_contact,total`、旧`modal_native.bin`の`Q_contact_adopted`が全区間/全モードで一致することを必ず確認（ビット一致か許容誤差かは項目ごとに明示）。
- `consistency.csv`の全side-interval検証合格、`D>=0`、有限値、採用反復対応、ファイルバイト数、非実行分ゼロを確認。保存VTUで四面体符号付き体積を再検証するなら既存`audit_vtu_geometry.py`を使い、検査した保存時刻のみ報告。
- 成功したら結果を解析。失敗したら`FAILED_QC`で終え、次のフルランは提案のみ（勝手に実施しない）。

## 7. 後処理で必ず出す比較・結果の読み方

入力は既存`kanaya_diagnostic/modal_native.bin`と今回追加`contact_component_native.bin`（同じrun由来）、前回12/6と比較用12/12の`report.md`/CSV。別runの`Q`と`Δq`を混ぜない。

1. 0.30–0.50秒の**左右それぞれ**`W_spring`、`W_damp`、`W_total_ref`、**左右合計**をmJで提示。`W_damp/W_total`は総和がゼロ付近なら`NA`。分母符号・正負の相殺を説明し、割合が100%を超えても勝手に丸めない。
2. 左右の第4・5・7・9・12・14モードと、`|W_contact_ref|`が大きい上位モードを比較。モードの固有周波数と元source mode番号を残す。**モード別仕事を「独立モードの物理エネルギー損失」と断言しない**（既報Gram非対角が0.06程度）。
3. 接触発生時刻、`max penetration`、`penDot`、接触圧（ばね/速度依存）、AP端/中央、`Amin`、代表点gap、左右変位を重ねる。**負の仕事が各周期に広く分布するか、少数の接触時刻に集中するか**を見る。周期分類が未確定なら`cycle`と断言せず時間窓で評価。
4. ばね項の仕事が大きく負の場合、接触面積/反力が大きい時点でどちらが変位したか、反復と位相に注意し「保存的な力でも有限窓の参考仕事は負になり得る」と明示。速度依存項の負仕事が主の場合も、「散逸候補の寄与」と表現して**絶対的な散逸エネルギーの確定値とは区別**。
5. **右だけでなく左との和**を必ず確認。作用反作用の仕事は各側で等しいとは限らない。**右の負値が大きい原因がspringとdampのどちらか**、左の負値低下がどちらに関係するかを説明する。
6. 12/12対称runについては**総接触参考仕事だけ**との背景比較に止める。対称runの`W_spring`/`W_damp`欄は`NOT_MEASURED`とし、「12/6の速度依存項が対称より大きい」とは主張しない。
7. 主報告では`(i) 観測した代数値`、`(ii) それに整合する機構候補`、`(iii) 未確定のエネルギー散逸・因果`を明確に分ける。

### 解釈の例（**仮説であり予想結果ではない**）

- **Hspring**：右の負の接触参考仕事の大部分がspring → 接触拘束・弾性エネルギー交換と整合。ただし有限窓のばね仕事は回復可能性もあるため、全量散逸とは言わない。
- **Hdamp**：右の負の仕事の大部分がdamp → 接触の閉鎖速度依存項が運動に負の仕事をすることと整合。ただし**これだけでは右の振幅が小さくなった原因とは断定できない**。
- **Hmixed**：spring/damp両者が重要、または相殺して総仕事が小さくなる → 成分別仕事と総和を併記し、無理に一方を「主因」にしない。
- **Hstage**：分解値と採用総荷重が一致しない → 収束/反復/診断専用呼び出し混入などの観測器問題を先に修正。機構解釈しない。

## 8. 成果物・ファイル構成

出力は固有の新run内の`kanaya_diagnostic/`と、その外側の**1つの解析フォルダ**にまとめる。物理ソルバからの通常出力は元runと比較するだけで編集しない。

| 生成物 | 必須内容 |
|---|---|
| `contact_component_native.bin` | 有効step×L/R×モード×spring/dampの一般化力。旧`modal_native.bin`と別 |
| `contact_component_metadata.json` | dtype、endian、shape、mode一覧、step除外、荷重の採用段階、単位、全検算閾値と最大残差 |
| `contact_component_summary.csv` | side/mode/window、ばね/速度項の仕事、総和、旧referenceとの差、品質 |
| `contact_component_regions.csv` | AP端/中央、採用反復、接触数、圧力・節点力の成分統計。可能なら地域別離散仕事 |
| `checks_contact_components.json` | OFF/ON、元run対照SHA、力/モード/仕事/節点一致、左右合計、QC、欠測 |
| `manifest.json` | 元run/新run identity、入力SHA、ソルバ/診断schema、窓、実行環境と不変更証跡 |
| `report.md` | まず3～5行の結論。左右・両側和、成分別主モード、AP/時間局在、整合性、限界、次の最小検証案 |
| `figures/` | 最大4枚：①左右の成分別仕事、②成分別累積仕事時間変化、③モード別仕事、④接触時刻/貫入/圧力/局所gap |

Python側の推奨場所：`tools/amplitude_asymmetry/analyze_contact_components.py`。既存`tools/amplitude_asymmetry/analyze_production.py`や`compare_symmetric_asymmetric.py`は**破壊的変更せず**利用できる集計を再利用する。

**CSVとレポートは作業runディレクトリの成果物であり、巨大な計算出力や生VTUをGitに追加しない**。追加コードと小さなテストだけを変更対象とする。実装完了時は`git diff --stat`と具体的な検証コマンドを報告。

## 9. 完了条件（Codexチェックリスト）

- [ ] 元12/6 runの実行条件・ハッシュ・前回の`W_contact_reference`（L/R）を再確認した。
- [ ] 既存`applyContactLoads`の採用荷重に対応するspring/dampの**節点別分解**ができた。反復/除外/`diagnosticOnly`を混入していない。
- [ ] `spring+damp=total` が **圧力、節点力、モード力、左右和、離散仕事**の全てで通った。
- [ ] 既存`modal_native.bin`と`dat`のschema・意味を壊していない。
- [ ] 単体テスト／短いOFF-ONスモークを実施し、**ソルバの数値結果が変わらない**ことを証明した。
- [ ] 12/6・900 Paの追加本番計算は最大**1回だけ**。出力を全て隔離し、旧runのSHA256と確認した。
- [ ] 0.30–0.50秒における**L/Rそれぞれと両側和**のspring/damp/totalをmJで報告。主モードと時間・APの偏りを示した。
- [ ] 物理的散逸と参照仕事を区別し、対称12/12のばね/速度項を**未測定**と明示した。
- [ ] 解析レポートで`Hspring/Hdamp/Hmixed/Hstage`のどれが支持されるか数値根拠で分類し、解釈の限界を述べた。

### Codexへの最終依頼

**12/6・Cover3/3・900 Paの接触仕事について、既存の総接触荷重の数値を一切変えず、実際の最終採用反復で使われた`f_spring`と`f_damp`由来の節点荷重を分離して記録してください。** 完走済み12/6の元runと全`dat`が一致するかをゲートにして、最大1回だけ新しい診断runを行い、各成分の`Q Δq`を同じ時間窓で積分してください。特に右の負の接触参考仕事約−1.206 mJがどちらに由来するかを、左・両側和・主要モード・AP/時間分布まで含めて検証し、`report.md`に観測と解釈を分けて報告してください。数値整合が失敗したら結論を保留し、追加本番計算は勝手に行わないでください。
