# Kanaya振動再現：既存1ランから原因を絞る3項目の後処理検証
> Codexへの実装・解析指示書 / 2026-10-21  
> 対象: 既存の `simulation-main-kanaya-diagnostic` 作業環境。**新たな連成計算は実行しない。**

## 0. 結論を先に：今回の目的

目的は Kanayaらの模型のような**周期的な声門閉鎖を伴う自励振動**が再現できない理由を、既存の **Ps=1000 Pa / L_vt=0 / 100 modes / nstep=20000 / dt=1e-5 s** のランの生データだけで、次の3論点に切り分けること。

1. **AP幅仮説** — 現行50断面の `Amin` は真の局所狭窄か、左右交線が存在する **共通AP積分幅の減少**による見かけ上の極小か。
2. **静的つり合い仮説** — 大きな平均開大（片側の平均流体開大力約+112 mN、Amin平均約6.62 mm²）は、既存モードの剛性・平均一般化力からどの程度説明できるか。どのモード・表面領域が平均変形に寄与しているか。
3. **自励仕事仮説** — モード仕事の評価窓合計が負でも、**平均変形の仕事と振動成分の仕事**を分離すると、どのモードにエネルギー供給があり、維持・成長・減衰のいずれが示唆されるか。

**達成水準**: 3項目それぞれに、定量値・図・実装上の検証・「支持/否定/未判定」・次に行う比較ランを最大1条件だけ記述する。いずれかが入力不足で不可能なら「不可能な理由と必要な入力」を明示。**後処理だけで因果の確定を装わない。**

## 1. 作業場所と絶対禁止事項

作業開始時に `pwd; git status --short; git branch --show-current; git rev-parse HEAD`、`AGENTS.md`、実データのパス・サイズ・更新時刻を記録し、現在の作業中のコード／ブランチを正とする。GitHub上の指示書ブランチとローカルworktreeでコミット・未コミット差分が異なる可能性がある。

想定ワークツリー: `/tmp/simulation-main-kanaya-diagnostic-20261021` またはこれに相当するCodexの既存分離worktree。ローカルに存在する場所を探す。**元の `/home/kajiyama/code/simulation-main` を変更しない**。既存worktreeでも未コミット変更、既存ファイルを保護する。

**絶対禁止**
- `./simulation`、本番再計算、80-stepスモークなど、時間積分を伴うソルバを**一切起動しない**。
- `src/simulation/`、`src/fluid/`、`src/forces/`、`src/structure/`、`include/` 内の物理・数値計算コード、`input/`、出力済みの `output/kanaya_diagnostic/report.md` およびCSV/binを変更しない。
- 接触条件、圧力回復幅、減衰、剛性、接続声道、nxsup、50断面ロジックの物理式を変更しない。
- `git reset`、checkoutによる未コミット変更破棄、古いソースへの巻き戻し、巨大生データのコピーやgit add、既存成果物の上書きはしない。
- 以前のレポートに書かれた「負圧比14.7倍」「正の変動仕事2.61 µJ」等を入力値としてハードコードしない。必ず元データから再計算し、既存集計と照合する。
- `NaN` の断面を0面積や閉鎖とみなさない。

**追加するコードは** `tools/kanaya_diagnostics/three_checks.py` を中心とした読み取り専用後処理。どうしても断面交線の再構成をPython単独で実装するのが不適切なら、独立した後処理用C++ helperを`tools/kanaya_diagnostics/`に作ってよいが、ソルバ本体の再コンパイル・変更・時間発展は行わない。大きな新規依存導入を避ける。標準Python、numpy、pandas、matplotlib、既存依存を優先。

## 2. 入力ファイルと時間・単位の規約

既存の `output/kanaya_diagnostic/` 内（場所を発見して利用）:
- `modal_native.bin`: **19999 interval × 2 sides × 100 modes × 7 float64 little-endian**。列は `Q_pos,Q_neg,Q_fluid,Q_contact_adopted,q_n,qdot_n,delta_q`。ファイルサイズとメタデータを照合してから `numpy.memmap(dtype="<f8",shape=(n,2,nModes,7))` として読む。一度に全配列をRAMへコピーしない。
- `run_metadata.json`（実dt、最終区間除外、モード数、仕事の定義、正規化情報）、`flow_geometry.csv`（20,000行のステップ同期Amin, xMin,xSep,圧力, Ug,Lg/Rv/Rk、probe等）、`modal_energy_summary.csv`、`modal_selected_timeseries.csv`、`contact_regions.csv`、`force_regions.csv`、`consistency.csv`。無いものがあれば`missing_inputs`として記録。
- `output/dat/area.dat`, `section_x.dat`, `gap_cubed.dat`, `pressure.dat`（提供物や当該本番worktreeの同一run）。各行: 先頭列 `step`、その後に**50断面**、通常 **5ステップ間隔**。異なるrunを混ぜないため`n=5k`について `flow_geometry.csv` と面積・最狭位置・剥離位置を照合。
- 元VTU/NAS/周波数ファイル（左右のsoft版）、`params_used.txt`、`manifest.txt`、必要なら本番 `result/deform_combined*.vtu`（**末尾5,000ステップのみ**、20ステップ間隔の可能性）。inputは読み取り専用。
- `kanaya_selected_modes.npz` は5モードの共有用派生データ。**幾何形状の厳密な再構成には全100モードの `modal_native.bin` が必要**。選択5モードだけで再構成した値を「本番の断面」と呼ばない。

**時間軸の重要なオフセット**
- `modal_native.bin` のレコード`n`は、流体荷重が使われる `t_n=n*dt` 時点の `q_n`, `qdot_n` と `Δq=q_(n+1)-q_n` を含む。**最終区間n=19999は記録されず、仕事積分から除外済み**。ソルバの総20,000ステップのうち19,999区間だけが対象。
- `area.dat` 等はステップ`n`の**更新前座標**で流路断面を記録する。
- `writeVTKCombined` は `state.uf2u()` 後、すなわちステップ`n`の**後の状態（t_(n+1)）**を保存する。ファイル名だけから時刻を決めずソースのカウンタを解読し、`q_n`との1ステップずれを合わせる。
- `State::mode2uf` の本体に既知の境界バグがあり最終更新がスキップされる。**除外済み最終区間を仕事・整合性・変形再構成評価に再投入しない。**
- 座標・モードを用いる再構成は `u[mm]=Σ φ_m q_m * 1e3`（`state.cpp`参照）。`x,y,z`はmm、`QΔq`は診断のJ相当の仕事。静的な`q`はmmではなく正規化された一般化座標。`f[Hz]`から`ω=2πf`。
- 以前の診断で選択モードのGram非対角要素最大絶対値約0.072。**厳密な全モード直交性やモード間の独立エネルギー保存を仮定しない**。質量行列の近似の限界を明記。

## 3. 検証A：AP有効幅と真のくびれを分ける【最優先】

### 問題の正確な定義

現行 `src/fluid/channelSections.cpp`（特に `intersectAtX`, `sectionSegments`, `positiveGapIntegralsYZ`）は、左右の各`j`にx平面交点があるときだけ、物理zの共通区間を分割積分する。各`j`の個別x交点がない区間は補間して埋めない。`GapIntegrals`内では
\[
A(x)=\int_{Z_{\rm common}(x)}\max(h(x,z),0)\,dz,\quad
G_3(x)=\int_{Z_{\rm common}(x)}\max(h(x,z),0)^3\,dz .
\]
`span`は**共通交線が存在するAP幅**（積分区間幅の合計）であり、**開口しているAP幅**とは違う。これを混同しない。

### 実装する後処理

1. `area.dat`・`section_x.dat`・`gap_cubed.dat` から各時刻の `A(x)`, `G3(x)`, `x`、無効断面を整理し、`argmin(A)`とその隣接2～3断面、入口側・出口側の変化を示す。末端NaNは閉鎖として数えない。
2. 全100モードの`q_n`と初期幾何・正規化モードから、選んだ代表時刻の**左右の変形後の節点座標**を再構成する。左面`geomL`、右は**y鏡映・モードuy符号反転、initialGapMm=0に伴う左右剛体並進**を`simulation.cpp`と同じ順序で実施する。固定系のx,zもモード変形を反映する。 `nxsup=42`, `nsurfl=64`, `nsurfz=69` はレポートの参考値であり、現行ソース・本番manifestと照合する。
3. 現行と同じ「50個の共通全体x区間」「AP行ごと内側y交点を選択」「zの左右共通区間だけを和」「負gapを0で切る」アルゴリズムで`A_recon`, `G3_recon`, `S_common`を**独立再評価**する。`S_open`も正gap部分のz幅を計算し、両者を別列にする。重複区間・穴あき区間を二重計数しない。1点だけの交点はAP幅0。
4. `A_recon`と`area.dat`の一致、`x_recon`と`section_x.dat`の一致を**最初の厳格なゲート**にする。面積一致まで到達しなければ`S`に物理的な結論を出さない。DATの印字丸め誤差を明示し、`flow_geometry.csv`の`Amin`も独立照合する。
5. 各断面で`h_mean_common=A/S_common`、`h_mean_open=A/S_open`（0幅ならNaN）を作り、`S_common`の減少率、`S_open`の減少率、平均gap、Aの減少率を比較する。「Aが小さい」原因が**狭いgap**か**少ない共通z被覆**かを定量化。参照幅は共通する正常断面の幅や初期AP範囲を明記し、`A/S`の意味を過大解釈しない。
6. 選択時刻は少なくとも`n=0`付近、`t≈0.10,0.125,0.15,0.175,0.195` s、および評価窓中のAmin極小/極大、xSep変位極値から重複を除いて**10～20フレーム**。高コストなら代表6フレーム以上に縮めて理由を記す。連成計算の再実行はしない。
7. 対応する最狭断面の`z`ごとの`yL,yR,h`曲線を少なくとも2時刻描く。見かけの幅減少と真のgap縮小を図で区別する。
8. 50断面のAがz被覆によって歪むと`Lg,Rv,Rk`にも波及する可能性を述べる。ただし**現行の`A`や流体式を書き換えない**。

**出力**: `ap_coverage.csv` (step,time,plane_index,x,A_observed,A_reconstructed,G3_observed,G3_reconstructed,common_span_mm,open_span_mm,mean_gap_common_mm,mean_gap_open_mm,valid,area_error,reason)、`figures/01_ap_coverage.png`、`figures/02_ap_gap_profile.png`。

**合否**: 数値の再現性を優先。AP被覆を調べた結果が「十分幅があり真のくびれ」「幅の減少に大きく依存」「再構成不一致で未判定」のどれかを数値根拠つきで判定。

## 4. 検証B：平均流体荷重とモードの準静的変形

### 方程式と注意

既存の各モードのモデルが
\[
\ddot q_m+2\zeta_m\omega_m\dot q_m+\omega_m^2q_m=Q_{m}^{fluid}+Q_m^{contact}
\]
の形であることを**現行のNewmark呼び出しとモード正規化から確認**する。違うならその実際の式を使う。

評価窓の準静的予測は
\[
\bar q_{m,\rm static}=\frac{\bar Q_m^{fluid}+\bar Q_m^{contact}}{\omega_m^2},
\quad
R_{m,\rm static}=\omega_m^2\bar q_m-(\bar Q_m^{fluid}+\bar Q_m^{contact}) .
\]
ただし非定常なら `<qddot>+2ζω<qdot>` も含めて解釈し、単純な`q_stat`一致を強制しない。`qddot`は`qdot_n`差分から近似可能だが、差分方式・終端処理・高周波ノイズを記録する。

### 実装手順

1. 全100モード／左右別に `q_n,Qfluid_n,Qcontact_adopted_n,qdot_n,delta_q`を取得する。`Qfluid = Qpos + Qneg`を再検証。既存`modal_energy_summary.csv`の積算`ΣQΔq`と再集計値が一致しなければ停止する。
2. 窓を`0–0.05`（ランプ）、`0.05–0.10`、`0.10–0.15`、`0.15–0.19999`、`0.10–0.19999`とする。実ランプ時間`0.05s`を確認し変更しない。各窓で`mean(Qfluid), mean(Qcontact), mean(q), mean(qdot), qdot_RMS_AC, approximate_mean(qddot), qstatic, static_residual`を保存。残差は`ω²\bar q`等の規模を使って正規化し、分母が小さいときは絶対値だけを示す。
3. 主にモード1,4,5,7,35、そして**平均`|Q|`または平均`|q|`上位5モード**を重点可視化するが、合計・ランキングは全100モードから算出する。
4. 再構成した平均変位`u_mean(x,z)=Σφ_m\bar q_m ×1e3 [mm]`のうち、声門開大方向y成分とx成分を別々に表示。**モニタ点、最狭近傍、AP中央**で各モードの寄与を列挙する。`monitor`と`near sep`が時間で節点切替される場合は「同じ材料点」扱いをせず、固定材料点の比較も実施する。
5. `<Qfluid>` と `<Qcontact>`が作る予測開大 `Σφ(Qbar/ω²)`の大きさを比較。接触反復から採用した`Qcontact`は予測時刻の荷重であり、左端流体荷重と同時刻とは断言しない。空間的な正圧広がりは`force_regions.csv`と照合。
6. `ω`の確認、モード正規化Gramの非直交性、線形モードの大変形適用限界を`limitations`として明記。「モデル内で説明される平均変形」と「物理的に正しい剛性・荷重」は別。

**出力**: `static_balance.csv`, `figures/03_static_modal_balance.png`。

**判定**: 開大は主に平均流体荷重と低次モードで説明可能か。それとも顕著な残留運動／高次モード・接触・時刻ずれがあり、準静的近似自体がまだ不適切か。

## 5. 検証C：振動成分の仕事と時間発展

### 最も重要な定義

これまでの`W=Σ Q_n Δq_n`は**平均変形の仕事を含む**。そこで各窓で厳密に成り立つ加法分解:
\[
W_m=\underbrace{\bar Q_m\sum_n\Delta q_{m,n}}_{W_{\rm mean-force}}
+\underbrace{\sum_n(Q_{m,n}-\bar Q_m)\Delta q_{m,n}}_{W_{\rm fluctuating-force}} .
\]
`bar Q`は同じ窓のサンプル平均。`W_fluctuating-force`は**厳密な自励振動仕事ではない**。変動する平均位置、ドリフト、複数周期の混合が含まれる。実装では上の足し算と生の`QΔq`の完全一致を必須にする。

さらに`q'=q-\bar q_\mathrm{slow}`,`Q'=Q-\bar Q_\mathrm{slow}`とした帯域・トレンド分離も**補助解析**として実施。移動平均/ローパスのカットオフを変えて（例20,30,40 Hz）、符号と順位が安定か調べる。基本周期が約60–130 Hzの混合で、窓長0.1秒は短い。カットオフに依存する結論は「不確実」とする。フィルター境界のアーティファクトを捨てるか感度として報告する。

### 実装手順

1. 窓を`0.05–0.10`, `0.10–0.15`, `0.15–0.19999`, `0.10–0.19999`に分け、左右×100モードで `Wtotal,Wmean_force,Wfluctuating_force,Wpos,Wneg,Wcontact_reference,D_total`を再計算。左端時刻`Q_n`×`delta_q`を厳守。`Wpos+Wneg=Wfluid`の恒等性も確認。
2. 減衰`D_total=Σ2ζω (qdot_n)^2 dt`と、参考の`D_ac=Σ2ζω(qdot_n-mean(qdot))² dt`を併記。**`W_fluctuating_force-D_ac`を厳密なエネルギー保存残差とは呼ばない**。質量正規化のGram限界を明示。
3. モード1,4,5,7,35と `|Wfluctuating_force|`上位5を重点表示。各窓で`q_AC_RMS`,`qdot_AC_RMS`,`Q_AC_RMS`,相関、`Q`と`qdot`の主要スペクトル位相を計算。ただし多峰性・帯域不一致・弱いコヒーレンスでは位相を一値に要約しない。
4. `0.10–0.15`と`0.15–0.19999`で**同じ長さ/同じ方法**による振動振幅・流体仕事率（各窓の仕事/窓長）を比較し、振幅が減衰・横ばい・成長か判定。窓長約50 msは45 Hzモードで約2周期なので統計の不確かさを示す。
5. 有意な狭帯域振動があれば、モード固有の連続ピーク/ゼロクロスによって**周期ごと**の`ΣQΔq`,`D`,振幅を計算。ただしピークが不安定なら無理に周期を抽出せず、短時間窓のみ提示。「130 Hzに近いスペクトルピーク」は閉鎖・自励の証明ではない。
6. 追加の原因解析: `flow_geometry.csv`のAmin, Ug, xSepと`modal_native`のモード運動の同時刻関係を示す。負圧による仕事`Wneg`が正であっても、そのモードがy閉鎖変形を持つかを検証Bの形状寄与と照合する。必要なら流れ方向で異なるx位置の局所開口の位相を、検証Aの**同じAP被覆を持つ**代表断面だけで評価し、疑似位相差を避ける。
7. 「平均開大で閉鎖しない」と「自励の成長率が小さい」を混同せず、**閉鎖不足量**の目安（Amin最低5.195 mm²、閉鎖閾値0.01 mm²）と動的Amin振幅・その時間変化を最後に対比。

**出力**: `dynamic_work.csv`, `figures/04_dynamic_work.png`, `figures/05_amplitude_and_gap.png`。

## 6. 共通のデータ検証とSTOP条件

本体に手を加える前に分かる問題を後処理で潰すため、次を**自動アサーション**にする。

- `run_metadata.json`で`complete=true`、`integrated_intervals=19999`、`excluded_intervals=1`、`nmode=100,100`、`dt_s=1e-5`等を確認。別ランなら数値を決め打ちせず整合する実測条件を採用。
- `consistency.csv`の非除外レコードが全て`valid=1`で残差が定義済みの許容範囲内。失敗なら力・仕事の物理解釈を停止。
- `modal_native.bin`の全ファイルバイト数=`nInterval*2*nMode*7*8`。`NaN`や`inf`は理由を追跡、勝手に0補間しない。
- `Qpos+Qneg=Qfluid`、`Wpos+Wneg=Wfluid`、各窓の仕事合計が`modal_energy_summary.csv`の保存済み`all/evaluation`などに丸め差内で一致。50ms分割窓の足し合わせも一致。
- `area.dat`と`section_x.dat`は行ごとに同じステップIDが並び、`n%5=0`である。出力の実桁数に応じた許容値を設定・報告。
- **断面再構成が既存`A(x)`に合わない場合**はAP被覆率仮説を肯定/否定しない。右鏡映、初期gap補正、qの時刻、ソースのモード正規化を重点調査。
- `modal_native.bin`がなく選択モードのみの場合: 静的・振動仕事は選択モードに限定し、**全100モードの結論を出さない**。AP形状再構成については「入力不足」と明記し推測で埋めない。
- もともとの`report.md`と今回の`postprocess_3checks/`に矛盾があれば、**採用した元ファイル/時刻/式**を記録して原因を説明してから判断する。
- 読み込み・図生成がメモリ不足なら`memmap`、選択時刻抽出、集約を利用。`modal_native.bin`全コピー・重い外部パッケージインストール・本番再計算はしない。

## 7. クリーンな成果物のツリー

コードの変更点は原則**1 Pythonスクリプト＋1テスト**、AP再構成の必要がある場合のみ独立helperを追加。出力は既存レポートを上書きせず専用フォルダ内に置く。

```text
simulation-main-kanaya-diagnostic/
├── docs/plans/
│   └── KANAYA_THREE_CHECKS_POSTPROCESS.md
├── tools/kanaya_diagnostics/
│   ├── analyze.py                       # 既存、変更しない
│   ├── three_checks.py                  # 新規。既存ランを読むだけ
│   └── geometry_replay.*                # 任意。独立した読取専用補助
├── tests/
│   └── test_kanaya_three_checks.py      # 新規、合成ケースのみ
└── output/kanaya_diagnostic/
    ├── modal_native.bin                 # 既存、変更しない
    ├── ...                              # 既存CSV/JSON/report.md
    └── postprocess_3checks/
        ├── report.md                    # 1ページ目に結論、仮説判定
        ├── manifest.json                # ソースとハッシュ、単位、条件、制約
        ├── ap_coverage.csv
        ├── static_balance.csv
        ├── dynamic_work.csv
        ├── validation.json              # assertion合否、再構成誤差
        └── figures/
            ├── 01_ap_coverage.png
            ├── 02_ap_gap_profile.png
            ├── 03_static_modal_balance.png
            ├── 04_dynamic_work.png
            └── 05_amplitude_and_gap.png
```

**出力が既に存在すれば上書きしない。** 既存runの各出力のSHA256を`manifest.json`に記録。再解析が必要な場合は明示的な新しいサブフォルダ名を使う。元`report.md`を変更しない。`git commit/push`は指示なしにしない。

## 8. テスト

`tests/test_kanaya_three_checks.py`に外部本番データ不要の合成例を用意:

- 全長17 mmの一定幅・一定gapの2面: `A=S*h`、`G3=S*h³`、`S_common=S_open`を再現。
- **AP交線の片側欠落**: `A`は減るが`h_mean_common`は一定で、幅不足を検出できる。
- **真のくびれ**: `S_common`一定でgapが減るため`A`が減る。
- **部分接触**: `S_common>S_open`, `A=∫max(h,0)dz`となる。
- 断面に交点がない/NaNのケース: 未定義を閉鎖と取り違えない。
- Q定数、qが単調増加: `Wtotal=QΔq`で変動仕事0。
- 零平均正弦Qと速度が同位相/逆位相: 変動仕事の符号が正/負で変わる。
- q変動0、平均開大のみ: 自励振動と判定しない。
- `q_n`とVTUの `q_(n+1)` を取り違えるケースで時刻照合が失敗する。
- `ω=2πf`と`ω²qbar≈Qbar`の静的人工ケース、非定常ケースの`qddot`補正。

いずれかが失敗した場合、後処理レポートに物理結論を記載しない。テストのためにソルバを実行しない。

## 9. 最終report.mdの書き方・意思決定

**冒頭に短く3結論**（下記の3点に実数・確信度を入れる）:
1. `Amin` の極小は、AP被覆が減ることによる偽のくびれか、平均gapの縮小か。その割合/時間帯/再構成精度。
2. 平均開大はどのモード（1,4,5等）の平均一般化力・準静的変形で説明できるか。実測平均開大との一致度。
3. 自励を支える振動仕事は後半でも正か。50ms前後半で振幅・仕事率は成長、横ばい、減衰のいずれか。

続いて以下の**判断表**を必ず埋める:

| 診断結果 | 次の1手（まだ実行しない） |
|---|---|
| 面積最小が共通AP幅の縮小に支配される | 断面積/剥離選択に対する**幾何アルゴリズムだけ**の修正案と、その副作用予測を提示 |
| 幅の問題は小さく、平均荷重で過大開大が説明できる | **現行材料・初期gap・荷重対象面**のうち識別力が高い1条件だけの比較実験を設計 |
| 平均開大は妥当だが、変動流体仕事が不足/負 | 流れ方向の変形位相と圧力回復近似を優先する比較案を提示 |
| 仕事は正、振幅が成長中 | まず**同条件の時間延長**が必要か、現データの振幅傾向から判断 |
| データ/再構成不一致 | 追加計算せず、何が不足・不一致かを具体的に報告 |

最終報告には下記だけを簡潔に載せる:
- 作業前後の`git status`と既存成果物が変更されていない証拠。
- 3診断の数値結果、図5枚以内、重要な限界、仮説の判定。
- **次に行う比較ランは最大1条件**。変更箇所と、変更すると何が増減するか・その識別理由を示すが、**このタスクでは実行しない**。
- `postprocess_3checks/report.md`への絶対パス、および成果物ツリー（深さ3）。

### 完了条件
追加の連成計算 **0回**、数値モデルの変更 **0件**、既存出力の上書き **0件**。既存の1ランから3仮説の支持・否定・未判定を明確にし、次の1条件を選ぶだけの定量的根拠が残ること。
