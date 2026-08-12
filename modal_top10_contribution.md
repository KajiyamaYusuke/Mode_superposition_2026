# 上位10モードの変位寄与率診断 実装仕様

## 目的

声帯の最終変位は、各固有モードのモード座標 \(q_m(t)\) とモード形状
\(\boldsymbol{\phi}_{m,i}\) の重ね合わせによって計算されている。

本改修では、現在の数値計算結果と既存出力を変更せず、診断用の新規ファイルを追加して、左右それぞれについて次を出力する。

- 表面全体の3次元変位に対する上位10モード
- 表面全体の開閉方向 \(u_y\) に対する上位10モード
- 監視点の開閉方向 \(u_y\) に対する上位10モード
- 各モードの大きさ寄与率
- モード間の強め合い・打ち消しを含む符号付き寄与率
- 全体のモード相殺の強さ

既存の以下の出力は変更しないこと。

- `modal_contribution.csv`
- `modal_dominant.csv`
- `displace.dat`
- `displace_xy.dat`
- その他すべての既存出力

新規診断ファイルとして、次だけを追加する。

```text
output/modal_top10.csv
```

---

## 現行コードでの変位合成

対象リポジトリ：

```text
KajiyamaYusuke/Mode_superposition_2026
```

対象ブランチ：

```text
main
```

### `src/state.cpp`

`State::mode2uf()` と `State::mode2ufSurface()` では、節点 \(i\) の変位を次のように計算している。

\[
\boldsymbol u_i(t)
=
10^3\sum_{m=1}^{N_m}
q_m(t)\boldsymbol\phi_{m,i}
\]

ここで \(10^3\) は m から mm への換算である。

コード上では、各モードについて次を加算している。

```cpp
predUx += modeAtPoint.ux * qi * 1.0e3;
predUy += modeAtPoint.uy * qi * 1.0e3;
predUz += modeAtPoint.uz * qi * 1.0e3;
```

したがって、モード \(m\) 単独に由来する節点変位は、

\[
\boldsymbol u_i^{(m)}
=
10^3 q_m\boldsymbol\phi_{m,i}
\]

であり、

\[
\boldsymbol u_i
=
\sum_m\boldsymbol u_i^{(m)}
\]

が成立する。

### 各関数の役割

- `State::mode2ufSurface()`
  - 接触反復中に表面節点だけを復元する。
- `State::mode2uf()`
  - 最終確定前に全節点を復元する。
- `State::uf2u()`
  - `qf → q`、`qfdot → qdot`、`predictedDisp → disp` として状態を確定する。

---

## 現在すでに存在するモード診断

### `src/simulation.cpp`

現在は次の2ファイルを出力している。

#### `modal_contribution.csv`

各時刻、左右、全モードについて次を出力する。

- モード番号
- 固有振動数
- \(q_m\)
- \(\dot q_m\)
- 監視点でのモード由来 \(u_x,u_y,u_z\)
- 表面上のモード由来 \(u_x,u_y,u_z\) のRMS
- 表面上の3成分合成RMS

現在の `surface_rms_norm_mm` は、概ね次である。

\[
r_m
=
|q_m|10^3
\sqrt{
\frac{1}{N_s}
\sum_{i\in S}
\|\boldsymbol\phi_{m,i}\|^2
}
\]

#### `modal_dominant.csv`

現在は次の最大モードを1つずつ記録する。

- 監視点 \(u_y\) の絶対値が最大のモード
- 表面3成分合成RMSが最大のモード

既存の計算を利用できるが、現在のRMSを単純に百分率化するだけでは、モード間の符号や位相による相殺を表せない。

---

## 寄与率の定義

単一の「寄与率」では意味が曖昧になるため、以下の2種類を必ず併記する。

## 1. 大きさ寄与率

表面上のモード \(m\) の変位ベクトルを、

\[
\boldsymbol U_m
=
\{
u_{x,i}^{(m)},
u_{y,i}^{(m)},
u_{z,i}^{(m)}
\}_{i\in S}
\]

とする。

そのノルムを、

\[
R_m=\|\boldsymbol U_m\|
\]

とし、大きさ寄与率を、

\[
\boxed{
A_m
=
\frac{R_m}{\sum_kR_k}
}
\]

と定義する。

### 性質

- \(0\le A_m\le1\)
- 全モードの合計は1
- 上位10モードの順位付けに使いやすい
- 符号やモード間の相殺は表さない

これは「最終変位の何％を作ったか」ではなく、

> 全モードの変位振幅総量に対する、そのモードの大きさの割合

である。

## 2. 符号付き投影寄与率

最終合成変位を、

\[
\boldsymbol U=\sum_k\boldsymbol U_k
\]

とする。

モード \(m\) の合成変位方向への寄与を、

\[
\boxed{
B_m
=
\frac{
\langle\boldsymbol U_m,\boldsymbol U\rangle
}{
\langle\boldsymbol U,\boldsymbol U\rangle
}
}
\]

と定義する。

### 性質

\[
\sum_mB_m=1
\]

が数値誤差の範囲で成立する。

- \(B_m>0\)：最終変位を強めている
- \(B_m<0\)：他モードを打ち消している
- \(B_m>1\)：そのモード単独の変位は最終変位より大きく、他モードに相殺されている

### 相殺例

\[
u^{(1)}=+2\ {\rm mm},
\qquad
u^{(2)}=-1\ {\rm mm}
\]

なら最終変位は、

\[
u=1\ {\rm mm}
\]

である。

大きさ寄与率は、

\[
A_1=66.7\%,
\qquad
A_2=33.3\%
\]

となる。

一方、符号付き投影寄与率は、

\[
B_1=200\%,
\qquad
B_2=-100\%
\]

となり、

\[
200\%-100\%=100\%
\]

である。

この定義により、mode 2が最終変位を打ち消していることを表現できる。

---

## 開閉方向だけの寄与率

3次元変位全体とは別に、開閉方向 \(u_y\) だけの寄与も計算する。

表面上の \(y\) 方向モード変位を、

\[
\boldsymbol U_{m,y}
=
\{u_{y,i}^{(m)}\}_{i\in S}
\]

とする。

### 開閉方向の大きさ寄与率

\[
\boxed{
A_{m,y}
=
\frac{\|\boldsymbol U_{m,y}\|}
{\sum_k\|\boldsymbol U_{k,y}\|}
}
\]

### 開閉方向の符号付き投影寄与率

\[
\boxed{
B_{m,y}
=
\frac{
\langle\boldsymbol U_{m,y},\boldsymbol U_y\rangle
}{
\langle\boldsymbol U_y,\boldsymbol U_y\rangle
}
}
\]

ここで、

\[
\boldsymbol U_y=\sum_k\boldsymbol U_{k,y}
\]

である。

---

## 監視点 \(u_y\) の寄与率

監視点では、

\[
u_y=\sum_mu_y^{(m)}
\]

である。

### 大きさ寄与率

\[
A_{m,\mathrm{probe}}
=
\frac{|u_y^{(m)}|}
{\sum_k|u_y^{(k)}|}
\]

### 符号付き寄与率

\[
B_{m,\mathrm{probe}}
=
\frac{u_y^{(m)}}{u_y}
\]

ただし、監視点の合成 \(u_y\) がゼロに近い場合、符号付き寄与率は発散する。ゼロ交差付近では `nan` を出すこと。

---

## 相殺係数

全体としてモード相殺がどの程度強いかを次で評価する。

\[
\boxed{
C
=
\frac{
\sum_m\|\boldsymbol U_m\|
}{
\left\|\sum_m\boldsymbol U_m\right\|
}
}
\]

### 解釈

- \(C\simeq1\)：各モードがほぼ同方向に加算されている
- \(C>1\)：モード間の相殺がある
- \(C\gg1\)：大きなモード変位同士が強く打ち消し合っている

`surface_xyz`、`surface_uy`、`probe_uy` ごとに計算する。

---

## 上位10モードの選び方

上位10モードは、符号付き投影率ではなく、大きさ \(R_m\) の降順で選択する。

符号付き投影率で順位付けすると、強い負の打消しモードを見落とす可能性がある。

したがって、

```cpp
sort by modalNorm descending
```

として上位10個を選び、それぞれについて次を出力する。

- 大きさ寄与率
- 符号付き投影寄与率
- 累積大きさ寄与率

評価範囲 `scope` ごとに独立して上位10を求める。

```text
surface_xyz
surface_uy
probe_uy
```

---

## 新規出力ファイル

ファイル名：

```text
output/modal_top10.csv
```

推奨ヘッダ：

```text
step,time_s,side,scope,rank,mode_index,frequency_hz,
q,modal_norm_mm,magnitude_ratio,signed_projection_ratio,
cumulative_magnitude_ratio,total_displacement_norm_mm,
cancellation_factor
```

実際のCSVヘッダは1行で出力する。

```cpp
fmodalTop10
    << "step,time_s,side,scope,rank,mode_index,frequency_hz,"
    << "q,modal_norm_mm,magnitude_ratio,signed_projection_ratio,"
    << "cumulative_magnitude_ratio,total_displacement_norm_mm,"
    << "cancellation_factor\n";
```

比率はCSV内部では `0–1` を基本とする。解析スクリプトやグラフで100倍して百分率表示する。

符号付き投影率については、負値や1を超える値を許容する。

---

## 変更対象

## 必須変更

### `src/simulation.cpp`

次を追加する。

1. `modal_top10.csv` のストリーム
2. CSVヘッダ
3. モード別変位の集計構造体
4. 上位10モードを計算する関数またはラムダ
5. 既存の `n % params.nwrite == 0` ブロック内で左右について呼び出す処理

## 原則として変更しないファイル

- `src/state.cpp`
- `src/modalProjector.cpp`
- Newmark積分部分
- 流体力計算
- 接触力計算
- 既存CSV出力

診断処理だけで実現できるため、運動方程式、モード座標、変位復元方法は変更しないこと。

コードが長くなりすぎる場合のみ、次を新設してもよい。

```text
src/ModalContributionDiagnostics.cpp
include/ModalContributionDiagnostics.h
```

ただし、最初の実装は `simulation.cpp` 内の局所関数またはラムダでもよい。

---

## 推奨データ構造

```cpp
struct ModalContributionEntry {
    int modeIndex = -1;
    double frequencyHz = 0.0;
    double q = 0.0;
    double modalNormMm = 0.0;
    double magnitudeRatio = 0.0;
    double signedProjectionRatio =
        std::numeric_limits<double>::quiet_NaN();
};
```

---

## 表面3次元変位の計算手順

対象は `geom.surfp[i][j]` に登録された表面節点とする。

同一節点の重複集計を避けること。`State::surfacePointIds` を外部から参照できない場合は、診断初期化時に表面節点IDの一意リストを作る。

### 手順1：合成変位を計算

各表面節点 \(i\) について、

```cpp
totalUx[i] = 0.0;
totalUy[i] = 0.0;
totalUz[i] = 0.0;
```

とし、全モードについて、

```cpp
const double scaleMm = state.q[m] * 1.0e3;
const auto& phi = modes.modes[m][pid];

totalUx[i] += scaleMm * phi.ux;
totalUy[i] += scaleMm * phi.uy;
totalUz[i] += scaleMm * phi.uz;
```

とする。

`state.disp - geom.points` から合成変位を求めてもよいが、符号付き寄与率の恒等式を確実に保つため、モード和から合成変位を作るほうが安全である。

### 手順2：モードごとのノルムと投影を集計

モード \(m\) ごとに各表面節点で、

```cpp
const double ux = scaleMm * phi.ux;
const double uy = scaleMm * phi.uy;
const double uz = scaleMm * phi.uz;
```

を計算する。

表面3次元変位について、

```cpp
modalNormSquared[m] += ux*ux + uy*uy + uz*uz;

projectionNumerator[m] +=
    ux * totalUx[s]
  + uy * totalUy[s]
  + uz * totalUz[s];
```

とする。

合成変位ノルムは、

```cpp
totalNormSquared +=
    totalUx[s]*totalUx[s]
  + totalUy[s]*totalUy[s]
  + totalUz[s]*totalUz[s];
```

である。

### 手順3：割合を計算

```cpp
modalNorm[m] = std::sqrt(modalNormSquared[m]);
sumModalNorm += modalNorm[m];
```

その後、

```cpp
magnitudeRatio[m] =
    modalNorm[m] / sumModalNorm;

signedProjectionRatio[m] =
    projectionNumerator[m] / totalNormSquared;
```

とする。

---

## 表面 \(u_y\) の計算手順

表面3次元変位と同時に集計する。

```cpp
modalUyNormSquared[m] += uy * uy;
uyProjectionNumerator[m] += uy * totalUy[s];
totalUyNormSquared += totalUy[s] * totalUy[s];
```

最終的に、

```cpp
modalUyNorm[m] =
    std::sqrt(modalUyNormSquared[m]);

uyMagnitudeRatio[m] =
    modalUyNorm[m] / sumModalUyNorm;

uySignedProjectionRatio[m] =
    uyProjectionNumerator[m] / totalUyNormSquared;
```

とする。

---

## 監視点 \(u_y\) の計算手順

```cpp
const double probeModeUy =
    state.q[m] * 1.0e3 * modes.modes[m][probeId].uy;
```

全モード和：

```cpp
probeTotalUy += probeModeUy;
```

大きさの和：

```cpp
sumAbsProbeModeUy += std::abs(probeModeUy);
```

割合：

```cpp
probeMagnitudeRatio[m] =
    std::abs(probeModeUy) / sumAbsProbeModeUy;

probeSignedRatio[m] =
    probeModeUy / probeTotalUy;
```

ただし `probeTotalUy` がゼロ付近なら `probeSignedRatio` は `nan` とする。

---

## ゼロ付近の処理

合成変位がほぼゼロの場合、符号付き投影率と相殺係数は不安定になる。

```cpp
constexpr double normEpsMm2 = 1.0e-24;
constexpr double scalarEpsMm = 1.0e-12;
```

を目安とし、

```cpp
if (totalNormSquared <= normEpsMm2) {
    signedProjectionRatio = NaN;
    cancellationFactor = NaN;
}
```

とする。

監視点については、

```cpp
if (std::abs(probeTotalUy) <= scalarEpsMm) {
    probeSignedRatio = NaN;
}
```

とする。

ゼロ変位時に0を出すと、

- 各モードの変位もゼロ
- 大きなモード同士が完全相殺

を区別できないため、`nan` が適切である。

---

## 並べ替え

モード数は100程度なので、通常の `std::sort` で十分である。

```cpp
std::sort(
    entries.begin(),
    entries.end(),
    [](const auto& a, const auto& b) {
        return a.modalNormMm > b.modalNormMm;
    }
);
```

出力件数は、

```cpp
const int topCount =
    std::min<int>(10, entries.size());
```

とする。

累積寄与率は順位順に、

```cpp
cumulative += entries[rank].magnitudeRatio;
```

として出力する。

---

## 呼び出し位置

`Simulation::run()` 内の既存部分：

```cpp
if (n % params.nwrite == 0) {
    writeModalDiagnostics(...);
}
```

と同じブロック内で呼ぶ。

```cpp
if (n % params.nwrite == 0) {
    writeModalDiagnostics(
        n, t, "L",
        mdataL, stateL,
        nearestIdxL,
        surfaceModeRmsL);

    writeModalDiagnostics(
        n, t, "R",
        mdataR, stateR,
        nearestIdxR,
        surfaceModeRmsR);

    writeTop10ModalContributions(
        n, t, "L",
        geomL, mdataL, stateL,
        nearestIdxL,
        fmodalTop10);

    writeTop10ModalContributions(
        n, t, "R",
        geomR, mdataR, stateR,
        nearestIdxR,
        fmodalTop10);
}
```

この位置の `state.q` は確定済みの \(t_n\) 状態であり、同時刻の `displace.dat`、`modal_contribution.csv` と整合する。

接触反復中の予測値 `qf` ではなく、必ず確定済みの `state.q` を使うこと。

---

## 計算結果を変えないための条件

診断追加によって次を変更してはならない。

- `state.q`
- `state.qf`
- `state.disp`
- `state.predictedDisp`
- 流体荷重
- 接触荷重
- 一般化力
- Newmark積分
- 接触反復回数
- 浮動小数点演算の順序を含む既存計算

診断処理は読み取り専用にする。

特に、既存の変位計算を診断関数の結果で置き換えないこと。

---

## 検証条件

実装後、最低限以下を自動または手動で確認する。

## 1. 既存結果の完全保持

同一パラメータで、診断追加前後の以下が一致すること。

- `displace.dat`
- `displace_xy.dat`
- `area.dat`
- `pressure.dat`
- `airflow_vt.dat`
- 接触診断
- 最終モード座標

理想的にはファイル単位で一致すること。少なくとも数値差が丸め誤差以下であること。

## 2. 大きさ寄与率の合計

各時刻、左右、scopeについて、

\[
\sum_m A_m=1
\]

が成立すること。

許容誤差例：

```text
abs(sum_magnitude_ratio - 1.0) < 1e-10
```

合成変位も全モード変位もゼロの場合は例外としてよい。

## 3. 符号付き投影寄与率の合計

合成変位ノルムがゼロでない場合、

\[
\sum_mB_m=1
\]

が成立すること。

許容誤差例：

```text
abs(sum_signed_projection_ratio - 1.0) < 1e-9
```

## 4. モード和と実変位の一致

監視点について、

\[
\sum_m
q_m\phi_{m,y}10^3
\]

と、

\[
\texttt{state.disp[probeId].uy}
-
\texttt{geom.points[probeId].y}
\]

が一致すること。

## 5. 上位10の順位

各scopeについて、

```text
rank 1 modal_norm >= rank 2 modal_norm >= ... >= rank 10 modal_norm
```

であること。

## 6. 人工的な相殺テスト

可能なら簡単な単体テストを追加する。

```text
mode 1 = +2
mode 2 = -1
total  = +1
```

のとき、

```text
magnitude ratio:
mode 1 = 2/3
mode 2 = 1/3

signed ratio:
mode 1 = 2
mode 2 = -1

signed ratio sum = 1
cancellation factor = 3
```

になること。

---

## 注意：RMSとエネルギー寄与率は別物

本診断は「変位への寄与率」であり、「振動エネルギーの寄与率」ではない。

質量正規化が適切なら、モードエネルギーは概ね、

\[
E_m
=
\frac12\dot q_m^2
+
\frac12\omega_m^2q_m^2
\]

で評価できる。

ただし、変位寄与とエネルギー寄与は一致しない。

- 低周波モードは、同じエネルギーでも変位が大きくなりやすい
- 高周波モードは、変位が小さくてもエネルギーが大きい場合がある

今回の `modal_top10.csv` は変位寄与専用とし、エネルギー割合を同じ列へ混在させないこと。

将来必要なら、別ファイルとして、

```text
modal_energy_top10.csv
```

を追加する。

---

## 現行コードで同時に確認された注意点

`src/simulation.cpp` のmanifest出力では、右側のモデルも `b12c3` と記録している。

一方、実際の右側読み込みは現時点で、

```text
M5_mode_T3_b2c2.vtu
M5_freq_T3_d2_b2c2.txt
```

となっている。

したがって、manifest上の右モデルと実際の入力が一致していない。

モード寄与率の実装とは独立した問題だが、実行条件の取り違えを防ぐため、manifestには実際に読み込む左右のファイルパスを記録すること。

可能ならファイルパスを文字列変数へまとめ、

```cpp
const fs::path leftModeFile = ...;
const fs::path rightModeFile = ...;
const fs::path leftFrequencyFile = ...;
const fs::path rightFrequencyFile = ...;
```

として、読み込みとmanifestで同じ変数を使う。

---

## 完了条件

以下をすべて満たしたら実装完了とする。

- `output/modal_top10.csv` が生成される
- 左右それぞれについて上位10モードが出力される
- `surface_xyz`、`surface_uy`、`probe_uy` が出力される
- 大きさ寄与率の全モード合計が1になる
- 符号付き投影寄与率の全モード合計が1になる
- 負の符号付き寄与率を保持する
- 合成変位ゼロ付近では `nan` を出す
- 既存のシミュレーション結果が変化しない
- 既存出力ファイルの形式を変更しない
- manifestが実際の左右入力ファイルと一致する

---

## 最終的な解釈

上位モードの評価では、次の2値を必ずセットで見る。

\[
\boxed{
A_m
=
\frac{\|\boldsymbol U_m\|}
{\sum_k\|\boldsymbol U_k\|}
}
\]

は、

> そのモードが持つ変位振幅の相対的な大きさ

を表す。

\[
\boxed{
B_m
=
\frac{
\langle\boldsymbol U_m,\boldsymbol U\rangle
}{
\langle\boldsymbol U,\boldsymbol U\rangle
}
}
\]

は、

> そのモードが最終合成変位を強めたか、打ち消したか

を表す。

上位10の順位は \(A_m\) で決め、物理的な強め合い・相殺の解釈には \(B_m\) と相殺係数 \(C\) を使うこと。
