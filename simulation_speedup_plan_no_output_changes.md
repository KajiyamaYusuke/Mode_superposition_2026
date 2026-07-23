# 現行計算内容を維持した高速化方針

更新日: 2026-07-23  
対象リポジトリ: `KajiyamaYusuke/Mode_superposition_2026`  
対象ブランチ: `main`

## 1. 目的

現在の声帯振動計算について、次の条件を維持したまま計算時間を短縮する。

- 流体モデルを変更しない
- 接触モデルを変更しない
- 使用モード数を勝手に減らさない
- 時間刻みを変更しない
- 接触反復の収束条件を変更しない
- 現在の出力内容・出力頻度を変更しない
- 得られる変位、流量、圧力、接触力、VTU形状を数値誤差範囲で一致させる

したがって、本書では「物理モデルの簡略化」ではなく、主に次を対象とする。

1. 同じ演算の重複削除
2. メモリアクセスの改善
3. 並列化
4. 行列–ベクトル積への整理
5. 接触反復中に必要な節点だけを更新する処理分離

---

# 2. 現在の計算時間

現在のTiming Summaryは次の通りである。

| 処理 | 時間 [ms] | 時間 [s] | 全体に対する割合 |
|---|---:|---:|---:|
| `calcArea` | 26,741.2 | 26.7 | 1.4% |
| `calcForce` | 33,467.9 | 33.5 | 1.7% |
| `f2mode` | 626,034 | 626.0 | 32.7% |
| `mode2uf` | 737,936 | 737.9 | 38.5% |
| `calcDis` | 135,095 | 135.1 | 7.0% |
| `output` | 357,837 | 357.8 | 18.7% |
| **合計** | **1,917,111** | **1,917.1** | **100%** |

最も大きい二つは、

\[
\texttt{f2mode}+\texttt{mode2uf}
=
1,363,970\ {\rm ms}
\]

であり、全計算時間の約71.1%を占める。

したがって、高速化の優先順位は次の通りである。

1. `mode2uf`
2. `f2mode`
3. 接触反復中の重複演算
4. `calcDis`
5. `calcArea`・`calcForce`

`calcArea`と`calcForce`は合計でも約3.1%であるため、最初に細かく最適化しても全体への効果は小さい。

---

# 3. 現行実装のボトルネック

## 3.1 `mode2uf`: 全節点 × 全モードの重ね合わせ

現在の `State::mode2uf()` は、全体積節点について全モードを重ね合わせている。

概略:

```cpp
for (int point = 0; point < nPoints; ++point) {
    for (int mode = 0; mode < nModes; ++mode) {
        predictedDisp[point] += modeShape[mode][point] * qf[mode];
        velocity[point]      += modeShape[mode][point] * qfdot[mode];
    }
}
```

計算量は片側あたり概ね、

\[
O(N_{\mathrm{point}}N_{\mathrm{mode}})
\]

である。

さらに、1時間ステップ内で、

- 左右の各接触反復
- 最終接触力を反映する再計算

のたびに呼ばれる。

接触判定と流路断面計算が実際に参照するのは主に内側表面節点であるため、接触反復中に全体積節点を毎回復元するのは大きな無駄になっている可能性が高い。

---

## 3.2 `f2mode`: 表面節点 × 全モードの逐次射影

現在の `ModalProjector::project()` は、各モードについて全表面節点を走査し、

\[
f_m
=
\sum_i
\boldsymbol F_i\cdot\boldsymbol\phi_{m,i}
\]

を計算している。

概略:

```cpp
for (int mode = 0; mode < nModes; ++mode) {
    for (int i = 0; i < nsurfl; ++i) {
        for (int j = 0; j < nsurfz; ++j) {
            force[mode] += load[i][j] * modeShape[mode][point];
        }
    }
}
```

現行コードでは、モード外側ループにOpenMP並列化が入っていない。

各モードの一般化力は独立に計算できるため、この部分は比較的安全に並列化できる。

---

## 3.3 モード配列のメモリアクセス方向が処理ごとに異なる

現在のモード配列は概ね、

```text
modes[mode][point]
```

の順で保持されている。

この配置は、

```text
for mode:
    for point:
```

と走査する `f2mode` には比較的向いている。

一方、`mode2uf` は、

```text
for point:
    for mode:
```

と走査するため、内側ループで離れたメモリ領域を読みやすい。

結果として、

- キャッシュミス
- メモリ帯域待ち
- SIMDベクトル化の阻害

が起こる可能性がある。

---

## 3.4 接触反復ごとに同じ固定量を再評価している

1時間ステップ内では流体荷重が固定され、接触荷重だけが反復更新される。

それでも毎反復で、

- 流体荷重と接触荷重を合成
- 全表面荷重を全モードへ射影
- 左右全節点変位を復元

している。

このうち流体一般化力は、その時間ステップ内では一定である。

したがって、

\[
f_m^{\mathrm{total}}
=
f_m^{\mathrm{fluid}}
+
f_m^{\mathrm{contact}}
\]

として分離すれば、毎反復で流体荷重を再射影する必要がない。

---

# 4. 優先度P0: 低リスクで先に行う高速化

## P0-1. `ModalProjector::project()` のモードループをOpenMP並列化

### 修正方針

左右それぞれについて、外側のモードループを並列化する。

```cpp
#pragma omp parallel for schedule(static)
for (int modeIndex = 0; modeIndex < leftModes.nModes; ++modeIndex) {
    double sum = 0.0;

    for (int i = 0; i < leftGeom.nsurfl; ++i) {
        for (int j = 0; j < leftGeom.nsurfz; ++j) {
            const int id = leftGeom.surfp[i][j];
            if (id < 0) continue;

            const auto& phi = leftModes.modes[modeIndex][id];

            sum += load.fxL[i][j] * phi.ux
                 + load.fyL[i][j] * phi.uy
                 + load.fzL[i][j] * phi.uz;
        }
    }

    forceL[modeIndex] = sum;
}
```

右側も同様にする。

### 安全性

各スレッドは異なる `forceL[modeIndex]` または `forceR[modeIndex]` へ書き込むため、モード間のデータ競合はない。

ループ内ではローカル変数 `sum` を使い、配列要素へ繰り返し加算しないほうがよい。

### 注意

左右のループをそれぞれ別の `parallel for` にすると、短いループではスレッド生成・同期コストが増える可能性がある。

次の二つを計測して比較する。

#### 方法A: 左右別parallel for

実装が単純で安全。

#### 方法B: 左右を一つのparallel regionに入れる

```cpp
#pragma omp parallel
{
    #pragma omp for nowait schedule(static)
    for (...) {
        // left
    }

    #pragma omp for schedule(static)
    for (...) {
        // right
    }
}
```

### 合格条件

- 逐次版と各モード一般化力が十分小さい丸め誤差内で一致する
- 左右の最大一般化力モード番号が一致する
- 変位・流量・接触結果が一致する
- `f2mode`時間が短縮する

---

## P0-2. 固有角振動数とNewmark定数を事前計算

現在、時間ステップ・接触反復・モードループ内で、

\[
\omega_m=2\pi f_m
\]

を毎回計算している。

またNewmark法では、各モードについて時間刻みと固有角振動数から決まる固定係数が存在する。

### 修正方針

初期化時に左右それぞれ、

```cpp
omega[m] = 2.0 * M_PI * frequency[m];
omega2[m] = omega[m] * omega[m];
```

を保存する。

減衰比、時間刻み、Newmarkの \(\beta,\gamma\) が計算中に一定なら、分母などの固定係数も事前計算できる。

例:

```cpp
struct NewmarkModeCoefficients {
    double omega;
    double omega2;
    double damping;
    double effectiveDenominator;
    // 必要な固定係数
};
```

### 効果

個々の削減量は小さいが、全モード・全時間ステップ・全接触反復で繰り返されるため、低コストで実施できる。

### 合格条件

- Newmark更新値が現行実装と丸め誤差内で一致
- 固有周波数・減衰を途中変更するコードがないことを確認

---

## P0-3. `force.assign()` の繰り返し再確保を避ける

現在の `ModalProjector::project()` は毎回、

```cpp
forceL.assign(leftModes.nModes, 0.0);
forceR.assign(rightModes.nModes, 0.0);
```

を実行する。

`assign()` はサイズによっては再確保や初期化コストを生む。

### 修正方針

初期化時に必要サイズを確保し、時間ループでは、

```cpp
std::fill(forceL.begin(), forceL.end(), 0.0);
std::fill(forceR.begin(), forceR.end(), 0.0);
```

を使う。

並列化後に各要素を必ず上書きするなら、ゼロ初期化自体も不要にできる。

```cpp
forceL[modeIndex] = sum;
```

### 合格条件

- ベクトルサイズが計算中に変化しない
- 各モード要素を必ず上書きする

---

## P0-4. OpenMPスレッド設定を固定して比較する

現在CMakeではOpenMPとRelease `-O3` が有効である。

実行時に次を設定する。

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

スレッド数は物理コア数を中心に比較する。

```bash
OMP_NUM_THREADS=2  ./build/simulation
OMP_NUM_THREADS=4  ./build/simulation
OMP_NUM_THREADS=6  ./build/simulation
OMP_NUM_THREADS=8  ./build/simulation
OMP_NUM_THREADS=12 ./build/simulation
```

### 注意

`mode2uf` はメモリ帯域律速になり得るため、スレッド数を増やしすぎると速くならない、または遅くなる可能性がある。

### 記録項目

- 全体時間
- `f2mode`
- `mode2uf`
- `calcDis`
- CPU使用率
- メモリ使用量

---

# 5. 優先度P1: 効果が大きい構造的高速化

## P1-1. 表面節点専用の `mode2ufSurface()` を追加

### 目的

接触反復中は、流体断面・流体力・接触判定に必要な節点だけを復元する。

全体積節点の復元は、現在と同じ出力を行う時刻や、全節点状態が必要な処理の直前だけ実施する。

出力内容・頻度は変更しない。

### 現行

```cpp
State::mode2uf(...)
```

が全節点を更新する。

### 推奨構造

```cpp
void State::mode2ufSurface(
    const Geometry& geom,
    const ModeData& modeData
);

void State::mode2ufAll(
    const Geometry& geom,
    const ModeData& modeData
);
```

### 表面節点リストの事前作成

初期化時に、`surfp[i][j]` から有効な節点IDを抽出し、重複を除く。

```cpp
std::vector<int> surfacePointIds;
```

重複除去例:

```cpp
std::vector<char> used(nPoints, false);

for (int i = 0; i < nsurfl; ++i) {
    for (int j = 0; j < nsurfz; ++j) {
        int pid = surfp[i][j];

        if (pid >= 0 && !used[pid]) {
            used[pid] = true;
            surfacePointIds.push_back(pid);
        }
    }
}
```

### 接触反復中

```cpp
stateL.mode2ufSurface(...);
stateR.mode2ufSurface(...);
```

### 全節点復元が必要な場合

現在のVTU出力内容・頻度を維持するため、VTUを書き込む直前に、

```cpp
stateL.mode2ufAll(...);
stateR.mode2ufAll(...);
```

を実行する。

ただし、`state.disp` の確定処理との関係を整理する必要がある。

### 重要な確認

実装前に、接触反復から次の出力までの間に、体積内部節点の `predictedDisp` または `disp` を読む処理がないことを確認する。

少なくとも次の処理は主に表面節点を使う。

- 流路断面
- 流体表面荷重
- 接触探索
- 監視点変位
- 表面モード診断

一方、VTUは全節点を必要とする。

### 注意: `uf2u()`との関係

現在の `uf2u()` は、

```cpp
disp = predictedDisp;
```

で全節点をコピーする。

表面だけ更新した状態でそのまま実行すると、内部節点は古い `predictedDisp` をコピーする。

したがって、次のどちらかに変更する。

#### 方法A: `uf2uSurface()`を追加

```cpp
void State::uf2uSurface() {
    q = qf;
    qdot = qfdot;
    qddot = qfddot;

    for (int pid : surfacePointIds) {
        disp[pid] = predictedDisp[pid];
    }
}
```

VTU出力前に全節点を復元して `disp` を更新する。

#### 方法B: 節点座標を毎時刻確定状態として保存しない

モード状態 \(q,\dot q,\ddot q\) を真の確定状態とし、節点座標は必要時に復元するキャッシュとみなす。

長期的には方法Bが整理しやすい。

### 期待効果

全節点数に対して表面節点数が小さいほど効果が大きい。

概算比:

\[
r=
\frac{N_{\mathrm{surface}}}{N_{\mathrm{all}}}
\]

接触反復中の`mode2uf`演算量は、おおよそこの比率まで減少する可能性がある。

実際には、

- OpenMPオーバーヘッド
- キャッシュ
- VTU出力前の全節点復元
- 左右別処理

があるため単純比例にはならない。

### 合格条件

- 同じ \(q_m\) から表面節点座標が従来版と一致
- 断面積、圧力、接触力が一致
- VTU出力時の全節点座標が一致
- `mode2uf`時間が明確に短縮

---

## P1-2. 表面モード形状を連続配列へコピー

現在の、

```text
modes[mode][globalPoint]
```

を時間ループ内で直接参照すると、表面節点IDが非連続な場合にメモリアクセスが散らばる。

### 修正方針

初期化時に、表面節点だけのモード形状を連続配列へコピーする。

候補:

```cpp
phiSurfaceX[mode * nSurface + surfaceIndex]
phiSurfaceY[mode * nSurface + surfaceIndex]
phiSurfaceZ[mode * nSurface + surfaceIndex]
```

これは`f2mode`の、

```text
for mode:
    for surfacePoint:
```

に向く。

`mode2ufSurface`向けに転置配列も用意できる。

```cpp
phiSurfaceX_T[surfaceIndex * nModes + mode]
phiSurfaceY_T[surfaceIndex * nModes + mode]
phiSurfaceZ_T[surfaceIndex * nModes + mode]
```

### メモリとの交換条件

転置配列を両方持つとメモリ量は増えるが、表面節点だけであれば全体モード配列に比べて小さい可能性が高い。

### 合格条件

- 元のモード配列からコピーした値が一致
- 初期化後にモード形状が変更されない
- 表面節点順序が荷重配列の順序と一致

---

## P1-3. 表面荷重を一次元ベクトルへ平坦化

現在の荷重は、

```text
fx[i][j], fy[i][j], fz[i][j]
```

で保持されている。

`f2mode`では三重ループの中で、

- `surfp[i][j]`
- 二重vector
- モード配列

を繰り返し参照する。

### 修正方針

初期化時に表面自由度の順序を固定する。

```text
surfaceIndex = 0 ... nSurface-1
```

各時刻で、

```cpp
loadVector[3*s + 0] = Fx;
loadVector[3*s + 1] = Fy;
loadVector[3*s + 2] = Fz;
```

へ変換する。

モード行列も同じ順序にする。

\[
\boldsymbol f_q
=
\Phi_{\mathrm{surface}}^\mathsf{T}
\boldsymbol F_{\mathrm{surface}}
\]

### 効果

- 分岐削減
- 二重vectorアクセス削減
- 連続メモリアクセス
- BLAS化しやすい

---

## P1-4. `f2mode`を流体成分と接触成分に分離

### 現状

各接触反復で、流体力と接触力を合計した表面荷重全体をモードへ射影する。

しかし1時間ステップ内では流体荷重は固定である。

### 修正方針

時間ステップ冒頭に一度だけ、

\[
f_m^{\mathrm{fluid}}
=
\Phi^\mathsf{T}F^{\mathrm{fluid}}
\]

を計算する。

接触反復では、

\[
f_m^{\mathrm{contact},k}
=
\Phi^\mathsf{T}F^{\mathrm{contact},k}
\]

だけを計算し、

\[
f_m^k
=
f_m^{\mathrm{fluid}}
+
f_m^{\mathrm{contact},k}
\]

とする。

### 必要な荷重管理

現在の`SurfaceLoad`を概念的に分ける。

```cpp
SurfaceLoad fluidLoad;
SurfaceLoad contactLoad;
```

または一般化力だけを分ける。

```cpp
std::vector<double> fluidModalForceL;
std::vector<double> fluidModalForceR;
std::vector<double> contactModalForceL;
std::vector<double> contactModalForceR;
```

### 効果

接触反復数が多いほど効果が大きい。

現在の完全荷重射影を毎反復行う必要がなくなる。

### 物理内容

線形射影なので、

\[
\Phi^\mathsf{T}
(F^{\mathrm{fluid}}+F^{\mathrm{contact}})
=
\Phi^\mathsf{T}F^{\mathrm{fluid}}
+
\Phi^\mathsf{T}F^{\mathrm{contact}}
\]

が厳密に成立する。

したがって、浮動小数点の加算順序による微小差を除けば、物理内容は変わらない。

### 合格条件

- 合算後の一般化力が従来版と一致
- 接触なし時の一般化力が完全一致
- 接触反復ごとの一般化力履歴が一致

---

# 6. 優先度P2: 行列–ベクトル積化

## P2-1. `mode2ufSurface`をGEMVとして実装

表面節点の変位復元は、

\[
\boldsymbol u_x=\Phi_x\boldsymbol q,
\]

\[
\boldsymbol u_y=\Phi_y\boldsymbol q,
\]

\[
\boldsymbol u_z=\Phi_z\boldsymbol q
\]

である。

速度も、

\[
\boldsymbol v_x=\Phi_x\dot{\boldsymbol q}
\]

のように同じ行列を使う。

これは密行列–ベクトル積である。

### 候補

- oneMKL
- Eigen
- Intel compilerの最適化を利用した自前連続配列ループ

現在`icpx`を使っているため、oneMKLとの相性がよい。

### 注意

小規模な行列ではBLAS呼出しコストが支配することがある。

そのため、

- 自前OpenMP/SIMD
- Eigen
- oneMKL

を同一ケースで比較する。

---

## P2-2. `f2mode`を転置GEMVとして実装

\[
\boldsymbol f_q
=
\Phi^\mathsf{T}\boldsymbol F
\]

をBLASのGEMVへ置き換える。

### データ構造

左右それぞれについて、

```text
PhiSurface:
rows    = 3 × nSurface
columns = nModes
```

またはその転置を連続配列で保持する。

### 注意

`f2mode`と`mode2uf`で最適な行列向きが逆になる。

BLASライブラリの転置指定を使えば、同じ行列から両方を計算できる。

---

## P2-3. `restrict`・SIMDを検討

自前ループを残す場合は、コンパイラが別名参照を疑わないように連続配列化する。

必要なら、

```cpp
#pragma omp simd reduction(+:sum)
```

を表面節点ループへ追加する。

ただし、二重vectorや構造体配列のままではSIMD化されにくい。

---

# 7. 接触反復の重複計算削減

## P1-5. 非接触かつ荷重変更なしの場合の最終再計算を省略

### 現状

接触反復が非接触で終了した場合でも、最終接触力との整合性を保証するため、

1. 再度`f2mode`
2. 再度Newmark
3. 再度`mode2uf`

を実行する。

この最終再計算は、接触力が変更された場合には必要である。

一方、その時間ステップで最初から接触力がゼロであり、反復前後で荷重が変化していない場合は、反復内で得た解と最終解が同じである。

### 安全な省略条件

単に、

```cpp
if (!contactFlag)
```

では不十分である。

接触が一度発生した後に消失した場合は、接触力除去後の再計算が必要だからである。

次をすべて満たす場合だけ省略する。

- その時間ステップで接触ペアが一度も検出されていない
- 反復開始時の接触荷重がゼロ
- 反復終了時の接触荷重がゼロ
- `max_force_diff == 0`
- 現在の`qf`が同じ流体荷重から計算されている

例:

```cpp
const bool canReuseTrialSolution =
    !contactEverDetectedThisStep &&
    contactLoadWasZeroAtStart &&
    contactLoadIsZeroAtEnd &&
    max_force_diff == 0.0;
```

### 合格条件

- 再計算あり版と`qf`, `qfdot`, `qfddot`が一致
- 非接触ケースのみに適用
- 接触消失ケースには適用しない

---

## P1-6. 接触反復内でNewmarkの固定項を再利用

各接触反復は同じ確定状態、

\[
q_n,\dot q_n,\ddot q_n
\]

から、異なる外力 \(f^{(k)}\) に対して \(q_{n+1}^{(k)}\) を計算している。

線形モード方程式と固定Newmark係数では、更新式を概念的に、

\[
q_{n+1}^{(k)}
=
q_{\mathrm{base}}
+
C_f f^{(k)}
\]

の形に整理できる。

### 修正方針

モードごとに、

- 確定状態に由来する固定項
- 外力に掛かる係数

を時間ステップ冒頭に計算する。

接触反復中は外力依存部分だけを更新する。

### 注意

現在の`newmarkStep()`式を正確に展開し、同じ数式であることを確認してから行う。

### 効果

接触反復数が多い場合に有効。

---

# 8. `calcDis`の高速化候補

`calcDis`は約7%であり、`f2mode`・`mode2uf`の次に検討する。

## 8.1 接触候補キャッシュの再利用範囲

現行コードでは、同一時間ステップの2回目以降の接触反復で候補トポロジーを再利用している。

これは有効な最適化である。

追加で検討できるのは、前時間ステップの候補を初期候補として使う方法である。

ただし大変形時には接触位置が急変する可能性があるため、

- 前回候補を先に試す
- その後に完全探索で漏れを検証

という段階的方式が安全である。

## 8.2 一時vectorの再利用

各`j`ループで作られる、

- `segL`
- `segR`
- `indexL`
- `indexR`

などについて、毎回の動的確保を避ける。

`reserve()`だけでなく、必要最大サイズのバッファをメンバとして保持する方法を検討する。

## 8.3 接触デバッグ無効時の分岐整理

本計算でデバッグファイルが無効なら、

- 最大値記録
- worst pair記録
- monitor pair記録

の一部をコンパイル時オプションで外すことができる。

ただし現在の出力内容を維持する必要がある場合、必要な診断は残す。

---

# 9. コンパイル設定

現在のRelease設定は、

```cmake
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
```

である。

## 推奨候補

使用CPUが実行環境で固定される場合、

```cmake
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -march=native -mtune=native"
)
```

を比較する。

Intel `icpx`では、

```text
-xHost
```

も候補になる。

ただし、別CPUへ実行ファイルを持ち運ぶ場合は互換性が低下する。

## 浮動小数点オプション

`-ffast-math`や同等オプションは、NaN、無限大、演算順序、丸めを変更する可能性がある。

接触・流体連成の分岐挙動へ影響し得るため、最初は使用しない。

まず、

- `-O3`
- CPU命令最適化
- OpenMP
- メモリ配置改善

だけで高速化する。

---

# 10. 高速化後の数値検証

高速化では計算内容を変えないことが目的である。

各段階で、同じ短時間ケースを旧版・新版で実行して比較する。

## 10.1 比較対象

- `displace.dat`
- `displace_xy.dat`
- `area.dat`
- `pressure.dat`
- `airflow_vt.dat`
- `separation.dat`
- `contact_iteration_debug.csv`
- `modal_contribution.csv`
- VTUの点座標

## 10.2 誤差指標

時系列 \(x_i\) と高速化版 \(\hat{x}_i\) に対して、

最大絶対誤差:

\[
E_{\infty}
=
\max_i|x_i-\hat{x}_i|
\]

相対RMS誤差:

\[
E_{\mathrm{RMS}}
=
\frac{
\sqrt{
\frac{1}{N}
\sum_i(x_i-\hat{x}_i)^2
}
}{
\max\left(
\sqrt{\frac{1}{N}\sum_i x_i^2},
\epsilon
\right)
}
\]

を計算する。

## 10.3 注意

OpenMP並列化やBLAS化では加算順序が変わるため、ビット単位で完全一致しない場合がある。

非線形接触系では微小丸め差が長時間後に増幅することがある。

したがって、最初は短時間ケースで、

- 一般化力
- 1ステップ更新
- 接触反復履歴

を直接比較する。

その後、長時間ケースで、

- 基本周波数
- 振幅
- 接触時間比
- 最大めり込み
- 平均流量
- 分岐状態

が一致するか確認する。

---

# 11. 期待できる短縮量

現在の全体時間は約31.95分である。

出力時間は変更しないため、最低でも約5.96分相当はそのまま残る。

## 11.1 `f2mode`と`mode2uf`が2倍高速化した場合

現在:

\[
T_{\mathrm{modal}}
=
626.0+737.9
=
1363.9\ {\rm s}
\]

2倍高速化後:

\[
681.95\ {\rm s}
\]

全体概算:

\[
1917.1-681.95
=
1235.2\ {\rm s}
\approx20.6\ {\rm min}
\]

## 11.2 `f2mode`と`mode2uf`が3倍高速化した場合

削減量:

\[
1363.9
-
\frac{1363.9}{3}
=
909.3\ {\rm s}
\]

全体概算:

\[
1917.1-909.3
=
1007.8\ {\rm s}
\approx16.8\ {\rm min}
\]

## 11.3 追加で接触反復の重複計算を削減した場合

非接触ステップや接触反復数によるが、さらに数分程度短縮する可能性がある。

実際の値は、

- 全節点数
- 表面節点数
- モード数
- 平均接触反復数
- CPUコア数
- メモリ帯域

に依存する。

---

# 12. 推奨実装順序

## 第1段階: 小変更・低リスク

1. `ModalProjector::project()` のモードループをOpenMP並列化
2. ローカル`sum`を使って一般化力を計算
3. ベクトルの再確保を避ける
4. \(\omega,\omega^2\) とNewmark固定係数を事前計算
5. OpenMPスレッド数とaffinityを比較
6. `-march=native`または`-xHost`を比較

## 第2段階: 表面専用復元

7. 表面節点IDの一意リストを作成
8. `mode2ufSurface()`を実装
9. `mode2ufAll()`と役割を分離
10. `uf2u()`の表面更新・全体更新を整理
11. 現在と同じVTU出力時には全節点を復元
12. 旧版との時系列・VTU比較

## 第3段階: 荷重射影の整理

13. 表面モード形状を連続配列化
14. 表面荷重を一次元化
15. 流体一般化力を時間ステップ内で再利用
16. 接触一般化力だけを反復更新
17. 非接触・荷重不変時の最終再計算省略

## 第4段階: BLAS化

18. `mode2ufSurface`をGEMV化
19. `f2mode`を転置GEMV化
20. 自前OpenMP版・Eigen版・oneMKL版を比較
21. 最速かつ数値結果が安定する実装を採用

---

# 13. 最初に実装すべき三項目

現状のTiming Summaryから、最も費用対効果が高いのは次の三つである。

## 1位: `mode2ufSurface()`の追加

接触反復中に全体積節点を復元しない。

## 2位: `f2mode`のOpenMP並列化

現行の逐次モード射影を並列化する。

## 3位: 流体一般化力の時間ステップ内再利用

接触反復では接触力の射影だけを更新する。

この三つは、流体モデル、接触則、時間積分条件、モード数、出力内容を変えずに、大きな短縮が期待できる。
