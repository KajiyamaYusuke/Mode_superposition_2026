# Codex実装指示：固定x断面のz節点交点判定を廃止し、区間積分に一本化

対象：`simulation-main`（`KajiyamaYusuke/Mode_superposition_2026`）  
対象コード：`src/fluid/channelSections.cpp`  
作成日：2026-10-20

## 目的・今回の範囲

`envelopeAtKnot()` の `three or more distinct intersections` 誤検出を、**断面積の幾何・物理式を変えずに**直す。

現象：
```text
[ChannelSections] time_s=0.01576 step=1576:
section=31 x_mm=10.550000:
three or more distinct intersections at z=0.231212
```
`envelopeAtKnot()` は `z=z_k` の節点座標で、`|b.z-a.z|<=cutTolerance` のほぼ水平な線分について両端の `y` を別々に `ys` へ追加している。ほぼ水平でも、検査zとの交点が厳密には1点の場合があり、水平線分の両端とほかの交点を合わせて3個以上と誤認できる。最小再現では同種の誤検出を確認済み。ただし上記停止が同じ原因かは未確定。

**ユーザーの合意方針：面積計算に必要なのは各z積分区間における左右のy座標であり、z節点で「交点の総数」を数える必要はない。** 区間内部で選んだ線分から区間両端のyを直接求め、そのままA/G3へ渡す。

今回はこのバグだけ修正する。固定AP中央帯への変更、AP被覆率90%緩和、流体抵抗・圧力式変更、減衰/接触変更、全面リファクタリング等は**実施しない**。

## 1. 実コードを最初に確認

- `AGENTS.md`、`git status --short`、対象ファイルの現状と `git diff` を読む。既存の未コミット変更を消さない。GitHub上の仕様書とローカルのコードが違う場合は**ローカルの現状を優先して読み取り、衝突する編集は行わない**。
- 現行実装の関連部分はおおむね次の通り（行番号は変化し得る）：
  - `SurfaceCuts::build()`：`x=const`で三角形表面を切断して`Contour`を作る。
  - `validateContour()`：真の自己交差・分岐・複数連結成分を検出。
  - `envelopeAtKnot()`：z節点で水平辺の両端などを列挙して`ys.size()>2`を判定。**修正対象**。
  - `envelopeEdge()`：z区間内部のmidで選択する内側包絡線の辺IDを返す。**こちらを残す**。
  - `yOnEdge()`：辺上でyを線形補間する。**こちらを使う**。
  - `positiveLinearIntegrals()`：符号付きギャップの正部分をA/G3へ積分。**変更しない**。
  - `ChannelSectionBuilder::build()` の`knots`処理で、`leftAtKnots`・`rightAtKnots`を求め、各区間端点の`yOnEdge`結果と比較している。**この節点由来の判定を撤去**。

## 2. 必須の変更

### A. 節点での点交差判定を面積計算経路から削除

`ChannelSectionBuilder::build()` 内の以下を除去：

- `leftAtKnots` / `rightAtKnots` とそれらへの`envelopeAtKnot(...)`呼び出し
- `leftAtKnots[j-1]`、`leftAtKnots[j]`、`rightAtKnots[j-1]`、`rightAtKnots[j]`との比較による `discontinuous inner envelope / branch switch` のエラー判定

`envelopeAtKnot()` が他に使われなくなる場合は関数自体も削除する。似た機能の別関数や複雑な点交差特例処理を追加しない。

### B. z区間の内点から積分用のy両端値を得る

従来の`knots`生成、AP被覆率算定、`validateContour()`、`envelopeEdge()`は基本的に維持。

`knots[j-1]=z0`, `knots[j]=z1` の各開区間（`z1>z0`）に対し：

```cpp
const double mid = 0.5 * (z0 + z1);
const int il = envelopeEdge(leftContour,  mid, maxLeft);
const int ir = envelopeEdge(rightContour, mid, !maxLeft);

if (il < 0 || ir < 0) {
    // 左右の一方でも輪郭が存在しない：欠損区間。橋渡ししない。
    continue;
}

// 選択した線分の片側極限（同じ辺の延長上の値）として端点を評価する。
// z0/z1が辺のz区間に含まれることを通常の丸め誤差の範囲で確認する。
const double yl0 = yOnEdge(leftContour,  il, z0);
const double yl1 = yOnEdge(leftContour,  il, z1);
const double yr0 = yOnEdge(rightContour, ir, z0);
const double yr1 = yOnEdge(rightContour, ir, z1);

const double h0 = sections.gapSign * (yr0 - yl0);
const double h1 = sections.gapSign * (yr1 - yl1);
positiveLinearIntegrals(h0, h1, z1 - z0, area, cubed);
```

実装上の条件：

- `envelopeEdge()` の内点での**実際の3本以上の異なる交差枝**のエラー判定は維持する。単純にすべての`ys.size()>2`を消してはいけない。
- `validateContour()` の自己交差・分岐・複数独立成分チェックも残す。
- `yOnEdge()` は `b.z-a.z` で除算するため、選択辺が`mid`を厳密にまたぎ、分母が0でないことを確認する。完全水平辺は区間内部の候補にならない。z区間長が丸め誤差程度で不適切なら、非有限値のまま続けず明示的に扱う。
- `z0,z1` が選択辺のz範囲を浮動小数丸め誤差を超えてはみ出す場合は、無制限の外挿をしない。`knots`の結合許容差を点検し、問題区間があればx・z・左右・辺ID・端点を添えて停止する。許容差を大幅に増やしてごまかさない。
- 完全水平の輪郭辺が`z=z_k`に存在しても、それはz積分上の**幅0の境界情報**であり、両端のyを「独立した2交点」と数えて停止しない。
- もし境界の左右で選択辺が変わっても、各開区間はそれぞれの辺の**片側極限**で積分する。z節点の`y`を唯一の値に無理に決めず、今回のバグを再導入するような「節点の全交点数」検査を加えない。
- 非有限y、ゼロ除算、有限z区間での本当の多重枝や自己交差等、幾何を定義できない場合は明確に停止する。両側の輪郭があるのに前回のyを使って補完しない。

### C. 数値積分式は変えない

```text
h(z)=gapSign*(y_right(z)-y_left(z))
A(x)=∫max(h(z),0) dz         [mm²]
G3(x)=∫max(h(z),0)^3 dz      [mm⁴]
```

`positiveLinearIntegrals()` のゼロ交差分割、正ギャップの積分、単位変換、流量・圧力・剥離設定、`Rv`の剥離位置までの積分を変更しない。

`x=const`三角形切断も、81断面（0〜10 mmを0.5mm、10〜13 mmを0.05mm）、AP被覆率90%・入口から連続した有効断面、`valid`判定も今回変更しない。

## 3. 短い検証（重い連成計算は禁止）

変更前後の`src/fluid/channelSections.cpp`だけをレビューする。新しいテストフレームワークや永続的な診断ファイルは作らない。可能なら`/tmp`等の一時スモークテストで次の条件を検査：

1. 通常の輪郭線（各zで1交点）から妥当な両端yが求まる。
2. `z=z_k`に水平辺が存在する場合でも、各z**開区間**で輪郭が一意なら3交点エラーにならない。
3. `|Δz|<=cutTolerance` のほぼ水平な辺による従来の最小再現ケースで偽の3交点エラーが消える。
4. 区間**内点**で真に3つ以上の異なる交点を持つ場合は、`envelopeEdge()`が従来どおり停止する。自己交差・分岐は`validateContour()`が停止する。
5. 単純な矩形ギャップ（幅`H`、AP長`L`）なら`A=H*L`、`G3=H³*L`となる。負ギャップ区間では正部分のみを積分する。
6. 欠損AP区間は補間せず、積分対象外とする。既存の被覆率判定は維持する。

通常のReleaseビルド：`cmake --build build -j`（未構成時のみ既存CMake手順）。ビルドディレクトリからの`./simulation ../input/param.txt`起動方法を変更しない。

**0.16秒等の通常フル計算、1000Pa計算、感度解析はCodexが勝手に実行しない。** 実データで`t=0.01576`の停止が解消されたかは、ユーザー側の再実行が必要ならその旨を明示する。

## 4. クリーンな変更管理

- 変更ファイルは原則 **`src/fluid/channelSections.cpp` の1つだけ**。
- `include/`、`input/`、`CMakeLists.txt`、`flowModel.cpp`、`forceCalculatorFluid.cpp`、構造/接触/モード計算は変更しない。
- 未コミット変更・生成済み出力・Git履歴を保護。`git clean/reset/stash` を勝手に実行しない。`output/`、`result/`へ新しいファイルを出さない。
- 余分なログを恒常的に追加しない。エラーメッセージ改善は本当に必要な範囲に限定。
- 変更後、`git diff --check` と該当ファイルのdiffを確認し、変更理由・ビルド結果・行った最小テスト・未検証点だけ短く報告する。
- **実際の輪郭多価やAP被覆90%による停止が残ったら、それを別件として報告する。判定を無断で緩めない。**

## 完了条件

1. 正常なz区間で、節点に重なる水平辺の両端を2つの交点と誤認して停止しない。
2. 区間の内側で選んだ実際の輪郭辺から`yl0,yl1,yr0,yr1`を直接得てA/G3へ積分する。
3. 真の多重交点・自己交差・分岐は黙殺しない。
4. 流体の物理式・モデル係数・AP被覆率・実行コマンド・出力ファイル形式に変更がない。
5. 変更は小さく、ビルドが通り、未コミットのユーザー作業を保護できている。

この文書は**Codexへの実装指示**であり、作成時点では物理コードを修正したものではない。
