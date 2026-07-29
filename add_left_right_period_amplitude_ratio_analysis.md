# 肺圧スイープ解析への左右周期比・振幅比の追加仕様

## 目的

現行の肺圧スイープおよび各条件の振動解析に、左右声帯の

- 周期比
- 振幅比

を追加する。

目的は、肺圧ごとに次の振動状態を自動抽出できるようにすることである。

1. 左右の代表変位の周期が約 2:1 または 1:2 となる条件
2. 左右の代表変位振幅が約 2:1 または 1:2 となる条件

既存の声門面積解析、Cycle Amplitude、Cycle Period、Poincaré Points、Return Map 等の出力は変更せず、追加解析として実装すること。

---

## 実装対象

リポジトリ内の現行スイープコードと、各圧力条件に対して解析図・CSVを作成している解析スクリプトを確認し、既存の処理フローに組み込む。

特に、現在の Poincaré Points で用いている左右代表変位

- `uyL_mm`
- `uyR_mm`

または同等の左右変位列を、そのまま周期比・振幅比の解析信号として再利用する。

列名や入力ファイル名は実装前に現行コードを確認し、既存仕様に合わせること。新しい固定列番号を無断で仮定しないこと。

---

# 1. 解析区間

## 1.1 過渡区間の除外

肺圧ランプおよび振動立ち上がりの影響を除外する。

既存解析に定常区間の開始時刻設定がある場合は、それを共通利用する。設定がない場合は、コマンドライン引数または設定値として次を追加する。

```python
ANALYSIS_START_TIME = 0.18  # [s]
ANALYSIS_END_TIME = None    # None の場合はデータ末尾
```

周期比・振幅比は、この解析区間内のデータだけから求める。

解析区間内の周期数が少なすぎる場合は、値を無理に算出せず `insufficient_data` とする。

推奨条件：

```python
MIN_PEAK_COUNT = 6
MIN_CYCLE_COUNT = 5
```

---

# 2. 入力信号の前処理

左右信号を

\[
y_L(t),\qquad y_R(t)
\]

とする。

## 2.1 平均値除去

各信号から解析区間内平均を引く。

\[
\tilde y_L(t)=y_L(t)-\overline{y_L}
\]

\[
\tilde y_R(t)=y_R(t)-\overline{y_R}
\]

周期・振幅解析には平均値除去後の信号を使用する。

## 2.2 単位

Poincaré 解析で mm に変換済みなら、その mm 単位を利用する。

変位が m で保存されている場合のみ、

```python
y_mm = y_m * 1.0e3
```

とする。

既存コードの単位を確認し、二重変換しないこと。

## 2.3 平滑化

数値ノイズによる偽ピークを避けるため、必要に応じて軽い平滑化を許可する。ただし、周期や振幅を変えるほど強いフィルタは使わない。

推奨：

- Savitzky–Golay filter
- 移動平均は短い窓に限定
- 平滑化なしも選択可能

例：

```python
SMOOTH_SIGNAL = True
SAVGOL_WINDOW = 11   # 奇数
SAVGOL_POLYORDER = 3
```

窓幅は基本周期の十分の一以下とする。サンプリング数に応じて自動調整し、不正な窓長を使わないこと。

---

# 3. 左右別ピーク検出

周期比と振幅比は、左右それぞれの代表変位のピークから求める。

```python
from scipy.signal import find_peaks
```

を使用する。

## 3.1 基本周波数の事前推定

既存の声門面積解析で基本周波数が得られている場合は、その値をピーク間隔条件の初期値に使う。

得られていない場合は、左右各信号の Welch PSD または FFT から支配周波数を推定する。

解析対象は、声帯振動として妥当な周波数範囲に限定する。

例：

```python
FREQ_MIN_HZ = 40.0
FREQ_MAX_HZ = 500.0
```

左右で周期が 2 倍異なる状態を探すため、両側に同じ周波数を強制しないこと。それぞれ独立に支配周波数を推定する。

## 3.2 `find_peaks` の条件

左右ごとに以下を設定する。

- `distance`
- `prominence`

`distance` は各側の推定周期から決める。

```python
distance_samples = int(0.6 * fs / dominant_frequency_hz)
```

`prominence` は信号振幅に対する相対値とする。

```python
prominence = PEAK_PROMINENCE_RATIO * np.ptp(signal)
```

推奨初期値：

```python
PEAK_PROMINENCE_RATIO = 0.10
```

固定値だけに依存せず、左右の振幅差が大きい場合でも小振幅側のピークを検出できるようにする。

## 3.3 ピーク検出結果の検証

次を確認する。

- ピーク時刻が単調増加している
- ピーク間隔が 0 以下でない
- 極端な外れ周期がない
- 解析区間端の不完全周期を除外する

外れ値除去には中央値を使用する。

\[
T_n=t_{n+1}-t_n
\]

として、例えば次の範囲外を除外する。

\[
0.5\,\mathrm{median}(T)
<
T_n
<
1.5\,\mathrm{median}(T)
\]

ただし、本研究では周期の不規則性そのものが重要であるため、外れ値を完全に捨てるのではなく、

- 生データ
- 外れ値除外後データ

の両方を保存する。

周期比の代表値には外れ値除外後を使い、生データの変動量もCSVに残す。

---

# 4. 周期比の定義

左右の有効周期列を

\[
T_{L,n},\qquad T_{R,n}
\]

とする。

代表周期は中央値を基本とする。

\[
\tilde T_L=\mathrm{median}(T_{L,n})
\]

\[
\tilde T_R=\mathrm{median}(T_{R,n})
\]

平均値も参考値として出力する。

## 4.1 符号付き周期比

左右の向きを残す比を

\[
R_T^{L/R}
=
\frac{\tilde T_L}{\tilde T_R}
\]

とする。

- \(R_T^{L/R}\approx2\)：左の周期が右の約2倍
- \(R_T^{L/R}\approx0.5\)：右の周期が左の約2倍
- \(R_T^{L/R}\approx1\)：左右の代表周期がほぼ同じ

## 4.2 大小関係を無視した周期比

スイープで候補を抽出しやすいように、

\[
R_T
=
\max\left(
\frac{\tilde T_L}{\tilde T_R},
\frac{\tilde T_R}{\tilde T_L}
\right)
\]

も出力する。

この値は常に \(R_T\ge1\) となる。

## 4.3 2:1 周期候補判定

初期設定：

```python
PERIOD_RATIO_TARGET = 2.0
PERIOD_RATIO_TOL = 0.20
```

候補条件：

\[
1.8 \le R_T \le 2.2
\]

ただし、比だけで誤判定しないため、次も同時に確認する。

- 左右とも有効周期数が `MIN_CYCLE_COUNT` 以上
- 左右の周期CVが過大でない
- PSDの支配周波数比も概ね 2:1 または 1:2

周期の変動係数：

\[
CV_{T,L}=
\frac{\sigma(T_L)}{\overline{T_L}}
\]

\[
CV_{T,R}=
\frac{\sigma(T_R)}{\overline{T_R}}
\]

候補判定用の上限例：

```python
MAX_PERIOD_CV_FOR_LOCKING = 0.20
```

`CV` が大きい場合は、`2to1_period_candidate` ではなく、

```text
ratio_near_2_but_irregular
```

と分類する。

---

# 5. 振幅比の定義

左右の各周期振幅は、連続するピークの間に存在する谷を利用して定義する。

左右それぞれについて、

\[
a_n =
y_{\mathrm{peak},n}
-
y_{\mathrm{trough},n}
\]

とする。

単純な正負対称振動を仮定して

\[
a_n = 2|y_{\mathrm{peak},n}|
\]

とはしないこと。左右で静的オフセットや波形非対称性が異なる可能性があるため、peak-to-trough を使用する。

## 5.1 周期振幅の求め方

各ピーク \(p_n\) に対し、前後の隣接ピークとの間の谷を検出する。

実装しやすい方法：

1. 正のピーク列を検出
2. 連続するピーク区間内の最小値を谷とする
3. 各周期の peak-to-trough 振幅を計算
4. 解析区間端の不完全周期は除外

例：

```python
cycle_amplitude = peak_value - local_trough_value
```

全振幅として peak-to-peak を使用したい場合は、設定で切り替え可能にしてもよいが、既定値は peak-to-trough とする。

## 5.2 代表振幅

左右の周期振幅中央値を

\[
\tilde a_L=\mathrm{median}(a_{L,n})
\]

\[
\tilde a_R=\mathrm{median}(a_{R,n})
\]

とする。

平均値、標準偏差、CVも保存する。

## 5.3 符号付き振幅比

\[
R_A^{L/R}
=
\frac{\tilde a_L}{\tilde a_R}
\]

- \(R_A^{L/R}\approx2\)：左振幅が右の約2倍
- \(R_A^{L/R}\approx0.5\)：右振幅が左の約2倍

## 5.4 大小関係を無視した振幅比

\[
R_A
=
\max\left(
\frac{\tilde a_L}{\tilde a_R},
\frac{\tilde a_R}{\tilde a_L}
\right)
\]

## 5.5 2:1 振幅候補判定

初期設定：

```python
AMPLITUDE_RATIO_TARGET = 2.0
AMPLITUDE_RATIO_TOL = 0.20
```

候補条件：

\[
1.8 \le R_A \le 2.2
\]

追加条件：

- 左右とも有効振幅周期数が `MIN_CYCLE_COUNT` 以上
- 小振幅側がノイズフロアより十分大きい
- 振幅比が単発の外れ周期だけで作られていない

小振幅側の代表振幅が極端に小さい場合は、比が発散するため候補から除外する。

```python
MIN_VALID_AMPLITUDE_MM = 1.0e-4
```

実際の変位スケールを確認して妥当な値へ調整すること。

---

# 6. 時間局所的な比の確認

左右で不規則性が強い条件では、全解析区間の代表値だけでなく、時間窓ごとの比も出す。

例：

```python
WINDOW_DURATION = 0.05  # [s]
WINDOW_OVERLAP = 0.5
```

各窓について、

- `period_ratio_unsigned`
- `amplitude_ratio_unsigned`

を計算する。

これにより、

- 全時間で約2倍
- 一部区間だけ約2倍
- 時間とともに比が遷移

を区別できる。

ただし、窓内周期数が不足する場合は `NaN` とする。

この時間局所解析は追加機能とし、まずは全定常区間の代表比を必須実装する。

---

# 7. 出力CSV

各圧力条件について、既存の解析結果ディレクトリへ次を保存する。

## 7.1 条件別詳細CSV

ファイル例：

```text
left_right_ratio_metrics.csv
```

1行のサマリーまたはkey-value形式ではなく、後で集計しやすい列形式とする。

推奨列：

```text
pressure_pa
analysis_start_s
analysis_end_s

left_peak_count
right_peak_count
left_cycle_count
right_cycle_count

left_period_median_s
right_period_median_s
left_period_mean_s
right_period_mean_s
left_period_std_s
right_period_std_s
left_period_cv
right_period_cv

period_ratio_left_over_right
period_ratio_unsigned
period_ratio_near_2
period_ratio_class

left_amplitude_median_mm
right_amplitude_median_mm
left_amplitude_mean_mm
right_amplitude_mean_mm
left_amplitude_std_mm
right_amplitude_std_mm
left_amplitude_cv
right_amplitude_cv

amplitude_ratio_left_over_right
amplitude_ratio_unsigned
amplitude_ratio_near_2
amplitude_ratio_class

left_dominant_frequency_hz
right_dominant_frequency_hz
frequency_ratio_left_over_right
frequency_ratio_unsigned

status
warning
```

`period_ratio_class` の例：

```text
same_period
left_period_about_twice_right
right_period_about_twice_left
ratio_near_2_but_irregular
other
insufficient_data
```

`amplitude_ratio_class` の例：

```text
same_scale
left_amplitude_about_twice_right
right_amplitude_about_twice_left
other
insufficient_data
```

## 7.2 周期ごとの詳細

ファイル例：

```text
left_cycle_metrics.csv
right_cycle_metrics.csv
```

推奨列：

```text
cycle_index
peak_time_s
period_raw_s
period_valid
peak_value_mm
trough_value_mm
amplitude_mm
```

左右別に保存する。

---

# 8. スイープ全体の集計CSV

現行スイープの各圧力条件の終了後、全結果を1つにまとめる。

ファイル例：

```text
pressure_sweep_left_right_ratios.csv
```

1圧力1行とし、最低限次を含める。

```text
pressure_pa
left_period_median_s
right_period_median_s
period_ratio_left_over_right
period_ratio_unsigned
period_ratio_class

left_amplitude_median_mm
right_amplitude_median_mm
amplitude_ratio_left_over_right
amplitude_ratio_unsigned
amplitude_ratio_class

left_period_cv
right_period_cv
left_amplitude_cv
right_amplitude_cv

left_dominant_frequency_hz
right_dominant_frequency_hz
frequency_ratio_unsigned

status
```

既存のスイープサマリーCSVがある場合は、新規CSVを増やすより既存CSVへ列追加してよい。ただし、既存列名・列順を壊さず末尾へ追加すること。

---

# 9. 追加図

各圧力条件について、最低限次の図を追加する。

## 9.1 左右変位重ね描き

ファイル例：

```text
left_right_displacement.png
```

- 横軸：Time [s]
- 縦軸：Displacement [mm]
- 左右信号を同一図に表示
- 解析区間を使用
- 左右ピーク位置をマーカー表示可能にする

この図は、ピーク数が2:1になっているかを目視確認するために必要。

## 9.2 左右周期列

ファイル例：

```text
left_right_cycle_period.png
```

- 横軸：左右それぞれの周期インデックス、またはピーク時刻
- 縦軸：Period [s]
- 左右を同一図へ表示

左右の周期系列は必ずしも同数でないため、インデックスを無理に揃えない。ピーク時刻を横軸にする方法を推奨する。

## 9.3 左右振幅列

ファイル例：

```text
left_right_cycle_amplitude.png
```

- 横軸：ピーク時刻または周期インデックス
- 縦軸：Amplitude [mm]
- 左右を同一図へ表示

## 9.4 圧力―比の分岐図

スイープ終了後に次を作る。

```text
pressure_vs_period_ratio.png
pressure_vs_amplitude_ratio.png
```

横軸：

```text
Pressure [Pa]
```

縦軸：

```text
Unsigned period ratio [-]
Unsigned amplitude ratio [-]
```

基準線として、

```text
ratio = 1
ratio = 2
```

を表示する。

候補点は色またはマーカーで区別してよいが、既存プロットのスタイルに合わせる。

---

# 10. 実装上の重要事項

## 10.1 左右信号を全体声門周期で分割しない

左右で周期が2倍近く異なる状態を探すため、左右の周期検出は独立に行う。

全体声門面積のイベント時刻を左右両方へ強制的に使用すると、2:1周期差を見逃すため禁止する。

## 10.2 FFTだけで周期比を決めない

非正弦波では高調波が主ピークになる場合がある。

周期比の主判定は時間領域のピーク間隔から行い、PSDの支配周波数比は整合性確認として使用する。

## 10.3 振幅比は同一時刻の瞬時値比ではない

左右の位相がずれている可能性があるため、

```python
abs(y_left(t)) / abs(y_right(t))
```

のような瞬時比を代表振幅比にしない。

左右それぞれの周期振幅中央値を比較する。

## 10.4 周期2の声門面積と左右周期差を混同しない

全体声門面積が大小交互になる周期2振動と、左右の基本周期が2:1になる状態は異なる。

今回追加する `period_ratio` は左右個別変位の代表周期比であり、全体声門面積の周期倍化判定とは別指標として扱う。

## 10.5 欠損・ゼロ除算

次の場合は `NaN` とし、例外でスイープ全体を停止させない。

- 入力ファイルなし
- 左右列なし
- ピーク不足
- 周期不足
- 振幅がゼロまたはノイズ以下
- 不正値・NaN・Infの混入

`status` と `warning` に理由を記録する。

---

# 11. コマンドライン引数または設定項目

現行解析コードの設定方式に合わせ、次を追加する。

```text
--analysis-start
--analysis-end
--left-column
--right-column
--peak-prominence-ratio
--period-ratio-target
--period-ratio-tolerance
--amplitude-ratio-target
--amplitude-ratio-tolerance
--min-peak-count
--min-valid-amplitude
```

既存の設定ファイルまたはスクリプト冒頭定数を使っている場合は、その形式に統一してよい。

デフォルト値を与え、従来の実行コマンドを壊さないこと。

---

# 12. 推奨関数構成

既存スクリプトの構造に合わせつつ、解析部分は再利用可能な関数へ分離する。

例：

```python
def preprocess_signal(time, signal, start_time, end_time):
    ...

def estimate_dominant_frequency(time, signal, freq_min, freq_max):
    ...

def detect_cycles(time, signal, fs, dominant_frequency, prominence_ratio):
    ...

def calculate_cycle_amplitudes(time, signal, peak_indices):
    ...

def summarize_periods(peak_times):
    ...

def summarize_amplitudes(amplitudes):
    ...

def calculate_unsigned_ratio(left_value, right_value):
    ...

def classify_period_ratio(
    ratio_left_over_right,
    unsigned_ratio,
    left_cv,
    right_cv,
    target=2.0,
    tolerance=0.2,
):
    ...

def classify_amplitude_ratio(
    ratio_left_over_right,
    unsigned_ratio,
    target=2.0,
    tolerance=0.2,
):
    ...

def analyze_left_right_ratios(time, left_signal, right_signal, config):
    ...
```

戻り値は辞書または dataclass とし、CSV出力とプロットの双方で利用する。

---

# 13. 検証

## 13.1 人工信号テスト

自動テスト基盤がない場合でも、解析関数単体を確認する短いテストスクリプトを追加する。

### 同周期・同振幅

```python
left  = np.sin(2*np.pi*100*t)
right = np.sin(2*np.pi*100*t)
```

期待値：

```text
period_ratio_unsigned ≈ 1
amplitude_ratio_unsigned ≈ 1
```

### 周期比 2:1

```python
left  = np.sin(2*np.pi*100*t)
right = np.sin(2*np.pi*50*t)
```

期待値：

```text
period_ratio_left_over_right ≈ 0.5
period_ratio_unsigned ≈ 2
```

### 振幅比 2:1

```python
left  = 2*np.sin(2*np.pi*100*t)
right = 1*np.sin(2*np.pi*100*t)
```

期待値：

```text
amplitude_ratio_left_over_right ≈ 2
amplitude_ratio_unsigned ≈ 2
```

### 周期比・振幅比とも2:1

```python
left  = 2*np.sin(2*np.pi*100*t)
right = 1*np.sin(2*np.pi*50*t)
```

期待値を確認する。

### ノイズ付き

上記に弱い乱数ノイズを加え、許容範囲内で結果が維持されることを確認する。

## 13.2 実データ検証

少なくとも次を確認する。

1. 既存の解析図・CSVが変更前と一致する
2. 左右変位プロットが現在の Poincaré 解析に使用している信号と一致する
3. ピークマーカーが目視上の振動ピークと対応する
4. 600 Pa、700 Pa、800 Paなど既存スイープ条件で比が異常値にならない
5. 入力欠損時もスイープ全体が停止しない

---

# 14. 完了条件

以下をすべて満たした時点で完了とする。

- 左右代表変位を独立にピーク検出できる
- 左右周期中央値・平均値・CVを算出できる
- 左右周期比を符号付き・unsignedの両方で出力できる
- 左右周期比約2の候補を分類できる
- 左右周期振幅を peak-to-trough で算出できる
- 左右振幅比を符号付き・unsignedの両方で出力できる
- 左右振幅比約2の候補を分類できる
- 各圧力条件の詳細CSVを出力できる
- スイープ全体の比一覧CSVを出力できる
- 圧力対周期比・振幅比の図を出力できる
- 既存解析結果を壊さない
- Releaseビルドおよび現行スイープ実行手順が維持される

---

# 15. Codexへの実装指示

1. まずリポジトリ内の現行スイープスクリプト、Poincaré解析、左右変位読み込み処理を確認する。
2. 現在の左右代表変位列と単位を特定する。
3. 上記仕様を既存設計に合わせて最小変更で追加する。
4. 重複コードを避け、周期・振幅解析を関数化する。
5. 従来の出力および数値計算部分には変更を加えない。
6. 実装後、変更ファイルと追加出力を一覧で報告する。
7. 人工信号テストと既存実データでの確認結果を報告する。
8. 周期比・振幅比の算出に使用した左右列名、解析区間、ピーク検出条件を明記する。
9. 仕様上判断が必要な点があれば、勝手に物理条件を変更せず、既存コードを優先して実装する。
