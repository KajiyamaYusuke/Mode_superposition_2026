import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt

# --- 1. データの読み込み ---
harea = np.loadtxt("../output/area.dat")
# steps列は20飛ばしになっているため時間計算には使いません
areas = harea[:, 1:]

# --- 2. 最小断面積の時系列 ---
row_min = np.min(areas, axis=1)

# --- 3. 正しい時間軸の作成 ---
# 1行が進むごとに 0.0002秒 経過すると仮定
dt = 0.0002
fs = 1 / dt  # 5000 Hz
N_total = len(row_min)
t = np.arange(N_total) * dt  # 行数基準で時間を再構築

# --- 4. 分析区間の抽出 ---
# データ後半の安定した区間 (0.15s - 0.30s) を使用
mask = (t >= 0.1) & (t <= 0.30)
y = row_min[mask]
t_seg = t[mask]

# データ点数の確認 (これ重要です！)
print(f"解析データ点数: {len(y)} 点") 
# → 正しくは750点前後になります。30〜40点しかない場合は設定ミスです。

# --- 5. DC成分除去 ---
y = y - np.mean(y)

# 目標の分解能（Hz）
target_resolution = 0.1

# 必要なデータ点数 = サンプリング周波数 / 目標分解能
# または、必要な時間長 / サンプリング間隔
# dt = 0.0002 [s]
N_padded = int(1 / (dt * target_resolution)) 
# 例: 1 / (0.0002 * 0.1) = 50000点

# n引数に大きな値を指定してFFTを実行
Y_padded = np.abs(fft(y, n=N_padded))
freqs_padded = fftfreq(N_padded, d=dt)
# --- 6. FFT ---


# 正の周波数のみ抽出
mask2 = freqs_padded > 0
freqs = freqs_padded[mask2]
Y = Y_padded[mask2]

# --- 7. ピーク検出 ---
peak_idx = np.argmax(Y)
peak_freq = freqs[peak_idx]
print(f"主要周波数: {peak_freq:.2f} Hz")

# --- 8. プロット ---
plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.plot(t_seg, y)
plt.title("Waveform (0.15s - 0.30s)")
plt.xlabel("Time [s]")
plt.ylabel("Min Area")
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(freqs, Y)
plt.xlim(0, 500)
plt.title(f"Spectrum (Peak: {peak_freq:.2f} Hz)")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Magnitude")
plt.grid()

plt.tight_layout()
plt.show()