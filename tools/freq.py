import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt
from scipy.fft import fft, fftfreq

# 1. データ読み込み
harea = np.loadtxt("../output/area.dat")
labels = harea[:, 0]
values = harea[:, 1:]

x = labels * 1e-5  # 秒
row_min = np.min(values, axis=1)

t_start = 0.3
t_end   = 0.5
idx = np.where((x >= t_start) & (x <= t_end))[0]

x_seg = x[idx]
row_seg = row_min[idx]

# 2. ノイズ除去（低域通過フィルタ例）
fs = 1 / (x_seg[1] - x_seg[0])  # サンプリング周波数
nyq = fs / 2
cutoff = 1000  # Hz, 音声なら 1 kHz 程度まで残す
b, a = butter(4, cutoff / nyq, btype='low')
row_seg_smooth = filtfilt(b, a, row_seg)

# 3. FFT
N = len(row_seg_smooth)
yf = fft(row_seg_smooth)
xf = fftfreq(N, d=1/fs)

# 正の周波数だけ
pos = xf > 0
xf = xf[pos]
yf = np.abs(yf[pos])

peak_freq = xf[np.argmax(yf)]
print(f"FFTでの主要周波数: {peak_freq:.2f} Hz")

# 4. ピーク検出（スムーズ化波形で）
peaks, _ = find_peaks(row_seg_smooth, prominence=0.001, distance=int(fs/peak_freq/2))
peak_times = x_seg[peaks]
periods = np.diff(peak_times)
avg_period = np.mean(periods)
frequency = 1 / avg_period

print(f"ピーク間平均周波数: {frequency:.2f} Hz")

# 5. 可視化
plt.figure(figsize=(12,5))
plt.plot(x_seg, row_seg, alpha=0.5, label="raw")
plt.plot(x_seg, row_seg_smooth, label="smoothed")
plt.plot(peak_times, row_seg_smooth[peaks], "ro", label="peaks")
plt.xlabel("time [s]")
plt.ylabel("area")
plt.title("Glottal area with peaks")
plt.grid(True)
plt.legend()
plt.show()

# FFTスペクトル確認
plt.figure(figsize=(10,4))
plt.plot(xf, yf)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.title("FFT (0.15-0.20 s)")
plt.grid(True)
plt.show()
