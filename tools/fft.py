import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt

# --- 1. load data ---
harea = np.loadtxt("../output/area.dat")
steps = harea[:, 0]
areas = harea[:, 1:]

# --- 2. 最小断面積の時系列を作る ---
row_min = np.min(areas, axis=1)

# --- 3. 時間軸 ---
dt = 0.0002                # 1ステップ=0.0002秒
fs = 1/dt                  # サンプリング周波数 (5000 Hz)
t = steps * dt             # 時間配列

# --- 4. 必要なら区間抽出 ---
mask = (t >= 0.15) & (t <= 0.30)
y = row_min[mask]
t_seg = t[mask]

# --- 5. DC成分除去 ---
y = y - np.mean(y)

# --- 6. FFT ---
N = len(y)
Y = np.abs(fft(y))
freqs = fftfreq(N, d=dt)

# 正の周波数のみ
mask2 = freqs > 0
freqs = freqs[mask2]
Y = Y[mask2]

# --- 7. 最大ピーク ---
peak_freq = freqs[np.argmax(Y)]
print("主要周波数:", peak_freq, "Hz")

# --- 8. debug plot ---
plt.plot(freqs, Y)
plt.xlim(0, 1000)  # 声帯のF0はだいたい0〜500Hz
plt.grid()
plt.show()
