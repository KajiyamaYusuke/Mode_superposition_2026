import numpy as np
import matplotlib.pyplot as plt

import scienceplots

plt.style.use(['science','ieee', 'no-latex'])

# =========================
# 設定
# =========================
filename = "../output/pressure_vt.dat"
sim_dt = 1.0e-5
output_interval = 5
dt = sim_dt * output_interval

# =========================
# データ読み込みと解析
# =========================
data = np.loadtxt(filename, comments='#')
pressure = data[:, 1]

# 時間切り出し (0.15s - 0.4s)
t_start = 0.1
t_end   = 0.4
start_idx = int(t_start / dt)
end_idx   = int(t_end / dt)

if end_idx > len(pressure):
    end_idx = len(pressure)

valid_pressure = pressure[start_idx:end_idx]

# DC成分除去 & 窓関数
valid_pressure = valid_pressure - np.mean(valid_pressure)
window = np.hanning(len(valid_pressure))
valid_pressure_windowed = valid_pressure * window

# FFT計算
N = len(valid_pressure)
freq = np.fft.rfftfreq(N, d=dt)
fft_val = np.fft.rfft(valid_pressure_windowed)

# 振幅スペクトル (Linear)
amplitude = np.abs(fft_val) / N * 2

# ★変更点1: 最大値を0dBにする正規化
# 最大振幅を見つける
max_linear_amp = np.max(amplitude)
# 最大値で割ってから対数をとる (+1e-12は0除算防止)
db_amplitude = 20 * np.log10(amplitude / max_linear_amp + 1e-12)

# --- F0検出ロジック (変更なし) ---
mask = freq > 20
masked_freq = freq[mask]
masked_amp = amplitude[mask] # ピーク検出はLinear振幅で行うのが安全

f0 = 0
if len(masked_amp) > 2:
    is_peak = (masked_amp[1:-1] > masked_amp[:-2]) & (masked_amp[1:-1] > masked_amp[2:])
    is_peak = np.r_[False, is_peak, False]
    
    peak_indices = np.where(is_peak)[0]
    peak_freqs = masked_freq[peak_indices]
    peak_amps = masked_amp[peak_indices]
    
    # 最大値の10% (-20dB) を閾値とする
    threshold = 0.1 * max_linear_amp 
    significant_peaks_idx = np.where(peak_amps > threshold)[0]
    
    if len(significant_peaks_idx) > 0:
        first_peak_idx = significant_peaks_idx[0]
        f0 = peak_freqs[first_peak_idx]

print(f"Detected Fundamental Frequency (F0): {f0:.2f} Hz")

# =========================
# プロット
# =========================
fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

# ★変更点2: 色を変更 (color='...')
# 色の例: 'steelblue', 'firebrick', 'darkgreen', 'navy', 'black', 'orange'
ax.plot(freq / 1000, db_amplitude, color='cornflowerblue',alpha = 0.8, label='Spectrum')

#ax.set_title(f"Normalized Frequency Domain (F0 = {f0:.2f} Hz)")
ax.set_xlabel("Frequency [kHz]", fontsize=20)
ax.set_ylabel("Sound Pressure Level [dB]", fontsize=20) # ラベルも変更

# 0dBが最大なので、上限を少し余裕を持たせて設定
ax.set_ylim(-140, 5) 
ax.set_xlim(0, 6)
ax.tick_params(direction='in', labelsize=14, top=True, right=True)
ax.grid(which='both', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("result_fft.png", dpi=300)
plt.show()