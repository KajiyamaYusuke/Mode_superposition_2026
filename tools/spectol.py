import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scienceplots

# 論文用のスタイル設定（ユーザー環境と同一）
plt.style.use(['science', 'ieee', 'no-latex'])

# =========================
# 設定
# =========================
filename = "../output/airflow_vt.dat"
sim_dt = 1.0e-5
output_interval = 5
dt = sim_dt * output_interval
fs = 1.0 / dt  # サンプリング周波数 (20,000 Hz)

# =========================
# データ読み込みと解析
# =========================
data = np.loadtxt(filename, comments='#')
pressure = data[:, 1]

# ★スペクトログラムは「変化の瞬間」を見るのが最大の目的なので、
# 立ち上がり（0秒）から全て描画することをおすすめします。
t_start = 0.0 
t_end   = 0.4
start_idx = int(t_start / dt)
end_idx   = int(t_end / dt)

if end_idx > len(pressure):
    end_idx = len(pressure)

valid_pressure = pressure[start_idx:end_idx]

# DC成分の除去
valid_pressure = valid_pressure - np.mean(valid_pressure)

# =========================
# スペクトログラム計算 (scipy.signal.spectrogram)
# =========================
# nperseg (窓幅): 時間分解能と周波数分解能のバランスを決める最重要パラメータ
# fs=20000Hzの場合、1024点で約0.05秒(約19.5Hz刻み)。声帯振動の解析に最適です。
nperseg = 1024 
noverlap = nperseg * 3 // 4  # 75%オーバーラップさせて画像を滑らかにする

f, t_spec, Sxx = signal.spectrogram(valid_pressure, fs=fs, window='hann',
                                    nperseg=nperseg, noverlap=noverlap)

# 時間軸を切り出し開始時間に合わせる
t_spec = t_spec + t_start

# ★振幅の0dB正規化（FFTコードと同じロジック）
# Sxxはパワースペクトル（振幅の2乗）なので、平方根をとってLinear振幅に戻す
Sxx_linear = np.sqrt(Sxx)
max_Sxx = np.max(Sxx_linear)
Sxx_db = 20 * np.log10(Sxx_linear / max_Sxx + 1e-12)

# =========================
# プロット
# =========================
fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

# pcolormeshで描画 (カラーマップは論文で映える 'magma' や 'viridis' がおすすめ)
# vmin=-80 で、-80dB以下の微小なノイズをカットしてコントラストを上げます
mesh = ax.pcolormesh(t_spec, f / 1000, Sxx_db, shading='gouraud', cmap='magma', vmin=-80, vmax=0)

ax.set_xlabel("Time [s]", fontsize=20)
ax.set_ylabel("Frequency [kHz]", fontsize=20)

# 注目したい周波数帯域（例: 0〜6kHz）に絞る
ax.set_ylim(0, 6)
ax.tick_params(direction='in', labelsize=14, top=True, right=True)

# カラーバーの追加
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label('Sound Pressure Level [dB]', fontsize=16)
cbar.ax.tick_params(labelsize=14)

plt.tight_layout()
# plt.savefig("result_spectrogram.png", dpi=300)
plt.show()