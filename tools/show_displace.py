import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scienceplots

# 論文用のスタイル設定
plt.style.use(['science', 'ieee', 'no-latex'])

# =========================
# ユーザー設定環境
# =========================
EXECUTABLE = "../build/simulation"  
PARAM_FILE = "../input/param.txt"     
OUTPUT_DATA = "../output/airflow_vt.dat"  
DISP_DATA = "../output/displace.dat"  # ★追加：変位データ


# 解析設定
sim_dt = 1.0e-5
output_interval = 5  
dt = sim_dt * output_interval
fs = 1.0 / dt
t_start = 0.05
t_end   = 0.5
flow_plot_start = 0.05
flow_plot_end = 0.1
pressure_val = 500

def save_displacement_waveform(pressure_val, steady_time, steady_x1l, steady_x1r):
    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)

    ax.plot(steady_time, steady_x1l, label='Left VF', color='#0072B2', linewidth=1.2)
    ax.plot(steady_time, steady_x1r, label='Right VF', color='#D55E00', linewidth=1.2, linestyle='--')
    #ax.set_title(f"Vocal Fold Displacement (Ps = {pressure_val} Pa)", fontsize=16)
    ax.set_xlabel("Time [s]", fontsize=26)
    ax.set_ylabel("Displacement [mm]", fontsize=26)
    ax.tick_params(direction='in', labelsize=20)
    #ax.set_ylim(-0.85, 3.3)
    ax.legend(loc='upper right', fontsize=18)
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    save_filename = f"../output/displacement_Ps_{pressure_val}Pa.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 声帯変位波形を保存しました: {save_filename}")


data_flow = np.loadtxt(OUTPUT_DATA, comments='#')
data_disp = np.loadtxt(DISP_DATA, comments='#')

# 万が一、書き込みのタイミングで行数がずれた場合の安全策
min_len = min(len(data_flow), len(data_disp))

# 時間軸と各データの抽出
time = np.arange(min_len) * dt
flow_data = data_flow[:min_len, 1]

# displace.datから変位を取得 (1列目=L, 2列目=R)
# 左声帯は符号が右声帯と逆なので、重ねて比較しやすいように反転する
x1l_data = -data_disp[:min_len, 1] 
x1r_data = data_disp[:min_len, 2] 

# 指定区間の切り出し
start_idx = int(t_start / dt)
end_idx   = min(int(t_end / dt), min_len)

valid_time = time[start_idx:end_idx]
valid_flow = flow_data[start_idx:end_idx]
valid_x1l  = x1l_data[start_idx:end_idx]
valid_x1r  = x1r_data[start_idx:end_idx]

valid_flow_ac = valid_flow - np.mean(valid_flow)

# スペクトログラム計算
nperseg = 1024 
noverlap = nperseg * 3 // 4 
f, t_spec, Sxx = signal.spectrogram(valid_flow_ac, fs=fs, window='hann',
                                    nperseg=nperseg, noverlap=noverlap)
t_spec = t_spec + t_start

Sxx_linear = np.sqrt(Sxx)
max_Sxx = np.max(Sxx_linear) if np.max(Sxx_linear) > 0 else 1.0
Sxx_db = 20 * np.log10(Sxx_linear / max_Sxx + 1e-12)

zoom_start = max(0, len(valid_time) - int(0.15 / dt)) # 最後の0.05秒間で解析
steady_time = valid_time[zoom_start:]
steady_flow = valid_flow[zoom_start:]
steady_x1l  = valid_x1l[zoom_start:]
steady_x1r  = valid_x1r[zoom_start:]

save_displacement_waveform(pressure_val, steady_time, steady_x1l, steady_x1r)