import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks
import math
import scienceplots

# 論文用のスタイル設定
plt.style.use(['science', 'ieee', 'no-latex'])

# =========================
# ユーザー設定環境
# =========================
EXECUTABLE = "/home/kajiyama/code/Mode_superposition_2026/build/simulation"  
PARAM_FILE = "/home/kajiyama/code/Mode_superposition_2026/input/param.txt"     
OUTPUT_DATA = "/home/kajiyama/code/Mode_superposition_2026/output/airflow_vt.dat"  
DISP_DATA = "/home/kajiyama/code/Mode_superposition_2026/output/displace.dat"  # ★追加：変位データ

# スイープする圧力のリスト (Pa)
pressure_list = [700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

# 解析設定
sim_dt = 1.0e-5
output_interval = 5  
dt = sim_dt * output_interval
fs = 1.0 / dt         
t_start = 0.0 
t_end   = 0.15        
flow_plot_start = 0.05
flow_plot_end = 0.1

# =========================
# 関数: paramファイルの圧力書き換え
# =========================
def update_param_file(filepath, new_pressure):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    for i in range(len(lines)):
        if "# subglottal pressure Ps (Pa)" in lines[i]:
            lines[i+1] = f"{new_pressure}\n"
            break
            
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(lines)
    print(f"[設定更新] 声門下圧を {new_pressure} Pa に変更しました。")

def save_flow_waveform(pressure_val):
    data_flow = np.loadtxt(OUTPUT_DATA, comments='#')
    time = data_flow[:, 0] * sim_dt
    flow = data_flow[:, 1]

    mask = (time >= flow_plot_start) & (time <= flow_plot_end)
    if not np.any(mask):
        print(
            f"[スキップ] 流量波形 {flow_plot_start:.1f}-{flow_plot_end:.1f} s のデータがありません。"
            f" 現在の最終時刻は {time[-1]:.5f} s です。"
        )
        return

    fig, ax = plt.subplots(figsize=(7, 3), dpi=150)
    ax.plot(time[mask], flow[mask], color='black', linewidth=1.0)
    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Airflow [l/s]", fontsize=14)
    ax.set_title(f"Airflow waveform (Ps = {pressure_val} Pa)", fontsize=15)
    ax.tick_params(direction='in')
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    save_filename = f"../output/airflow_waveform_{flow_plot_start:.1f}_{flow_plot_end:.1f}s_Ps_{pressure_val}Pa.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 流量波形を保存しました: {save_filename}")

# =========================
# 関数: シミュレーション結果の可視化とピーク解析
# =========================
def analyze_and_plot(pressure_val):
    # 2つのデータを読み込む
    data_flow = np.loadtxt(OUTPUT_DATA, comments='#')
    data_disp = np.loadtxt(DISP_DATA, comments='#')
    
    # 万が一、書き込みのタイミングで行数がずれた場合の安全策
    min_len = min(len(data_flow), len(data_disp))
    
    # 時間軸と各データの抽出
    time = np.arange(min_len) * dt
    flow_data = data_flow[:min_len, 1]
    
    # displace.datから変位を取得 (1列目=L, 2列目=R)
    x1l_data = data_disp[:min_len, 1] 
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

    # =========================
    # 極大値（ピーク）の検出と比率計算
    # =========================
    zoom_start = max(0, len(valid_time) - int(0.05 / dt)) # 最後の0.05秒間で解析
    steady_time = valid_time[zoom_start:]
    steady_flow = valid_flow[zoom_start:]
    steady_x1l  = valid_x1l[zoom_start:]
    steady_x1r  = valid_x1r[zoom_start:]

    # ノイズを拾わないよう、波形の振幅の10%をピーク検出の閾値（prominence）に設定
    prom_l = max((np.max(steady_x1l) - np.min(steady_x1l)) * 0.1, 1e-8)
    prom_r = max((np.max(steady_x1r) - np.min(steady_x1r)) * 0.1, 1e-8)

    peaks_l, _ = find_peaks(steady_x1l, prominence=prom_l)
    peaks_r, _ = find_peaks(-steady_x1r, prominence=prom_r)

    num_l = len(peaks_l)
    num_r = len(peaks_r)
    divisor = math.gcd(num_l, num_r) if (num_l > 0 and num_r > 0) else 1
    ratio_l = num_l // divisor
    ratio_r = num_r // divisor
    ratio_str = f"{ratio_l}:{ratio_r}"
    
    print(f"  -> 極大値検出: 左={num_l}回, 右={num_r}回 (比率 {ratio_str})")

    # =========================
    # 描画 (3つのサブプロット)
    # =========================
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), dpi=150, gridspec_kw={'height_ratios': [1, 1, 1.5]})
    
    # 上段: 変位の波形とピーク位置
    ax1.plot(steady_time, steady_x1l, label='Left VF', color='blue', linewidth=1.0, alpha=0.7)
    ax1.plot(steady_time, steady_x1r, label='Right VF', color='red', linewidth=1.0, alpha=0.7)
    # ピーク位置にバツ印を打つ
    ax1.plot(steady_time[peaks_l], steady_x1l[peaks_l], "x", color='blue', markersize=6)
    ax1.plot(steady_time[peaks_r], steady_x1r[peaks_r], "x", color='red', markersize=6)
    ax1.set_title(f"Ratio of Maxima = {ratio_str} (Ps = {pressure_val} Pa)", fontsize=16)
    ax1.set_ylabel("Displacement", fontsize=14)
    ax1.tick_params(direction='in')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5)

    # 中段: 気流波形
    ax2.plot(steady_time, steady_flow, color='black', linewidth=1.0)
    ax2.set_ylabel("Flow/Pressure", fontsize=14)
    ax2.tick_params(direction='in')
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # 下段: スペクトログラム
    mesh = ax3.pcolormesh(t_spec, f / 1000, Sxx_db, shading='gouraud', cmap='magma', vmin=-80, vmax=0)
    ax3.set_xlabel("Time [s]", fontsize=16)
    ax3.set_ylabel("Frequency [kHz]", fontsize=16)
    ax3.set_ylim(0, 4) 
    ax3.tick_params(direction='in')

    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label('SPL [dB]', fontsize=14)

    plt.tight_layout()
    
    # 画像の保存
    save_filename = f"../output/analysis_Ps_{pressure_val}Pa.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 画像を保存しました: {save_filename}\n")


# =========================
# メインのスイープ実行ループ
# =========================
if __name__ == "__main__":
    for p in pressure_list:
        print(f"========== Start Simulation for Ps = {p} Pa ==========")
        update_param_file(PARAM_FILE, p)
        
        print("[実行中] C++ソルバーを計算しています...")
        # cwdを付与してルートで実行させる
        result = subprocess.run([EXECUTABLE], capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"[エラー] シミュレーションが異常終了しました (Ps={p}Pa)")
            print(result.stderr)
            continue
            
        analyze_and_plot(p)
        save_flow_waveform(p)
        
    print("========== 全てのスイープ解析が完了しました！ ==========")
