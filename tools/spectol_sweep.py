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
EXECUTABLE = "/home/kajiyama/code/Mode_superposition_2026/build/simulation"  
PARAM_FILE = "/home/kajiyama/code/Mode_superposition_2026/input/param.txt"     
OUTPUT_DATA = "/home/kajiyama/code/Mode_superposition_2026/output/airflow_vt.dat"  
DISP_DATA = "/home/kajiyama/code/Mode_superposition_2026/output/displace.dat"  # ★追加：変位データ

# スイープする圧力のリスト (Pa)
pressure_list = [600, 700, 800, 900, 1000, 1100, 1200, 1300]

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
# 関数: シミュレーション結果の可視化
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

    zoom_start = max(0, len(valid_time) - int(0.05 / dt)) # 最後の0.05秒間で解析
    steady_time = valid_time[zoom_start:]
    steady_flow = valid_flow[zoom_start:]
    steady_x1l  = valid_x1l[zoom_start:]
    steady_x1r  = valid_x1r[zoom_start:]

    save_displacement_waveform(pressure_val, steady_time, steady_x1l, steady_x1r)

    # =========================
    # 描画 (2つのサブプロット)
    # =========================
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), dpi=150, gridspec_kw={'height_ratios': [1, 1.5]})

    # 上段: 気流波形
    ax1.plot(steady_time, steady_flow, color='black', linewidth=1.0)
    ax1.set_title(f"Airflow and Spectrogram (Ps = {pressure_val} Pa)", fontsize=16)
    ax1.set_ylabel("Flow/Pressure", fontsize=14)
    ax1.tick_params(direction='in')
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # 下段: スペクトログラム
    mesh = ax2.pcolormesh(t_spec, f / 1000, Sxx_db, shading='gouraud', cmap='magma', vmin=-80, vmax=0)
    ax2.set_xlabel("Time [s]", fontsize=16)
    ax2.set_ylabel("Frequency [kHz]", fontsize=16)
    ax2.set_ylim(0, 4) 
    ax2.tick_params(direction='in')

    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label('SPL [dB]', fontsize=14)

    plt.tight_layout()
    
    # 画像の保存
    save_filename = f"../output/analysis_Ps_{pressure_val}Pa.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 画像を保存しました: {save_filename}\n")


def save_displacement_waveform(pressure_val, steady_time, steady_x1l, steady_x1r):
    fig, ax = plt.subplots(figsize=(8, 3), dpi=150)

    ax.plot(steady_time, steady_x1l, label='Left VF', color='#0072B2', linewidth=1.2)
    ax.plot(steady_time, steady_x1r, label='Right VF', color='#D55E00', linewidth=1.2, linestyle='--')
    ax.set_title(f"Vocal Fold Displacement (Ps = {pressure_val} Pa)", fontsize=16)
    ax.set_xlabel("Time [s]", fontsize=18)
    ax.set_ylabel("Displacement", fontsize=18)
    ax.set_ylim(-0.85, 3.3)
    ax.tick_params(direction='in')
    ax.legend(loc='upper right', fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    save_filename = f"../output/displacement_waveform_Ps_{pressure_val}Pa.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 声帯変位波形を保存しました: {save_filename}")


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
