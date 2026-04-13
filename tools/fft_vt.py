import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# =========================
# 設定
# =========================
filename = "../output/pressure_vt.dat"  # ファイルパス
sim_dt = 1.0e-5
output_interval = 5
dt = sim_dt * output_interval

# =========================
# データ処理
# =========================
try:
    data = np.loadtxt(filename, comments='#')
    pressure = data[:, 1]
    
    t_start = 0.15  # 開始時間 (秒)
    t_end   = 0.4   # 終了時間 (秒)

    # 時間をインデックス（配列の何番目か）に変換
    start_idx = int(t_start / dt)
    end_idx   = int(t_end / dt)

    # 配列の範囲外エラーを防ぐための安全策
    if end_idx > len(pressure):
        end_idx = len(pressure)

    # スライス機能で指定範囲を抽出
    valid_pressure = pressure[start_idx:end_idx]
    
    # DC除去 & 窓関数
    valid_pressure = valid_pressure - np.mean(valid_pressure)
    window = np.hanning(len(valid_pressure))
    
    # FFT計算
    N = len(valid_pressure)
    freq = np.fft.rfftfreq(N, d=dt)
    fft_val = np.fft.rfft(valid_pressure * window)
    
    # dB変換 (Sound Pressure Level)
    amplitude = np.abs(fft_val) / N * 2
    # 0除算防止
    db_amplitude = 20 * np.log10(amplitude + 1e-12 / 20e-6)
    
    # =================================================
    # ★追加：H1, H2, H2k の特定と指標計算
    # =================================================
    
    # ピーク検出 (50Hz以上)
    search_mask = freq > 50.0
    masked_db = db_amplitude[search_mask]
    masked_freq = freq[search_mask]
    
    # 最大値から20dB以内にある山をすべて候補とする
    max_db = np.max(masked_db)
    peaks, _ = find_peaks(masked_db, height=max_db - 20.0, distance=5)
    
    H1_freq, H1_amp = 0, -100
    H2_freq, H2_amp = 0, -100
    H2k_freq, H2k_amp = 0, -100 # 追加: 2kHz付近のピーク
    
    H1_H2 = 0
    H1_H2k = 0 # 追加
    
    if len(peaks) > 0:
        # --- H1 (第1高調波 = 基本周波数 F0) ---
        h1_idx = peaks[0]
        H1_freq = masked_freq[h1_idx]
        H1_amp = masked_db[h1_idx]
        
        # --- H2 (第2高調波) ---
        target_f2 = H1_freq * 2.0
        # ±20Hzの範囲で探す
        h2_mask = (masked_freq > target_f2 - 20.0) & (masked_freq < target_f2 + 20.0)
        
        if np.any(h2_mask):
            band_freqs = masked_freq[h2_mask]
            band_dbs = masked_db[h2_mask]
            local_max_idx = np.argmax(band_dbs)
            H2_freq = band_freqs[local_max_idx]
            H2_amp = band_dbs[local_max_idx]
            H1_H2 = H1_amp - H2_amp

        # --- H2k (2000Hzに最も近い高調波) ---
        # 2000Hzを中心に ±(F0/2) の範囲で最大のピークを探すのが安全
        target_f2k = 2000.0
        search_width = H1_freq * 0.5 # F0の半分くらいの幅で探す
        
        h2k_mask = (masked_freq > target_f2k - search_width) & (masked_freq < target_f2k + search_width)
        
        # もし範囲内に強いピークがなければ、純粋に2000Hzに一番近いピークを探す
        if np.any(h2k_mask):
            band_freqs = masked_freq[h2k_mask]
            band_dbs = masked_db[h2k_mask]
            local_max_idx = np.argmax(band_dbs)
            H2k_freq = band_freqs[local_max_idx]
            H2k_amp = band_dbs[local_max_idx]
        else:
            # 2000Hz周辺にめぼしいピークがない場合、最も近い周波数の値を拾う
            idx_closest = np.argmin(np.abs(masked_freq - 2000.0))
            H2k_freq = masked_freq[idx_closest]
            H2k_amp = masked_db[idx_closest]

        H1_H2k = H1_amp - H2k_amp

        print(f"H1 (F0): {H1_freq:.2f} Hz, {H1_amp:.2f} dB")
        if H2_freq > 0:
            print(f"H2 (2f): {H2_freq:.2f} Hz, {H2_amp:.2f} dB")
        print(f"H2k (~2kHz): {H2k_freq:.2f} Hz, {H2k_amp:.2f} dB")
        print("-" * 30)
        print(f"H1-H2  : {H1_H2:.2f} dB")
        print(f"H1-H2k : {H1_H2k:.2f} dB") # これが求める指標

    else:
        print("No peaks found.")

    # =========================
    # プロット
    # =========================
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # 1. 時間波形
    axes[0].plot(data[:, 0], data[:, 1])
    axes[0].set_title("Time Domain")
    axes[0].set_xlabel("Step")
    axes[0].set_ylabel("Pressure [Pa]")
    axes[0].grid()
    
    # 2. FFT結果
    axes[1].plot(freq, db_amplitude, label='Spectrum')
    axes[1].set_title(f"Freq Domain | H1-H2={H1_H2:.1f}dB, H1-H2k={H1_H2k:.1f}dB")
    axes[1].set_xlabel("Frequency [Hz]")
    axes[1].set_ylabel("SPL [dB]")
    axes[1].set_xlim(0, 3000) # 2kHzが見えるように
    axes[1].grid()
    
    # マーカー表示
    if H1_freq > 0:
        axes[1].plot(H1_freq, H1_amp, 'ro', label='H1')
    if H2_freq > 0:
        axes[1].plot(H2_freq, H2_amp, 'bo', label='H2')
    if H2k_freq > 0:
        axes[1].plot(H2k_freq, H2k_amp, 'go', label='H2k')
        axes[1].annotate(f"H2k", (H2k_freq, H2k_amp), textcoords="offset points", xytext=(0,10), ha='center')

    axes[1].legend()
    plt.tight_layout()
    plt.show()

except Exception as e:
    print(f"Error: {e}")