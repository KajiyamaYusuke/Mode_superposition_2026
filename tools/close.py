import numpy as np
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt

# --- 1. データの読み込み ---
harea = np.loadtxt("../output/airflow_vt.dat")
areas = harea[:, 1:]

# --- 2. 最小断面積の時系列 ---
row_min = np.min(areas, axis=1)

# --- 3. 正しい時間軸の作成 ---
dt = 0.00005
N_total = len(row_min)
t = np.arange(N_total) * dt 

# --- 4. 分析区間の抽出 ---
mask = (t >= 0.1) & (t <= 0.5)
y = row_min[mask]      # ★ここではまだ生の断面積（正の値）のままにする
t_seg = t[mask]

print(f"解析データ点数: {len(y)} 点") 

# --- (削除) 5. DC成分除去 ---
# y = y - np.mean(y)  <-- ★★★ この行を削除、またはコメントアウト ★★★

# --- 6. FFT (周波数特定用) ---
target_resolution = 0.1
N_padded = int(1 / (dt * target_resolution)) 

# ★FFTにかけるときだけ平均を引いて、DC成分(0Hz)のピークを抑制する
y_for_fft = y - np.mean(y) 
Y_padded = np.abs(fft(y_for_fft, n=N_padded))
freqs_padded = fftfreq(N_padded, d=dt)

mask2 = freqs_padded > 0
freqs = freqs_padded[mask2]
Y = Y_padded[mask2]

# --- 7. ピーク検出 ---
peak_idx = np.argmax(Y)
peak_freq = freqs[peak_idx]

# --- 9. 閉鎖判定 ---
# 単位に注意してください。もし断面積が m^2 なら 1e-6 は 1mm^2 です。
# 完全に閉じるモデルなら 1e-10 程度でも良い場合があります。
area_eps = 1e-6 

# ★ y は平均を引いていないので、正しく「0に近いかどうか」で判定される
is_closed = y <= area_eps

T_closed_total = np.sum(is_closed) * dt
# print(f"閉鎖時間（総和）: {T_closed_total:.6f} s") # 必要なら表示

T_period = 1.0 / peak_freq
print(f"周期: {T_period:.6f} s")
print(f"基本周波数: {peak_freq:.2f} Hz")

# 閉鎖区間の開始・終了インデックスを検出
diff = np.diff(is_closed.astype(int))
start_idx = np.where(diff == 1)[0] + 1
end_idx   = np.where(diff == -1)[0] + 1

# 端点補正
if is_closed[0]:
    start_idx = np.insert(start_idx, 0, 0)
if is_closed[-1]:
    end_idx = np.append(end_idx, len(is_closed))

CQ_list = []

# エラー回避: 閉鎖区間が検出されなかった場合の処理
if len(start_idx) == 0:
    print("閉鎖区間が検出されませんでした。area_eps の値を見直してください。")
else:
    for s, e in zip(start_idx, end_idx):
        closed_time = (e - s) * dt
        CQ = closed_time / T_period
        
        # 異常値（周期より長い閉鎖など）を除外する場合
        if CQ <= 1.0: 
            CQ_list.append(CQ)

    CQ_list = np.array(CQ_list)
    print(f"平均 Closed Quotient: {np.mean(CQ_list):.3f}")
    # print(f"Closed Quotient 分散: {np.std(CQ_list):.3f}")

# --- プロット確認 ---
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(t_seg, y, label="Min Area (Raw)")
# 閾値を赤い点線で表示
plt.axhline(y=area_eps, color='orange', linestyle='--', label=f"Threshold ({area_eps})")

plt.fill_between(
    t_seg, y.min(), y.max(),
    where=is_closed,
    color="red", alpha=0.3,
    label="Closed Detected"
)
plt.legend()
plt.ylabel("Area")
plt.title(f"Analysis Range ({0.1}-{0.14}s)")
plt.show()