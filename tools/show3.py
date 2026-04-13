import numpy as np
import matplotlib.pyplot as plt
import scienceplots

# --- ファイル名を指定 ---
file1 = "../output/displace_c5.dat"   # 1点目
file2 = "../output/displace2_c5.dat"  # 2点目

# file1 = "../output/displace.dat"   # 1点目
# file2 = "../output/displace2.dat"  # 2点目

# --- データ読み込み ---
try:
    data1 = np.loadtxt(file1)
    data2 = np.loadtxt(file2)
except OSError:
    print("エラー: ファイルが見つかりません。")
    exit()

# --- データ整理 ---
# 時間軸 (両方のファイルで同じと仮定)
raw_time = data1[:, 0]

start_time = 0.15
mask = raw_time >= start_time

# マスクを適用してデータを切り出し
time = raw_time[mask]

# Point 1 (x=10.0)
p1_lat = -data1[mask, 1] # Lateral (Y) - 符号反転
p1_ver = data1[mask, 2]  # Vertical (X)

# Point 2 (x=7.2)
p2_lat = -data2[mask, 1] # Lateral (Y) - 符号反転
p2_ver = data2[mask, 2]  # Vertical (X)

# --- プロット設定 (共通) ---
plt.style.use(['science','muted', 'ieee', 'no-latex'])


# ==========================================
# 2. 垂直変位 (Vertical) のグラフ
# ==========================================
plt.figure(figsize=(4, 4))
plt.plot(p1_lat, p1_ver, label=r'upper point',color='cornflowerblue', linestyle='-')
plt.plot(p2_lat, p2_ver, label=r'lower point',  color='tab:orange', linestyle='-')

plt.xlabel(r'lateral displacement [mm]', fontsize = 12)
plt.ylabel(r'Vertical displacement [mm]',fontsize = 12)
plt.legend(loc='upper left', frameon=True, fontsize=8)
# plt.axis('equal')

plt.xlim(0.0, 1)

plt.ylim(1.6, 2.4)

plt.tight_layout()
plt.savefig('../output/displacement.png', dpi=300, bbox_inches='tight')
plt.close()

print("Saved: displacement.png")