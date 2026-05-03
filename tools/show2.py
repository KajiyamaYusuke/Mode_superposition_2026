import numpy as np
import matplotlib.pyplot as plt

# --- ファイル名を指定 ---
filename = "../output/displace.dat"

# --- データ読み込み ---
# 空白区切りなので delimiter=' '（自動判別でもOK）
data = np.loadtxt(filename)

# 1列目: 時間, 2列目: 変位
time = data[:, 0]
dispY = data[:, 1]
dispX = data[:, 2]

# --- グラフ描画 ---
plt.figure(figsize=(6,4))
plt.plot(time, dispY, linewidth=2, color='g', label='left')
plt.plot(time, dispX, linewidth=2, color='orange', label='right')
plt.xlabel("Time [s]")
plt.ylabel("Displacement [mm]")
plt.title("Node displacement over time")
plt.legend(loc='best', fontsize=10, frameon=True)
plt.grid(True)
plt.tight_layout()
plt.show()
