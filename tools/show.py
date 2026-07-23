import numpy as np
import matplotlib.pyplot as plt
import scienceplots

# 論文用のスタイル設定
plt.style.use(['science', 'ieee', 'no-latex'])

# 1. datファイルを読み込む（空白区切りやタブ区切りを想定）
harea = np.loadtxt("../runs/param_20260722_190210/area.dat")

labels = harea[:, 0]        # ← 1列目（横軸のラベル）
values = harea[:, 1:]  

# 2. 各行の最小値を求める
row_min = np.min(values, axis=1)  # 各行の最小値を1次元配列で返す

x = labels * 1e-5
# mask = x >= 0.1
# x = x[mask]
# row_min = row_min[mask]

plt.figure(figsize=(12,4), dpi=100)
plt.plot(x, row_min, linestyle='-', color='#0072B2')  # ← マーカーなし線だけ
plt.xlabel('Time [s]', fontsize=18)
plt.ylabel('area [mm^2]', fontsize=18)
#plt.title('glottal area')
# plt.ylim(-0.3,8.2)
plt.grid(True)
plt.tight_layout()
plt.show()
