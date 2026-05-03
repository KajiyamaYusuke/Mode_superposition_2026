import numpy as np
import matplotlib.pyplot as plt
# 1. datファイルを読み込む（空白区切りやタブ区切りを想定）
harea = np.loadtxt("../output/area.dat")

labels = harea[:, 0]        # ← 1列目（横軸のラベル）
values = harea[:, 1:]  

# 2. 各行の最小値を求める
row_min = np.min(values, axis=1)  # 各行の最小値を1次元配列で返す

x = labels * 1e-5

plt.figure(figsize=(8,4))
plt.plot(x, row_min, linestyle='-', color='g')  # ← マーカーなし線だけ
plt.xlabel('time')
plt.ylabel('area')
plt.title('glottal area')
# plt.ylim(-0.3,8.2)
plt.grid(True)
plt.tight_layout()
plt.show()
