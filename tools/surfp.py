

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# ==== データ読み込み ====
# ファイル名を適宜変更してください
filename = "../output/surfp_output.csv"

# コメント行 (#) をスキップして読み込む
df = pd.read_csv(filename, comment='#', header=None)
df.columns = ["i", "j", "node_id", "x", "y", "z"]

# ==== 可視化 ====
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# i の値で色を変える
sc = ax.scatter(df['x'], df['y'], df['z'],
                c=df['j'], cmap='viridis', s=40)

# カラーバーを追加
cb = plt.colorbar(sc, ax=ax, label='i index')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Node Distribution (colored by i, marker by j)')
plt.tight_layout()
plt.show()

