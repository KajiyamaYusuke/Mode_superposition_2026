
import numpy as np

# row_min: 各ステップの最小面積
# t: 各サンプルの時間（秒）
# 例: t = labels * 1e-5

harea = np.loadtxt("../output/area.dat")
labels = harea[:, 0]
values = harea[:, 1:]

row_min = np.min(values, axis=1)

t = labels * 1e-5

# 0.05秒以降のデータを抽出
indices = np.where(t > 0.05)[0]
T_sub = row_min[indices]
N_sub = len(T_sub)

# 幅5の移動平均で局所平均との差の絶対値を計算
diff = []
for i in range(N_sub):
    start = max(0, i-2)
    end   = min(N_sub, i+3)  # i+2まで含む
    local_avg = np.mean(T_sub[start:end])
    diff.append(np.abs(T_sub[i] - local_avg))

diff = np.array(diff)

# 正規化ジッター
jitter = np.mean(diff) / np.mean(T_sub)
print("0.05秒以降のジッター:", jitter)
