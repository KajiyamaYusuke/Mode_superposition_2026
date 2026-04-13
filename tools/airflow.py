import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','ieee', 'no-latex'])

# 1. datファイルを読み込む
harea = np.loadtxt("../output/airflow_vt.dat")

labels = harea[:, 0]        # 1列目（時間ステップなど）
values_m3 = harea[:, 1]     # 2列目（流量 [m^3/s]）

# ★ここで単位変換します
# [m^3/s] を 1,000,000倍して [ml/s] に変換
values_ml = values_m3 * 1e6 

x = labels * 1e-5
start_time = 0.205
duration   = 0.03

plt.figure(figsize=(8,2), dpi=100)

# 変換後の values_ml をプロット
plt.plot(x, values_ml, linestyle='-', color='cornflowerblue', alpha = 0.8)

plt.xlabel('Time [s]', fontsize=15)
plt.ylabel('Airflow [ml/s]', fontsize=15)

plt.grid(True)

# 表示範囲の設定（ml/s になったので桁が変わります）
# 元が 0.001 m^3/s だったら 1000 ml/s です
plt.xlim(start_time, start_time + duration)
plt.ylim(0, 700) # ※データの振幅に合わせて調整してください

plt.tight_layout()
# plt.legend() # labelを設定していないのでコメントアウトしました
plt.savefig("result_airflow.png", dpi=300)
plt.show()