
import numpy as np

def compare_resistance(x_mm, area_mm2, gap_cubed_mm4,
                       rho=1.13, mu=1.83e-5, ap_length_mm=17.0):
    # SI単位に変換
    x = np.asarray(x_mm) * 1e-3
    A = np.asarray(area_mm2) * 1e-6
    G = np.asarray(gap_cubed_mm4) * 1e-12

    # 開放状態で、全断面が有効な場合を比較
    if np.any(A <= 0) or np.any(G <= 0):
        raise ValueError("閉鎖・無効断面は別途処理が必要です")

    Amin = np.min(A[1:-1])
    length = x[-1] - x[0]
    ap_length = ap_length_mm * 1e-3

    # 論文公開版：最小断面積による集中近似
    L_old = 0.5 * rho * length / Amin
    Rv_old = 12 * mu * ap_length**2 * length / Amin**3

    # 現在のC++版：50断面での分布積分
    dx = np.diff(x)
    A_avg = (A[:-1] + A[1:]) / 2
    G_avg = (G[:-1] + G[1:]) / 2

    L_new = 0.5 * rho * np.sum(dx / A_avg)
    Rv_new = 12 * mu * np.sum(dx / G_avg)

    return {
        "Lg_ratio": L_new / L_old,
        "Rv_ratio": Rv_new / Rv_old,
        "Lg_old": L_old,
        "Lg_new": L_new,
        "Rv_old": Rv_old,
        "Rv_new": Rv_new,
    }


from pathlib import Path
import matplotlib.pyplot as plt

# C++シミュレーションの出力先
data_dir = Path("output/dat")
dt = 1e-5  # C++の計算時間刻み [s]

# データ読み込み
area_data = np.loadtxt(data_dir / "/home/kajiyama/code/simulation-main/output/dat/area.dat")
x_data = np.loadtxt(data_dir / "/home/kajiyama/code/simulation-main/output/dat/section_x.dat")
g3_data = np.loadtxt(data_dir / "/home/kajiyama/code/simulation-main/output/dat/gap_cubed.dat")

# ステップ番号の一致を確認
assert np.array_equal(area_data[:, 0], x_data[:, 0])
assert np.array_equal(area_data[:, 0], g3_data[:, 0])

steps = area_data[:, 0]
times = steps * dt

# 1列目を除いた50断面のデータ
areas = area_data[:, 1:]
xs = x_data[:, 1:]
g3s = g3_data[:, 1:]

Lg_ratio = []
Rv_ratio = []

for x, A, G in zip(xs, areas, g3s):
    try:
        result = compare_resistance(
            x, A, G,
            rho=1.13,
            mu=1.83e-5,
            ap_length_mm=17.0
        )
        Lg_ratio.append(result["Lg_ratio"])
        Rv_ratio.append(result["Rv_ratio"])

    except ValueError:
        # 声門閉鎖時などは比較から除外
        Lg_ratio.append(np.nan)
        Rv_ratio.append(np.nan)

# 時系列グラフ
plt.figure(figsize=(10, 5))
plt.plot(times, Lg_ratio, label="Lg: Current / Original")
plt.plot(times, Rv_ratio, label="Rv: Current / Original")

plt.axhline(1.0, color="gray", linestyle="--")
plt.xlabel("Time [s]")
plt.ylabel("Resistance coefficient ratio")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
