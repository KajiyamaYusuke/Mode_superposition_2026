import os
import subprocess
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import optuna

EXECUTABLE = "/home/kajiyama/code/Mode_superposition_2026/build/simulation"
PARAM_FILE = "/home/kajiyama/code/Mode_superposition_2026/input/param.txt"
DISP_FILE = "/home/kajiyama/code/Mode_superposition_2026/output/displace.dat"
AREA_FILE = "/home/kajiyama/code/Mode_superposition_2026/output/area.dat"
DEBUG_FILE = "/home/kajiyama/code/Mode_superposition_2026/output/debug_step_summary.csv"

TARGET_F0 = 196.0
TARGET_OQ = 0.5
TARGET_DISP_MM = 1.5
TARGET_GAW_MM2 = 10.0

SIM_DT = 1.0e-5
OUTPUT_INTERVAL = 5
DT = SIM_DT * OUTPUT_INTERVAL


def update_param_file(kc1, kc2, kc3, zeta=0.1, ps=800.0):
    with open(PARAM_FILE, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    numeric_count = 0

    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            new_lines.append(line)
            continue

        numeric_count += 1

        # param.txt 構造に対応
        if numeric_count == 6:      # zeta
            new_lines.append(f"{zeta}\n")
        elif numeric_count == 7:    # kc1 kc2 kc3
            new_lines.append(f"{kc1} {kc2} {kc3}\n")
        elif numeric_count == 15:   # ps
            new_lines.append(f"{ps}\n")
        else:
            new_lines.append(line)

    with open(PARAM_FILE, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


def run_sim():
    # shell=True でも動くならそのままでOK
    result = subprocess.run(EXECUTABLE, shell=True)
    return result.returncode


def load_disp_signal():
    data = np.loadtxt(DISP_FILE, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    t = data[:, 0]          # simulation.cpp 側ですでに n*dt [s]
    disp_l = data[:, 1]     # 左側 uy 変位
    return t, disp_l


def load_gaw_proxy_signal():
    # area.dat は各時刻に対して harea[i] の列
    data = np.loadtxt(AREA_FILE, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    step = data[:, 0]
    areas = data[:, 1:]

    # 厳密な GAW ではなく proxy
    gaw_proxy = np.max(areas, axis=1)

    t = step * DT
    return t, gaw_proxy


def compute_f0(signal, dt):
    x = signal - np.mean(signal)
    if len(x) < 16:
        return np.nan
    if np.allclose(np.std(x), 0.0):
        return np.nan

    freqs = np.fft.rfftfreq(len(x), d=dt)
    spec = np.abs(np.fft.rfft(x))

    # DC を除外
    if len(spec) < 2:
        return np.nan

    idx = np.argmax(spec[1:]) + 1
    return freqs[idx]


def compute_oq(gaw_proxy):
    x = gaw_proxy - np.min(gaw_proxy)
    peak = np.max(x)
    if peak <= 1e-12:
        return 0.0

    x = x / peak
    open_mask = x > 0.1
    return np.mean(open_mask)


def read_diverged():
    if not os.path.exists(DEBUG_FILE):
        return False

    try:
        df = pd.read_csv(DEBUG_FILE)
    except Exception:
        return True

    if len(df) == 0 or "diverged" not in df.columns:
        return False

    return int(df.iloc[-1]["diverged"]) == 1


def base_loss(f0, oq, gaw_mm2, disp_mm):
    return (
        ((f0 - TARGET_F0) / TARGET_F0) ** 2
        + (oq - TARGET_OQ) ** 2
        + (np.log(max(gaw_mm2, 1e-6) / TARGET_GAW_MM2)) ** 2
        + (np.log(max(disp_mm, 1e-6) / TARGET_DISP_MM)) ** 2
    )


def objective(trial):
    # ---- 探索範囲 ----
    kc1 = trial.suggest_float("kc1", 1e-2, 1e2, log=True)
    kc2 = trial.suggest_float("kc2", 1e-1, 1e4, log=True)
    kc3 = trial.suggest_float("kc3", 1e-1, 1e4, log=True)

    update_param_file(kc1, kc2, kc3, zeta=0.1, ps=800.0)

    rc = run_sim()
    if rc != 0:
        return 1e6

    if read_diverged():
        return 1e6

    try:
        t_disp, disp = load_disp_signal()
        t_gaw, gaw = load_gaw_proxy_signal()
    except Exception:
        return 1e6

    # ファイルが極端に短いときも失敗扱い
    if len(disp) < 20 or len(gaw) < 20:
        return 1e6

    # 過渡を捨てる
    start_d = len(disp) // 2
    start_g = len(gaw) // 2
    disp = disp[start_d:]
    gaw = gaw[start_g:]

    if len(disp) < 20 or len(gaw) < 20:
        return 1e6

    # 指標計算
    f0 = compute_f0(gaw, DT)
    oq = compute_oq(gaw)

    # 単位はあなたの現状解釈に合わせる
    disp_mm = np.max(np.abs(disp))
    gaw_mm2 = np.max(gaw)

    if not np.isfinite(f0):
        return 1e6

    penalty = 0.0

    # ---- 失敗ケース ----
    # 全く開かない
    if oq < 0.02:
        penalty += 300.0

    # 開口が極端に小さい
    if gaw_mm2 < 0.5:
        penalty += 200.0

    # 振幅が極端に小さい
    if disp_mm < 0.1:
        penalty += 120.0

    # 異常周波数
    if f0 < 80.0 or f0 > 400.0:
        penalty += 250.0

    # 周期性不足
    peaks, _ = find_peaks(gaw)
    if len(peaks) < 3:
        penalty += 180.0

    # 完全に異常なら即失格
    if oq <= 0.0 or oq >= 1.0 or f0 < 50.0 or f0 > 1000.0:
        print(
            f"[PENALTY-HARD] "
            f"kc1={kc1:.3g}, kc2={kc2:.3g}, kc3={kc3:.3g} | "
            f"f0={f0:.1f}, OQ={oq:.3f}, GAW={gaw_mm2:.3g}, DISP={disp_mm:.3g}"
        )
        return 1e6

    value = base_loss(f0, oq, gaw_mm2, disp_mm) + penalty

    print(
        f"[TRIAL] kc1={kc1:.3g}, kc2={kc2:.3g}, kc3={kc3:.3g} | "
        f"f0={f0:.1f}, OQ={oq:.3f}, GAW={gaw_mm2:.3g}, DISP={disp_mm:.3g}, "
        f"loss={value:.4f}"
    )

    return value


if __name__ == "__main__":
    study = optuna.create_study(direction="minimize")

    # 既知の良さそうな点を最初に評価
    study.enqueue_trial({"kc1": 5.0, "kc2": 6000.0, "kc3": 5000.0})

    study.optimize(objective, n_trials=10)

    print("Best trial:")
    print(study.best_trial.params)
    print("Best value:", study.best_trial.value)