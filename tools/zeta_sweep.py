#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib import cm
from matplotlib.colors import Normalize


# ============================================================
# Plot style
# ============================================================

def setup_plot_style():
    """
    Use SciencePlots if available.
    Fall back to standard Matplotlib style if unavailable.
    """
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "ieee", "no-latex"])
        print("[INFO] Using SciencePlots style: science, ieee, no-latex")
    except Exception as e:
        print(f"[WARN] SciencePlots is unavailable. Falling back to default style.")
        print(f"[WARN] Reason: {e}")
        plt.style.use("default")


# ============================================================
# Path settings
# ============================================================

PROJECT_DIR = Path("/home/kajiyama/code/Mode_superposition_2026")

EXECUTABLE = PROJECT_DIR / "build" / "simulation"
PARAM_FILE = PROJECT_DIR / "input" / "param.txt"
DISP_DATA = PROJECT_DIR / "output" / "displace.dat"

SAVE_DIR = Path("./damping_sweep")


# ============================================================
# Sweep condition
# ============================================================

zetaR = 0.10

zetaL_list = [
    0.025,
    0.05,
    0.075,
    0.10,
    0.15,
    0.20,
    0.30,
    0.50,
]


# ============================================================
# Analysis settings
# ============================================================

steady_start = 0.05

# Waveform and phase portrait display window
plot_tmin = 0.05
plot_tmax = 0.08


# ============================================================
# Utility functions
# ============================================================

def check_paths():
    """
    Check required files before running the sweep.
    """
    if not EXECUTABLE.exists():
        raise FileNotFoundError(f"Executable not found: {EXECUTABLE}")

    if not os.access(EXECUTABLE, os.X_OK):
        raise PermissionError(
            f"Executable exists but is not executable: {EXECUTABLE}\n"
            f"Run: chmod +x {EXECUTABLE}"
        )

    if not PARAM_FILE.exists():
        raise FileNotFoundError(f"Parameter file not found: {PARAM_FILE}")

    output_dir = DISP_DATA.parent
    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory not found: {output_dir}")


def update_zeta(param_file, zetaL, zetaR):
    """
    Update zetaL and zetaR in param.txt.

    Expected format:

        # damping coefficient zetaL
        0.1

        # damping coefficient zetaR
        0.1
    """
    param_file = Path(param_file)

    with param_file.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    found_zetaL = False
    found_zetaR = False

    for i, line in enumerate(lines):
        if "# damping coefficient zetaL" in line:
            if i + 1 >= len(lines):
                raise ValueError(
                    f"Found zetaL label but no value line after it in {param_file}"
                )
            lines[i + 1] = f"{zetaL:.10g}\n"
            found_zetaL = True

        if "# damping coefficient zetaR" in line:
            if i + 1 >= len(lines):
                raise ValueError(
                    f"Found zetaR label but no value line after it in {param_file}"
                )
            lines[i + 1] = f"{zetaR:.10g}\n"
            found_zetaR = True

    if not found_zetaL:
        raise ValueError(
            f'Could not find "# damping coefficient zetaL" in {param_file}'
        )

    if not found_zetaR:
        raise ValueError(
            f'Could not find "# damping coefficient zetaR" in {param_file}'
        )

    with param_file.open("w", encoding="utf-8") as f:
        f.writelines(lines)


def run_simulation():
    """
    Run the external simulation.

    The simulation appears to use relative paths such as:
        ../input/...
        ../output/...

    Therefore, it should be executed from the build directory.
    """
    subprocess.run(
        [str(EXECUTABLE)],
        cwd=str(EXECUTABLE.parent),
        check=True,
    )


def load_displacement(disp_file):
    """
    Load displacement data.

    Expected columns:
        column 0: time
        column 1: left displacement
        column 2: right displacement
    """
    disp_file = Path(disp_file)

    if not disp_file.exists():
        raise FileNotFoundError(f"Displacement file not found: {disp_file}")

    data = np.loadtxt(disp_file, comments="#")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 3:
        raise ValueError(
            f"{disp_file} must have at least 3 columns: time, uL, uR. "
            f"Actual shape: {data.shape}"
        )

    t = data[:, 0]
    uL = -data[:, 1]
    uR = data[:, 2]

    return t, uL, uR


def calc_amplitude(sig):
    """
    Half peak-to-peak amplitude.
    """
    sig = np.asarray(sig)

    if sig.size == 0:
        return np.nan

    return 0.5 * (np.nanmax(sig) - np.nanmin(sig))


def calc_f0(time, signal, min_peak_distance_samples=20):
    """
    Estimate fundamental frequency from peak intervals.

    Parameters
    ----------
    time : array_like
        Time array.
    signal : array_like
        Signal array.
    min_peak_distance_samples : int
        Minimum peak distance in samples for scipy.signal.find_peaks.

    Returns
    -------
    float
        Estimated frequency in Hz. Returns np.nan if estimation fails.
    """
    time = np.asarray(time)
    signal = np.asarray(signal)

    if time.size < 3 or signal.size < 3:
        return np.nan

    if time.size != signal.size:
        return np.nan

    signal = signal - np.nanmean(signal)

    if not np.all(np.isfinite(signal)):
        return np.nan

    if np.nanmax(signal) - np.nanmin(signal) <= 0.0:
        return np.nan

    prominence = 0.05 * np.ptp(signal)

    peaks, _ = find_peaks(
        signal,
        distance=min_peak_distance_samples,
        prominence=prominence,
    )

    if len(peaks) < 2:
        return np.nan

    periods = np.diff(time[peaks])

    periods = periods[np.isfinite(periods)]
    periods = periods[periods > 0.0]

    if periods.size == 0:
        return np.nan

    return 1.0 / np.mean(periods)


def asymmetry_index(ampL, ampR):
    """
    Asymmetry index:
        (ampL - ampR) / (ampL + ampR)
    """
    denom = ampL + ampR

    if not np.isfinite(denom) or denom == 0.0:
        return np.nan

    return (ampL - ampR) / denom


def safe_ratio(numerator, denominator):
    """
    Safe division. Returns np.nan if denominator is zero or invalid.
    """
    if not np.isfinite(denominator) or denominator == 0.0:
        return np.nan

    return numerator / denominator


def make_zeta_filename_value(zeta):
    """
    Convert zeta value to file-safe string.
    Example:
        0.025 -> 0p025
    """
    return f"{zeta:.5g}".replace(".", "p")


def get_time_window(t, tmin, tmax):
    """
    Return boolean mask for requested time window.
    If no points exist in the requested window, return all True as fallback.
    """
    idx = (t > tmin) & (t < tmax)

    if np.count_nonzero(idx) < 2:
        print(
            f"[WARN] Not enough data in time window {tmin} < t < {tmax}. "
            f"Using full time range for plotting instead."
        )
        idx = np.ones_like(t, dtype=bool)

    return idx


# ============================================================
# Plot helpers
# ============================================================

def format_axis(ax):
    ax.tick_params(
        direction="in",
        which="both",
        top=True,
        right=True,
    )

    ax.grid(
        True,
        linestyle="--",
        alpha=0.35,
        linewidth=0.6,
    )

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


def save_summary_figure(results, zetaR, save_dir):
    """
    Summary figure.

    results columns:
        0 zetaL
        1 ampL
        2 ampR
        3 amp_ratio
        4 asymmetry
        5 F0
    """
    save_dir = Path(save_dir)

    if results.size == 0:
        print("[WARN] No results available. Summary figure was not created.")
        return

    x = results[:, 0] / zetaR

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(4.2, 6.0),
        dpi=200,
        sharex=True,
    )

    # Amplitude ratio
    axes[0].plot(
        x,
        results[:, 3],
        "o-",
        color="#0072B2",
        lw=1.6,
        ms=4,
        markerfacecolor="white",
        markeredgewidth=1.1,
    )

    axes[0].axhline(
        1.0,
        color="black",
        ls="--",
        lw=0.8,
        alpha=0.7,
    )

    axes[0].set_ylabel(r"$A_L/A_R$")
    axes[0].set_title("Effect of left-side damping", pad=8)

    # Asymmetry index
    axes[1].plot(
        x,
        results[:, 4],
        "s-",
        color="#D55E00",
        lw=1.6,
        ms=4,
        markerfacecolor="white",
        markeredgewidth=1.1,
    )

    axes[1].axhline(
        0.0,
        color="black",
        ls="--",
        lw=0.8,
        alpha=0.7,
    )

    axes[1].set_ylabel("Asymmetry")

    # F0
    axes[2].plot(
        x,
        results[:, 5],
        "^-",
        color="#009E73",
        lw=1.6,
        ms=4,
        markerfacecolor="white",
        markeredgewidth=1.1,
    )

    axes[2].set_ylabel(r"$F_0$ [Hz]")
    axes[2].set_xlabel(r"$\zeta_L/\zeta_R$")

    if np.all(x > 0.0):
        axes[2].set_xscale("log")
    else:
        print("[WARN] Non-positive x value found. Log scale was not applied.")

    for ax in axes:
        format_axis(ax)

    fig.align_ylabels(axes)
    fig.tight_layout()

    png_path = save_dir / "summary_damping_sweep.png"
    pdf_path = save_dir / "summary_damping_sweep.pdf"

    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)

    plt.close(fig)

    print(f"[INFO] Saved: {png_path}")
    print(f"[INFO] Saved: {pdf_path}")


def save_waveform_stack(waveform_data, zetaR, save_dir):
    """
    Save stacked waveform plots.

    waveform_data:
        list of (zetaL, t, uL, uR)
    """
    save_dir = Path(save_dir)

    if len(waveform_data) == 0:
        print("[WARN] No waveform data available. Waveform figure was not created.")
        return

    n = len(waveform_data)

    fig, axes = plt.subplots(
        n,
        1,
        figsize=(7.0, max(2.0, 1.0 * n)),
        dpi=200,
        sharex=True,
    )

    if n == 1:
        axes = [axes]

    ratios = np.array([z / zetaR for z, _, _, _ in waveform_data])
    norm = Normalize(vmin=np.min(ratios), vmax=np.max(ratios))
    cmap = cm.viridis

    for ax, data in zip(axes, waveform_data):
        zetaL, t, uL, uR = data

        idx = get_time_window(t, plot_tmin, plot_tmax)

        tt = t[idx]
        uLs = uL[idx]
        uRs = uR[idx]

        if tt.size >= 1:
            tt = (tt - tt[0]) * 1000.0

        uLs = uLs - np.mean(uLs)
        uRs = uRs - np.mean(uRs)

        c = cmap(norm(zetaL / zetaR))

        ax.plot(
            tt,
            uLs,
            color="#0072B2",
            lw=1.0,
            label="Left",
        )

        ax.plot(
            tt,
            uRs,
            color="#D55E00",
            lw=1.0,
            ls="--",
            label="Right",
        )

        ax.text(
            0.01,
            0.78,
            rf"$\zeta_L/\zeta_R={zetaL / zetaR:.2f}$",
            transform=ax.transAxes,
            fontsize=8,
            color=c,
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
                pad=1.5,
            ),
        )

        ax.set_ylabel("Disp.")
        format_axis(ax)

    axes[0].legend(
        loc="upper right",
        frameon=True,
        fontsize=8,
    )

    axes[-1].set_xlabel("Time [ms]")

    fig.tight_layout()

    png_path = save_dir / "waveform_stack.png"
    pdf_path = save_dir / "waveform_stack.pdf"

    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)

    plt.close(fig)

    print(f"[INFO] Saved: {png_path}")
    print(f"[INFO] Saved: {pdf_path}")


def save_phase_portrait_stack(waveform_data, zetaR, save_dir):
    """
    Save phase portraits: Left displacement vs Right displacement.
    """
    save_dir = Path(save_dir)

    if len(waveform_data) == 0:
        print("[WARN] No waveform data available. Phase portrait was not created.")
        return

    n = len(waveform_data)

    fig, axes = plt.subplots(
        1,
        n,
        figsize=(max(4.0, 1.25 * n), 1.8),
        dpi=200,
        sharex=False,
        sharey=False,
    )

    if n == 1:
        axes = [axes]

    ratios = np.array([z / zetaR for z, _, _, _ in waveform_data])
    norm = Normalize(vmin=np.min(ratios), vmax=np.max(ratios))
    cmap = cm.viridis

    for ax, data in zip(axes, waveform_data):
        zetaL, t, uL, uR = data

        idx = get_time_window(t, plot_tmin, plot_tmax)

        uLs = uL[idx] - np.mean(uL[idx])
        uRs = uR[idx] - np.mean(uR[idx])

        c = cmap(norm(zetaL / zetaR))

        ax.plot(
            uLs,
            uRs,
            color=c,
            lw=0.9,
        )

        ax.set_title(
            rf"{zetaL / zetaR:.2f}",
            fontsize=8,
        )

        ax.tick_params(
            direction="in",
            labelsize=7,
            top=True,
            right=True,
        )

        ax.grid(
            True,
            linestyle="--",
            alpha=0.25,
            linewidth=0.5,
        )

        ax.set_xlabel("Left")

    axes[0].set_ylabel("Right disp.")

    fig.suptitle(
        r"Phase portraits by $\zeta_L/\zeta_R$",
        y=1.08,
        fontsize=10,
    )

    fig.tight_layout()

    png_path = save_dir / "phase_portraits.png"
    pdf_path = save_dir / "phase_portraits.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")

    plt.close(fig)

    print(f"[INFO] Saved: {png_path}")
    print(f"[INFO] Saved: {pdf_path}")


# ============================================================
# Main sweep
# ============================================================

def main():
    setup_plot_style()

    SAVE_DIR.mkdir(parents=True, exist_ok=True)

    check_paths()

    results = []
    waveform_data = []

    for zetaL in zetaL_list:
        print("=" * 60)
        print(f"zetaL = {zetaL:.5g}, zetaR = {zetaR:.5g}")
        print("=" * 60)

        # ----------------------------------------------------
        # Update parameter file
        # ----------------------------------------------------
        update_zeta(
            PARAM_FILE,
            zetaL,
            zetaR,
        )

        # ----------------------------------------------------
        # Run simulation
        # ----------------------------------------------------
        run_simulation()

        # ----------------------------------------------------
        # Load displacement data
        # ----------------------------------------------------
        t, uL, uR = load_displacement(DISP_DATA)

        # ----------------------------------------------------
        # Save raw displacement data for this zetaL
        # ----------------------------------------------------
        zeta_name = make_zeta_filename_value(zetaL)
        save_name = SAVE_DIR / f"disp_zetaL_{zeta_name}.dat"

        shutil.copy(DISP_DATA, save_name)
        print(f"[INFO] Copied displacement data to: {save_name}")

        # ----------------------------------------------------
        # Extract steady-state region
        # ----------------------------------------------------
        idx = t > steady_start

        if np.count_nonzero(idx) < 2:
            print(
                f"[WARN] Not enough data after steady_start={steady_start}. "
                f"Skipping metrics for zetaL={zetaL:.5g}."
            )
            continue

        ts = t[idx]
        uLs = uL[idx]
        uRs = uR[idx]

        # ----------------------------------------------------
        # Metrics
        # ----------------------------------------------------
        ampL = calc_amplitude(uLs)
        ampR = calc_amplitude(uRs)

        amp_ratio = safe_ratio(ampL, ampR)

        asym = asymmetry_index(
            ampL,
            ampR,
        )

        f0 = calc_f0(
            ts,
            uLs,
        )

        results.append([
            zetaL,
            ampL,
            ampR,
            amp_ratio,
            asym,
            f0,
        ])

        waveform_data.append(
            (zetaL, t, uL, uR)
        )

        print(f"ampL       = {ampL:.6g}")
        print(f"ampR       = {ampR:.6g}")
        print(f"amp ratio  = {amp_ratio:.6g}")
        print(f"asym index = {asym:.6g}")
        print(f"F0         = {f0:.6g} Hz")

    # --------------------------------------------------------
    # Save CSV
    # --------------------------------------------------------
    if len(results) == 0:
        print("[ERROR] No valid results were obtained.")
        return

    results = np.array(results, dtype=float)

    csv_path = SAVE_DIR / "damping_sweep.csv"

    np.savetxt(
        csv_path,
        results,
        delimiter=",",
        header="zetaL,ampL,ampR,amp_ratio,asymmetry,F0",
        comments="",
    )

    print(f"[INFO] Saved: {csv_path}")

    # --------------------------------------------------------
    # Visualization
    # --------------------------------------------------------
    save_summary_figure(
        results,
        zetaR,
        SAVE_DIR,
    )

    save_waveform_stack(
        waveform_data,
        zetaR,
        SAVE_DIR,
    )

    save_phase_portrait_stack(
        waveform_data,
        zetaR,
        SAVE_DIR,
    )

    print("=" * 60)
    print("Damping sweep finished successfully.")
    print(f"Output directory: {SAVE_DIR.resolve()}")
    print("=" * 60)


if __name__ == "__main__":
    main()