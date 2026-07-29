#!/usr/bin/env python3
"""Analyze periodicity and modulation in one simulation result directory."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dynamics_common import (
    cycle_table,
    detect_events,
    interpolate_at,
    json_value,
    load_min_area,
    load_optional_signals,
    local_variation,
    read_dt,
    select_interval,
    smooth_signal,
    spectrum_metrics,
    threshold_crossings,
    write_json,
)
from left_right_ratios import RatioConfig, analyze_left_right_ratios


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Quantify vocal-fold vibration dynamics.")
    parser.add_argument("--run-dir", type=Path, default=Path("output"))
    parser.add_argument("--dt", type=float)
    parser.add_argument("--start-time", "--analysis-start", dest="start_time", type=float)
    parser.add_argument("--end-time", "--analysis-end", dest="end_time", type=float)
    parser.add_argument("--last-cycles", type=int)
    parser.add_argument(
        "--event-method",
        choices=["rising_threshold", "peak", "trough", "hilbert"],
        default="rising_threshold",
    )
    parser.add_argument("--section-level", type=float, default=0.5)
    parser.add_argument("--direction", choices=["rising", "falling"], default="rising")
    parser.add_argument("--min-peak-distance", type=float)
    parser.add_argument("--prominence", type=float)
    parser.add_argument(
        "--smooth", choices=["none", "savgol", "moving_average"], default="savgol"
    )
    parser.add_argument("--smooth-window", type=int, default=7)
    parser.add_argument("--smooth-polyorder", type=int, default=2)
    parser.add_argument("--minimum-cycles", type=int, default=20)
    parser.add_argument("--left-column", default="uyL_mm")
    parser.add_argument("--right-column", default="uyR_mm")
    parser.add_argument("--peak-prominence-ratio", type=float, default=0.10)
    parser.add_argument("--period-ratio-target", type=float, default=2.0)
    parser.add_argument("--period-ratio-tolerance", type=float, default=0.20)
    parser.add_argument("--amplitude-ratio-target", type=float, default=2.0)
    parser.add_argument("--amplitude-ratio-tolerance", type=float, default=0.20)
    parser.add_argument("--min-peak-count", type=int, default=6)
    parser.add_argument("--min-valid-amplitude", type=float, default=1.0e-4)
    parser.add_argument("--allow-nan", action="store_true")
    parser.add_argument("--save-plots", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def classification(period_cv: float, amplitude_cv: float, cycles: pd.DataFrame) -> str:
    if len(cycles) < 5:
        return "insufficient_points"
    amplitudes = cycles["amplitude_mm2"].to_numpy()
    if period_cv < 0.01 and amplitude_cv < 0.01:
        return "period_1_candidate"
    if period_cv < 0.02 and amplitude_cv > max(0.03, 3 * period_cv):
        even, odd = amplitudes[::2], amplitudes[1::2]
        separation = abs(np.mean(even) - np.mean(odd))
        spread = np.std(amplitudes)
        if spread > 0 and separation > 1.5 * spread:
            return "period_2_candidate"
        return "amplitude_modulation_candidate"
    if period_cv < 0.08:
        return "quasiperiodic_or_multiperiod_candidate"
    return "scattered_or_chaotic_candidate"


def save_figures(
    run_dir: Path,
    time: np.ndarray,
    raw: np.ndarray,
    filtered: np.ndarray,
    event_times: np.ndarray,
    cycles: pd.DataFrame,
    spectrum: pd.DataFrame,
    poincare: pd.DataFrame,
    mode_energy: pd.DataFrame | None,
    show: bool,
) -> None:
    figures = run_dir / "figures"
    figures.mkdir(exist_ok=True)

    fig, axis = plt.subplots(figsize=(10, 4))
    axis.plot(time, raw, label="minimum area", linewidth=1)
    axis.plot(time, filtered, label="detection signal", linewidth=1, alpha=0.8)
    if len(event_times):
        axis.scatter(event_times, interpolate_at(time, raw, event_times), s=18, color="red", label="events")
    axis.set(xlabel="Time [s]", ylabel="Minimum area [mm²]", title="Detected Cycle Events")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(figures / "time_series.png", dpi=200)

    fig, axis = plt.subplots(figsize=(9, 4))
    if not cycles.empty:
        axis.plot(cycles["event_time_s"], cycles["period_s"], marker=".")
    axis.set(xlabel="Time [s]", ylabel="Period [s]", title="Cycle Period")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "cycle_period.png", dpi=200)

    fig, axis = plt.subplots(figsize=(9, 4))
    if not cycles.empty:
        axis.plot(cycles["cycle_index"], cycles["amplitude_mm2"], marker=".")
    axis.set(xlabel="Cycle index", ylabel="Amplitude [mm²]", title="Cycle Amplitude")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "amplitude_sequence.png", dpi=200)

    fig, axis = plt.subplots(figsize=(9, 4))
    axis.semilogy(spectrum["frequency_hz"], np.maximum(spectrum["fft_amplitude"], 1e-30))
    axis.set(xlabel="Frequency [Hz]", ylabel="FFT amplitude", title="Minimum-area Spectrum")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "spectrum.png", dpi=200)

    fig, axis = plt.subplots(figsize=(5, 5))
    amplitudes = cycles["amplitude_mm2"].to_numpy()
    if len(amplitudes) > 1:
        axis.scatter(amplitudes[:-1], amplitudes[1:], s=20)
    axis.set(xlabel=r"$a_n$ [mm²]", ylabel=r"$a_{n+1}$ [mm²]", title="Amplitude Return Map")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "return_map.png", dpi=200)

    fig, axis = plt.subplots(figsize=(5, 5))
    candidates = [name for name in poincare.columns if name != "crossing_index"]
    x_name, y_name = ("uyL_mm", "uyR_mm") if {"uyL_mm", "uyR_mm"} <= set(candidates) else ("time_s", "darea_dt_mm2_s")
    axis.scatter(poincare[x_name], poincare[y_name], c=poincare["crossing_index"], s=22)
    axis.set(xlabel=x_name, ylabel=y_name, title="Poincaré Points")
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figures / "poincare.png", dpi=200)

    if mode_energy is not None and not mode_energy.empty:
        fig, axis = plt.subplots(figsize=(10, 5))
        for (side, mode), group in mode_energy.groupby(["side", "mode_index"]):
            if int(mode) >= 5:
                continue
            axis.plot(
                group["time_s"],
                group["energy_like"],
                linewidth=0.9,
                label=f"{side} mode {int(mode) + 1}",
            )
        axis.set(
            xlabel="Time [s]",
            ylabel="Energy-like quantity",
            title="Leading Modal Energy-like Quantities",
        )
        axis.set_yscale("log")
        axis.grid(alpha=0.3)
        axis.legend(ncol=2, fontsize=8)
        fig.tight_layout()
        fig.savefig(figures / "mode_energy.png", dpi=200)

    if show:
        plt.show()
    else:
        plt.close("all")


def main() -> None:
    args = arguments()
    run_dir = args.run_dir.resolve()
    dt = read_dt(run_dir, args.dt)
    full_time, full_area, _ = load_min_area(run_dir, dt, args.allow_nan)
    time, area = select_interval(full_time, full_area, args.start_time, args.end_time)
    filtered = smooth_signal(area, args.smooth, args.smooth_window, args.smooth_polyorder)

    event_times, _ = detect_events(
        time,
        filtered,
        args.event_method,
        args.section_level,
        args.direction,
        args.min_peak_distance,
        args.prominence,
    )
    cycles = cycle_table(time, area, event_times)
    if args.last_cycles:
        cycles = cycles.tail(args.last_cycles).reset_index(drop=True)
        if not cycles.empty:
            event_times = event_times[event_times >= cycles["event_time_s"].iloc[0]]
    if len(cycles) < args.minimum_cycles:
        warnings.warn(
            f"only {len(cycles)} cycles detected; at least {args.minimum_cycles} are recommended",
            stacklevel=2,
        )

    periods = cycles["period_s"].to_numpy() if not cycles.empty else np.array([])
    amplitudes = cycles["amplitude_mm2"].to_numpy() if not cycles.empty else np.array([])
    period_mean = float(np.mean(periods)) if len(periods) else np.nan
    amplitude_mean = float(np.mean(amplitudes)) if len(amplitudes) else np.nan
    period_cv = float(np.std(periods) / period_mean) if period_mean else np.nan
    amplitude_cv = float(np.std(amplitudes) / amplitude_mean) if amplitude_mean else np.nan
    spectrum, spectral_metrics = spectrum_metrics(time, area)

    crossing_times, crossing_indices = threshold_crossings(
        time, filtered, args.section_level, args.direction
    )
    darea = np.gradient(area, time)
    poincare = load_optional_signals(run_dir, dt, crossing_times)
    poincare.insert(0, "crossing_index", np.arange(len(poincare)))
    poincare["area_mm2"] = interpolate_at(time, area, crossing_times)
    poincare["darea_dt_mm2_s"] = interpolate_at(time, darea, crossing_times)

    return_map = pd.DataFrame(
        {
            "cycle_index": cycles["cycle_index"].iloc[:-1].to_numpy() if len(cycles) > 1 else [],
            "amplitude_n_mm2": amplitudes[:-1],
            "amplitude_next_mm2": amplitudes[1:],
            "period_n_s": periods[:-1],
            "period_next_s": periods[1:],
        }
    )
    mode_energy: pd.DataFrame | None = None
    mode_energy_summary = pd.DataFrame()
    modal_path = run_dir / "modal_contribution.csv"
    if modal_path.is_file():
        modal = pd.read_csv(modal_path)
        required = {"time_s", "side", "mode_index", "frequency_hz", "q", "qdot"}
        if required <= set(modal.columns):
            modal = modal[
                (modal["time_s"] >= time[0]) & (modal["time_s"] <= time[-1])
            ].copy()
            omega = 2.0 * np.pi * modal["frequency_hz"].to_numpy()
            modal["energy_like"] = 0.5 * modal["qdot"].to_numpy() ** 2 + 0.5 * (
                omega * modal["q"].to_numpy()
            ) ** 2
            mode_energy = modal[
                ["time_s", "side", "mode_index", "frequency_hz", "energy_like"]
            ]
            mode_energy_summary = (
                mode_energy.groupby(["side", "mode_index", "frequency_hz"])[
                    "energy_like"
                ]
                .agg(["mean", "std"])
                .reset_index()
                .rename(columns={"mean": "mean_energy_like", "std": "std_energy_like"})
            )
            mode_energy_summary["energy_cv"] = (
                mode_energy_summary["std_energy_like"]
                / mode_energy_summary["mean_energy_like"]
            )
            mode_energy.to_csv(run_dir / "mode_energy_timeseries.csv", index=False)
            mode_energy_summary.to_csv(
                run_dir / "mode_energy_summary.csv", index=False
            )
            energy_wide = mode_energy.assign(
                mode_label=lambda frame: frame["side"]
                + "_mode_"
                + (frame["mode_index"] + 1).astype(int).astype(str)
            ).pivot_table(
                index="time_s", columns="mode_label", values="energy_like"
            )
            energy_wide.corr().to_csv(run_dir / "mode_energy_correlation.csv")
    metrics = {
        "analysis_start_s": float(time[0]),
        "analysis_end_s": float(time[-1]),
        "dt_s": dt,
        "event_method": args.event_method,
        "section_level": args.section_level,
        "number_of_cycles": len(cycles),
        "mean_period_s": period_mean,
        "f0_from_cycles_hz": 1.0 / period_mean if period_mean else np.nan,
        "period_std_s": float(np.std(periods)) if len(periods) else np.nan,
        "period_cv": period_cv,
        "jitter_local": local_variation(periods),
        "mean_amplitude_mm2": amplitude_mean,
        "amplitude_cv": amplitude_cv,
        "shimmer_local": local_variation(amplitudes),
        "mean_max_area_mm2": float(cycles["peak_area_mm2"].mean()) if len(cycles) else np.nan,
        "mean_min_area_mm2": float(cycles["trough_area_mm2"].mean()) if len(cycles) else np.nan,
        **spectral_metrics,
        "poincare_class": classification(period_cv, amplitude_cv, cycles),
        "classification_note": "heuristic classification; not a proof of chaos",
        "status": "completed" if len(cycles) >= args.minimum_cycles else "insufficient_cycles",
    }
    write_json(run_dir / "metrics.json", metrics)
    cycles.to_csv(run_dir / "cycle_metrics.csv", index=False)
    poincare.to_csv(run_dir / "poincare_points.csv", index=False)
    return_map.to_csv(run_dir / "return_map.csv", index=False)
    spectrum.to_csv(run_dir / "spectrum.csv", index=False)
    if args.save_plots:
        save_figures(
            run_dir,
            time,
            area,
            filtered,
            event_times,
            cycles,
            spectrum,
            poincare,
            mode_energy,
            args.show,
        )
    ratio_config = RatioConfig(
        left_column=args.left_column,
        right_column=args.right_column,
        peak_prominence_ratio=args.peak_prominence_ratio,
        period_ratio_target=args.period_ratio_target,
        period_ratio_tolerance=args.period_ratio_tolerance,
        amplitude_ratio_target=args.amplitude_ratio_target,
        amplitude_ratio_tolerance=args.amplitude_ratio_tolerance,
        min_peak_count=args.min_peak_count,
        min_valid_amplitude=args.min_valid_amplitude,
        smooth=args.smooth,
        smooth_window=max(args.smooth_window, 11),
        smooth_polyorder=max(args.smooth_polyorder, 3),
    )
    ratio_metrics = analyze_left_right_ratios(
        run_dir,
        args.start_time,
        args.end_time,
        ratio_config,
        save_plots=args.save_plots,
    )
    print(f"Analyzed {run_dir}: {len(cycles)} cycles, {metrics['poincare_class']}")
    print(f"period_cv={json_value(period_cv)}, amplitude_cv={json_value(amplitude_cv)}")
    print(
        "left/right ratios: "
        f"period={json_value(ratio_metrics.get('period_ratio_unsigned', np.nan))}, "
        f"amplitude={json_value(ratio_metrics.get('amplitude_ratio_unsigned', np.nan))}, "
        f"status={ratio_metrics['status']}"
    )


if __name__ == "__main__":
    main()
