#!/usr/bin/env python3
"""Export and visualize irregular-vibration indicators without classifying states."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal


SIGNALS = ("x_left_mm", "x_right_mm", "glottal_area_min_mm2")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate left/right and glottal irregularity indicators."
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=PROJECT_ROOT / "output",
        help="simulation result directory (default: repository output directory)",
    )
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--prominence-ratio", type=float, default=0.05)
    parser.add_argument("--min-frequency", type=float, default=40.0)
    parser.add_argument("--max-frequency", type=float, default=1000.0)
    parser.add_argument("--max-repeat-lag", type=int, default=10)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def load_timeseries(run_dir: Path, start: float | None, end: float | None) -> pd.DataFrame:
    path = run_dir / "irregularity_timeseries.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} is missing. Rebuild and rerun the simulation to create it: "
            "cmake --build build -j && (cd build && ./simulation ../input/param.txt)"
        )
    frame = pd.read_csv(path)
    required = {"step", "time_s", *SIGNALS, "flow_rate_m3_s"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["time_s", *SIGNALS])
    if start is not None:
        frame = frame[frame["time_s"] >= start]
    if end is not None:
        frame = frame[frame["time_s"] <= end]
    frame = frame.reset_index(drop=True)
    if len(frame) < 8 or np.any(np.diff(frame["time_s"]) <= 0):
        raise ValueError("analysis interval needs at least 8 strictly ordered samples")
    return frame


def peak_table(
    time: np.ndarray, values: np.ndarray, ratio: float, max_frequency: float | None = None
) -> pd.DataFrame:
    centered = values - np.mean(values)
    prominence = ratio * np.ptp(centered)
    dt = float(np.median(np.diff(time)))
    distance = (
        max(1, int(0.5 / (max_frequency * dt)))
        if max_frequency is not None and max_frequency > 0
        else None
    )
    peaks, properties = signal.find_peaks(
        centered, prominence=prominence if prominence > 0 else None, distance=distance
    )
    rows = []
    for index, peak in enumerate(peaks):
        previous = peaks[index - 1] if index else None
        period = time[peak] - time[previous] if previous is not None else np.nan
        if previous is None:
            amplitude = np.nan
            trough_value = np.nan
        else:
            trough = previous + int(np.argmin(values[previous : peak + 1]))
            trough_value = values[trough]
            amplitude = values[peak] - trough_value
        rows.append(
            {
                "peak_index": index,
                "sample_index": int(peak),
                "peak_time_s": float(time[peak]),
                "peak_value": float(values[peak]),
                "trough_value": float(trough_value),
                "peak_to_trough_amplitude": float(amplitude),
                "period_s": float(period),
                "prominence": float(properties["prominences"][index]),
            }
        )
    return pd.DataFrame(rows)


def stats(values: np.ndarray) -> tuple[float, float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    mean = float(np.mean(values))
    std = float(np.std(values))
    return mean, std, std / abs(mean) if mean else np.nan


def repeat_errors(values: np.ndarray, max_lag: int) -> list[dict[str, float]]:
    values = values[np.isfinite(values)]
    scale = np.ptp(values)
    result = []
    for lag in range(1, max_lag + 1):
        error = float(np.mean(np.abs(values[lag:] - values[:-lag]))) if len(values) > lag else np.nan
        result.append(
            {"lag_cycles": lag, "mean_absolute_error": error,
             "normalized_error": error / scale if scale > 0 else np.nan}
        )
    return result


def spectrum(time: np.ndarray, values: np.ndarray) -> pd.DataFrame:
    dt = float(np.median(np.diff(time)))
    window = np.hanning(len(values))
    amplitude = np.abs(np.fft.rfft((values - np.mean(values)) * window))
    amplitude *= 2.0 / max(float(np.sum(window)), 1.0)
    return pd.DataFrame(
        {"frequency_hz": np.fft.rfftfreq(len(values), dt), "amplitude": amplitude}
    )


def dominant_frequency(table: pd.DataFrame, low: float, high: float) -> float:
    selected = table[(table.frequency_hz >= low) & (table.frequency_hz <= high)]
    return float(selected.loc[selected.amplitude.idxmax(), "frequency_hz"]) if len(selected) else np.nan


def phase_data(time: np.ndarray, left: np.ndarray, right: np.ndarray) -> tuple[pd.DataFrame, dict[str, float]]:
    left_phase = np.unwrap(np.angle(signal.hilbert(left - np.mean(left))))
    right_phase = np.unwrap(np.angle(signal.hilbert(right - np.mean(right))))
    unwrapped = right_phase - left_phase
    wrapped = np.angle(np.exp(1j * unwrapped))
    slope, intercept = np.polyfit(time, unwrapped, 1)
    slip_coordinate = (unwrapped - unwrapped[0]) / (2.0 * np.pi)
    slip_count = int(np.sum(np.abs(np.diff(np.floor(slip_coordinate + 0.5)))))
    duration = time[-1] - time[0]
    table = pd.DataFrame(
        {"time_s": time, "phase_left_rad": left_phase, "phase_right_rad": right_phase,
         "phase_diff_wrapped_rad": wrapped, "phase_diff_unwrapped_rad": unwrapped}
    )
    phase_coherence = min(max(abs(np.mean(np.exp(1j * wrapped))), 1e-15), 1.0)
    metrics = {
        "mean_phase_diff_rad": float(np.angle(np.mean(np.exp(1j * wrapped)))),
        "std_phase_diff_rad": float(np.sqrt(-2.0 * np.log(phase_coherence))),
        "phase_diff_range_rad": float(np.ptp(wrapped)),
        "phase_drift_rate_rad_s": float(slope),
        "phase_drift_intercept_rad": float(intercept),
        "phase_slip_count": slip_count,
        "phase_slip_rate_hz": slip_count / duration if duration > 0 else np.nan,
    }
    return table, metrics


def autocorrelation(time: np.ndarray, values: np.ndarray, name: str) -> pd.DataFrame:
    centered = values - np.mean(values)
    corr = signal.correlate(centered, centered, mode="full", method="fft")[len(values) - 1 :]
    if corr[0] != 0:
        corr = corr / corr[0]
    return pd.DataFrame({"lag_s": np.arange(len(corr)) * np.median(np.diff(time)), name: corr})


def save_plots(run_dir: Path, frame: pd.DataFrame, peaks: dict[str, pd.DataFrame],
               phases: pd.DataFrame, spectra: dict[str, pd.DataFrame],
               repeats: pd.DataFrame, correlations: pd.DataFrame,
               max_frequency: float, show: bool) -> None:
    figure_dir = run_dir / "figures" / "irregularity"
    figure_dir.mkdir(parents=True, exist_ok=True)
    time = frame.time_s.to_numpy()

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for axis, name in zip(axes, SIGNALS):
        axis.plot(time, frame[name], linewidth=0.8)
        table = peaks[name]
        if len(table):
            axis.scatter(table.peak_time_s, table.peak_value, s=10, color="tab:red")
        axis.set_ylabel(name)
        axis.grid(alpha=0.25)
    axes[-1].set_xlabel("Time [s]")
    fig.tight_layout(); fig.savefig(figure_dir / "timeseries_and_peaks.png", dpi=180)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for name, color in (("x_left_mm", "tab:blue"), ("x_right_mm", "tab:orange")):
        axes[0].plot(peaks[name].peak_time_s, peaks[name].period_s, ".-", label=name, color=color)
        axes[1].plot(peaks[name].peak_time_s, peaks[name].peak_to_trough_amplitude, ".-", label=name, color=color)
    axes[0].set_ylabel("Period [s]"); axes[1].set_ylabel("Amplitude [mm]")
    axes[1].set_xlabel("Peak time [s]")
    for axis in axes: axis.grid(alpha=0.25); axis.legend()
    fig.tight_layout(); fig.savefig(figure_dir / "cycle_sequences.png", dpi=180)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(frame.x_left_mm, frame.x_right_mm, linewidth=0.6)
    axes[0].set(xlabel="x_left [mm]", ylabel="x_right [mm]", title="Left-right trajectory")
    axes[1].plot(phases.time_s, phases.phase_diff_unwrapped_rad, linewidth=0.8)
    axes[1].set(xlabel="Time [s]", ylabel="Unwrapped phase difference [rad]", title="Phase drift and slips")
    for axis in axes: axis.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(figure_dir / "phase_relationship.png", dpi=180)

    fig, axis = plt.subplots(figsize=(10, 5))
    for name, table in spectra.items():
        axis.semilogy(table.frequency_hz, np.maximum(table.amplitude, 1e-15), label=name)
    axis.set(xlabel="Frequency [Hz]", ylabel="Amplitude", title="Spectra")
    axis.set_xlim(0.0, max_frequency)
    axis.grid(alpha=0.25); axis.legend(); fig.tight_layout()
    fig.savefig(figure_dir / "spectra.png", dpi=180)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for (name, quantity), group in repeats.groupby(["signal", "quantity"]):
        axes[0].plot(group.lag_cycles, group.normalized_error, ".-", label=f"{name} {quantity}")
    axes[0].set(xlabel="Lag [cycles]", ylabel="Normalized repeat error", title="Cycle repeatability")
    for name in SIGNALS:
        axes[1].plot(correlations.lag_s, correlations[name], label=name)
    axes[1].set(xlabel="Lag [s]", ylabel="Autocorrelation", title="Autocorrelation")
    for axis in axes: axis.grid(alpha=0.25); axis.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(figure_dir / "repeatability_autocorrelation.png", dpi=180)
    if show: plt.show()
    else: plt.close("all")


def main() -> None:
    args = arguments()
    run_dir = args.run_dir.resolve()
    frame = load_timeseries(run_dir, args.start_time, args.end_time)
    time = frame.time_s.to_numpy()
    peaks = {
        name: peak_table(
            time, frame[name].to_numpy(), args.prominence_ratio, args.max_frequency
        )
        for name in SIGNALS
    }
    for name, table in peaks.items():
        table.to_csv(run_dir / f"irregularity_peaks_{name}.csv", index=False)

    phases, phase_metrics = phase_data(time, frame.x_left_mm.to_numpy(), frame.x_right_mm.to_numpy())
    phases.to_csv(run_dir / "irregularity_phase.csv", index=False)
    spectra = {name: spectrum(time, frame[name].to_numpy()) for name in SIGNALS}
    pd.concat(spectra, names=["signal"]).reset_index(level=0).to_csv(
        run_dir / "irregularity_spectra.csv", index=False
    )

    repeat_rows = []
    metrics: dict[str, float | int] = {
        "analysis_start_s": float(time[0]), "analysis_end_s": float(time[-1]),
        "sample_interval_s": float(np.median(np.diff(time))), **phase_metrics,
    }
    for name in SIGNALS:
        table = peaks[name]
        periods = table.period_s.to_numpy()
        amplitudes = table.peak_to_trough_amplitude.to_numpy()
        for label, values in (("period", periods), ("amplitude", amplitudes)):
            mean, std, cv = stats(values)
            metrics[f"{name}_{label}_mean"] = mean
            metrics[f"{name}_{label}_std"] = std
            metrics[f"{name}_{label}_cv"] = cv
            for row in repeat_errors(values, args.max_repeat_lag):
                repeat_rows.append({"signal": name, "quantity": label, **row})
        metrics[f"{name}_peak_count"] = len(table)
        metrics[f"{name}_dominant_frequency_hz"] = dominant_frequency(
            spectra[name], args.min_frequency, args.max_frequency
        )
    left_count, right_count = metrics["x_left_mm_peak_count"], metrics["x_right_mm_peak_count"]
    metrics["peak_count_right_over_left"] = right_count / left_count if left_count else np.nan

    repeats = pd.DataFrame(repeat_rows)
    repeats.to_csv(run_dir / "irregularity_repeatability.csv", index=False)
    correlations = pd.DataFrame({"lag_s": autocorrelation(time, frame[SIGNALS[0]].to_numpy(), SIGNALS[0]).lag_s})
    for name in SIGNALS:
        correlations[name] = autocorrelation(time, frame[name].to_numpy(), name)[name]
    correlations.to_csv(run_dir / "irregularity_autocorrelation.csv", index=False)
    json_metrics = {
        key: (None if isinstance(value, float) and not np.isfinite(value) else value)
        for key, value in metrics.items()
    }
    (run_dir / "irregularity_metrics.json").write_text(
        json.dumps(json_metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    pd.DataFrame([metrics]).to_csv(run_dir / "irregularity_metrics.csv", index=False)
    save_plots(
        run_dir, frame, peaks, phases, spectra, repeats, correlations,
        args.max_frequency, args.show
    )
    print(f"Wrote irregularity indicators and figures to {run_dir}")


if __name__ == "__main__":
    main()
