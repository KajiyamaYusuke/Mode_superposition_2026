#!/usr/bin/env python3
"""Independent left/right displacement cycle and ratio analysis."""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal


@dataclass
class RatioConfig:
    left_column: str = "uyL_mm"
    right_column: str = "uyR_mm"
    peak_prominence_ratio: float = 0.10
    period_ratio_target: float = 2.0
    period_ratio_tolerance: float = 0.20
    amplitude_ratio_target: float = 2.0
    amplitude_ratio_tolerance: float = 0.20
    min_peak_count: int = 6
    min_cycle_count: int = 5
    min_valid_amplitude: float = 1.0e-4
    max_period_cv_for_locking: float = 0.20
    freq_min_hz: float = 40.0
    freq_max_hz: float = 500.0
    smooth: str = "savgol"
    smooth_window: int = 11
    smooth_polyorder: int = 3


SUMMARY_COLUMNS = [
    "pressure_pa",
    "analysis_start_s",
    "analysis_end_s",
    "left_column",
    "right_column",
    "left_peak_count",
    "right_peak_count",
    "left_cycle_count",
    "right_cycle_count",
    "left_period_median_s",
    "right_period_median_s",
    "left_period_mean_s",
    "right_period_mean_s",
    "left_period_std_s",
    "right_period_std_s",
    "left_period_cv",
    "right_period_cv",
    "period_ratio_left_over_right",
    "period_ratio_unsigned",
    "period_ratio_near_2",
    "period_ratio_class",
    "left_amplitude_median_mm",
    "right_amplitude_median_mm",
    "left_amplitude_mean_mm",
    "right_amplitude_mean_mm",
    "left_amplitude_std_mm",
    "right_amplitude_std_mm",
    "left_amplitude_cv",
    "right_amplitude_cv",
    "amplitude_ratio_left_over_right",
    "amplitude_ratio_unsigned",
    "amplitude_ratio_near_2",
    "amplitude_ratio_class",
    "left_dominant_frequency_hz",
    "right_dominant_frequency_hz",
    "frequency_ratio_left_over_right",
    "frequency_ratio_unsigned",
    "status",
    "warning",
]


def _header_columns(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            tokens = line[1:].strip().split()
            if len(tokens) >= 2:
                return [re.sub(r"\[.*?\]", "", token) for token in tokens]
    raise ValueError(f"{path} has no column header")


def load_displacements(
    run_dir: Path, left_column: str, right_column: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = run_dir / "displace.dat"
    if not path.is_file():
        raise FileNotFoundError(path)
    columns = _header_columns(path)
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(columns):
        raise ValueError(
            f"{path}: header has {len(columns)} names but data has {data.shape[1]} columns"
        )
    aliases = {"time": "time_s"}
    columns = [aliases.get(name, name) for name in columns]
    missing = [name for name in ("time_s", left_column, right_column) if name not in columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}; available: {columns}")
    selected = data[:, [columns.index("time_s"), columns.index(left_column), columns.index(right_column)]]
    if not np.all(np.isfinite(selected)):
        raise ValueError(f"{path} contains NaN or Inf in selected columns")
    if np.any(np.diff(selected[:, 0]) <= 0):
        raise ValueError(f"{path} time is not strictly increasing")
    return selected[:, 0], selected[:, 1], selected[:, 2]


def preprocess_signal(
    time: np.ndarray,
    values: np.ndarray,
    start_time: float | None,
    end_time: float | None,
    config: RatioConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    start = time[0] if start_time is None else start_time
    end = time[-1] if end_time is None else end_time
    mask = (time >= start) & (time <= end)
    if np.count_nonzero(mask) < 5:
        raise ValueError("left/right analysis interval contains fewer than five samples")
    selected_time = time[mask]
    centered = values[mask] - np.mean(values[mask])
    filtered = centered.copy()
    if config.smooth != "none":
        window = min(config.smooth_window, len(centered) if len(centered) % 2 else len(centered) - 1)
        if window % 2 == 0:
            window -= 1
        if window >= 3:
            if config.smooth == "savgol":
                filtered = signal.savgol_filter(
                    centered, window, min(config.smooth_polyorder, window - 1)
                )
            elif config.smooth == "moving_average":
                filtered = np.convolve(centered, np.ones(window) / window, mode="same")
            else:
                raise ValueError(f"unknown left/right smoothing method: {config.smooth}")
    return selected_time, centered, filtered


def estimate_dominant_frequency(
    time: np.ndarray, values: np.ndarray, freq_min_hz: float, freq_max_hz: float
) -> float:
    sample_dt = float(np.median(np.diff(time)))
    frequencies, psd = signal.welch(
        values, fs=1.0 / sample_dt, nperseg=min(4096, len(values))
    )
    mask = (frequencies >= freq_min_hz) & (frequencies <= freq_max_hz)
    if not np.any(mask) or np.max(psd[mask]) <= 0:
        return np.nan
    indices = np.flatnonzero(mask)
    return float(frequencies[indices[np.argmax(psd[mask])]])


def detect_cycles(
    time: np.ndarray,
    raw_values: np.ndarray,
    filtered_values: np.ndarray,
    dominant_frequency: float,
    config: RatioConfig,
) -> tuple[np.ndarray, pd.DataFrame]:
    sample_dt = float(np.median(np.diff(time)))
    distance = (
        max(1, int(0.6 / (dominant_frequency * sample_dt)))
        if np.isfinite(dominant_frequency) and dominant_frequency > 0
        else 1
    )
    prominence = config.peak_prominence_ratio * float(np.ptp(filtered_values))
    peaks, _ = signal.find_peaks(
        filtered_values,
        distance=distance,
        prominence=prominence if prominence > 0 else None,
    )
    periods = np.diff(time[peaks])
    median_period = float(np.median(periods)) if len(periods) else np.nan
    valid = (
        (periods > 0.5 * median_period) & (periods < 1.5 * median_period)
        if np.isfinite(median_period)
        else np.zeros(len(periods), dtype=bool)
    )
    rows: list[dict[str, Any]] = []
    for index in range(max(0, len(peaks) - 1)):
        start, stop = int(peaks[index]), int(peaks[index + 1])
        trough_offset = int(np.argmin(raw_values[start : stop + 1]))
        trough_index = start + trough_offset
        rows.append(
            {
                "cycle_index": index,
                "peak_time_s": float(time[start]),
                "period_raw_s": float(periods[index]),
                "period_valid": bool(valid[index]),
                "peak_value_mm": float(raw_values[start]),
                "trough_value_mm": float(raw_values[trough_index]),
                "amplitude_mm": float(raw_values[start] - raw_values[trough_index]),
            }
        )
    return peaks, pd.DataFrame(rows)


def _stats(values: np.ndarray) -> dict[str, float]:
    if not len(values):
        return {"median": np.nan, "mean": np.nan, "std": np.nan, "cv": np.nan}
    mean = float(np.mean(values))
    return {
        "median": float(np.median(values)),
        "mean": mean,
        "std": float(np.std(values)),
        "cv": float(np.std(values) / mean) if mean != 0 else np.nan,
    }


def calculate_unsigned_ratio(left: float, right: float) -> tuple[float, float]:
    if not np.isfinite(left) or not np.isfinite(right) or left <= 0 or right <= 0:
        return np.nan, np.nan
    signed = left / right
    return float(signed), float(max(signed, 1.0 / signed))


def classify_period_ratio(
    signed: float,
    unsigned: float,
    left_cv: float,
    right_cv: float,
    frequency_unsigned: float,
    config: RatioConfig,
) -> tuple[bool, str]:
    if not np.all(np.isfinite([signed, unsigned, left_cv, right_cv])):
        return False, "insufficient_data"
    near_two = abs(unsigned - config.period_ratio_target) <= config.period_ratio_tolerance
    if 0.8 <= signed <= 1.2:
        return False, "same_period"
    if near_two:
        frequency_consistent = (
            np.isfinite(frequency_unsigned)
            and abs(frequency_unsigned - config.period_ratio_target)
            <= config.period_ratio_tolerance
        )
        if (
            left_cv > config.max_period_cv_for_locking
            or right_cv > config.max_period_cv_for_locking
            or not frequency_consistent
        ):
            return True, "ratio_near_2_but_irregular"
        return (
            True,
            "left_period_about_twice_right"
            if signed > 1
            else "right_period_about_twice_left",
        )
    return False, "other"


def classify_amplitude_ratio(
    signed: float, unsigned: float, config: RatioConfig
) -> tuple[bool, str]:
    if not np.isfinite(signed) or not np.isfinite(unsigned):
        return False, "insufficient_data"
    if 0.8 <= signed <= 1.2:
        return False, "same_scale"
    near_two = abs(unsigned - config.amplitude_ratio_target) <= config.amplitude_ratio_tolerance
    if near_two:
        return (
            True,
            "left_amplitude_about_twice_right"
            if signed > 1
            else "right_amplitude_about_twice_left",
        )
    return False, "other"


def _read_pressure(run_dir: Path) -> float:
    path = run_dir / "params_used.txt"
    if not path.is_file():
        return np.nan
    lines = path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if "subglottal pressure Ps" in line:
            for value in lines[index + 1 :]:
                value = value.strip()
                if value and not value.startswith("#"):
                    return float(value.split()[0])
    return np.nan


def _ratio_figures(
    run_dir: Path,
    time: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    left_peaks: np.ndarray,
    right_peaks: np.ndarray,
    left_cycles: pd.DataFrame,
    right_cycles: pd.DataFrame,
) -> None:
    figures = run_dir / "figures"
    figures.mkdir(exist_ok=True)
    fig, axis = plt.subplots(figsize=(10, 4))
    axis.plot(time, left, label="Left", linewidth=1)
    axis.plot(time, right, label="Right", linewidth=1)
    axis.scatter(time[left_peaks], left[left_peaks], s=16, marker="o")
    axis.scatter(time[right_peaks], right[right_peaks], s=18, marker="x")
    axis.set(xlabel="Time [s]", ylabel="Displacement [mm]", title="Left/Right Displacement Peaks")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(figures / "left_right_displacement.png", dpi=200)

    fig, axis = plt.subplots(figsize=(9, 4))
    for cycles, label in ((left_cycles, "Left"), (right_cycles, "Right")):
        if not cycles.empty:
            axis.plot(cycles["peak_time_s"], cycles["period_raw_s"], marker=".", label=label)
    axis.set(xlabel="Peak time [s]", ylabel="Period [s]", title="Left/Right Cycle Period")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(figures / "left_right_cycle_period.png", dpi=200)

    fig, axis = plt.subplots(figsize=(9, 4))
    for cycles, label in ((left_cycles, "Left"), (right_cycles, "Right")):
        if not cycles.empty:
            axis.plot(cycles["peak_time_s"], cycles["amplitude_mm"], marker=".", label=label)
    axis.set(xlabel="Peak time [s]", ylabel="Amplitude [mm]", title="Left/Right Cycle Amplitude")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(figures / "left_right_cycle_amplitude.png", dpi=200)
    plt.close("all")


def analyze_left_right_ratios(
    run_dir: Path,
    start_time: float | None,
    end_time: float | None,
    config: RatioConfig,
    save_plots: bool = True,
) -> dict[str, Any]:
    warning_messages: list[str] = []
    try:
        full_time, full_left, full_right = load_displacements(
            run_dir, config.left_column, config.right_column
        )
        time, left, left_filtered = preprocess_signal(
            full_time, full_left, start_time, end_time, config
        )
        right_time, right, right_filtered = preprocess_signal(
            full_time, full_right, start_time, end_time, config
        )
        if not np.array_equal(time, right_time):
            raise ValueError("left and right analysis time axes differ")
        left_frequency = estimate_dominant_frequency(
            time, left_filtered, config.freq_min_hz, config.freq_max_hz
        )
        right_frequency = estimate_dominant_frequency(
            time, right_filtered, config.freq_min_hz, config.freq_max_hz
        )
        left_peaks, left_cycles = detect_cycles(
            time, left, left_filtered, left_frequency, config
        )
        right_peaks, right_cycles = detect_cycles(
            time, right, right_filtered, right_frequency, config
        )
        left_cycles.to_csv(run_dir / "left_cycle_metrics.csv", index=False)
        right_cycles.to_csv(run_dir / "right_cycle_metrics.csv", index=False)

        left_periods = (
            left_cycles.loc[left_cycles["period_valid"], "period_raw_s"].to_numpy()
            if not left_cycles.empty
            else np.array([])
        )
        right_periods = (
            right_cycles.loc[right_cycles["period_valid"], "period_raw_s"].to_numpy()
            if not right_cycles.empty
            else np.array([])
        )
        left_amplitudes = left_cycles["amplitude_mm"].to_numpy() if not left_cycles.empty else np.array([])
        right_amplitudes = right_cycles["amplitude_mm"].to_numpy() if not right_cycles.empty else np.array([])
        left_period = _stats(left_periods)
        right_period = _stats(right_periods)
        left_amplitude = _stats(left_amplitudes)
        right_amplitude = _stats(right_amplitudes)
        period_signed, period_unsigned = calculate_unsigned_ratio(
            left_period["median"], right_period["median"]
        )
        amplitude_values_valid = (
            left_amplitude["median"] >= config.min_valid_amplitude
            and right_amplitude["median"] >= config.min_valid_amplitude
        )
        amplitude_signed, amplitude_unsigned = (
            calculate_unsigned_ratio(left_amplitude["median"], right_amplitude["median"])
            if amplitude_values_valid
            else (np.nan, np.nan)
        )
        frequency_signed, frequency_unsigned = calculate_unsigned_ratio(
            left_frequency, right_frequency
        )
        enough = (
            len(left_peaks) >= config.min_peak_count
            and len(right_peaks) >= config.min_peak_count
            and len(left_periods) >= config.min_cycle_count
            and len(right_periods) >= config.min_cycle_count
        )
        if not enough:
            warning_messages.append("too few valid left/right peaks or cycles")
        if not amplitude_values_valid:
            warning_messages.append("one or both representative amplitudes are below the noise floor")
        period_near, period_class = (
            classify_period_ratio(
                period_signed,
                period_unsigned,
                left_period["cv"],
                right_period["cv"],
                frequency_unsigned,
                config,
            )
            if enough
            else (False, "insufficient_data")
        )
        amplitude_near, amplitude_class = (
            classify_amplitude_ratio(amplitude_signed, amplitude_unsigned, config)
            if enough and amplitude_values_valid
            else (False, "insufficient_data")
        )
        summary = {
            "pressure_pa": _read_pressure(run_dir),
            "analysis_start_s": float(time[0]),
            "analysis_end_s": float(time[-1]),
            "left_column": config.left_column,
            "right_column": config.right_column,
            "left_peak_count": len(left_peaks),
            "right_peak_count": len(right_peaks),
            "left_cycle_count": len(left_periods),
            "right_cycle_count": len(right_periods),
            **{f"left_period_{key}_s" if key != "cv" else "left_period_cv": value for key, value in left_period.items()},
            **{f"right_period_{key}_s" if key != "cv" else "right_period_cv": value for key, value in right_period.items()},
            "period_ratio_left_over_right": period_signed,
            "period_ratio_unsigned": period_unsigned,
            "period_ratio_near_2": period_near,
            "period_ratio_class": period_class,
            **{f"left_amplitude_{key}_mm" if key != "cv" else "left_amplitude_cv": value for key, value in left_amplitude.items()},
            **{f"right_amplitude_{key}_mm" if key != "cv" else "right_amplitude_cv": value for key, value in right_amplitude.items()},
            "amplitude_ratio_left_over_right": amplitude_signed,
            "amplitude_ratio_unsigned": amplitude_unsigned,
            "amplitude_ratio_near_2": amplitude_near,
            "amplitude_ratio_class": amplitude_class,
            "left_dominant_frequency_hz": left_frequency,
            "right_dominant_frequency_hz": right_frequency,
            "frequency_ratio_left_over_right": frequency_signed,
            "frequency_ratio_unsigned": frequency_unsigned,
            "status": "completed" if enough else "insufficient_data",
            "warning": "; ".join(warning_messages),
        }
        if save_plots:
            _ratio_figures(
                run_dir,
                time,
                left,
                right,
                left_peaks,
                right_peaks,
                left_cycles,
                right_cycles,
            )
    except Exception as error:
        summary = {
            "pressure_pa": _read_pressure(run_dir),
            "analysis_start_s": start_time,
            "analysis_end_s": end_time,
            "left_column": config.left_column,
            "right_column": config.right_column,
            "status": "insufficient_data",
            "warning": str(error),
        }
        warnings.warn(f"left/right ratio analysis skipped: {error}", stacklevel=2)
    normalized = {column: summary.get(column, np.nan) for column in SUMMARY_COLUMNS}
    normalized["status"] = summary["status"]
    normalized["warning"] = summary["warning"]
    pd.DataFrame([normalized]).to_csv(
        run_dir / "left_right_ratio_metrics.csv", index=False
    )
    return normalized


def write_sweep_ratio_outputs(sweep_dir: Path) -> None:
    rows: list[pd.DataFrame] = []
    for path in sorted(sweep_dir.glob("pressure_*/left_right_ratio_metrics.csv")):
        rows.append(pd.read_csv(path))
    if not rows:
        return
    data = pd.concat(rows, ignore_index=True)
    output = sweep_dir / "pressure_sweep_left_right_ratios.csv"
    data.to_csv(output, index=False)
    figures = sweep_dir / "figures"
    figures.mkdir(exist_ok=True)
    plot_specs = (
        ("period_ratio_unsigned", "Unsigned period ratio [-]", "pressure_vs_period_ratio.png"),
        ("amplitude_ratio_unsigned", "Unsigned amplitude ratio [-]", "pressure_vs_amplitude_ratio.png"),
    )
    for column, ylabel, filename in plot_specs:
        if column not in data:
            continue
        valid = data.dropna(subset=["pressure_pa", column])
        if valid.empty:
            continue
        fig, axis = plt.subplots(figsize=(8, 5))
        axis.scatter(valid["pressure_pa"], valid[column], s=28)
        axis.axhline(1.0, color="gray", linestyle="--", linewidth=1, label="ratio = 1")
        axis.axhline(2.0, color="red", linestyle="--", linewidth=1, label="ratio = 2")
        axis.set(xlabel="Pressure [Pa]", ylabel=ylabel)
        axis.grid(alpha=0.3)
        axis.legend()
        fig.tight_layout()
        fig.savefig(figures / filename, dpi=220)
        plt.close(fig)
