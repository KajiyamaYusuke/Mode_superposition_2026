#!/usr/bin/env python3
"""Shared loading and signal-processing utilities for dynamics analysis."""

from __future__ import annotations

import json
import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import signal


def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw in path.read_text(encoding="utf-8").splitlines():
        if "=" in raw and not raw.lstrip().startswith("#"):
            key, value = raw.split("=", 1)
            values[key.strip()] = value.strip().strip('"')
    return values


def read_dt(run_dir: Path, explicit_dt: float | None = None) -> float:
    if explicit_dt is not None:
        if explicit_dt <= 0:
            raise ValueError("--dt must be positive")
        return explicit_dt
    manifest = read_key_value_file(run_dir / "manifest.txt")
    if "dt_s" in manifest:
        return float(manifest["dt_s"])
    params = run_dir / "params_used.txt"
    if params.is_file():
        lines = params.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            if "time step for time integration dt" in line:
                for value in lines[index + 1 :]:
                    value = value.strip()
                    if value and not value.startswith("#"):
                        return float(value.split()[0])
    raise ValueError("dt is unavailable; pass --dt or provide manifest.txt/params_used.txt")


def load_numeric(path: Path, minimum_columns: int = 2) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"required input does not exist: {path}")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] < minimum_columns:
        raise ValueError(f"{path} has {data.shape[1]} columns; expected at least {minimum_columns}")
    return np.asarray(data, dtype=float)


def validate_time_series(
    time: np.ndarray, values: np.ndarray, name: str, allow_nan: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    if len(time) != len(values) or len(time) < 3:
        raise ValueError(f"{name} must contain at least three matching samples")
    finite = np.isfinite(time) & np.all(
        np.isfinite(values.reshape(len(time), -1)), axis=1
    )
    if not np.all(finite):
        bad = np.flatnonzero(~finite).tolist()
        if not allow_nan:
            raise ValueError(f"{name} contains NaN/Inf at rows {bad[:10]}")
        warnings.warn(f"{name}: removing {len(bad)} non-finite rows", stacklevel=2)
        time, values = time[finite], values[finite]
    if np.any(np.diff(time) <= 0):
        raise ValueError(f"{name} time/step values are not strictly increasing")
    return time, values


def load_min_area(
    run_dir: Path, dt: float, allow_nan: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = load_numeric(run_dir / "area.dat")
    steps, sections = validate_time_series(data[:, 0], data[:, 1:], "area.dat", allow_nan)
    if np.any(sections < 0):
        warnings.warn("area.dat contains negative area values; values were not clipped", stacklevel=2)
    return steps * dt, np.min(sections, axis=1), steps


def select_interval(
    time: np.ndarray,
    values: np.ndarray,
    start_time: float | None,
    end_time: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    start = time[0] if start_time is None else start_time
    end = time[-1] if end_time is None else end_time
    if start >= end:
        raise ValueError("analysis start time must be smaller than end time")
    mask = (time >= start) & (time <= end)
    if np.count_nonzero(mask) < 3:
        raise ValueError("selected analysis interval contains fewer than three samples")
    return time[mask], values[mask]


def smooth_signal(
    values: np.ndarray,
    method: str = "savgol",
    window: int = 7,
    polyorder: int = 2,
) -> np.ndarray:
    if method == "none":
        return values.copy()
    if window < 3:
        raise ValueError("smoothing window must be at least 3")
    if window % 2 == 0:
        window += 1
    window = min(window, len(values) if len(values) % 2 else len(values) - 1)
    if window < 3:
        return values.copy()
    if method == "savgol":
        return signal.savgol_filter(values, window, min(polyorder, window - 1))
    if method == "moving_average":
        return np.convolve(values, np.ones(window) / window, mode="same")
    raise ValueError(f"unknown smoothing method: {method}")


def threshold_crossings(
    time: np.ndarray,
    values: np.ndarray,
    level: float = 0.5,
    direction: str = "rising",
) -> tuple[np.ndarray, np.ndarray]:
    low, high = np.percentile(values, [5, 95])
    if not high > low:
        return np.array([]), np.array([], dtype=int)
    threshold = low + level * (high - low)
    if direction == "rising":
        indices = np.flatnonzero((values[:-1] < threshold) & (values[1:] >= threshold))
    elif direction == "falling":
        indices = np.flatnonzero((values[:-1] > threshold) & (values[1:] <= threshold))
    else:
        raise ValueError("direction must be rising or falling")
    fractions = (threshold - values[indices]) / (values[indices + 1] - values[indices])
    return time[indices] + fractions * (time[indices + 1] - time[indices]), indices


def detect_events(
    time: np.ndarray,
    values: np.ndarray,
    method: str,
    section_level: float = 0.5,
    direction: str = "rising",
    min_peak_distance: float | None = None,
    prominence: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    sampling = float(np.median(np.diff(time)))
    distance = None
    if min_peak_distance is not None:
        distance = max(1, int(round(min_peak_distance / sampling)))
    if method in {"peak", "maximum"}:
        indices, _ = signal.find_peaks(values, distance=distance, prominence=prominence)
        return time[indices], indices
    if method in {"trough", "minimum"}:
        indices, _ = signal.find_peaks(-values, distance=distance, prominence=prominence)
        return time[indices], indices
    if method == "rising_threshold":
        return threshold_crossings(time, values, section_level, direction)
    if method == "hilbert":
        phase = np.unwrap(np.angle(signal.hilbert(values - np.mean(values))))
        turns = np.floor((phase - phase[0]) / (2 * np.pi))
        indices = np.flatnonzero(np.diff(turns) > 0)
        return time[indices + 1], indices + 1
    raise ValueError(f"unknown event method: {method}")


def cycle_table(
    time: np.ndarray, values: np.ndarray, event_times: np.ndarray
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for cycle, (start, end) in enumerate(zip(event_times[:-1], event_times[1:])):
        mask = (time >= start) & (time < end)
        if not np.any(mask):
            continue
        segment_t, segment = time[mask], values[mask]
        imax, imin = int(np.argmax(segment)), int(np.argmin(segment))
        maximum, minimum = float(segment[imax]), float(segment[imin])
        rows.append(
            {
                "cycle_index": cycle,
                "event_time_s": start,
                "period_s": end - start,
                "peak_time_s": segment_t[imax],
                "peak_area_mm2": maximum,
                "trough_time_s": segment_t[imin],
                "trough_area_mm2": minimum,
                "amplitude_mm2": maximum - minimum,
            }
        )
    return pd.DataFrame(rows)


def local_variation(values: np.ndarray) -> float:
    mean = float(np.mean(values)) if len(values) else np.nan
    return float(np.mean(np.abs(np.diff(values))) / mean) if len(values) > 1 and mean != 0 else np.nan


def spectrum_metrics(time: np.ndarray, values: np.ndarray) -> tuple[pd.DataFrame, dict[str, Any]]:
    dt = float(np.median(np.diff(time)))
    centered = signal.detrend(values)
    frequencies = np.fft.rfftfreq(len(centered), dt)
    amplitude = 2.0 * np.abs(np.fft.rfft(centered)) / len(centered)
    nperseg = min(4096, len(centered))
    welch_f, psd = signal.welch(centered, fs=1.0 / dt, nperseg=nperseg)
    positive = welch_f > 0
    peak_indices, _ = signal.find_peaks(psd[positive])
    candidates = np.flatnonzero(positive)[peak_indices]
    ranked = candidates[np.argsort(psd[candidates])[::-1]] if len(candidates) else np.array([], dtype=int)
    probability = psd / np.sum(psd) if np.sum(psd) > 0 else np.zeros_like(psd)
    nonzero = probability > 0
    entropy = -np.sum(probability[nonzero] * np.log(probability[nonzero]))
    entropy /= np.log(len(probability)) if len(probability) > 1 else 1.0
    metrics = {
        "spectral_entropy": float(entropy),
        "dominant_frequency_1_hz": float(welch_f[ranked[0]]) if len(ranked) else np.nan,
        "dominant_frequency_2_hz": float(welch_f[ranked[1]]) if len(ranked) > 1 else np.nan,
        "spectral_peak_count": int(len(candidates)),
    }
    return pd.DataFrame({"frequency_hz": frequencies, "fft_amplitude": amplitude}), metrics


def interpolate_at(source_time: np.ndarray, values: np.ndarray, target_time: np.ndarray) -> np.ndarray:
    return np.interp(target_time, source_time, values, left=np.nan, right=np.nan)


def load_optional_signals(run_dir: Path, dt: float, target_time: np.ndarray) -> pd.DataFrame:
    result = pd.DataFrame({"time_s": target_time})
    specs = {
        "airflow_m3_s": ("airflow_vt.dat", 1, True),
        "outlet_pressure_pa": ("pressure_vt.dat", 1, True),
        "uyL_mm": ("displace.dat", 1, False),
        "uyR_mm": ("displace.dat", 2, False),
    }
    cache: dict[str, np.ndarray] = {}
    for name, (filename, column, first_is_step) in specs.items():
        path = run_dir / filename
        if not path.is_file():
            continue
        data = cache.setdefault(filename, load_numeric(path, column + 1))
        source_time = data[:, 0] * dt if first_is_step else data[:, 0]
        result[name] = interpolate_at(source_time, data[:, column], target_time)
    modal_path = run_dir / "modal_contribution.csv"
    if modal_path.is_file():
        modal = pd.read_csv(modal_path)
        for side in ("L", "R"):
            for mode in sorted(modal["mode_index"].unique())[:2]:
                selected = modal[(modal["side"] == side) & (modal["mode_index"] == mode)]
                if selected.empty:
                    continue
                suffix = int(mode) + 1
                result[f"q{side}_{suffix}"] = interpolate_at(
                    selected["time_s"].to_numpy(), selected["q"].to_numpy(), target_time
                )
                result[f"qdot{side}_{suffix}"] = interpolate_at(
                    selected["time_s"].to_numpy(), selected["qdot"].to_numpy(), target_time
                )
    return result


def json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_value(item) for item in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(json_value(data), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def replace_pressure_parameter(text: str, pressure_pa: float) -> str:
    lines = text.splitlines()
    pattern = re.compile(r"subglottal pressure Ps", re.IGNORECASE)
    for index, line in enumerate(lines):
        if pattern.search(line):
            for value_index in range(index + 1, len(lines)):
                stripped = lines[value_index].strip()
                if stripped and not stripped.startswith("#"):
                    lines[value_index] = f"{pressure_pa:g}"
                    return "\n".join(lines) + "\n"
    raise ValueError("could not locate the value after '# subglottal pressure Ps (Pa)'")
