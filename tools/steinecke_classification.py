"""Steinecke--Herzel locking classification for left/right displacement."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks, get_window


@dataclass
class PeakSeries:
    indices: np.ndarray
    times: np.ndarray
    values: np.ndarray


def extract_peaks(
    time: np.ndarray,
    values: np.ndarray,
    *,
    prominence_fraction: float = 0.05,
    min_distance_time: float | None = None,
) -> PeakSeries:
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if time.ndim != 1 or values.ndim != 1 or len(time) != len(values):
        raise ValueError("time and values must be matching one-dimensional arrays")
    if len(time) < 3 or not np.all(np.diff(time) > 0.0):
        raise ValueError("time must contain at least three strictly increasing values")
    signal_range = float(np.ptp(values))
    if not np.isfinite(signal_range) or signal_range <= 0.0:
        return PeakSeries(np.array([], dtype=int), np.array([]), np.array([]))
    kwargs: dict[str, float | int] = {
        "prominence": prominence_fraction * signal_range
    }
    if min_distance_time is not None:
        sample_interval = float(np.median(np.diff(time)))
        kwargs["distance"] = max(1, int(round(min_distance_time / sample_interval)))
    indices, _ = find_peaks(values, **kwargs)
    return PeakSeries(indices, time[indices], values[indices])


def detect_peak_period(
    peak_times: np.ndarray,
    peak_values: np.ndarray,
    *,
    signal_scale: float,
    max_period: int = 12,
    min_repeats: int = 3,
    amplitude_rtol: float = 0.02,
    interval_rtol: float = 0.02,
) -> tuple[int | None, dict[int, tuple[float, float]]]:
    times = np.asarray(peak_times, dtype=float)
    values = np.asarray(peak_values, dtype=float)
    errors: dict[int, tuple[float, float]] = {}
    if len(times) != len(values) or len(values) < min_repeats:
        return None, errors
    amplitude_scale = max(abs(float(signal_scale)), np.finfo(float).eps)
    median_interval = float(np.median(np.diff(times)))
    if median_interval <= 0.0:
        return None, errors
    intervals = np.r_[np.nan, np.diff(times)]
    for period in range(1, max_period + 1):
        if len(values) < (min_repeats + 1) * period:
            break
        amplitude_error = float(
            np.sqrt(np.mean((values[period:] - values[:-period]) ** 2))
            / amplitude_scale
        )
        valid = np.arange(period, len(values))
        valid = valid[(valid - period) > 0]
        if not len(valid):
            continue
        interval_error = float(
            np.sqrt(np.mean((intervals[valid] - intervals[valid - period]) ** 2))
            / median_interval
        )
        errors[period] = (amplitude_error, interval_error)
        if amplitude_error <= amplitude_rtol and interval_error <= interval_rtol:
            return period, errors
    return None, errors


def next_maximum_map(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    return (values[:-1], values[1:]) if len(values) >= 2 else (np.array([]), np.array([]))


def one_sided_spectrum(
    time: np.ndarray, values: np.ndarray, *, window: str = "hann"
) -> tuple[np.ndarray, np.ndarray]:
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    sample_interval = float(np.median(np.diff(time)))
    weights = get_window(window, len(values), fftbins=True)
    transformed = np.fft.rfft((values - np.mean(values)) * weights)
    return np.fft.rfftfreq(len(values), sample_interval), np.abs(transformed)


def map_cluster_count(values: np.ndarray, tolerance: float) -> int:
    """Count compact next-maximum-map point clouds with a deterministic rule."""
    map_x, map_y = next_maximum_map(values)
    if not len(map_x):
        return 0
    points = np.column_stack((map_x, map_y))
    clusters: list[list[np.ndarray]] = []
    for point in points:
        for cluster in clusters:
            center = np.mean(cluster, axis=0)
            if np.linalg.norm(point - center) <= tolerance:
                cluster.append(point)
                break
        else:
            clusters.append([point])
    return len(clusters)


def classify_locking(
    time: np.ndarray,
    x1l: np.ndarray,
    x1r: np.ndarray,
    *,
    transient_end: float,
    prominence_fraction: float = 0.05,
    min_distance_time: float | None = None,
    max_period: int = 12,
    min_peaks: int = 8,
    amplitude_rtol: float = 0.02,
    interval_rtol: float = 0.02,
    minimum_signal_range: float = 1.0e-8,
    count_agreement: float = 0.8,
    map_tolerance_fraction: float = 0.02,
) -> dict:
    time = np.asarray(time, dtype=float)
    x1l = np.asarray(x1l, dtype=float)
    x1r = np.asarray(x1r, dtype=float)
    if not (len(time) == len(x1l) == len(x1r)):
        raise ValueError("time, x1l, and x1r lengths must match")
    use = time >= transient_end
    analysis_time, left, right = time[use], x1l[use], x1r[use]
    result = {
        "label": "insufficient_data",
        "analysis_time": analysis_time,
        "left_signal": left,
        "right_signal": right,
        "thresholds": {
            "prominence_fraction": prominence_fraction,
            "max_period": max_period,
            "min_peaks": min_peaks,
            "amplitude_rtol": amplitude_rtol,
            "interval_rtol": interval_rtol,
            "minimum_signal_range": minimum_signal_range,
            "count_agreement_required": count_agreement,
            "map_tolerance_fraction": map_tolerance_fraction,
        },
    }
    if len(analysis_time) < 3:
        return result
    range_left, range_right = float(np.ptp(left)), float(np.ptp(right))
    result.update(signal_range_left=range_left, signal_range_right=range_right)
    if min(range_left, range_right) < minimum_signal_range:
        result["label"] = "no_or_insufficient_oscillation"
        return result
    peaks_left = extract_peaks(
        analysis_time, left, prominence_fraction=prominence_fraction,
        min_distance_time=min_distance_time,
    )
    peaks_right = extract_peaks(
        analysis_time, right, prominence_fraction=prominence_fraction,
        min_distance_time=min_distance_time,
    )
    result.update(peaks_left=peaks_left, peaks_right=peaks_right)
    if min(len(peaks_left.times), len(peaks_right.times)) < min_peaks:
        return result
    period_right, errors_right = detect_peak_period(
        peaks_right.times, peaks_right.values, signal_scale=range_right,
        max_period=max_period, amplitude_rtol=amplitude_rtol,
        interval_rtol=interval_rtol,
    )
    period_left, errors_left = detect_peak_period(
        peaks_left.times, peaks_left.values, signal_scale=range_left,
        max_period=max_period, amplitude_rtol=amplitude_rtol,
        interval_rtol=interval_rtol,
    )
    map_x, map_y = next_maximum_map(peaks_right.values)
    cluster_count = map_cluster_count(
        peaks_right.values,
        max(map_tolerance_fraction * range_right, np.finfo(float).eps),
    )
    result.update(
        period_right=period_right, period_left=period_left,
        period_errors_right=errors_right, period_errors_left=errors_left,
        next_maximum_x=map_x, next_maximum_y=map_y,
        map_cluster_count=cluster_count,
    )
    if period_right is None:
        result["label"] = "aperiodic_or_high_order"
        return result
    count_pairs, total_periods, windows = [], [], []
    for index in range(len(peaks_right.times) - period_right):
        start, end = peaks_right.times[index], peaks_right.times[index + period_right]
        n_right = int(np.sum((peaks_right.times >= start) & (peaks_right.times < end)))
        n_left = int(np.sum((peaks_left.times >= start) & (peaks_left.times < end)))
        count_pairs.append((n_right, n_left))
        total_periods.append(float(end - start))
        windows.append((start, end, n_right, n_left))
    pairs = np.asarray(count_pairs, dtype=int)
    unique_pairs, occurrences = np.unique(pairs, axis=0, return_counts=True)
    best = int(np.argmax(occurrences))
    n_right, n_left = map(int, unique_pairs[best])
    agreement = float(occurrences[best] / len(pairs))
    base_label = f"{n_right}:{n_left}"
    consistency = []
    if agreement < count_agreement:
        consistency.append("unstable_peak_counts")
    if cluster_count != period_right:
        consistency.append("next_maximum_cluster_mismatch")
    if period_left is not None and period_left != n_left:
        consistency.append("left_period_mismatch")
    label = base_label if not consistency else "uncertain"
    result.update(
        label=label, base_label=base_label, n_right=n_right, n_left=n_left,
        T_total=float(np.median(total_periods)), count_pairs=pairs,
        count_windows=np.asarray(windows), count_agreement=agreement,
        consistency_issues=consistency,
    )
    return result
