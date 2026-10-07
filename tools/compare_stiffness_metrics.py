#!/usr/bin/env python3
"""Compare cycle metrics for multiple vocal-fold stiffness simulation runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dynamics_common import read_key_value_file
from left_right_ratios import (
    RatioConfig,
    detect_cycles,
    estimate_dominant_frequency,
    load_displacements,
    preprocess_signal,
)
from output_layout import output_path, result_path
from steinecke_classification import classify_locking


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare jitter, shimmer, frequency, and amplitude by stiffness condition."
    )
    parser.add_argument(
        "runs", nargs="*", metavar="LABEL=RUN_DIR",
        help="Condition and simulation directory; default: current=output",
    )
    parser.add_argument("--duration", type=float, default=0.15)
    parser.add_argument("--prominence-fraction", type=float, default=0.05)
    parser.add_argument("--minimum-peaks", type=int, default=6)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "output")
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def parse_runs(specifications: list[str]) -> list[tuple[str, Path]]:
    if not specifications:
        return [("current", PROJECT_ROOT / "output")]
    runs = []
    labels = set()
    for specification in specifications:
        if "=" not in specification:
            raise ValueError(f"run must have LABEL=RUN_DIR form: {specification}")
        label, directory = specification.split("=", 1)
        label = label.strip()
        if not label or label in labels:
            raise ValueError(f"run labels must be non-empty and unique: {label!r}")
        labels.add(label)
        path = Path(directory).expanduser()
        runs.append((label, path if path.is_absolute() else (Path.cwd() / path).resolve()))
    return runs


def local_variation(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if len(values) < 2 or not np.all(np.isfinite(values)):
        return np.nan
    mean = float(np.mean(values))
    return float(np.mean(np.abs(np.diff(values))) / mean) if mean > 0.0 else np.nan


def distribution_stats(values: np.ndarray, prefix: str) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            f"{prefix}_mean": np.nan, f"{prefix}_median": np.nan,
            f"{prefix}_std": np.nan, f"{prefix}_cv": np.nan,
        }
    mean = float(np.mean(values))
    std = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
    return {
        f"{prefix}_mean": mean,
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_std": std,
        f"{prefix}_cv": std / mean if mean != 0.0 else np.nan,
    }


def bootstrap_ratio(
    left: np.ndarray, right: np.ndarray, *, samples: int, rng: np.random.Generator
) -> dict[str, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left[np.isfinite(left) & (left > 0.0)]
    right = right[np.isfinite(right) & (right > 0.0)]
    empty = {"ratio": np.nan, "bootstrap_std": np.nan, "ci95_low": np.nan, "ci95_high": np.nan}
    if not len(left) or not len(right):
        return empty
    ratio = float(np.median(left) / np.median(right))
    if samples < 2:
        return {**empty, "ratio": ratio}
    left_draws = rng.choice(left, size=(samples, len(left)), replace=True)
    right_draws = rng.choice(right, size=(samples, len(right)), replace=True)
    ratios = np.median(left_draws, axis=1) / np.median(right_draws, axis=1)
    return {
        "ratio": ratio,
        "bootstrap_std": float(np.std(ratios, ddof=1)),
        "ci95_low": float(np.percentile(ratios, 2.5)),
        "ci95_high": float(np.percentile(ratios, 97.5)),
    }


def _cycle_arrays(cycles: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    if cycles.empty:
        return np.array([]), np.array([])
    valid = cycles["period_valid"].astype(bool).to_numpy()
    periods = cycles.loc[valid, "period_raw_s"].to_numpy(dtype=float)
    amplitudes = cycles.loc[valid, "amplitude_mm"].to_numpy(dtype=float)
    return periods, amplitudes


def analyze_run(
    label: str, run_dir: Path, *, duration: float, prominence_fraction: float,
    minimum_peaks: int, bootstrap_samples: int, seed: int,
) -> tuple[dict, pd.DataFrame]:
    full_time, full_left, full_right = load_displacements(run_dir, "uyL_mm", "uyR_mm")
    analysis_end = float(full_time[-1])
    analysis_start = max(float(full_time[0]), analysis_end - duration)
    config = RatioConfig(peak_prominence_ratio=prominence_fraction)
    time, left, left_filtered = preprocess_signal(
        full_time, full_left, analysis_start, analysis_end, config
    )
    _, right, right_filtered = preprocess_signal(
        full_time, full_right, analysis_start, analysis_end, config
    )
    left_fft_frequency = estimate_dominant_frequency(
        time, left_filtered, config.freq_min_hz, config.freq_max_hz
    )
    right_fft_frequency = estimate_dominant_frequency(
        time, right_filtered, config.freq_min_hz, config.freq_max_hz
    )
    left_peaks, left_cycles = detect_cycles(
        time, left, left_filtered, left_fft_frequency, config
    )
    right_peaks, right_cycles = detect_cycles(
        time, right, right_filtered, right_fft_frequency, config
    )
    left_periods, left_amplitudes = _cycle_arrays(left_cycles)
    right_periods, right_amplitudes = _cycle_arrays(right_cycles)
    left_frequencies = 1.0 / left_periods if len(left_periods) else np.array([])
    right_frequencies = 1.0 / right_periods if len(right_periods) else np.array([])

    rng = np.random.default_rng(seed)
    frequency_ratio = bootstrap_ratio(
        left_frequencies, right_frequencies, samples=bootstrap_samples, rng=rng
    )
    amplitude_ratio = bootstrap_ratio(
        left_amplitudes, right_amplitudes, samples=bootstrap_samples, rng=rng
    )
    locking = classify_locking(
        time, -left, right, transient_end=analysis_start,
        prominence_fraction=prominence_fraction, min_peaks=minimum_peaks,
    )
    manifest = read_key_value_file(result_path(run_dir, "manifest.txt"))
    enough = (
        len(left_peaks) >= minimum_peaks
        and len(right_peaks) >= minimum_peaks
        and len(left_periods) >= 2
        and len(right_periods) >= 2
    )
    summary = {
        "condition": label,
        "run_directory": str(run_dir),
        "analysis_start_s": analysis_start,
        "analysis_end_s": analysis_end,
        "analysis_duration_s": analysis_end - analysis_start,
        "left_peak_count": len(left_peaks),
        "right_peak_count": len(right_peaks),
        "left_valid_cycle_count": len(left_periods),
        "right_valid_cycle_count": len(right_periods),
        **distribution_stats(left_periods, "left_period_s"),
        **distribution_stats(right_periods, "right_period_s"),
        **distribution_stats(left_frequencies, "left_frequency_hz"),
        **distribution_stats(right_frequencies, "right_frequency_hz"),
        "frequency_ratio_left_over_right": frequency_ratio["ratio"],
        "frequency_ratio_bootstrap_std": frequency_ratio["bootstrap_std"],
        "frequency_ratio_ci95_low": frequency_ratio["ci95_low"],
        "frequency_ratio_ci95_high": frequency_ratio["ci95_high"],
        **distribution_stats(left_amplitudes, "left_amplitude_mm"),
        **distribution_stats(right_amplitudes, "right_amplitude_mm"),
        "amplitude_ratio_left_over_right": amplitude_ratio["ratio"],
        "amplitude_ratio_bootstrap_std": amplitude_ratio["bootstrap_std"],
        "amplitude_ratio_ci95_low": amplitude_ratio["ci95_low"],
        "amplitude_ratio_ci95_high": amplitude_ratio["ci95_high"],
        "left_jitter_local": local_variation(left_periods),
        "right_jitter_local": local_variation(right_periods),
        "left_shimmer_local": local_variation(left_amplitudes),
        "right_shimmer_local": local_variation(right_amplitudes),
        "locking_label": locking["label"],
        "locking_count_label": locking.get("base_label", ""),
        "locking_count_agreement": locking.get("count_agreement", np.nan),
        "kc1": manifest.get("kc1", ""),
        "kc2": manifest.get("kc2", ""),
        "kc3": manifest.get("kc3", ""),
        "zetaL": manifest.get("zetaL", ""),
        "zetaR": manifest.get("zetaR", ""),
        "left_mode_vtu": manifest.get("left_mode_vtu", ""),
        "right_mode_vtu": manifest.get("right_mode_vtu", ""),
        "prominence_fraction": prominence_fraction,
        "bootstrap_samples": bootstrap_samples,
        "bootstrap_seed": seed,
        "status": "completed" if enough else "insufficient_cycles",
    }
    cycle_tables = []
    for side, cycles in (("L", left_cycles), ("R", right_cycles)):
        table = cycles.copy()
        table.insert(0, "side", side)
        table.insert(0, "condition", label)
        cycle_tables.append(table)
    return summary, pd.concat(cycle_tables, ignore_index=True)


def plot_comparison(summary: pd.DataFrame, destination: Path) -> None:
    labels = summary["condition"].astype(str).to_numpy()
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), dpi=150)
    axes[0, 0].plot(x, 100.0 * summary["left_jitter_local"], "o-", label="Left")
    axes[0, 0].plot(x, 100.0 * summary["right_jitter_local"], "s--", label="Right")
    axes[0, 0].set_ylabel("Local jitter [%]")
    axes[0, 1].plot(x, 100.0 * summary["left_shimmer_local"], "o-", label="Left")
    axes[0, 1].plot(x, 100.0 * summary["right_shimmer_local"], "s--", label="Right")
    axes[0, 1].set_ylabel("Local shimmer [%]")
    for axis, value, error, ylabel in (
        (axes[1, 0], "frequency_ratio_left_over_right", "frequency_ratio_bootstrap_std", "Frequency ratio L/R"),
        (axes[1, 1], "amplitude_ratio_left_over_right", "amplitude_ratio_bootstrap_std", "Amplitude ratio L/R"),
    ):
        axis.errorbar(x, summary[value], yerr=summary[error], fmt="o-", capsize=4)
        axis.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
        axis.set_ylabel(ylabel)
    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=25, ha="right")
        axis.grid(alpha=0.3)
    axes[0, 0].legend()
    axes[0, 1].legend()
    fig.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=250, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = arguments()
    if args.duration <= 0.0 or args.bootstrap_samples < 0:
        raise ValueError("--duration must be positive and --bootstrap-samples non-negative")
    summaries, cycles = [], []
    for index, (label, run_dir) in enumerate(parse_runs(args.runs)):
        summary, cycle_table = analyze_run(
            label, run_dir, duration=args.duration,
            prominence_fraction=args.prominence_fraction,
            minimum_peaks=args.minimum_peaks,
            bootstrap_samples=args.bootstrap_samples, seed=args.seed + index,
        )
        summaries.append(summary)
        cycles.append(cycle_table)
    summary_frame = pd.DataFrame(summaries)
    summary_path = output_path(args.output_dir, "stiffness_comparison_metrics.csv")
    cycles_path = output_path(args.output_dir, "stiffness_cycle_metrics.csv")
    figure_path = output_path(args.output_dir, "stiffness_comparison.png")
    summary_frame.to_csv(summary_path, index=False)
    pd.concat(cycles, ignore_index=True).to_csv(cycles_path, index=False)
    plot_comparison(summary_frame, figure_path)
    metadata = {
        "ratio_uncertainty": "independent bootstrap of left/right cycle medians",
        "jitter": "mean(abs(diff(period))) / mean(period)",
        "shimmer": "mean(abs(diff(peak-to-trough amplitude))) / mean(amplitude)",
    }
    output_path(args.output_dir, "stiffness_comparison_definition.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {summary_path}")
    print(f"Wrote {cycles_path}")
    print(f"Wrote {figure_path}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
