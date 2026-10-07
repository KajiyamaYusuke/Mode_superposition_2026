#!/usr/bin/env python3
"""Analyze one simulation run and classify its self-oscillation response."""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal

from output_layout import output_path, result_path


CLASSIFICATIONS = {
    "DECAYING_TO_STATIC", "GROWING_OSCILLATION", "SUSTAINED_PERIODIC",
    "SUSTAINED_IRREGULAR", "DIVERGED_OR_INVALID", "INSUFFICIENT_DATA",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--analysis-start", type=float, default=0.15)
    parser.add_argument("--analysis-end", type=float, default=0.50)
    parser.add_argument("--tail-duration", type=float, default=0.10)
    parser.add_argument("--band-low-hz", type=float, default=20.0)
    parser.add_argument("--band-high-hz", type=float, default=300.0)
    parser.add_argument("--growth-threshold-per-s", type=float, default=1.0)
    parser.add_argument("--decay-threshold-per-s", type=float, default=-1.0)
    parser.add_argument("--static-rms-threshold-mm", type=float, default=1.0e-4)
    parser.add_argument("--periodic-envelope-cv-threshold", type=float, default=0.15)
    parser.add_argument("--envelope-fit-trim-fraction", type=float, default=0.30,
                        help="fraction trimmed from each Hilbert-envelope edge")
    parser.add_argument("--closure-area-threshold-mm2", type=float, default=1.0e-6)
    parser.add_argument("--time-step-relative-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--minimum-samples", type=int, default=64)
    return parser.parse_args()


def finite_or_none(value: Any) -> Any:
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def read_dat(path: Path, names: list[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=names)


def load_optional_csv(run_dir: Path, name: str, notes: list[str]) -> pd.DataFrame:
    path = result_path(run_dir, name)
    if not path.is_file():
        notes.append(f"missing optional file: {name}")
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception as exception:  # keep partial analysis available
        notes.append(f"could not read {name}: {exception}")
        return pd.DataFrame()


def window(frame: pd.DataFrame, start: float, end: float) -> pd.DataFrame:
    if frame.empty or "time_s" not in frame:
        return frame.iloc[0:0]
    return frame[(frame["time_s"] >= start) & (frame["time_s"] <= end)].copy()


def statistics(values: np.ndarray, prefix: str) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {}
    mean = float(np.mean(values))
    return {
        f"{prefix}_mean": mean,
        f"{prefix}_max_amplitude": float(np.max(np.abs(values - mean))),
        f"{prefix}_peak_to_peak": float(np.ptp(values)),
        f"{prefix}_rms": float(np.sqrt(np.mean(values ** 2))),
        f"{prefix}_fluctuation_rms": float(np.sqrt(np.mean((values - mean) ** 2))),
    }


def dominant_and_growth(
    time: np.ndarray,
    values: np.ndarray,
    band_low: float,
    band_high: float,
    trim_fraction: float,
) -> tuple[float, float, float, np.ndarray, np.ndarray]:
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(time) & np.isfinite(values)
    time, values = time[finite], values[finite]
    if len(time) < 8:
        return np.nan, np.nan, np.nan, time, np.full_like(time, np.nan)
    dt = float(np.median(np.diff(time)))
    fs = 1.0 / dt
    detrended = signal.detrend(values, type="linear")
    frequency = np.fft.rfftfreq(len(detrended), dt)
    amplitude = np.abs(np.fft.rfft(detrended * signal.windows.hann(len(detrended))))
    usable = (frequency >= band_low) & (frequency <= min(band_high, 0.49 * fs))
    dominant = float(frequency[np.flatnonzero(usable)[np.argmax(amplitude[usable])]]) \
        if np.any(usable) else np.nan
    low = max(band_low / (0.5 * fs), 1.0e-6)
    high = min(band_high / (0.5 * fs), 0.999999)
    filtered = detrended
    if low < high and len(values) > 24:
        sos = signal.butter(3, [low, high], btype="bandpass", output="sos")
        try:
            filtered = signal.sosfiltfilt(sos, detrended)
        except ValueError:
            filtered = detrended
    envelope = np.abs(signal.hilbert(filtered))
    floor = max(float(np.nanmax(envelope)) * 1.0e-8, np.finfo(float).tiny)
    span = time[-1] - time[0]
    fit_start = time[0] + trim_fraction * span
    fit_end = time[-1] - trim_fraction * span
    valid = (np.isfinite(envelope) & (envelope > floor)
             & (time >= fit_start) & (time <= fit_end))
    sigma = float(np.polyfit(time[valid], np.log(envelope[valid]), 1)[0]) \
        if np.count_nonzero(valid) >= 8 else np.nan
    envelope_mean = float(np.mean(envelope[valid])) if np.any(valid) else np.nan
    envelope_cv = float(np.std(envelope[valid]) / envelope_mean) \
        if envelope_mean > 0 else np.nan
    return dominant, sigma, envelope_cv, time, envelope


def check_time(frame: pd.DataFrame, name: str, tolerance: float, notes: list[str]) -> None:
    if frame.empty or "time_s" not in frame or len(frame) < 3:
        return
    time = frame["time_s"].to_numpy(dtype=float)
    differences = np.diff(time)
    median = float(np.median(differences))
    if median <= 0 or np.any(differences <= 0):
        notes.append(f"{name}: non-increasing time values")
    elif np.max(np.abs(differences - median)) > tolerance * median:
        notes.append(f"{name}: inconsistent time step")


def main() -> None:
    args = parse_arguments()
    if not 0.0 <= args.envelope_fit_trim_fraction < 0.5:
        raise ValueError("--envelope-fit-trim-fraction must be in [0, 0.5)")
    if args.analysis_end <= args.analysis_start:
        raise ValueError("--analysis-end must be greater than --analysis-start")
    run_dir = args.run_directory.resolve()
    notes: list[str] = []
    invalid = False

    displacement_path = result_path(run_dir, "displace.dat")
    if displacement_path.is_file():
        displacement = read_dat(displacement_path, ["time_s", "left_mm", "right_mm"])
    else:
        displacement = pd.DataFrame()
        notes.append("missing required primary signal: displace.dat")
    irregularity = load_optional_csv(run_dir, "irregularity_timeseries.csv", notes)
    flow = load_optional_csv(run_dir, "flow_diagnostics_v2.csv", notes)
    flow_source_name = "flow_diagnostics_v2.csv"
    energy = load_optional_csv(run_dir, "energy_summary.csv", notes)
    area_history = pd.DataFrame()
    area_path = result_path(run_dir, "area.dat")
    if area_path.is_file():
        try:
            raw_area = pd.read_csv(area_path, sep=r"\s+", comment="#", header=None)
            if raw_area.shape[1] >= 4 and len(raw_area) == len(displacement):
                # Columns 1 and -1 are boundary planes, matching the solver's
                # physical constriction search over sections.area[1:-1].
                area_history = pd.DataFrame({
                    "time_s": displacement["time_s"].to_numpy(dtype=float),
                    "min_area_mm2": np.min(
                        raw_area.iloc[:, 2:-1].to_numpy(dtype=float), axis=1
                    ),
                })
            elif not raw_area.empty:
                notes.append("area.dat length does not match displace.dat")
        except Exception as exception:
            notes.append(f"could not read area.dat: {exception}")
    if flow.empty:
        legacy_flow_path = result_path(run_dir, "debug_flow_detail.csv")
        if legacy_flow_path.is_file():
            try:
                legacy_flow = pd.read_csv(legacy_flow_path, comment="#")
                legacy_flow = legacy_flow.rename(columns={
                    "time": "time_s", "minA_mm2": "min_area_mm2",
                    "currentPg": "current_pg_pa", "currentUg": "flow_rate_m3_s",
                    "hasNonFinite": "has_nonfinite",
                })
                required = {"time_s", "min_area_mm2", "current_pg_pa", "flow_rate_m3_s"}
                if required.issubset(legacy_flow.columns):
                    flow = legacy_flow
                    flow_source_name = "debug_flow_detail.csv"
                    notes.append("using legacy debug_flow_detail.csv fallback")
            except Exception as exception:
                notes.append(f"could not read legacy debug_flow_detail.csv: {exception}")

    for name, frame in (("displace.dat", displacement),
                        ("irregularity_timeseries.csv", irregularity),
                        (flow_source_name, flow),
                        ("energy_summary.csv", energy)):
        check_time(frame, name, args.time_step_relative_tolerance, notes)
        if not frame.empty:
            numeric = frame.select_dtypes(include=[np.number])
            if not np.all(np.isfinite(numeric.to_numpy(dtype=float))):
                invalid = True
                notes.append(f"{name}: NaN or Inf detected")

    selected_displacement = window(displacement, args.analysis_start, args.analysis_end)
    selected_flow = window(flow, args.analysis_start, args.analysis_end)
    selected_irregularity = window(irregularity, args.analysis_start, args.analysis_end)
    selected_energy = window(energy, args.analysis_start, args.analysis_end)
    selected_area_history = window(area_history, args.analysis_start, args.analysis_end)
    summary: dict[str, Any] = {
        "run_directory": str(run_dir),
        "analysis_start_s": args.analysis_start,
        "analysis_end_s": args.analysis_end,
        "tail_duration_s": args.tail_duration,
    }

    if not selected_displacement.empty:
        summary.update(statistics(selected_displacement["left_mm"].to_numpy(), "left_displacement_mm"))
        summary.update(statistics(selected_displacement["right_mm"].to_numpy(), "right_displacement_mm"))
        symmetry = selected_displacement["left_mm"].to_numpy() + selected_displacement["right_mm"].to_numpy()
        summary["left_right_symmetry_error_rms_mm"] = float(np.sqrt(np.mean(symmetry ** 2)))
        combined = 0.5 * (
            selected_displacement["right_mm"].to_numpy()
            - selected_displacement["left_mm"].to_numpy()
        )
        dominant, sigma, envelope_cv, envelope_time, envelope = dominant_and_growth(
            selected_displacement["time_s"].to_numpy(), combined,
            args.band_low_hz, args.band_high_hz, args.envelope_fit_trim_fraction,
        )
        summary["dominant_frequency_hz"] = dominant
        summary["sigma_per_s"] = sigma
        summary["zeta_effective"] = (
            -sigma / (2.0 * np.pi * dominant)
            if np.isfinite(sigma) and np.isfinite(dominant) and dominant > 0 else np.nan
        )
        summary["envelope_cv"] = envelope_cv
        tail_start = max(args.analysis_start, args.analysis_end - args.tail_duration)
        tail = selected_displacement[selected_displacement["time_s"] >= tail_start]
        tail_signal = 0.5 * (tail["right_mm"].to_numpy() - tail["left_mm"].to_numpy())
        summary["tail_oscillation_rms_mm"] = float(np.std(signal.detrend(tail_signal))) \
            if len(tail_signal) >= 3 else np.nan
    else:
        dominant = sigma = envelope_cv = np.nan
        envelope_time = envelope = np.array([])

    if not selected_flow.empty:
        area = selected_flow["min_area_mm2"].to_numpy(dtype=float)
        summary["min_area_mean_mm2"] = float(np.mean(area))
        summary["min_area_rms_mm2"] = float(np.sqrt(np.mean(area ** 2)))
        summary["min_area_fluctuation_rms_mm2"] = float(np.std(area))
        summary["min_area_min_mm2"] = float(np.min(area))
        summary["min_area_max_mm2"] = float(np.max(area))
        pg = selected_flow["current_pg_pa"].to_numpy(dtype=float)
        summary["current_pg_mean_pa"] = float(np.mean(pg))
        summary["current_pg_rms_pa"] = float(np.sqrt(np.mean(pg ** 2)))
        summary["current_pg_fluctuation_rms_pa"] = float(np.std(pg))
        tail_start = max(args.analysis_start, args.analysis_end - args.tail_duration)
        tail_flow = selected_flow[selected_flow["time_s"] >= tail_start]
        summary["tail_current_pg_mean_pa"] = (
            float(np.mean(tail_flow["current_pg_pa"])) if not tail_flow.empty else np.nan
        )
        if "input_pressure_pa" in selected_flow:
            input_pressure = float(np.mean(selected_flow["input_pressure_pa"]))
            summary["current_pg_to_input_pressure_ratio"] = (
                summary["current_pg_mean_pa"] / input_pressure if input_pressure else np.nan
            )
        flow_rate = selected_flow["flow_rate_m3_s"].to_numpy(dtype=float)
        summary["flow_rate_mean_m3_s"] = float(np.mean(flow_rate))
        summary["flow_rate_rms_m3_s"] = float(np.sqrt(np.mean(flow_rate ** 2)))
        summary["flow_rate_fluctuation_rms_m3_s"] = float(np.std(flow_rate))
        summary["closure_rate"] = float(np.mean(area <= args.closure_area_threshold_mm2))
        if "has_nonfinite" in selected_flow and selected_flow["has_nonfinite"].astype(bool).any():
            invalid = True
            notes.append("flow diagnostics reported a non-finite state")
    elif not selected_area_history.empty or not selected_irregularity.empty:
        area = (selected_area_history["min_area_mm2"].to_numpy(dtype=float)
                if not selected_area_history.empty
                else selected_irregularity["glottal_area_min_mm2"].to_numpy(dtype=float))
        summary["min_area_mean_mm2"] = float(np.mean(area))
        summary["min_area_rms_mm2"] = float(np.sqrt(np.mean(area ** 2)))
        summary["min_area_fluctuation_rms_mm2"] = float(np.std(area))
        summary["min_area_min_mm2"] = float(np.min(area))
        summary["min_area_max_mm2"] = float(np.max(area))
        summary["closure_rate"] = float(np.mean(area <= args.closure_area_threshold_mm2))
        if not selected_irregularity.empty and "flow_rate_m3_s" in selected_irregularity:
            flow_rate = selected_irregularity["flow_rate_m3_s"].to_numpy(dtype=float)
            summary["flow_rate_mean_m3_s"] = float(np.mean(flow_rate))
            summary["flow_rate_rms_m3_s"] = float(np.sqrt(np.mean(flow_rate ** 2)))
            summary["flow_rate_fluctuation_rms_m3_s"] = float(np.std(flow_rate))

    contact_source = selected_energy if "contact_flag" in selected_energy else selected_irregularity
    if not contact_source.empty and "contact_flag" in contact_source:
        summary["contact_rate"] = float(np.mean(contact_source["contact_flag"].astype(bool)))

    if not selected_energy.empty:
        for column, key in (
            ("fluid_power_total", "mean_fluid_power"),
            ("damping_power_total", "mean_damping_power"),
        ):
            if column in selected_energy:
                summary[key] = float(np.mean(selected_energy[column]))
        time_values = selected_energy["time_s"].to_numpy(dtype=float)
        if len(time_values) >= 2:
            summary["interval_fluid_work"] = float(np.trapezoid(
                selected_energy["fluid_power_total"], time_values))
            summary["interval_damping_loss"] = float(np.trapezoid(
                selected_energy["damping_power_total"], time_values))

    if not displacement.empty and not irregularity.empty and len(displacement) != len(irregularity):
        notes.append("left/right displacement and irregularity output lengths differ")

    sample_count = len(selected_displacement)
    tail_rms = float(summary.get("tail_oscillation_rms_mm", np.nan))
    if invalid:
        classification = "DIVERGED_OR_INVALID"
        reason = "non-finite or invalid numeric state was detected"
    elif sample_count < args.minimum_samples or not np.isfinite(sigma):
        classification = "INSUFFICIENT_DATA"
        reason = f"only {sample_count} usable displacement samples or no growth fit"
    elif tail_rms < args.static_rms_threshold_mm or sigma <= args.decay_threshold_per_s:
        classification = "DECAYING_TO_STATIC"
        reason = f"tail_rms={tail_rms:.6g} mm, sigma={sigma:.6g} 1/s"
    elif sigma >= args.growth_threshold_per_s:
        classification = "GROWING_OSCILLATION"
        reason = f"sigma={sigma:.6g} 1/s exceeds growth threshold"
    elif np.isfinite(envelope_cv) and envelope_cv <= args.periodic_envelope_cv_threshold:
        classification = "SUSTAINED_PERIODIC"
        reason = f"near-zero growth and envelope_cv={envelope_cv:.6g}"
    else:
        classification = "SUSTAINED_IRREGULAR"
        reason = f"near-zero growth and envelope_cv={envelope_cv:.6g}"
    assert classification in CLASSIFICATIONS
    summary["classification"] = classification
    summary["classification_reason"] = reason
    summary["warnings"] = notes

    json_path = output_path(run_dir, "analysis_summary.json")
    json_path.write_text(
        json.dumps({key: finite_or_none(value) for key, value in summary.items()},
                   indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    csv_path = output_path(run_dir, "analysis_summary.csv")
    csv_summary = {key: ("; ".join(value) if isinstance(value, list) else finite_or_none(value))
                   for key, value in summary.items()}
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_summary))
        writer.writeheader()
        writer.writerow(csv_summary)

    figure, axes = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
    if not displacement.empty:
        axes[0, 0].plot(displacement["time_s"], displacement["left_mm"], label="left")
        axes[0, 0].plot(displacement["time_s"], displacement["right_mm"], label="right")
        axes[0, 0].legend()
    axes[0, 0].set(title="Probe displacement", ylabel="mm")
    if not flow.empty:
        axes[0, 1].plot(flow["time_s"], flow["min_area_mm2"])
        axes[1, 0].plot(flow["time_s"], flow["current_pg_pa"])
        axes[1, 1].plot(flow["time_s"], flow["flow_rate_m3_s"])
    elif not area_history.empty:
        axes[0, 1].plot(area_history["time_s"], area_history["min_area_mm2"])
        if not irregularity.empty and "flow_rate_m3_s" in irregularity:
            axes[1, 1].plot(irregularity["time_s"], irregularity["flow_rate_m3_s"])
    axes[0, 1].set(title="Minimum area", ylabel="mm$^2$")
    axes[1, 0].set(title="Glottal entry pressure", ylabel="Pa")
    axes[1, 1].set(title="Flow rate", ylabel="m$^3$/s")
    if len(envelope_time):
        axes[2, 0].semilogy(envelope_time, np.maximum(envelope, np.finfo(float).tiny))
    axes[2, 0].set(title="Band-passed amplitude envelope", ylabel="mm")
    if not energy.empty:
        axes[2, 1].plot(energy["time_s"], energy["fluid_power_total"], label="fluid")
        axes[2, 1].plot(energy["time_s"], energy["damping_power_total"], label="damping")
        axes[2, 1].legend()
    axes[2, 1].set(title="Power (mass-normalized coordinate units)")
    for axis in axes.flat:
        axis.set_xlabel("time [s]")
        axis.grid(alpha=0.25)
    figure.suptitle(f"{classification}: {reason}")
    figure.savefig(output_path(run_dir, "diagnostic_overview.png"), dpi=160)
    plt.close(figure)

    for note in notes:
        warnings.warn(note)
    print(json.dumps({"classification": classification, "reason": reason,
                      "summary": str(json_path)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
