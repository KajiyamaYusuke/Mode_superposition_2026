#!/usr/bin/env python3
"""Run the minimal live/freeze-fluid modal-work exit test.

The large modal CSV files are read in chunks and only the requested modes are
retained. When the solver's internal-dt perturbation integrals are present,
they are used as the primary work estimate. Decimated integrations remain in
the output as an uncertainty and consistency check.
"""

from __future__ import annotations

import argparse
import json
import math
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from output_layout import result_path


MODAL_COLUMNS = [
    "step", "load_time_s", "state_time_s", "side", "mode_index_1based",
    "frequency_hz", "q", "qdot", "fluid_force", "total_force",
    "applied_force", "diagnostic_load_mode",
]
NATIVE_COLUMNS = [
    "step", "load_time_s", "state_time_s", "side", "mode_index_1based",
    "frequency_hz", "dq", "dv", "pfluid_pert_generalized_per_s",
    "wfluid_pert_generalized", "wdamp_pert_generalized", "event_step_excluded",
    "integration_method", "units",
]
KEY = ["step", "side", "mode_index_1based"]


@dataclass(frozen=True)
class RunInfo:
    label: str
    path: Path
    manifest: dict[str, str]
    identity: str
    load_mode: str


def arguments(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("live_run", type=Path)
    parser.add_argument("freeze_fluid_run", type=Path)
    parser.add_argument("--modes", default="1,4,5",
                        help="comma-separated one-based source-mode indices")
    parser.add_argument("--window", nargs=2, type=float, action="append",
                        metavar=("START", "END"),
                        help="analysis window; pass twice for primary and confirmation windows")
    parser.add_argument("--reference-window", nargs=2, type=float,
                        metavar=("START", "END"))
    parser.add_argument("--switch-time", type=float,
                        help="override the freeze/perturbation switch time")
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--representative-side", choices=("L", "R"), default="L")
    parser.add_argument("--chunksize", type=int, default=200_000)
    return parser.parse_args(argv)


def parse_modes(text: str) -> list[int]:
    modes = sorted({int(value.strip()) for value in text.split(",") if value.strip()})
    if not modes or any(mode < 1 for mode in modes):
        raise ValueError("--modes must contain positive one-based indices")
    return modes


def parse_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().lower()] = value.strip().strip('"')
    return values


def first_existing(run: Path, names: Iterable[str]) -> Path | None:
    for name in names:
        path = result_path(run, name)
        if path.is_file():
            return path
    return None


def run_info(label: str, run: Path) -> RunInfo:
    run = run.resolve()
    manifest_path = first_existing(run, ("manifest.txt",))
    manifest = parse_key_values(manifest_path) if manifest_path else {}
    identity = str(run)
    json_path = run / "run_manifest.json"
    if json_path.is_file():
        data = json.loads(json_path.read_text(encoding="utf-8"))
        identity = f"{data.get('case', label)}@{data.get('git_commit', 'unknown')}"
    return RunInfo(label, run, manifest, identity,
                   manifest.get("diagnostic_load_mode", "unknown").lower())


def require_columns(path: Path, required: Iterable[str]) -> None:
    columns = set(pd.read_csv(path, nrows=0).columns)
    missing = sorted(set(required) - columns)
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")


def read_selected(path: Path, columns: list[str], modes: list[int], chunksize: int,
                  max_time: float | None = None) -> pd.DataFrame:
    require_columns(path, columns)
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=columns, chunksize=chunksize):
        selected = chunk[chunk["mode_index_1based"].isin(modes)]
        if max_time is not None:
            selected = selected[selected["state_time_s"] <= np.nextafter(max_time, np.inf)]
        if not selected.empty:
            pieces.append(selected.copy())
    if not pieces:
        raise ValueError(f"{path}: no rows found for modes {modes}")
    result = pd.concat(pieces, ignore_index=True)
    result["mode_index_1based"] = result["mode_index_1based"].astype(int)
    return result.sort_values(["side", "mode_index_1based", "state_time_s"])


def final_csv_time(path: Path, time_column: str = "state_time_s") -> float:
    """Read the final CSV record without scanning a multi-gigabyte file."""
    with path.open("rb") as stream:
        header = stream.readline().decode("utf-8").rstrip("\r\n").split(",")
        try:
            time_index = header.index(time_column)
        except ValueError as error:
            raise ValueError(f"{path}: missing {time_column}") from error
        stream.seek(0, 2)
        position = stream.tell() - 1
        while position >= 0:
            stream.seek(position)
            if stream.read(1) not in {b"\n", b"\r"}:
                break
            position -= 1
        end = position + 1
        while position >= 0:
            stream.seek(position)
            if stream.read(1) == b"\n":
                break
            position -= 1
        stream.seek(position + 1)
        row = stream.read(end - position - 1).decode("utf-8").split(",")
    return float(row[time_index])


def applied_fluid_force(frame: pd.DataFrame) -> pd.Series:
    """Separate applied fluid from computed contact for live/freeze-fluid."""
    return frame["applied_force"] - (frame["total_force"] - frame["fluid_force"])


def closed_mask(values: pd.Series, start: float, end: float) -> pd.Series:
    return values.between(np.nextafter(start, -np.inf), np.nextafter(end, np.inf),
                          inclusive="both")


def manifest_float(info: RunInfo, key: str) -> float | None:
    value = info.manifest.get(key.lower())
    try:
        return float(value) if value is not None else None
    except ValueError:
        return None


def choose_times(args: argparse.Namespace, live: RunInfo, freeze: RunInfo,
                 available_end: float) -> tuple[float, tuple[float, float], list[tuple[float, float]]]:
    switch = args.switch_time
    if switch is None:
        switch = next((value for value in (
            manifest_float(freeze, "perturbation_time_sec"),
            manifest_float(live, "perturbation_time_sec"),
        ) if value is not None), None)
    if switch is None:
        raise ValueError("switch time is absent from the manifests; pass --switch-time")

    if args.reference_window:
        reference = tuple(args.reference_window)
    else:
        start = manifest_float(freeze, "diagnostic_reference_window_start_sec")
        end = manifest_float(freeze, "diagnostic_reference_window_end_sec")
        if start is None or end is None:
            raise ValueError("reference window is absent from the manifest; pass --reference-window")
        reference = (start, end)

    if args.window:
        windows = [tuple(window) for window in args.window]
    else:
        # Derive both windows relative to the recorded switch; there is no
        # absolute freeze-time default.
        guard = min(0.02, max(0.0, (available_end - switch) / 4.0))
        usable = available_end - (switch + guard)
        width = min(0.08, usable / 2.0)
        if width <= 0.0:
            raise ValueError("not enough post-switch data to construct two windows")
        first = switch + guard
        windows = [(first, first + width), (first + width, first + 2.0 * width)]

    if not reference[0] < reference[1] <= switch:
        raise ValueError("reference window must end no later than the switch time")
    if len(windows) > 2:
        raise ValueError("at most two analysis windows are supported")
    for start, end in windows:
        if not switch < start < end <= available_end + 1.0e-12:
            raise ValueError(f"invalid post-switch window {(start, end)}")
    return switch, reference, windows


def attach_native(modal: pd.DataFrame, native: pd.DataFrame) -> pd.DataFrame:
    keep = [column for column in NATIVE_COLUMNS if column not in
            {"load_time_s", "state_time_s", "frequency_hz"}]
    # Native perturbation rows start at the switch. Preserve pre-switch modal
    # rows because they define the reference window and the A/B identity check.
    return modal.merge(native[keep], on=KEY, how="left", validate="one_to_one")


def endpoint_row(group: pd.DataFrame, time: float) -> pd.Series:
    index = int(np.argmin(np.abs(group["state_time_s"].to_numpy() - time)))
    return group.iloc[index]


def trapz(values: np.ndarray, time: np.ndarray) -> float:
    if len(time) < 2:
        return math.nan
    return float(np.trapezoid(values, time))


def fit_growth(group: pd.DataFrame, start: float, end: float,
               omega: float, q_ref: float, v_ref: float) -> dict[str, float]:
    data = group[closed_mask(group["state_time_s"], start, end)]
    if len(data) < 16:
        return {"growth_rate_per_s": math.nan, "growth_fit_r2": math.nan,
                "growth_fit_log_rmse": math.nan, "growth_sample_count": len(data)}
    if {"dq", "dv"}.issubset(data.columns) and data[["dq", "dv"]].notna().all().all():
        dq = data["dq"].to_numpy(dtype=float)
        dv = data["dv"].to_numpy(dtype=float)
    else:
        dq = data["q"].to_numpy(dtype=float) - q_ref
        dv = data["qdot"].to_numpy(dtype=float) - v_ref
    amplitude = np.sqrt(dq * dq + np.square(dv / omega))
    time = data["state_time_s"].to_numpy(dtype=float)
    floor = max(float(np.nanmax(amplitude)) * 1.0e-8, np.finfo(float).tiny)
    valid = np.isfinite(amplitude) & (amplitude > floor)
    if np.count_nonzero(valid) < 8:
        return {"growth_rate_per_s": math.nan, "growth_fit_r2": math.nan,
                "growth_fit_log_rmse": math.nan,
                "growth_sample_count": int(np.count_nonzero(valid))}
    x, y = time[valid], np.log(amplitude[valid])
    coefficients = np.polyfit(x, y, 1)
    prediction = np.polyval(coefficients, x)
    residual = y - prediction
    total = float(np.sum(np.square(y - np.mean(y))))
    r2 = 1.0 - float(np.sum(np.square(residual))) / total if total > 0.0 else math.nan
    return {"growth_rate_per_s": float(coefficients[0]), "growth_fit_r2": r2,
            "growth_fit_log_rmse": float(np.sqrt(np.mean(np.square(residual)))),
            "growth_sample_count": int(len(x))}


def integrate_decimated(group: pd.DataFrame, start: float, end: float,
                        q_ref_fluid: float, damping_coefficient: float) -> dict[str, float]:
    data = group[closed_mask(group["state_time_s"], start, end)].copy()
    if len(data) < 2:
        raise ValueError(f"window [{start}, {end}] contains fewer than two samples")
    time = data["state_time_s"].to_numpy(dtype=float)
    q = data["q"].to_numpy(dtype=float)
    velocity = data["qdot"].to_numpy(dtype=float)
    force = data["Q_fluid_applied"].to_numpy(dtype=float)
    raw_power = trapz(force * velocity, time)
    raw_displacement = float(np.sum(0.5 * (force[:-1] + force[1:]) * np.diff(q)))
    static = q_ref_fluid * float(q[-1] - q[0])
    damping = trapz(damping_coefficient * velocity * velocity, time)
    sample_dt = np.diff(time)
    return {
        "sample_start_s": float(time[0]), "sample_end_s": float(time[-1]),
        "sample_count": int(len(data)), "sample_interval_median_s": float(np.median(sample_dt)),
        "q_start": float(q[0]), "q_end": float(q[-1]),
        "W_raw_decimated_trapezoid": raw_power,
        "W_raw_decimated_displacement": raw_displacement,
        "W_static": static,
        "W_dynamic_decimated_trapezoid": raw_power - static,
        "W_dynamic_decimated_displacement": raw_displacement - static,
        "D_decimated_trapezoid": damping,
    }


def integrate_native(group: pd.DataFrame, start: float, end: float,
                     q_ref_fluid: float, v_ref: float,
                     damping_coefficient: float) -> dict[str, float]:
    first, last = endpoint_row(group, start), endpoint_row(group, end)
    duration = float(last["state_time_s"] - first["state_time_s"])
    delta_q = float(last["q"] - first["q"])
    native_dynamic_centered = float(
        last["wfluid_pert_generalized"] - first["wfluid_pert_generalized"])
    native_damping_centered = float(
        last["wdamp_pert_generalized"] - first["wdamp_pert_generalized"])
    data = group[closed_mask(group["state_time_s"], float(first["state_time_s"]),
                             float(last["state_time_s"]))]
    force_delta_integral = trapz(
        data["Q_fluid_applied"].to_numpy(dtype=float) - q_ref_fluid,
        data["state_time_s"].to_numpy(dtype=float),
    )
    dynamic = native_dynamic_centered + v_ref * force_delta_integral
    damping = (native_damping_centered + 2.0 * damping_coefficient * v_ref * delta_q
               - damping_coefficient * v_ref * v_ref * duration)
    static = q_ref_fluid * delta_q
    return {
        "native_start_s": float(first["state_time_s"]),
        "native_end_s": float(last["state_time_s"]),
        "W_dynamic_internal": dynamic,
        "W_raw_internal": dynamic + static,
        "W_static_internal": static,
        "D_internal": damping,
        "native_centered_fluid_work": native_dynamic_centered,
        "native_centered_damping": native_damping_centered,
        "native_velocity_reference_correction": v_ref * force_delta_integral,
    }


def classify(w_dynamic: float, damping: float, work_error: float,
             damping_error: float, sensitivity_low: float,
             sensitivity_high: float) -> str:
    if not all(np.isfinite([w_dynamic, damping, work_error, damping_error,
                            sensitivity_low, sensitivity_high])):
        return "NEAR_ZERO_OR_INCONCLUSIVE"
    low = min(w_dynamic - work_error, sensitivity_low)
    high = max(w_dynamic + work_error, sensitivity_high)
    d_low = max(0.0, damping - damping_error)
    d_high = damping + damping_error
    scale = max(abs(w_dynamic), abs(damping), abs(low), abs(high), 1.0)
    numerical_floor = 64.0 * np.finfo(float).eps * scale
    if damping <= max(damping_error, numerical_floor):
        return "NEAR_ZERO_OR_INCONCLUSIVE"
    if high < -numerical_floor:
        return "FLUID_DAMPING"
    if low > numerical_floor and high < d_low:
        return "EXCITING_INSUFFICIENT"
    if low >= d_high:
        return "INPUT_AT_LEAST_DAMPING"
    return "NEAR_ZERO_OR_INCONCLUSIVE"


def reference_table(live: pd.DataFrame, freeze: pd.DataFrame,
                    reference: tuple[float, float], switch: float) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for key, live_group in live.groupby(["side", "mode_index_1based"], sort=True):
        freeze_group = freeze[(freeze["side"] == key[0]) &
                              (freeze["mode_index_1based"] == key[1])]
        live_ref = live_group[closed_mask(live_group["load_time_s"], *reference)]
        frozen = freeze_group[freeze_group["load_time_s"] >= switch]
        if live_ref.empty or frozen.empty:
            raise ValueError(f"missing reference/frozen samples for {key}")
        qfrozen = frozen["Q_fluid_applied"].to_numpy(dtype=float)
        row: dict[str, float | int | str] = {
            "side": key[0], "mode_index_1based": int(key[1]),
            "Q_fluid_ref": float(np.median(qfrozen)),
            "Q_fluid_ref_basis": "paired_freeze_actual_applied_value",
            "Q_fluid_ref_sensitivity_low": float(live_ref["Q_fluid_applied"].quantile(0.05)),
            "Q_fluid_ref_sensitivity_high": float(live_ref["Q_fluid_applied"].quantile(0.95)),
            "freeze_applied_fluid_range": float(np.ptp(qfrozen)),
            "q_ref": float(live_ref["q"].mean()),
            "v_ref": float(live_ref["qdot"].mean()),
        }
        if {"dq", "dv"}.issubset(live_group.columns):
            native_post = live_group[(live_group["state_time_s"] >= switch)
                                     & live_group["dq"].notna() & live_group["dv"].notna()]
            if not native_post.empty:
                row["q_ref"] = float(np.median(native_post["q"] - native_post["dq"]))
                row["v_ref"] = float(np.median(native_post["qdot"] - native_post["dv"]))
                row["state_reference_basis"] = "solver_internal_dt_reference_from_dq_dv"
        row.setdefault("state_reference_basis", "decimated_reference_window_mean")
        rows.append(row)
    return pd.DataFrame(rows)


def analyze_case(case: str, frame: pd.DataFrame, refs: pd.DataFrame,
                 info: RunInfo, windows: list[tuple[float, float]], dt: float | None) -> list[dict]:
    rows: list[dict] = []
    ref_index = refs.set_index(["side", "mode_index_1based"])
    for (side, mode), group in frame.groupby(["side", "mode_index_1based"], sort=True):
        group = group.sort_values("state_time_s")
        ref = ref_index.loc[(side, mode)]
        frequency = float(group["frequency_hz"].median())
        omega = 2.0 * np.pi * frequency
        zeta_key = "zetal" if side == "L" else "zetar"
        zeta = manifest_float(info, zeta_key)
        if zeta is None:
            raise ValueError(f"{info.path}: {zeta_key} is absent from manifest")
        damping_coefficient = 2.0 * zeta * omega
        for window_number, (start, end) in enumerate(windows, 1):
            decimated = integrate_decimated(group, start, end, float(ref.Q_fluid_ref),
                                             damping_coefficient)
            has_native = {"wfluid_pert_generalized", "wdamp_pert_generalized"}.issubset(group)
            native = (integrate_native(group, start, end, float(ref.Q_fluid_ref),
                                       float(ref.v_ref), damping_coefficient)
                      if has_native and group["wfluid_pert_generalized"].notna().any() else {})
            w_primary = native.get("W_dynamic_internal",
                                   decimated["W_dynamic_decimated_displacement"])
            d_primary = native.get("D_internal", decimated["D_decimated_trapezoid"])
            work_candidates = [decimated["W_dynamic_decimated_trapezoid"],
                               decimated["W_dynamic_decimated_displacement"]]
            work_error = max(abs(w_primary - value) for value in work_candidates)
            damping_error = abs(d_primary - decimated["D_decimated_trapezoid"])
            delta_q = decimated["q_end"] - decimated["q_start"]
            if case == "freeze_fluid":
                # Its actual applied force is the reference itself; varying an
                # unrelated candidate equilibrium load would defeat this
                # constant-force identity check.
                sensitivity = [w_primary, w_primary]
            else:
                sensitivity = sorted([
                    w_primary - (float(ref.Q_fluid_ref_sensitivity_low)
                                 - float(ref.Q_fluid_ref)) * delta_q,
                    w_primary - (float(ref.Q_fluid_ref_sensitivity_high)
                                 - float(ref.Q_fluid_ref)) * delta_q,
                ])
            classification = classify(w_primary, d_primary, work_error, damping_error,
                                      sensitivity[0], sensitivity[1])
            growth = fit_growth(group, start, end, omega, float(ref.q_ref), float(ref.v_ref))
            raw_primary = native.get("W_raw_internal",
                                     decimated["W_raw_decimated_displacement"])
            static_primary = native.get("W_static_internal", decimated["W_static"])
            closure = raw_primary - w_primary - static_primary
            flags: list[str] = []
            if not native:
                flags.append("DECIMATED_INTEGRATION_ONLY")
            if abs(float(ref.freeze_applied_fluid_range)) > max(
                    1.0e-12, 1.0e-9 * abs(float(ref.Q_fluid_ref))):
                flags.append("FREEZE_FORCE_NOT_CONSTANT")
            if classification == "NEAR_ZERO_OR_INCONCLUSIVE":
                flags.append("CLASSIFICATION_UNCERTAIN")
            if work_error >= max(abs(w_primary), np.finfo(float).tiny):
                flags.append("INTEGRATION_SENSITIVE")
            if not np.isfinite(growth["growth_fit_r2"]) or growth["growth_fit_r2"] < 0.5:
                flags.append("GROWTH_FIT_POOR")
            source_methods = (group["integration_method"].dropna().astype(str).unique()
                              if "integration_method" in group else [])
            row = {
                "case": case, "side": side, "mode_index_1based": int(mode),
                "window_index": window_number, "window_requested_start_s": start,
                "window_requested_end_s": end, "frequency_hz": frequency,
                "omega_rad_s": omega, "zeta": zeta, "c_i": damping_coefficient,
                "dt_s": dt, "run_identity": info.identity,
                "integration_method": ("solver_internal_dt_cumulative_with_reference_correction"
                                       if native else "decimated_displacement_and_trapezoid"),
                "source_internal_integration_method": (
                    ";".join(source_methods) if len(source_methods) else "not_available"),
                **ref.to_dict(), **decimated, **native,
                "W_raw": raw_primary, "W_static": static_primary,
                "W_dynamic": w_primary, "D": d_primary,
                "W_dynamic_over_D": w_primary / d_primary if d_primary > 0.0 else math.nan,
                "work_integration_error_estimate": work_error,
                "damping_integration_error_estimate": damping_error,
                "W_dynamic_reference_sensitivity_low": sensitivity[0],
                "W_dynamic_reference_sensitivity_high": sensitivity[1],
                "raw_dynamic_static_closure_error": closure,
                "classification": classification, "quality_flags": ";".join(flags) or "OK",
                **growth,
            }
            rows.append(row)
    return rows


def pre_switch_check(live: pd.DataFrame, freeze: pd.DataFrame, switch: float) -> dict[str, float]:
    merged = live[live["state_time_s"] < switch].merge(
        freeze[freeze["state_time_s"] < switch], on=KEY, suffixes=("_live", "_freeze"),
        validate="one_to_one")
    if merged.empty:
        return {"pre_switch_max_abs_q_difference": math.nan,
                "pre_switch_max_abs_qdot_difference": math.nan}
    return {
        "pre_switch_max_abs_q_difference": float(
            np.max(np.abs(merged["q_live"] - merged["q_freeze"]))),
        "pre_switch_max_abs_qdot_difference": float(
            np.max(np.abs(merged["qdot_live"] - merged["qdot_freeze"]))),
    }


def mark_window_stability(summary: pd.DataFrame) -> None:
    summary["window_stable_classification"] = True
    live = summary[summary["case"] == "live"]
    for key, group in live.groupby(["side", "mode_index_1based"]):
        stable = group["classification"].nunique() == 1
        mask = ((summary["case"] == "live") & (summary["side"] == key[0]) &
                (summary["mode_index_1based"] == key[1]))
        summary.loc[mask, "window_stable_classification"] = stable
        if not stable:
            summary.loc[mask, "quality_flags"] = summary.loc[mask, "quality_flags"].map(
                lambda value: ((value + ";") if value != "OK" else "") + "WINDOW_DEPENDENT")


def overall_conclusion(summary: pd.DataFrame) -> str:
    live = summary[(summary["case"] == "live") & (summary["window_index"] == 1)]
    stable = live[live["window_stable_classification"]]
    certain = set(stable.loc[stable["classification"] != "NEAR_ZERO_OR_INCONCLUSIVE",
                             "classification"])
    if len(certain) > 1:
        return "MIXED_MODAL_EFFECT"
    if len(certain) == 1:
        return next(iter(certain))
    return "NEAR_ZERO_OR_INCONCLUSIVE"


def response_plot(path: Path, frames: dict[str, pd.DataFrame], modes: list[int],
                  side: str, windows: list[tuple[float, float]], switch: float) -> None:
    fig, axes = plt.subplots(len(modes), 2, figsize=(11, 2.8 * len(modes)), sharex=True,
                             constrained_layout=True, squeeze=False)
    for row, mode in enumerate(modes):
        for label, frame in frames.items():
            data = frame[(frame["side"] == side) &
                         (frame["mode_index_1based"] == mode)]
            axes[row, 0].plot(data["state_time_s"], data["q"], label=label, linewidth=1.0)
            axes[row, 1].plot(data["state_time_s"], data["qdot"], label=label, linewidth=1.0)
        axes[row, 0].set_ylabel(f"mode {mode}\nq")
        axes[row, 1].set_ylabel(f"mode {mode}\nqdot")
        for axis in axes[row]:
            axis.axvline(switch, color="black", linewidth=0.7, linestyle=":")
            for start, end in windows:
                axis.axvspan(start, end, color="grey", alpha=0.08)
            axis.grid(alpha=0.2)
    axes[0, 0].legend()
    axes[0, 1].legend()
    axes[-1, 0].set_xlabel("state time [s]")
    axes[-1, 1].set_xlabel("state time [s]")
    fig.suptitle(f"Modal response: representative side {side}")
    fig.savefig(path, dpi=170)
    plt.close(fig)


def work_plot(path: Path, live: pd.DataFrame, refs: pd.DataFrame, modes: list[int],
              side: str, window: tuple[float, float], info: RunInfo) -> None:
    fig, axes = plt.subplots(len(modes), 1, figsize=(9, 2.7 * len(modes)), sharex=True,
                             constrained_layout=True, squeeze=False)
    ref_index = refs.set_index(["side", "mode_index_1based"])
    for row, mode in enumerate(modes):
        data = live[(live["side"] == side) & (live["mode_index_1based"] == mode)].copy()
        data = data[closed_mask(data["state_time_s"], *window)]
        axis = axes[row, 0]
        if data.empty:
            continue
        ref = ref_index.loc[(side, mode)]
        frequency = float(data["frequency_hz"].median())
        zeta = manifest_float(info, "zetal" if side == "L" else "zetar")
        assert zeta is not None
        coefficient = 2.0 * zeta * 2.0 * np.pi * frequency
        time = data["state_time_s"].to_numpy(dtype=float)
        if {"wfluid_pert_generalized", "wdamp_pert_generalized"}.issubset(data):
            force_integral = np.zeros(len(data))
            force_delta = data["Q_fluid_applied"].to_numpy(dtype=float) - float(ref.Q_fluid_ref)
            if len(data) > 1:
                force_integral[1:] = np.cumsum(
                    0.5 * np.diff(time) * (force_delta[:-1] + force_delta[1:]))
            w_dynamic = (data["wfluid_pert_generalized"].to_numpy(dtype=float)
                         - float(data["wfluid_pert_generalized"].iloc[0])
                         + float(ref.v_ref) * force_integral)
            delta_q = data["q"].to_numpy(dtype=float) - float(data["q"].iloc[0])
            duration = time - time[0]
            damping = (data["wdamp_pert_generalized"].to_numpy(dtype=float)
                       - float(data["wdamp_pert_generalized"].iloc[0])
                       + 2.0 * coefficient * float(ref.v_ref) * delta_q
                       - coefficient * float(ref.v_ref) ** 2 * duration)
        else:
            velocity = data["qdot"].to_numpy(dtype=float)
            force = data["Q_fluid_applied"].to_numpy(dtype=float)
            w_dynamic = np.zeros(len(data))
            damping = np.zeros(len(data))
            if len(data) > 1:
                w_dynamic[1:] = np.cumsum(0.5 * np.diff(time) * (
                    (force[:-1] - float(ref.Q_fluid_ref)) * velocity[:-1]
                    + (force[1:] - float(ref.Q_fluid_ref)) * velocity[1:]))
                damping[1:] = np.cumsum(0.5 * np.diff(time) * coefficient * (
                    velocity[:-1] ** 2 + velocity[1:] ** 2))
        axis.plot(time, w_dynamic, label="cumulative W_dynamic")
        axis.plot(time, damping, label="cumulative D", linestyle="--")
        axis.set_ylabel(f"mode {mode}\ngeneralized")
        axis.grid(alpha=0.2)
        axis.legend(loc="best")
    axes[-1, 0].set_xlabel("state time [s]")
    fig.suptitle(f"Live modal work: representative side {side}")
    fig.savefig(path, dpi=170)
    plt.close(fig)


def format_number(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return "NA" if not np.isfinite(number) else f"{number:.4e}"


def report_text(summary: pd.DataFrame, live: RunInfo, freeze: RunInfo,
                switch: float, reference: tuple[float, float],
                windows: list[tuple[float, float]], checks: dict[str, float],
                conclusion: str, command: str) -> str:
    live_primary = summary[(summary["case"] == "live") & (summary["window_index"] == 1)]
    freeze_rows = summary[summary["case"] == "freeze_fluid"]
    freeze_dynamic_max = float(freeze_rows["W_dynamic"].abs().max())
    freeze_tolerance_max = float(freeze_rows["work_integration_error_estimate"].max())
    side_pivot = summary.pivot_table(
        index=["case", "mode_index_1based", "window_index"], columns="side",
        values=["W_dynamic", "D"], aggfunc="first")
    side_differences = []
    for quantity in ("W_dynamic", "D"):
        if (quantity, "L") in side_pivot and (quantity, "R") in side_pivot:
            side_differences.append(
                float((side_pivot[(quantity, "L")] - side_pivot[(quantity, "R")]).abs().max()))
    max_side_difference = max(side_differences, default=math.nan)
    source_methods = ", ".join(sorted(summary["source_internal_integration_method"].unique()))
    lines = [
        "# Modal exit test report", "",
        f"- Overall result: **{conclusion}**",
        f"- live run: `{live.identity}` (`{live.path}`)",
        f"- freeze_fluid run: `{freeze.identity}` (`{freeze.path}`)",
        f"- switch/perturbation time: {switch:.9g} s",
        f"- reference window: [{reference[0]:.9g}, {reference[1]:.9g}] s",
        f"- analysis windows: " + ", ".join(
            f"[{start:.9g}, {end:.9g}] s" for start, end in windows),
        "- units: generalized mass-normalized units (not asserted as J/W)", "",
        "## Method", "",
        "The solver's existing `perturbation_energy.csv` internal-dt cumulative integrals are "
        "used as the primary estimate. Its centered velocity is corrected back to qdot; "
        "decimated `modal_force_split.csv` displacement and power quadratures provide the "
        "integration-error estimate. Applied fluid force is separated as "
        "`applied_force - (total_force - fluid_force)`. The paired freeze run's actual "
        "constant applied value is Q_fluid_ref; the 5th--95th percentile of the live "
        "reference window is retained as a sensitivity range, not asserted to be a unique "
        "physical equilibrium load.", "",
        f"The reused solver integration is recorded as `{source_methods}`. It evaluates load "
        "at t and state velocity at t+dt; therefore its difference from the decimated "
        "displacement/power quadratures is retained as uncertainty rather than hidden.", "",
        "## Primary live classifications", "",
        "| Side | Mode | W_dynamic | D | W_dynamic/D | Growth [1/s] | Classification | Stable |",
        "|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for _, row in live_primary.iterrows():
        lines.append(
            f"| {row.side} | {int(row.mode_index_1based)} | {format_number(row.W_dynamic)} | "
            f"{format_number(row.D)} | {format_number(row.W_dynamic_over_D)} | "
            f"{format_number(row.growth_rate_per_s)} | {row.classification} | "
            f"{bool(row.window_stable_classification)} |")
    lines.extend(["", "## Validation", "",
                  f"- Pre-switch max |live-freeze| in q: "
                  f"{format_number(checks['pre_switch_max_abs_q_difference'])}",
                  f"- Pre-switch max |live-freeze| in qdot: "
                  f"{format_number(checks['pre_switch_max_abs_qdot_difference'])}",
                  f"- Maximum freeze applied-fluid range: "
                  f"{format_number(summary['freeze_applied_fluid_range'].max())}",
                  f"- Maximum |freeze_fluid W_dynamic|: {format_number(freeze_dynamic_max)} "
                  f"(maximum integration estimate {format_number(freeze_tolerance_max)}; "
                  f"{'PASS' if freeze_dynamic_max <= freeze_tolerance_max else 'FAIL'})",
                  f"- Maximum left-right difference in W_dynamic or D: "
                  f"{format_number(max_side_difference)}",
                  f"- Minimum D: {format_number(summary['D'].min())}",
                  f"- Maximum |raw - dynamic - static|: "
                  f"{format_number(summary['raw_dynamic_static_closure_error'].abs().max())}",
                  "", "## Interpretation", ""])
    if conclusion == "MIXED_MODAL_EFFECT":
        lines.append("The resolved modes do not share one fluid-effect classification. This is "
                     "evidence of mode-dependent fluid action, not proof of direct intermodal "
                     "energy transfer.")
    elif conclusion == "NEAR_ZERO_OR_INCONCLUSIVE":
        lines.append("The available windows, reference sensitivity, and integration estimates do "
                     "not support a robust overall sign classification. A primary-window result "
                     "may still be informative when explicitly marked unstable. No broad rerun or "
                     "coefficient sweep is implied; if a strict sign is required, the next minimal "
                     "check is a two-case selected-mode accumulator using "
                     "`Q_n * (q_(n+1) - q_n)`.")
    else:
        lines.append(f"The resolved primary-window modes support `{conclusion}` where their "
                     "quality flags and confirmation-window classifications remain stable.")
    lines.extend(["", "## Reproduction", "", "```bash", command, "```", ""])
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = arguments(argv)
    modes = parse_modes(args.modes)
    if args.chunksize < 1:
        raise ValueError("--chunksize must be positive")
    live_info = run_info("live", args.live_run)
    freeze_info = run_info("freeze_fluid", args.freeze_fluid_run)
    if live_info.load_mode not in {"live", "unknown"}:
        raise ValueError(f"live run reports diagnostic_load_mode={live_info.load_mode}")
    if freeze_info.load_mode not in {"freeze_fluid", "unknown"}:
        raise ValueError("comparison run is not freeze_fluid")

    modal_paths = {
        "live": result_path(live_info.path, "modal_force_split.csv"),
        "freeze_fluid": result_path(freeze_info.path, "modal_force_split.csv"),
    }
    for path in modal_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    available_end = min(final_csv_time(path) for path in modal_paths.values())
    switch, reference, windows = choose_times(args, live_info, freeze_info, available_end)
    max_time = max(end for _, end in windows)

    modal = {label: read_selected(path, MODAL_COLUMNS, modes, args.chunksize, max_time)
             for label, path in modal_paths.items()}
    for frame in modal.values():
        frame["Q_fluid_applied"] = applied_fluid_force(frame)

    native_available = True
    for label, info in (("live", live_info), ("freeze_fluid", freeze_info)):
        path = result_path(info.path, "perturbation_energy.csv")
        if not path.is_file() or path.stat().st_size == 0:
            native_available = False
            break
        native = read_selected(path, NATIVE_COLUMNS, modes, args.chunksize, max_time)
        modal[label] = attach_native(modal[label], native)

    refs = reference_table(modal["live"], modal["freeze_fluid"], reference, switch)
    checks = pre_switch_check(modal["live"], modal["freeze_fluid"], switch)
    rows: list[dict] = []
    for label, info in (("live", live_info), ("freeze_fluid", freeze_info)):
        rows.extend(analyze_case(label, modal[label], refs, info, windows,
                                 manifest_float(info, "dt_s")))
    summary = pd.DataFrame(rows)
    mark_window_stability(summary)
    conclusion = overall_conclusion(summary)

    output = (args.output_directory.resolve() if args.output_directory else
              live_info.path.parent / "modal_exit_test")
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output / "modal_exit_summary.csv", index=False)
    response_plot(output / "modal_exit_response.png", modal, modes,
                  args.representative_side, windows, switch)
    work_plot(output / "modal_exit_work.png", modal["live"], refs, modes,
              args.representative_side, windows[0], live_info)

    command_parts = [sys.executable, str(Path(__file__).resolve()), str(live_info.path),
                     str(freeze_info.path), "--modes", ",".join(map(str, modes)),
                     "--output-directory", str(output)]
    for start, end in windows:
        command_parts.extend(["--window", f"{start:.17g}", f"{end:.17g}"])
    command_parts.extend(["--reference-window", f"{reference[0]:.17g}",
                          f"{reference[1]:.17g}", "--switch-time", f"{switch:.17g}"])
    command = " ".join(shlex.quote(part) for part in command_parts)
    (output / "modal_exit_report.md").write_text(
        report_text(summary, live_info, freeze_info, switch, reference, windows,
                    checks, conclusion, command), encoding="utf-8")
    metadata = {
        "native_internal_dt_used": native_available,
        "switch_time_s": switch, "reference_window_s": reference,
        "analysis_windows_s": windows, "modes": modes, "checks": checks,
        "overall_conclusion": conclusion,
    }
    (output / "modal_exit_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote modal exit test to {output}")
    print(f"Overall result: {conclusion}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
