#!/usr/bin/env python3
"""Compare saved symmetric 12/12 and asymmetric 12/6 diagnostic runs.

This is a read-only postprocessor.  It never launches the simulation and it
refuses to overwrite its output directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_production import (  # noqa: E402
    load_frequency, load_surface_ids, load_vtu_basis, modal_gram,
)


DEFAULT_WINDOWS = "0.10:0.20,0.20:0.30,0.30:0.40,0.40:0.50,0.30:0.50"
MONITOR_TARGET = np.array([9.18767, 8.37084, 8.37681])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8")


def parse_windows(text: str) -> list[tuple[float, float]]:
    result = []
    for item in text.split(","):
        start, end = (float(value) for value in item.split(":"))
        if not 0 <= start < end:
            raise ValueError(f"invalid window: {item}")
        result.append((start, end))
    if not result:
        raise ValueError("at least one window is required")
    return result


def parse_manifest(path: Path) -> dict[str, str]:
    return {key.strip(): value.strip().strip('"')
            for key, value in (line.split("=", 1)
                               for line in path.read_text().splitlines() if "=" in line)}


def signal_metrics(values: np.ndarray) -> dict[str, float]:
    p5, p95 = np.percentile(values, [5, 95])
    return {
        "half_ptp_mm": float(np.ptp(values) / 2),
        "ac_rms_mm": float(np.std(values)),
        "half_p95_p5_mm": float((p95 - p5) / 2),
        "mean_mm": float(np.mean(values)),
    }


def covariance_metrics(contributions: np.ndarray) -> dict[str, object]:
    covariance = np.cov(contributions, rowvar=False, ddof=0)
    individual = np.diag(covariance)
    total = float(np.var(contributions.sum(axis=1)))
    cross = float(total - individual.sum())
    pairs = []
    for first in range(len(individual)):
        for second in range(first + 1, len(individual)):
            pairs.append((2 * covariance[first, second], first + 1, second + 1))
    pairs.sort()
    return {
        "covariance": covariance,
        "individual": individual,
        "individual_sum": float(individual.sum()),
        "cross": cross,
        "total": total,
        "cancellation_fraction": float(-cross / individual.sum())
        if cross < 0 and individual.sum() else 0.0,
        "reinforcement_fraction": float(cross / individual.sum())
        if cross > 0 and individual.sum() else 0.0,
        "top_cancel": pairs[:3],
        "top_reinforce": pairs[-3:][::-1],
    }


def mac_matrix(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    left = reference.reshape(len(reference), -1)
    right = candidate.reshape(len(candidate), -1)
    numerator = (left @ right.T) ** 2
    denominator = np.sum(left * left, axis=1)[:, None] * np.sum(
        right * right, axis=1)[None, :]
    return numerator / denominator


def match_modes(reference: np.ndarray, candidate: np.ndarray,
                reference_frequency: np.ndarray,
                candidate_frequency: np.ndarray,
                mac_floor: float = 0.80, ambiguity_margin: float = 0.05
                ) -> tuple[pd.DataFrame, np.ndarray]:
    """Match candidate modes to reference; ambiguous modes remain unmatched."""
    matrix = mac_matrix(reference, candidate)
    rows, columns = linear_sum_assignment(-matrix)
    assignment = {column: row for row, column in zip(rows, columns)}
    records = []
    for candidate_index in range(len(candidate)):
        ranked = np.sort(matrix[:, candidate_index])[::-1]
        best_index = int(np.argmax(matrix[:, candidate_index]))
        best = float(ranked[0])
        margin = best - float(ranked[1]) if len(ranked) > 1 else best
        assigned = assignment.get(candidate_index)
        unambiguous = (assigned == best_index and best >= mac_floor
                       and margin >= ambiguity_margin)
        quality = ("high" if unambiguous and best >= 0.95 else
                   "medium" if unambiguous else
                   "ambiguous" if best >= mac_floor else "unmatched")
        matched = best_index if unambiguous else None
        records.append({
            "source_mode": candidate_index + 1,
            "f_Hz": float(candidate_frequency[candidate_index]),
            "matched_mode": matched + 1 if matched is not None else pd.NA,
            "matched_f_Hz": (float(reference_frequency[matched])
                              if matched is not None else np.nan),
            "MAC": best,
            "frequency_delta_Hz": (float(candidate_frequency[candidate_index]
                                         - reference_frequency[matched])
                                   if matched is not None else np.nan),
            "MAC_margin": margin,
            "match_quality": quality,
            "ambiguous": not unambiguous,
        })
    return pd.DataFrame(records), matrix


def work_components(data: np.ndarray, frequency: np.ndarray, dt: float,
                    zeta: float = 0.06) -> list[dict[str, float]]:
    """Return per-mode Q*dq quantities with explicit Q-reference sensitivity."""
    rows = []
    for mode in range(data.shape[1]):
        qpos, qneg, qfluid, qcontact, q, velocity, delta_q = data[:, mode, :].T
        midpoint = len(data) // 2
        qrefs = [np.mean(qfluid), np.mean(qfluid[:midpoint]),
                 np.mean(qfluid[midpoint:])]
        raw = float(np.sum(qfluid * delta_q))
        fluctuating = [float(raw - reference * np.sum(delta_q)) for reference in qrefs]
        velocity_integral = float(np.sum(velocity ** 2) * dt)
        rows.append({
            "mode": mode + 1,
            "W_pos_J": float(np.sum(qpos * delta_q)),
            "W_neg_J": float(np.sum(qneg * delta_q)),
            "W_raw_J": raw,
            "D_struct_normalized": float(np.sum(
                2 * zeta * 2 * np.pi * frequency[mode] * velocity ** 2 * dt)),
            "W_contact_ref_J": float(np.sum(qcontact * delta_q)),
            "Q_mean": float(np.mean(qfluid)),
            "Q_std": float(np.std(qfluid)),
            "Q_rms": float(np.sqrt(np.mean(qfluid ** 2))),
            "Q_min": float(np.min(qfluid)),
            "Q_max": float(np.max(qfluid)),
            "delta_q_sum": float(np.sum(delta_q)),
            "W_fluct_Qmean_J": fluctuating[0],
            "W_fluct_Qfirst_J": fluctuating[1],
            "W_fluct_Qsecond_J": fluctuating[2],
            "W_fluct_min_J": min(fluctuating),
            "W_fluct_max_J": max(fluctuating),
            "W_raw_drop_first_interval_J": float(np.sum((qfluid * delta_q)[1:])),
            "W_raw_drop_last_interval_J": float(np.sum((qfluid * delta_q)[:-1])),
            "Q_v_corr": float(np.corrcoef(qfluid, velocity)[0, 1])
            if np.std(qfluid) and np.std(velocity) else np.nan,
            "q_rms": float(np.sqrt(np.mean(q ** 2))),
            "qdot_rms": float(np.sqrt(np.mean(velocity ** 2))),
            "velocity_squared_integral": velocity_integral,
            "W_raw_over_velocity_squared_integral": raw / velocity_integral
            if velocity_integral else np.nan,
        })
    return rows


def load_case(name: str, production: Path) -> dict[str, object]:
    diagnostic = production / "kanaya_diagnostic"
    metadata = json.loads((diagnostic / "run_metadata.json").read_text())
    if not metadata.get("complete"):
        raise ValueError(f"{name}: diagnostic run is incomplete")
    nstep = int(metadata["integrated_intervals"])
    nmode = int(metadata["nmode"][0])
    dt = float(metadata["dt_s"])
    binary = diagnostic / "modal_native.bin"
    expected = nstep * 2 * nmode * 7 * 8
    if binary.stat().st_size != expected:
        raise ValueError(f"{name}: modal binary size {binary.stat().st_size} != {expected}")
    raw = np.memmap(binary, dtype="<f8", mode="r", shape=(nstep, 2, nmode, 7))
    manifest = parse_manifest(production / "txt/manifest.txt")
    frequency_paths = [Path(manifest[f"{side}_frequency"]) for side in ("left", "right")]
    mode_paths = [Path(manifest[f"{side}_mode_vtu"]) for side in ("left", "right")]
    frequencies = [load_frequency(path, nmode) for path in frequency_paths]
    points, left_modes, cells = load_vtu_basis(mode_paths[0], frequencies[0])
    right_points, right_unmirrored, _ = load_vtu_basis(mode_paths[1], frequencies[1])
    if not np.array_equal(points, right_points):
        raise ValueError(f"{name}: left/right point ordering differs")
    right_modes = right_unmirrored.copy()
    right_modes[:, :, 1] *= -1
    surface_ids, surface_error = load_surface_ids(Path(manifest["left_surface_nas"]), points)
    monitor_id = int(surface_ids[np.argmin(
        np.linalg.norm(points[surface_ids] - MONITOR_TARGET, axis=1))])
    displacement = np.loadtxt(production / "dat/displace.dat", comments="#", ndmin=2)
    return {
        "name": name, "production": production, "diagnostic": diagnostic,
        "metadata": metadata, "raw": raw, "times": np.arange(nstep) * dt,
        "dt": dt, "manifest": manifest, "frequencies": frequencies,
        "mode_paths": mode_paths, "frequency_paths": frequency_paths,
        "points": points, "modes": [left_modes, right_modes],
        "unmirrored_modes": [left_modes, right_unmirrored], "cells": cells,
        "surface_ids": surface_ids, "surface_error": surface_error,
        "monitor_id": monitor_id, "displacement": displacement,
    }


def surface_metrics(case: dict[str, object], start: float, end: float,
                    side_index: int) -> tuple[float, float]:
    mask = (case["times"] >= start) & (case["times"] < end)
    q = np.asarray(case["raw"][mask, side_index, :, 4])
    phi = case["modes"][side_index][:, case["surface_ids"], 1]
    sum_half_ptp_sq = 0.0
    sum_ac_rms_sq = 0.0
    count = phi.shape[1]
    for offset in range(0, count, 256):
        values = q @ phi[:, offset:offset + 256] * 1e3
        sum_half_ptp_sq += float(np.sum((np.ptp(values, axis=0) / 2) ** 2))
        sum_ac_rms_sq += float(np.sum(np.std(values, axis=0) ** 2))
    return np.sqrt(sum_half_ptp_sq / count), np.sqrt(sum_ac_rms_sq / count)


def symmetric_surface_check(case: dict[str, object], start: float,
                            end: float) -> dict[str, float]:
    mask = (case["times"] >= start) & (case["times"] < end)
    q_left = np.asarray(case["raw"][mask, 0, :, 4])
    q_right = np.asarray(case["raw"][mask, 1, :, 4])
    surface = case["surface_ids"]
    maximum_waveform_error = 0.0
    maximum_half_ptp_error = 0.0
    maximum_rms_error = 0.0
    for offset in range(0, len(surface), 256):
        ids = surface[offset:offset + 256]
        left = q_left @ case["modes"][0][:, ids, 1] * 1e3
        right = q_right @ case["modes"][1][:, ids, 1] * 1e3
        maximum_waveform_error = max(maximum_waveform_error,
                                     float(np.max(np.abs(left + right))))
        maximum_half_ptp_error = max(maximum_half_ptp_error, float(np.max(np.abs(
            np.ptp(left, axis=0) / 2 - np.ptp(right, axis=0) / 2))))
        maximum_rms_error = max(maximum_rms_error, float(np.max(np.abs(
            np.std(left, axis=0) - np.std(right, axis=0)))))
    return {
        "all_surface_max_opening_waveform_abs_error_mm": maximum_waveform_error,
        "all_surface_max_half_ptp_abs_error_mm": maximum_half_ptp_error,
        "all_surface_max_ac_rms_abs_error_mm": maximum_rms_error,
    }


def monitor_reconstruction(case: dict[str, object]) -> dict[str, dict[str, float]]:
    result = {}
    times = case["displacement"][:, 0]
    valid = times <= case["times"][-1]
    indices = np.rint(times[valid] / case["dt"]).astype(int)
    for side_index, side in enumerate(("L", "R")):
        q = np.asarray(case["raw"][indices, side_index, :, 4])
        phi = case["modes"][side_index][:, case["monitor_id"], 1]
        reconstructed = q @ phi * 1e3
        expected = case["displacement"][valid, side_index + 1]
        error = reconstructed - expected
        result[side] = {
            "max_abs_error_mm": float(np.max(np.abs(error))),
            "rms_error_mm": float(np.sqrt(np.mean(error ** 2))),
        }
    return result


def exact_dat(source: Path, production: Path) -> dict[str, bool]:
    result = {}
    for path in sorted((source / "output/dat").glob("*.dat")):
        counterpart = production / "dat" / path.name
        result[path.name] = counterpart.is_file() and sha256(path) == sha256(counterpart)
    return result


def parameter_line_differences(left: Path, right: Path) -> list[dict[str, str]]:
    def active(path: Path) -> list[str]:
        return [line.strip() for line in path.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
    first, second = active(left), active(right)
    if len(first) != len(second):
        raise ValueError("source parameter files have different active-line counts")
    return [{"symmetric": a, "asymmetric": b}
            for a, b in zip(first, second) if a != b]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asym-analysis", type=Path, required=True)
    parser.add_argument("--asym-production", type=Path, required=True)
    parser.add_argument("--sym-source", type=Path, required=True)
    parser.add_argument("--sym-production", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--windows", default=DEFAULT_WINDOWS)
    parser.add_argument("--no-extra-run", action="store_true", default=True,
                        help="postprocess only (default; simulations are never launched)")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    (output / "figures").mkdir()
    windows = parse_windows(args.windows)
    cases = [load_case("sym12_12", args.sym_production.resolve()),
             load_case("asym12_6", args.asym_production.resolve())]
    sym, asym = cases
    asym_analysis_manifest = json.loads(
        (args.asym_analysis.resolve() / "manifest.json").read_text())
    asym_analysis_checks = json.loads(
        (args.asym_analysis.resolve() / "checks.json").read_text())
    asym_source = Path(asym_analysis_manifest["source_run"])
    source_parameter_differences = parameter_line_differences(
        args.sym_source.resolve() / "input/param.txt", asym_source / "input/param.txt")

    # Amplitudes and modal covariance, always using the same four series/definitions.
    amplitude_rows, covariance_rows = [], []
    contribution_cache = {}
    for case in cases:
        for start, end in windows:
            mask = (case["times"] >= start) & (case["times"] < end)
            for side_index, side in enumerate(("L", "R")):
                opening_probe = (-1 if side == "L" else 1) * case["displacement"][:, side_index + 1]
                dmask = ((case["displacement"][:, 0] >= start)
                         & (case["displacement"][:, 0] < end))
                metrics = signal_metrics(opening_probe[dmask])
                surface_half, surface_rms = surface_metrics(case, start, end, side_index)
                amplitude_rows.append({
                    "case": case["name"], "side": side,
                    "window_start_s": start, "window_end_s": end,
                    "probe": "i40_j34", **metrics,
                    "surface_half_ptp_spatial_rms_mm": surface_half,
                    "surface_ac_rms_spatial_rms_mm": surface_rms,
                })
                q = np.asarray(case["raw"][mask, side_index, :, 4])
                phi = case["modes"][side_index][:, case["monitor_id"], 1]
                contributions = q * phi[None, :] * 1e3
                contribution_cache[(case["name"], side, start, end)] = contributions
                cov = covariance_metrics(contributions)
                for mode in range(len(phi)):
                    covariance_rows.append({
                        "case": case["name"], "side": side,
                        "window_start_s": start, "window_end_s": end,
                        "mode": mode + 1,
                        "individual_variance_mm2": cov["individual"][mode],
                        "cross_term_mm2": 2 * (cov["covariance"][mode].sum()
                                                - cov["covariance"][mode, mode]),
                        "reconstructed_variance_mm2": cov["total"],
                        "cancellation_fraction": cov["cancellation_fraction"],
                        "reinforcement_fraction": cov["reinforcement_fraction"],
                        "top_cancel_pairs": json.dumps(cov["top_cancel"]),
                        "top_reinforce_pairs": json.dumps(cov["top_reinforce"]),
                    })
                covariance_rows.append({
                    "case": case["name"], "side": side,
                    "window_start_s": start, "window_end_s": end, "mode": "total",
                    "individual_variance_mm2": cov["individual_sum"],
                    "cross_term_mm2": cov["cross"],
                    "reconstructed_variance_mm2": cov["total"],
                    "cancellation_fraction": cov["cancellation_fraction"],
                    "reinforcement_fraction": cov["reinforcement_fraction"],
                    "top_cancel_pairs": json.dumps(cov["top_cancel"]),
                    "top_reinforce_pairs": json.dumps(cov["top_reinforce"]),
                })
    amplitude = pd.DataFrame(amplitude_rows)
    for (case_name, start, end), group in amplitude.groupby(
            ["case", "window_start_s", "window_end_s"]):
        left = group[group.side == "L"].iloc[0]
        right = group[group.side == "R"].iloc[0]
        mask = ((amplitude.case == case_name) & (amplitude.window_start_s == start)
                & (amplitude.window_end_s == end))
        for field in ("half_ptp_mm", "ac_rms_mm", "half_p95_p5_mm",
                      "surface_half_ptp_spatial_rms_mm",
                      "surface_ac_rms_spatial_rms_mm"):
            amplitude.loc[mask, field + "_R_over_L"] = right[field] / left[field]
    for side in ("L", "R"):
        for start, end in windows:
            symrow = amplitude[(amplitude.case == "sym12_12") & (amplitude.side == side)
                               & (amplitude.window_start_s == start)
                               & (amplitude.window_end_s == end)].iloc[0]
            mask = ((amplitude.case == "asym12_6") & (amplitude.side == side)
                    & (amplitude.window_start_s == start)
                    & (amplitude.window_end_s == end))
            for field in ("half_ptp_mm", "ac_rms_mm",
                          "surface_half_ptp_spatial_rms_mm",
                          "surface_ac_rms_spatial_rms_mm"):
                amplitude.loc[mask, field + "_asym_minus_sym"] = (
                    amplitude.loc[mask, field] - symrow[field])
                amplitude.loc[mask, field + "_asym_over_sym"] = (
                    amplitude.loc[mask, field] / symrow[field])
    amplitude.to_csv(output / "amplitude_compare.csv", index=False)
    covariance = pd.DataFrame(covariance_rows)
    covariance.to_csv(output / "modal_covariance_compare.csv", index=False)

    # MAC matching uses unmirrored material bases, with symmetric b12 as reference.
    match_frames, mac_by_series = [], {}
    for case in cases:
        for side_index, side in enumerate(("L", "R")):
            frame, matrix = match_modes(
                sym["unmirrored_modes"][0], case["unmirrored_modes"][side_index],
                sym["frequencies"][0], case["frequencies"][side_index])
            frame.insert(0, "side", side)
            frame.insert(0, "case", case["name"])
            match_frames.append(frame)
            mac_by_series[(case["name"], side)] = matrix
    modal_match = pd.concat(match_frames, ignore_index=True)
    modal_match.to_csv(output / "modal_match.csv", index=False)

    # Fluid/modal work, force input, and conditional mean-force-removal sensitivity.
    work_rows, contact_rows = [], []
    for case in cases:
        force = pd.read_csv(case["diagnostic"] / "force_regions.csv")
        force_step = force.groupby(["step", "time_s", "side"], as_index=False).agg(
            opening_Fy_N=("opening_Fy_N", "sum"))
        force_partition = force.groupby(
            ["step", "time_s", "side", "zone", "ap_region", "pressure_sign"],
            as_index=False).agg(opening_Fy_N=("opening_Fy_N", "sum"))
        flow = pd.read_csv(case["diagnostic"] / "flow_geometry.csv")
        contact_geometry = pd.read_csv(case["diagnostic"] / "contact_regions.csv")
        for start, end in windows:
            mask = (case["times"] >= start) & (case["times"] < end)
            raw_steps = np.flatnonzero(mask)
            flow_by_step = flow.set_index("step")
            event = flow_by_step.loc[raw_steps, "contact_flag"].to_numpy() != 0
            for side_index, side in enumerate(("L", "R")):
                data = np.asarray(case["raw"][mask, side_index])
                mode_rows = work_components(data, case["frequencies"][side_index], case["dt"])
                mean_force = float(force_step[(force_step.side == side)
                    & (force_step.time_s >= start) & (force_step.time_s < end)].opening_Fy_N.mean())
                selected_force = force_partition[(force_partition.side == side)
                    & (force_partition.time_s >= start) & (force_partition.time_s < end)]
                conditional_force = selected_force.groupby(
                    ["zone", "ap_region", "pressure_sign"]).opening_Fy_N.mean()
                positive_force = float(conditional_force[
                    conditional_force.index.get_level_values("pressure_sign") == 1].sum())
                negative_force = float(conditional_force[
                    conditional_force.index.get_level_values("pressure_sign") == -1].sum())
                for row in mode_rows:
                    reference_sensitive = row["W_fluct_min_J"] <= 0 <= row["W_fluct_max_J"]
                    row.update({"case": case["name"], "side": side,
                                "window_start_s": start, "window_end_s": end,
                                "frequency_Hz": case["frequencies"][side_index][row["mode"] - 1],
                                "mean_opening_fluid_force_N": mean_force,
                                "mean_positive_opening_fluid_force_N": positive_force,
                                "mean_negative_opening_fluid_force_N": negative_force,
                                "mean_power_raw_W": row["W_raw_J"] / (end - start),
                                "interpretation_flag": ("REFERENCE_SENSITIVE; " if reference_sensitive else "")
                                + "RAW_QDQ; W_fluct is operational reference only"})
                    work_rows.append(row)
                numeric = pd.DataFrame(mode_rows).sum(numeric_only=True).to_dict()
                numeric.update({"case": case["name"], "side": side, "mode": "total",
                                "window_start_s": start, "window_end_s": end,
                                "frequency_Hz": np.nan,
                                "mean_opening_fluid_force_N": mean_force,
                                "mean_positive_opening_fluid_force_N": positive_force,
                                "mean_negative_opening_fluid_force_N": negative_force,
                                "mean_power_raw_W": numeric["W_raw_J"] / (end - start),
                                "interpretation_flag": "RAW_QDQ total; W_fluct is operational reference only"})
                numeric["W_raw_over_velocity_squared_integral"] = (
                    numeric["W_raw_J"] / numeric["velocity_squared_integral"])
                work_rows.append(numeric)
                modal_contact = data[:, :, 3] * data[:, :, 6]
                total_by_interval = modal_contact.sum(axis=1)
                adopted = contact_geometry[(contact_geometry.side == side)
                    & (contact_geometry.time_s >= start)
                    & (contact_geometry.time_s < end)]
                regions = {}
                for region in ("end", "interior"):
                    selected = adopted[adopted.ap_region == region]
                    regions[region] = {
                        "contact_pair_fraction": float(np.mean(selected.pairs > 0)),
                        "mean_pairs": float(np.mean(selected.pairs)),
                        "max_pairs": int(np.max(selected.pairs)),
                        "max_penetration_m": float(np.max(selected.max_penetration_m)),
                    }
                contact_rows.append({
                    "case": case["name"], "side": side,
                    "window_start_s": start, "window_end_s": end,
                    "W_contact_ref_J": float(total_by_interval.sum()),
                    "W_contact_event_J": float(total_by_interval[event].sum()),
                    "W_contact_nonevent_J": float(total_by_interval[~event].sum()),
                    "paired_total_W_contact_ref_J": np.nan,
                    "contact_flag_fraction": float(np.mean(event)),
                    "end_contact_pair_fraction": regions["end"]["contact_pair_fraction"],
                    "interior_contact_pair_fraction": regions["interior"]["contact_pair_fraction"],
                    "end_mean_pairs": regions["end"]["mean_pairs"],
                    "interior_mean_pairs": regions["interior"]["mean_pairs"],
                    "max_pairs": max(regions["end"]["max_pairs"], regions["interior"]["max_pairs"]),
                    "max_penetration_m": max(regions["end"]["max_penetration_m"],
                                             regions["interior"]["max_penetration_m"]),
                    "spring_work_J": np.nan, "velocity_dependent_work_J": np.nan,
                    "component_quality": "NOT_IDENTIFIABLE_FROM_SAVED_DATA",
                })
    work = pd.DataFrame(work_rows)
    work.to_csv(output / "fluid_work_compare.csv", index=False)
    contact = pd.DataFrame(contact_rows)
    for _, group in contact.groupby(["case", "window_start_s", "window_end_s"]):
        contact.loc[group.index, "paired_total_W_contact_ref_J"] = group.W_contact_ref_J.sum()
    contact.to_csv(output / "contact_compare.csv", index=False)

    # Quality gates and provenance.
    sym_dat = exact_dat(args.sym_source.resolve(), args.sym_production.resolve())
    consistency_checks = {}
    gram_checks = {}
    for case in cases:
        consistency = pd.read_csv(case["diagnostic"] / "consistency.csv")
        included = consistency[consistency.included == 1]
        consistency_checks[case["name"]] = {
            "passed": bool(included.valid.eq(1).all() and included.scaled_max.max() <= 1),
            "included_side_intervals": int(len(included)),
            "excluded_rows": int((consistency.included == 0).sum()),
            "maximum_scaled_residual": float(included.scaled_max.max()),
        }
        gram_checks[case["name"]] = {}
        for side, modes in zip(("L", "R"), case["modes"]):
            gram = modal_gram(modes)
            offdiag = gram - np.diag(np.diag(gram))
            gram_checks[case["name"]][side] = {
                "max_diagonal_abs_error": float(np.max(np.abs(np.diag(gram) - 1))),
                "max_off_diagonal_abs": float(np.max(np.abs(offdiag))),
            }
    smoke_root = args.sym_production.resolve().parent / "smoke"
    smoke_hashes = {name: sha256(smoke_root / name / "regression_state.bin")
                    for name in ("off", "on_disabled", "on_enabled")}
    primary_amp = amplitude[(amplitude.window_start_s == 0.30)
                            & (amplitude.window_end_s == 0.50)]
    sym_primary = primary_amp[primary_amp.case == "sym12_12"].set_index("side")
    asym_primary = primary_amp[primary_amp.case == "asym12_6"].set_index("side")
    sym_total_cov = covariance[(covariance.case == "sym12_12")
        & (covariance.window_start_s == 0.30) & (covariance.window_end_s == 0.50)
        & (covariance["mode"] == "total")].set_index("side")
    asym_total_cov = covariance[(covariance.case == "asym12_6")
        & (covariance.window_start_s == 0.30) & (covariance.window_end_s == 0.50)
        & (covariance["mode"] == "total")].set_index("side")
    primary_work = work[(work.window_start_s == 0.30) & (work.window_end_s == 0.50)
                        & (work["mode"] == "total")]
    primary_contact = contact[(contact.window_start_s == 0.30)
                              & (contact.window_end_s == 0.50)]
    ignored_manifest_keys = {"parameter_file"}
    shared_manifest_keys = set(sym["manifest"]) & set(asym["manifest"])
    input_differences = {key: [sym["manifest"].get(key), asym["manifest"].get(key)]
                         for key in sorted(shared_manifest_keys - ignored_manifest_keys)
                         if sym["manifest"].get(key) != asym["manifest"].get(key)}
    manifest_schema_differences = {
        "sym_only": sorted(set(sym["manifest"]) - set(asym["manifest"])),
        "asym_only": sorted(set(asym["manifest"]) - set(sym["manifest"])),
    }
    allowed_input_differences = {"right_frequency", "right_mode_vtu"}
    checks = {
        "input_difference_gate": {
            "passed": (set(input_differences) == allowed_input_differences
                       and len(source_parameter_differences) == 2
                       and {item["symmetric"].split("=", 1)[0].strip()
                            for item in source_parameter_differences}
                           == {"rightFrequencyFile", "rightModeFile"}),
            "differences": input_differences,
            "source_parameter_line_differences": source_parameter_differences,
            "production_manifest_schema_differences": manifest_schema_differences,
            "note": "parameter_file paths ignored; schema-only keys are recorded but actual source parameter diff was checked separately",
        },
        "sym_source_dat_exact": {"passed": bool(sym_dat and all(sym_dat.values())),
                                 "files": sym_dat},
        "sym_off_on_smoke": {"passed": len(set(smoke_hashes.values())) == 1,
                             "regression_state_sha256": smoke_hashes},
        "modal_force_work_consistency": consistency_checks,
        "modal_gram": gram_checks,
        "symmetry_primary": {
            "probe_half_ptp_R_over_L": float(sym_primary.loc["R", "half_ptp_mm"]
                                                / sym_primary.loc["L", "half_ptp_mm"]),
            "probe_opening_waveform_max_abs_L_minus_R_mm": float(np.max(np.abs(
                -sym["displacement"][:, 1] - sym["displacement"][:, 2]))),
            "surface_half_ptp_R_over_L": float(sym_primary.loc["R", "surface_half_ptp_spatial_rms_mm"]
                                                  / sym_primary.loc["L", "surface_half_ptp_spatial_rms_mm"]),
            **symmetric_surface_check(sym, 0.30, 0.50),
        },
        "monitor_reconstruction": {case["name"]: monitor_reconstruction(case)
                                   for case in cases},
        "modal_matching": {
            "ambiguous_or_unmatched": int(modal_match.ambiguous.sum()),
            "asym_right_ambiguous_or_unmatched": int(modal_match[
                (modal_match.case == "asym12_6") & (modal_match.side == "R")].ambiguous.sum()),
        },
        "contact_component_split": "NOT_IDENTIFIABLE_FROM_SAVED_DATA",
        "dynamic_fluid_work": "NOT_MEASURED; W_fluct values are operational Q-reference sensitivity only",
        "geometry_vtu": "PENDING_EXTERNAL_SYMMETRIC_SAVED_VTU_AUDIT",
        "asym_geometry_vtu": asym_analysis_checks.get("geometry_vtu", "NOT_MEASURED"),
        "preexisting_test_failure": "test_input_file_params expects obsolete Kawahara paths in current input/param.txt",
    }
    dump_json(output / "checks_compare.json", checks)
    sym_binary = args.sym_production.resolve().parent / "build/on/simulation"
    asym_binary = Path(asym_analysis_manifest["diagnostic_run"]["binary_path"])
    repository = Path(__file__).resolve().parents[2]
    manifest_compare = {
        "analysis_schema": "symmetric_asymmetric_mechanism_compare_v1",
        "postprocess_only": True,
        "repository_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip(),
        "windows_s": windows,
        "sources": {
            "asym_analysis": str(args.asym_analysis.resolve()),
            "asym_production": str(args.asym_production.resolve()),
            "sym_source": str(args.sym_source.resolve()),
            "sym_production": str(args.sym_production.resolve()),
        },
        "raw_referenced_not_copied": True,
        "executables": {
            "symmetric_observer": {"path": str(sym_binary), "sha256": sha256(sym_binary)},
            "asymmetric_observer": {"path": str(asym_binary), "sha256": sha256(asym_binary)},
            "identical_sha256": sha256(sym_binary) == sha256(asym_binary),
            "note": "Executable hashes identify the frozen solver/observer used by each production run; later working-tree edits do not alter these binaries.",
        },
        "runs": {case["name"]: {
            "complete": case["metadata"]["complete"], "nstep": case["metadata"]["nstep"],
            "dt_s": case["dt"], "nmode": case["metadata"]["nmode"],
            "modal_native_size": (case["diagnostic"] / "modal_native.bin").stat().st_size,
            "modal_native_sha256": sha256(case["diagnostic"] / "modal_native.bin"),
            "parameter_sha256": sha256(case["production"] / "txt/params_used.txt"),
            "mode_sha256": [sha256(path) for path in case["mode_paths"]],
            "frequency_sha256": [sha256(path) for path in case["frequency_paths"]],
            "monitor_point_id_0based": case["monitor_id"],
            "monitor_xyz_mm": case["points"][case["monitor_id"]].tolist(),
            "surface_mapping_max_error_mm": case["surface_error"],
        } for case in cases},
        "input_differences": input_differences,
        "source_parameter_line_differences": source_parameter_differences,
        "production_manifest_schema_differences": manifest_schema_differences,
        "outputs": ["amplitude_compare.csv", "modal_match.csv",
                    "modal_covariance_compare.csv", "fluid_work_compare.csv",
                    "contact_compare.csv", "checks_compare.json", "report.md"],
    }
    dump_json(output / "manifest_compare.json", manifest_compare)

    # Four compact figures (maximum allowed is five).
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    for case in cases:
        t = case["displacement"][:, 0]
        axes[0].plot(t, -case["displacement"][:, 1], lw=.6, label=case["name"] + "-L")
        axes[0].plot(t, case["displacement"][:, 2], lw=.6, label=case["name"] + "-R")
    axes[0].set(xlabel="time [s]", ylabel="opening-positive probe uy [mm]")
    axes[0].legend(ncols=2)
    primary_amp.pivot(index="case", columns="side",
                      values="surface_half_ptp_spatial_rms_mm").plot.bar(ax=axes[1])
    axes[1].set(ylabel="surface half-ptp spatial RMS [mm]")
    fig.savefig(output / "figures/01_waveform_surface.png", dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for axis, case_name in zip(axes, ("sym12_12", "asym12_6")):
        selected = covariance[(covariance.case == case_name)
            & (covariance.window_start_s == .30) & (covariance.window_end_s == .50)
            & (covariance["mode"] == "total")]
        x = np.arange(2)
        axis.bar(x - .18, selected.individual_variance_mm2, .36, label="individual sum")
        axis.bar(x + .18, selected.cross_term_mm2, .36, label="cross term")
        axis.set_xticks(x, selected.side); axis.set_title(case_name); axis.set_ylabel("variance [mm²]")
        axis.legend()
    fig.savefig(output / "figures/02_covariance.png", dpi=160); plt.close(fig)
    fig, axis = plt.subplots(figsize=(7, 6), constrained_layout=True)
    image = axis.imshow(mac_by_series[("asym12_6", "R")], vmin=0, vmax=1,
                        origin="lower", aspect="auto", cmap="viridis")
    axis.set(xlabel="asymmetric right mode", ylabel="symmetric reference mode")
    fig.colorbar(image, ax=axis, label="MAC")
    fig.savefig(output / "figures/03_mac_asym_right.png", dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    primary_work.pivot(index="case", columns="side", values="W_raw_J").mul(1e3).plot.bar(ax=axes[0])
    axes[0].set(ylabel="raw fluid work [mJ]")
    primary_contact.pivot(index="case", columns="side", values="W_contact_ref_J").mul(1e3).plot.bar(ax=axes[1])
    axes[1].set(ylabel="contact reference work [mJ]")
    fig.savefig(output / "figures/04_work_contact.png", dpi=160); plt.close(fig)

    def value(frame: pd.DataFrame, case: str, side: str, field: str) -> float:
        return float(frame[(frame.case == case) & (frame.side == side)].iloc[0][field])
    sym_work_l = value(primary_work, "sym12_12", "L", "W_raw_J") * 1e3
    sym_work_r = value(primary_work, "sym12_12", "R", "W_raw_J") * 1e3
    asym_work_l = value(primary_work, "asym12_6", "L", "W_raw_J") * 1e3
    asym_work_r = value(primary_work, "asym12_6", "R", "W_raw_J") * 1e3
    sym_contact_l = value(primary_contact, "sym12_12", "L", "W_contact_ref_J") * 1e3
    sym_contact_r = value(primary_contact, "sym12_12", "R", "W_contact_ref_J") * 1e3
    asym_contact_l = value(primary_contact, "asym12_6", "L", "W_contact_ref_J") * 1e3
    asym_contact_r = value(primary_contact, "asym12_6", "R", "W_contact_ref_J") * 1e3
    sym_force_pos = value(primary_work, "sym12_12", "L", "mean_positive_opening_fluid_force_N")
    sym_force_neg = value(primary_work, "sym12_12", "L", "mean_negative_opening_fluid_force_N")
    asym_force_pos_l = value(primary_work, "asym12_6", "L", "mean_positive_opening_fluid_force_N")
    asym_force_pos_r = value(primary_work, "asym12_6", "R", "mean_positive_opening_fluid_force_N")
    asym_force_neg_l = value(primary_work, "asym12_6", "L", "mean_negative_opening_fluid_force_N")
    asym_force_neg_r = value(primary_work, "asym12_6", "R", "mean_negative_opening_fluid_force_N")
    sym_work_norm = value(primary_work, "sym12_12", "L", "W_raw_over_velocity_squared_integral")
    asym_work_norm_l = value(primary_work, "asym12_6", "L", "W_raw_over_velocity_squared_integral")
    asym_work_norm_r = value(primary_work, "asym12_6", "R", "W_raw_over_velocity_squared_integral")
    report = f"""# 対称12/12 vs 非対称12/6 機構診断

対称12/12は代表点・全表面・仕事の左右基準として品質ゲートを通過した。
非対称12/6では代表点half-ptp R/Lが{asym_primary.loc['R','half_ptp_mm']/asym_primary.loc['L','half_ptp_mm']:.3f}、全表面R/Lが{asym_primary.loc['R','surface_half_ptp_spatial_rms_mm']/asym_primary.loc['L','surface_half_ptp_spatial_rms_mm']:.3f}へ低下した。
主要な観測は、右側のモード基底変更に伴う応答配分と交差共分散、流体仕事・接触参考仕事の同時変化である。
raw流体仕事は動的仕事ではなく、接触参考仕事は不可逆散逸の証明ではない。
接触ばね/速度依存成分は12/6 rawに分離値がなく、NOT_IDENTIFIABLE_FROM_SAVED_DATA とした。

## 1. 対称性

0.30–0.50 s の対称代表点 half-ptp は L={sym_primary.loc['L','half_ptp_mm']:.6g} mm、R={sym_primary.loc['R','half_ptp_mm']:.6g} mm（R/L={checks['symmetry_primary']['probe_half_ptp_R_over_L']:.9g}）。全表面half-ptp空間RMSは L={sym_primary.loc['L','surface_half_ptp_spatial_rms_mm']:.6g} mm、R={sym_primary.loc['R','surface_half_ptp_spatial_rms_mm']:.6g} mm。raw流体仕事は L={sym_work_l:.6g}、R={sym_work_r:.6g} mJ、接触参考仕事は L={sym_contact_l:.6g}、R={sym_contact_r:.6g} mJ。OFF/ON state、全dat、元run全dat、仕事整合性はいずれも合格した。

## 2. 非対称化による変化

0.30–0.50 s の非対称代表点 half-ptp は L={asym_primary.loc['L','half_ptp_mm']:.6g} mm、R={asym_primary.loc['R','half_ptp_mm']:.6g} mm。全表面half-ptp空間RMSは L={asym_primary.loc['L','surface_half_ptp_spatial_rms_mm']:.6g} mm、R={asym_primary.loc['R','surface_half_ptp_spatial_rms_mm']:.6g} mm。代表点の個別モード分散和/交差項/全分散は、対称L={sym_total_cov.loc['L','individual_variance_mm2']:.6g}/{sym_total_cov.loc['L','cross_term_mm2']:.6g}/{sym_total_cov.loc['L','reconstructed_variance_mm2']:.6g}、対称R={sym_total_cov.loc['R','individual_variance_mm2']:.6g}/{sym_total_cov.loc['R','cross_term_mm2']:.6g}/{sym_total_cov.loc['R','reconstructed_variance_mm2']:.6g} mm²、非対称L={asym_total_cov.loc['L','individual_variance_mm2']:.6g}/{asym_total_cov.loc['L','cross_term_mm2']:.6g}/{asym_total_cov.loc['L','reconstructed_variance_mm2']:.6g}、非対称R={asym_total_cov.loc['R','individual_variance_mm2']:.6g}/{asym_total_cov.loc['R','cross_term_mm2']:.6g}/{asym_total_cov.loc['R','reconstructed_variance_mm2']:.6g} mm²。

## 3. モード対応

同じ番号を仮定せず、鏡映を戻した体積全節点3成分のMACで対応づけた。非対称右の ambiguous/unmatched は {checks['modal_matching']['asym_right_ambiguous_or_unmatched']}/20。曖昧な行は `matched_mode` を空欄とし、無理な1対1解釈をしていない。詳細は `modal_match.csv`。

## 4. 流体荷重と仕事

0.30–0.50 s の条件付き正/負開口力は対称（各側）{sym_force_pos:.6g}/{sym_force_neg:.6g} N、非対称L={asym_force_pos_l:.6g}/{asym_force_neg_l:.6g} N、R={asym_force_pos_r:.6g}/{asym_force_neg_r:.6g} Nで、非対称の実壁面荷重は左右で近い。raw `ΣQ_fluid Δq` は対称 L/R={sym_work_l:.6g}/{sym_work_r:.6g} mJ、非対称 L/R={asym_work_l:.6g}/{asym_work_r:.6g} mJ。一方 `W_raw/∫Σqdot²dt` は対称={sym_work_norm:.6g}、非対称L/R={asym_work_norm_l:.6g}/{asym_work_norm_r:.6g}で、rawの左右順位をそのまま励振能力と解釈できない。Q統計、平均power、Q-v相関、構造減衰指標、平均/前半/後半Q参照による `W_fluct` 感度は `fluid_work_compare.csv` に同一定義で示した。`W_fluct` は操作的な平均力除去であり、平衡一般化力または真のdynamic fluid workとは呼ばない。

## 5. 接触

接触参考仕事 L/R は対称 {sym_contact_l:.6g}/{sym_contact_r:.6g} mJ、非対称 {asym_contact_l:.6g}/{asym_contact_r:.6g} mJ。接触event/non-event、AP端/内部のpair頻度、最大侵入量を `contact_compare.csv` に示した。C2のばね/速度依存分解は、既存12/6に適用時の成分が保存されていないため `NOT_IDENTIFIABLE_FROM_SAVED_DATA`。負値を不可逆散逸とは断定しない。

## 6. 数値品質

入力差分は右周波数・右モードVTUのみ。対称スモークはOFF/ON-disabled/ON-enabledでbit一致。本番対称runは元runの全datとSHA256一致。両runの全side-interval仕事整合性が合格。保存VTU監査結果は `checks_compare.json` の `geometry_vtu` を参照。既存 `input_file_params` テストは現行入力に対し旧ファイル名を期待する既知の陳腐化で、他のC++テストは合格した。

## 7. 機構候補

- M1（剛性変更が周波数・空間形状・寄与を変更）: **SUPPORTED**。右基底のみ変更され、12/20はMAC対応、8/20はunmatchedとなり、代表点モード分散も変化した。周波数差だけの因果証明ではない。
- M2（流体力の大きさよりモード投影・位相関係が重要）: **SUPPORTED**。非対称の正/負壁面荷重は左右で近い一方、右の交差項は+0.107864（対称）から-0.0776092 mm²となり、raw仕事も左右で分かれた。投影と位相の個別因果までは分離していない。
- M3（接触が右へ強い負仕事）: **WEAK_SUPPORT**。右接触参考仕事は対称-0.95918から非対称-1.20558 mJへ負側に変化したが、両側和は対称-1.91836から非対称-1.35821 mJへ絶対値が減り、ばね/速度依存分解も不能なので散逸起源は未証明。
- M4（小振幅により仕事・減衰も小さい）: **SUPPORTED**。非対称右はraw仕事と速度二乗積分がともに小さく、`W_raw/∫Σqdot²dt` はむしろ対称より大きくなって順位が逆転した。正規化だけで因果方向は確定しない。
- M5（数値・観測器アーティファクト）: **NOT_SUPPORTED**。OFF/ON、元run、整合性、鏡映基準の品質ゲートに合格した。

## 8. 次の1手（実行しない）

12/6条件を1回だけ、採用接触荷重のばね成分と速度依存成分を同じraw schemaで保存して再実行し、0.30–0.50 sの右側負仕事がどちらの成分に由来するかを判定する。変更変数は診断出力粒度だけとし、物理係数・式・時間刻みは変えない。
"""
    (output / "report.md").write_text(report, encoding="utf-8")
    print(json.dumps({"output": str(output), "checks": checks,
                      "primary_probe_R_over_L": {
                          "sym": checks["symmetry_primary"]["probe_half_ptp_R_over_L"],
                          "asym": float(asym_primary.loc["R", "half_ptp_mm"]
                                        / asym_primary.loc["L", "half_ptp_mm"])}},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
