#!/usr/bin/env python3
"""Analyze one completed opt-in diagnostic run for left/right asymmetry."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


WINDOWS = ((0.10, 0.20), (0.20, 0.30), (0.30, 0.40),
           (0.40, 0.50), (0.30, 0.50))
MASS_KG = 1.07155257e-3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def load_frequency(path: Path, count: int) -> np.ndarray:
    values = np.loadtxt(path, comments="%", ndmin=2)
    if len(values) < count:
        raise ValueError(f"{path}: fewer than {count} modes")
    return values[:count, 0]


def load_vtu_basis(path: Path, frequencies: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Return internal xyz points, normalized internal modes, and cell arrays."""
    count = len(frequencies)
    components: dict[str, list[np.ndarray | None]] = {
        axis: [None] * count for axis in "XYZ"
    }
    points_raw: np.ndarray | None = None
    cells: dict[str, np.ndarray] = {}
    pattern = re.compile(r"Displacement_field,_([XYZ])-?component_@_([\d.eE+-]+)_Hz")
    stack: list[str] = []
    for event, element in ET.iterparse(path, events=("start", "end")):
        tag = element.tag.split("}")[-1]
        if event == "start":
            stack.append(tag)
            continue
        if tag != "DataArray":
            stack.pop()
            element.clear()
            continue
        name = element.attrib.get("Name", "")
        match = pattern.fullmatch(name)
        if match:
            axis, text_frequency = match.groups()
            frequency = float(text_frequency)
            index = int(np.argmin(np.abs(frequencies - frequency)))
            tolerance = max(1.0e-3, 5.1e-5 * frequencies[index])
            if abs(frequencies[index] - frequency) <= tolerance and components[axis][index] is None:
                components[axis][index] = np.fromstring(element.text or "", sep=" ")
        elif name == "Points" or (len(stack) >= 2 and stack[-2] == "Points"):
            points_raw = np.fromstring(element.text or "", sep=" ").reshape(-1, 3)
        elif name in {"connectivity", "offsets", "types"}:
            cells[name] = np.fromstring(element.text or "", sep=" ", dtype=np.int64)
        element.clear()
        stack.pop()
    if points_raw is None:
        raise ValueError(f"{path}: missing Points")
    if any(value is None for axis in "XYZ" for value in components[axis]):
        raise ValueError(f"{path}: missing one of the selected modal components")
    npoints = len(points_raw)
    raw = np.stack(
        [np.stack([components[axis][mode] for axis in "XYZ"], axis=1)
         for mode in range(count)]
    )
    if raw.shape != (count, npoints, 3):
        raise ValueError(f"{path}: unexpected modal shape {raw.shape}")
    scale = 1.0 / np.sqrt((MASS_KG / npoints) * np.sum(raw * raw, axis=(1, 2)))
    raw *= scale[:, None, None]
    # COMSOL raw X/Y/Z -> solver internal uy/uz/ux, respectively.
    modes = raw[:, :, [2, 0, 1]]
    points = points_raw[:, [2, 0, 1]]
    return points, modes, cells


def load_surface_ids(path: Path, points: np.ndarray) -> tuple[np.ndarray, float]:
    nodes: dict[int, tuple[float, float, float]] = {}
    quads: list[list[int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = [field.strip() for field in line.split(",")]
        if fields[0] == "GRID":
            gx, gy, gz = (float(value) for value in fields[3:6])
            nodes[int(fields[1])] = (gz, gx, gy)
        elif fields[0] == "CQUAD4":
            quads.append([int(value) for value in fields[3:7]])
    used = sorted({node for quad in quads for node in quad})
    coordinates = np.asarray([nodes[node] for node in used])
    distances, indices = cKDTree(points).query(coordinates)
    if np.max(distances) >= 1.0e-3:
        raise ValueError(f"NAS/VTU surface mapping exceeds tolerance: {np.max(distances)} mm")
    if len(np.unique(indices)) != len(indices):
        raise ValueError("NAS surface mapping is not one-to-one")
    return indices.astype(np.int64), float(np.max(distances))


def modal_gram(modes: np.ndarray) -> np.ndarray:
    flattened = modes.reshape(len(modes), -1)
    return (MASS_KG / modes.shape[1]) * (flattened @ flattened.T)


def summarize_monitor(
    raw: np.ndarray,
    times: np.ndarray,
    side_index: int,
    side: str,
    phi_y: np.ndarray,
    frequencies: np.ndarray,
    displacement: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, float], np.ndarray]:
    mask = (times >= 0.30) & (times < 0.50)
    q = np.asarray(raw[mask, side_index, :, 4])
    contributions = q * phi_y[None, :] * 1.0e3
    total = contributions.sum(axis=1)
    covariance = np.cov(contributions, rowvar=False, ddof=0)
    total_variance = float(np.var(total))
    individual_variance = np.diag(covariance)
    covariance_total = float(total_variance - np.sum(individual_variance))
    rows = []
    for mode in range(len(frequencies)):
        values = contributions[:, mode]
        pair_covariance = covariance[mode, :].copy()
        pair_covariance[mode] = np.nan
        strongest_cancel = int(np.nanargmin(pair_covariance))
        strongest_reinforce = int(np.nanargmax(pair_covariance))
        rows.append(
            {
                "side": side,
                "mode": mode + 1,
                "source_mode": mode + 1,
                "frequency_hz": frequencies[mode],
                "probe_point_id_0based": None,
                "probe_phi_y_normalized": phi_y[mode],
                "mean_contribution_mm": np.mean(values),
                "half_ptp_contribution_mm": np.ptp(values) / 2.0,
                "ac_rms_contribution_mm": np.std(values),
                "variance_contribution_mm2": individual_variance[mode],
                "twice_covariance_with_other_modes_mm2": 2.0 * (
                    np.sum(covariance[mode, :]) - covariance[mode, mode]
                ),
                "strongest_canceling_partner_mode": strongest_cancel + 1,
                "strongest_canceling_pair_covariance_mm2": pair_covariance[strongest_cancel],
                "strongest_reinforcing_partner_mode": strongest_reinforce + 1,
                "strongest_reinforcing_pair_covariance_mm2": pair_covariance[strongest_reinforce],
                "all_mode_variance_sum_mm2": np.sum(individual_variance),
                "all_cross_covariance_mm2": covariance_total,
                "reconstructed_total_variance_mm2": total_variance,
            }
        )
    frame = pd.DataFrame(rows)
    frame["rms_rank"] = frame["ac_rms_contribution_mm"].rank(
        ascending=False, method="min"
    ).astype(int)
    saved_mask = (displacement[:, 0] >= 0.30) & (displacement[:, 0] < 0.50)
    saved = displacement[saved_mask]
    indices = np.rint(saved[:, 0] / (times[1] - times[0])).astype(int)
    reconstructed = np.sum(
        np.asarray(raw[indices, side_index, :, 4]) * phi_y[None, :] * 1.0e3,
        axis=1,
    )
    expected = saved[:, 1 if side == "L" else 2]
    error = reconstructed - expected
    checks = {
        "max_abs_error_mm": float(np.max(np.abs(error))),
        "rms_error_mm": float(np.sqrt(np.mean(error * error))),
        "max_abs_error_over_saved_ptp": float(np.max(np.abs(error)) / np.ptp(expected)),
        "sum_individual_variance_mm2": float(np.sum(individual_variance)),
        "cross_covariance_mm2": covariance_total,
        "total_variance_mm2": total_variance,
        "cancellation_fraction": float(-covariance_total / np.sum(individual_variance))
        if covariance_total < 0 else 0.0,
    }
    return frame, checks, contributions


def surface_amplitudes(
    raw: np.ndarray,
    times: np.ndarray,
    surface_ids: np.ndarray,
    points: np.ndarray,
    modes_by_side: list[np.ndarray],
) -> tuple[pd.DataFrame, dict[str, object]]:
    mask = (times >= 0.30) & (times < 0.50)
    rows_by_side: dict[str, pd.DataFrame] = {}
    z_values = points[surface_ids, 2]
    z_levels = np.unique(np.round(z_values, 6))
    z_index = np.argmin(np.abs(z_values[:, None] - z_levels[None, :]), axis=1)
    ap_region = np.where((z_index < 3) | (z_index >= len(z_levels) - 3), "end", "interior")
    for side_index, side in enumerate(("L", "R")):
        q = np.asarray(raw[mask, side_index, :, 4])
        phi = modes_by_side[side_index][:, surface_ids, 1]
        mean = np.zeros(len(surface_ids))
        rms = np.zeros(len(surface_ids))
        half_ptp = np.zeros(len(surface_ids))
        for start in range(0, len(surface_ids), 256):
            stop = min(start + 256, len(surface_ids))
            displacement = q @ phi[:, start:stop] * 1.0e3
            mean[start:stop] = np.mean(displacement, axis=0)
            rms[start:stop] = np.std(displacement, axis=0)
            half_ptp[start:stop] = np.ptp(displacement, axis=0) / 2.0
        rows_by_side[side] = pd.DataFrame(
            {
                "side": side,
                "point_id_0based": surface_ids,
                "x_mm": points[surface_ids, 0],
                "y_initial_mm": points[surface_ids, 1],
                "z_mm": z_values,
                "ap_level_index": z_index,
                "ap_region": ap_region,
                "mean_uy_mm": mean,
                "ac_rms_uy_mm": rms,
                "half_ptp_uy_mm": half_ptp,
            }
        )
    left = rows_by_side["L"]
    right = rows_by_side["R"]
    ratio_ptp = right.half_ptp_uy_mm.to_numpy() / np.maximum(
        left.half_ptp_uy_mm.to_numpy(), 1.0e-15
    )
    ratio_rms = right.ac_rms_uy_mm.to_numpy() / np.maximum(
        left.ac_rms_uy_mm.to_numpy(), 1.0e-15
    )
    for frame in rows_by_side.values():
        frame["paired_half_ptp_r_over_l"] = ratio_ptp
        frame["paired_ac_rms_r_over_l"] = ratio_rms
    combined = pd.concat([rows_by_side["L"], rows_by_side["R"]], ignore_index=True)
    summary: dict[str, object] = {
        "point_count_per_side": len(surface_ids),
        "half_ptp_spatial_rms_mm": {
            side: float(np.sqrt(np.mean(frame.half_ptp_uy_mm.to_numpy() ** 2)))
            for side, frame in rows_by_side.items()
        },
        "ac_rms_spatial_rms_mm": {
            side: float(np.sqrt(np.mean(frame.ac_rms_uy_mm.to_numpy() ** 2)))
            for side, frame in rows_by_side.items()
        },
        "median_paired_half_ptp_r_over_l": float(np.median(ratio_ptp)),
        "p05_paired_half_ptp_r_over_l": float(np.percentile(ratio_ptp, 5)),
        "p95_paired_half_ptp_r_over_l": float(np.percentile(ratio_ptp, 95)),
        "fraction_points_right_larger_half_ptp": float(np.mean(ratio_ptp > 1.0)),
        "global_half_ptp_spatial_rms_r_over_l": float(
            np.sqrt(np.mean(right.half_ptp_uy_mm.to_numpy() ** 2))
            / np.sqrt(np.mean(left.half_ptp_uy_mm.to_numpy() ** 2))
        ),
        "global_ac_rms_spatial_rms_r_over_l": float(
            np.sqrt(np.mean(right.ac_rms_uy_mm.to_numpy() ** 2))
            / np.sqrt(np.mean(left.ac_rms_uy_mm.to_numpy() ** 2))
        ),
    }
    return combined, summary


def modal_work(raw: np.ndarray, times: np.ndarray, frequencies: list[np.ndarray]) -> pd.DataFrame:
    rows = []
    dt = times[1] - times[0]
    for start, end in WINDOWS:
        mask = (times >= start) & (times < end)
        for side_index, side in enumerate(("L", "R")):
            data = np.asarray(raw[mask, side_index, :, :])
            for mode in range(data.shape[1]):
                qpos, qneg, qfluid, qcontact, q, velocity, delta_q = data[:, mode, :].T
                rows.append(
                    {
                        "side": side,
                        "mode": mode + 1,
                        "source_mode": mode + 1,
                        "frequency_hz": frequencies[side_index][mode],
                        "window_start_s": start,
                        "window_end_s": end,
                        "interval_count": len(data),
                        "W_pos_J": np.sum(qpos * delta_q),
                        "W_neg_J": np.sum(qneg * delta_q),
                        "W_fluid_raw_J": np.sum(qfluid * delta_q),
                        "Q_reference": np.nan,
                        "W_static_J": np.nan,
                        "W_fluid_dynamic_J": np.nan,
                        "dynamic_work_status": "not_computed_no_justified_equilibrium_reference",
                        "D_struct_normalized": np.sum(
                            2.0 * 0.06 * 2.0 * np.pi * frequencies[side_index][mode]
                            * velocity * velocity * dt
                        ),
                        "W_contact_reference_J": np.sum(qcontact * delta_q),
                        "q_rms": np.sqrt(np.mean(q * q)),
                        "qdot_rms": np.sqrt(np.mean(velocity * velocity)),
                        "delta_q_sum": np.sum(delta_q),
                    }
                )
    return pd.DataFrame(rows)


def contact_geometry(production: Path) -> pd.DataFrame:
    diagnostic = production / "kanaya_diagnostic"
    flow = pd.read_csv(diagnostic / "flow_geometry.csv")
    contact = pd.read_csv(diagnostic / "contact_regions.csv")
    rows = []
    for start, end in WINDOWS:
        selected_flow = flow[(flow.time_s >= start) & (flow.time_s < end)]
        for side in ("L", "R"):
            for region in ("end", "interior"):
                selected = contact[
                    (contact.time_s >= start)
                    & (contact.time_s < end)
                    & (contact.side == side)
                    & (contact.ap_region == region)
                ]
                rows.append(
                    {
                        "window_start_s": start,
                        "window_end_s": end,
                        "side": side,
                        "ap_region": region,
                        "samples": len(selected),
                        "contact_pair_fraction": np.mean(selected.pairs > 0),
                        "mean_pairs": np.mean(selected.pairs),
                        "max_pairs": np.max(selected.pairs),
                        "max_penetration_m": np.max(selected.max_penetration_m),
                        "mean_opening_contact_force_N": np.mean(selected.opening_Fy_N),
                        "Amin_min_mm2": np.min(selected_flow.Amin_mm2),
                        "Amin_mean_mm2": np.mean(selected_flow.Amin_mm2),
                        "area_closure_fraction_le_0p01mm2": np.mean(selected_flow.Amin_mm2 <= 0.01),
                        "monitor_gap_min_mm": np.min(selected_flow.monitor_gap_mm),
                        "monitor_gap_max_mm": np.max(selected_flow.monitor_gap_mm),
                        "mean_delta_gap_mm": np.mean(selected_flow.mean_delta_gap_mm),
                        "contact_flag_fraction": np.mean(selected_flow.contact_flag != 0),
                        "flow_source": str(diagnostic / "flow_geometry.csv"),
                        "contact_source": str(diagnostic / "contact_regions.csv"),
                    }
                )
    return pd.DataFrame(rows)


def plot_surface(path: Path, surface: pd.DataFrame, monitor_z: float) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for side, group in surface.groupby("side"):
        summary = group.groupby("ap_level_index").agg(
            z_mm=("z_mm", "mean"),
            median_half_ptp=("half_ptp_uy_mm", "median"),
            p95_half_ptp=("half_ptp_uy_mm", lambda value: np.percentile(value, 95)),
        )
        axes[0].plot(summary.z_mm, summary.median_half_ptp, label=f"{side} median")
        axes[0].plot(summary.z_mm, summary.p95_half_ptp, linestyle="--", label=f"{side} P95")
    axes[0].axvline(monitor_z, color="black", alpha=0.4, label="monitor z")
    axes[0].set(xlabel="AP z [mm]", ylabel="surface half-ptp uy [mm]")
    axes[0].legend(ncols=2)
    paired = surface[surface.side == "L"]
    axes[1].scatter(
        paired.z_mm,
        paired.paired_half_ptp_r_over_l,
        s=5,
        alpha=0.35,
    )
    axes[1].axhline(1.0, color="black", linewidth=0.8)
    axes[1].set(xlabel="AP z [mm]", ylabel="paired surface half-ptp R/L", ylim=(0, 2.5))
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_modal(path: Path, modal: pd.DataFrame) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for axis, side in zip(axes, ("L", "R")):
        selected = modal[modal.side == side].sort_values("ac_rms_contribution_mm", ascending=False)
        colors = ["tab:red" if value < 0 else "tab:blue" for value in selected.twice_covariance_with_other_modes_mm2]
        axis.bar(selected["mode"].astype(str), selected.ac_rms_contribution_mm, color=colors)
        axis.set(title=side, xlabel="mode (ranked by RMS order in table)", ylabel="probe modal AC RMS [mm]")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def plot_work(path: Path, work: pd.DataFrame, contact: pd.DataFrame) -> None:
    selected = work[(work.window_start_s == 0.30) & (work.window_end_s == 0.50)]
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    totals = selected.groupby("side")[["W_pos_J", "W_neg_J", "W_fluid_raw_J", "D_struct_normalized", "W_contact_reference_J"]].sum()
    totals.plot.bar(ax=axes[0])
    axes[0].set(ylabel="integrated algebraic value", title="0.30-0.50 s all modes")
    c = contact[(contact.window_start_s == 0.30) & (contact.window_end_s == 0.50)]
    labels = c.side + "-" + c.ap_region
    axes[1].bar(labels, c.contact_pair_fraction)
    axes[1].set(ylabel="fraction with adopted contact pairs", title="contact by AP region", ylim=(0, 1))
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument("--production-run", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    arguments = parser.parse_args()
    source = arguments.source_run.resolve()
    production = arguments.production_run.resolve()
    artifact = arguments.artifact_dir.resolve()
    diagnostic = production / "kanaya_diagnostic"
    metadata = json.loads((diagnostic / "run_metadata.json").read_text(encoding="utf-8"))
    nstep = int(metadata["integrated_intervals"])
    nmode = metadata["nmode"][0]
    dt = float(metadata["dt_s"])
    binary = diagnostic / "modal_native.bin"
    expected_bytes = nstep * 2 * nmode * 7 * 8
    if binary.stat().st_size != expected_bytes:
        raise ValueError("modal_native.bin size disagrees with metadata")
    raw = np.memmap(binary, dtype="<f8", mode="r", shape=(nstep, 2, nmode, 7))
    times = np.arange(nstep) * dt

    manifest_text = (production / "txt/manifest.txt").read_text(encoding="utf-8")
    manifest = dict(
        (key.strip(), value.strip())
        for key, value in (line.split("=", 1) for line in manifest_text.splitlines() if "=" in line)
    )
    frequency_paths = [Path(manifest["left_frequency"]), Path(manifest["right_frequency"])]
    mode_paths = [Path(manifest["left_mode_vtu"]), Path(manifest["right_mode_vtu"])]
    frequencies = [load_frequency(path, nmode) for path in frequency_paths]
    points_left, modes_left, cells = load_vtu_basis(mode_paths[0], frequencies[0])
    points_right, modes_right, _ = load_vtu_basis(mode_paths[1], frequencies[1])
    if not np.array_equal(points_left, points_right):
        raise ValueError("left/right geometry point clouds differ")
    modes_right[:, :, 1] *= -1.0
    surface_ids, surface_mapping_error = load_surface_ids(
        Path(manifest["left_surface_nas"]), points_left
    )
    monitor_target = np.array([9.18767, 8.37084, 8.37681])
    monitor_id = int(surface_ids[np.argmin(np.linalg.norm(points_left[surface_ids] - monitor_target, axis=1))])
    monitor_error = float(np.linalg.norm(points_left[monitor_id] - monitor_target))
    if monitor_error > 1.0e-4:
        raise ValueError("could not recover archived monitor point ID")

    displacement = np.loadtxt(production / "dat/displace.dat", comments="#")
    modal_frames = []
    modal_checks = {}
    contributions = []
    for side_index, (side, modes) in enumerate((("L", modes_left), ("R", modes_right))):
        frame, check, contribution = summarize_monitor(
            raw, times, side_index, side, modes[:, monitor_id, 1], frequencies[side_index], displacement
        )
        frame["probe_point_id_0based"] = monitor_id
        modal_frames.append(frame)
        modal_checks[side] = check
        contributions.append(contribution)
    modal = pd.concat(modal_frames, ignore_index=True)
    modal.to_csv(artifact / "modal_contribution.csv", index=False)

    surface, surface_summary = surface_amplitudes(
        raw, times, surface_ids, points_left, [modes_left, modes_right]
    )
    surface.to_csv(artifact / "surface_amplitude.csv", index=False)
    work = modal_work(raw, times, frequencies)
    work.to_csv(artifact / "modal_work.csv", index=False)
    contact = contact_geometry(production)
    contact.to_csv(artifact / "contact_geometry.csv", index=False)

    consistency = pd.read_csv(diagnostic / "consistency.csv")
    included = consistency[consistency.included == 1]
    gram = {"L": modal_gram(modes_left), "R": modal_gram(modes_right)}
    gram_checks = {}
    for side, values in gram.items():
        off_diagonal = values - np.diag(np.diag(values))
        gram_checks[side] = {
            "max_diagonal_abs_error": float(np.max(np.abs(np.diag(values) - 1.0))),
            "max_off_diagonal_abs": float(np.max(np.abs(off_diagonal))),
        }

    original_dat = source / "output/dat"
    production_dat = production / "dat"
    exact_outputs = {}
    for path in sorted(original_dat.glob("*.dat")):
        counterpart = production_dat / path.name
        exact_outputs[path.name] = counterpart.is_file() and sha256(path) == sha256(counterpart)
    smoke = artifact / "smoke"
    smoke_hashes = {
        name: sha256(smoke / name / "regression_state.bin")
        for name in ("off", "on_disabled", "on_enabled")
    }
    checks = json.loads((artifact / "checks.json").read_text(encoding="utf-8"))
    checks["data_gaps"] = [
        "W_fluid_dynamic is unmeasured because no independently justified equilibrium Q_reference is available.",
        "Fundamental frequency/cycle amplitudes remain unconfirmed for the irregular archived waveform.",
        "VTU geometry is saved only for 0.45-0.4998 s; the online AP-order check covers all steps.",
    ]
    checks.update(
        {
            "off_on_consistency": {
                "passed": len(set(smoke_hashes.values())) == 1,
                "regression_state_sha256": smoke_hashes,
                "full_run_dat_exact_by_sha256": exact_outputs,
            },
            "modal_force_work_consistency": {
                "passed": bool(
                    len(included) == 2 * nstep
                    and included.valid.eq(1).all()
                    and included.scaled_max.max() <= 1.0
                ),
                "included_side_intervals": len(included),
                "excluded_side_intervals": int(np.sum(consistency.included == 0)),
                "maximum_scaled_residual": float(included.scaled_max.max()),
                "maximum_force_abs_N": float(included.force_abs.max()),
                "maximum_Q_abs": float(included.Q_abs.max()),
                "maximum_work_split_abs_J": float(included.work_split_abs.max()),
                "maximum_nodal_modal_work_abs_J": float(included.nodal_work_abs.max()),
                "maximum_contact_Q_abs": float(included.contact_Q_abs.max()),
            },
            "monitor_reconstruction": modal_checks,
            "surface_amplitude": surface_summary,
            "surface_mapping": {
                "point_count": len(surface_ids),
                "maximum_nas_vtu_distance_mm": surface_mapping_error,
                "monitor_point_id_0based": monitor_id,
                "monitor_coordinate_error_mm": monitor_error,
            },
            "modal_gram": gram_checks,
            "dynamic_work": "not computed: no independently justified equilibrium Q_reference",
            "geometry_vtu": checks.get(
                "geometry_vtu", "pending separate saved-frame signed-volume audit"
            ),
        }
    )
    dump_json(artifact / "checks.json", checks)

    analysis_manifest_path = artifact / "manifest.json"
    analysis_manifest = json.loads(analysis_manifest_path.read_text(encoding="utf-8"))
    analysis_manifest["monitor"]["l_point_id_0based"] = monitor_id
    analysis_manifest["monitor"]["r_point_id_0based"] = monitor_id
    analysis_manifest["monitor"]["point_id_note"] = (
        "Recovered by one-to-one NAS-to-VTU mapping and archived monitor coordinates."
    )
    analysis_manifest["diagnostic_run"] = {
        "path": str(production),
        "observer_compile_option": "ENABLE_KANAYA_DIAGNOSTICS=ON",
        "observer_runtime_option": "SIMULATION_KANAYA_DIAGNOSTICS=1",
        "binary_path": str(artifact / "build/on/simulation"),
        "binary_sha256": sha256(artifact / "build/on/simulation"),
        "modal_native_path": str(binary),
        "modal_native_sha256": sha256(binary),
        "run_metadata_path": str(diagnostic / "run_metadata.json"),
        "integrated_intervals": nstep,
        "excluded_final_intervals": int(metadata["excluded_intervals"]),
        "raw_data_copied": False,
        "full_run_dat_exact_to_archived": bool(all(exact_outputs.values())),
    }
    analysis_manifest["analysis_quality"].update(
        {
            "modal_native_available": True,
            "off_on_smoke_exact": len(set(smoke_hashes.values())) == 1,
            "full_run_dat_exact_to_archived": bool(all(exact_outputs.values())),
            "modal_nodal_consistency_passed": checks["modal_force_work_consistency"]["passed"],
            "monitor_reconstruction_passed": all(
                value["max_abs_error_mm"] <= 5.1e-6 for value in modal_checks.values()
            ),
            "dynamic_work_status": "unmeasured_no_justified_equilibrium_reference",
        }
    )
    analysis_manifest["derived_artifacts"] = [
        "amplitude_by_window.csv",
        "modal_contribution.csv",
        "surface_amplitude.csv",
        "modal_work.csv",
        "contact_geometry.csv",
        "fluid_load_summary.csv",
        "geometry_audit.csv",
        "checks.json",
        "report.md",
        "figures/01_waveform_window_ratios.png",
        "figures/02_surface_amplitude.png",
        "figures/03_modal_contribution.png",
        "figures/04_work_contact.png",
    ]
    dump_json(analysis_manifest_path, analysis_manifest)

    plot_surface(artifact / "figures/02_surface_amplitude.png", surface, points_left[monitor_id, 2])
    plot_modal(artifact / "figures/03_modal_contribution.png", modal)
    plot_work(artifact / "figures/04_work_contact.png", work, contact)
    print(json.dumps({
        "monitor_id": monitor_id,
        "monitor_reconstruction": modal_checks,
        "surface": surface_summary,
        "consistency_max": float(included.scaled_max.max()),
        "gram": gram_checks,
    }, indent=2))


if __name__ == "__main__":
    main()
