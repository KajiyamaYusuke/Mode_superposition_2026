#!/usr/bin/env python3
"""Postprocess left/right monitor displacement without changing the solver.

The script deliberately treats an existing run as read-only.  It writes all
derived files to an explicitly supplied output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_WINDOWS = ((0.10, 0.20), (0.20, 0.30), (0.30, 0.40),
                   (0.40, 0.50), (0.30, 0.50))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def parse_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip('"')
    return values


def parse_monitor(log_path: Path) -> dict[str, object]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    grid = re.search(r"Monitor surface point \(i,j\)=\((\d+),\s*(\d+)\)", text)
    result: dict[str, object] = {
        "target_x_mm": 9.2,
        "target_z_mm": 8.5,
        "surface_i": int(grid.group(1)) if grid else None,
        "surface_j": int(grid.group(2)) if grid else None,
    }
    for side in ("L", "R"):
        found = re.search(
            rf"Monitor Node {side} idx=([-+\d.eE]+),\s*([-+\d.eE]+),\s*([-+\d.eE]+)",
            text,
        )
        result[f"{side.lower()}_point_id_0based"] = None
        result[f"{side.lower()}_initial_xyz_mm"] = (
            [float(found.group(i)) for i in range(1, 4)] if found else None
        )
    result["point_id_note"] = (
        "The archived log prints coordinates after the label 'idx' but does not print the node ID."
    )
    return result


def git_value(root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", *arguments], cwd=root, text=True, stderr=subprocess.DEVNULL
    ).strip()


def load_displacement(path: Path) -> np.ndarray:
    values = np.loadtxt(path, comments="#", ndmin=2)
    if values.shape[1] != 3:
        raise ValueError(f"{path}: expected time, uyL, uyR; got shape {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{path}: non-finite monitor displacement")
    if np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError(f"{path}: time is not strictly increasing")
    return values


def signal_metrics(values: np.ndarray) -> dict[str, float]:
    centered = values - np.mean(values)
    p5, p95 = np.percentile(values, [5, 95])
    return {
        "half_ptp_mm": float(np.ptp(values) / 2.0),
        "mean_mm": float(np.mean(values)),
        "median_mm": float(np.median(values)),
        "std_mm": float(np.std(values, ddof=0)),
        "ac_rms_mm": float(np.sqrt(np.mean(centered * centered))),
        "half_p95_p5_mm": float((p95 - p5) / 2.0),
    }


def ratio(numerator: float, denominator: float) -> float | None:
    floor = 1.0e-12
    return float(numerator / denominator) if abs(denominator) > floor else None


def analyze_windows(data: np.ndarray) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    ratios: list[dict[str, object]] = []
    for start, end in DEFAULT_WINDOWS:
        mask = (data[:, 0] >= start) & (data[:, 0] <= end)
        if np.count_nonzero(mask) < 5:
            raise ValueError(f"window {start:.2f}-{end:.2f} s has fewer than five samples")
        metrics = {
            "L": signal_metrics(data[mask, 1]),
            "R": signal_metrics(data[mask, 2]),
        }
        window = f"{start:.2f}-{end:.2f}"
        ratio_values = {
            "half_ptp_r_over_l": ratio(metrics["R"]["half_ptp_mm"], metrics["L"]["half_ptp_mm"]),
            "ac_rms_r_over_l": ratio(metrics["R"]["ac_rms_mm"], metrics["L"]["ac_rms_mm"]),
            "half_p95_p5_r_over_l": ratio(
                metrics["R"]["half_p95_p5_mm"], metrics["L"]["half_p95_p5_mm"]
            ),
        }
        quality = "ok"
        if min(metrics[side]["half_ptp_mm"] for side in ("L", "R")) < 1.0e-4:
            quality = "insufficient_signal"
        for side in ("L", "R"):
            rows.append(
                {
                    "side": side,
                    "window_start_s": start,
                    "window_end_s": end,
                    "window": window,
                    "probe": "simulation_monitor_i40_j34",
                    **metrics[side],
                    **ratio_values,
                    "cycle_amplitude_mm": None,
                    "cycle_quality": "not_reported_archived_f0_unconfirmed",
                    "quality_flag": quality,
                }
            )
        ratios.append({"window": window, **ratio_values})
    return rows, ratios


def write_amplitude_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(path: Path, data: np.ndarray, ratios: list[dict[str, object]]) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    axes[0].plot(data[:, 0], -data[:, 1], linewidth=0.8, label="L opening-positive (-uyL)")
    axes[0].plot(data[:, 0], data[:, 2], linewidth=0.8, label="R opening-positive (+uyR)")
    axes[0].axvspan(0.30, 0.50, alpha=0.12, color="black", label="primary window")
    axes[0].set(xlabel="time [s]", ylabel="monitor displacement [mm]")
    axes[0].legend(loc="upper left", ncols=3)
    labels = [entry["window"] for entry in ratios]
    positions = np.arange(len(labels))
    axes[1].plot(positions, [entry["half_ptp_r_over_l"] for entry in ratios], "o-", label="half-ptp")
    axes[1].plot(positions, [entry["ac_rms_r_over_l"] for entry in ratios], "s-", label="AC RMS")
    axes[1].plot(
        positions,
        [entry["half_p95_p5_r_over_l"] for entry in ratios],
        "^-",
        label="half(P95-P5)",
    )
    axes[1].set_xticks(positions, labels, rotation=20)
    axes[1].set(xlabel="evaluation window [s]", ylabel="R/L ratio")
    axes[1].legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()

    run_dir = arguments.run_dir.resolve()
    output_dir = arguments.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True)
    (output_dir / "figures").mkdir()

    repository = Path(__file__).resolve().parents[2]
    displacement_path = run_dir / "output/dat/displace.dat"
    manifest_path = run_dir / "output/txt/manifest.txt"
    parameter_path = run_dir / "output/txt/params_used.txt"
    status_path = run_dir / "status.json"
    log_path = run_dir / "stdout.log"
    required = [displacement_path, manifest_path, parameter_path, status_path, log_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing archived run files: {missing}")

    data = load_displacement(displacement_path)
    rows, ratios = analyze_windows(data)
    primary = next(entry for entry in ratios if entry["window"] == "0.30-0.50")
    expected = 0.3023554621751491
    observed = primary["half_ptp_r_over_l"]
    difference = abs(float(observed) - expected)
    stage1_passed = difference <= 5.0e-7
    write_amplitude_csv(output_dir / "amplitude_by_window.csv", rows)
    make_figure(output_dir / "figures/01_waveform_window_ratios.png", data, ratios)

    archived_manifest = parse_key_values(manifest_path)
    status = json.loads(status_path.read_text(encoding="utf-8"))
    input_paths = {
        "parameter_file": parameter_path,
        "left_mode_vtu": Path(archived_manifest["left_mode_vtu"]),
        "right_mode_vtu": Path(archived_manifest["right_mode_vtu"]),
        "left_frequency": Path(archived_manifest["left_frequency"]),
        "right_frequency": Path(archived_manifest["right_frequency"]),
        "left_surface_nas": Path(archived_manifest["left_surface_nas"]),
        "right_surface_nas": Path(archived_manifest["right_surface_nas"]),
    }
    hashes = {
        name: {"path": str(path), "sha256": sha256(path)}
        for name, path in input_paths.items()
    }
    diagnostics = run_dir / "output/kanaya_diagnostic"
    diagnostic_files = sorted(path.name for path in diagnostics.glob("*")) if diagnostics.is_dir() else []
    result_files = sorted((run_dir / "output/result").glob("deform_combined*.vtu"))
    manifest = {
        "analysis": "right_body_6kpa_900pa_asymmetry",
        "analysis_created_date": "2026-10-29",
        "repository": {
            "path": str(repository),
            "branch_at_analysis": git_value(repository, "branch", "--show-current"),
            "head_at_analysis": git_value(repository, "rev-parse", "HEAD"),
            "status_at_analysis": git_value(repository, "status", "--short"),
            "archived_run_branch": "codex/time-dependent-pressure-sweep-20261025",
            "archived_run_tracked_diff_sha256": "e4d5b31f900ea5873c77f76a4d00a4703de63a489ab190a030411ab55f829d90",
            "archived_run_executable_sha256": "647629a07828efbd1c3743aabfdc9c4f21ca17344fd6dbf2f4f3cbf7f9e6fd3b",
        },
        "source_run": str(run_dir),
        "run_status": status.get("status"),
        "return_code": status.get("return_code"),
        "elapsed_s": status.get("elapsed_s"),
        "conditions": {
            key: archived_manifest.get(key)
            for key in (
                "nstep", "dt_s", "nmode", "mode_selection", "zetaL", "zetaR",
                "kc1", "kc2", "kc3", "ncont", "target_initial_minimum_gap_mm",
                "pressure_ramp_time_sec", "pressure_schedule_final_pa",
                "glottal_method_profile", "write_deformation_vtu",
            )
        },
        "monitor": parse_monitor(log_path),
        "inputs": hashes,
        "outputs": {
            "displace": {"path": str(displacement_path), "sha256": sha256(displacement_path)},
            "displacement_rows": int(len(data)),
            "displacement_interval_s": float(np.median(np.diff(data[:, 0]))),
            "deformation_vtu_count": len(result_files),
            "deformation_vtu_time_scope": "last 0.05 s only (solver saves when nstep-n <= 5000)",
            "kanaya_diagnostic_present": diagnostics.is_dir(),
            "kanaya_diagnostic_files": diagnostic_files,
        },
        "evaluation_windows_s": [list(window) for window in DEFAULT_WINDOWS],
        "units": {"time": "s", "monitor_displacement": "mm", "ratios": "R/L"},
        "analysis_quality": {
            "stage1_ratio_reproduced": stage1_passed,
            "archived_fundamental_frequency": "unconfirmed/candidate only",
            "surface_scope_available": "0.45-0.50 s saved VTU only",
            "modal_native_available": diagnostics.is_dir() and (diagnostics / "modal_native.bin").is_file(),
        },
    }
    write_json(output_dir / "manifest.json", manifest)
    checks = {
        "stage1": {
            "passed": stage1_passed,
            "expected_half_ptp_r_over_l": expected,
            "observed_half_ptp_r_over_l": observed,
            "absolute_difference": difference,
            "tolerance": 5.0e-7,
            "finite_monitor_data": bool(np.all(np.isfinite(data))),
            "strictly_increasing_time": bool(np.all(np.diff(data[:, 0]) > 0)),
            "last_saved_time_s": float(data[-1, 0]),
        },
        "modal_force_work_consistency": "not_measured_no_target_kanaya_diagnostic",
        "off_on_consistency": "not_measured",
        "data_gaps": [
            "Archived target run has no modal_native.bin or Kanaya diagnostic CSVs.",
            "Archived deformation VTUs cover only the final 0.05 s, not the full 0.30-0.50 s window.",
            "Archived log omits monitor point IDs although it records surface (i,j) and coordinates.",
        ],
    }
    write_json(output_dir / "checks.json", checks)
    if not stage1_passed:
        raise RuntimeError("stage-1 amplitude ratio did not reproduce; stop before H1")
    print(json.dumps(checks["stage1"], indent=2))


if __name__ == "__main__":
    main()
