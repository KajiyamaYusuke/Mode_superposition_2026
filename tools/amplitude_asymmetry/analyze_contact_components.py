#!/usr/bin/env python3
"""Analyze adopted contact spring/damping loads without modifying a run.

The component binary is paired with the legacy seven-field modal binary from
the same run.  Work is therefore evaluated with exactly the legacy ``delta_q``
field and is reported as algebraic reference work, not irreversible loss.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np


WINDOWS = (
    (0.10, 0.20),
    (0.20, 0.30),
    (0.30, 0.40),
    (0.40, 0.50),
    (0.30, 0.50),
)
SIDES = ("L", "R")
COMPONENTS = ("spring", "damp")
BASELINE_PRIMARY_MJ = {"L": -0.15263216245777325,
                       "R": -1.2055820694547033}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8")


def window_label(start: float, end: float) -> str:
    return f"[{start:.2f},{end:.2f})"


def interval_mask(count: int, dt: float, start: float, end: float) -> np.ndarray:
    """Select observer intervals by their recorded left-end time."""
    time = np.arange(count, dtype=np.float64) * dt
    return (time >= start) & (time < end)


def parse_key_values(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip().strip('"')
    return result


def load_frequency(path: Path, count: int) -> np.ndarray:
    data = np.loadtxt(path, comments="%", ndmin=2)
    if len(data) < count:
        raise ValueError(f"{path}: only {len(data)} modes, expected {count}")
    return np.asarray(data[:count, 0], dtype=np.float64)


def exact_tree_hashes(root: Path, relative_dir: str) -> dict[str, str]:
    directory = root / relative_dir
    return {str(path.relative_to(directory)): sha256(path)
            for path in sorted(directory.glob("*")) if path.is_file()}


def summarize_work(modal: np.ndarray, component: np.ndarray, dt: float,
                   frequencies: list[np.ndarray]) -> tuple[list[dict], dict]:
    count, _, nmode, _ = modal.shape
    if component.shape != (count, 2, nmode, 2):
        raise ValueError(f"component shape mismatch: {component.shape}")
    delta_q = np.asarray(modal[..., 6])
    component_increment = np.asarray(component) * delta_q[..., None]
    reference_increment = np.asarray(modal[..., 3]) * delta_q
    decomposition = component_increment.sum(axis=-1) - reference_increment
    rows: list[dict] = []
    for start, end in WINDOWS:
        mask = interval_mask(count, dt, start, end)
        if not np.any(mask):
            raise ValueError(f"empty analysis window {start}:{end}")
        for side_index, side in enumerate(SIDES):
            for mode_index in range(nmode):
                increments = component_increment[mask, side_index, mode_index, :]
                reference = reference_increment[mask, side_index, mode_index]
                values = increments.sum(axis=0)
                row = {
                    "window": window_label(start, end),
                    "window_start_s": start,
                    "window_end_s": end,
                    "side": side,
                    "mode": mode_index + 1,
                    "source_mode": mode_index + 1,
                    "frequency_hz": float(frequencies[side_index][mode_index]),
                    "interval_count": int(mask.sum()),
                    "W_spring_J": float(values[0]),
                    "W_damp_J": float(values[1]),
                    "W_total_J": float(values.sum()),
                    "W_reference_J": float(reference.sum()),
                    "residual_J": float(values.sum() - reference.sum()),
                    "W_spring_positive_J": float(increments[:, 0].clip(min=0).sum()),
                    "W_spring_negative_J": float(increments[:, 0].clip(max=0).sum()),
                    "W_damp_positive_J": float(increments[:, 1].clip(min=0).sum()),
                    "W_damp_negative_J": float(increments[:, 1].clip(max=0).sum()),
                    "quality": "PASS",
                }
                rows.append(row)
            selected = [row for row in rows
                        if row["window"] == window_label(start, end)
                        and row["side"] == side and row["mode"] != "ALL"]
            total = {key: sum(float(row[key]) for row in selected)
                     for key in ("W_spring_J", "W_damp_J", "W_total_J",
                                 "W_reference_J", "residual_J",
                                 "W_spring_positive_J", "W_spring_negative_J",
                                 "W_damp_positive_J", "W_damp_negative_J")}
            rows.append({
                "window": window_label(start, end),
                "window_start_s": start,
                "window_end_s": end,
                "side": side,
                "mode": "ALL",
                "source_mode": "ALL",
                "frequency_hz": "",
                "interval_count": int(mask.sum()),
                **total,
                "quality": "PASS",
            })
    checks = {
        "max_abs_interval_mode_decomposition_residual_J":
            float(np.max(np.abs(decomposition))),
        "all_component_values_finite": bool(np.all(np.isfinite(component))),
        "all_reference_values_finite": bool(np.all(np.isfinite(reference_increment))),
    }
    return rows, checks


def load_regions(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                         encoding="utf-8")
    if data.ndim == 0:
        data = np.asarray([data], dtype=data.dtype)
    return {name: np.asarray(data[name]) for name in data.dtype.names or ()}


def summarize_regions(regions: dict[str, np.ndarray], dt: float) -> list[dict]:
    rows: list[dict] = []
    for start, end in WINDOWS:
        window = (regions["time_s"] >= start) & (regions["time_s"] < end)
        for side in SIDES:
            for ap_region in ("end", "interior", "ALL"):
                mask = window & (regions["side"] == side)
                if ap_region != "ALL":
                    mask &= regions["ap_region"] == ap_region
                pairs = np.asarray(regions["pairs"][mask], dtype=np.float64)
                pair_total = float(pairs.sum())
                def weighted_mean(field: str) -> float:
                    values = np.asarray(regions[field][mask], dtype=np.float64)
                    return float(np.sum(values * pairs) / pair_total) if pair_total else 0.0
                def finite_max(field: str) -> float:
                    values = np.asarray(regions[field][mask], dtype=np.float64)
                    values = values[np.isfinite(values)]
                    return float(np.max(values)) if len(values) else 0.0
                rows.append({
                    "window": window_label(start, end),
                    "window_start_s": start,
                    "window_end_s": end,
                    "side": side,
                    "ap_region": ap_region,
                    "row_count": int(mask.sum()),
                    "pair_samples": int(pair_total),
                    "excluded_pair_samples": int(np.sum(regions["excluded_pairs"][mask])),
                    "spring_pressure_pair_mean_Pa": weighted_mean("spring_pressure_mean_Pa"),
                    "spring_pressure_max_Pa": finite_max("spring_pressure_max_Pa"),
                    "damp_pressure_pair_mean_Pa": weighted_mean("damp_pressure_mean_Pa"),
                    "damp_pressure_max_Pa": finite_max("damp_pressure_max_Pa"),
                    "spring_work_nodal_J": float(np.sum(regions["spring_work_nodal_J"][mask])),
                    "damp_work_nodal_J": float(np.sum(regions["damp_work_nodal_J"][mask])),
                    "max_penetration_m": finite_max("max_penetration_m"),
                    "max_penetration_rate_m_s": finite_max("max_penetration_rate_m_s"),
                    "damp_pair_power_candidate_integral_J": float(
                        np.sum(regions["damp_pair_power_candidate_W"][mask]) * dt),
                })
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing empty CSV {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def primary_totals(rows: list[dict]) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for side in SIDES:
        row = next(item for item in rows if item["window"] == "[0.30,0.50)"
                   and item["side"] == side and item["mode"] == "ALL")
        result[side] = {key: float(row[key]) for key in
                        ("W_spring_J", "W_damp_J", "W_total_J", "W_reference_J")}
    result["both"] = {key: sum(result[side][key] for side in SIDES)
                      for key in result["L"]}
    return result


def make_figures(output: Path, modal: np.ndarray, component: np.ndarray,
                 dt: float, work_rows: list[dict], region_rows: list[dict]) -> None:
    figures = output / "figures"
    figures.mkdir()
    primary = primary_totals(work_rows)
    x = np.arange(3)
    spring = [primary[key]["W_spring_J"] * 1e3 for key in ("L", "R", "both")]
    damp = [primary[key]["W_damp_J"] * 1e3 for key in ("L", "R", "both")]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x, spring, label="spring")
    ax.bar(x, damp, bottom=spring, label="velocity-dependent")
    ax.set(xticks=x, xticklabels=("L", "R", "L+R"), ylabel="reference work [mJ]",
           title="Adopted contact components, [0.30, 0.50) s")
    ax.axhline(0, color="black", linewidth=.7); ax.legend(); fig.tight_layout()
    fig.savefig(figures / "01_component_work.png", dpi=160); plt.close(fig)

    increments = np.asarray(component) * np.asarray(modal[..., 6])[..., None]
    time = np.arange(len(increments)) * dt
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for side_index, side in enumerate(SIDES):
        cumulative = np.cumsum(increments[:, side_index].sum(axis=1), axis=0) * 1e3
        axes[side_index].plot(time, cumulative[:, 0], label="spring")
        axes[side_index].plot(time, cumulative[:, 1], label="velocity-dependent")
        axes[side_index].plot(time, cumulative.sum(axis=1), label="total", linewidth=1)
        axes[side_index].set_ylabel(f"{side} cumulative [mJ]"); axes[side_index].legend()
    axes[-1].set_xlabel("interval left-end time [s]"); fig.tight_layout()
    fig.savefig(figures / "02_cumulative_work.png", dpi=160); plt.close(fig)

    selected = [row for row in work_rows if row["window"] == "[0.30,0.50)"
                and row["mode"] != "ALL"]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    modes = np.arange(1, modal.shape[2] + 1)
    for ax, side in zip(axes, SIDES):
        side_rows = [row for row in selected if row["side"] == side]
        ax.bar(modes - .2, [row["W_spring_J"] * 1e3 for row in side_rows], .4,
               label="spring")
        ax.bar(modes + .2, [row["W_damp_J"] * 1e3 for row in side_rows], .4,
               label="velocity-dependent")
        ax.set_ylabel(f"{side} [mJ]"); ax.legend()
    axes[-1].set(xlabel="active/source mode", xticks=modes); fig.tight_layout()
    fig.savefig(figures / "03_modal_work.png", dpi=160); plt.close(fig)

    selected_regions = [row for row in region_rows if row["window"] == "[0.30,0.50)"
                        and row["ap_region"] != "ALL"]
    labels = [f"{side}-{region}" for side in SIDES for region in ("end", "interior")]
    lookup = {(row["side"], row["ap_region"]): row for row in selected_regions}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(np.arange(4)-.2, [lookup[tuple(label.split("-"))]["spring_work_nodal_J"]*1e3 for label in labels], .4, label="spring")
    axes[0].bar(np.arange(4)+.2, [lookup[tuple(label.split("-"))]["damp_work_nodal_J"]*1e3 for label in labels], .4, label="velocity-dependent")
    axes[0].set(ylabel="nodal reference work [mJ]", xticks=np.arange(4), xticklabels=labels); axes[0].legend()
    axes[1].bar(np.arange(4)-.2, [lookup[tuple(label.split("-"))]["spring_pressure_max_Pa"] for label in labels], .4, label="spring")
    axes[1].bar(np.arange(4)+.2, [lookup[tuple(label.split("-"))]["damp_pressure_max_Pa"] for label in labels], .4, label="velocity-dependent")
    axes[1].set(ylabel="peak pressure [Pa]", xticks=np.arange(4), xticklabels=labels); axes[1].legend()
    fig.tight_layout(); fig.savefig(figures / "04_ap_regions.png", dpi=160); plt.close(fig)


def classify(primary: dict[str, dict[str, float]]) -> tuple[str, str]:
    spring = abs(primary["R"]["W_spring_J"])
    damp = abs(primary["R"]["W_damp_J"])
    if max(spring, damp) == 0:
        return "Hstage", "right-side component work is identically zero"
    fraction = spring / (spring + damp)
    if fraction >= 0.75:
        return "Hspring", f"spring absolute-share proxy={fraction:.3f}"
    if fraction <= 0.25:
        return "Hdamp", f"velocity-dependent absolute-share proxy={1-fraction:.3f}"
    return "Hmixed", f"spring absolute-share proxy={fraction:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-run", type=Path, required=True,
                        help="archived original output directory")
    parser.add_argument("--previous-analysis", type=Path, required=True)
    parser.add_argument("--production-run", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--smoke-root", type=Path, required=True)
    parser.add_argument("--solver-binary", type=Path, required=True)
    args = parser.parse_args()
    if args.artifact_dir.exists():
        raise FileExistsError(f"refusing to overwrite {args.artifact_dir}")

    diagnostic = args.production_run / "kanaya_diagnostic"
    metadata = json.loads((diagnostic / "run_metadata.json").read_text())
    component_metadata = json.loads((diagnostic / "contact_component_metadata.json").read_text())
    if not metadata.get("complete") or not component_metadata.get("complete"):
        raise ValueError("production diagnostic is incomplete")
    count = int(metadata["integrated_intervals"])
    nmode = int(metadata["nmode"][0])
    dt = float(metadata["dt_s"])
    modal_path = diagnostic / "modal_native.bin"
    component_path = diagnostic / "contact_component_native.bin"
    if modal_path.stat().st_size != count * 2 * nmode * 7 * 8:
        raise ValueError("legacy modal binary size mismatch")
    if component_path.stat().st_size != count * 2 * nmode * 2 * 8:
        raise ValueError("component binary size mismatch")
    modal = np.memmap(modal_path, dtype="<f8", mode="r", shape=(count, 2, nmode, 7))
    component = np.memmap(component_path, dtype="<f8", mode="r", shape=(count, 2, nmode, 2))
    manifest_text = parse_key_values(args.production_run / "txt/manifest.txt")
    frequency_paths = [Path(manifest_text[f"{name}_frequency"])
                       for name in ("left", "right")]
    frequencies = [load_frequency(path, nmode) for path in frequency_paths]
    work_rows, work_checks = summarize_work(modal, component, dt, frequencies)
    regions = load_regions(diagnostic / "contact_component_regions.csv")
    region_rows = summarize_regions(regions, dt)
    primary = primary_totals(work_rows)

    previous_production = args.previous_analysis / "production"
    previous_modal = previous_production / "kanaya_diagnostic/modal_native.bin"
    source_dat = exact_tree_hashes(args.source_run, "dat")
    current_dat = exact_tree_hashes(args.production_run, "dat")
    smoke_hashes = {name: {
        "state": sha256(args.smoke_root / name / "regression_state.bin"),
        "dat": exact_tree_hashes(args.smoke_root / name, "dat"),
        "result": exact_tree_hashes(args.smoke_root / name, "result"),
    } for name in ("off", "component_disabled", "component_enabled")}
    smoke_exact = smoke_hashes["off"] == smoke_hashes["component_disabled"] == smoke_hashes["component_enabled"]
    residual_keys = [key for key in component_metadata
                     if key.startswith("max_") and "_residual_" in key]
    baseline_error = {side: primary[side]["W_reference_J"] * 1e3 - BASELINE_PRIMARY_MJ[side]
                      for side in SIDES}
    region_primary = {(row["side"], row["ap_region"]): row for row in region_rows
                      if row["window"] == "[0.30,0.50)" and row["ap_region"] == "ALL"}
    region_modal_residual = {side: (
        region_primary[(side, "ALL")]["spring_work_nodal_J"]
        + region_primary[(side, "ALL")]["damp_work_nodal_J"]
        - primary[side]["W_reference_J"]
    ) for side in SIDES}
    checks = {
        "overall_qc": "PASS",
        "production_complete": True,
        "trajectory_identical_all_dat": source_dat == current_dat,
        "trajectory_dat_source": source_dat,
        "trajectory_dat_production": current_dat,
        "legacy_modal_binary_identical_to_previous_diagnostic":
            sha256(modal_path) == sha256(previous_modal),
        "legacy_modal_binary_sha256": sha256(modal_path),
        "previous_legacy_modal_binary_sha256": sha256(previous_modal),
        "smoke_off_component_disabled_component_enabled_exact": smoke_exact,
        "smoke_hashes": smoke_hashes,
        "component_metadata_max_residuals": {key: component_metadata[key]
                                             for key in residual_keys},
        "work": work_checks,
        "primary_reference_baseline_error_mJ": baseline_error,
        "primary_region_nodal_minus_modal_total_J": region_modal_residual,
        "component_binary_expected_bytes": count * 2 * nmode * 2 * 8,
        "component_binary_actual_bytes": component_path.stat().st_size,
        "symmetric_12_12_components": "NOT_MEASURED",
    }
    pass_conditions = [
        checks["trajectory_identical_all_dat"],
        checks["legacy_modal_binary_identical_to_previous_diagnostic"],
        smoke_exact,
        work_checks["all_component_values_finite"],
        work_checks["all_reference_values_finite"],
        work_checks["max_abs_interval_mode_decomposition_residual_J"] <= 1e-12,
        max(abs(value) for value in baseline_error.values()) <= 1e-9,
        max(abs(value) for value in region_modal_residual.values()) <= 1e-12,
        all(float(component_metadata[key]) <= 1e-8 for key in residual_keys),
    ]
    if not all(pass_conditions):
        checks["overall_qc"] = "FAILED_QC"

    args.artifact_dir.mkdir(parents=True)
    write_csv(args.artifact_dir / "contact_component_summary.csv", work_rows)
    write_csv(args.artifact_dir / "contact_component_region_summary.csv", region_rows)
    write_json(args.artifact_dir / "checks_contact_components.json", checks)
    try:
        git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
        git_status = subprocess.check_output(["git", "status", "--short"], text=True).splitlines()
    except (OSError, subprocess.CalledProcessError):
        git_sha, git_status = "UNKNOWN", []
    input_paths = {key: Path(value) for key, value in manifest_text.items()
                   if key in {"left_mode_vtu", "right_mode_vtu", "left_frequency",
                              "right_frequency", "left_surface_nas", "right_surface_nas"}}
    windows = [{"label": window_label(start, end), "start_s": start, "end_s": end,
                "interval_count": int(interval_mask(count, dt, start, end).sum())}
               for start, end in WINDOWS]
    analysis_manifest = {
        "schema": "contact-component-analysis-v1",
        "git_head": git_sha,
        "git_status_short": git_status,
        "source_run": str(args.source_run.resolve()),
        "previous_analysis": str(args.previous_analysis.resolve()),
        "production_run": str(args.production_run.resolve()),
        "solver_binary": {"path": str(args.solver_binary.resolve()),
                          "sha256": sha256(args.solver_binary)},
        "diagnostic_binaries": {
            "modal_native": {"sha256": sha256(modal_path), "shape": [count, 2, nmode, 7]},
            "contact_component_native": {"sha256": sha256(component_path),
                                         "shape": [count, 2, nmode, 2]},
        },
        "input_files": {key: {"path": str(path), "sha256": sha256(path)}
                        for key, path in input_paths.items()},
        "parameter_file": {"path": manifest_text.get("parameter_file"),
                           "sha256": sha256(Path(manifest_text["parameter_file"]))},
        "conditions": {key: manifest_text.get(key) for key in
                       ("nstep", "dt_s", "nmode", "zetaL", "zetaR", "kc1", "kc2",
                        "kc3", "contact_reference_frequency_hz", "pressure_ramp_time_sec",
                        "pressure_schedule_final_pa", "target_initial_minimum_gap_mm")},
        "runtime": {"OMP_NUM_THREADS": "4", "OMP_PROC_BIND": "close",
                    "OMP_PLACES": "cores", "SIMULATION_KANAYA_DIAGNOSTICS": "1",
                    "SIMULATION_CONTACT_COMPONENTS": "1", "SIMULATION_ACOUSTIC_TRACE": "1"},
        "windows": windows,
        "interval_convention": "left-end time; final solver step excluded (nstep-1 intervals)",
        "contact_stage": component_metadata["contact_stage"],
        "work_definition": component_metadata["work_definition"],
        "primary_totals_J": primary,
        "symmetric_12_12_reference": {"L_mJ": -0.95918, "R_mJ": -0.95918,
                                      "components": "NOT_MEASURED"},
        "qc": checks["overall_qc"],
        "artifacts": ["contact_component_summary.csv", "contact_component_region_summary.csv",
                      "checks_contact_components.json", "manifest.json", "report.md", "figures/"],
    }
    write_json(args.artifact_dir / "manifest.json", analysis_manifest)
    make_figures(args.artifact_dir, modal, component, dt, work_rows, region_rows)

    hypothesis, reason = classify(primary) if checks["overall_qc"] == "PASS" else ("Hstage", "QC failed")
    primary_mode_rows = [row for row in work_rows if row["window"] == "[0.30,0.50)"
                         and row["mode"] != "ALL"]
    top = sorted(primary_mode_rows, key=lambda row: abs(row["W_total_J"]), reverse=True)[:8]
    def format_line(key: str) -> str:
        values = primary[key]
        damp_fraction = values["W_damp_J"] / values["W_total_J"] if values["W_total_J"] else float("nan")
        return (f"| {key} | {values['W_spring_J']*1e3:.9g} | "
                f"{values['W_damp_J']*1e3:.9g} | {values['W_total_J']*1e3:.9g} | "
                f"{damp_fraction:.6g} |")
    trajectory_statement = (
        "既存9 datの軌道と旧7-field modal binaryはともにバイト一致した。"
        if checks["trajectory_identical_all_dat"]
        and checks["legacy_modal_binary_identical_to_previous_diagnostic"]
        else "既存runとのbit一致ゲートは不合格であり、旧runそのものの成分分解とは認定しない。"
    )
    report = f"""# 12/6・900 Pa 接触ばね／速度依存項の仕事分解

QCは **{checks['overall_qc']}**。{trajectory_statement}
主窓 `[0.30,0.50) s` の右総接触参考仕事は {primary['R']['W_total_J']*1e3:.9g} mJで、ばね {primary['R']['W_spring_J']*1e3:.9g} mJ、速度依存 {primary['R']['W_damp_J']*1e3:.9g} mJに分解された。
QC不合格時の分類は **{hypothesis}**（{reason}）として機構結論を保留する。これは採用予測反復の荷重による代数的参照仕事であり、不可逆散逸量や振幅低下の因果を確定するものではない。
対称12/12の総接触参考仕事は既報 L=R=-0.95918 mJだが、その成分別内訳は **NOT_MEASURED** である。

## 主結果

| side | spring [mJ] | velocity-dependent [mJ] | total [mJ] | damp/total |
|---|---:|---:|---:|---:|
{format_line('L')}
{format_line('R')}
{format_line('both')}

正負の割合は代数和を分母にするため100%を超え得る。詳細は `contact_component_summary.csv` に正・負増分を別記した。

## 主なモード（|total|上位）

| side | mode | frequency [Hz] | spring [mJ] | velocity-dependent [mJ] | total [mJ] |
|---|---:|---:|---:|---:|---:|
"""
    for row in top:
        report += (f"| {row['side']} | {row['mode']} | {row['frequency_hz']:.7g} | "
                   f"{row['W_spring_J']*1e3:.7g} | {row['W_damp_J']*1e3:.7g} | "
                   f"{row['W_total_J']*1e3:.7g} |\n")
    report += f"""

指定モード4, 5, 7, 9, 12, 14を含む全モード値はCSVに収録した。モード基底は完全直交ではないため、各モード値を独立した物理エネルギー損失とは解釈しない。

## AP・時刻分布と整合性

AP端／中央の節点仕事、接触pair数、ばね／速度依存圧、最大貫入量、閉鎖速度候補、局所速度依存power候補は `contact_component_region_summary.csv`、全時系列はproductionの `contact_component_regions.csv` にある。4枚の図は成分総量、累積時系列、モード別、AP別を示す。

- 全interval・side・modeの最大分解残差: {work_checks['max_abs_interval_mode_decomposition_residual_J']:.3e} J
- LのAP節点和−modal総仕事: {region_modal_residual['L']:.3e} J
- RのAP節点和−modal総仕事: {region_modal_residual['R']:.3e} J
- observer内部の圧力・節点力・モード力・仕事・作用反作用検査: metadata記録値を `checks_contact_components.json` に転記

## 観測・機構候補・限界

観測した新runの代数値は上表のとおりだが、`FAILED_QC` の場合は Hspring/Hdamp/Hmixed の判定を行わない。速度依存項の局所pair power積分は補助量で、最終Newmark変位増分を使う本体の `Q Δq` と同一ではない。接触面積・法線・対応点・採用予測反復が時間変化するため、ばね仕事も単純な保存ポテンシャル差ではない。受動観測1条件なので「接触成分が右振幅低下を引き起こした」とは断定しない。

次の最小検証案は、今回の1 runを再利用して採用反復と最終変位増分の位相差を接触イベントごとに可視化することであり、新たな物理アブレーションは別タスクとする。
"""
    (args.artifact_dir / "report.md").write_text(report, encoding="utf-8")
    if checks["overall_qc"] != "PASS":
        raise RuntimeError("analysis artifacts written, but QC failed")


if __name__ == "__main__":
    main()
