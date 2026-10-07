#!/usr/bin/env python3
"""Run the two-stage pressure/alignment and initial-gap diagnostic sweep."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from output_layout import output_path, result_path


ROOT = Path(__file__).resolve().parent.parent
INPUT_PATH_KEYS = {
    "leftfrequencyfile", "rightfrequencyfile", "leftmodefile", "rightmodefile",
    "leftsurfacenasfile", "rightsurfacenasfile", "leftsurfacefile", "rightsurfacefile",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-param", type=Path, default=ROOT / "input" / "param.txt")
    parser.add_argument("--executable", type=Path, default=ROOT / "build" / "simulation")
    parser.add_argument("--working-directory", type=Path, default=ROOT / "build")
    parser.add_argument("--sweep-directory", type=Path)
    parser.add_argument("--stage", choices=("all", "pressure", "gap"), default="all")
    parser.add_argument("--pressures", type=float, nargs="+", default=[2900, 3300, 3600, 4000])
    parser.add_argument("--initial-gaps", type=float, nargs="+", default=[0.0, -0.1, -0.2])
    parser.add_argument("--target-current-pg", type=float, default=2900.0)
    parser.add_argument("--selected-pressure", type=float,
                        help="required for --stage gap when no pressure summary is available")
    parser.add_argument("--zeta", type=float, default=0.06)
    parser.add_argument("--pressure-duration", type=float, default=0.30)
    parser.add_argument("--gap-duration", type=float, default=0.50)
    parser.add_argument("--pressure-ramp-time", type=float, default=0.10)
    parser.add_argument("--perturbation-time", type=float, default=0.20)
    parser.add_argument("--perturbation-mode", type=int, default=5)
    parser.add_argument("--perturbation-amplitude-mm", type=float, default=0.001)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--omp-threads", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def normalize_key(value: str) -> str:
    return value.strip().lower().replace("_", "").replace("-", "")


def data_line_indices(lines: list[str]) -> list[int]:
    return [index for index, line in enumerate(lines)
            if line.strip() and not line.lstrip().startswith("#")]


def named_line_map(lines: list[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        result[normalize_key(stripped.split("=", 1)[0])] = index
    return result


def prepare_parameter_text(
    base_path: Path,
    pressure: float,
    zeta: float,
    duration: float,
    initial_gap: float,
    perturbation: bool,
    args: argparse.Namespace,
) -> str:
    lines = base_path.read_text(encoding="utf-8").splitlines()
    data_indices = data_line_indices(lines)
    if len(data_indices) < 23:
        raise ValueError("base parameter file has fewer than 23 positional values")
    dt = float(lines[data_indices[4]].strip())
    nstep = int(round(duration / dt))
    lines[data_indices[2]] = str(nstep)
    lines[data_indices[5]] = f"{zeta:.17g}"
    lines[data_indices[6]] = f"{zeta:.17g}"
    lines[data_indices[15]] = f"{pressure:.17g}"

    named = named_line_map(lines)
    # A copied parameter file resolves relative model paths from its new
    # location, so freeze existing input references to absolute paths.
    for key in INPUT_PATH_KEYS:
        if key not in named:
            continue
        index = named[key]
        lhs, rhs = lines[index].split("=", 1)
        source = Path(rhs.strip())
        if not source.is_absolute():
            source = (base_path.parent / source).resolve()
        lines[index] = f"{lhs.strip()} = {source}"

    values = {
        "initialGapMm": f"{initial_gap:.17g}",
        "pressureRampTimeSec": f"{args.pressure_ramp_time:.17g}",
        "diagnosticOutputIntervalSteps": "5",
        "perturbationEnabled": "1" if perturbation else "0",
        "perturbationTimeSec": f"{args.perturbation_time:.17g}",
        "perturbationModeIndex": str(args.perturbation_mode),
        "perturbationAmplitudeAtProbeMm": f"{args.perturbation_amplitude_mm:.17g}",
        "perturbationPattern": "symmetric_opening",
    }
    named = named_line_map(lines)
    for key, value in values.items():
        normalized = normalize_key(key)
        if normalized in named:
            lines[named[normalized]] = f"{key}={value}"
        else:
            lines.append(f"{key}={value}")
    return "\n".join(lines) + "\n"


def git_commit() -> str:
    process = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    )
    return process.stdout.strip() if process.returncode == 0 else "unknown"


def condition_name(stage: str, pressure: float, gap: float) -> str:
    pressure_text = f"{pressure:g}".replace(".", "p").replace("-", "m")
    gap_text = f"{gap:g}".replace(".", "p").replace("-", "m")
    return f"{stage}_ps_{pressure_text}_gap_{gap_text}"


def run_condition(condition: dict[str, Any], args: argparse.Namespace,
                  sweep_dir: Path, commit: str) -> dict[str, Any]:
    run_dir = sweep_dir / condition_name(
        condition["stage"], condition["input_pressure_pa"], condition["initial_gap_mm"]
    )
    run_dir.mkdir(parents=False, exist_ok=False)
    parameter_path = run_dir / "param.txt"
    parameter_path.write_text(
        prepare_parameter_text(
            args.base_param.resolve(), condition["input_pressure_pa"], args.zeta,
            condition["duration_s"], condition["initial_gap_mm"],
            condition["perturbation_enabled"], args,
        ),
        encoding="utf-8",
    )
    command = [str(args.executable.resolve()), str(parameter_path)]
    started = datetime.now(timezone.utc).isoformat()
    environment = os.environ.copy()
    environment["SIMULATION_RUN_DIR"] = str(run_dir)
    if args.omp_threads is not None:
        environment["OMP_NUM_THREADS"] = str(args.omp_threads)
    with (run_dir / "stdout.log").open("w", encoding="utf-8") as stdout, \
            (run_dir / "stderr.log").open("w", encoding="utf-8") as stderr:
        process = subprocess.run(
            command, cwd=args.working_directory.resolve(), env=environment,
            stdout=stdout, stderr=stderr, check=False,
        )
    manifest = {
        **condition,
        "parameter_file": str(parameter_path),
        "command": command,
        "git_commit": commit,
        "started_at_utc": started,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": process.returncode,
    }
    if process.returncode == 0:
        analysis_start = min(0.15, max(0.0, condition["duration_s"] - 0.10))
        analysis_command = [
            sys.executable, str(ROOT / "tools" / "analyze_self_oscillation.py"), str(run_dir),
            "--analysis-start", str(analysis_start),
            "--analysis-end", str(condition["duration_s"]),
            "--tail-duration", "0.10",
        ]
        analysis = subprocess.run(
            analysis_command, cwd=ROOT, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        manifest["analysis_command"] = analysis_command
        manifest["analysis_exit_code"] = analysis.returncode
        (run_dir / "analysis_stdout.log").write_text(analysis.stdout, encoding="utf-8")
        (run_dir / "analysis_stderr.log").write_text(analysis.stderr, encoding="utf-8")
    (run_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    summary_path = result_path(run_dir, "analysis_summary.json")
    analysis_summary = json.loads(summary_path.read_text(encoding="utf-8")) \
        if summary_path.is_file() else {}
    return {**condition, "run_directory": str(run_dir), "exit_code": process.returncode,
            **analysis_summary}


def execute(conditions: list[dict[str, Any]], args: argparse.Namespace,
            sweep_dir: Path, commit: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_condition, item, args, sweep_dir, commit)
                   for item in conditions]
        for future in concurrent.futures.as_completed(futures):
            try:
                row = future.result()
            except Exception as exception:
                row = {"exit_code": -1, "runner_error": str(exception)}
            rows.append(row)
            print(json.dumps(row, ensure_ascii=False))
    return rows


def write_outputs(sweep_dir: Path, rows: list[dict[str, Any]], target: float) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (sweep_dir / "sweep_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    valid = [row for row in rows if row.get("tail_current_pg_mean_pa") is not None]
    if not valid:
        return
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    pressures = [row["input_pressure_pa"] for row in valid]
    current_pg = [row["tail_current_pg_mean_pa"] for row in valid]
    axes[0].scatter(pressures, current_pg)
    axes[0].axhline(target, color="black", linestyle="--", linewidth=0.8)
    axes[0].set(xlabel="input pressure [Pa]", ylabel="tail mean currentPg [Pa]")
    metric = [row.get("sigma_per_s", float("nan")) for row in valid]
    axes[1].scatter(current_pg, metric)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set(xlabel="tail mean currentPg [Pa]", ylabel="sigma [1/s]")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.savefig(output_path(sweep_dir, "diagnostic_sweep_comparison.png"), dpi=160)
    plt.close(figure)


def main() -> None:
    args = arguments()
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.jobs > 1:
        print("warning: each simulation uses OpenMP; --jobs may oversubscribe the machine",
              file=sys.stderr)
    if not args.base_param.is_file():
        raise FileNotFoundError(args.base_param)
    if not args.dry_run and not args.executable.is_file():
        raise FileNotFoundError(args.executable)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep_dir = (args.sweep_directory.resolve() if args.sweep_directory
                 else ROOT / "diagnostic_runs" / f"diagnostic_sweep_{timestamp}")
    if sweep_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing sweep directory: {sweep_dir}")

    pressure_conditions = [
        {"stage": "pressure", "input_pressure_pa": pressure, "initial_gap_mm": 0.0,
         "duration_s": args.pressure_duration, "perturbation_enabled": False}
        for pressure in args.pressures
    ]
    if args.dry_run:
        for condition in pressure_conditions if args.stage != "gap" else []:
            print(json.dumps(condition, ensure_ascii=False))
        selected = args.selected_pressure
        if args.stage in {"all", "gap"}:
            print(json.dumps({"stage": "gap", "input_pressure_pa": selected or "SELECTED_BY_currentPg",
                              "initial_gaps_mm": args.initial_gaps,
                              "duration_s": args.gap_duration,
                              "perturbation_enabled": True}, ensure_ascii=False))
        return

    sweep_dir.mkdir(parents=True)
    commit = git_commit()
    rows: list[dict[str, Any]] = []
    pressure_rows: list[dict[str, Any]] = []
    if args.stage in {"all", "pressure"}:
        pressure_rows = execute(pressure_conditions, args, sweep_dir, commit)
        rows.extend(pressure_rows)

    selected_pressure = args.selected_pressure
    if args.stage == "all":
        candidates = [row for row in pressure_rows
                      if row.get("exit_code") == 0
                      and row.get("tail_current_pg_mean_pa") is not None]
        if candidates:
            selected_pressure = min(
                candidates,
                key=lambda row: abs(
                    float(row["tail_current_pg_mean_pa"]) - args.target_current_pg
                ),
            )["input_pressure_pa"]
        else:
            print("warning: no successful pressure case; skipping gap stage", file=sys.stderr)
    if args.stage == "gap" and selected_pressure is None:
        raise ValueError("--stage gap requires --selected-pressure")
    if args.stage in {"all", "gap"} and selected_pressure is not None:
        gap_conditions = [
            {"stage": "gap", "input_pressure_pa": selected_pressure,
             "initial_gap_mm": gap, "duration_s": args.gap_duration,
             "perturbation_enabled": True}
            for gap in args.initial_gaps
        ]
        rows.extend(execute(gap_conditions, args, sweep_dir, commit))
    write_outputs(sweep_dir, rows, args.target_current_pg)
    failures = sum(row.get("exit_code", -1) != 0 for row in rows)
    print(f"completed {len(rows)} conditions with {failures} failures: {sweep_dir}")


if __name__ == "__main__":
    main()
