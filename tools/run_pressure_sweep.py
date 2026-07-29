#!/usr/bin/env python3
"""Run independent subglottal-pressure cases and archive every result."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dynamics_common import replace_pressure_parameter, write_json
from left_right_ratios import write_sweep_ratio_outputs


ROOT = Path(__file__).resolve().parent.parent
ARCHIVE_FILES = [
    "manifest.txt",
    "area.dat",
    "airflow_vt.dat",
    "displace.dat",
    "displace_xy.dat",
    "pressure.dat",
    "pressure_vt.dat",
    "modal_contribution.csv",
    "modal_dominant.csv",
    "contact_iteration_debug.csv",
]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run an independent pressure sweep.")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--base-param", type=Path)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--working-directory", type=Path)
    parser.add_argument("--start", type=float)
    parser.add_argument("--stop", type=float)
    parser.add_argument("--step", type=float)
    parser.add_argument("--pressures", nargs="+", type=float)
    parser.add_argument("--analysis-start", type=float)
    parser.add_argument("--analysis-end", type=float)
    parser.add_argument("--last-cycles", type=int)
    parser.add_argument("--left-column", default="uyL_mm")
    parser.add_argument("--right-column", default="uyR_mm")
    parser.add_argument("--peak-prominence-ratio", type=float, default=0.10)
    parser.add_argument("--period-ratio-target", type=float, default=2.0)
    parser.add_argument("--period-ratio-tolerance", type=float, default=0.20)
    parser.add_argument("--amplitude-ratio-target", type=float, default=2.0)
    parser.add_argument("--amplitude-ratio-tolerance", type=float, default=0.20)
    parser.add_argument("--min-peak-count", type=int, default=6)
    parser.add_argument("--min-valid-amplitude", type=float, default=1.0e-4)
    parser.add_argument("--archive-root", type=Path)
    parser.add_argument("--sweep-dir", type=Path, help="existing/new exact sweep directory")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rerun-completed", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def nested(config: dict[str, Any], *keys: str, default: Any = None) -> Any:
    value: Any = config
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def resolve(path: Path | str, base: Path = ROOT) -> Path:
    value = Path(path)
    return value if value.is_absolute() else (base / value).resolve()


def pressure_values(args: argparse.Namespace, config: dict[str, Any]) -> list[float]:
    if args.pressures:
        return args.pressures
    start = args.start if args.start is not None else nested(config, "pressure", "start_pa")
    stop = args.stop if args.stop is not None else nested(config, "pressure", "stop_pa")
    step = args.step if args.step is not None else nested(config, "pressure", "step_pa")
    if start is None or stop is None or step is None or step == 0:
        raise ValueError("provide --pressures or valid --start/--stop/--step")
    if (stop - start) * step < 0:
        raise ValueError("pressure step points away from stop")
    count = int(np.floor((stop - start) / step + 1e-10)) + 1
    return [float(start + index * step) for index in range(count)]


def pressure_name(pressure: float) -> str:
    if pressure.is_integer():
        return f"pressure_{int(pressure):05d}"
    return "pressure_" + f"{pressure:.6f}".rstrip("0").rstrip(".").replace(".", "p")


def write_summary(sweep_dir: Path) -> None:
    rows: list[dict[str, Any]] = []
    bifurcation: list[dict[str, Any]] = []
    for status_path in sorted(sweep_dir.glob("pressure_*/status.json")):
        status = json.loads(status_path.read_text(encoding="utf-8"))
        run_dir = status_path.parent
        metrics_path = run_dir / "metrics.json"
        metrics = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.is_file() else {}
        ratio_path = run_dir / "left_right_ratio_metrics.csv"
        ratio: dict[str, Any] = {}
        if ratio_path.is_file():
            ratio_row = pd.read_csv(ratio_path).iloc[0].to_dict()
            ratio_status = ratio_row.pop("status", "")
            ratio_warning = ratio_row.pop("warning", "")
            for duplicate in ("pressure_pa", "analysis_start_s", "analysis_end_s"):
                ratio_row.pop(duplicate, None)
            ratio = {
                **ratio_row,
                "left_right_status": ratio_status,
                "left_right_warning": ratio_warning,
            }
        rows.append(
            {
                "pressure_pa": status["pressure_pa"],
                "status": status["status"],
                "return_code": status.get("return_code"),
                **metrics,
                "run_directory": str(run_dir),
                "error_message": status.get("error_message", ""),
                **ratio,
            }
        )
        cycles_path = run_dir / "cycle_metrics.csv"
        if cycles_path.is_file():
            with cycles_path.open(newline="", encoding="utf-8") as stream:
                for cycle in csv.DictReader(stream):
                    bifurcation.append({"pressure_pa": status["pressure_pa"], **cycle})
    if rows:
        keys = list(dict.fromkeys(key for row in rows for key in row))
        with (sweep_dir / "sweep_summary.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys)
            writer.writeheader()
            writer.writerows(rows)
    if bifurcation:
        keys = list(dict.fromkeys(key for row in bifurcation for key in row))
        with (sweep_dir / "bifurcation_points.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys)
            writer.writeheader()
            writer.writerows(bifurcation)


def main() -> None:
    args = arguments()
    config: dict[str, Any] = {}
    if args.config:
        config = json.loads(args.config.read_text(encoding="utf-8"))
    base_param = resolve(args.base_param or config.get("base_param_file", "input/param.txt"))
    executable = resolve(args.executable or config.get("executable", "build/simulation"))
    workdir = resolve(args.working_directory or config.get("working_directory", "build"))
    archive_root = resolve(args.archive_root or config.get("archive_root", "analysis_runs"))
    pressures = pressure_values(args, config)
    analysis_start = args.analysis_start
    if analysis_start is None:
        analysis_start = nested(config, "analysis", "start_time_s")
    analysis_end = args.analysis_end
    if analysis_end is None:
        analysis_end = nested(config, "analysis", "end_time_s")
    last_cycles = args.last_cycles
    if last_cycles is None:
        last_cycles = nested(config, "analysis", "last_cycles")
    if not base_param.is_file():
        raise FileNotFoundError(base_param)
    if not args.dry_run and not executable.is_file():
        raise FileNotFoundError(executable)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep_dir = resolve(args.sweep_dir) if args.sweep_dir else archive_root / f"pressure_sweep_{timestamp}"
    sweep_dir.mkdir(parents=True, exist_ok=True)
    effective_config = {
        "executable": str(executable),
        "base_param_file": str(base_param),
        "working_directory": str(workdir),
        "pressures_pa": pressures,
        "analysis_start_s": analysis_start,
        "analysis_end_s": analysis_end,
        "last_cycles": last_cycles,
        "left_column": args.left_column,
        "right_column": args.right_column,
        "peak_prominence_ratio": args.peak_prominence_ratio,
        "period_ratio_target": args.period_ratio_target,
        "period_ratio_tolerance": args.period_ratio_tolerance,
        "amplitude_ratio_target": args.amplitude_ratio_target,
        "amplitude_ratio_tolerance": args.amplitude_ratio_tolerance,
        "min_peak_count": args.min_peak_count,
        "min_valid_amplitude_mm": args.min_valid_amplitude,
        "independent_initial_conditions": True,
    }
    write_json(sweep_dir / "sweep_config.json", effective_config)
    base_text = base_param.read_text(encoding="utf-8")
    environment = os.environ.copy()
    configured_environment = config.get("environment", {})
    environment.update({str(key): str(value) for key, value in configured_environment.items()})
    if args.threads:
        environment["OMP_NUM_THREADS"] = str(args.threads)

    for pressure in pressures:
        run_dir = sweep_dir / pressure_name(pressure)
        status_path = run_dir / "status.json"
        if status_path.is_file() and args.resume and not args.rerun_completed:
            old = json.loads(status_path.read_text(encoding="utf-8"))
            if old.get("status") in {"completed", "insufficient_cycles"}:
                print(f"Skipping completed {pressure:g} Pa")
                continue
        run_dir.mkdir(parents=True, exist_ok=True)
        param_dir = run_dir / "input"
        param_dir.mkdir(exist_ok=True)
        generated_param = param_dir / "param.txt"
        generated_param.write_text(
            replace_pressure_parameter(base_text, pressure), encoding="utf-8"
        )
        if args.dry_run:
            write_json(status_path, {"pressure_pa": pressure, "status": "pending", "dry_run": True})
            print(f"Dry run prepared {pressure:g} Pa: {generated_param}")
            continue

        output_dir = run_dir / "output"
        started = time.monotonic()
        write_json(status_path, {"pressure_pa": pressure, "status": "running"})
        with (run_dir / "stdout.log").open("w", encoding="utf-8") as stdout, (
            run_dir / "stderr.log"
        ).open("w", encoding="utf-8") as stderr:
            process = subprocess.run(
                [str(executable), str(generated_param)],
                cwd=workdir,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        if process.returncode != 0:
            status = {
                "pressure_pa": pressure,
                "status": "simulation_failed",
                "return_code": process.returncode,
                "elapsed_time_s": time.monotonic() - started,
                "error_message": "simulation returned a non-zero status",
            }
            write_json(status_path, status)
            if args.fail_fast:
                break
            continue
        for name in ARCHIVE_FILES:
            source = output_dir / name
            if source.is_file():
                shutil.copy2(source, run_dir / name)
        shutil.copy2(generated_param, run_dir / "params_used.txt")
        analysis_command = [
            sys.executable,
            str(ROOT / "tools" / "analyze_dynamics.py"),
            "--run-dir",
            str(run_dir),
        ]
        if analysis_start is not None:
            analysis_command += ["--start-time", str(analysis_start)]
        if analysis_end is not None:
            analysis_command += ["--end-time", str(analysis_end)]
        if last_cycles is not None:
            analysis_command += ["--last-cycles", str(last_cycles)]
        analysis_command += [
            "--left-column",
            args.left_column,
            "--right-column",
            args.right_column,
            "--peak-prominence-ratio",
            str(args.peak_prominence_ratio),
            "--period-ratio-target",
            str(args.period_ratio_target),
            "--period-ratio-tolerance",
            str(args.period_ratio_tolerance),
            "--amplitude-ratio-target",
            str(args.amplitude_ratio_target),
            "--amplitude-ratio-tolerance",
            str(args.amplitude_ratio_tolerance),
            "--min-peak-count",
            str(args.min_peak_count),
            "--min-valid-amplitude",
            str(args.min_valid_amplitude),
        ]
        analysis = subprocess.run(analysis_command, env=environment, check=False)
        metrics = {}
        if (run_dir / "metrics.json").is_file():
            metrics = json.loads((run_dir / "metrics.json").read_text(encoding="utf-8"))
        final_status = metrics.get("status", "completed") if analysis.returncode == 0 else "analysis_failed"
        write_json(
            status_path,
            {
                "pressure_pa": pressure,
                "status": final_status,
                "return_code": analysis.returncode,
                "elapsed_time_s": time.monotonic() - started,
                "error_message": "" if analysis.returncode == 0 else "analysis failed",
            },
        )
        write_summary(sweep_dir)
        write_sweep_ratio_outputs(sweep_dir)
        if analysis.returncode != 0 and args.fail_fast:
            break
    write_summary(sweep_dir)
    write_sweep_ratio_outputs(sweep_dir)
    print(f"Sweep directory: {sweep_dir}")


if __name__ == "__main__":
    main()
