#!/usr/bin/env python3
"""Sweep right-fold body stiffness model files and compare cycle metrics."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILE_KEYS = (
    "leftFrequencyFile", "rightFrequencyFile",
    "leftModeFile", "rightModeFile",
    "leftSurfaceNasFile", "rightSurfaceNasFile",
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep right vocal-fold body stiffness files and run compare_stiffness_metrics."
    )
    parser.add_argument("--base-param", type=Path, default=ROOT / "input" / "param.txt")
    parser.add_argument("--executable", type=Path, default=ROOT / "build" / "simulation")
    parser.add_argument("--working-directory", type=Path, default=ROOT / "build")
    parser.add_argument("--archive-root", type=Path, default=ROOT / "analysis_runs")
    parser.add_argument("--sweep-dir", type=Path)
    parser.add_argument("--b-values", nargs="+", type=int, default=[2, 4, 6, 8, 10, 12])
    parser.add_argument("--c-value", type=int, default=2)
    parser.add_argument("--duration", type=float, default=0.15)
    parser.add_argument("--prominence-fraction", type=float, default=0.05)
    parser.add_argument("--minimum-peaks", type=int, default=6)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rerun-completed", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def resolve(path: Path, base: Path = ROOT) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def named_value(text: str, key: str) -> str:
    match = re.search(rf"(?mi)^\s*{re.escape(key)}\s*=\s*(.*?)\s*$", text)
    if not match or not match.group(1):
        raise ValueError(f"{key} is missing from the base parameter file")
    return match.group(1)


def replace_named_value(text: str, key: str, value: str | Path) -> str:
    pattern = re.compile(rf"(?mi)^(\s*{re.escape(key)}\s*=\s*).*$")
    replaced, count = pattern.subn(lambda match: match.group(1) + str(value), text)
    if count != 1:
        raise ValueError(f"expected exactly one {key} entry, found {count}")
    return replaced


def make_parameter_text(base_text: str, base_dir: Path, body: int, c_value: int) -> tuple[str, Path, Path]:
    text = base_text
    # Generated params live below analysis_runs, so make every model input
    # absolute before moving the parameter file away from input/.
    for key in FILE_KEYS:
        configured = Path(named_value(text, key)).expanduser()
        absolute = configured if configured.is_absolute() else (base_dir / configured).resolve()
        text = replace_named_value(text, key, absolute)

    frequency = (base_dir / f"M5_test/M5_freq_T3_d2_b{body}c{c_value}.txt").resolve()
    mode = (base_dir / f"M5_test/M5_mode_T3_b{body}c{c_value}.vtu").resolve()
    text = replace_named_value(text, "rightFrequencyFile", frequency)
    text = replace_named_value(text, "rightModeFile", mode)
    return text, frequency, mode


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    args = arguments()
    base_param = resolve(args.base_param)
    executable = resolve(args.executable)
    workdir = resolve(args.working_directory)
    archive_root = resolve(args.archive_root)
    if not base_param.is_file():
        raise FileNotFoundError(base_param)
    if not args.dry_run and not executable.is_file():
        raise FileNotFoundError(executable)
    if args.duration <= 0.0 or args.bootstrap_samples < 0:
        raise ValueError("duration must be positive and bootstrap samples non-negative")
    if len(set(args.b_values)) != len(args.b_values):
        raise ValueError("--b-values must not contain duplicates")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep_dir = resolve(args.sweep_dir) if args.sweep_dir else archive_root / f"right_stiffness_sweep_{timestamp}"
    sweep_dir.mkdir(parents=True, exist_ok=True)
    base_text = base_param.read_text(encoding="utf-8")
    environment = os.environ.copy()
    if args.threads is not None:
        environment["OMP_NUM_THREADS"] = str(args.threads)

    write_json(sweep_dir / "sweep_config.json", {
        "base_param": str(base_param), "executable": str(executable),
        "working_directory": str(workdir), "right_body_values": args.b_values,
        "c_value": args.c_value, "analysis_duration_s": args.duration,
        "independent_initial_conditions": True,
    })
    completed: list[tuple[str, Path]] = []
    for body in args.b_values:
        label = f"right_b{body}c{args.c_value}"
        condition_dir = sweep_dir / label
        status_path = condition_dir / "status.json"
        output_dir = condition_dir / "output"
        if status_path.is_file() and args.resume and not args.rerun_completed:
            old = json.loads(status_path.read_text(encoding="utf-8"))
            if old.get("status") == "completed":
                print(f"Skipping completed {label}")
                completed.append((label, output_dir))
                continue

        param_dir = condition_dir / "input"
        param_dir.mkdir(parents=True, exist_ok=True)
        text, frequency, mode = make_parameter_text(
            base_text, base_param.parent, body, args.c_value
        )
        missing = [path for path in (frequency, mode) if not path.is_file()]
        if missing:
            write_json(status_path, {
                "condition": label, "status": "missing_input",
                "missing_files": [str(path) for path in missing],
            })
            message = f"{label}: missing " + ", ".join(str(path) for path in missing)
            if args.fail_fast:
                raise FileNotFoundError(message)
            print(message, file=sys.stderr)
            continue
        generated_param = param_dir / "param.txt"
        generated_param.write_text(text, encoding="utf-8")
        if args.dry_run:
            write_json(status_path, {"condition": label, "status": "pending", "dry_run": True})
            print(f"Prepared {label}: {generated_param}")
            continue

        started = time.monotonic()
        write_json(status_path, {"condition": label, "status": "running"})
        with (condition_dir / "stdout.log").open("w", encoding="utf-8") as stdout, \
             (condition_dir / "stderr.log").open("w", encoding="utf-8") as stderr:
            process = subprocess.run(
                [str(executable), str(generated_param)], cwd=workdir,
                env=environment, stdout=stdout, stderr=stderr, check=False,
            )
        status = "completed" if process.returncode == 0 else "simulation_failed"
        write_json(status_path, {
            "condition": label, "right_body_value": body, "c_value": args.c_value,
            "status": status, "return_code": process.returncode,
            "elapsed_time_s": time.monotonic() - started,
            "right_frequency_file": str(frequency), "right_mode_file": str(mode),
        })
        if process.returncode == 0:
            completed.append((label, output_dir))
        else:
            print(f"{label}: simulation failed; see {condition_dir / 'stderr.log'}", file=sys.stderr)
            if args.fail_fast:
                break

    if not args.dry_run and completed:
        command = [
            sys.executable, str(ROOT / "tools" / "compare_stiffness_metrics.py"),
            *(f"{label}={output}" for label, output in completed),
            "--duration", str(args.duration),
            "--prominence-fraction", str(args.prominence_fraction),
            "--minimum-peaks", str(args.minimum_peaks),
            "--bootstrap-samples", str(args.bootstrap_samples),
            "--seed", str(args.seed),
            "--output-dir", str(sweep_dir),
        ]
        comparison = subprocess.run(command, env=environment, check=False)
        if comparison.returncode != 0:
            raise RuntimeError("compare_stiffness_metrics.py failed")
    print(f"Sweep directory: {sweep_dir}")


if __name__ == "__main__":
    main()
