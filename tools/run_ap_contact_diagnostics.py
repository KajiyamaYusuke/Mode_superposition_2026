#!/usr/bin/env python3
"""Run isolated A/B AP-contact and live/frozen-load diagnostic cases."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PATH_KEYS = {
    "leftfrequencyfile", "rightfrequencyfile", "leftmodefile", "rightmodefile",
    "leftsurfacenasfile", "rightsurfacenasfile", "leftsurfacefile", "rightsurfacefile",
    "fixednodeidsfile",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-param", type=Path, default=ROOT / "input" / "param.txt")
    parser.add_argument("--executable", type=Path, default=ROOT / "build" / "simulation")
    parser.add_argument("--working-directory", type=Path, default=ROOT / "build")
    parser.add_argument("--run-directory", type=Path)
    parser.add_argument("--stage", choices=("ab", "freeze", "all"), default="ab")
    parser.add_argument("--duration", type=float, default=0.70)
    parser.add_argument("--reference-window", type=float, nargs=2, default=(0.25, 0.30))
    parser.add_argument("--perturbation-time", type=float, default=0.35)
    parser.add_argument("--excluded-rows", type=int, default=3)
    parser.add_argument("--detail-window", type=float, nargs=2, default=(0.349, 0.351))
    parser.add_argument("--omp-threads", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def normalize(value: str) -> str:
    return value.strip().lower().replace("_", "").replace("-", "")


def parameter_text(base: Path, values: dict[str, str], duration: float) -> str:
    lines = base.read_text(encoding="utf-8").splitlines()
    data = [index for index, line in enumerate(lines)
            if line.strip() and not line.lstrip().startswith("#")]
    if len(data) < 23:
        raise ValueError("base parameter file has fewer than 23 positional values")
    dt = float(lines[data[4]].strip())
    lines[data[2]] = str(int(round(duration / dt)))
    named: dict[str, int] = {}
    for index, line in enumerate(lines):
        if line.strip() and not line.lstrip().startswith("#") and "=" in line:
            named[normalize(line.split("=", 1)[0])] = index
    for key in PATH_KEYS:
        if key not in named:
            continue
        index = named[key]
        lhs, rhs = lines[index].split("=", 1)
        path = Path(rhs.strip())
        if not path.is_absolute():
            path = (base.parent / path).resolve()
        lines[index] = f"{lhs.strip()} = {path}"
    named = {normalize(line.split("=", 1)[0]): index for index, line in enumerate(lines)
             if line.strip() and not line.lstrip().startswith("#") and "=" in line}
    for key, value in values.items():
        normalized = normalize(key)
        if normalized in named:
            lines[named[normalized]] = f"{key} = {value}"
        else:
            lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def cases(args: argparse.Namespace) -> list[tuple[str, dict[str, str]]]:
    common = {
        "apContactDiagnosticEnabled": "1",
        "apContactAnalysisRowsPerEnd": str(args.excluded_rows),
        "apContactAnalysisDistanceMm": "0",
        "apContactExcludedDistanceMm": "0",
        "contactDetailStartSec": f"{args.detail_window[0]:.17g}",
        "contactDetailEndSec": f"{args.detail_window[1]:.17g}",
        "diagnosticReferenceWindowStartSec": f"{args.reference_window[0]:.17g}",
        "diagnosticReferenceWindowEndSec": f"{args.reference_window[1]:.17g}",
        "perturbationEnabled": "1",
        "perturbationTimeSec": f"{args.perturbation_time:.17g}",
        "diagnosticOutputIntervalSteps": "5",
    }
    result: list[tuple[str, dict[str, str]]] = []
    if args.stage in {"ab", "all"}:
        result.extend([
            ("A_live_contact", {**common, "apContactExcludedRowsPerEnd": "0",
                                "diagnosticLoadMode": "live"}),
            ("B_exclude_ap_ends", {**common,
                                   "apContactExcludedRowsPerEnd": str(args.excluded_rows),
                                   "diagnosticLoadMode": "live"}),
        ])
    if args.stage in {"freeze", "all"}:
        for mode in ("live", "freeze_contact", "freeze_fluid", "freeze_all"):
            result.append((f"A_{mode}", {**common, "apContactExcludedRowsPerEnd": "0",
                                         "diagnosticLoadMode": mode}))
    # De-duplicate A_live when --stage all.
    return list(dict(result).items())


def git_commit() -> str:
    process = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False)
    return process.stdout.strip() if process.returncode == 0 else "unknown"


def main() -> None:
    args = arguments()
    if not args.reference_window[0] < args.reference_window[1] <= args.perturbation_time:
        raise ValueError("reference window must end no later than perturbation time")
    if args.excluded_rows < 0:
        raise ValueError("--excluded-rows must be >= 0")
    selected = cases(args)
    if args.dry_run:
        print(json.dumps({"duration_s": args.duration, "cases": [name for name, _ in selected]},
                         indent=2, ensure_ascii=False))
        return
    if not args.base_param.is_file() or not args.executable.is_file():
        raise FileNotFoundError("base parameter file or simulation executable is missing")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = (args.run_directory.resolve() if args.run_directory
            else ROOT / "diagnostic_runs" / f"ap_contact_{stamp}")
    if root.exists():
        raise FileExistsError(f"refusing to overwrite {root}")
    root.mkdir(parents=True)
    commit = git_commit()
    summary: list[dict[str, object]] = []
    for name, values in selected:
        run = root / name
        run.mkdir()
        param = run / "param.txt"
        param.write_text(parameter_text(args.base_param.resolve(), values, args.duration),
                         encoding="utf-8")
        command = [str(args.executable.resolve()), str(param)]
        environment = os.environ.copy()
        environment["SIMULATION_RUN_DIR"] = str(run)
        if args.omp_threads is not None:
            environment["OMP_NUM_THREADS"] = str(args.omp_threads)
        started = datetime.now(timezone.utc).isoformat()
        with (run / "stdout.log").open("w", encoding="utf-8") as stdout, \
             (run / "stderr.log").open("w", encoding="utf-8") as stderr:
            process = subprocess.run(command, cwd=args.working_directory.resolve(),
                                     env=environment, stdout=stdout, stderr=stderr, check=False)
        manifest = {"case": name, "command": command, "git_commit": commit,
                    "parameter_file": str(param), "started_at_utc": started,
                    "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                    "exit_code": process.returncode, "settings": values}
        (run / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        summary.append({"case": name, "run_directory": str(run),
                        "exit_code": process.returncode})
        print(f"{name}: exit {process.returncode}")
    (root / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    failures = [item for item in summary if item["exit_code"] != 0]
    if failures:
        (root / "failures.json").write_text(
            json.dumps(failures, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(root)


if __name__ == "__main__":
    main()
