#!/usr/bin/env python3
"""Build an isolated channel-section diagnostic executable; never edit src/.

Build only:
    python3 tools/run_channel_section_diagnostic.py
Explicitly build and run with the current input (may take a full run):
    python3 tools/run_channel_section_diagnostic.py --run

All generated sources, builds and results stay in a fresh directory outside
the repository (default: /tmp). Only a failed interval-intersection check
produces diagnostic.json; ordinary output is captured in run.log. The original
exception is rethrown, and no geometric or physical acceptance rule is relaxed.
Instrumentation may affect compiler optimization: this diagnostic executable
is not a byte-identical or numerically certified production build.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent.parent
CHANNEL = Path("src/fluid/channelSections.cpp")
MARKER = "CHANNEL_SECTION_DIAGNOSTIC_JSON "


def instrument(source: str) -> str:
    """Insert failure-only reporting into a disposable source copy."""
    definition = "\nvoid positiveLinearIntegrals("
    call = ("                const int il = envelopeEdge(l, mid, maxLeft), "
            "ir = envelopeEdge(r, mid, !maxLeft);")
    if source.count(definition) != 1 or source.count(call) != 1:
        raise RuntimeError("channelSections.cpp changed: cannot safely install the diagnostic hook")
    replacement = """                int il, ir;
                const char* failedSide = "left";
                try {
                    il = envelopeEdge(l, mid, maxLeft);
                    failedSide = "right";
                    ir = envelopeEdge(r, mid, !maxLeft);
                } catch (const std::exception& error) {
                    try {
                        reportChannelSectionFailure(l, r, failedSide, k, sections.x[k],
                                                    z0, z1, mid, maxLeft, error.what());
                    } catch (...) {
                        std::cerr << "[diagnostic] Failed to serialize section details\\n";
                    }
                    throw;
                }"""
    source = source.replace(definition, '\n#include "../../channel_section_failure.inc"\n' + definition)
    return source.replace(call, replacement)


def prepare(work: Path, source: str) -> Path:
    staged = work / "source"
    staged.mkdir()
    for name in ("src", "include", "tests"):
        shutil.copytree(ROOT / name, staged / name)
    shutil.copy2(ROOT / "CMakeLists.txt", staged / "CMakeLists.txt")
    shutil.copy2(ROOT / "tools/diagnostics/channel_section_failure.inc",
                 staged / "channel_section_failure.inc")
    (staged / CHANNEL).write_text(source, encoding="utf-8")
    return staged


def run_diagnostic(work: Path, parameter: Path, source_hash: str) -> int:
    env = os.environ.copy()
    env["SIMULATION_RUN_DIR"] = str(work / "run")
    log_path = work / "run.log"
    report = None
    time_context = None
    with log_path.open("x", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(work / "build/simulation"), str(parameter)], cwd=work / "build",
            env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1,
        )
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                if line.startswith(MARKER):
                    try:
                        report = json.loads(line[len(MARKER):])
                    except json.JSONDecodeError:
                        print("[diagnostic] Invalid report JSON; raw details are in run.log", file=sys.stderr)
                else:
                    print(line, end="", flush=True)
                match = re.search(r"\[ChannelSections\] time_s=(\S+) step=(\d+)", line)
                if match:
                    time_context = {"time_s": float(match[1]), "step": int(match[2])}
            returncode = process.wait()
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            process.stdout.close()
    if report is not None:
        used_parameter = work / "run/txt/params_used.txt"
        parameter_snapshot = used_parameter if used_parameter.is_file() else parameter
        report.update(time_context or {})
        report.update({
            "simulation_exit_code": returncode,
            "unmodified_channel_source_sha256": source_hash,
            "parameter_file": str(parameter),
            "parameter_sha256": hashlib.sha256(parameter_snapshot.read_bytes()).hexdigest(),
            "openmp_environment": {key: env.get(key) for key in (
                "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES")},
            "caution": "Isolated instrumented build; compiler optimization may differ from production.",
        })
        report_path = work / "diagnostic.json"
        with report_path.open("x", encoding="utf-8") as output:
            json.dump(report, output, indent=2, ensure_ascii=False, allow_nan=False)
            output.write("\n")
        print(f"[diagnostic] Failure details: {report_path}")
    else:
        print("[diagnostic] No interval-intersection report captured; inspect run.log.")
    print(f"[diagnostic] Complete log: {log_path}")
    return returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--param", type=Path, default=ROOT / "input/param.txt")
    parser.add_argument("--work-dir", type=Path, help="new, nonexistent directory outside this repository")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--run", action="store_true", help="explicitly run the isolated simulation after building")
    args = parser.parse_args()
    parameter = args.param.resolve()
    if not parameter.is_file():
        parser.error(f"parameter file does not exist: {parameter}")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    original = (ROOT / CHANNEL).read_text(encoding="utf-8")
    source = instrument(original)  # Validate anchors before creating anything.
    if args.work_dir is None:
        work = Path(tempfile.mkdtemp(prefix="channel-section-diagnostic-", dir="/tmp"))
    else:
        work = args.work_dir.resolve()
        if work == ROOT or ROOT in work.parents:
            parser.error("--work-dir must be outside the repository")
        if work.exists():
            parser.error("--work-dir already exists; existing files will not be overwritten")
        work.mkdir(parents=True, exist_ok=False)
    print(f"[diagnostic] Isolated workspace (retained): {work}", flush=True)
    staged = prepare(work, source)
    subprocess.run(["cmake", "-S", str(staged), "-B", str(work / "build"),
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run(["cmake", "--build", str(work / "build"), "--target", "simulation",
                    "-j", str(args.jobs)], check=True)
    if not args.run:
        print("[diagnostic] Build succeeded. Simulation NOT run; use --run explicitly to build and run.")
        return 0
    # The build cwd also isolates legacy ../output/{csv,dat} geometry exports.
    return run_diagnostic(work, parameter, hashlib.sha256(original.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"[diagnostic] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("\n[diagnostic] Interrupted; isolated workspace retained.", file=sys.stderr)
        raise SystemExit(130)
