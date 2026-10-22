#!/usr/bin/env python3
"""Isolated Kanaya diagnostic validation, single-run launcher and analysis.

No parameter/model tuning, package installation, or writes into the original
worktree. All generated artifacts are under this worktree's diagnostic output.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/kanaya_diagnostic"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def preservation_audit(baseline):
    original = Path("/home/kajiyama/code/simulation-main")
    snapshot = {}
    for folder in ["input", "output", "result", "build", "tools/diagnostics"]:
        rows = []
        for base, directories, files in os.walk(original / folder):
            directories.sort()
            for name in sorted(files):
                path = Path(base) / name
                stat = path.stat()
                rows.append([str(path.relative_to(original)), stat.st_size, stat.st_mtime_ns])
        snapshot[folder] = {"count": len(rows), "sha256_metadata": hashlib.sha256(json.dumps(rows).encode()).hexdigest()}
    for key, command in [("status", ["status", "--short"]), ("diff", ["diff", "--binary"]),
                         ("cached", ["diff", "--cached", "--binary"]), ("head", ["rev-parse", "HEAD"])]:
        value = subprocess.check_output(["git"]+command, cwd=original, text=True)
        snapshot[key] = value.strip() if key == "head" else value
    if snapshot != baseline:
        raise RuntimeError("original worktree preservation audit failed")
    return {"passed": True, "method": "original Git HEAD/status/diff/index and per-file path/size/mtime metadata aggregate", "snapshot": snapshot}


def environment(enabled):
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="4", OMP_PROC_BIND="close", OMP_PLACES="cores",
               SIMULATION_KANAYA_DIAGNOSTICS="1" if enabled else "0")
    return env


def regression(path):
    arrays = {}
    with path.open("rb") as stream:
        while label := stream.read(32):
            count = struct.unpack("<Q", stream.read(8))[0]
            values = np.frombuffer(stream.read(count * 8), dtype="<f8").copy()
            if len(values) != count or np.isinf(values).any():
                raise ValueError("invalid regression data")
            arrays[label.rstrip(b"\0").decode()] = values
    return arrays


def difference(a, b):
    if a.shape != b.shape or np.isinf(a).any() or np.isinf(b).any() or not np.array_equal(np.isnan(a), np.isnan(b)):
        raise ValueError("shape/nonfinite regression mismatch")
    missing = int(np.isnan(a).sum())
    mask = np.isfinite(a)
    original_a, original_b = a, b
    a, b = a[mask], b[mask]
    delta = np.abs(a - b)
    return {"max_abs": float(delta.max(initial=0)),
            "max_rel": float((delta / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)).max(initial=0)),
            "matching_undefined_values": missing,
            "exact": bool(np.array_equal(original_a, original_b, equal_nan=True))}


def smoke(compare_only=False, open_gap=False):
    validation = OUT / "validation"
    validation.mkdir(parents=True, exist_ok=True)
    cases = [("off-build", ROOT / "build/off/kanaya_smoke_driver", False),
             ("on-disabled", ROOT / "build/on/test_kanaya_diagnostics", False),
             ("on-enabled", ROOT / "build/on/test_kanaya_diagnostics", True)]
    prefix = "smoke-open-" if open_gap else "smoke-"
    for name, binary, enabled in ([] if compare_only else cases):
        run = validation / (prefix + name)
        if run.exists():
            raise RuntimeError(f"refusing to overwrite {run}")
        run.mkdir()
        env = environment(enabled)
        env["SIMULATION_RUN_DIR"] = str(run)
        print("Starting 80-step smoke", name, flush=True)
        with (run / "run.log").open("x") as log:
            subprocess.run([str(binary), "--smoke", str(ROOT / "input" / ("kanaya_smoke_open.txt" if open_gap else "kanaya_smoke.txt"))],
                           cwd=ROOT / "build", env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    base = validation / (prefix + "off-build")
    comparisons = {}
    for name, _, _ in cases[1:]:
        run = validation / (prefix + name)
        a, b = regression(base / "regression_state.bin"), regression(run / "regression_state.bin")
        if a.keys() != b.keys():
            raise ValueError("regression field mismatch")
        results = {key: difference(a[key], b[key]) for key in a}
        for path in sorted((base / "dat").glob("*.dat")):
            results["dat/" + path.name] = difference(np.loadtxt(path), np.loadtxt(run / "dat" / path.name))
        reference_files = sorted((base / "result").glob("*.vtu"))
        actual_files = sorted((run / "result").glob("*.vtu"))
        if [p.name for p in reference_files] != [p.name for p in actual_files]:
            raise ValueError("VTU schedule mismatch")
        for path in reference_files:
            ta, tb = ET.parse(path), ET.parse(run / "result" / path.name)
            pa, pb = ta.find(".//Piece"), tb.find(".//Piece")
            if pa.attrib != pb.attrib:
                raise ValueError("VTU mesh shape mismatch")
            for key in ["connectivity", "offsets", "types", "Points"]:
                va = np.fromstring(ta.find(f'.//DataArray[@Name="{key}"]').text, sep=" ")
                vb = np.fromstring(tb.find(f'.//DataArray[@Name="{key}"]').text, sep=" ")
                results[path.name + "/" + key] = difference(va, vb)
        patterns = [r"total contact iterations\s*:\s*(\d+)", r"max icont used\s*:\s*(\d+)",
                    r"steps no contact break\s*:\s*(\d+)", r"steps converged break\s*:\s*(\d+)",
                    r"steps reached ncont\s*:\s*(\d+)"]
        logs = [(base / "run.log").read_text(), (run / "run.log").read_text()]
        counts = [[int(re.search(p, log)[1]) for p in patterns] for log in logs]
        if counts[0] != counts[1]:
            raise ValueError("contact iteration mismatch")
        comparisons[name] = {"fields": results, "contact_summary": counts[0],
                             "all_exact": all(v["exact"] for v in results.values())}
    consistency = pd.read_csv(validation / (prefix + "on-enabled") / "kanaya_diagnostic/consistency.csv")
    included = consistency[consistency.included == 1]
    if len(included) != 2 * 79 or not included.valid.eq(1).all() or included.scaled_max.max() > 1:
        raise ValueError("smoke mandatory consistency gate failed")
    report = {"steps": 80, "integrated_intervals": 79, "excluded_intervals": 1,
              "conditions": {"threads": 4, "placement": "close/cores", "dt_s": 1e-5,
                             "nmode": 100, "ps_pa": 1000, "zetaL": .06, "zetaR": .06,
                             "initialGapMm": .1 if open_gap else 0},
              "comparisons": comparisons,
              "max_consistency_scaled_residual": float(included.scaled_max.max()),
              "passed": all(v["all_exact"] for v in comparisons.values())}
    report["max_abs_flow_m3_s"] = float(np.abs(regression(base / "regression_state.bin")["Ug_history"]).max())
    if open_gap and report["max_abs_flow_m3_s"] == 0:
        raise RuntimeError("additional open-flow smoke did not exercise nonzero flow")
    dump(validation / ("smoke_open_comparison.json" if open_gap else "smoke_comparison.json"), report)
    print(json.dumps({k: v for k, v in report.items() if k != "comparisons"}, indent=2), flush=True)
    if not report["passed"]:
        raise RuntimeError("numerical changes detected; no tolerance relaxation; full run forbidden")


def production():
    validation = OUT / "validation"
    if (OUT / "production_attempt.json").exists():
        raise RuntimeError("production attempt already recorded; repeat full runs forbidden")
    gate = json.loads((validation / "smoke_comparison.json").read_text())
    open_gate = json.loads((validation / "smoke_open_comparison.json").read_text())
    if not gate["passed"] or not open_gate["passed"]:
        raise RuntimeError("smoke gate not passed")
    audit = {}
    for name in ["src/simulation/simulation.cpp", "src/fluid/flowModel.cpp", "src/forces/forceCalculatorFluid.cpp",
                 "src/forces/forceCalculatorContact.cpp", "include/forces/ForceCalculator.h"]:
        original = subprocess.check_output(["git", "show", "HEAD:" + name], cwd=ROOT, text=True)
        stripped = re.sub(r"#ifdef ENABLE_KANAYA_DIAGNOSTICS\n.*?#endif\n", "", (ROOT / name).read_text(), flags=re.S)
        if stripped != original:
            raise RuntimeError("non-observational source change: " + name)
        audit[name] = "OFF source identical to current HEAD after removing diagnostic #ifdef blocks"
    for name in ["src/structure/state.cpp", "src/structure/rungeStep.cpp", "src/fluid/channelSections.cpp",
                 "src/fluid/separationSelection.cpp", "src/forces/modalProjector.cpp"]:
        if (ROOT / name).read_bytes() != subprocess.check_output(["git", "show", "HEAD:" + name], cwd=ROOT):
            raise RuntimeError("protected source changed: " + name)
        audit[name] = "byte identical to HEAD"
    for obj in (ROOT / "build/on/CMakeFiles/simulation.dir/src").rglob("*.o"):
        if obj.name == "main.cpp.o":
            continue
        driver = Path(str(obj).replace("simulation.dir", "test_kanaya_diagnostics.dir"))
        if sha(obj) != sha(driver):
            raise RuntimeError("smoke/production object mismatch: " + str(obj))
    audit["smoke_vs_production_solver_objects"] = "byte-identical"
    tests = {}
    for config in ["off", "on"]:
        tests[config] = subprocess.check_output(["ctest", "--test-dir", str(ROOT / "build" / config),
                                                 "--output-on-failure"], cwd=ROOT, text=True)
    if not (ROOT / "build/on/profile_contact_check.o").exists():
        raise RuntimeError("PROFILE_CALCDIS compatibility compile missing")
    audit["profile_contact_compile"] = "PROFILE_CALCDIS and ENABLE_KANAYA_DIAGNOSTICS simultaneous compile passed"
    dump(validation / "source_and_test_audit.json", {"sources": audit, "ctest": tests})
    preflight = json.loads((OUT / "manifest.json").read_text())
    if sha(ROOT / "input/kanaya_param.txt") != preflight["original_input_sha256"]:
        raise RuntimeError("production parameters differ from protected original working-tree input")
    # Exclusive marker: even a failed production attempt is never retried here.
    attempt = OUT / "production_attempt.json"
    with attempt.open("x") as stream:
        json.dump({"started_unix": time.time(), "full_run_attempts": 1, "status": "running"}, stream)
    binary = ROOT / "build/on/simulation"
    param = ROOT / "input/kanaya_param.txt"
    inputs = [param] + [ROOT / "input/M5_test" / name for name in
                       ["M5_freq_kawahara_mesh7_soft.txt", "M5_mode_kawahara_mesh7_soft.vtu", "M5_surface_kawahara_mesh7_soft.nas"]]
    provenance = {"git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  "git_branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
                  "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
                  "binary": str(binary), "binary_sha256": sha(binary),
                  "inputs": [{"path": str(p), "resolved_path": str(p.resolve()), "sha256": sha(p)} for p in inputs],
                  "environment": {k: v for k, v in environment(True).items() if k.startswith("OMP_") or "KANAYA" in k},
                  "build": "Intel icpx 2025.1.0; C++17; -O3; OpenMP; Release; ENABLE_KANAYA_DIAGNOSTICS=ON",
                  "solver_source_diff": subprocess.check_output(["git", "diff"], cwd=ROOT, text=True)}
    dump(OUT / "validation/production_provenance.json", provenance)
    env = environment(True)
    env["SIMULATION_RUN_DIR"] = str(ROOT / "output")
    start = time.time()
    with (OUT / "production.log").open("x") as log:
        code = subprocess.run([str(binary), str(param)], cwd=ROOT / "build", env=env,
                              stdout=log, stderr=subprocess.STDOUT).returncode
    dump(attempt, {"started_unix": start, "finished_unix": time.time(), "elapsed_s": time.time() - start,
                   "full_run_attempts": 1, "exit_code": code, "status": "complete" if code == 0 else "failed"})
    if code:
        (OUT / "report.md").write_text(
            "# Kanaya診断：本番ランで停止\n\n"
            "結論1：本番試行は1回実行したが、計算または必須整合性ゲートで停止した。\n"
            "結論2：物理解釈は行わない。H1〜H5は未判定。\n"
            "結論3：最終区間除外の例外を他の不一致へ拡張せず、再実行も行わない。\n\n"
            f"終了コード：{code}\n\n```text\n" + (OUT / "production.log").read_text()[-6000:] + "\n```\n")
        raise RuntimeError(f"production failed ({code}); no retry allowed; inspect production.log")
    print("Single production run completed", flush=True)


def analyze():
    os.environ.setdefault("MPLCONFIGDIR", str(OUT / "validation/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    previous_manifest = json.loads((OUT / "manifest.json").read_text())
    preservation = preservation_audit(previous_manifest["preservation_snapshot"])

    meta = json.loads((OUT / "run_metadata.json").read_text())
    attempt = json.loads((OUT / "production_attempt.json").read_text())
    if not meta["complete"] or attempt["exit_code"] != 0 or attempt["full_run_attempts"] != 1:
        raise RuntimeError("a successful single production run is required")
    n, dt = meta["integrated_intervals"], meta["dt_s"]
    if n != 19999 or meta["excluded_intervals"] != 1 or meta["excluded_step"] != 19999:
        raise ValueError("production interval exclusion/accounting mismatch")
    consistency = pd.read_csv(OUT / "consistency.csv")
    included = consistency[consistency.included == 1]
    excluded = consistency[consistency.included == 0]
    if len(included) != 2*n or not included.valid.eq(1).all() or included.scaled_max.max() > 1:
        raise RuntimeError("mandatory consistency gate failed; no physical interpretation")
    if len(excluded) != 2 or not excluded.step.eq(19999).all():
        raise RuntimeError("unexpected excluded consistency samples")
    flow = pd.read_csv(OUT / "flow_geometry.csv")
    if len(flow) != 20000 or not np.array_equal(flow.step.to_numpy(), np.arange(20000)):
        raise ValueError("missing/reordered flow observation")
    forces = pd.read_csv(OUT / "force_regions.csv")
    contacts = pd.read_csv(OUT / "contact_regions.csv")
    if contacts.step.max() >= n or len(contacts) != 4*n:
        raise ValueError("last contact interval was not fully excluded")
    summary = pd.read_csv(OUT / "modal_energy_summary.csv")
    nm = meta["nmode"]
    if nm != [100, 100]:
        raise ValueError("unexpected modal dimension")
    raw = np.memmap(OUT / "modal_native.bin", dtype="<f8", mode="r", shape=(n, 2, 100, 7))
    if (OUT / "modal_native.bin").stat().st_size != n*2*100*7*8 or not np.isfinite(raw).all():
        raise ValueError("invalid native modal data")
    q_split_scaled = np.abs(raw[..., 0]+raw[..., 1]-raw[..., 2]) / (1e-12+1e-8*(np.abs(raw[..., 0])+np.abs(raw[..., 1])))
    q_split_scaled_max = float(q_split_scaled.max())
    if q_split_scaled_max > 1:
        raise RuntimeError("individual-mode Q decomposition failed; no physical interpretation")
    del q_split_scaled
    times = np.arange(n)*dt
    evaluation = times >= .1
    records, phase = [], []
    selected = {}
    for side, name in enumerate(["L", "R"]):
        all_modes = summary[(summary.side == name) & (summary.window == "all")]
        top = all_modes.assign(abs_work=all_modes.W_fluid.abs()).nlargest(5, "abs_work")["mode"].to_list()
        selected[name] = sorted(set([1, 4, 5] + top))
        for mode in range(100):
            data = np.asarray(raw[:, side, mode, :])
            target = all_modes[all_modes["mode"] == mode+1].iloc[0]
            for column, expected in [(0, target.W_pos), (1, target.W_neg), (2, target.W_fluid), (3, target.W_contact_reference)]:
                measured = float(np.sum(data[:, column]*data[:, 6]))
                if abs(measured-expected) > 1e-12+1e-8*abs(expected):
                    raise ValueError("native modal history/summary mismatch")
            if mode+1 not in selected[name]:
                continue
            q, v, qforce = data[evaluation, 4], data[evaluation, 5], data[evaluation, 2]
            fq, fv, ff = [np.fft.rfft((x-x.mean())*np.hanning(len(x))) for x in [q, v, qforce]]
            peak = 1+int(np.argmax(np.abs(fq[1:])))
            hz = float(np.fft.rfftfreq(len(q), dt)[peak])
            phase_deg = float(np.angle(ff[peak]*np.conjugate(fv[peak]), deg=True)) if abs(fv[peak]*ff[peak]) > 0 else None
            correlation = float(np.corrcoef(v, qforce)[0, 1]) if np.std(v)>0 and np.std(qforce)>0 else None
            phase.append({"side": name, "mode": mode+1, "dominant_q_peak_Hz": hz,
                          "Q_vs_qdot_phase_deg": phase_deg, "Q_qdot_correlation": correlation})
            records.append(pd.DataFrame({"step": np.arange(n), "time_s": times, "side": name, "mode": mode+1,
                                         "Q_pos": data[:, 0], "Q_neg": data[:, 1], "Q_fluid": data[:, 2],
                                         "Q_contact_adopted_predicted_stage": data[:, 3], "q_n": data[:, 4],
                                         "qdot_n": data[:, 5], "delta_q": data[:, 6],
                                         "power_left_endpoint": data[:, 2]*data[:, 5],
                                         "xSep_mm": flow.xSep_mm.to_numpy()[:n]}))
    pd.concat(records, ignore_index=True).to_csv(OUT / "modal_selected_timeseries.csv", index=False)
    ev = flow[(flow.interval_included == 1) & (flow.time_s >= .1)]
    ef = forces[forces.time_s >= .1]
    ec = contacts[contacts.time_s >= .1]
    valid_area = ev.Amin_mm2[np.isfinite(ev.Amin_mm2)]
    if len(valid_area) != len(ev):
        raise ValueError("undefined minimum area: do not substitute closure")
    force_stats = [r for r in meta["force_statistics"] if r["window"] == 3]
    pos = next(r for r in force_stats if r["side"] == "L" and r["pressure_sign"] == 1)
    neg = next(r for r in force_stats if r["side"] == "L" and r["pressure_sign"] == -1)
    negative = ef[(ef.pressure_sign == -1) & (ef.side == "L")]
    end_contacts = ec[(ec.side == "L") & (ec.ap_region == "end")]
    interior_contacts = ec[(ec.side == "L") & (ec.ap_region == "interior")]
    pair_total = int(end_contacts.pairs.sum()+interior_contacts.pairs.sum())
    end_fraction = float(end_contacts.pairs.sum()/pair_total) if pair_total else None
    evaluation_summary = summary[summary.window == "evaluation"]
    left_summary = evaluation_summary[evaluation_summary.side == "L"]
    modal_totals = {key: float(left_summary[key].sum()) for key in
                    ["W_pos", "W_neg", "W_fluid", "D_struct", "W_contact_reference"]}
    work_decomposition = {}
    for column, label in [(0, "positive"), (1, "negative"), (2, "fluid")]:
        data = np.asarray(raw[evaluation, 0, :, :])
        total = float(np.sum(data[..., column]*data[..., 6]))
        bias = float(np.sum(np.mean(data[..., column], axis=0)*np.sum(data[..., 6], axis=0)))
        work_decomposition[label] = {"total": total, "mean_force_endpoint_work": bias,
                                     "force_fluctuation_work": total-bias}
    wrong_normal = int(forces.wrong_normal_count.sum())
    peaks = np.abs(np.fft.rfft((ev.Ug_m3_s.to_numpy()-ev.Ug_m3_s.mean())*np.hanning(len(ev))))
    frequencies = np.fft.rfftfreq(len(ev), dt)
    maxima = [k for k in range(1, len(peaks)-1) if peaks[k]>peaks[k-1] and peaks[k]>=peaks[k+1]]
    top_peaks = sorted(maxima, key=lambda k: peaks[k], reverse=True)[:5]
    max_columns = ["force_abs", "Q_abs", "work_split_abs", "nodal_work_abs", "contact_Q_abs", "scaled_max"]
    gram = meta["gram"]
    gram_off = max(abs(row[3]) for row in gram if row[1] != row[2])
    gram_diag = max(abs(row[3]-1) for row in gram if row[1] == row[2])
    centroid = None
    if len(negative) and negative.opening_Fy_N.abs().sum()>0:
        centroid = float(np.average(.5*(negative.x_min_mm+negative.x_max_mm), weights=negative.opening_Fy_N.abs()))
    metrics = {"evaluation_window_s": [.1, meta["work_interval_end_s"]], "evaluation_intervals": len(ev),
               "mean_Ug_ml_s": float(ev.Ug_m3_s.mean()*1e6), "mean_Amin_mm2": float(valid_area.mean()),
               "min_Amin_mm2": float(valid_area.min()), "area_closure_fraction": float((valid_area <= .01).mean()),
               "mean_Pg_Pa": float(ev.Pg_Pa.mean()), "mean_Psep_Pa": float(ev.Psep_Pa.mean()),
               "mean_Pd0_Pa": float(ev.Pd0_Pa.mean()), "negative_force_regions_sampled": len(negative),
               "negative_x_min_mm": float(negative.x_min_mm.min()) if len(negative) else None,
               "negative_x_max_mm": float(negative.x_max_mm.max()) if len(negative) else None,
               "negative_pressure_min_Pa": float(negative.p_min_Pa.min()) if len(negative) else None,
               "negative_force_centroid_x_mm_binned": centroid, "force_statistics_evaluation": force_stats,
               "wrong_normal_samples": wrong_normal, "contact_AP_end_pair_fraction": end_fraction,
               "contact_step_fraction": float(ev.contact_flag.mean()),
               "contact_end_active_fraction": float((end_contacts.pairs>0).mean()),
               "contact_interior_active_fraction": float((interior_contacts.pairs>0).mean()),
               "mean_contact_iterations": float(ev.contact_iterations.mean()),
               "contact_end_mean_opening_N": float(end_contacts.opening_Fy_N.mean()),
               "contact_interior_mean_opening_N": float(interior_contacts.opening_Fy_N.mean()),
               "maximum_contact_penetration_m": float(ec.max_penetration_m.max()),
               "mean_delta_gap_mm": float(ev.mean_delta_gap_mm.mean()),
               "monitor_uyL_mean_mm": float(ev.uyL_mm.mean()), "monitor_uyL_ac_rms_mm": float(ev.uyL_mm.std(ddof=0)),
               "monitor_gap_min_mm": float(ev.monitor_gap_mm.min()),
               "xSep_range_mm": [float(ev.xSep_mm.min()), float(ev.xSep_mm.max())],
               "mean_xSep_mm": float(ev.xSep_mm.mean()), "mean_xBlendEnd_mm": float(ev.xBlendEnd_mm.mean()),
               "mean_max_delta_gap_x_mm": float(ev.max_delta_gap_x_mm.mean()),
               "mean_monitor_x_mm": float(ev.monitor_x_mm.mean()),
               "mean_sep_node_x_mm": float(ev.sep_node_x_mm.mean()),
               "sep_node_uyL_ac_rms_mm": float(ev.sep_node_uyL_mm.std(ddof=0)),
               "Lg_range": [float(ev.Lg.min()), float(ev.Lg.max())],
               "Rv_range": [float(ev.Rv.min()), float(ev.Rv.max())],
               "Rk_range": [float(ev.Rk.min()), float(ev.Rk.max())],
               "selected_gram_max_off_diagonal": gram_off, "selected_gram_max_diagonal_error": gram_diag,
               "left_normalized_modal_totals": modal_totals, "selected_modes": selected,
               "phase": phase, "Ug_spectral_peaks_Hz_not_f0": [float(frequencies[k]) for k in top_peaks],
               "consistency_maxima": {key: float(included[key].max()) for key in max_columns}}
    metrics["individual_mode_Q_split_max_scaled_residual"] = q_split_scaled_max
    metrics["left_evaluation_mean_force_work_decomposition"] = work_decomposition
    metrics["positive_negative_opening_force_magnitude_ratio"] = pos["mean_opening_N"]/abs(neg["mean_opening_N"])
    metrics["monitor_minus_xSep_mean_mm"] = float((ev.monitor_x_mm-ev.xSep_mm).mean())
    metrics["maximum_opening_minus_xSep_mean_mm"] = float((ev.max_delta_gap_x_mm-ev.xSep_mm).mean())
    metrics["monitor_beyond_blend_fraction"] = float((ev.monitor_x_mm>=ev.xBlendEnd_mm).mean())
    metrics["negative_force_by_zone_N_decimated"] = (negative.groupby("zone").opening_Fy_N.sum()/ef.step.nunique()).to_dict()
    metrics["negative_force_by_AP_N_decimated"] = (negative.groupby("ap_region").opening_Fy_N.sum()/ef.step.nunique()).to_dict()
    metrics["contact_pair_positions"] = {}
    for ap, group in ec[ec.side == "L"].groupby("ap_region"):
        if group.pairs.sum():
            metrics["contact_pair_positions"][ap] = {"pair_weighted_mean_x_mm": float(np.average(group.x_mean_mm, weights=group.pairs)),
                                                      "pair_weighted_mean_z_mm": float(np.average(group.z_mean_mm, weights=group.pairs))}
    # Reuse the user's archived run read-only. This is an additional comparison,
    # not a second full computation, and does not replace the OFF/ON gate.
    original = Path("/home/kajiyama/code/simulation-main")
    reference = {"provenance_limit": "archived parameters match; complete historical binary/build/thread provenance unavailable"}
    if sha(original / "output/txt/params_used.txt") == sha(ROOT / "input/kanaya_param.txt"):
        reference["parameter_snapshot_identical"] = True
        reference["fields"] = {}
        for path in sorted((ROOT / "output/dat").glob("*.dat")):
            old_path = original / "output/dat" / path.name
            if old_path.exists():
                try:
                    reference["fields"][path.name] = difference(np.loadtxt(old_path), np.loadtxt(path))
                except ValueError as error:
                    reference["fields"][path.name] = {"exact": False, "comparable": False, "reason": str(error)}
        old_flow = np.loadtxt(original / "output/dat/airflow_vt.dat")
        old_area = np.loadtxt(original / "output/dat/area.dat")
        mask = (old_flow[:, 0] >= 10000) & (old_flow[:, 0] < 20000)
        reference["evaluation_decimated_mean_flow_ml_s"] = float(old_flow[mask, 1].mean()*1e6)
        reference["evaluation_decimated_mean_Amin_mm2"] = float(np.nanmin(old_area[mask, 2:-1], axis=1).mean())
        reference["all_compared_scalar_outputs_exact"] = all(v["exact"] for v in reference["fields"].values())
        reference["sampled_vtu_byte_comparison"] = {}
        for k in [0, 50, 100, 249]:
            name = f"deform_combined{k:04d}.vtu"
            ref, current = original / "result" / name, ROOT / "output/result" / name
            reference["sampled_vtu_byte_comparison"][name] = sha(ref) == sha(current)
    else:
        reference["parameter_snapshot_identical"] = False
    metrics["archived_reference_comparison"] = reference
    dump(OUT / "analysis_metrics.json", metrics)

    figures = OUT / "figures"
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, side in zip(axes, ["L", "R"]):
        f = ef[ef.side == side]
        denominator = f.step.nunique()
        for pressure_sign, label in [(1, "positive pressure"), (-1, "negative pressure")]:
            v = f[f.pressure_sign == pressure_sign].groupby("x_bin_mm").opening_Fy_N.sum()/denominator
            ax.plot(v.index, v.values, marker=".", label=label)
        ax.axhline(0, color="k", lw=.5)
        ax.axvline(ev.xSep_mm.mean(), ls="--", color="gray", label="mean xSep")
        ax.set(title=side, xlabel="deformed physical x bin [mm]", ylabel="mean opening force [N]")
        ax.legend()
    fig.tight_layout(); fig.savefig(figures / "01_forces_by_x.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    for ax, side in zip(axes, ["L", "R"]):
        modes = selected[side]
        table = evaluation_summary[evaluation_summary.side == side].set_index("mode").loc[modes]
        x = np.arange(len(modes))
        for offset, key, label in [(-.3, "W_pos", "positive pressure W"), (-.1, "W_neg", "negative pressure W"),
                                   (.1, "D_struct", "structural D (normalized)"), (.3, "W_contact_reference", "contact reference W")]:
            ax.bar(x+offset, table[key], width=.19, label=label)
        ax.set_xticks(x, modes); ax.axhline(0, color="k", lw=.5)
        ax.set(title=side, xlabel="mode (1-based)", ylabel="work / normalized modal units")
        ax.legend(ncol=2)
    fig.tight_layout(); fig.savefig(figures / "02_modal_work.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 6))
    plot_modes = list(dict.fromkeys([1, 4, 5] + selected["L"]))[:4]
    for ax, mode in zip(axes.flat, plot_modes):
        qforce = np.array(raw[:, 0, mode-1, 2]); v = np.array(raw[:, 0, mode-1, 5])
        for data, label in [(qforce, "Q_fluid,n"), (v, "qdot_n")]:
            ac = data[evaluation]-data[evaluation].mean()
            ax.plot(times[evaluation], ac/max(np.max(np.abs(ac)), 1e-300), label=label, lw=.7)
        ax.set(title=f"L mode {mode}: AC normalized", xlabel="t_n [s]"); ax.legend()
    fig.tight_layout(); fig.savefig(figures / "03_modal_phase.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes[0, 0].plot(flow.time_s, flow.Amin_mm2); axes[0, 0].axhline(.01, color="r", ls="--")
    axes[0, 0].set(ylabel="Amin [mm2]", xlabel="time [s]")
    axes[0, 1].plot(flow.time_s, flow.uyL_mm, label="monitor uyL")
    axes[0, 1].plot(flow.time_s, flow.sep_node_uyL_mm, label="near-separation uyL")
    axes[0, 1].plot(flow.time_s, flow.mean_delta_gap_mm, label="mean gap increase")
    axes[0, 1].set(ylabel="[mm]", xlabel="time [s]"); axes[0, 1].legend()
    for ap in ["end", "interior"]:
        c = contacts[(contacts.side == "L") & (contacts.ap_region == ap)]
        axes[1, 0].plot(c.time_s, c.pairs, label=ap, lw=.7)
    axes[1, 0].set(ylabel="adopted contact pairs", xlabel="interval t_n [s]"); axes[1, 0].legend()
    axes[1, 1].plot(flow.time_s, flow.Pg_Pa, label="Pg")
    axes[1, 1].plot(flow.time_s, flow.Psep_Pa, label="Psep")
    axes[1, 1].plot(flow.time_s, flow.Pd0_Pa, label="Pd0")
    axes[1, 1].set(ylabel="pressure [Pa]", xlabel="time [s]"); axes[1, 1].legend()
    fig.tight_layout(); fig.savefig(figures / "04_contact_and_gap.png", dpi=150); plt.close(fig)

    def f(value):
        return "未測定" if value is None else f"{value:.6g}"
    total_fluid = modal_totals["W_fluid"]
    rows = ["|side|mode|f [Hz]|W+|W−|Wfluid|Dstruct（正規化）|Wcontact参考|",
            "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for side in ["L", "R"]:
        for mode in selected[side]:
            r = evaluation_summary[(evaluation_summary.side == side) & (evaluation_summary["mode"] == mode)].iloc[0]
            rows.append(f"|{side}|{mode}|{r.f_Hz:.6g}|{r.W_pos:.6g}|{r.W_neg:.6g}|{r.W_fluid:.6g}|{r.D_struct:.6g}|{r.W_contact_reference:.6g}|")
    force_rows = ["|side|圧力符号|平均開大方向合力 [N]|平均作用面積 [mm²]|最小/最大合力 [N]|",
                  "|---|---|---:|---:|---|"]
    for r in force_stats:
        force_rows.append(f"|{r['side']}|{r['pressure_sign']:+d}|{r['mean_opening_N']:.6g}|{r['mean_area_m2']*1e6:.6g}|{r['min_opening_N']:.6g} / {r['max_opening_N']:.6g}|")
    phase_rows = ["|side|mode|qの主要ピーク [Hz]|Qfluid対qdotの位相差 [deg]|AC相関|",
                  "|---|---:|---:|---:|---:|"]
    for r in phase:
        phase_rows.append(f"|{r['side']}|{r['mode']}|{f(r['dominant_q_peak_Hz'])}|{f(r['Q_vs_qdot_phase_deg'])}|{f(r['Q_qdot_correlation'])}|")
    h1 = "支持（観測上）" if pos["mean_opening_N"] > abs(neg["mean_opening_N"]) else "単純な正圧優勢は支持されない"
    h2 = "検査範囲では支持されない" if wrong_normal == 0 else "符号不整合の疑いあり"
    end_text = "未測定" if end_fraction is None else f"{100*end_fraction:.6g}%"
    next_case = "AP端の最初/最後3列の接触力だけを除く1条件比較（apContactDiagnosticEnabled=1、apContactExcludedRowsPerEnd=3）。これは1つの介入条件であり、他の物理パラメータは固定する。" if end_fraction is not None and end_fraction > .5 else \
        "flowBlendLengthMmだけを0.5→0.25 mmにする1条件比較。他条件を固定し、負圧の壁面作用域とモード仕事の変化でH1/H5を識別する。現行parserはこの名前付きキーに未対応なので、実施時は既存の末尾位置オプションを0/40/0.25/1（forceDirection/contactReferenceFrequencyHz/flowBlendLengthMm/flowSeparationAreaRatio）として他の有効値を維持する。今回は未実行。"
    report = f"""# Kanaya流体力・モード仕事診断

結論1：評価窓の左正圧開大合力 {pos['mean_opening_N']:.6g} N、負圧開大合力 {neg['mean_opening_N']:.6g} N。H1は{h1}。
結論2：平均Ug {metrics['mean_Ug_ml_s']:.6g} mL/s、平均Amin {metrics['mean_Amin_mm2']:.6g} mm²、面積閉鎖率 {metrics['area_closure_fraction']:.6g}。周期的閉鎖の再現とは判定しない。
結論3：19,999区間の必須整合性を検証。最終1区間だけ承認どおり除外し、本体の20,000ステップは通常実行した。

## 実行・検証

- worktree: {ROOT}、HEAD: df024238b39d610da42df8faa9556dcc18a98880。
- 作業ブランチ: codex/kanaya-force-energy-implementation-20261021。元ブランチcodex/kanaya-force-energy-plan-20261021は維持。commit/pushは未実施。
- 既存ファイルの変更はCMakeLists.txt、include/forces/ForceCalculator.h、src/fluid/flowModel.cpp、src/forces/forceCalculatorFluid.cpp、src/forces/forceCalculatorContact.cpp、src/simulation/simulation.cppの6件（設定と計測フックのみ）。
- 追加はinclude/diagnostics/KanayaDiagnostics.h、src/diagnostics/KanayaDiagnostics.cpp、tools/kanaya_diagnostics/analyze.py、tests/test_kanaya_diagnostics.cppの4件。
- Release / Intel icpx 2025.1.0 / -O3 / OpenMP、4 threads、close/cores。
- 診断ビルドオプションと実行フラグは既定OFF。ONはENABLE_KANAYA_DIAGNOSTICS=ONとSIMULATION_KANAYA_DIAGNOSTICS=1の両方が必要。
- CTest：OFF 3/3、ON 4/4。合成定数断面・鏡映正負圧・実節点仕事・符号反転・閾値超過検出・最終区間だけの除外を検証。PROFILE_CALCDISとの同時コンパイルも成功。
- 80ステップの閉鎖初期形状、および初期間隙0.1 mmの開流路で、各々OFFビルド／ONビルド無効／ON有効を比較。高精度流量・モード状態・全節点座標・最終採用接触力・候補キャッシュ・通常DAT・VTU・接触反復統計を検証。詳細はvalidation/*comparison.json。
- 本番試行は1回のみ。nstep=20000、dt=1e-5、100 modes、Ps=1000 Pa、zetaL/R=0.06、kc1=5、kc2=6000、kc3=5000、contactReferenceFrequencyHz=40、初期間隙0、L_vt=0。
- 本番の壁時計時間は {attempt['elapsed_s']:.3f} s。既存タイマーは診断処理の一部や初期化を含まないため、これだけから診断オーバーヘッド率は算出しない。
- 既存出力の入力スナップショット一致：{reference['parameter_snapshot_identical']}。既存スカラー出力との追加比較はanalysis_metrics.json/archived_reference_comparison。過去バイナリ・ビルド・スレッドの完全証跡はないため、同条件OFF/ON比較とは区別する。
- 元の未コミット入力を専用input/kanaya_param.txtへコピー。元入力・結果・buildは不変更。VTU/NAS/周波数は元データへのリンクを読取り用途に限定。
- 元作業ツリーの最終保護監査は合格。Git HEAD/status/diff/indexと、input/output/result/build/既存診断のファイル一覧・サイズ・更新時刻の集約値が作業開始時と一致。
- state.cpp、Newmark、投影、断面・剥離式は変更なし。変更した既存ソースも診断#ifdefを除去するとHEADと同一。スモークと本番のソルバオブジェクトはバイト一致。

## 最終区間の明示除外と時間対応

除外は n=19999 の1区間 [0.19999,0.20000] s。既存State境界判定が再構成をスキップしモードだけ確定する不整合は未修正。
仕事と節点/モード検証の対象は n=0..19998、積分範囲は [0,0.19999] s。最後の力の左端時刻は0.19998 s。
最終区間は接触仕事、接触位置、接触ペア・力統計からも除外。計算本体の最終ステップは削減・変更していない。
全流体係数行は20,000件保持しinterval_included列で識別する。すべての物理解釈・評価統計は除外後のデータを使う。
流体仕事はFfluid,nとq_(n+1)-q_n、追加パワーはQ_n*qdot_n。
接触は最終採用された予測反復荷重であり、流体と同時刻の荷重とは扱わない。接触仕事は参考離散仕事で、散逸とは呼ばない。
非最終区間の不一致は一切例外扱いせず、相対1e-8 + 絶対1e-12を超えたら停止する実装。
再現した既存バグはvalidation/state_endpoint.cpp。別変更での修正・検証待ちであり、今回コミットしない。

## 内部整合性

- 非最終区間39,998件（左右）のvalid=1。最大規格化残差 {metrics['consistency_maxima']['scaled_max']:.6g}（合格限界1）。
- F加法分解 最大絶対差 {metrics['consistency_maxima']['force_abs']:.6g} N。
- Q加法分解 最大絶対差 {metrics['consistency_maxima']['Q_abs']:.6g}（正規化一般化力）。
- 全有効区間×左右×100モードを個別のQスケールで再検査した最大規格化残差 {q_split_scaled_max:.6g}（限界1）。
- W加法分解 最大絶対差 {metrics['consistency_maxima']['work_split_abs']:.6g}。
- 実際の節点変位差（mm→m）対モード仕事の最大絶対差 {metrics['consistency_maxima']['nodal_work_abs']:.6g}。
- 最終採用Qfluid+Qcontact対実ソルバfiの最大絶対差 {metrics['consistency_maxima']['contact_Q_abs']:.6g}。
- 全時間および各時間窓で正負別・合計の累積節点/モード仕事も検査済み（run_metadata.json/cumulative_work_checks）。
- 小値の絶対下限は丸め誤差・座標差分の相殺に対する指示書の基準。ゼロや欠損を物理閉鎖へ置換せず、NaN位置の一致を比較した。

## H1/H2：負圧が実際に作る力

評価窓は [0.1,0.19999] s、{len(ev)}区間。下表の平均と極値は毎区間集計、空間分布は20ステップ毎。

{chr(10).join(force_rows)}

左の負圧作用域（記録サンプル）はx={f(metrics['negative_x_min_mm'])}〜{f(metrics['negative_x_max_mm'])} mm、最低割当圧力 {f(metrics['negative_pressure_min_Pa'])} Pa。
作用面積は実際に適用された面積ベクトル|dA|の和であり、別の表面積推定を追加していない。
負圧の力重み付き作用中心は約 {f(centroid)} mm（0.5 mmビン内の範囲中点による近似）。AP端/中央、upstream/blend/downstreamの内訳はforce_regions.csv。
左負圧合力の記録サンプル平均は、上流 {metrics['negative_force_by_zone_N_decimated'].get('upstream', float('nan')):.6g} N、ブレンド帯 {metrics['negative_force_by_zone_N_decimated'].get('blend', float('nan')):.6g} N、AP端 {metrics['negative_force_by_AP_N_decimated'].get('end', float('nan')):.6g} N、中央 {metrics['negative_force_by_AP_N_decimated'].get('interior', float('nan')):.6g} N。
各節点で実際に使ったp、面積ベクトル、最終f=p*dAを観測。法線・圧力に対する開大符号違反は記録サンプル {wrong_normal} 件。
したがってH2は{h2}。ただし対象表面の物理的妥当性や流路近似まで証明したものではない。

![forces](figures/01_forces_by_x.png)

## H3/H4：モード仕事・減衰・接触

評価窓の主要モード（1始まり、全100モードを積分し1/4/5と全時間|Wfluid|上位5を表示）：

{chr(10).join(rows)}

左の正規化モード合計：W+={modal_totals['W_pos']:.6g}、W−={modal_totals['W_neg']:.6g}、Wfluid={total_fluid:.6g}、Dstruct={modal_totals['D_struct']:.6g}、接触参考W={modal_totals['W_contact_reference']:.6g}。
正圧開大力/負圧閉鎖力の平均絶対値比は {metrics['positive_negative_opening_force_magnitude_ratio']:.6g}。
圧力の符号から仕事の符号は決めていない。QΔqを実測し、瞬間パワーも左端時刻に固定した。
補助的な代数分解 W=Q平均×ΣΔq+Σ(Q−Q平均)Δq では、左の平均力による端点変化の仕事 {work_decomposition['fluid']['mean_force_endpoint_work']:.6g}、変動力の仕事 {work_decomposition['fluid']['force_fluctuation_work']:.6g}。
従って正味Wfluidが負でも「流体から振動への供給が全くない」とは結論しない。これは平均開大の変化と変動を分ける参考分解であり、新しいエネルギー保存則ではない。
選択モードGram行列の最大非対角絶対値 {gram_off:.6g}、対角の1からの最大偏差 {gram_diag:.6g}。
節点等分質量による各モード正規化はモード間直交を保証しない。従って合計を厳密な独立モードの物理エネルギー収支と断定しない。
H3はモード別の仕事・位相の観測まで判定可能だが、「位相が原因で閉鎖を阻害する」という因果は未判定。位相の数値はanalysis_metrics.json/phase、波形はmodal_selected_timeseries.csv。

{chr(10).join(phase_rows)}

位相は評価窓のq主要ACピークにおけるQfluidとqdot_nのFFT位相差、相関は同じ評価窓のゼロ遅れAC相関。
短い非定常窓での参考値であり、微小振幅モードや複数ピークから単一の自励周波数を断定しない。

接触のある区間割合 {metrics['contact_step_fraction']:.6g}、平均反復数 {metrics['mean_contact_iterations']:.6g}。
接触合力は実際にfdisへ分配する最終採用ペア反力の領域別和、接触位置はそのペアの算術平均。
Qcontactと接触参考仕事は最終採用されたprevFdisの節点荷重から直接投影しており、最後のdiagnosticOnly再評価では置換しない。
左接触ペアのAP端割合 {end_text}、端の接触あり割合 {metrics['contact_end_active_fraction']:.6g}、中央 {metrics['contact_interior_active_fraction']:.6g}。
左の接触開大合力平均は端 {metrics['contact_end_mean_opening_N']:.6g} N、中央 {metrics['contact_interior_mean_opening_N']:.6g} N。
接触ペアの位置（左、ペア数重みの平均x）はanalysis_metrics.json/contact_pair_positionsに記録。
AP端は両端を合算するため、z_meanが中央付近でも中央接触を意味しない。個々のjごとの位置・極値は今回保存していない。
端だけが常時接触しているわけではなく、中央の接触の有無も同時に記録して区別した。
最大めり込み {metrics['maximum_contact_penetration_m']:.6g} m。AP端はj<3またはj>=nsurfz-3。境界j自体を接触ループが除く点にも注意。
H4の接触局在と減衰量は観測済みだが、阻害要因であるかは単一ランでは未判定。ばねの蓄積/返還を散逸とはせず、D_contactは作成していない。

![work](figures/02_modal_work.png)
![phase](figures/03_modal_phase.png)

## H5：幾何・流量・平均開大・位置のずれ

- 平均Ug {metrics['mean_Ug_ml_s']:.6g} mL/s、Amin平均 {metrics['mean_Amin_mm2']:.6g}・最小 {metrics['min_Amin_mm2']:.6g} mm²。閉鎖率 {metrics['area_closure_fraction']:.6g}（有限のAmin<=0.01 mm²）。CQや確定f0は算出していない。
- 平均Pg {metrics['mean_Pg_Pa']:.6g} Pa、Psep {metrics['mean_Psep_Pa']:.6g} Pa、Pd0 {metrics['mean_Pd0_Pa']:.6g} Pa。
- xSep範囲 {metrics['xSep_range_mm']} mm、平均xSep {metrics['mean_xSep_mm']:.6g}・xBlendEnd {metrics['mean_xBlendEnd_mm']:.6g} mm。
- 監視点x平均 {metrics['mean_monitor_x_mm']:.6g} mm、最大増分開大点x平均 {metrics['mean_max_delta_gap_x_mm']:.6g} mm、剥離近傍プローブx平均 {metrics['mean_sep_node_x_mm']:.6g} mm。
- 平均位置差は監視点−剥離点 {metrics['monitor_minus_xSep_mean_mm']:.6g} mm、最大増分開大点−剥離点 {metrics['maximum_opening_minus_xSep_mean_mm']:.6g} mm。監視点がブレンド終端以上の区間割合は {metrics['monitor_beyond_blend_fraction']:.6g}。この条件では実節点圧力は現行境界値Pd0である。
- 平均ギャップ増分 {metrics['mean_delta_gap_mm']:.6g} mm、監視点uyL平均 {metrics['monitor_uyL_mean_mm']:.6g} mm・AC RMS {metrics['monitor_uyL_ac_rms_mm']:.6g} mm、剥離近傍uyL AC RMS {metrics['sep_node_uyL_ac_rms_mm']:.6g} mm。
- 剥離近傍プローブは各ステップでxSep,z=8.5 mmに最も近い節点を選び直す。固定材料点ではなく、RMSや波形には節点切替も含まれる。
- Ugスペクトルの上位局所ピーク {metrics['Ug_spectral_peaks_Hz_not_f0']} Hz。約0.1秒の有限評価窓、非閉鎖波形のピークであり、実験の133 Hzとの一致や周期的閉鎖の成立とは解釈しない。
- Lg範囲 {metrics['Lg_range']}、Rv範囲 {metrics['Rv_range']}、Rk範囲 {metrics['Rk_range']}。実ソルバで使用した局所変数を採取し、再積分を実係数と呼んでいない。閉鎖枝で未使用の係数はNaN+coefficients_active=0。
- Lg/Rvの終点は共通の選択剥離点。剥離後のA/G3を含めない。50断面は変形後の左右全体x範囲に配置し、AP閾値はない。部分AP断面が最狭位置を左右する近似は残る。
- pressure.datの剥離後保存値はPd0、実際の節点圧力はPsep→Pd0の0.5 mmブレンド。負圧分類には保存配列でなく実際の割当値を使用した。
- 実ランプは {meta['actual_ramp_time_s']} s。省略オプションの表示値0.1秒でなく、FlowModelのlegacy分岐による0.05秒。窓は0–0.05 ramp、0.05–0.10 transient、0.10–0.19999 evaluation。
- L_vt=0では現行回路がPd[i]=0、Ud[i]=Ugを代入することを確認しただけで、診断側から書き換えていない。PsepとPd0の不一致も無修正。

H5の局在・位置ずれは観測可能だが、幾何または圧力回復則が原因という因果は未判定。

![contact and geometry](figures/04_contact_and_gap.png)

## 仮説の順位・次の一条件

以下は追加検証の優先度であり、単一ランから確定した原因の順位ではない。

1. H1：{h1}。上記の正負合力と平均開大を根拠とする。
2. H5：負圧作用域・剥離・監視点の位置差と平均開大を定量化。因果の分離は次の1条件比較へ持ち越す。
3. H3：モードごとの仕事・位相を測定。閉鎖阻害の因果は未判定。
4. H4：AP端だけでなく中央にも接触を確認。接触局在と減衰は定量化したが、阻害の因果は未判定。
5. H2：{h2}。適用した荷重の符号と加法分解、投影/再構成を検査した範囲に限る。

次の比較案（未実行・別承認が必要）：{next_case}
単一観測ランでは原因の完全確定はできない。物理式、係数、並列順序、出力桁数は変更していない。

## 出力の定義と追加ファイル

必須CSV6件に加え、flow_geometry.csvは不足していた実係数・時刻同期プローブを1か所に集約したもの。
modal_native.binは全100モードの7列を各有効区間で一度だけ保存し、上位モードを再計算なしに抽出する固定長内部ファイル。
run_metadata.jsonには毎区間集計の極値・正負力平均、Gram、累積仕事検査を記録。analysis_metrics.jsonはレポートの機械可読な数値。
validation/には前回の境界バグ再現、ビルド/テスト監査、短時間比較、本番の入力・バイナリ識別子を集約。
本番の通常出力は作業worktreeのoutput/dat・wav・result以下だけ。元作業ツリーへ書き戻していない。
"""
    (OUT / "report.md").write_text(report)
    old = json.loads((OUT / "manifest.json").read_text())
    old.update(meta)
    old.update(status="complete", diagnostics_implemented=True, project_build_performed=True,
               input_copy_created=True, production_full_runs=1, no_physical_interpretation=False,
               production_attempt=attempt,
               validation=json.loads((OUT / "validation/source_and_test_audit.json").read_text()),
               production_provenance=json.loads((OUT / "validation/production_provenance.json").read_text()),
               current_original_worktree_preservation="verified unchanged", final_preservation_audit=preservation,
               diagnostic_build_default="OFF", diagnostic_runtime_default="OFF",
               source_ramp_reason="omitted pressureRampTimeSec uses historical 0.05 s; verified from source and runtime",
               excluded_diagnostics=["modal and nodal work", "nodal/modal consistency", "contact work", "contact pairs/forces/positions", "derived statistics"],
               integrated_interval_indices=[0, 19998], excluded_interval_indices=[19999],
               excluded_interval_time_s=[19999*dt, 20000*dt],
               actual_work_integrated_time_s=[0, 19999*dt], solver_steps_executed=20000,
               csv_stride={"force_regions": 20, "modal_selected_timeseries": 1, "contact_regions": 1, "consistency": 1, "flow_geometry": 1},
               units={"time": "s", "position": "mm", "nodal_velocity": "mm/s at t_n", "force": "N",
                      "pressure": "Pa", "applied_area": "m^2 (norm of actual applied area vector)",
                      "nodal_work": "N*m (J)", "modal_coordinates": "sampled-mass-normalized; each mode diagonal mass norm=1",
                      "modal_work": "algebraically equals nodal SI work; independent physical modal energy not guaranteed",
                      "Lg": "Pa*s^2/m^3 (=kg/m^4)", "Rv": "Pa*s/m^3 (=kg/(m^4*s))",
                      "Rk": "Pa*s^2/m^6 (=kg/m^7)"},
               vtu_timing="250 files: file index k corresponds to solver n=15000+20*k, coordinates after commitment at t_(n+1); original writer unchanged",
               source_hash_snapshot_stage="post-run; core source unmodified since validated build",
               source_sha256={str(p.relative_to(ROOT)): sha(p) for folder in ["src", "include", "tests"]
                              for p in sorted((ROOT / folder).rglob("*")) if p.is_file()},
               analysis_script_sha256=sha(Path(__file__)),
               analysis_dependencies={"numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
               selected_modes=selected, analysis_metrics="analysis_metrics.json")
    dump(OUT / "manifest.json", old)
    print(json.dumps(metrics, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["smoke", "production", "analyze"])
    parser.add_argument("--compare-only", action="store_true", help="compare completed smoke outputs without rerunning")
    parser.add_argument("--open-gap", action="store_true", help="additional open-flow smoke, initialGapMm=0.1")
    options = parser.parse_args()
    if options.command == "smoke":
        smoke(options.compare_only, options.open_gap)
    else:
        {"production": production, "analyze": analyze}[options.command]()
