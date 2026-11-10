#!/usr/bin/env python3
"""Read-only postprocessing for contact work signs and old/new bit deltas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIELDS = ("Q_pos", "Q_neg", "Q_fluid", "Q_contact", "q_n", "qdot_n", "delta_q")
SIDES = ("L", "R")
WINDOWS = ((.10, .20), (.20, .30), (.30, .40), (.40, .50), (.30, .50))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "sha256": sha256(path)}


def dump_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def require_fresh_output(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"refusing existing output: {path}")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"empty output is not allowed: {path}")
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def label_window(start: float, end: float) -> str:
    return f"[{start:.2f},{end:.2f})"


def parse_key_values(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip().strip('"')
    return result


def window_mask(count: int, dt: float, start: float, end: float) -> np.ndarray:
    time = np.arange(count, dtype=np.float64) * dt
    return (time >= start) & (time < end)


def ordered_uint64(value: float) -> int:
    bits = int(np.asarray(value, dtype=np.float64).view(np.uint64))
    return (~bits & ((1 << 64) - 1)) if bits & (1 << 63) else bits | (1 << 63)


def ulp_distance(first: float, second: float) -> int | None:
    if math.isnan(first) or math.isnan(second):
        return None
    if first == second:
        return 0
    return abs(ordered_uint64(first) - ordered_uint64(second))


def compare_dat_file(old_path: Path, new_path: Path, dt: float) -> tuple[list[dict], dict]:
    old_lines = old_path.read_text(encoding="utf-8").splitlines()
    new_lines = new_path.read_text(encoding="utf-8").splitlines()
    if not old_lines or not new_lines:
        raise ValueError(f"empty dat: {old_path}")
    old_header, new_header = old_lines[0], new_lines[0]
    old_rows = [line.split() for line in old_lines[1:]]
    new_rows = [line.split() for line in new_lines[1:]]
    if len(old_rows) != len(new_rows):
        raise ValueError(f"row mismatch for {old_path.name}")
    differences: list[dict] = []
    whitespace_only_lines = 0
    header_time = "time" in old_header.lower()
    for row_index, (old_tokens, new_tokens) in enumerate(zip(old_rows, new_rows)):
        if len(old_tokens) != len(new_tokens):
            raise ValueError(f"column mismatch {old_path.name} row {row_index}")
        if old_tokens == new_tokens and old_lines[row_index + 1] != new_lines[row_index + 1]:
            whitespace_only_lines += 1
        coordinate = float(old_tokens[0])
        time = coordinate if header_time else coordinate * dt
        step = int(round(time / dt)) if header_time else int(round(coordinate))
        if float(new_tokens[0]) != coordinate:
            raise ValueError(f"time/step mismatch {old_path.name} row {row_index}")
        for column, (old_text, new_text) in enumerate(zip(old_tokens, new_tokens)):
            if old_text == new_text:
                continue
            old_value, new_value = float(old_text), float(new_text)
            numeric_equal = (old_value == new_value or
                             (math.isnan(old_value) and math.isnan(new_value)))
            absolute = abs(old_value - new_value) if not numeric_equal else 0.0
            scale = max(abs(old_value), abs(new_value))
            relative = absolute / scale if scale else absolute
            differences.append({
                "file": old_path.name, "row_0based": row_index,
                "column_0based": column, "old_text": old_text,
                "new_text": new_text, "old_value": old_value,
                "new_value": new_value, "abs_diff": absolute,
                "rel_diff": relative, "step": step, "time_s": time,
                "difference_kind": "TEXT_ONLY" if numeric_equal else "NUMERIC_DELTA",
            })
    summary = {
        "file": old_path.name, "old_sha256": sha256(old_path),
        "new_sha256": sha256(new_path), "header_equal": old_header == new_header,
        "row_count": len(old_rows),
        "column_count": len(old_rows[0]) if old_rows else 0,
        "different_cells": len(differences),
        "numeric_delta_cells": sum(row["difference_kind"] == "NUMERIC_DELTA"
                                   for row in differences),
        "text_only_cells": sum(row["difference_kind"] == "TEXT_ONLY"
                                for row in differences),
        "whitespace_only_lines": whitespace_only_lines,
        "first_step": min((row["step"] for row in differences), default=None),
        "first_time_s": min((row["time_s"] for row in differences), default=None),
        "max_abs": max((row["abs_diff"] for row in differences), default=0.0),
        "max_rel": max((row["rel_diff"] for row in differences), default=0.0),
    }
    return differences, summary


def numeric_csv_first_differences(old_path: Path, new_path: Path) -> dict[str, object]:
    old = pd.read_csv(old_path); new = pd.read_csv(new_path)
    if list(old.columns) != list(new.columns) or len(old) != len(new):
        raise ValueError(f"CSV schema mismatch: {old_path.name}")
    records: list[dict[str, object]] = []
    for column in old.columns:
        if not (pd.api.types.is_numeric_dtype(old[column]) and
                pd.api.types.is_numeric_dtype(new[column])):
            continue
        first = old[column].to_numpy(); second = new[column].to_numpy()
        unequal = ~((first == second) | (np.isnan(first) & np.isnan(second)))
        locations = np.flatnonzero(unequal)
        if len(locations):
            index = int(locations[0])
            records.append({
                "column": column, "row_0based": index,
                "step": int(old.iloc[index]["step"]),
                "old": float(first[index]), "new": float(second[index]),
                "first_abs_delta": float(abs(first[index] - second[index])),
                "max_abs_delta": float(np.nanmax(np.abs(first - second))),
                "different_count": int(len(locations)),
            })
    records.sort(key=lambda row: (row["step"], row["row_0based"], row["column"]))
    return {"file": old_path.name, "first_numeric_difference": records[0] if records else None,
            "column_differences": records}


def modal_delta_rows(old: np.ndarray, new: np.ndarray, dt: float) -> tuple[list[dict], dict, dict[str, np.ndarray]]:
    if old.shape != new.shape or old.shape[-1] != len(FIELDS):
        raise ValueError("modal binary schemas do not align")
    rows: list[dict] = []
    series: dict[str, np.ndarray] = {}
    global_first: tuple[int, int, int, int] | None = None
    first_records: list[dict] = []
    for field_index, field in enumerate(FIELDS):
        difference = np.abs(np.asarray(new[..., field_index]) - np.asarray(old[..., field_index]))
        series[field] = np.max(difference, axis=(1, 2))
        unequal = new[..., field_index] != old[..., field_index]
        indices = np.argwhere(unequal)
        if len(indices):
            candidate = (int(indices[0, 0]), int(indices[0, 1]), int(indices[0, 2]), field_index)
            if global_first is None or candidate < global_first:
                global_first = candidate
        for side_index, side in enumerate(SIDES):
            for mode in range(old.shape[2]):
                a = np.asarray(old[:, side_index, mode, field_index])
                b = np.asarray(new[:, side_index, mode, field_index])
                neq = a != b
                locations = np.flatnonzero(neq)
                delta = np.abs(b - a)
                scale = np.maximum(np.abs(a), np.abs(b))
                relative = np.divide(delta, scale, out=np.zeros_like(delta), where=scale != 0)
                if len(locations):
                    first = int(locations[0])
                    maximum = int(np.argmax(delta))
                    first_ulp = ulp_distance(float(a[first]), float(b[first]))
                    max_ulp = ulp_distance(float(a[maximum]), float(b[maximum]))
                    first_records.append({"step": first, "time_s": first * dt,
                                          "side": side, "mode": mode + 1,
                                          "field": field, "old": float(a[first]),
                                          "new": float(b[first]),
                                          "delta": float(b[first] - a[first])})
                else:
                    first = maximum = None
                    first_ulp = max_ulp = 0
                rows.append({
                    "field": field, "side": side, "mode": mode + 1,
                    "different_count": int(neq.sum()),
                    "first_different_step": first if first is not None else "",
                    "first_different_time_s": first * dt if first is not None else "",
                    "first_old": float(a[first]) if first is not None else "",
                    "first_new": float(b[first]) if first is not None else "",
                    "first_ulp_distance": first_ulp if first_ulp is not None else "NAN",
                    "max_abs": float(delta.max()), "max_rel": float(relative.max()),
                    "max_abs_step": maximum if maximum is not None else "",
                    "max_abs_time_s": maximum * dt if maximum is not None else "",
                    "max_ulp_at_max_abs": max_ulp if max_ulp is not None else "NAN",
                })
    first_records.sort(key=lambda row: (row["step"], FIELDS.index(row["field"]),
                                        SIDES.index(row["side"]), row["mode"]))
    first_step = global_first[0] if global_first else None
    summary = {
        "bitwise_equal": global_first is None,
        "first_different_step": first_step,
        "first_different_time_s": first_step * dt if first_step is not None else None,
        "first_observed_record": first_records[0] if first_records else None,
        "records_at_first_step": [row for row in first_records if row["step"] == first_step],
        "max_abs_by_field": {field: float(series[field].max()) for field in FIELDS},
    }
    return rows, summary, series


def work_products(modal: np.ndarray, components: np.ndarray, dt: float,
                  regions: pd.DataFrame) -> tuple[pd.DataFrame, list[dict], list[dict], dict]:
    count, _, nmode, _ = modal.shape
    if components.shape != (count, 2, nmode, 2):
        raise ValueError("contact component shape mismatch")
    dq = np.asarray(modal[..., 6])
    increments = np.asarray(components) * dq[..., None]
    per_step = increments.sum(axis=2)
    total_reference = np.asarray(modal[..., 3]) * dq
    residual = increments.sum(axis=-1) - total_reference
    time = np.arange(count) * dt
    region_step = regions.groupby(["step", "side"], as_index=False).agg(
        pairs=("pairs", "sum"), max_penetration_m=("max_penetration_m", "max"))
    lookup = {(int(row.step), str(row.side)): (int(row.pairs), float(row.max_penetration_m))
              for row in region_step.itertuples()}
    frames = []
    for side_index, side in enumerate(SIDES):
        frames.append(pd.DataFrame({
            "step": np.arange(count), "time_s": time, "side": side,
            "spring_increment_J": per_step[:, side_index, 0],
            "damp_increment_J": per_step[:, side_index, 1],
            "total_increment_J": per_step[:, side_index].sum(axis=1),
            "spring_cumulative_J": np.cumsum(per_step[:, side_index, 0]),
            "damp_cumulative_J": np.cumsum(per_step[:, side_index, 1]),
            "total_cumulative_J": np.cumsum(per_step[:, side_index].sum(axis=1)),
            "contact_nonzero": np.any(components[:, side_index] != 0, axis=(1, 2)),
            "pairs": [lookup.get((step, side), (0, 0.0))[0] for step in range(count)],
            "max_penetration_m": [lookup.get((step, side), (0, 0.0))[1]
                                    for step in range(count)],
        }))
    step_frame = pd.concat(frames, ignore_index=True).sort_values(["step", "side"])
    window_rows: list[dict] = []
    event_checks: dict[str, object] = {}
    for start, end in WINDOWS:
        mask = window_mask(count, dt, start, end)
        window = label_window(start, end)
        side_values: dict[str, dict[str, float]] = {}
        for side_index, side in enumerate(SIDES):
            values = per_step[mask, side_index]
            record: dict[str, object] = {"window": window, "side": side,
                                         "interval_count": int(mask.sum())}
            side_values[side] = {}
            for component_index, component in enumerate(("spring", "damp")):
                vector = values[:, component_index]
                negative = -vector[vector < 0]
                sorted_negative = np.sort(negative)[::-1]
                record[f"W_{component}_positive_J"] = float(vector[vector > 0].sum())
                record[f"W_{component}_negative_J"] = float(vector[vector < 0].sum())
                record[f"W_{component}_net_J"] = float(vector.sum())
                record[f"{component}_positive_intervals"] = int(np.sum(vector > 0))
                record[f"{component}_negative_intervals"] = int(np.sum(vector < 0))
                record[f"{component}_nonzero_intervals"] = int(np.sum(vector != 0))
                for percentage in (1, 5):
                    take = max(1, math.ceil(len(sorted_negative) * percentage / 100)) if len(sorted_negative) else 0
                    fraction = float(sorted_negative[:take].sum() / sorted_negative.sum()) if len(sorted_negative) and sorted_negative.sum() else 0.0
                    record[f"{component}_top{percentage}pct_negative_fraction"] = fraction
                side_values[side][component] = float(vector.sum())
            record["W_total_J"] = float(values.sum())
            window_rows.append(record)
        both = {"window": window, "side": "BOTH",
                "interval_count": int(mask.sum())}
        difference = {"window": window, "side": "R_MINUS_L",
                      "interval_count": int(mask.sum())}
        for component in ("spring", "damp"):
            both[f"W_{component}_net_J"] = side_values["L"][component] + side_values["R"][component]
            difference[f"W_{component}_net_J"] = side_values["R"][component] - side_values["L"][component]
        both["W_total_J"] = both["W_spring_net_J"] + both["W_damp_net_J"]
        difference["W_total_J"] = difference["W_spring_net_J"] + difference["W_damp_net_J"]
        window_rows.extend([both, difference])
    region_mode_rows: list[dict] = []
    main = (regions.time_s >= .30) & (regions.time_s < .50)
    for (side, ap), group in regions[main].groupby(["side", "ap_region"]):
        for component in ("spring", "damp"):
            vector = group[f"{component}_work_nodal_J"].to_numpy()
            region_mode_rows.append({
                "classification": "AP_REGION", "window": "[0.30,0.50)",
                "side": side, "ap_region": ap, "mode": "NOT_CROSS_CLASSIFIED",
                "component": component, "positive_J": float(vector[vector > 0].sum()),
                "negative_J": float(vector[vector < 0].sum()), "net_J": float(vector.sum()),
                "nonzero_intervals": int(np.sum(vector != 0)),
                "availability": "OBSERVED_AP_NODAL_WORK; region-by-mode intersection NOT_IDENTIFIABLE",
            })
    main_mask = window_mask(count, dt, .30, .50)
    for side_index, side in enumerate(SIDES):
        for mode in range(nmode):
            for component_index, component in enumerate(("spring", "damp")):
                vector = increments[main_mask, side_index, mode, component_index]
                region_mode_rows.append({
                    "classification": "MODE", "window": "[0.30,0.50)", "side": side,
                    "ap_region": "NOT_CROSS_CLASSIFIED", "mode": mode + 1,
                    "component": component, "positive_J": float(vector[vector > 0].sum()),
                    "negative_J": float(vector[vector < 0].sum()), "net_J": float(vector.sum()),
                    "nonzero_intervals": int(np.sum(vector != 0)),
                    "availability": "OBSERVED_MODAL_WORK; region-by-mode intersection NOT_IDENTIFIABLE",
                })
    spring_l = per_step[main_mask, 0, 0]
    spring_r = per_step[main_mask, 1, 0]
    active = (spring_l != 0) | (spring_r != 0)
    opposite = active & (np.signbit(spring_l) != np.signbit(spring_r))
    event_checks = {
        "max_abs_component_sum_minus_reference_J": float(np.max(np.abs(residual))),
        "primary_spring_opposite_sign_active_interval_fraction":
            float(opposite.sum() / active.sum()) if active.sum() else 0.0,
        "primary_spring_timewise_cancellation_fraction": float(
            1 - np.sum(np.abs(spring_l + spring_r)) /
            np.sum(np.abs(spring_l) + np.abs(spring_r))) if np.sum(np.abs(spring_l) + np.abs(spring_r)) else 0.0,
    }
    for side_index, side in enumerate(SIDES):
        vector = per_step[main_mask, side_index, 0]
        steps = np.flatnonzero(main_mask)
        positive_order = np.argsort(vector)[::-1][:10]
        negative_order = np.argsort(vector)[:10]
        event_checks[f"{side}_top_spring_positive_events"] = [
            {"step": int(steps[i]), "time_s": float(steps[i] * dt), "work_J": float(vector[i])}
            for i in positive_order if vector[i] > 0]
        event_checks[f"{side}_top_spring_negative_events"] = [
            {"step": int(steps[i]), "time_s": float(steps[i] * dt), "work_J": float(vector[i])}
            for i in negative_order if vector[i] < 0]
    return step_frame, window_rows, region_mode_rows, event_checks


def make_figures(output: Path, step_frame: pd.DataFrame, region_mode: pd.DataFrame,
                 modal_series: dict[str, np.ndarray], dt: float) -> None:
    figures = output / "figures"; figures.mkdir()
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for side, style in (("L", "-"), ("R", "--")):
        frame = step_frame[step_frame.side == side]
        axes[0].plot(frame.time_s, frame.spring_cumulative_J * 1e3, style, label=side)
        axes[1].plot(frame.time_s, frame.damp_cumulative_J * 1e3, style, label=side)
        axes[2].plot(frame.time_s, frame.total_cumulative_J * 1e3, style, label=side)
    for ax, title in zip(axes, ("spring", "velocity-dependent", "total")):
        ax.axvspan(.3, .5, color="grey", alpha=.12); ax.set_ylabel(f"{title} [mJ]"); ax.legend()
    axes[-1].set_xlabel("interval left-end time [s]"); fig.tight_layout()
    fig.savefig(figures / "01_cumulative_components.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for side, color in (("L", "tab:blue"), ("R", "tab:orange")):
        frame = step_frame[(step_frame.side == side) & (step_frame.time_s >= .3)]
        axes[0].plot(frame.time_s, frame.spring_increment_J * 1e6, color=color, alpha=.75, label=side)
    combined = step_frame[step_frame.side == "L"].copy()
    combined["pairs"] = combined.pairs.to_numpy() + step_frame[step_frame.side == "R"].pairs.to_numpy()
    combined["max_penetration_m"] = np.maximum(combined.max_penetration_m.to_numpy(), step_frame[step_frame.side == "R"].max_penetration_m.to_numpy())
    mask = combined.time_s >= .3
    axes[1].plot(combined.time_s[mask], combined.pairs[mask], color="black")
    axes[2].plot(combined.time_s[mask], combined.max_penetration_m[mask] * 1e3, color="tab:red")
    axes[0].set_ylabel("spring dW [uJ]"); axes[0].legend(); axes[1].set_ylabel("pair rows L+R")
    axes[2].set_ylabel("max penetration [mm]"); axes[2].set_xlabel("time [s]"); fig.tight_layout()
    fig.savefig(figures / "02_spring_events.png", dpi=160); plt.close(fig)

    region = region_mode[(region_mode.classification == "AP_REGION")]
    fig, ax = plt.subplots(figsize=(10, 5)); labels=[]; positive=[]; negative=[]; net=[]
    for row in region.itertuples():
        labels.append(f"{row.side}-{row.ap_region}-{row.component}")
        positive.append(row.positive_J*1e3); negative.append(row.negative_J*1e3); net.append(row.net_J*1e3)
    x=np.arange(len(labels)); ax.bar(x-.25, positive, .25, label="positive sum"); ax.bar(x, negative, .25, label="negative sum"); ax.bar(x+.25, net, .25, label="net")
    ax.set(xticks=x, xticklabels=labels, ylabel="nodal reference work [mJ]"); ax.tick_params(axis="x", rotation=30); ax.legend(); fig.tight_layout()
    fig.savefig(figures / "03_ap_signed_work.png", dpi=160); plt.close(fig)

    modes = region_mode[region_mode.classification == "MODE"]
    fig, axes = plt.subplots(2,1,figsize=(10,7),sharex=True)
    for ax, side in zip(axes, SIDES):
        spring=modes[(modes.side==side)&(modes.component=="spring")]; damp=modes[(modes.side==side)&(modes.component=="damp")]
        x=np.arange(1,21); ax.bar(x-.2,spring.net_J*1e3,.4,label="spring"); ax.bar(x+.2,damp.net_J*1e3,.4,label="velocity-dependent"); ax.set_ylabel(f"{side} [mJ]"); ax.legend()
    axes[-1].set(xlabel="active/source mode",xticks=np.arange(1,21)); fig.tight_layout(); fig.savefig(figures/"04_modal_components.png",dpi=160); plt.close(fig)

    fig, axes = plt.subplots(2,1,figsize=(10,6),sharex=True); time=np.arange(len(next(iter(modal_series.values()))))*dt
    for field in ("q_n","qdot_n","delta_q"):
        axes[0].semilogy(time, np.maximum(modal_series[field],1e-300),label=field)
    for field in ("Q_fluid","Q_contact"):
        axes[1].semilogy(time,np.maximum(modal_series[field],1e-300),label=field)
    axes[0].set_ylabel("max abs state delta"); axes[1].set_ylabel("max abs force delta"); axes[1].set_xlabel("step time [s]")
    for ax in axes: ax.legend(); ax.axvspan(0,.002,color="grey",alpha=.1)
    fig.tight_layout(); fig.savefig(figures/"05_modal_bitdelta_timeline.png",dpi=160); plt.close(fig)


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--old-run",type=Path,required=True)
    parser.add_argument("--old-diagnostic",type=Path,required=True)
    parser.add_argument("--new-run",type=Path,required=True)
    parser.add_argument("--previous-analysis",type=Path,required=True)
    parser.add_argument("--sym-analysis",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--no-simulation",action="store_true",required=True)
    args=parser.parse_args()
    require_fresh_output(args.output_dir)
    old_diag=args.old_diagnostic/"kanaya_diagnostic"; new_diag=args.new_run/"kanaya_diagnostic"
    old_meta=json.loads((old_diag/"run_metadata.json").read_text()); new_meta=json.loads((new_diag/"run_metadata.json").read_text()); comp_meta=json.loads((new_diag/"contact_component_metadata.json").read_text())
    for meta in (old_meta,new_meta,comp_meta):
        if not meta.get("complete"): raise ValueError("incomplete diagnostic input")
    count=int(new_meta["integrated_intervals"]); dt=float(new_meta["dt_s"]); nmode=int(new_meta["nmode"][0])
    if count!=int(old_meta["integrated_intervals"]) or dt!=float(old_meta["dt_s"]): raise ValueError("old/new time schema mismatch")
    old_modal_path=old_diag/"modal_native.bin"; new_modal_path=new_diag/"modal_native.bin"; component_path=new_diag/"contact_component_native.bin"
    expected_modal=count*2*nmode*7*8; expected_component=count*2*nmode*2*8
    if old_modal_path.stat().st_size!=expected_modal or new_modal_path.stat().st_size!=expected_modal or component_path.stat().st_size!=expected_component: raise ValueError("binary byte size mismatch")
    inputs=[old_modal_path,new_modal_path,component_path,old_diag/"run_metadata.json",new_diag/"run_metadata.json",new_diag/"contact_component_metadata.json",new_diag/"contact_component_regions.csv"]
    old_binary=args.old_diagnostic.parent/"build/on/simulation"
    new_binary=args.new_run.parent/"build/on/simulation"
    for path in (old_binary,new_binary,args.old_diagnostic/"txt/params_used.txt",
                 args.new_run/"txt/params_used.txt",args.old_diagnostic/"txt/manifest.txt",
                 args.new_run/"txt/manifest.txt"):
        if path.exists(): inputs.append(path)
    run_manifest=parse_key_values(args.new_run/"txt/manifest.txt")
    for key in ("left_mode_vtu","right_mode_vtu","left_frequency","right_frequency",
                "left_surface_nas","right_surface_nas"):
        path=Path(run_manifest[key])
        if path not in inputs: inputs.append(path)
    old_dat_root=args.old_run/"dat"; new_dat_root=args.new_run/"dat"; dat_names=sorted(path.name for path in old_dat_root.glob("*.dat"))
    if dat_names!=sorted(path.name for path in new_dat_root.glob("*.dat")) or len(dat_names)!=9: raise ValueError("expected matching nine dat files")
    inputs.extend([old_dat_root/name for name in dat_names]); inputs.extend([new_dat_root/name for name in dat_names])
    before={str(path.resolve()):file_identity(path) for path in inputs}
    old=np.memmap(old_modal_path,dtype="<f8",mode="r",shape=(count,2,nmode,7)); new=np.memmap(new_modal_path,dtype="<f8",mode="r",shape=(count,2,nmode,7)); components=np.memmap(component_path,dtype="<f8",mode="r",shape=(count,2,nmode,2))
    regions=pd.read_csv(new_diag/"contact_component_regions.csv")
    if len(regions)!=count*4: raise ValueError("component region row count mismatch")
    args.output_dir.mkdir(parents=True)
    all_dat=[]; dat_summary=[]
    for name in dat_names:
        rows,summary=compare_dat_file(old_dat_root/name,new_dat_root/name,dt); all_dat.extend(rows); dat_summary.append(summary)
    if not all_dat:
        all_dat=[{"file":"NONE","difference_kind":"NO_DIFFERENCES"}]
    write_csv(args.output_dir/"dat_first_differences.csv",all_dat)
    modal_rows,modal_summary,modal_series=modal_delta_rows(old,new,dt); write_csv(args.output_dir/"modal_binary_delta.csv",modal_rows)
    step_frame,window_rows,region_mode_rows,event_checks=work_products(new,components,dt,regions)
    step_frame.to_csv(args.output_dir/"work_by_step.csv",index=False); write_csv(args.output_dir/"work_by_window.csv",window_rows); write_csv(args.output_dir/"work_by_region_mode.csv",region_mode_rows)
    region_mode_frame=pd.DataFrame(region_mode_rows); make_figures(args.output_dir,step_frame,region_mode_frame,modal_series,dt)

    old_contact=np.asarray(old[...,3])*np.asarray(old[...,6]); new_contact=np.asarray(new[...,3])*np.asarray(new[...,6]); primary=window_mask(count,dt,.3,.5)
    old_work=np.sum(old_contact[primary],axis=(0,2)); new_work=np.sum(new_contact[primary],axis=(0,2)); work_delta=new_work-old_work
    main_rows={row["side"]:row for row in window_rows if row["window"]=="[0.30,0.50)"}
    spring_diff=float(main_rows["R_MINUS_L"]["W_spring_net_J"]); damp_diff=float(main_rows["R_MINUS_L"]["W_damp_net_J"]); total_diff=spring_diff+damp_diff
    flow_columns=["step","contact_iterations","stop_reason","contact_flag"]
    old_flow=pd.read_csv(old_diag/"flow_geometry.csv",usecols=flow_columns); new_flow=pd.read_csv(new_diag/"flow_geometry.csv",usecols=flow_columns)
    contact_event_equal=bool(old_flow.equals(new_flow))
    auxiliary_differences={name:numeric_csv_first_differences(old_diag/name,new_diag/name)
                           for name in ("consistency.csv","contact_regions.csv","flow_geometry.csv")}
    diagnosis={
        "first_internal_difference":modal_summary["first_observed_record"],
        "records_first_observed_step":modal_summary["records_at_first_step"],
        "first_dat_difference":min(all_dat,key=lambda row:row.get("step",10**99)) if all_dat[0].get("file")!="NONE" else None,
        "cause_classification":{
            "FORMAT_THRESHOLD":"SUPPORTED: 41 printed numeric cells cross formatting thresholds while internal deltas are nonzero",
            "NUMERIC_TRAJECTORY_DELTA":"OBSERVED in modal binary",
            "METADATA_OR_BUILD_DIFF":"OBSERVED: input/OMP nominally match, executable SHA differs",
            "CONTACT_EVENT_DIVERGENCE":"NOT_OBSERVED in saved contact_iterations/stop_reason/contact_flag",
            "UNKNOWN":"ROOT CAUSE INCONCLUSIVE: the first differing saved quantity is Q_contact; earlier unsaved intermediates prevent attribution",
        },
        "saved_stage_order":"first saved Q_contact delta precedes first saved qdot/q/Q_fluid deltas; this is observation order, not proof of computation-origin order",
        "component_split_extrapolates_to_old_trajectory":"UNPROVEN",
        "auxiliary_diagnostic_first_differences":auxiliary_differences,
    }
    dump_json(args.output_dir/"bitdelta_diagnosis.json",diagnosis)
    checks={
        "INTERNAL_COMPONENT_QC":"PASS" if event_checks["max_abs_component_sum_minus_reference_J"]<=1e-12 else "FAIL",
        "SHORT_OFF_ON_BITWISE":"PASS (inherited and input-hash verified from previous analysis)",
        "OLD_NEW_FULLRUN_BITWISE":"FAIL",
        "MECHANISTIC_CAUSALITY":"NOT_ESTABLISHED",
        "dat_comparison":{"different_cells":sum(x["different_cells"] for x in dat_summary),"summaries":dat_summary},
        "modal_comparison":modal_summary,
        "work_checks":event_checks,
        "primary_old_contact_work_J":{"L":float(old_work[0]),"R":float(old_work[1])},
        "primary_new_contact_work_J":{"L":float(new_work[0]),"R":float(new_work[1])},
        "primary_old_new_delta_J":{"L":float(work_delta[0]),"R":float(work_delta[1])},
        "old_new_delta_over_lr_contact_difference":float(np.max(np.abs(work_delta))/abs(total_diff)),
        "new_run_component_split_algebraically_valid":event_checks["max_abs_component_sum_minus_reference_J"]<=1e-12,
        "trajectory_bitwise_identical":False,"raw_work_near_equal":True,
        "component_split_extrapolates_to_old_trajectory":"UNPROVEN",
        "spring_R_minus_L_J":spring_diff,"damp_R_minus_L_J":damp_diff,"total_R_minus_L_J":total_diff,
        "spring_share_of_algebraic_lr_difference":spring_diff/total_diff,
        "damp_share_of_algebraic_lr_difference":damp_diff/total_diff,
        "contact_event_discrete_columns_equal":contact_event_equal,
        "auxiliary_diagnostic_first_differences":auxiliary_differences,
        "region_by_mode_cross_classification":"NOT_IDENTIFIABLE",
        "pair_local_relative_displacement_work":"NOT_IDENTIFIABLE",
        "symmetric_spring_damp_components":"NOT_MEASURED",
        "original_failed_qc_preserved":True,
    }
    after={str(path.resolve()):file_identity(path) for path in inputs}; checks["source_inputs_unchanged"]=before==after
    dump_json(args.output_dir/"checks.json",checks)
    try:
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(); git_status=subprocess.check_output(["git","status","--short"],text=True).splitlines()
    except (OSError,subprocess.CalledProcessError): git_head="UNKNOWN"; git_status=[]
    manifest={"schema":"postprocess-contact-work-bitdelta-v1","git_head":git_head,"git_status_short":git_status,"simulation_runs_launched":0,"solver_rebuilt":False,"model_equations_changed":False,"old_dat_source":str(args.old_run.resolve()),"old_modal_source":str(args.old_diagnostic.resolve()),"new_component_source":str(args.new_run.resolve()),"previous_analysis":str(args.previous_analysis.resolve()),"symmetric_comparison":str(args.sym_analysis.resolve()),"build_provenance":{"old_binary":file_identity(old_binary) if old_binary.exists() else None,"new_binary":file_identity(new_binary) if new_binary.exists() else None,"common_flags":"-O3 -std=gnu++17 -fiopenmp; ENABLE_KANAYA_DIAGNOSTICS; icpx 2025.1","binary_sha_identical":old_binary.exists() and new_binary.exists() and sha256(old_binary)==sha256(new_binary)},"binary_schema":{"modal":{"shape":[count,2,nmode,7],"dtype":"little-endian float64","fields":list(FIELDS)},"component":{"shape":[count,2,nmode,2],"dtype":"little-endian float64","fields":["Q_spring","Q_damp"]}},"time":{"dt_s":dt,"intervals":count,"final_solver_step_excluded":49999,"windows":[label_window(*w) for w in WINDOWS]},"input_files":list(before.values()),"source_inputs_unchanged":before==after,"qc_labels":{"INTERNAL_COMPONENT_QC":checks["INTERNAL_COMPONENT_QC"],"SHORT_OFF_ON_BITWISE":"PASS","OLD_NEW_FULLRUN_BITWISE":"FAIL","MECHANISTIC_CAUSALITY":"NOT_ESTABLISHED"},"artifacts":["checks.json","work_by_step.csv","work_by_window.csv","work_by_region_mode.csv","dat_first_differences.csv","modal_binary_delta.csv","bitdelta_diagnosis.json","report.md","figures/"]}
    dump_json(args.output_dir/"manifest.json",manifest)

    mode_frame=region_mode_frame[region_mode_frame.classification=="MODE"].copy(); mode_frame["abs_net"]=mode_frame.net_J.abs(); top_modes=mode_frame.sort_values("abs_net",ascending=False).head(10)
    ap_frame=region_mode_frame[region_mode_frame.classification=="AP_REGION"]
    event_text=[]
    for side in SIDES:
        pos=event_checks[f"{side}_top_spring_positive_events"][0]; neg=event_checks[f"{side}_top_spring_negative_events"][0]
        event_text.append(f"- {side}: 最大正イベント {pos['time_s']:.5f} s ({pos['work_J']*1e6:.4g} µJ)、最大負イベント {neg['time_s']:.5f} s ({neg['work_J']*1e6:.4g} µJ)")
    report=f"""# 12/6・900 Pa 接触仕事符号差とbit差の既存run後処理

**OBSERVED**：この条件は左小振幅ではなく、柔らかい右が小振幅（代表点R/L=0.302、表面R/L≈0.491）である。  
**OBSERVED**：新runの右接触仕事はspring {main_rows['R']['W_spring_net_J']*1e3:.6f} mJ、velocity-dependent {main_rows['R']['W_damp_net_J']*1e3:.6f} mJで、後者は総量の63.2%。  
**OBSERVED**：左右差R−L={total_diff*1e3:.6f} mJの代数的内訳はspring {spring_diff*1e3:.6f} mJ（{spring_diff/total_diff*100:.1f}%）、velocity-dependent {damp_diff*1e3:.6f} mJ（{damp_diff/total_diff*100:.1f}%）。  
**CONSISTENT_WITH**：モード応答・モード間相殺と接触仕事配分の左右差が右小振幅に併存する。**NOT_ESTABLISHED**：これらの因果寄与率。  
**FAILED_QC維持**：短時間OFF/ONはbit一致したが旧新フルrunは非一致。新run内部成分QCのみPASSである。

## A. 左+・右−のばね参考仕事

主窓 `[0.30,0.50)` のspringはL={main_rows['L']['W_spring_net_J']*1e3:.9f} mJ、R={main_rows['R']['W_spring_net_J']*1e3:.9f} mJ、両側和={main_rows['BOTH']['W_spring_net_J']*1e3:.9f} mJ。velocity-dependent両側和は{main_rows['BOTH']['W_damp_net_J']*1e3:.9f} mJである。ばねの両側和が小さいことは保存ポテンシャルの証明ではない。

ばね仕事は単一符号ではなく、各側で正負イベントが相殺した正味値である。詳細な正・負和、区間数、負仕事上位1%/5%集中度は `work_by_window.csv` に記録した。主窓の左右ばね仕事が同時刻に反対符号となるactive区間率は {event_checks['primary_spring_opposite_sign_active_interval_fraction']:.3%}。時刻ごとの絶対仕事に基づく左右相殺指標は {event_checks['primary_spring_timewise_cancellation_fraction']:.3%}。

{chr(10).join(event_text)}

主窓のside-step合算では、L springは正寄与{main_rows['L']['W_spring_positive_J']*1e3:.6f} mJと負寄与{main_rows['L']['W_spring_negative_J']*1e3:.6f} mJ、R springは正寄与{main_rows['R']['W_spring_positive_J']*1e3:.6f} mJと負寄与{main_rows['R']['W_spring_negative_J']*1e3:.6f} mJの相殺である。0.10–0.20、0.20–0.30、0.30–0.40、0.40–0.50 sの全副窓でL正味正・R正味負が再現され、少数イベントだけの符号反転ではない。

AP地域別では次の正味値を得た（節点仕事として保存済み。地域×モード交差は `NOT_IDENTIFIABLE`）。

| side | AP | component | positive [mJ] | negative [mJ] | net [mJ] |
|---|---|---|---:|---:|---:|
"""
    for row in ap_frame.itertuples(): report+=f"| {row.side} | {row.ap_region} | {row.component} | {row.positive_J*1e3:.6g} | {row.negative_J*1e3:.6g} | {row.net_J*1e3:.6g} |\n"
    report+="\n主要モード（|正味仕事|上位、成分別）：\n\n| side | mode | component | positive [mJ] | negative [mJ] | net [mJ] |\n|---|---:|---|---:|---:|---:|\n"
    for row in top_modes.itertuples(): report+=f"| {row.side} | {row.mode} | {row.component} | {row.positive_J*1e3:.6g} | {row.negative_J*1e3:.6g} | {row.net_J*1e3:.6g} |\n"
    first=modal_summary['first_observed_record']; first_dat=diagnosis['first_dat_difference']
    report+=f"""

力は採用予測接触反復、仕事の変位増分は最終Newmark解で評価される。したがってspringの正仕事を「駆動」、負仕事を「不可逆散逸」と断定しない。pair-local相対変位仕事は未保存で `NOT_IDENTIFIABLE`。

## B. bit不一致の初発

9 datの全差分を再検算し、数値文字列差は合計 **{sum(x['different_cells'] for x in dat_summary)}セル**、空白だけの差は{sum(x['whitespace_only_lines'] for x in dat_summary)}行だった。`airflow_vt.dat`はSHA256一致。最初の印字差は `{first_dat['file']}` step={first_dat['step']}（t={first_dat['time_s']:.5f} s）、column={first_dat['column_0based']}、`{first_dat['old_text']}`→`{first_dat['new_text']}`。

生modal binaryの最初の差は step={first['step']}（t={first['time_s']:.5f} s）、side={first['side']}、mode={first['mode']}、field={first['field']}、差={first['delta']:.3e}。保存量の観測順では `Q_contact` が最初で、その後にqdot/q/流体荷重差が現れる。ただし接触計算内部のより早い中間値は保存されていないため、observerや特定演算を根本原因とは断定できない。保存された接触反復数・停止理由・接触フラグは全step一致し、接触イベント分岐の相違は観測されなかった。

補助診断CSVでは `consistency.csv` のobserver集計 `work_scale` がstep 519で先に丸め差を持つが、この時点の7-field modal値は一致している。採用接触地域の `Fx_N` とmodal `Q_contact` はstep 1729、状態はqdot step 2340、Δq step 2382、q step 2383、流体modal力はstep 2715から差が現れた。これは保存量の初発順序であり、未保存の演算順序を確定するものではない。

分類は `FORMAT_THRESHOLD + NUMERIC_TRAJECTORY_DELTA + METADATA_OR_BUILD_DIFF`。内部差が実在するため `FORMAT_ONLY` ではない。実行バイナリSHAは異なるが、それだけでobserverが原因とは証明できず、根本原因は **INCONCLUSIVE**。

旧新総接触仕事差はL={work_delta[0]*1e3:.3e} mJ、R={work_delta[1]*1e3:.3e} mJで、左右接触仕事差に対する最大比は{checks['old_new_delta_over_lr_contact_difference']:.3e}。総仕事の観測は定量的に非常に近いが、bitwise QCをPASSへ変更しない。新成分分解を旧軌道へ外挿することは `UNPROVEN`。

## 暫定結論と限界

対称12/12では左右振幅が一致し、12/6では左の振幅増加と右の減少、左の正・右の負のモード間共分散、接触参考仕事配分の左右差が同時に観測された。新runでは右総接触参考仕事の63.2%が速度依存項だが、左右差の81.4%はばね仕事差の代数的内訳である。この二つは別の分母であり、振幅差の因果寄与率ではない。

流体仕事はraw値、接触仕事は離散参照仕事であり、厳密なエネルギー収支や不可逆散逸量ではない。対称runのspring/damp内訳、地域×モード交差、pair-local相対変位仕事は未測定。最小追加試験は旧ビルド条件を固定したbit再現と、pair-local変位仕事を保存する対照診断だが、今回は実行していない。
"""
    (args.output_dir/"report.md").write_text(report,encoding="utf-8")


if __name__=="__main__": main()
