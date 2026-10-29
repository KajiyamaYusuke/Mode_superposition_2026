#!/usr/bin/env python3
"""Summarize actual applied nodal fluid loads by pressure sign and region."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


WINDOWS = ((0.10, 0.20), (0.20, 0.30), (0.30, 0.40),
           (0.40, 0.50), (0.30, 0.50))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    source = pd.read_csv(arguments.input)
    rows = []
    keys = ["step", "time_s", "side", "zone", "ap_region", "pressure_sign"]
    per_step = source.groupby(keys, as_index=False).agg(
        Fx_N=("Fx_N", "sum"),
        Fy_N=("Fy_N", "sum"),
        Fz_N=("Fz_N", "sum"),
        opening_Fy_N=("opening_Fy_N", "sum"),
        applied_node_count=("applied_node_count", "sum"),
        area_m2=("area_m2", "sum"),
        wrong_normal_count=("wrong_normal_count", "sum"),
    )
    for start, end in WINDOWS:
        selected = per_step[(per_step.time_s >= start) & (per_step.time_s < end)]
        for group_keys, group in selected.groupby(
            ["side", "zone", "ap_region", "pressure_sign"], sort=True
        ):
            side, zone, ap_region, pressure_sign = group_keys
            rows.append(
                {
                    "window_start_s": start,
                    "window_end_s": end,
                    "side": side,
                    "zone": zone,
                    "ap_region": ap_region,
                    "pressure_sign": pressure_sign,
                    "sampled_steps": len(group),
                    "mean_Fx_N": group.Fx_N.mean(),
                    "mean_Fy_N": group.Fy_N.mean(),
                    "mean_Fz_N": group.Fz_N.mean(),
                    "mean_opening_Fy_N": group.opening_Fy_N.mean(),
                    "rms_opening_Fy_N": np.sqrt(np.mean(group.opening_Fy_N ** 2)),
                    "min_opening_Fy_N": group.opening_Fy_N.min(),
                    "max_opening_Fy_N": group.opening_Fy_N.max(),
                    "mean_applied_node_count": group.applied_node_count.mean(),
                    "mean_area_m2": group.area_m2.mean(),
                    "wrong_normal_count_total": int(group.wrong_normal_count.sum()),
                    "source_stride_steps": 20,
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(arguments.output, index=False)
    primary = result[(result.window_start_s == 0.30) & (result.window_end_s == 0.50)]
    print(
        primary.groupby(["side", "pressure_sign"])["mean_opening_Fy_N"]
        .sum()
        .to_string()
    )


if __name__ == "__main__":
    main()
