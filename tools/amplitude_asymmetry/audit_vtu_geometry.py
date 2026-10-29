#!/usr/bin/env python3
"""Audit every saved deformation VTU for tetrahedral sign changes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd


def input_points(path: Path) -> np.ndarray:
    inside_points = False
    inside_array = False
    values: list[tuple[float, float, float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if "<Points>" in line:
                inside_points = True
                continue
            if not inside_points:
                continue
            if "<DataArray" in line:
                inside_array = True
                continue
            if "</DataArray>" in line:
                break
            if inside_array:
                tokens = line.split()
                for start in range(0, len(tokens), 3):
                    if start + 2 < len(tokens):
                        raw_x, raw_y, raw_z = map(float, tokens[start : start + 3])
                        values.append((raw_z, raw_x, raw_y))
    result = np.asarray(values)
    if result.ndim != 2 or result.shape[1] != 3:
        raise ValueError(f"failed to read input points from {path}")
    return result


def output_points(path: Path) -> np.ndarray:
    inside_points = False
    inside_array = False
    values: list[tuple[float, float, float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if "<Points>" in line:
                inside_points = True
                continue
            if not inside_points:
                continue
            if "<DataArray" in line:
                inside_array = True
                continue
            if "</DataArray>" in line:
                break
            if inside_array:
                tokens = line.split()
                for start in range(0, len(tokens), 3):
                    if start + 2 < len(tokens):
                        values.append(tuple(map(float, tokens[start : start + 3])))
    return np.asarray(values)


def output_cells(path: Path) -> tuple[np.ndarray, np.ndarray]:
    tree = ET.parse(path)
    connectivity = np.fromstring(
        tree.find('.//DataArray[@Name="connectivity"]').text or "", sep=" ", dtype=np.int64
    )
    offsets = np.fromstring(
        tree.find('.//DataArray[@Name="offsets"]').text or "", sep=" ", dtype=np.int64
    )
    types = np.fromstring(
        tree.find('.//DataArray[@Name="types"]').text or "", sep=" ", dtype=np.int64
    )
    starts = np.r_[0, offsets[:-1]]
    tetra_cells = np.flatnonzero(types == 10)
    if not np.all(offsets[tetra_cells] - starts[tetra_cells] == 4):
        raise ValueError("VTK tetra connectivity is not four-node")
    indices = np.stack(
        [connectivity[starts[tetra_cells] + offset] for offset in range(4)], axis=1
    )
    return indices, types


def signed_six_volume(points: np.ndarray, tetra: np.ndarray) -> np.ndarray:
    origin = points[tetra[:, 0]]
    a = points[tetra[:, 1]] - origin
    b = points[tetra[:, 2]] - origin
    c = points[tetra[:, 3]] - origin
    return np.einsum("ij,ij->i", np.cross(a, b), c)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-mode", type=Path, required=True)
    parser.add_argument("--vtu-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checks", type=Path, required=True)
    arguments = parser.parse_args()
    files = sorted(arguments.vtu_dir.glob("deform_combined*.vtu"))
    if not files:
        raise FileNotFoundError("no deformation VTUs")
    initial_left = input_points(arguments.input_mode)
    initial_right = initial_left.copy()
    initial_right[:, 1] = 2.0 * 8.4 - initial_right[:, 1]
    initial = np.vstack([initial_left, initial_right])
    tetra, types = output_cells(files[0])
    reference = signed_six_volume(initial, tetra)
    valid_reference = np.abs(reference) > 1.0e-15
    if not np.all(valid_reference):
        raise ValueError("initial mesh contains zero-volume tetrahedra")
    rows = []
    for index, path in enumerate(files):
        points = output_points(path)
        if points.shape != initial.shape:
            raise ValueError(f"{path}: point count mismatch")
        current = signed_six_volume(points, tetra)
        ratio = current / reference
        displacement = np.linalg.norm(points - initial, axis=1)
        rows.append(
            {
                "frame": index,
                "time_s": 0.45 + index * 20.0e-5,
                "file": str(path),
                "tetra_count": len(tetra),
                "sign_change_count": int(np.sum(ratio <= 0.0)),
                "minimum_signed_volume_ratio": float(np.min(ratio)),
                "p01_signed_volume_ratio": float(np.percentile(ratio, 1)),
                "maximum_nodal_displacement_mm": float(np.max(displacement)),
                "nonfinite_points": int(np.sum(~np.isfinite(points))),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(arguments.output, index=False)
    checks = json.loads(arguments.checks.read_text(encoding="utf-8"))
    checks["geometry_vtu"] = {
        "passed": bool(
            frame.sign_change_count.eq(0).all()
            and frame.nonfinite_points.eq(0).all()
        ),
        "saved_time_scope_s": [float(frame.time_s.min()), float(frame.time_s.max())],
        "saved_frame_count": len(frame),
        "tetra_count_bilateral": int(frame.tetra_count.iloc[0]),
        "maximum_sign_change_count": int(frame.sign_change_count.max()),
        "minimum_signed_volume_ratio": float(frame.minimum_signed_volume_ratio.min()),
        "maximum_nodal_displacement_mm": float(frame.maximum_nodal_displacement_mm.max()),
        "limitation": "VTU audit covers saved 0.45-0.4998 s only; online channel AP-order checks covered all solver steps.",
    }
    arguments.checks.write_text(
        json.dumps(checks, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(checks["geometry_vtu"], indent=2))


if __name__ == "__main__":
    main()
