#!/usr/bin/env python3
"""Plot arbitrary variables from generated Poincaré crossing points."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a Poincaré section from one run.")
    parser.add_argument("--run-dir", type=Path, default=Path("output"))
    parser.add_argument("--section-signal", default="min_area", help="metadata label")
    parser.add_argument("--section-level", type=float, default=0.5, help="metadata label")
    parser.add_argument("--direction", choices=["rising", "falling"], default="rising")
    parser.add_argument("--x", dest="x_column")
    parser.add_argument("--y", dest="y_column")
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--last-points", type=int)
    parser.add_argument("--connect", action="store_true")
    parser.add_argument("--annotate-order", action="store_true")
    parser.add_argument("--list-columns", action="store_true")
    parser.add_argument("--save", type=Path)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    path = args.run_dir / "poincare_points.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} does not exist; run tools/analyze_dynamics.py first"
        )
    points = pd.read_csv(path)
    print("Available columns:", ", ".join(points.columns))
    if args.list_columns:
        return
    if args.start_time is not None:
        points = points[points["time_s"] >= args.start_time]
    if args.end_time is not None:
        points = points[points["time_s"] <= args.end_time]
    if args.last_points:
        points = points.tail(args.last_points)
    if points.empty:
        raise ValueError("no Poincaré points remain after filtering")

    available = set(points.columns)
    x_column = args.x_column or ("qL_1" if "qL_1" in available else "uyL_mm")
    y_column = args.y_column or ("qdotL_1" if "qdotL_1" in available else "uyR_mm")
    missing = [name for name in (x_column, y_column) if name not in available]
    if missing:
        raise ValueError(f"unknown columns {missing}; available: {sorted(available)}")

    fig, axis = plt.subplots(figsize=(6, 6))
    axis.scatter(points[x_column], points[y_column], c=points["crossing_index"], s=28)
    if args.connect:
        axis.plot(points[x_column], points[y_column], linewidth=0.7, alpha=0.5)
    if args.annotate_order:
        for _, row in points.iterrows():
            axis.annotate(
                str(int(row["crossing_index"])),
                (row[x_column], row[y_column]),
                fontsize=7,
            )
    axis.set(
        xlabel=x_column,
        ylabel=y_column,
        title=(
            f"Poincaré section: {args.section_signal}, "
            f"{args.direction} at {args.section_level:g}"
        ),
    )
    axis.grid(alpha=0.3)
    fig.tight_layout()
    output = args.save or args.run_dir / "figures" / f"poincare_{x_column}_{y_column}.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=250, bbox_inches="tight")
    print(f"Saved Poincaré plot: {output}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
