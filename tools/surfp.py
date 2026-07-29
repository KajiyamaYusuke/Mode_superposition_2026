#!/usr/bin/env python3
"""Visualize the surface map and the displacement/contact monitor point."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
MONITOR_TARGET_X = 9.6
MONITOR_TARGET_Z = 8.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot surfp nodes and highlight the nearest monitor point."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "output" / "surfp_output.csv",
        help="surface-point CSV (default: output/surfp_output.csv)",
    )
    parser.add_argument("--target-x", type=float, default=MONITOR_TARGET_X)
    parser.add_argument("--target-z", type=float, default=MONITOR_TARGET_Z)
    parser.add_argument(
        "--nxsup",
        type=int,
        default=46,
        help="number of i rows searched by Simulation (default: 46)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="save the figure instead of only displaying it",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="do not open an interactive plot window",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.nxsup <= 0:
        raise ValueError("--nxsup must be positive")
    if not args.input.is_file():
        raise FileNotFoundError(f"surface-point file does not exist: {args.input}")

    columns = ["i", "j", "node_id", "x", "y", "z"]
    points = pd.read_csv(args.input, comment="#", header=None, names=columns)
    if points.empty:
        raise ValueError(f"surface-point file contains no data: {args.input}")

    # Match Simulation::run(): search only i < nxsup and minimize distance
    # in the x-z plane. Invalid surfp entries are absent from this CSV.
    candidates = points[points["i"] < args.nxsup].copy()
    if candidates.empty:
        raise ValueError("no monitor candidates remain after applying --nxsup")
    candidates["distance2"] = (
        (candidates["x"] - args.target_x) ** 2
        + (candidates["z"] - args.target_z) ** 2
    )
    monitor = candidates.loc[candidates["distance2"].idxmin()]

    fig = plt.figure(figsize=(9, 7))
    axis = fig.add_subplot(111, projection="3d")
    scatter = axis.scatter(
        points["x"],
        points["y"],
        points["z"],
        c=points["j"],
        cmap="viridis",
        s=16,
        alpha=0.75,
    )
    fig.colorbar(scatter, ax=axis, label="j index", shrink=0.75)

    axis.scatter(
        monitor["x"],
        monitor["y"],
        monitor["z"],
        color="red",
        marker="*",
        s=220,
        edgecolor="black",
        linewidth=0.8,
        label="Selected monitor node",
        depthshade=False,
    )
    axis.scatter(
        args.target_x,
        monitor["y"],
        args.target_z,
        color="black",
        marker="x",
        s=90,
        linewidth=2,
        label="Monitor target (x-z)",
        depthshade=False,
    )
    axis.plot(
        [args.target_x, monitor["x"]],
        [monitor["y"], monitor["y"]],
        [args.target_z, monitor["z"]],
        color="red",
        linestyle="--",
        linewidth=1,
    )
    axis.text(
        monitor["x"],
        monitor["y"],
        monitor["z"],
        f"  (i,j)=({int(monitor['i'])},{int(monitor['j'])})"
        f"\n  node={int(monitor['node_id'])}",
        color="red",
    )

    axis.set_xlabel("X [mm]")
    axis.set_ylabel("Y [mm]")
    axis.set_zlabel("Z [mm]")
    axis.set_title("Surface Nodes and Monitor Point")
    axis.legend(loc="upper left")
    fig.tight_layout()

    print(
        "Monitor point: "
        f"target=({args.target_x:g}, {args.target_z:g}), "
        f"(i,j)=({int(monitor['i'])}, {int(monitor['j'])}), "
        f"node={int(monitor['node_id'])}, "
        f"xyz=({monitor['x']:g}, {monitor['y']:g}, {monitor['z']:g})"
    )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=300, bbox_inches="tight")
        print(f"Saved surface plot: {args.output}")
    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
