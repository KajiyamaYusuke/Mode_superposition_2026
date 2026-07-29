#!/usr/bin/env python3
"""Plot left/right x-y displacement histories from displace_xy.dat."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize x and y vocal-fold displacements."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "output" / "displace_xy.dat",
        help="input data file (default: output/displace_xy.dat)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=PROJECT_ROOT / "output" / "displacement_xy.png",
        help="image output path (default: output/displacement_xy.png)",
    )
    parser.add_argument("--start", type=float, help="start time in seconds")
    parser.add_argument("--end", type=float, help="end time in seconds")
    parser.add_argument(
        "--show", action="store_true", help="also display the graph interactively"
    )
    return parser.parse_args()


def load_displacement(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"input file does not exist: {path}")

    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != 5:
        raise ValueError(
            f"expected 5 columns (time, uxL, uyL, uxR, uyR), got {data.shape[1]}"
        )
    if not np.all(np.isfinite(data)):
        raise ValueError(f"input contains non-finite values: {path}")
    return data


def main() -> None:
    args = parse_args()
    if args.start is not None and args.end is not None and args.start >= args.end:
        raise ValueError("--start must be smaller than --end")

    data = load_displacement(args.input)
    time = data[:, 0]
    mask = np.ones(time.shape, dtype=bool)
    if args.start is not None:
        mask &= time >= args.start
    if args.end is not None:
        mask &= time <= args.end
    if not np.any(mask):
        raise ValueError("the selected time range contains no samples")

    time = time[mask]
    ux_left, uy_left, ux_right, uy_right = (data[mask, i] for i in range(1, 5))

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    series = (
        (axes[0], ux_left, ux_right, "x displacement [mm]"),
        (axes[1], uy_left, uy_right, "y displacement [mm]"),
    )
    for axis, left, right, ylabel in series:
        axis.plot(time, left, label="Left vocal fold", linewidth=1.2)
        axis.plot(time, right, label="Right vocal fold", linewidth=1.2)
        axis.set_ylabel(ylabel)
        axis.grid(True, linestyle="--", alpha=0.4)
        axis.legend()

    axes[1].set_xlabel("Time [s]")
    fig.suptitle("Vocal Fold Displacement")
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    #fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Saved displacement plot: {args.output}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
