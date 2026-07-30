#!/usr/bin/env python3
"""Plot the minimum glottal area over the final part of a simulation."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401  (registers the SciencePlots styles)


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "output" / "area.dat"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the minimum area, showing the final 0.05 s by default."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input area.dat path (default: output/area.dat).",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.05,
        help="Displayed duration [s] (default: 0.05).",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=None,
        help="Start time [s]. If omitted, show the final --duration seconds.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0:
        raise ValueError("--duration must be greater than zero.")

    harea = np.loadtxt(args.input, ndmin=2)
    if harea.shape[1] < 2:
        raise ValueError(f"{args.input} must contain at least two columns.")

    # As in show.py, column 1 is the step number and the remaining columns
    # contain area values. One simulation step corresponds to 1e-5 s.
    time = harea[:, 0] * 1.0e-5
    row_min = np.min(harea[:, 1:], axis=1)

    data_end_time = np.max(time)
    if args.start is None:
        start_time = data_end_time - args.duration
        end_time = data_end_time
    else:
        start_time = args.start
        end_time = start_time + args.duration

    mask = (time >= start_time) & (time <= end_time)
    if not np.any(mask):
        raise ValueError(
            f"No data in the requested interval {start_time:g}--{end_time:g} s "
            f"(data range: {np.min(time):g}--{data_end_time:g} s)."
        )

    plt.style.use(["science", "ieee", "no-latex"])
    plt.figure(figsize=(12, 4), dpi=100)
    plt.plot(time[mask], row_min[mask], linestyle="-", color="#0072B2")
    plt.xlabel("Time [s]", fontsize=24)
    plt.ylabel("area [mm^2]", fontsize=24)
    plt.tick_params(labelsize=20)
    plt.ylim(-0.5, 15)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
