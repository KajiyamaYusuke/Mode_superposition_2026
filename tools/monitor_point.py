import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot monitor_point.csv.")
    parser.add_argument(
        "--input",
        default="../output_pair_offset/monitor_point.csv",
        help="Path to monitor_point.csv.",
    )
    parser.add_argument("--start", type=float, default=None, help="Start time after scaling.")
    parser.add_argument("--end", type=float, default=None, help="End time after scaling.")
    parser.add_argument(
        "--time-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to the CSV time column.",
    )
    parser.add_argument(
        "--time-offset",
        type=float,
        default=0.0,
        help="Offset added after scaling the CSV time column.",
    )
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Save PNGs with this prefix instead of showing interactive windows.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df["plot_time"] = df["time"] * args.time_scale + args.time_offset

    if args.start is not None:
        df = df[df["plot_time"] >= args.start]
    if args.end is not None:
        df = df[df["plot_time"] <= args.end]
    if df.empty:
        raise RuntimeError("Selected time range is empty.")

    def plot_pair(columns, ylabel, suffix):
        plt.figure()
        for column in columns:
            plt.plot(df["plot_time"], df[column], label=column)
        plt.legend()
        plt.xlabel("time [s]")
        plt.ylabel(ylabel)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if args.out_prefix:
            out_path = Path(f"{args.out_prefix}_{suffix}.png")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(out_path, dpi=200)
            print(f"Saved: {out_path}")

    plot_pair(["xDispL", "xDispR"], "x displacement [mm]", "x")
    plot_pair(["zDispL", "zDispR"], "z displacement [mm]", "z")

    if not args.out_prefix:
        plt.show()


if __name__ == "__main__":
    main()
