import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot a value column from pair offset summary CSV."
    )
    parser.add_argument(
        "--summary",
        default="output_pair_offset/pair_offset_summary.csv",
        help="Path to pair offset summary CSV.",
    )
    parser.add_argument(
        "--value-column",
        default="max_abs_dx_deformation",
        help="CSV column to plot.",
    )
    parser.add_argument(
        "--i-column",
        default="max_def_i",
        help="CSV column for the peak i annotation.",
    )
    parser.add_argument(
        "--j-column",
        default="max_def_j",
        help="CSV column for the peak j annotation.",
    )
    parser.add_argument("--dt", type=float, default=1.0e-5, help="Time step [s].")
    parser.add_argument("--start", type=float, default=0.05, help="Start time [s].")
    parser.add_argument("--end", type=float, default=0.10, help="End time [s].")
    parser.add_argument(
        "--out",
        default="output_pair_offset/max_abs_dx_deformation_0p05_0p10s.png",
        help="Output PNG path.",
    )
    args = parser.parse_args()

    times = []
    values = []
    max_i = []
    max_j = []

    with open(args.summary, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            step = int(row["step"])
            time = step * args.dt
            if args.start <= time <= args.end:
                times.append(time)
                values.append(float(row[args.value_column]))
                max_i.append(int(row[args.i_column]))
                max_j.append(int(row[args.j_column]))

    if not times:
        raise RuntimeError(
            f"No rows found in {args.summary} for {args.start} <= t <= {args.end}."
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    peak_idx = max(range(len(values)), key=values.__getitem__)

    plt.figure(figsize=(8, 4.5))
    plt.plot(times, values, marker="o", markersize=3, linewidth=1.5)
    plt.scatter([times[peak_idx]], [values[peak_idx]], color="crimson", zorder=3)
    plt.annotate(
        f"max={values[peak_idx]:.6g}\n"
        f"t={times[peak_idx]:.5f}s\n"
        f"(i,j)=({max_i[peak_idx]},{max_j[peak_idx]})",
        xy=(times[peak_idx], values[peak_idx]),
        xytext=(10, 15),
        textcoords="offset points",
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "crimson"},
    )
    plt.xlabel("time [s]")
    plt.ylabel(f"{args.value_column} [mm]")
    #plt.title("Correspondence offset in flow direction")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)

    print(f"Saved: {out_path}")
    print(
        "Peak: "
        f"t={times[peak_idx]:.8f} s, "
        f"{args.value_column}={values[peak_idx]:.12g}, "
        f"(i,j)=({max_i[peak_idx]},{max_j[peak_idx]})"
    )


if __name__ == "__main__":
    main()
