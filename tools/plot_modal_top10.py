#!/usr/bin/env python3
"""Plot the three largest time-averaged modal displacement contributions."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots  # noqa: F401  (registers the SciencePlots styles)
from matplotlib.ticker import PercentFormatter


REPO_ROOT = Path(__file__).resolve().parents[1]
SCOPE_TITLES = {
    "surface_uy": "opening-direction",
    "surface_xyz": "total-displacement",
}
SIDES = ("L", "R")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the three largest mean modal displacement contributions."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=REPO_ROOT / "output" / "modal_top10.csv",
        help="Input CSV (default: output/modal_top10.csv).",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=None,
        help="Averaging interval start [s].",
    )
    parser.add_argument(
        "--end",
        type=float,
        default=None,
        help="Averaging interval end [s].",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.05,
        help="Interval duration [s] when start or end is omitted (default: 0.05).",
    )
    parser.add_argument(
        "--scope",
        choices=tuple(SCOPE_TITLES),
        default="surface_uy",
        help="Contribution scope (default: surface_uy).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "output",
        help="Directory for PNG files (default: output).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show figures interactively after saving.",
    )
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_csv(path)
    required = {
        "step",
        "time_s",
        "side",
        "scope",
        "rank",
        "mode_index",
        "frequency_hz",
        "magnitude_ratio",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing CSV columns: {', '.join(sorted(missing))}")
    for column in required.difference({"side", "scope"}):
        data[column] = pd.to_numeric(data[column], errors="coerce")
    if data.empty:
        raise ValueError(f"No rows found in {path}")
    return data


def averaging_interval(
    data: pd.DataFrame,
    start: float | None,
    end: float | None,
    duration: float,
) -> tuple[float, float]:
    if duration <= 0.0:
        raise ValueError("--duration must be greater than zero.")
    data_start = float(data["time_s"].min())
    data_end = float(data["time_s"].max())

    if start is None and end is None:
        end = data_end
        start = end - duration
    elif start is None:
        start = end - duration
    elif end is None:
        end = start + duration

    assert start is not None and end is not None
    if end < start:
        raise ValueError("--end must be greater than or equal to --start.")
    if end < data_start or start > data_end:
        raise ValueError(
            f"Requested interval {start:g}--{end:g} s does not overlap "
            f"data range {data_start:g}--{data_end:g} s."
        )
    return start, end


def average_contribution(
    interval_data: pd.DataFrame,
    side: str,
    scope: str,
) -> tuple[pd.DataFrame, int]:
    rows = interval_data[
        (interval_data["side"] == side)
        & (interval_data["scope"] == scope)
    ].copy()
    times = np.sort(rows["time_s"].dropna().unique())
    if len(times) == 0:
        raise ValueError(f"No {scope} data for side={side}")

    # modal_top10.csv only stores the instantaneous top 10. A mode absent at
    # one output time is therefore treated as zero contribution at that time.
    frequency = rows.groupby("mode_index")["frequency_hz"].first()
    magnitude = rows.groupby("mode_index")["magnitude_ratio"].sum() / len(times)

    average = pd.DataFrame(
        {
            "frequency_hz": frequency,
            "magnitude_ratio": magnitude,
        }
    )
    average.index.name = "mode_index"
    average = average.sort_values(
        ["magnitude_ratio", "frequency_hz"],
        ascending=[False, True],
    ).head(3)
    return average, len(times)


def labels(rows: pd.DataFrame) -> list[str]:
    return [
        f"Mode {int(mode)+1}  ({frequency:.1f} Hz)"
        for mode, frequency in zip(rows.index, rows["frequency_hz"])
    ]


def plot_side(
    interval_data: pd.DataFrame,
    side: str,
    scope: str,
    start: float,
    end: float,
    output: Path,
) -> plt.Figure:
    average, sample_count = average_contribution(
        interval_data, side, scope
    )
    plot_rows = average.sort_values("magnitude_ratio", ascending=True)
    magnitude = plot_rows["magnitude_ratio"].to_numpy(dtype=float)
    y = np.arange(len(plot_rows))

    fig, axis = plt.subplots(figsize=(7, 3.6), constrained_layout=True)
    bars = axis.barh(y, magnitude, color="#0072B2")
    axis.set_yticks(y, labels(plot_rows), fontsize=20)
    axis.xaxis.set_major_formatter(PercentFormatter(1.0))
    axis.set_xlabel("Time-averaged modal contribution ratio [%]", fontsize="20")
    # axis.set_title(
    #     f"Side {side}: {SCOPE_TITLES[scope]} modal contribution\n"
    #     f"t = {start:.6f}--{end:.6f} s, samples = {sample_count}"
    # )
    # axis.bar_label(bars, labels=[f"{value:.1%}" for value in magnitude], padding=4)
    axis.set_xlim(0.0, max(magnitude) * 1.18 if len(magnitude) else 1.0)
    axis.grid(axis="x", alpha=0.3)
    axis.tick_params(labelsize=16)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    return fig


def main() -> None:
    args = parse_args()
    data = load_data(args.input)
    start, end = averaging_interval(data, args.start, args.end, args.duration)
    interval_data = data[
        (data["time_s"] >= start) & (data["time_s"] <= end)
    ]
    if interval_data.empty:
        raise ValueError(f"No output samples in interval {start:g}--{end:g} s.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.style.use(["science", "ieee", "no-latex"])

    figures = []
    for side in SIDES:
        output = args.output_dir / f"modal_top3_average_{args.scope}_{side}.png"
        figures.append(
            plot_side(interval_data, side, args.scope, start, end, output)
        )
        print(f"Saved: {output}")
    print(f"Averaging interval: {start:.12g}--{end:.12g} s")

    if args.show:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


if __name__ == "__main__":
    main()
