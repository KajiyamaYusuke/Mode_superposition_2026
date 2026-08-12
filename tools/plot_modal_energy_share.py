#!/usr/bin/env python3
"""Plot the time history of modal-energy shares for the leading modes."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot modal-energy shares calculated from q and qdot."
    )
    parser.add_argument(
        "--run-dir", type=Path, default=PROJECT_ROOT / "output",
        help="simulation output directory (default: repository output)",
    )
    parser.add_argument("--start-time", type=float)
    parser.add_argument("--end-time", type=float)
    parser.add_argument("--top-count", type=int, default=5)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def calculate_energy_share(
    modal: pd.DataFrame, start_time: float | None, end_time: float | None
) -> pd.DataFrame:
    required = {"time_s", "side", "mode_index", "frequency_hz", "q", "qdot"}
    missing = required - set(modal.columns)
    if missing:
        raise ValueError(f"modal_contribution.csv is missing columns: {sorted(missing)}")

    selected = modal.copy()
    if start_time is not None:
        selected = selected[selected["time_s"] >= start_time]
    if end_time is not None:
        selected = selected[selected["time_s"] <= end_time]
    if selected.empty:
        raise ValueError("the selected analysis interval contains no modal samples")

    numeric = ["time_s", "mode_index", "frequency_hz", "q", "qdot"]
    if not np.all(np.isfinite(selected[numeric].to_numpy(dtype=float))):
        raise ValueError("modal_contribution.csv contains NaN or Inf")

    omega = 2.0 * np.pi * selected["frequency_hz"].to_numpy(dtype=float)
    q = selected["q"].to_numpy(dtype=float)
    qdot = selected["qdot"].to_numpy(dtype=float)
    selected["kinetic_energy_j"] = 0.5 * qdot * qdot
    selected["potential_energy_j"] = 0.5 * (omega * q) ** 2
    selected["modal_energy_j"] = (
        selected["kinetic_energy_j"] + selected["potential_energy_j"]
    )
    total = selected.groupby(["time_s", "side"])["modal_energy_j"].transform("sum")
    selected["energy_share"] = np.where(
        total > np.finfo(float).tiny, selected["modal_energy_j"] / total, 0.0
    )
    selected["energy_share_percent"] = 100.0 * selected["energy_share"]
    return selected


def summarize_modes(energy: pd.DataFrame) -> pd.DataFrame:
    summary = (
        energy.groupby(["side", "mode_index", "frequency_hz"], as_index=False)
        .agg(
            mean_energy_j=("modal_energy_j", "mean"),
            mean_energy_share_percent=("energy_share_percent", "mean"),
            max_energy_share_percent=("energy_share_percent", "max"),
        )
    )
    summary["mean_energy_rank"] = (
        summary.groupby("side")["mean_energy_j"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    return summary.sort_values(["side", "mean_energy_rank"])


def plot_energy_share(
    energy: pd.DataFrame, summary: pd.DataFrame, top_count: int, output: Path
) -> None:
    sides = [side for side in ("L", "R") if side in set(energy["side"])]
    if not sides:
        raise ValueError("no L or R modal records were found")

    fig, axes = plt.subplots(
        len(sides), 1, figsize=(12, 4.2 * len(sides)), sharex=True, squeeze=False
    )
    for axis, side in zip(axes[:, 0], sides):
        leaders = summary[
            (summary["side"] == side) & (summary["mean_energy_rank"] <= top_count)
        ]
        side_data = energy[energy["side"] == side]
        for row in leaders.itertuples(index=False):
            mode_data = side_data[side_data["mode_index"] == row.mode_index]
            axis.plot(
                mode_data["time_s"], mode_data["energy_share_percent"],
                linewidth=1.1,
                label=(f"Mode {int(row.mode_index) + 1} "
                       f"({row.frequency_hz:.1f} Hz, mean {row.mean_energy_share_percent:.1f}%)"),
            )
        axis.set_ylabel("Energy share [%]")
        axis.set_title(f"{side} vocal fold: top {len(leaders)} modes by mean energy")
        axis.set_ylim(0.0, 100.0)
        axis.grid(alpha=0.3)
        axis.legend(loc="upper right", fontsize=8, ncol=2)
    axes[-1, 0].set_xlabel("Time [s]")
    fig.suptitle("Modal energy share relative to all retained modes")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)


def main() -> None:
    args = arguments()
    if args.top_count < 1:
        raise ValueError("--top-count must be at least 1")
    run_dir = args.run_dir.resolve()
    source = run_dir / "modal_contribution.csv"
    if not source.is_file():
        raise FileNotFoundError(
            f"{source} is missing; run the simulation with modal output enabled"
        )

    energy = calculate_energy_share(
        pd.read_csv(source), args.start_time, args.end_time
    )
    summary = summarize_modes(energy)
    energy.to_csv(run_dir / "modal_energy_share_timeseries.csv", index=False)
    summary.to_csv(run_dir / "modal_energy_share_summary.csv", index=False)
    output = run_dir / "figures" / "modal_energy_share_top5.png"
    plot_energy_share(energy, summary, args.top_count, output)
    print(f"Wrote {output}")
    print(f"Wrote {run_dir / 'modal_energy_share_timeseries.csv'}")
    print(f"Wrote {run_dir / 'modal_energy_share_summary.csv'}")
    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
