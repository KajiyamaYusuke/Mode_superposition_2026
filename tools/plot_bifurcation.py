#!/usr/bin/env python3
"""Create pressure bifurcation and summary plots from a sweep directory."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


POINT_QUANTITIES = {
    "peak_area_mm2",
    "trough_area_mm2",
    "period_s",
    "amplitude_mm2",
}
SUMMARY_LABELS = {
    "period_cv": "Period CV",
    "amplitude_cv": "Amplitude CV",
    "f0_from_cycles_hz": "Fundamental frequency [Hz]",
    "spectral_entropy": "Spectral entropy",
    "contact_fraction": "Contact fraction",
    "period_ratio_unsigned": "Unsigned left/right period ratio [-]",
    "amplitude_ratio_unsigned": "Unsigned left/right amplitude ratio [-]",
    "period_ratio_left_over_right": "Left/right period ratio [-]",
    "amplitude_ratio_left_over_right": "Left/right amplitude ratio [-]",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a pressure-sweep bifurcation diagram.")
    parser.add_argument("--sweep-dir", type=Path, required=True)
    parser.add_argument(
        "--quantity",
        default="peak_area_mm2",
        choices=sorted(POINT_QUANTITIES | set(SUMMARY_LABELS)),
    )
    parser.add_argument("--color-by")
    parser.add_argument("--pressure-min", type=float)
    parser.add_argument("--pressure-max", type=float)
    parser.add_argument("--save", type=Path)
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    sweep_dir = args.sweep_dir.resolve()
    point_file = sweep_dir / "bifurcation_points.csv"
    summary_file = sweep_dir / "sweep_summary.csv"
    source_file = point_file if args.quantity in POINT_QUANTITIES else summary_file
    if not source_file.is_file():
        raise FileNotFoundError(source_file)
    data = pd.read_csv(source_file)
    data["pressure_pa"] = pd.to_numeric(data["pressure_pa"], errors="coerce")
    data[args.quantity] = pd.to_numeric(data[args.quantity], errors="coerce")
    data = data.dropna(subset=["pressure_pa", args.quantity])
    if args.pressure_min is not None:
        data = data[data["pressure_pa"] >= args.pressure_min]
    if args.pressure_max is not None:
        data = data[data["pressure_pa"] <= args.pressure_max]
    if data.empty:
        raise ValueError("no valid points remain for the requested plot")

    fig, axis = plt.subplots(figsize=(8, 5))
    color = None
    if args.color_by:
        if args.color_by not in data:
            raise ValueError(f"--color-by column does not exist: {args.color_by}")
        color = pd.to_numeric(data[args.color_by], errors="coerce")
    scatter = axis.scatter(
        data["pressure_pa"],
        data[args.quantity],
        c=color,
        s=10 if args.quantity in POINT_QUANTITIES else 28,
        alpha=0.75,
    )
    if color is not None:
        fig.colorbar(scatter, ax=axis, label=args.color_by)
    axis.set_xlabel("Subglottal pressure Ps [Pa]")
    axis.set_ylabel(SUMMARY_LABELS.get(args.quantity, args.quantity))
    axis.set_title(f"Pressure Bifurcation: {args.quantity}")
    axis.grid(alpha=0.25)
    fig.tight_layout()
    output = args.save or sweep_dir / "figures" / f"bifurcation_{args.quantity}.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=250, bbox_inches="tight")
    print(f"Plotted {len(data)} points; saved: {output}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
