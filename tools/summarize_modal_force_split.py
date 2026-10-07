#!/usr/bin/env python3
"""Build compact modal-force summaries without loading the source CSV at once."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


KEY_COLUMNS = ["side", "mode_index_1based"]
CONTACT_COLUMNS = [
    "contact_force_ap_low",
    "contact_force_ap_high",
    "contact_force_interior",
]
EQUILIBRIUM_COLUMNS = [
    "q",
    "fluid_force",
    *CONTACT_COLUMNS,
    "contact_force_total",
    "total_force",
    "applied_force",
]
LOAD_COLUMNS = [column for column in EQUILIBRIUM_COLUMNS if column != "q"]
CHECK_COLUMNS = [
    "balance_residual",
    "load_decomposition_error",
    "computed_load_decomposition_error",
    "applied_minus_total_force",
]
REQUIRED_COLUMNS = {
    "step",
    "load_time_s",
    "state_time_s",
    "side",
    "mode_index_1based",
    "frequency_hz",
    "q",
    "fluid_force",
    *CONTACT_COLUMNS,
    "total_force",
    "applied_force",
    "equation_residual",
    "computed_contact_ap_low",
    "computed_contact_ap_high",
    "computed_contact_interior",
    "diagnostic_load_mode",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--selected", type=Path)
    parser.add_argument("--equilibrium-window", nargs=2, type=float, default=(0.25, 0.30))
    parser.add_argument("--post-window", nargs=2, type=float, default=(0.37, 0.60))
    parser.add_argument("--top-count", type=int, default=3)
    parser.add_argument("--required-mode", type=int, default=5)
    parser.add_argument("--output-interval-s", type=float, default=2.0e-4)
    parser.add_argument("--chunksize", type=int, default=200_000)
    return parser.parse_args()


def add_derived_columns(frame: pd.DataFrame) -> None:
    frame["contact_force_total"] = frame[CONTACT_COLUMNS].sum(axis=1)
    frame["balance_residual"] = frame["equation_residual"]
    frame["load_decomposition_error"] = (
        frame["total_force"] - frame["fluid_force"] - frame["contact_force_total"]
    )
    frame["computed_load_decomposition_error"] = (
        frame["total_force"]
        - frame["fluid_force"]
        - frame["computed_contact_ap_low"]
        - frame["computed_contact_ap_high"]
        - frame["computed_contact_interior"]
    )
    frame["applied_minus_total_force"] = frame["applied_force"] - frame["total_force"]


def accumulate(
    accumulator: pd.DataFrame | None,
    update: pd.DataFrame,
) -> pd.DataFrame:
    if accumulator is None:
        return update
    return accumulator.add(update, fill_value=0.0)


def grouped_sum(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame.groupby(KEY_COLUMNS, observed=True, sort=False)[columns].sum()


def grouped_count(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame.groupby(KEY_COLUMNS, observed=True, sort=False)[columns].count()


def grouped_max_abs(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    values = frame[columns].abs().copy()
    values[KEY_COLUMNS] = frame[KEY_COLUMNS]
    return values.groupby(KEY_COLUMNS, observed=True, sort=False)[columns].max()


def prefixed_moments(
    output: pd.DataFrame,
    sums: pd.DataFrame,
    squares: pd.DataFrame,
    counts: pd.DataFrame,
    maxima: pd.DataFrame,
    prefix: str,
) -> None:
    means = sums / counts
    rms = np.sqrt((squares / counts).clip(lower=0.0))
    for column in sums.columns:
        output[f"{prefix}_{column}_mean"] = means[column]
        output[f"{prefix}_{column}_rms"] = rms[column]
        output[f"{prefix}_{column}_max_abs"] = maxima[column]


def validate_window(name: str, values: tuple[float, float] | list[float]) -> None:
    if len(values) != 2 or values[0] > values[1]:
        raise ValueError(f"invalid {name}: {values}")


def closed_window_mask(values: pd.Series, start: float, end: float) -> pd.Series:
    """Include decimal CSV endpoints despite their nearest binary representation."""
    lower = np.nextafter(start, -np.inf)
    upper = np.nextafter(end, np.inf)
    return values.between(lower, upper, inclusive="both")


def main() -> None:
    args = arguments()
    validate_window("equilibrium window", args.equilibrium_window)
    validate_window("post window", args.post_window)
    if args.top_count < 1 or args.output_interval_s <= 0.0 or args.chunksize < 1:
        raise ValueError("top-count, output-interval-s, and chunksize must be positive")
    if not args.input_csv.is_file():
        raise FileNotFoundError(args.input_csv)

    output_directory = args.input_csv.parent
    summary_path = args.summary or output_directory / "modal_force_summary.csv"
    selected_path = args.selected or output_directory / "modal_force_selected.csv"

    header = pd.read_csv(args.input_csv, nrows=0)
    missing = sorted(REQUIRED_COLUMNS - set(header.columns))
    if missing:
        raise ValueError(f"missing required columns: {', '.join(missing)}")
    source_columns = list(header.columns)

    equilibrium_sums = equilibrium_counts = None
    check_sums = check_squares = check_counts = check_maxima = None
    all_check_sums = all_check_squares = all_check_counts = all_check_maxima = None
    frequencies: dict[tuple[str, int], float] = {}
    load_modes: dict[tuple[str, int], set[str]] = {}
    timing_rows: list[pd.DataFrame] = []

    eq_start, eq_end = args.equilibrium_window
    for chunk_number, chunk in enumerate(pd.read_csv(args.input_csv, chunksize=args.chunksize)):
        add_derived_columns(chunk)
        if chunk_number == 0:
            timing_rows.append(chunk[["step", "load_time_s"]].drop_duplicates().head(1000))

        metadata = chunk.groupby(KEY_COLUMNS, observed=True, sort=False)
        for key, value in metadata["frequency_hz"].first().items():
            frequencies[(str(key[0]), int(key[1]))] = float(value)
        for key, values in metadata["diagnostic_load_mode"]:
            normalized = (str(key[0]), int(key[1]))
            load_modes.setdefault(normalized, set()).update(values.astype(str).unique())

        all_check_sums = accumulate(all_check_sums, grouped_sum(chunk, CHECK_COLUMNS))
        all_check_squares = accumulate(
            all_check_squares, grouped_sum(chunk.assign(**{
                column: chunk[column] * chunk[column] for column in CHECK_COLUMNS
            }), CHECK_COLUMNS)
        )
        all_check_counts = accumulate(all_check_counts, grouped_count(chunk, CHECK_COLUMNS))
        maxima = grouped_max_abs(chunk, CHECK_COLUMNS)
        all_check_maxima = (
            maxima if all_check_maxima is None
            else all_check_maxima.combine(maxima, np.fmax)
        )

        equilibrium = chunk[closed_window_mask(chunk["load_time_s"], eq_start, eq_end)]
        if equilibrium.empty:
            continue
        equilibrium_sums = accumulate(
            equilibrium_sums, grouped_sum(equilibrium, EQUILIBRIUM_COLUMNS)
        )
        equilibrium_counts = accumulate(
            equilibrium_counts, grouped_count(equilibrium, EQUILIBRIUM_COLUMNS)
        )
        check_sums = accumulate(check_sums, grouped_sum(equilibrium, CHECK_COLUMNS))
        check_squares = accumulate(
            check_squares, grouped_sum(equilibrium.assign(**{
                column: equilibrium[column] * equilibrium[column] for column in CHECK_COLUMNS
            }), CHECK_COLUMNS)
        )
        check_counts = accumulate(check_counts, grouped_count(equilibrium, CHECK_COLUMNS))
        maxima = grouped_max_abs(equilibrium, CHECK_COLUMNS)
        check_maxima = maxima if check_maxima is None else check_maxima.combine(maxima, np.fmax)

    if equilibrium_sums is None or equilibrium_counts is None:
        raise ValueError("the equilibrium window contains no samples")
    assert check_sums is not None and check_squares is not None
    assert check_counts is not None and check_maxima is not None
    assert all_check_sums is not None and all_check_squares is not None
    assert all_check_counts is not None and all_check_maxima is not None

    equilibrium_means = equilibrium_sums / equilibrium_counts
    post_delta_squares = post_counts = None
    post_start, post_end = args.post_window
    for chunk in pd.read_csv(args.input_csv, chunksize=args.chunksize):
        add_derived_columns(chunk)
        post = chunk[closed_window_mask(chunk["load_time_s"], post_start, post_end)].copy()
        if post.empty:
            continue
        index = pd.MultiIndex.from_frame(post[KEY_COLUMNS])
        delta_squares = post[KEY_COLUMNS].copy()
        for column in LOAD_COLUMNS:
            baseline = equilibrium_means[column].reindex(index).to_numpy()
            delta_squares[column] = np.square(post[column].to_numpy() - baseline)
        post_delta_squares = accumulate(
            post_delta_squares, grouped_sum(delta_squares, LOAD_COLUMNS)
        )
        post_counts = accumulate(post_counts, grouped_count(delta_squares, LOAD_COLUMNS))

    if post_delta_squares is None or post_counts is None:
        raise ValueError("the post-perturbation window contains no samples")

    index = equilibrium_means.index.sort_values()
    summary = pd.DataFrame(index=index)
    summary["frequency_hz"] = [frequencies[(str(side), int(mode))] for side, mode in index]
    summary["diagnostic_load_mode"] = [
        ";".join(sorted(load_modes[(str(side), int(mode))])) for side, mode in index
    ]
    summary["equilibrium_window_start_s"] = eq_start
    summary["equilibrium_window_end_s"] = eq_end
    summary["post_window_start_s"] = post_start
    summary["post_window_end_s"] = post_end
    summary["equilibrium_sample_count"] = equilibrium_counts["q"].reindex(index).astype(int)
    summary["post_sample_count"] = post_counts[LOAD_COLUMNS[0]].reindex(index).astype(int)
    for column in EQUILIBRIUM_COLUMNS:
        summary[f"equilibrium_{column}_mean"] = equilibrium_means[column].reindex(index)
    for column in LOAD_COLUMNS:
        summary[f"post_{column}_delta_rms"] = np.sqrt(
            post_delta_squares[column].reindex(index) / post_counts[column].reindex(index)
        )
    prefixed_moments(
        summary,
        check_sums.reindex(index),
        check_squares.reindex(index),
        check_counts.reindex(index),
        check_maxima.reindex(index),
        "equilibrium",
    )
    prefixed_moments(
        summary,
        all_check_sums.reindex(index),
        all_check_squares.reindex(index),
        all_check_counts.reindex(index),
        all_check_maxima.reindex(index),
        "all_samples",
    )

    summary["equilibrium_q_abs_rank_within_side"] = (
        summary["equilibrium_q_mean"].abs().groupby(level="side").rank(method="first", ascending=False)
    ).astype(int)
    summary["post_contact_delta_rms_rank_within_side"] = (
        summary["post_contact_force_total_delta_rms"]
        .groupby(level="side")
        .rank(method="first", ascending=False)
    ).astype(int)
    summary["selected_for_timeseries"] = False
    summary["selection_reason"] = "not_selected"
    selected_keys: set[tuple[str, int]] = set()
    reasons: dict[tuple[str, int], list[str]] = {}
    for key, row in summary.iterrows():
        normalized = (str(key[0]), int(key[1]))
        key_reasons: list[str] = []
        if int(key[1]) == args.required_mode:
            key_reasons.append(f"required_mode_{args.required_mode}")
        if row["equilibrium_q_abs_rank_within_side"] <= args.top_count:
            key_reasons.append(f"top_{args.top_count}_equilibrium_abs_q")
        if row["post_contact_delta_rms_rank_within_side"] <= args.top_count:
            key_reasons.append(f"top_{args.top_count}_contact_delta_rms")
        if key_reasons:
            selected_keys.add(normalized)
            reasons[normalized] = key_reasons
            summary.loc[key, "selected_for_timeseries"] = True
            summary.loc[key, "selection_reason"] = ";".join(key_reasons)

    timing = pd.concat(timing_rows, ignore_index=True).drop_duplicates().sort_values("step")
    if len(timing) < 2:
        raise ValueError("cannot infer the source output interval")
    step_differences = np.diff(timing["step"].to_numpy(dtype=np.int64))
    time_differences = np.diff(timing["load_time_s"].to_numpy(dtype=float))
    step_increment = int(round(float(np.median(step_differences[step_differences > 0]))))
    source_interval = float(np.median(time_differences[time_differences > 0.0]))
    output_stride = max(1, int(round(args.output_interval_s / source_interval)))
    actual_output_interval = output_stride * source_interval
    first_step = int(timing["step"].iloc[0])
    summary["selected_output_interval_target_s"] = args.output_interval_s
    summary["selected_output_interval_actual_s"] = actual_output_interval
    summary["selected_source_sample_stride"] = output_stride

    summary.reset_index().to_csv(summary_path, index=False, float_format="%.15e")

    metadata_columns = [
        "selection_reason",
        "equilibrium_q_abs_rank_within_side",
        "post_contact_delta_rms_rank_within_side",
    ]
    wrote_header = False
    for chunk in pd.read_csv(args.input_csv, chunksize=args.chunksize):
        normalized_keys = list(zip(chunk["side"].astype(str), chunk["mode_index_1based"].astype(int)))
        selected_mask = np.fromiter((key in selected_keys for key in normalized_keys), dtype=bool)
        ordinals = np.rint((chunk["step"].to_numpy() - first_step) / step_increment).astype(np.int64)
        sampled = chunk[selected_mask & (ordinals % output_stride == 0)].copy()
        if sampled.empty:
            continue
        sampled_keys = list(zip(
            sampled["side"].astype(str), sampled["mode_index_1based"].astype(int)
        ))
        sampled["selection_reason"] = [";".join(reasons[key]) for key in sampled_keys]
        sampled_index = pd.MultiIndex.from_tuples(sampled_keys, names=KEY_COLUMNS)
        sampled["equilibrium_q_abs_rank_within_side"] = (
            summary["equilibrium_q_abs_rank_within_side"].reindex(sampled_index).to_numpy()
        )
        sampled["post_contact_delta_rms_rank_within_side"] = (
            summary["post_contact_delta_rms_rank_within_side"].reindex(sampled_index).to_numpy()
        )
        sampled[source_columns + metadata_columns].to_csv(
            selected_path,
            mode="a" if wrote_header else "w",
            header=not wrote_header,
            index=False,
            float_format="%.15e",
        )
        wrote_header = True
    if not wrote_header:
        raise ValueError("no rows were written to the selected time series")

    print(f"wrote {summary_path} ({len(summary)} modal rows)")
    print(f"wrote {selected_path} ({len(selected_keys)} selected side/mode pairs)")
    print(f"selected output interval: {actual_output_interval:.9g} s")


if __name__ == "__main__":
    main()
