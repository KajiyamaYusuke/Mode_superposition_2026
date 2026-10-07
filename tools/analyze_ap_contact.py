#!/usr/bin/env python3
"""Compare AP-contact diagnostic runs and build hypothesis-oriented outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal

from output_layout import output_path, result_path


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_a", type=Path, help="baseline/live run")
    parser.add_argument("run_b", type=Path, help="comparison run")
    parser.add_argument("--reference-window", type=float, nargs=2, default=(0.25, 0.30))
    parser.add_argument("--perturbation-time", type=float, default=0.35)
    parser.add_argument("--fit-window", type=float, nargs=2, default=(0.37, 0.60))
    parser.add_argument("--output-directory", type=Path)
    return parser.parse_args()


def csv_file(run: Path, name: str) -> Path:
    return result_path(run, name)


def load(run: Path, name: str, required: bool = True) -> pd.DataFrame:
    path = csv_file(run, name)
    if not path.is_file():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path)


def fit_growth(time: np.ndarray, values: np.ndarray, start: float, end: float) -> dict[str, float]:
    mask = np.isfinite(time) & np.isfinite(values) & (time >= start) & (time <= end)
    time, values = time[mask], values[mask]
    if len(time) < 16:
        return {"sigma_per_s": np.nan, "frequency_hz": np.nan, "fit_rmse": np.nan,
                "amplitude_floor": np.nan}
    detrended = signal.detrend(values)
    envelope = np.abs(signal.hilbert(detrended))
    floor = max(np.max(envelope) * 1.0e-8, np.finfo(float).tiny)
    trim = max(1, int(0.1 * len(time)))
    fit = np.arange(len(time)) >= trim
    fit &= np.arange(len(time)) < len(time) - trim
    fit &= envelope > floor
    if np.count_nonzero(fit) < 8:
        sigma = rmse = np.nan
    else:
        coefficients = np.polyfit(time[fit], np.log(envelope[fit]), 1)
        predicted = np.polyval(coefficients, time[fit])
        sigma = float(coefficients[0])
        rmse = float(np.sqrt(np.mean((np.log(envelope[fit]) - predicted) ** 2)))
    dt = float(np.median(np.diff(time)))
    frequency = np.fft.rfftfreq(len(time), dt)
    amplitude = np.abs(np.fft.rfft(detrended * signal.windows.hann(len(time))))
    usable = (frequency > 0.0) & (frequency <= 0.49 / dt)
    dominant = float(frequency[usable][np.argmax(amplitude[usable])]) if np.any(usable) else np.nan
    return {"sigma_per_s": sigma, "frequency_hz": dominant, "fit_rmse": rmse,
            "amplitude_floor": float(floor)}


def trapezoid_cumulative(time: np.ndarray, power: np.ndarray) -> np.ndarray:
    result = np.zeros(len(time))
    if len(time) > 1:
        result[1:] = np.cumsum(0.5 * np.diff(time) * (power[:-1] + power[1:]))
    return result


def energy_outputs(run: Path, modal: pd.DataFrame, reference: tuple[float, float],
                   perturbation_time: float) -> pd.DataFrame:
    native = load(run, "perturbation_energy.csv", required=False)
    if not native.empty and "integration_method" in native.columns:
        return native
    rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, float | str | int]] = []
    for (side, mode), group in modal.groupby(["side", "mode_index_1based"]):
        group = group.sort_values("state_time_s").copy()
        base = group[(group.load_time_s >= reference[0]) & (group.load_time_s <= reference[1])]
        post = group[group.load_time_s > perturbation_time].copy()
        if base.empty or post.empty:
            continue
        qbar = float(base.q.mean())
        vbar = float(base.qdot.mean())
        ffluid = float(base.fluid_force.mean())
        contact_columns = ["contact_force_ap_low", "contact_force_ap_high", "contact_force_interior"]
        fcontact = float(base[contact_columns].sum(axis=1).mean())
        omega = 2.0 * np.pi * float(post.frequency_hz.iloc[0])
        dq = post.q.to_numpy() - qbar
        dv = post.qdot.to_numpy() - vbar
        df_fluid = post.fluid_force.to_numpy() - ffluid
        contact = post[contact_columns].sum(axis=1).to_numpy()
        df_contact = contact - fcontact
        time = post.state_time_s.to_numpy()
        epert = 0.5 * (dv * dv + omega * omega * dq * dq)
        pfluid = df_fluid * dv
        pcontact = df_contact * dv
        damping_coefficient = np.divide(
            post.structural_damping_force.to_numpy(), post.qdot.to_numpy(),
            out=np.zeros(len(post)), where=np.abs(post.qdot.to_numpy()) > 1.0e-30)
        pdamp = damping_coefficient * dv * dv
        wfluid = trapezoid_cumulative(time, pfluid)
        wcontact = trapezoid_cumulative(time, pcontact)
        wdamp = trapezoid_cumulative(time, pdamp)
        balance = epert - epert[0] - wfluid - wcontact + wdamp
        post["dq"] = dq
        post["dv"] = dv
        post["epert_generalized"] = epert
        post["pfluid_pert_generalized_per_s"] = pfluid
        post["pcontact_pert_generalized_per_s"] = pcontact
        post["pdamp_pert_generalized_per_s"] = pdamp
        post["wfluid_pert_generalized"] = wfluid
        post["wcontact_pert_generalized"] = wcontact
        post["wdamp_pert_generalized"] = wdamp
        post["balance_residual_generalized"] = balance
        rows.append(post)
        summary_rows.append({
            "side": side, "mode_index_1based": int(mode),
            "q_bar": qbar, "v_bar": vbar, "fluid_force_bar": ffluid,
            "contact_force_bar": fcontact, "final_epert_generalized": float(epert[-1]),
            "final_wfluid_generalized": float(wfluid[-1]),
            "final_wcontact_generalized": float(wcontact[-1]),
            "final_wdamp_generalized": float(wdamp[-1]),
            "final_balance_residual_generalized": float(balance[-1]),
        })
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    out_dir = run / "csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "perturbation_energy.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(out_dir / "perturbation_energy_summary.csv", index=False)
    return result


def summarize_run(label: str, run: Path, reference: tuple[float, float],
                  perturbation: float, fit_window: tuple[float, float]) -> tuple[list[dict], pd.DataFrame]:
    modal = load(run, "modal_force_split.csv")
    flow = load(run, "flow_diagnostics_v2.csv")
    contact = load(run, "contact_region_summary.csv", required=False)
    reference_flow = flow[(flow.time_s >= reference[0]) & (flow.time_s <= reference[1])]
    common = {
        "case": label,
        "run_directory": str(run.resolve()),
        "equilibrium_min_area_mm2": float(reference_flow.min_area_mm2.mean()),
        "equilibrium_pg_pa": float(reference_flow.current_pg_pa.mean()),
        "equilibrium_flow_m3_s": float(reference_flow.flow_rate_m3_s.mean()),
        "separation_fallback_fraction": float(reference_flow.separation_used_fallback.mean()),
        "separation_index_std": float(reference_flow.separation_index.std(ddof=0)),
        "max_equation_residual": float(np.max(np.abs(modal.equation_residual))),
    }
    if not contact.empty:
        ref_contact = contact[(contact.time_s >= reference[0]) & (contact.time_s <= reference[1])]
        for region, values in ref_contact.groupby("region"):
            common[f"{region}_mean_pair_force_n"] = float(values.sum_pair_force_magnitude_n.mean())
            common[f"{region}_max_penetration_mm"] = float(values.max_penetration_mm.max())
            common[f"{region}_negative_gap_fraction"] = float(values.negative_gap_fraction.mean())
    rows: list[dict] = []
    mid = 0.5 * (fit_window[0] + fit_window[1])
    for (side, mode), group in modal.groupby(["side", "mode_index_1based"]):
        group = group.sort_values("state_time_s")
        first = fit_growth(group.state_time_s.to_numpy(), group.q.to_numpy(), *fit_window)
        second_window = (fit_window[0] + 0.1 * (fit_window[1] - fit_window[0]), fit_window[1])
        second = fit_growth(group.state_time_s.to_numpy(), group.q.to_numpy(), *second_window)
        status = "OK"
        if np.isfinite(first["sigma_per_s"]) and np.isfinite(second["sigma_per_s"]):
            if np.sign(first["sigma_per_s"]) != np.sign(second["sigma_per_s"]):
                status = "INCONCLUSIVE_WINDOW_DEPENDENT"
        rows.append({**common, "side": side, "mode_index_1based": int(mode),
                     **first, "sigma_alt_window_per_s": second["sigma_per_s"],
                     "growth_assessment": status, "fit_window_mid_s": mid})
    return rows, energy_outputs(run, modal, reference, perturbation)


def plot_outputs(out: Path, runs: list[tuple[str, Path]], summary: pd.DataFrame,
                 energies: dict[str, pd.DataFrame], perturbation: float) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for label, run in runs:
        detail = load(run, "contact_pairs_detail.csv", required=False)
        if detail.empty:
            region = load(run, "contact_region_summary.csv", required=False)
            if not region.empty:
                grouped = region.groupby("region").max_penetration_mm.max()
                ax.scatter([label] * len(grouped), grouped, label=label)
        else:
            active = detail[detail.excluded == 0]
            ax.scatter(active.x_mm, active.ap_j_0based, s=8,
                       c=active.penetration_mm, cmap="viridis", alpha=0.7, label=label)
    ax.set(xlabel="x [mm] (or case when detail logging is disabled)", ylabel="AP j index",
           title="AP contact location / penetration")
    ax.grid(alpha=0.2)
    fig.savefig(out / "contact_map.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for label, group in summary.groupby("case"):
        x = np.arange(len(group))
        ax.scatter(x, group.sigma_per_s, label=label)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set(xlabel="side/mode row", ylabel="growth rate sigma [1/s]", title="Modal growth")
    ax.legend(); ax.grid(alpha=0.2)
    fig.savefig(out / "modal_growth.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for label, energy in energies.items():
        if energy.empty:
            continue
        total = energy.groupby("state_time_s")[["wfluid_pert_generalized",
            "wcontact_pert_generalized", "wdamp_pert_generalized"]].sum()
        ax.plot(total.index, total.wfluid_pert_generalized, label=f"{label}: fluid")
        ax.plot(total.index, total.wcontact_pert_generalized, linestyle="--",
                label=f"{label}: contact")
    ax.axvline(perturbation, color="black", linewidth=0.8)
    ax.set(xlabel="time [s]", ylabel="cumulative generalized work",
           title="Perturbation work (not asserted as J)")
    ax.legend(ncol=2, fontsize=8); ax.grid(alpha=0.2)
    fig.savefig(out / "perturbation_work.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(len(runs), 1, figsize=(9, 3.5 * len(runs)), squeeze=False,
                             constrained_layout=True)
    for axis, (label, run) in zip(axes[:, 0], runs):
        area_path, pressure_path = result_path(run, "area.dat"), result_path(run, "pressure.dat")
        if not area_path.is_file() or not pressure_path.is_file():
            continue
        area = np.loadtxt(area_path, comments="#", ndmin=2)
        pressure = np.loadtxt(pressure_path, comments="#", ndmin=2)
        flow = load(run, "flow_diagnostics_v2.csv", required=False)
        perturb_step = (float(flow.loc[flow.time_s >= perturbation, "step"].iloc[0])
                        if not flow.empty and np.any(flow.time_s >= perturbation)
                        else area[0, 0])
        candidates = np.flatnonzero(area[:, 0] >= perturb_step)
        if not len(candidates):
            candidates = np.arange(len(area))
        pressure_axis = axis.twinx()
        for phase, index in enumerate(np.linspace(candidates[0], candidates[-1], 4).astype(int)):
            color = f"C{phase}"
            axis.plot(np.arange(area.shape[1] - 1), area[index, 1:],
                      color=color, label=f"A step {int(area[index, 0])}")
            pressure_index = int(np.argmin(np.abs(pressure[:, 0] - area[index, 0])))
            pressure_axis.plot(np.arange(pressure.shape[1] - 1), pressure[pressure_index, 1:],
                               color=color, linestyle="--", alpha=0.65,
                               label=f"p step {int(pressure[pressure_index, 0])}")
        axis.set(title=f"{label}: four sampled flow-area profiles", xlabel="flow section (0-based)",
                 ylabel="area [mm²]")
        pressure_axis.set_ylabel("pressure [Pa]")
        handles, labels = axis.get_legend_handles_labels()
        handles_p, labels_p = pressure_axis.get_legend_handles_labels()
        axis.legend(handles + handles_p, labels + labels_p, fontsize=7, ncol=2)
        axis.grid(alpha=0.2)
    fig.savefig(out / "flow_phase_profiles.png", dpi=170)
    plt.close(fig)


def main() -> None:
    args = arguments()
    reference = tuple(args.reference_window)
    fit_window = tuple(args.fit_window)
    if not reference[0] < reference[1] <= args.perturbation_time < fit_window[0] < fit_window[1]:
        raise ValueError("require reference_start < reference_end <= perturbation < fit_start < fit_end")
    run_a, run_b = args.run_a.resolve(), args.run_b.resolve()
    out = (args.output_directory.resolve() if args.output_directory
           else run_b.parent / f"comparison_{run_a.name}_vs_{run_b.name}")
    out.mkdir(parents=True, exist_ok=True)
    rows_a, energy_a = summarize_run("A", run_a, reference, args.perturbation_time, fit_window)
    rows_b, energy_b = summarize_run("B", run_b, reference, args.perturbation_time, fit_window)
    summary = pd.DataFrame(rows_a + rows_b)
    modal_a = load(run_a, "modal_force_split.csv")
    modal_b = load(run_b, "modal_force_split.csv")
    keys = ["step", "side", "mode_index_1based"]
    prefix_a = modal_a[modal_a.load_time_s < args.perturbation_time]
    prefix_b = modal_b[modal_b.load_time_s < args.perturbation_time]
    prefix = prefix_a.merge(prefix_b, on=keys, suffixes=("_a", "_b"))
    state_columns = ("q", "qdot", "qddot", "fluid_force", "total_force")
    prefix_error = (max(
        float(np.max(np.abs(prefix[f"{name}_a"] - prefix[f"{name}_b"])))
        for name in state_columns) if not prefix.empty else np.nan)
    summary["branch_prefix_max_abs_state_load_difference"] = prefix_error
    summary.to_csv(out / "comparison_summary.csv", index=False)
    plot_outputs(out, [("A", run_a), ("B", run_b)], summary,
                 {"A": energy_a, "B": energy_b}, args.perturbation_time)
    inconclusive = bool((summary.growth_assessment != "OK").any())
    native_energy = all(
        not energy.empty and "integration_method" in energy.columns
        and energy.integration_method.astype(str).str.startswith("internal_dt").all()
        for energy in (energy_a, energy_b)
    )
    report = ["# AP端部接触・自己発振診断", "",
              f"- A: `{run_a}`", f"- B: `{run_b}`", "",
              "## 判定", "",
              "- H1: A/B の平衡面積、AP端部接触モード力を比較してください。",
              "- H2: live/freeze_contact の sigma 差が必要です。A/B がその組合せでなければ未確定です。",
              "- H3: gap_field.csv、contact_map.png、fixed_boundary_mode_audit.csv を照合してください。",
              "- H4: freeze_fluid と差分流体仕事が必要です。",
              "- H5: separation_fallback_fraction と flow_phase_profiles.png を確認してください。",
              "",
              f"総合フラグ: {'INCONCLUSIVE（fit窓依存あり）' if inconclusive else '数値表を比較可能'}", "",
              f"分岐前の構造・荷重最大差: `{prefix_error:.6e}`", "",
              "エネルギーは質量正規化のSI整合を未確認のため generalized units であり、J/Wとは断定しません。",
              ("仕事はソルバ内部dtで台形積分されています。"
               if native_energy else
               "旧出力のため仕事をmodal_force_split.csvの間隔で再積分しており、内部dt積分より粗い点に注意してください。")]
    (out / "hypothesis_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"output_directory": str(out), "rows": len(summary)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
