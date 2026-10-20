#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
声門圧力分布：現行C++の断面出力と、声帯表面へ適用する圧力関数を比較する。

必要:
  numpy, matplotlib
  area.dat, section_x.dat, pressure.dat, separation.dat

推奨:
  gap_cubed.dat
  airflow_vt.dat  (これがあると論文式(9)による独立計算も実行する)

使用例:
  python pressure_compare.py --data-dir output/dat --out-dir pressure_comparison
  python pressure_compare.py --data-dir output/dat --steps 40000,50000,60000
  python pressure_compare.py --data-dir . --steps 0,5,10  # ユーザー提供の *_1016 ファイルにも対応

重要:
- 同じ声門形状と同じ Ug、同じ入口圧力で比較する「式の比較」であり、
  論文の計算を連成系ごと再現したものではない。
- pressure.dat は断面圧力 psurf[] であり、実際の各表面節点の圧力ではない。
  本コードでは、現行C++の pressureAtSurfaceX() と同じ連続圧力関数を再構成。
  FEM節点ごとの荷重積分には表面節点座標・面積ベクトル等が別途必要。
- 1Dモデルの圧力はAP方向一様と仮定する。
- 論文本文と公開Fortranの剥離判定/下流処理は完全には一致しない。
  そこで論文式による「剥離点まで」の曲線と、剥離後の仮定を分けて出力する。
"""
from __future__ import annotations

import argparse
import csv
import json
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


AREA_CLOSED_MM2 = 0.01  # C++: 1e-8 [m^2]
G3_CUTOFF_MM4 = 1e-12  # current forceCalculatorFluid.cpp
AREA_BERNOULLI_FLOOR_MM2 = 1e-3  # current forceCalculatorFluid.cpp


@dataclass
class SepRecord:
    step: int
    index: int
    x_sep_mm: float
    x_blend_end_mm: float
    p_sep_pa: float
    amin_mm2: float
    target_area_mm2: float
    sep_area_mm2: float
    fallback: bool


def read_table(path: Path, expected_ncols: int | None = None) -> dict[int, np.ndarray]:
    """#コメント付き「ステップ番号＋数値列」。行数の異なるファイルをステップで突き合わせる。"""
    arr = np.loadtxt(path, comments="#", ndmin=2)
    if arr.size == 0:
        raise ValueError(f"Empty file: {path}")
    if arr.shape[1] < 2:
        raise ValueError(f"Need step + data columns: {path}")
    if expected_ncols is not None and arr.shape[1] != expected_ncols:
        raise ValueError(f"{path}: expected {expected_ncols} columns, got {arr.shape[1]}")
    steps = arr[:, 0].astype(int)
    if not np.all(arr[:, 0] == steps):
        raise ValueError(f"Non-integer steps: {path}")
    if len(set(steps)) != len(steps):
        raise ValueError(f"Duplicate step IDs in {path}")
    return {int(k): row[1:].copy() for k, row in zip(steps, arr)}


def locate_file(directory: Path, stem: str, required: bool = True) -> Path | None:
    """標準名を優先。ユーザー提供の area_1016.dat なども自動判定。"""
    for ext in (".dat", ".txt"):
        p = directory / (stem + ext)
        if p.is_file():
            return p
    candidates = sorted(
        p for ext in (".dat", ".txt")
        for p in directory.glob(stem + "_*" + ext)
        if p.is_file()
    )
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple files for {stem}: {[str(p) for p in candidates]}. "
            f"Use a directory containing one set of output files."
        )
    if required:
        raise FileNotFoundError(f"{stem}.dat or {stem}.txt not found under {directory}")
    return None


def parse_separation(table: dict[int, np.ndarray]) -> dict[int, SepRecord]:
    result = {}
    for step, row in table.items():
        if len(row) != 8:
            raise ValueError("separation.dat must contain step + 8 columns")
        result[step] = SepRecord(
            step=step, index=int(row[0]), x_sep_mm=float(row[1]),
            x_blend_end_mm=float(row[2]), p_sep_pa=float(row[3]),
            amin_mm2=float(row[4]), target_area_mm2=float(row[5]),
            sep_area_mm2=float(row[6]), fallback=bool(int(row[7])),
        )
    return result


def reconstruct_applied_pressure(
    x_query_mm: np.ndarray,
    x_stations_mm: np.ndarray,
    p_stations_pa: np.ndarray,
    sep: SepRecord,
) -> np.ndarray:
    """
    現行C++ pressureAtSurfaceX() の構成式を再現。
    x<=x_sep: 断面圧力線形補間
    x_sep<x<x_blend: p_sepからp_downへsmoothstep
    x>=x_blend: p_down
    """
    xq = np.asarray(x_query_mm, dtype=float)
    x = np.asarray(x_stations_mm, dtype=float)
    p = np.asarray(p_stations_pa, dtype=float)
    if not (np.all(np.diff(x) > 0) and np.isfinite(x).all()):
        raise ValueError("section_x must be strictly increasing and finite")
    k = sep.index
    if k <= 0 or k >= len(x):
        raise ValueError(f"Invalid separation index {k}")
    p_down = p[-1]  # C++: all stations i>sepIndex are downstreamPressure
    p_before = np.interp(np.minimum(xq, sep.x_sep_mm), x[:k+1], p[:k+1])
    out = np.where(xq <= sep.x_sep_mm, p_before, p_down)
    mask = (xq > sep.x_sep_mm) & (xq < sep.x_blend_end_mm)
    blend_length = sep.x_blend_end_mm - sep.x_sep_mm
    if np.any(mask) and blend_length > 1e-12:
        s = np.clip((xq[mask] - sep.x_sep_mm) / blend_length, 0, 1)
        smooth = s*s*(3 - 2*s)
        out[mask] = (1-smooth)*sep.p_sep_pa + smooth*p_down
    return out


def paper_bern_visc_pressure(
    x_mm: np.ndarray,
    area_mm2: np.ndarray,
    ug_m3s: float,
    p0_pa: float,
    rho: float,
    mu: float,
    ap_length_mm: float,
    sep_index: int,
    height_scheme: str = "fortran",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Yoshinaga & Zhang 2025 Eq.(9), 各項の分解。
    height_scheme="fortran": 公開 calcForce.f90 の隣接 A[i-1],A[i] 平均。
    height_scheme="paper":   本文の heq=(A[i]+A[i+1])/(2 L_AP) を直接採用。
    下流はNaN（仮定を勝手に追加しない）。
    """
    n = len(x_mm)
    p = np.full(n, np.nan)
    dp_bern = np.full(n, np.nan)
    dp_visc = np.full(n, np.nan)
    if ug_m3s < 0:
        warnings.warn("Ug<0: original paper/Fortran approximation is intended for forward flow")
    p[0] = p0_pa
    if not (0 < sep_index < n):
        return p, dp_bern, dp_visc
    if np.min(area_mm2[:sep_index+1]) <= AREA_CLOSED_MM2:
        # 閉鎖時にはC++式・旧式とも別分岐。狭窄での0割回避。
        return p, dp_bern, dp_visc
    U = ug_m3s * 1e6  # For A in mm2: U/A has velocity units of m/s
    for i in range(1, sep_index+1):
        dx = x_mm[i] - x_mm[i-1]
        if dx <= 0:
            raise ValueError("non-positive section spacing")
        a0, a1 = area_mm2[i-1], area_mm2[i]
        bern = 0.5 * rho * U**2 * (1/a0**2 - 1/a1**2)
        if height_scheme == "fortran":
            a_avg = 0.5*(a0+a1)
        elif height_scheme == "paper":
            a_avg = 0.5*(a1 + area_mm2[min(i+1, n-1)])
        else:
            raise ValueError("height_scheme must be fortran or paper")
        h_eq_mm = a_avg / ap_length_mm
        # μ [Pa s] * Δx[mm] * U[mm³/s] / (L_AP[mm] h³[mm³])
        # multiply 1e3: mm³/s and mm lengths to SI
        visc = 12 * mu * dx * U / (ap_length_mm * h_eq_mm**3) * 1e3
        dp_bern[i] = bern
        dp_visc[i] = -visc
        p[i] = p[i-1] + bern - visc
    return p, dp_bern, dp_visc


def current_bern_visc_pressure(
    x_mm: np.ndarray,
    area_mm2: np.ndarray,
    g3_mm4: np.ndarray,
    ug_m3s: float,
    p0_pa: float,
    rho: float,
    mu: float,
    sep_index: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """forceCalculatorFluid.cpp の圧力積算・上昇時クランプを再現（開放時）。"""
    n = len(x_mm)
    p = np.full(n, np.nan)
    dp_bern = np.full(n, np.nan)
    dp_visc = np.full(n, np.nan)
    p[0] = p0_pa
    if np.min(area_mm2[1:-1]) <= AREA_CLOSED_MM2:
        return p, dp_bern, dp_visc
    U = ug_m3s * 1e6
    for i in range(1, sep_index+1):
        a0 = max(area_mm2[i-1], AREA_BERNOULLI_FLOOR_MM2)
        a1 = max(area_mm2[i], AREA_BERNOULLI_FLOOR_MM2)
        dx_mm = abs(x_mm[i]-x_mm[i-1])
        g3 = 0.5*(g3_mm4[i-1]+g3_mm4[i])
        bern = 0.5*rho*U**2*(1/a0**2-1/a1**2)
        visc = 12*mu*dx_mm*U/g3*1e3 if g3 > G3_CUTOFF_MM4 else 0.0
        dp_bern[i], dp_visc[i] = bern, -visc
        p[i] = p[i-1] + bern - visc
        if p[i] > p[i-1]:
            # C++: 圧力上昇時は直前の値に固定して残りの付着流区間を埋める
            p[i:sep_index+1] = p[i-1]
            break
    return p, dp_bern, dp_visc


def read_all(data_dir: Path):
    stems = ("area", "section_x", "pressure", "separation")
    paths = {stem: locate_file(data_dir, stem) for stem in stems}
    paths["gap_cubed"] = locate_file(data_dir, "gap_cubed", required=False)
    paths["airflow_vt"] = locate_file(data_dir, "airflow_vt", required=False)

    tables = {stem: read_table(path) for stem, path in paths.items() if path is not None}
    nsections = len(next(iter(tables["area"].values())))
    for name in ("section_x", "pressure"):
        for step, row in tables[name].items():
            if len(row) != nsections:
                raise ValueError(f"{name}: step {step} has {len(row)} stations; expected {nsections}")
    if "gap_cubed" in tables:
        for step, row in tables["gap_cubed"].items():
            if len(row) != nsections:
                raise ValueError(f"gap_cubed: step {step} has wrong station count")
    for stem in ("area", "section_x", "pressure"):
        if not tables[stem]:
            raise ValueError(f"No rows for {stem}")
    separators = parse_separation(tables["separation"])
    core_steps = set(tables["area"]) & set(tables["section_x"]) & set(tables["pressure"]) & set(separators)
    if not core_steps:
        raise ValueError("No common step among area, section_x, pressure, separation")
    return paths, tables, separators, sorted(core_steps), nsections


def prepare_snapshot(step, tables, separators, nsections, rho, mu, ap_mm,
                     paper_height_scheme):
    area = tables["area"][step]
    x = tables["section_x"][step]
    raw = tables["pressure"][step]
    sep = separators[step]
    if not (0 < sep.index < nsections-1):
        raise ValueError(f"Invalid separation index={sep.index}, step={step}")
    grid_x = np.linspace(x[0], x[-1], 800)
    p_applied = reconstruct_applied_pressure(grid_x, x, raw, sep)
    p_applied_stations = reconstruct_applied_pressure(x, x, raw, sep)
    p_down = float(raw[-1])

    paper = None
    current_recomputed = None
    if "airflow_vt" in tables and step in tables["airflow_vt"]:
        ug = float(tables["airflow_vt"][step][0])
        paper = paper_bern_visc_pressure(
            x, area, ug, float(raw[0]), rho, mu, ap_mm, sep.index,
            paper_height_scheme
        )
        if "gap_cubed" in tables and step in tables["gap_cubed"]:
            current_recomputed = current_bern_visc_pressure(
                x, area, tables["gap_cubed"][step], ug,
                float(raw[0]), rho, mu, sep.index
            )
    else:
        ug = None

    return dict(
        step=step, x=x, area=area, raw=raw, sep=sep, grid_x=grid_x,
        p_applied=p_applied, p_applied_stations=p_applied_stations,
        p_down=p_down, paper=paper, current_recomputed=current_recomputed,
        ug=ug,
    )


def save_snapshot_plot(out_path, snap, dt):
    x, area, raw, sep = snap["x"], snap["area"], snap["raw"], snap["sep"]
    xx, app = snap["grid_x"], snap["p_applied"]
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 8.0), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1.8]})
    ax = axes[0]
    ax.plot(x, area, "o-", markersize=2, lw=1.4, label="A(x) [mm2]")
    ax.axvline(sep.x_sep_mm, color="tab:red", ls="--", lw=1.4, label="x_sep")
    ax.axvline(sep.x_blend_end_mm, color="tab:orange", ls=":", lw=1.3,
               label="end of blend")
    ax.set_ylabel("Glottal area [mm2]")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    ax = axes[1]
    ax.plot(x, raw, "o", ms=3, color="tab:blue",
            label="pressure.dat: section pressure psurf")
    ax.plot(xx, app, "-", color="tab:orange", lw=2,
            label="Actual C++ wall-pressure law (reconstructed)")
    abrupt = np.where(xx <= sep.x_sep_mm, np.interp(np.minimum(xx, sep.x_sep_mm),
                                                   x[:sep.index+1], raw[:sep.index+1]),
                      snap["p_down"])
    ax.plot(xx, abrupt, ":", color="0.4", lw=1.5,
            label="Abrupt downstream-pressure cutoff (comparison)")
    if snap["paper"] is not None:
        paper = snap["paper"][0]
        ax.plot(x, paper, "s-", ms=3, lw=1.8, color="tab:green",
                label="Paper Eq.(9): attached-flow pressure")
    if snap["current_recomputed"] is not None:
        current_rebuilt = snap["current_recomputed"][0]
        ax.plot(x, current_rebuilt, "--", lw=1.5, color="tab:purple",
                label="Current formula recalculated from Ug/A/gap3")
    ax.axhline(0, color="0.45", lw=1)
    ax.axvline(sep.x_sep_mm, color="tab:red", ls="--", lw=1.4)
    ax.axvline(sep.x_blend_end_mm, color="tab:orange", ls=":", lw=1.3)
    ax.set_ylabel("Wall pressure (gauge) [Pa]")
    ax.set_xlabel("Flow-direction x [mm]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    suptitle = f"step={snap['step']:,}, t={snap['step']*dt:.6f} s"
    if snap["ug"] is not None:
        suptitle += f", Ug={snap['ug']*1e6:.3f} cm3/s"
    if sep.amin_mm2 <= AREA_CLOSED_MM2:
        suptitle += "  [NEAR-CLOSED]"
    fig.suptitle(suptitle)
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def save_csv(out_path, snap):
    x, area, raw, sep = snap["x"], snap["area"], snap["raw"], snap["sep"]
    app = snap["p_applied_stations"]
    p_paper = snap["paper"][0] if snap["paper"] else np.full(len(x), np.nan)
    p_current = (snap["current_recomputed"][0] if snap["current_recomputed"]
                 else np.full(len(x), np.nan))
    rows = zip(range(len(x)), x, area, raw, app, p_paper, p_current)
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section_index_0based", "x_mm", "area_mm2",
                         "pressure_dat_pa", "current_applied_wall_pressure_pa",
                         "paper_eq9_presep_pa", "current_formula_presep_pa"])
        writer.writerows(rows)


def save_overview(out_path, common_steps, tables, separators, dt):
    steps = np.asarray(common_steps, dtype=int)
    times = steps * dt
    minareas = np.array([separators[s].amin_mm2 for s in steps])
    psep = np.array([separators[s].p_sep_pa for s in steps])
    pdown = np.array([tables["pressure"][s][-1] for s in steps])
    pmin = np.array([np.min(tables["pressure"][s]) for s in steps])
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(times, minareas)
    axes[0].axhline(AREA_CLOSED_MM2, ls="--", color="0.4", label="closure threshold")
    axes[0].set_ylabel("Min area [mm2]")
    axes[0].legend(fontsize=9)
    axes[1].plot(times, psep, label="p_sep (separation.dat)")
    axes[1].plot(times, pdown, label="p_down (last pressure station)")
    axes[1].plot(times, pmin, label="min pressure.dat")
    axes[1].axhline(0, color="0.5", lw=.8)
    axes[1].set_ylabel("Pressure [Pa]")
    axes[1].legend(fontsize=9)
    axes[2].plot(times, [separators[s].x_sep_mm for s in steps],
                 label="x_sep")
    axes[2].plot(times, [separators[s].x_blend_end_mm for s in steps],
                 label="x_blend_end")
    axes[2].set_ylabel("Flow-direction x [mm]")
    axes[2].set_xlabel("Time [s]")
    axes[2].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.3)
    fig.suptitle("Simulation pressure and separation overview")
    fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def closest_step(target, valid_steps):
    return min(valid_steps, key=lambda st: abs(st-target))


def main():
    p = argparse.ArgumentParser(
        description="現行声門圧力の可視化とYoshinaga & Zhang(2025)式(9)の同一形状比較"
    )
    p.add_argument("--data-dir", type=Path, default=Path("output/dat"),
                   help="*.dat出力があるディレクトリ。*_1016.txt なども検出可能")
    p.add_argument("--out-dir", type=Path, default=Path("pressure_comparison"))
    p.add_argument("--steps", default="latest",
                   help="例: latest または 40000,50000,60000。近い保存ステップを選ぶ")
    p.add_argument("--times", default=None,
                   help="例: 0.8,0.802,0.804 [s]。指定すれば--stepsを上書き")
    p.add_argument("--dt", type=float, default=1.0e-5, help="時間刻み [s]")
    p.add_argument("--rho", type=float, default=1.13, help="空気密度 [kg/m3]")
    p.add_argument("--mu", type=float, default=1.83e-5, help="粘度 [Pa s]")
    p.add_argument("--ap-mm", type=float, default=17.0, help="声門AP長さ [mm]")
    p.add_argument("--paper-height", choices=["fortran", "paper"],
                   default="fortran", help="式(9)粘性項heqの添字規則")
    args = p.parse_args()

    if args.dt <= 0 or args.rho <= 0 or args.mu <= 0 or args.ap_mm <= 0:
        p.error("dt, rho, mu, ap-mm must all be positive")

    paths, tables, sep, common_steps, nsec = read_all(args.data_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print("INPUT:")
    for name, path in paths.items():
        print(f"  {name:<14}: {path if path else '[missing / optional]'}")
    print(f"\nCOMMON STEPS: {len(common_steps)} | {common_steps[0]} .. {common_steps[-1]}")
    print(f"SECTIONS: {nsec}")
    if len(common_steps) <= 3:
        print("NOTE: these appear to be short sample files. Full output is needed to assess phonation.")

    if args.times is not None:
        try:
            requested = [round(float(v.strip()) / args.dt)
                         for v in args.times.split(",") if v.strip()]
        except ValueError as err:
            raise ValueError("--times must be comma-separated seconds") from err
        if not requested:
            raise ValueError("--times is empty")
    elif args.steps.lower() == "latest":
        requested = [common_steps[-1]]
    else:
        try:
            requested = [int(x.strip()) for x in args.steps.split(",") if x.strip()]
        except ValueError as err:
            raise ValueError("--steps must be latest or comma-separated integers") from err
        if not requested:
            raise ValueError("--steps is empty")
    selected = []
    for st in requested:
        matched = closest_step(st, common_steps)
        if matched != st:
            print(f"NOTE: requested step {st} -> closest common step {matched}")
        if matched not in selected:
            selected.append(matched)

    summaries = []
    for step in selected:
        snap = prepare_snapshot(step, tables, sep, nsec, args.rho, args.mu, args.ap_mm,
                                args.paper_height)
        prefix = f"step_{step:07d}"
        save_snapshot_plot(args.out_dir / f"{prefix}.png", snap, args.dt)
        save_csv(args.out_dir / f"{prefix}.csv", snap)
        info = asdict(snap["sep"])
        info.update({
            "time_s": step*args.dt,
            "p_down_pa": snap["p_down"],
            "p_station_min_pa": float(np.min(snap["raw"])),
            "p_station_max_pa": float(np.max(snap["raw"])),
            "p_wall_reconstructed_min_pa": float(np.min(snap["p_applied"])),
            "p_wall_reconstructed_max_pa": float(np.max(snap["p_applied"])),
            "flow_rate_m3_s": snap["ug"],
            "paper_comparison_available": snap["paper"] is not None,
        })
        summaries.append(info)
        print(f"\nstep={step}, t={step*args.dt:.6f}s, minArea={info['amin_mm2']:.6g}mm2, "
              f"psep={info['p_sep_pa']:.6g}Pa, pdown={info['p_down_pa']:.6g}Pa")
        if snap["paper"] is None:
            print("  Paper Eq.(9) NOT recalculated: airflow_vt.dat not supplied or step absent")
        if info["amin_mm2"] <= AREA_CLOSED_MM2:
            print("  NEAR-CLOSED: do not interpret this snapshot as attached-flow phonation")

    save_overview(args.out_dir / "overview.png", common_steps, tables, sep, args.dt)
    manifest = {
        "input_files": {k: str(v) if v is not None else None for k, v in paths.items()},
        "n_common_steps": len(common_steps), "n_sections": nsec,
        "parameters": {
            "dt_s": args.dt, "rho_kg_m3": args.rho, "mu_pa_s": args.mu,
            "ap_length_mm": args.ap_mm, "paper_height": args.paper_height,
            "closure_threshold_mm2": AREA_CLOSED_MM2
        },
        "snapshots": summaries,
        "warnings": [
            "The pressure.dat file contains sections, not per-node loads.",
            "The reconstructed C++ wall pressure is a continuous function of x; "
            "node positions/normals are not supplied.",
            "Paper Eq.(9) is evaluated under the C++ flow and inlet pressure; "
            "this is NOT a complete coupled paper-model simulation.",
            "Paper post-separation constant-pressure choice is not fitted or identified here.",
            "C++ per-interval valid-section flags were not exported; pressure.dat "
            "is the reliable existing output for current pressure.",
        ],
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"\nOUTPUT: {args.out_dir.resolve()}")
    print("  overview.png, step_XXXXXXX.png, step_XXXXXXX.csv, summary.json")


if __name__ == "__main__":
    main()
