# Dynamics Analysis Tools

These tools classify simulation behavior without changing the C++ physical,
contact, fluid, or integration models. Classifications are exploratory
heuristics: `scattered_or_chaotic_candidate` is not proof of chaos.

## Installation

```bash
python3 -m pip install -r requirements-analysis.txt
```

## Analyze One Run

Run the simulator from `build/`, then analyze the latest output:

```bash
cd build
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./simulation ../input/param.txt
cd ..
python3 tools/analyze_dynamics.py \
  --run-dir output --start-time 0.15 --last-cycles 50
```

## Output layout

Generated results are grouped by file type instead of being written directly
under `output/`:

- `output/dat/`: solver time histories (`area.dat`, `displace.dat`, pressures, flow)
- `output/csv/`: tabular diagnostics and analysis results
- `output/png/`: plots (irregularity plots use `output/png/irregularity/`)
- `output/txt/`: `manifest.txt`, `params_used.txt`, and text logs
- `output/wav/`: synthesized and reference audio
- `output/json/`: analysis summaries

Analysis commands still take `--run-dir output`; maintained analysis tools find
the appropriate subdirectory automatically and can also read legacy flat runs.

For contact-coefficient tuning, compare `csv/contact_iteration_debug.csv`,
`csv/contact_coupling_debug.csv`, and `csv/contact_final_state.csv`. The last
file measures penetration after the final contact load has been applied.

## Steinecke--Herzel sweep classification

`python3 tools/spectol_sweep.py` now classifies the final 0.15 s of each
successful pressure/damping case using the unreduced right:left peak-count
label (`1:1`, `2:2`, `5:8`, etc.). It writes the sweep table to
`output/csv/steinecke_classification_summary.csv`, detailed peak/count-window/
period-error/map/spectrum tables under `output/csv/classification/`, and one
four-panel verification plot per condition under `output/png/classification/`.
`uncertain` and candidate labels are intentional: finite-period locking is not
asserted when peak-count, next-maximum-map, left-period, or subharmonic checks
conflict.

## Stiffness-condition comparison

Calculate left/right displacement jitter, shimmer, frequency/amplitude ratios,
per-side standard deviations, and bootstrap uncertainty for several completed
runs in one command:

```bash
python3 tools/compare_stiffness_metrics.py \
  baseline=analysis_runs/baseline \
  stiff_left=analysis_runs/stiff_left \
  stiff_right=analysis_runs/stiff_right
```

With no run arguments it analyzes the current `output/`. Results are written
to `output/csv/stiffness_comparison_metrics.csv`, cycle-level values to
`output/csv/stiffness_cycle_metrics.csv`, and a four-panel comparison to
`output/png/stiffness_comparison.png`. The default analysis interval is the
last 0.15 s; change it with `--duration`.

The analyzer reads `dt` from `manifest.txt` unless `--dt` is supplied. It
interprets the first column of `area.dat` as a step number, produces
`metrics.json`, `cycle_metrics.csv`, `poincare_points.csv`, `return_map.csv`,
and diagnostic figures under `png/`. Always inspect
`png/time_series.png` to confirm that detected events match the waveform.

The same command independently detects peaks in the `uyL_mm` and `uyR_mm`
columns of `displace.dat`. It writes:

- `left_right_ratio_metrics.csv`: representative periods, amplitudes, ratios,
  dominant frequencies, and 2:1 candidate classes.
- `left_cycle_metrics.csv` and `right_cycle_metrics.csv`: raw periods,
  validity flags, peaks, troughs, and peak-to-trough amplitudes.
- `png/left_right_displacement.png`,
  `left_right_cycle_period.png`, and `left_right_cycle_amplitude.png`.

Defaults can be adjusted without breaking the existing command:

```bash
python3 tools/analyze_dynamics.py --run-dir output \
  --analysis-start 0.18 \
  --peak-prominence-ratio 0.10 \
  --period-ratio-target 2.0 --period-ratio-tolerance 0.20 \
  --amplitude-ratio-target 2.0 --amplitude-ratio-tolerance 0.20
```

The unsigned ratios are always at least one. The signed `left_over_right`
columns retain which side has the longer period or larger amplitude.

Plot any available Poincaré variables:

```bash
python3 tools/plot_poincare.py --run-dir output --list-columns
python3 tools/plot_poincare.py --run-dir output \
  --x qL_1 --y qdotL_1 --connect
```

## Pressure Sweep

Each pressure starts from the same initial state. The base parameter file is
never modified, and every condition is archived separately.

```bash
python3 tools/run_pressure_sweep.py \
  --config tools/pressure_sweep_config.json \
  --start 500 --stop 3000 --step 100 --analysis-start 0.15
```

Use `--dry-run` to inspect generated parameters, `--sweep-dir` to resume an
existing directory, and `--rerun-completed` to replace completed cases.
Failures are recorded in `status.json`, `stdout.log`, and `stderr.log`.
The sweep also appends left/right fields to `sweep_summary.csv`, writes
`pressure_sweep_left_right_ratios.csv`, and creates
`figures/pressure_vs_period_ratio.png` and
`figures/pressure_vs_amplitude_ratio.png`. Horizontal reference lines mark
ratios 1 and 2.

```bash
python3 tools/plot_bifurcation.py \
  --sweep-dir analysis_runs/pressure_sweep_YYYYMMDD_HHMMSS \
  --quantity peak_area_mm2
```

Peak/trough diagrams plot every retained cycle, not only averages.

## Interpretation and Validation

Remove pressure-ramp and startup transients and retain at least 20 cycles for
classification (50+ for bifurcation diagrams). Repeat candidates with smaller
`dt`, different output resolution, mode count, contact iteration limit, and
fixed OpenMP settings. A feature that disappears under these checks may be a
sampling artifact, unresolved transient, or numerical instability.

Run artificial-signal tests with:

```bash
python3 -m unittest discover -s tools/tests -v
```

## Minimal modal-work exit test

Compare one `live` run with the matching `freeze_fluid` run for modes 1, 4,
and 5:

```bash
python3 tools/analyze_modal_exit_test.py \
  diagnostic_runs/example/A_live \
  diagnostic_runs/example/A_freeze_fluid \
  --modes 1,4,5
```

The freeze/perturbation time, reference window, time step, damping ratios, and
run identity are read from each run's manifests. Use `--switch-time`,
`--reference-window START END`, and one or two `--window START END` options to
override them. With no `--window`, two adjacent post-switch windows are derived
from the recorded switch and common data duration; no absolute freeze time is
embedded in the tool.

The tool reads `modal_force_split.csv` and `perturbation_energy.csv` in chunks,
retaining only the selected modes. It separates applied fluid force from the
total applied load, reuses the solver's internal-dt cumulative work when
available, and compares it against two integrations of the decimated data.
Outputs are `modal_exit_summary.csv`, `modal_exit_response.png`,
`modal_exit_work.png`, and `modal_exit_report.md` under a sibling
`modal_exit_test/` directory unless `--output-directory` is supplied. Work is
reported in generalized mass-normalized units, not labelled J/W.

## Irregular-vibration indicator export

The simulator writes `output/irregularity_timeseries.csv` at the same cadence
as `area.dat`. It contains synchronized time, left/right monitor displacement,
minimum and maximum glottal area, flow rate, outlet pressure, and contact state.
After removing the startup transient, calculate the indicators described in
`vocal_fold_irregularity_metrics.md` with:

```bash
python3 tools/analyze_irregularity.py --run-dir output --start-time 0.15
```

`--run-dir` may be omitted. The default is the repository's `output/`
directory even when the command is launched from inside `tools/`.

This analyzer intentionally does not assign Type A--D. It saves one-row
summaries (`irregularity_metrics.csv` and `.json`) plus peak sequences, Hilbert
phase, spectra, cycle-lag repeat errors, and autocorrelation as separate CSVs.
The figures are written under `output/figures/irregularity/`. Adjust
`--prominence-ratio` after checking the peak markers in
`timeseries_and_peaks.png`; the default is 5% of each signal's range.

## Modal-energy share history

Use the saved mass-normalized modal coordinates to compare how the structural
energy is distributed among modes:

```bash
python3 tools/plot_modal_energy_share.py --run-dir output --start-time 0.15
```

For each side, the script ranks modes by their time-averaged energy-like
quantity and plots the leading five energy shares relative to all retained
modes. Both sides appear in one figure,
`output/figures/modal_energy_share_top5.png`. The complete time history and
ranking are saved as `modal_energy_share_timeseries.csv` and
`modal_energy_share_summary.csv`. Use `--top-count` to change the number of
displayed modes.

## Self-oscillation diagnostics

The solver writes `csv/flow_diagnostics_v2.csv`, `csv/modal_energy.csv`, and
`csv/energy_summary.csv` at `diagnosticOutputIntervalSteps` (default: 5).
`modal_energy.csv` separates the contact-free fluid generalized force from the
total force used by the integrator. Its force, power, work, and energy values
are in **mass-normalized coordinate units**; they must not be interpreted as
SI joules until the mode normalization has been independently verified.

Optional trailing `key=value` parameters are:

```text
initialGapMm=0.0
pressureRampTimeSec=0.10
diagnosticOutputIntervalSteps=5
perturbationEnabled=0
perturbationTimeSec=0.20
perturbationModeIndex=5
perturbationAmplitudeAtProbeMm=0.001
perturbationPattern=symmetric_opening
```

The repository's historical implementation actually used a 0.05 s ramp even
though the diagnostic specification described it as 0.10 s. To preserve
bitwise behavior for legacy parameter files, omission of
`pressureRampTimeSec` retains 0.05 s; new diagnostic and sweep inputs write
`pressureRampTimeSec=0.10` explicitly.

A negative `initialGapMm` is diagnostic geometric precompression, not a
physical adduction force. The perturbation mode is one-based. The requested
probe displacement increment is applied once in opposite opening directions.

Analyze a completed run with:

```bash
python tools/analyze_self_oscillation.py output \
  --analysis-start 0.15 --analysis-end 0.50 --tail-duration 0.10
```

This creates `json/analysis_summary.json`, `csv/analysis_summary.csv`, and
`png/diagnostic_overview.png`. Missing optional inputs are reported as
warnings and do not prevent the available signals from being analyzed.

Run the pressure-alignment and initial-gap sweep with:

```bash
python tools/run_diagnostic_sweep.py --dry-run
python tools/run_diagnostic_sweep.py --jobs 1 --omp-threads 4
```

The first stage runs 2900, 3300, 3600, and 4000 Pa and selects the case whose
tail-average `currentPg` is closest to 2900 Pa. The second stage runs 0.0,
-0.1, and -0.2 mm at that pressure with the perturbation enabled. Each
condition has an independent directory and records its input, command, Git
commit, UTC timestamps, and exit status. Use `--jobs` cautiously: every
process also uses OpenMP threads.

## AP端部接触診断

通常計算を変えずにログだけ有効化する最小設定は次のとおりです。

```text
apContactDiagnosticEnabled=1
apContactAnalysisRowsPerEnd=3
apContactAnalysisDistanceMm=0
apContactExcludedRowsPerEnd=0
apContactExcludedDistanceMm=0
contactDetailStartSec=0.349
contactDetailEndSec=0.351
diagnosticLoadMode=live
diagnosticReferenceWindowStartSec=0.25
diagnosticReferenceWindowEndSec=0.30
```

行数指定と距離指定は同時に使えません。`apContactExcludedRowsPerEnd=3`
は `j<3` と `j>=N_AP-3` で検出された接触をログに残したまま、適用反力を
0にします。詳細区間をともに負値にするとペア詳細ログを無効化します。

A/Bケースと荷重凍結ケースは既存 `output/` を上書きせずに実行できます。

```bash
python tools/run_ap_contact_diagnostics.py --dry-run
python tools/run_ap_contact_diagnostics.py --stage ab --omp-threads 4
python tools/run_ap_contact_diagnostics.py --stage freeze --omp-threads 4
```

`--stage freeze` は `live`, `freeze_contact`, `freeze_fluid`, `freeze_all`
を同じ入力から決定的に再実行します。凍結値は基準窓の非ゼロ平均荷重で、
切替は擾乱時刻です。各ケースの `run_manifest.json` と最上位の
`run_summary.json`（失敗時は `failures.json`）を確認してください。

2ケースの比較例です。

```bash
python tools/analyze_ap_contact.py \
  diagnostic_runs/ap_contact_YYYYMMDD_HHMMSS/A_live_contact \
  diagnostic_runs/ap_contact_YYYYMMDD_HHMMSS/B_exclude_ap_ends \
  --reference-window 0.25 0.30 \
  --perturbation-time 0.35 --fit-window 0.37 0.60
```

主な出力の見方:

- `contact_region_summary.csv`: `ap_low/ap_high` と `interior` のめり込み・反力を比較します。`negative_gap_fraction` は候補ペアのサンプル数比です。
- `contact_pairs_detail.csv`: 接触した `j` と左右セグメント、マスクされたペア (`excluded=1`) を確認します。
- `modal_force_split.csv`: 流体・AP低端・AP高端・内部接触力と、積分へ渡した `applied_force` を比較します。荷重時刻は `t`、状態時刻は `t+dt` です。
- `gap_field.csv`: 初期、基準窓終端、擾乱直後の符号付きy間隙です。対応表面格子節点の平均 `(x,z)` を使うことを `sampling_method` に明記しています。
- `flow_diagnostics_v2.csv`: 剥離理由とフォールバック率を確認します。
- `comparison_summary.csv`: 成長率、周波数、fit誤差、別fit窓の結果を比較します。符号が窓依存なら `INCONCLUSIVE_WINDOW_DEPENDENT` です。
- `perturbation_energy*.csv`: 平衡差分仕事と収支です。仕事は出力間隔ではなく内部dtで累積し、擾乱注入ステップを除外します。SI整合を確認するまではJ/Wではなく generalized units と解釈します。

`hypothesis_report.md` はH1〜H5に必要な比較を案内します。接触除外で発振しても
正式モデルへ直ちに採用せず、形状・境界・接触探索の妥当性を確認してください。
固定節点IDファイル（0始まりVTU点ID）を指定すると
`fixed_boundary_mode_audit.csv` に各モードの固定節点最大値と全体最大値の比を
出力します。
