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

The analyzer reads `dt` from `manifest.txt` unless `--dt` is supplied. It
interprets the first column of `area.dat` as a step number, produces
`metrics.json`, `cycle_metrics.csv`, `poincare_points.csv`, `return_map.csv`,
and diagnostic figures under `figures/`. Always inspect
`figures/time_series.png` to confirm that detected events match the waveform.

The same command independently detects peaks in the `uyL_mm` and `uyR_mm`
columns of `displace.dat`. It writes:

- `left_right_ratio_metrics.csv`: representative periods, amplitudes, ratios,
  dominant frequencies, and 2:1 candidate classes.
- `left_cycle_metrics.csv` and `right_cycle_metrics.csv`: raw periods,
  validity flags, peaks, troughs, and peak-to-trough amplitudes.
- `figures/left_right_displacement.png`,
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
