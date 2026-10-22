# Kanaya legacy-fluid ablation

The default is still `modern`. Append named options to a **copied** parameter
file; do not replace production input files.

```ini
glottalMethodProfile = legacy_fluid_all
```

This switches material-row area, minimum-area flow coefficients, wall
viscosity, pressure-recovery stopping, downstream wall pressure, and nodal
pressure loading. It does **not** switch contact, Newmark, damping, pressure
ramp, inlet/subglottal acoustics, mass or modes. It is not a complete replay of
the historical coupled solver.

For a selective experiment use `custom` and explicitly list the six flags:

```ini
glottalMethodProfile = custom
legacyArea = 1
legacyFlowCoefficients = 0
legacyWallViscosity = 1
legacySeparation = 1
legacyPostSeparation = 1
legacySurfacePressureLoad = 1
```

`legacySurfacePressureLoad=1` requires `legacyArea=1`. Changing back to modern
area therefore requires changing mapping too. Legacy material rows are
restricted to matching topology and mirror-symmetric geometry; asymmetric
motion raises an error rather than silently applying a symmetric model.
`legacyContact=1` and `legacyIntegrator=1` are explicitly unsupported and
rejected. Contradictory flags on `modern`/`legacy_fluid_all` are also rejected.

The user-confirmed reference is `KajiyamaYusuke/Mode_superposition`, commit
`58cd92281aa6a36c732baa680d7817acbc98ddb2`. The harness compiles original
method bodies extracted verbatim from archived sources. Frozen-flow tests
check area, pressure, all three nodal force components, modal force/work and
Newmark. Source snapshots, hashes and limitations accompany the results.
The public Fortran is a secondary reference, not a substitute for this C++.

Additional optional controls, unchanged by default:

```ini
writeDeformationVtu = 0  # suppress only deformation VTU; numeric invariance tested
vocalTractLengthM = 0.175  # explicit separate acoustic experiment, NOT primary R0/R1
```

Build with Intel icpx / Release / OpenMP and diagnostics enabled:

```bash
cmake -S . -B build/switchable -DCMAKE_BUILD_TYPE=Release -DENABLE_KANAYA_DIAGNOSTICS=ON
cmake --build build/switchable -j
ctest --test-dir build/switchable --output-on-failure
```

`run_matrix.py` uses absolute model paths and independent working/output
directories; existing cases are never overwritten. It is tailored to this
machine's compiler and the protected input template. Inspect those paths
before reuse. It records parameters, source/input/binary/output hashes,
effective flags, acoustic conditions, affinity and resource usage. The
pristine baseline binary must be built **before** implementation changes;
rebuilding it afterwards would invalidate that reference.

The complete workflow is intentionally gated and staged:

1. `run_matrix.py prepare`, then `run_matrix.py baseline-build`, before edits.
2. Implement; `run_matrix.py build`; `validate.py fixed`, `parser`, `modern`,
   `open-modern` (the regression cases must first have been generated).
3. `benchmark.py`, then `benchmark_extended.py`: same 2000-step conditions,
   single-job and concurrent-job measurements, nonoversubscribed affinity,
   trajectory hash checks.
4. `production.py primary`; `compare_runs.py R0_modern R1_legacy_fluid_all --tag primary`.
5. Inspect primary results. Only then choose `production.py weak` (two separate
   damping/acoustic comparisons) or `production.py fluid` (R2–R6). Do not run
   both branches blindly. The conditional decision belongs in the report.
   `production.py targeted` selects just R3/R6 when both primary runs sustain
   open oscillations and a small mechanism comparison is more informative
   than the entire R2–R6 matrix. Failed directories are retained; retries
   must use a new case name.
6. `production.py baseline` optionally makes a full-length pristine modern
   reference. `run_matrix.py audit` checks the protected worktrees at the end.

Analysis excludes the known final geometry-skipped interval. Equal 4999-step
windows are compared: steps 10000–14998 and 15000–19998. Closure uses
`Amin <= 0.01 mm²`; CQ and closure-derived f0 are withheld unless repeated
closure has at least three periods and period CV <= 0.15. FFT peaks alone do
not establish f0 or self-excitation. Work is recomputed from native binary
data as `Q_n * (q_(n+1)-q_n)`, split exactly into mean and fluctuation work.
Mirror, sign-split, summary consistency and excluded-interval gates must pass.
Structural damping work is a modal reference, not a complete mechanical
energy balance for the nonorthogonal modal basis and numerical integrator.

For hybrid physical-x mapping with legacy abrupt zero downstream pressure,
coordinates within 32 machine epsilons of the pressure endpoint use the
upstream endpoint pressure. This is a floating-point equality tie rule, not
a physical smoothing length. A one-ulp endpoint test and a point 1e-10 mm
downstream guard the distinction. The failed pre-fix R6 attempt is excluded
from physical conclusions. `validate.py boundary` rechecks the user-source
oracle and modern/all-legacy coupled invariance after this hybrid-only fix.

Machine-specific results are in `output/kanaya_ablation/comparison_report.md`;
generated runs and figures are ignored by Git.
