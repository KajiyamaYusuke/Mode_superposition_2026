# R0 residual aperture + single S3 upstream-circuit experiment

The fetched specification is
`docs/plans/KANAYA_R0_RESIDUAL_GAP_AND_S3_LEGACY_UPSTREAM_20261021.md`.
Its Git blob is `22e0ee741950d54bbe705c03116e33664c4dd5e9`, fetched from
`codex/kanaya-force-energy-implementation-20261021` in
`KajiyamaYusuke/Mode_superposition_2026`. No branch switch, merge, stash,
commit or push is performed.

Append to a **copy of the actual S2 params**, not production `input/param.txt`:

```ini
upstreamAcousticProfile = legacy_main
```

The default `modern` preserves existing upstream constants and equations.
`legacy_main` chooses inlet radius .15 m / length .50 m, subglottal length
.25 m / three sections, and `Lu=rho*dx/A`. The subglottal radius, vocal
tract, structural integration, damping, contact, ramp and all six legacy
fluid choices are inherited from S2. Do not replace the existing `2*Lu`
flow-update denominator: that is also what the original C++ uses.

The profile is validated; unknown strings are rejected. A coefficient CSV
and manifest contain the initialized `Lui,Cui,Lu,Cu,R2,La,Ca,Lr,Rr`.
`SIMULATION_ACOUSTIC_TRACE=1` adds compact per-step Uu/Pu, physical boundary
pressures and a post-step NR residual/term breakdown. It does not change
the production NR loop. NR iteration count -1 means **not instrumented**.
The original loop uses 1e-9 Pa; independent postchecking allows 1e-8 Pa
roundoff. Closed-flow rows have NR_active=0 and are not NR balance samples.

An initial inside-NR-loop instrumentation attempt changed tiny-flow last
bits (up to 2.52e-29 m3/s) and was rejected before production. The final
observer is outside the kernel, whose body is byte-identical to the
archived previous source. S2/default/explicit-modern physical outputs were
then bitwise verified. Failed/reviewed short-test artifacts are retained.

Build in a NEW directory; never rebuild the previous reference binary:

```bash
cmake -S . -B build/r0_s3_20261022 -DCMAKE_BUILD_TYPE=Release -DENABLE_KANAYA_DIAGNOSTICS=ON
cmake --build build/r0_s3_20261022 -j
ctest --test-dir build/r0_s3_20261022 --output-on-failure
```

Workflow scripts are machine-specific and collision-checked:

- `next_common.py`: prepares inventories, copied S2 fixtures, CPU selection
  and provenance. Actual data paths must be checked before reuse.
- `validate_next.py smoke`: short old/new S2, explicit-modern, open modern,
  closed/open S3, original initialization oracle and parser checks. A failed
  gate forbids production. Existing matching short results are reusable.
- `validate_next.py production`: exactly one 20,000-step S3, eight threads
  pinned to distinct currently exposed cores, no deformation VTU.
- `analyze_r0_residual_aperture.py`: existing all-mode q geometry replay,
  with the unchanged production C++ intersections compiled into an OFFLINE
  driver. It never invokes the coupled run loop. A separate old VTU checks
  the q_n versus updated q_(n+1) time convention.
- `compare_s2_s3.py`: same old 4999/4999-step windows and final-interval
  exclusion, numerical gates, work decomposition and figures. Old S2 Uu/Pu
  are reconstructed algebraically from STORED Ug, not a new NR/structural/
  contact solve; recorded Pg/Pd0 and short actual Uu/Pu validate this replay.
- `finish_next.py audit`: checks protected inputs/results, source scope,
  kernel identity, execution count, tests and final provenance.
- `finish_next.py archive`: appends the completed source/results to the
  persistent archive; never overwrites its previous snapshot.

Generated outputs are in `output/kanaya_ablation_next/`. Prior data remain in
`output/kanaya_ablation/` and must not be overwritten. The persistent archive
is `/home/kajiyama/code/kanaya_validation_archive/20261022_r0_s3/`.
Generated heavy VTUs are omitted from the archive with explicit manifests;
all-mode histories and original model inputs are retained. Archived params
retain original absolute paths; use the archive input mapping to resolve
them on another machine. Intel runtime/library paths are also machine-specific.
