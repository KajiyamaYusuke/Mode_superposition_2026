# Repository Guidelines

## Project Structure & Module Organization

This is a C++17 vocal-fold mode-superposition simulation.

- `src/main.cpp`: executable entry point.
- `src/simulation/`: simulation orchestration and input parameters.
- `src/structure/`: geometry, modal data, structural state, and integration. The existing `rungeStep.cpp` filename is retained; it currently includes Newmark integration.
- `src/fluid/`: channel sections, flow, and separation selection.
- `src/forces/`: fluid/contact loads, surface loads, and modal force projection.
- `include/`: matching headers and shared data structures in the same four categories: `simulation/`, `structure/`, `fluid/`, and `forces/`.
- `docs/cleanup/`: cleanup instructions, historical validation reports, and directory reorganization reports.
- `docs/plans/`: research and implementation specifications and plans.
- `input/`: parameter files, modal VTU data, frequencies, and surface meshes required at runtime.
- `tools/`: Python scripts for waveform, pressure, displacement, and spectrum analysis.
- `tests/`: existing C++ test sources; location unchanged.
- `output/`: latest scalar, CSV, and WAV results; runs overwrite this workspace.
- `result/`: generated VTU deformation output.
- `nonperiodic_vibration_analysis_plan.md`: existing reference to dynamics analysis and pressure-sweep requirements; the file is not present in this working tree and has not been replaced by another document.

Keep generated files out of commits. `build/`, `output/`, `result/`, and common analysis artifacts are ignored.

## Build, Test, and Development Commands

The build expects Intel `icpx` and OpenMP:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cd build && ./simulation ../input/param.txt
```

Enable detailed contact profiling with:

```bash
cmake -S . -B build -DPROFILE_CALCDIS=ON
```

Set reproducible OpenMP placement when benchmarking, for example:

```bash
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores ./simulation ../input/param.txt
```

Run from `build/`; several model paths are resolved relative to the working directory.

## Coding Style & Naming Conventions

Use four-space indentation and C++17 features. Follow existing naming: classes in `PascalCase`, methods and variables in `camelCase`, and files such as `forceCalculator.cpp` paired with `ForceCalculator.h`. Keep OpenMP loops deterministic with `schedule(static)`. Avoid `-ffast-math` or algebraic rewrites that can change floating-point behavior, contact branching, or physical equations.

## Testing Guidelines

There is currently no automated test framework. Every change must at least compile in Release mode. For numerical changes, compare baseline and modified runs for displacement, flow, pressure, contact force, and VTU coordinates. Record timing separately from correctness checks. Use shortened parameter files for smoke tests, but never commit reduced production settings.

## Commit & Pull Request Guidelines

History uses short descriptive commits but has no strict convention. Prefer concise imperative subjects, for example `Parallelize modal force projection`. Keep commits focused. Pull requests should describe the affected model stage, list build and validation commands, report numerical tolerances and timing changes, and explicitly note whether outputs or equations changed. Link relevant issues and include plots only when they clarify numerical behavior.
