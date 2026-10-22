#!/usr/bin/env python3
"""Export existing R0 q_n to compressed volume VTU/PVD; no coupled simulation."""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy


ROOT = Path(__file__).resolve().parents[2]
CASE = Path(os.environ.get("KANAYA_R0_CASE", ROOT / "output/kanaya_ablation/R0_modern")).resolve()
ARCHIVE = Path(os.environ.get("KANAYA_ARCHIVE_ROOT", ROOT / "output/kanaya_archive")).resolve()
BUILD = Path(os.environ.get("KANAYA_BUILD_DIR", ROOT / "build/diagnostic")).resolve()
SPEC = Path(os.environ.get("KANAYA_VTU_SPEC_DIR", ROOT / "output/kanaya_ablation_next/r0_vtu_export_spec")).resolve()
REFERENCE = Path(os.environ.get("KANAYA_NEXT_REFERENCE_ROOT", ROOT / "output/kanaya_ablation_next")).resolve()
DESTINATION = ARCHIVE / "visualization/R0_0p12_0p17"
_CACHE = {}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def archived_reference_record(saved, path, digest):
    # A relocated case keeps its content identity, not its old /tmp pathname.
    direct = saved.get(str(path))
    if direct is not None and direct["sha256"] == digest:
        return direct
    candidates = [row for row in saved.values()
                  if row["sha256"] == digest and Path(row["source"]).name == Path(path).name]
    if not candidates:
        raise AssertionError(f"no archived content identity for {path}")
    return candidates[0]


def dump(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def snapshot(roots):
    entries = {}
    for root in roots:
        for path in root.rglob("*"):
            if path.is_file() and ".git" not in path.parts:
                stat = path.stat()
                entries[str(path)] = [stat.st_size, stat.st_mtime_ns]
    return entries


def git_state(root):
    def call(*args):
        return subprocess.check_output(["git", "--no-optional-locks", "-C", str(root), *args], text=True)
    return dict(head=call("rev-parse", "HEAD"), branch=call("branch", "--show-current"),
                status=call("status", "--porcelain=v1", "--untracked-files=all"),
                diff_sha256=hashlib.sha256(call("diff", "--binary").encode()).hexdigest(),
                staged_sha256=hashlib.sha256(call("diff", "--cached", "--binary").encode()).hexdigest())


def reconstruct(step):
    if step in _CACHE:
        return _CACHE[step]
    displacements = []
    for side, modes in enumerate(_MODES):
        displacement = np.zeros((_COUNTS[side], 3), dtype=np.float64)
        # Same mode order and multiplication order as State::mode2uf.
        # Vectorize across nodes, not across the summation dimension.
        for m in range(20):
            displacement += modes[m] * _NATIVE[step, side, m, 4] * 1e3
        displacements.append(displacement)
    displacement = np.concatenate(displacements)
    return _REFERENCE + displacement, displacement


def named_array(values, name):
    array = numpy_to_vtk(np.ascontiguousarray(values), deep=True)
    array.SetName(name)
    return array


def make_grid(step, points, displacement):
    grid = vtk.vtkUnstructuredGrid()
    coords = vtk.vtkPoints()
    coords.SetData(named_array(points, "Points"))
    grid.SetPoints(coords)
    cells = vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(_OFFSETS, deep=True),
                  numpy_to_vtkIdTypeArray(_CONNECTIVITY, deep=True))
    grid.SetCells(named_array(_TYPES, "types"), cells)
    data = grid.GetPointData()
    data.AddArray(named_array(_FOLDS, "fold_id"))
    data.AddArray(named_array(np.full(len(points), step * _DT), "time_s"))
    data.AddArray(named_array(_POINT_IDS, "original_point_id"))
    data.SetVectors(named_array(displacement, "displacement_mm"))
    grid.GetCellData().AddArray(named_array(_CELL_FOLDS, "fold_id"))
    grid.GetFieldData().AddArray(named_array(np.array([step * _DT]), "time_s"))
    grid.GetFieldData().AddArray(named_array(np.array([step], dtype=np.int64), "step"))
    return grid


def worker_init():
    vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(1)
    vtk.vtkSMPTools.Initialize(1)


def write_frame(task):
    step, filename, kind = task
    path = Path(filename)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    started = time.perf_counter()
    points, displacement = reconstruct(step)
    if not np.isfinite(points).all() or not np.isfinite(displacement).all():
        raise ValueError("nonfinite reconstructed coordinates")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(make_grid(step, points, displacement))
    # Installed VTK fails to reread large appended files (raw AND base64).
    # Inline compressed binary passes the actual-volume reader gate.
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetCompressionLevel(1)
    if writer.Write() != 1 or writer.GetErrorCode() != 0:
        raise RuntimeError(f"VTU write failed: {path}")
    return dict(kind=kind, file=str(path.relative_to(_OUTPUT)), step=step, time_s=step * _DT,
                points=len(points), cells=len(_TYPES), bytes=path.stat().st_size,
                max_displacement_mm=float(np.linalg.norm(displacement, axis=1).max()),
                write_elapsed_s=time.perf_counter() - started, sha256=sha(path))


def parallel_write(tasks, workers):
    started = time.perf_counter()
    with mp.get_context("fork").Pool(workers, initializer=worker_init) as pool:
        rows = pool.map(write_frame, tasks, chunksize=1)
    return rows, time.perf_counter() - started


def write_pvd(path, rows):
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = ET.SubElement(root, "Collection")
    for row in sorted(rows, key=lambda row: row["step"]):
        ET.SubElement(collection, "DataSet", timestep=f"{row['time_s']:.8f}", group="", part="0", file=row["file"])
    with path.open("xb") as stream:
        ET.ElementTree(root).write(stream, encoding="utf-8", xml_declaration=True)


def verify_frame(row):
    path = _OUTPUT / row["file"]
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert not reader.GetErrorCode()
    assert grid.GetNumberOfPoints() == len(_REFERENCE) and grid.GetNumberOfCells() == len(_TYPES)
    points = vtk_to_numpy(grid.GetPoints().GetData())
    expected, displacement = reconstruct(row["step"])
    np.testing.assert_array_equal(points, expected)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetCells().GetConnectivityArray()), _CONNECTIVITY)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetCells().GetOffsetsArray()), _OFFSETS)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetCellTypesArray()), _TYPES)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetPointData().GetArray("fold_id")), _FOLDS)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetCellData().GetArray("fold_id")), _CELL_FOLDS)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetPointData().GetArray("original_point_id")), _POINT_IDS)
    np.testing.assert_array_equal(vtk_to_numpy(grid.GetPointData().GetArray("displacement_mm")), displacement)
    assert np.isfinite(points).all()
    assert grid.GetPointData().GetArray("time_s").GetRange() == (row["time_s"], row["time_s"])
    assert sha(path) == row["sha256"]


def prepare(output, inputs):
    directory = output / "preparation"
    directory.mkdir()
    work = directory / "work"
    work.mkdir()
    # Link the already-built unchanged production objects. Do not rebuild or
    # modify the existing build tree, CMake files or simulation executable.
    words = shlex.split((BUILD / "CMakeFiles/kanaya_section_replay.dir/link.txt").read_text())
    objects = [BUILD / word for word in words if word.endswith(".o") and "/src/" in word]
    includes = [ROOT / "include" / folder for folder in ["", "simulation", "structure", "fluid", "forces", "diagnostics"]]
    command = [words[0], "-O3", "-std=gnu++17", "-fiopenmp", "-DENABLE_KANAYA_DIAGNOSTICS"]
    command += [f"-I{p}" for p in includes]
    command += [str(ROOT / "tools/kanaya_ablation/r0_vtu_basis.cpp"), *map(str, objects), "-o", str(directory / "r0_vtu_basis")]
    command += [word for word in words if word.startswith("/") and word.endswith((".so", ".a"))]
    inputs.update({str(path): sha(path) for path in objects})
    dump(directory / "build_command.json", command)
    with (directory / "build.log").open("x") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="4", OMP_PROC_BIND="close", OMP_PLACES="cores",
               SIMULATION_RUN_DIR=str(directory / "initialization_output"))
    env.pop("SIMULATION_ACOUSTIC_TRACE", None)
    with (directory / "initialize.log").open("x") as log:
        subprocess.run([str(directory / "r0_vtu_basis"), str(CASE / "params.txt"),
                        str(CASE / "output/kanaya_diagnostic/modal_native.bin"), str(directory)],
                       cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    return directory


def main():
    global _OUTPUT, _REFERENCE, _MODES, _NATIVE, _COUNTS, _DT, _OFFSETS, _CONNECTIVITY, _TYPES
    global _FOLDS, _POINT_IDS, _CELL_FOLDS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DESTINATION)
    args = parser.parse_args()
    _OUTPUT = args.output.resolve()
    if _OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {_OUTPUT}")
    assert sys.byteorder == "little"
    started = time.perf_counter()
    original_roots = [Path(p).resolve() for p in os.environ.get("KANAYA_PROTECTED_ROOTS", "").split(os.pathsep) if p]
    original_roots = [p for p in original_roots if p != ROOT]
    git_before = {str(p): git_state(p) for p in [ROOT, *original_roots]}
    before = snapshot([ROOT, ARCHIVE])
    manifest = json.loads((CASE / "manifest.json").read_text())
    native_meta_path = CASE / "output/kanaya_diagnostic/run_metadata.json"
    native_meta = json.loads(native_meta_path.read_text())
    assert manifest["name"] == "R0_modern" and manifest["effective_parameters"]["glottal_method_profile"] == "modern"
    assert native_meta["nmode"] == [20, 20] and native_meta["integrated_intervals"] == 19999
    assert native_meta["modal_binary_fields"][4] == "q_n" and "step, side L/R" in native_meta["modal_binary_order"]
    source_paths = [CASE / "manifest.json", CASE / "params.txt", native_meta_path,
                    CASE / "output/kanaya_diagnostic/modal_native.bin", CASE / "output/kanaya_diagnostic/flow_geometry.csv",
                    REFERENCE / "r0_aperture/validation.json",
                    REFERENCE / "r0_aperture/selected_events.csv",
                    ARCHIVE / "archive_manifest.json", SPEC / "plan.md", SPEC / "plan_provenance.json"]
    inputs = {str(path): sha(path) for path in source_paths}
    for path, expected in manifest["model_inputs"].items():
        assert sha(path) == expected["sha256"]
        inputs[path] = expected["sha256"]
    for relative in ["src/structure/geometry.cpp", "src/structure/modedata.cpp", "src/structure/state.cpp", "src/fluid/channelSections.cpp"]:
        assert sha(ROOT / relative) == manifest["source_sha256"][relative]
        inputs[str(ROOT / relative)] = sha(ROOT / relative)
    archive_manifest = json.loads((ARCHIVE / "archive_manifest.json").read_text())
    saved = {row["source"]: row for row in archive_manifest["copied"]}
    for path in [CASE / "params.txt", native_meta_path, CASE / "output/kanaya_diagnostic/modal_native.bin", *map(Path, manifest["model_inputs"])]:
        entry = archived_reference_record(saved, path, inputs[str(path)])
        assert sha(ARCHIVE / entry["archive"]) == inputs[str(path)]
    plan = (SPEC / "plan.md").read_bytes()
    plan_source = json.loads((SPEC / "plan_provenance.json").read_text())
    assert hashlib.sha1(b"blob " + str(len(plan)).encode() + b"\0" + plan).hexdigest() == plan_source["git_blob_sha"]
    _OUTPUT.mkdir(parents=True, exist_ok=False)
    dump(_OUTPUT / "preservation_before.json", dict(files=before, git=git_before))
    shutil.copy2(SPEC / "plan.md", _OUTPUT / "plan.md")
    shutil.copy2(SPEC / "plan_provenance.json", _OUTPUT / "plan_provenance.json")
    (_OUTPUT / "tools").mkdir()
    for name in ["export_r0_vtu.py", "r0_vtu_basis.cpp", "validate_r0_pvd.py", "test_r0_vtu_export.py"]:
        path = ROOT / "tools/kanaya_ablation" / name
        inputs[str(path)] = sha(path)
        shutil.copy2(path, _OUTPUT / "tools" / name)
    with (_OUTPUT / "unit_tests.log").open("x") as log:
        subprocess.run([sys.executable, "-B", str(ROOT / "tools/kanaya_ablation/test_r0_vtu_export.py")],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    print("Inputs and archived R0 identity verified; preparing exact production-normalized basis", flush=True)
    preparation_started = time.perf_counter()
    preparation = prepare(_OUTPUT, inputs)
    preparation_elapsed = time.perf_counter() - preparation_started
    basis = json.loads((preparation / "basis.json").read_text())
    _COUNTS = basis["points_per_fold"]
    _DT = basis["dt_s"]
    assert _DT == 1e-5 == native_meta["dt_s"]
    _REFERENCE = np.fromfile(preparation / "reference.bin", dtype="<f8").reshape(-1, 3)
    all_modes = np.fromfile(preparation / "modes.bin", dtype="<f8")
    _MODES = []
    start = 0
    for count in _COUNTS:
        stop = start + 20 * count * 3
        _MODES.append(all_modes[start:stop].reshape(20, count, 3))
        start = stop
    assert start == len(all_modes)
    _OFFSETS = np.r_[np.int64(0), np.fromfile(preparation / "offsets.bin", dtype="<i8")]
    _CONNECTIVITY = np.fromfile(preparation / "connectivity.bin", dtype="<i8")
    _TYPES = np.fromfile(preparation / "types.bin", dtype=np.uint8)
    _FOLDS = np.repeat(np.arange(2, dtype=np.uint8), _COUNTS)
    _CELL_FOLDS = np.repeat(np.arange(2, dtype=np.uint8), basis["cells_per_fold"])
    _POINT_IDS = np.concatenate([np.arange(n, dtype=np.int64) for n in _COUNTS])
    native_path = CASE / "output/kanaya_diagnostic/modal_native.bin"
    assert native_path.stat().st_size == 19999 * 2 * 20 * 7 * 8
    _NATIVE = np.memmap(native_path, dtype="<f8", mode="r", shape=(19999, 2, 20, 7))
    np.testing.assert_allclose(_NATIVE[1:, :, :, 4] - _NATIVE[:-1, :, :, 4], _NATIVE[:-1, :, :, 6], rtol=1e-12, atol=1e-22)
    steps = list(range(12000, 17001, 50))
    keys = [12480, 12792]
    flow = pd.read_csv(CASE / "output/kanaya_diagnostic/flow_geometry.csv", float_precision="round_trip").set_index("step")
    for step in steps + keys:
        assert step in flow.index and abs(flow.loc[step, "time_s"] - step * _DT) < 1e-14
    control_errors = {}
    for step in [12000, 12480, 12792, 15001, 17000]:
        points, _ = reconstruct(step)
        control = np.fromfile(preparation / f"control_{step}.bin", dtype="<f8").reshape(-1, 3)
        np.testing.assert_allclose(points, control, rtol=0, atol=1e-11)
        control_errors[str(step)] = float(np.max(abs(points - control)))
    pristine = Path(os.environ.get("KANAYA_PRISTINE_VTU", CASE.parent / "validation/pristine_modern_full/output/result/deform_combined0000.vtu"))
    inputs[str(pristine)] = sha(pristine)
    reference_reader = vtk.vtkXMLUnstructuredGridReader()
    reference_reader.SetFileName(str(pristine))
    reference_reader.Update()
    reference_grid = reference_reader.GetOutput()
    reference_points = vtk_to_numpy(reference_grid.GetPoints().GetData())
    points, _ = reconstruct(15001)
    np.testing.assert_allclose(points, reference_points, rtol=0, atol=1e-11)
    pristine_error = float(np.max(abs(points - reference_points)))
    np.testing.assert_array_equal(vtk_to_numpy(reference_grid.GetCells().GetConnectivityArray()), _CONNECTIVITY)
    np.testing.assert_array_equal(vtk_to_numpy(reference_grid.GetCells().GetOffsetsArray()), _OFFSETS)
    np.testing.assert_array_equal(vtk_to_numpy(reference_grid.GetCellTypesArray()), _TYPES)
    old_control = REFERENCE / "r0_aperture/control_coordinates.csv"
    inputs[str(old_control)] = sha(old_control)
    old_points = pd.read_csv(old_control, float_precision="round_trip")[["x_mm", "y_mm", "z_mm"]].to_numpy()
    np.testing.assert_allclose(points, old_points, rtol=0, atol=1e-11)
    checks = pd.read_csv(preparation / "geometry_checks.csv", float_precision="round_trip")
    np.testing.assert_allclose(checks.Amin_mm2, flow.loc[checks.step, "Amin_mm2"], rtol=0, atol=1e-10)
    np.testing.assert_allclose(checks.x_min_mm, flow.loc[checks.step, "xMin_mm"], rtol=0, atol=1e-10)
    minimum = checks[checks.step == 12792].iloc[0]
    assert minimum.min_index == 32 and abs(minimum.Amin_mm2 - .075628082917) < 1e-10
    print(f"Geometry gates passed; {len(_REFERENCE)} points, {len(_TYPES)} cells; benchmarking 1/2/4 writers", flush=True)
    probe = write_frame((12792, str(preparation / "compressed_format_probe.vtu"), "format_probe"))
    verify_frame(probe)
    for step in steps[::14][:8]:
        _CACHE[step] = reconstruct(step)
    bench = []
    for repeat, counts in enumerate([[1, 2, 4], [4, 2, 1]], 1):
        for workers in counts:
            target = _OUTPUT / "benchmarks" / f"repeat{repeat}_workers{workers}"
            target.mkdir(parents=True)
            tasks = [(step, str(target / f"deform_combined_n{step:06d}.vtu"), "benchmark") for step in sorted(_CACHE)]
            rows, elapsed = parallel_write(tasks, workers)
            bench.append(dict(repeat=repeat, workers=workers, frames=len(rows), elapsed_s=elapsed,
                              bytes=sum(row["bytes"] for row in rows)))
            print(f"writers={workers}, repeat={repeat}: {elapsed:.3f} s for {len(rows)} frames", flush=True)
    benchmark = pd.DataFrame(bench)
    benchmark.to_csv(_OUTPUT / "writer_benchmark.csv", index=False)
    selected = int(benchmark.groupby("workers").elapsed_s.mean().idxmin())
    projected_bytes = benchmark.bytes.mean() / len(_CACHE) * 103
    if projected_bytes > 2 * 2**30:
        raise RuntimeError("projected volume series exceeds 2 GiB; review surface-only option without overwriting this attempt")
    for name in ["frames", "keyframes"]:
        (_OUTPUT / name).mkdir()
    tasks = [(step, str(_OUTPUT / "frames" / f"deform_combined_n{step:06d}.vtu"), "frame") for step in steps]
    tasks += [(step, str(_OUTPUT / "keyframes" / f"deform_combined_n{step:06d}.vtu"), "keyframe") for step in keys]
    rows, export_elapsed = parallel_write(tasks, selected)
    frame_rows = [row for row in rows if row["kind"] == "frame"]
    key_rows = [row for row in rows if row["kind"] == "keyframe"]
    write_pvd(_OUTPUT / "R0_0p12_0p17.pvd", frame_rows)
    write_pvd(_OUTPUT / "R0_keyframes.pvd", key_rows)
    pd.DataFrame(rows).to_csv(_OUTPUT / "frame_index.csv", index=False)
    print(f"Exported 101 volume frames + 2 keyframes with {selected} writers; validating all files", flush=True)
    validation_started = time.perf_counter()
    for row in rows:
        verify_frame(row)
    for filename, expected in [("R0_0p12_0p17.pvd", steps), ("R0_keyframes.pvd", keys)]:
        entries = ET.parse(_OUTPUT / filename).findall("./Collection/DataSet")
        np.testing.assert_allclose([float(row.attrib["timestep"]) for row in entries], np.array(expected) * _DT, rtol=0, atol=1e-14)
        assert all((_OUTPUT / row.attrib["file"]).is_file() for row in entries)
    with (_OUTPUT / "paraview_reader.log").open("x") as log:
        subprocess.run(["pvpython", str(ROOT / "tools/kanaya_ablation/validate_r0_pvd.py"), str(_OUTPUT)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    validation_elapsed = time.perf_counter() - validation_started
    for path, digest in inputs.items():
        assert sha(path) == digest, f"input changed: {path}"
    for path, stat in before.items():
        p = Path(path)
        assert p.is_file() and [p.stat().st_size, p.stat().st_mtime_ns] == stat, f"existing file changed: {p}"
    assert git_before == {str(p): git_state(p) for p in [ROOT, *original_roots]}
    total_bytes = sum(row["bytes"] for row in rows)
    metadata = dict(passed=True, case="R0_modern", output=str(_OUTPUT), coupled_steps_executed=0,
                    mode_count_per_fold=[20, 20], native_shape=[19999, 2, 20, 7], native_coordinate_field="q_n",
                    state_timing="q_n at n*dt, before step-n integration; pristine step15000 output corresponds to q_15001",
                    interval_s=[.120, .170], frame_interval_s=.0005, frame_count=len(frame_rows), keyframe_count=len(key_rows),
                    key_steps=keys, points_per_fold=_COUNTS, cells_per_fold=basis["cells_per_fold"],
                    geometry="complete bilateral volume, original cell connectivity; no surface extraction or decimation",
                    point_units="mm", displacement_scale=1., compression="VTK inline base64 binary, zlib level1, Float64",
                    format_probe=probe, failed_appended_attempt=str(_OUTPUT.with_name(_OUTPUT.name + "_incomplete_appended_readback_failure")),
                    inputs_sha256=inputs, plan=plan_source, cpp_control_max_abs_error_mm=control_errors,
                    pristine_q15001_max_abs_error_mm=pristine_error,
                    prior_aperture_control_max_abs_error_mm=float(np.max(abs(points - old_points))),
                    minimum_keyframe_section=dict(step=12792, section_index=int(minimum.min_index),
                                                  x_mm=float(minimum.x_min_mm), Amin_mm2=float(minimum.Amin_mm2)),
                    all_103_readbacks_exact=True, topology_matches_pristine=True, finite_coordinates=True,
                    paraview_headless_reader_passed=True, previous_files_unchanged=len(before), previous_git_states_unchanged=True,
                    selected_writer_processes=selected, benchmark_selection="mean of two 8-frame trials, forward/reverse worker order",
                    worker_VTK_and_BLAS_threads=1, volume_VTU_bytes=total_bytes,
                    benchmark_VTU_bytes=int(benchmark.bytes.sum()), preparation_elapsed_s=preparation_elapsed,
                    export_elapsed_s=export_elapsed, validation_elapsed_s=validation_elapsed,
                    total_elapsed_s=time.perf_counter() - started,
                    parent_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    max_single_child_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    memory_note="RSS peaks are per process, not summed simultaneous memory; fork shares read-only modal basis",
                    vtk_version=vtk.vtkVersion.GetVTKVersion(), numpy_version=np.__version__)
    dump(_OUTPUT / "metadata.json", metadata)
    readme = f"""# R0 0.120〜0.170 s：ParaView用VTUアニメーション

保存先：`{_OUTPUT}`

## 開く方法

1. ParaViewの File → Open で **`R0_0p12_0p17.pvd`** を選び、Apply。
2. Representation を Surface または Surface With Edges にする。Reset Cameraで全体を表示する。
3. Coloring を `fold_id`（0=left、1=right）にすると左右を区別できる。別々に表示したい場合はThresholdを2つ作り、Cell Dataのfold_idをそれぞれ0だけ／1だけ選ぶ。
4. Playで101時刻を順に再生する。必要なら View → Animation View で時間範囲0.12〜0.17 s、モードSnap To TimeStepsを選ぶ。再生の壁時計速度はGUIの設定であり物理時間ではない。
5. 最大開口付近0.12480 sと最小開口付近0.12792 sの実ステップは、別の **`R0_keyframes.pvd`** で開く。これらは等間隔の主時系列に混ぜていない。

Pointsは**変形済み絶対座標[mm]・変形倍率1倍**。Warp By Vectorを重ねると二重変形になるため不要。`displacement_mm`は初期gap調整後の参照形状からの物理変位で、点順序は左右各元メッシュのまま。点の`original_point_id`は各fold内の0始まり元節点ID。点とセルの`fold_id`、点とFieldDataの`time_s`、FieldDataの`step`を付加した。

## 内容と正確性

- R0_modernの保存済みmodal_native.binから左右全20モードのq_nだけを使用。連成計算・流量更新・接触更新・時間積分は0回。101枚の間隔は50step=0.00050 s、末端0.17000 sも実在。
- **体積メッシュ全体を101枚＋キー2枚**保存。表面のみの簡略版ではなく、元のセル・接続関係を保持。計{len(_REFERENCE):,}点、{len(_TYPES):,}セル／枚。圧縮バイナリVTU合計{total_bytes / 2**20:.2f} MiB。表示をSurfaceにしても内部体積セルはファイルに残る。
- 元のC++読込・モード順・質量正規化・右鏡映・初期gap処理を変更せず、初期化を各foldにつき一度だけ利用。NumPyは節点方向だけベクトル化し、20モードの加算順を保持した。5状態を変更していないState::mode2ufと照合した（最大差{max(control_errors.values()):.4e} mm）。
- 既存の正規出力（n=15000更新後）とq_15001の全節点座標差は最大{pristine_error:.4e} mm。旧R0空間再生の全座標とも照合。n=12792の同じ50断面計算はAmin={minimum.Amin_mm2:.12g} mm²、x={minimum.x_min_mm:.12g} mmで一致。3D体積メッシュそのものと最狭断面積は別量である。
- 全103枚をVTKで再読込し、点・接続・セル種類・fold・時刻・変位・有限性を検査。ParaView PVDReaderでも主時系列101時刻とキー2時刻、および代表フレームの実読込に合格した。
- 元コード・入力・R0/S2/S3出力・以前のarchiveを変更しない。検査対象の既存{len(before):,}ファイルはサイズ・mtime不変、使用入力はSHA256不変、4worktreeのGit状態不変。新規出力だけを保存した。

## 速度と記録

書込プロセス1/2/4を同じ8枚で2回ずつ比較し、平均壁時計が最速の**{selected}プロセス**を採用。各プロセス内のVTK/BLASは1スレッド。`writer_benchmark.csv`に全実測を記録。比較試験の計48 VTUは`benchmarks/`に別保存（{int(benchmark.bytes.sum()) / 2**20:.2f} MiB）し、PVDの再生対象に含めない。

初期化・Releaseツール準備{preparation_elapsed:.3f} s、103枚出力{export_elapsed:.3f} s、再読込検証{validation_elapsed:.3f} s、後処理全体{metadata['total_elapsed_s']:.3f} s（metadata書込直前まで）。メモリは`metadata.json`の親／単一子プロセスRSSを参照。巨大な元VTUモード場を各フレームで再読込していない。

`frame_index.csv`：全103枚のファイル・時刻・step・点数・セル数・サイズ・最大変位・SHA256。
`metadata.json`：入力SHA・時間対応・検証誤差・使用方式。
`preparation/`：正規化済みモード、参照形状、元接続、C++照合データ・ログ。
`tools/`：今回の専用後処理ツール。`plan.md`：取得した指示原文。

書式互換性：この環境では大きなappended形式（raw/base64とも）のVTU読込が失敗したため、完成版は圧縮を維持した**inline binary（base64）**を採用。完成版は実体積メッシュ1枚の先行再読込と全103枚の再読込を通した。最初のappended試行は隣の`R0_0p12_0p17_incomplete_appended_readback_failure/`に未完試行として保持し、完成品・物理比較・再生対象に含めない。失敗は書式読込であり連成計算はどちらの試行でも0回。

再実行は元worktreeで `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B tools/kanaya_ablation/export_r0_vtu.py --output <未使用の新規パス>`。既存出力先は拒否する。元buildのReleaseオブジェクト、記載した入力絶対パス、Intelランタイム・VTK/Python・pvpythonが必要。全dt出力や再連成計算は実装していない。
"""
    with (_OUTPUT / "README.md").open("x") as stream:
        stream.write(readme)
    print(f"PASS: {_OUTPUT}; volume VTU {total_bytes / 2**20:.2f} MiB; total {metadata['total_elapsed_s']:.3f} s", flush=True)


if __name__ == "__main__":
    main()
