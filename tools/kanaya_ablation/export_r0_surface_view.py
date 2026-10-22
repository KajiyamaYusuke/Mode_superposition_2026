#!/usr/bin/env python3
"""Lightweight boundary display derived from verified R0 volume VTU files."""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from export_r0_vtu import DESTINATION, dump, sha, write_pvd, worker_init


def read_volume(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfPoints() == 90666 and grid.GetNumberOfCells() == 496536
    return grid


def convert(task):
    source_row, output_file, kind = task
    path = Path(output_file)
    if path.exists():
        raise FileExistsError(path)
    source = _ROOT / source_row['file']
    assert sha(source) == source_row['sha256']
    volume = read_volume(source)
    points = vtk_to_numpy(volume.GetPoints().GetData())[_SOURCE_IDS].copy()
    surface = vtk.vtkPolyData()
    surface.DeepCopy(_SURFACE)
    coords = vtk.vtkPoints()
    coords.SetData(numpy_to_vtk(points, deep=True))
    surface.SetPoints(coords)
    data = surface.GetPointData()
    data.Initialize()
    for name in ['fold_id', 'time_s', 'original_point_id', 'displacement_mm']:
        values = vtk_to_numpy(volume.GetPointData().GetArray(name))[_SOURCE_IDS].copy()
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        data.AddArray(array)
    ids = numpy_to_vtk(_SOURCE_IDS, deep=True)
    ids.SetName('source_combined_point_id')
    data.AddArray(ids)
    data.SetActiveVectors('displacement_mm')
    surface.GetFieldData().DeepCopy(volume.GetFieldData())
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(surface)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.SetCompressionLevel(1)
    assert writer.Write() == 1
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    replay = reader.GetOutput()
    assert replay.GetNumberOfCells() == _SURFACE.GetNumberOfCells()
    np.testing.assert_array_equal(vtk_to_numpy(replay.GetPoints().GetData()), points)
    np.testing.assert_array_equal(vtk_to_numpy(replay.GetPolys().GetConnectivityArray()),
                                  vtk_to_numpy(_SURFACE.GetPolys().GetConnectivityArray()))
    for name in ['fold_id', 'time_s', 'original_point_id', 'displacement_mm', 'source_combined_point_id']:
        np.testing.assert_array_equal(vtk_to_numpy(replay.GetPointData().GetArray(name)),
                                      vtk_to_numpy(surface.GetPointData().GetArray(name)))
    assert np.isfinite(points).all() and replay.GetPointData().GetArray('fold_id').GetRange() == (0., 1.)
    return dict(kind=kind, file=str(path.relative_to(_OUTPUT)), step=source_row['step'], time_s=source_row['time_s'],
                points=len(points), cells=surface.GetNumberOfCells(), bytes=path.stat().st_size, sha256=sha(path),
                max_displacement_mm=float(np.linalg.norm(vtk_to_numpy(data.GetArray('displacement_mm')), axis=1).max()),
                source_file=source_row['file'], source_sha256=source_row['sha256'])


def run(tasks, workers):
    started = time.perf_counter()
    with mp.get_context('fork').Pool(workers, initializer=worker_init) as pool:
        rows = pool.map(convert, tasks, chunksize=1)
    return rows, time.perf_counter() - started


def main():
    global _ROOT, _OUTPUT, _SURFACE, _SOURCE_IDS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--volume-root', type=Path, default=DESTINATION)
    args = parser.parse_args()
    _ROOT = args.volume_root.resolve()
    _OUTPUT = _ROOT / 'surface'
    if _OUTPUT.exists():
        raise FileExistsError(_OUTPUT)
    started = time.perf_counter()
    metadata = json.loads((_ROOT / 'metadata.json').read_text())
    assert metadata['passed'] and metadata['all_103_readbacks_exact'] and metadata['coupled_steps_executed'] == 0
    index = pd.read_csv(_ROOT / 'frame_index.csv', float_precision='round_trip')
    source_rows = index[index.kind == 'frame'].to_dict('records')
    assert len(source_rows) == 101
    source_stat = {str(_ROOT / row['file']): [(_ROOT / row['file']).stat().st_size,
                                           (_ROOT / row['file']).stat().st_mtime_ns] for row in source_rows}
    volume = read_volume(_ROOT / source_rows[0]['file'])
    extract = vtk.vtkDataSetSurfaceFilter()
    extract.SetInputData(volume)
    extract.PassThroughPointIdsOn()
    extract.Update()
    _SURFACE = vtk.vtkPolyData()
    _SURFACE.DeepCopy(extract.GetOutput())
    _SOURCE_IDS = vtk_to_numpy(_SURFACE.GetPointData().GetArray('vtkOriginalPointIds')).copy()
    assert len(np.unique(_SOURCE_IDS)) == len(_SOURCE_IDS)
    _OUTPUT.mkdir()
    benchmarks = []
    for workers in [1, 2, 4]:
        folder = _OUTPUT / f'benchmark_workers{workers}'
        folder.mkdir()
        tasks = [(r, str(folder / f"surface_n{r['step']:06d}.vtp"), 'benchmark') for r in source_rows[::14][:8]]
        rows, elapsed = run(tasks, workers)
        benchmarks.append(dict(workers=workers, elapsed_s=elapsed, frames=len(rows), bytes=sum(r['bytes'] for r in rows)))
        print(f'surface workers={workers}: {elapsed:.3f} s', flush=True)
    benchmark = pd.DataFrame(benchmarks)
    benchmark.to_csv(_OUTPUT / 'writer_benchmark.csv', index=False)
    selected = int(benchmark.loc[benchmark.elapsed_s.idxmin(), 'workers'])
    (_OUTPUT / 'frames').mkdir()
    rows, elapsed = run([(r, str(_OUTPUT / 'frames' / f"surface_n{r['step']:06d}.vtp"), 'frame') for r in source_rows], selected)
    write_pvd(_OUTPUT / 'R0_surface_0p12_0p17.pvd', rows)
    pd.DataFrame(rows).to_csv(_OUTPUT / 'frame_index.csv', index=False)
    pv_code = '''import json,sys
from pathlib import Path
from paraview import servermanager
from paraview.simple import PVDReader,UpdatePipeline
p=Path(sys.argv[1]);r=PVDReader(FileName=str(p/'R0_surface_0p12_0p17.pvd'));t=list(r.TimestepValues)
assert len(t)==101 and abs(t[0]-.12)<1e-14 and abs(t[-1]-.17)<1e-14
for a,b in zip(t,t[1:]): assert abs(b-a-.0005)<1e-14
for i in [0,50,100]:
 UpdatePipeline(time=t[i],proxy=r);g=servermanager.Fetch(r)
 assert g.GetNumberOfPoints()==15358 and g.GetNumberOfCells()==30708
 assert g.GetPointData().GetArray('fold_id').GetRange()==(0.,1.)
 assert all(abs(v-t[i])<1e-14 for v in g.GetPointData().GetArray('time_s').GetRange())
with (p/'paraview_reader_validation.json').open('x') as f:json.dump(dict(passed=True,time_count=len(t),points=15358,cells=30708),f)
print('PASS ParaView lightweight surface PVD reader')
'''
    with (_OUTPUT / 'paraview_reader.log').open('x') as log:
        subprocess.run(['pvpython', '-c', pv_code, str(_OUTPUT)], stdout=log, stderr=subprocess.STDOUT, check=True)
    for path, stat in source_stat.items():
        p = Path(path)
        assert [p.stat().st_size, p.stat().st_mtime_ns] == stat
    total_bytes = sum(row['bytes'] for row in rows)
    dump(_OUTPUT / 'metadata.json', dict(passed=True, coupled_steps_executed=0, frames=101,
        points_per_frame=_SURFACE.GetNumberOfPoints(), cells_per_frame=_SURFACE.GetNumberOfCells(),
        source='verified full volume VTU files; identical displacement, no resynthesis',
        point_mapping='source_combined_point_id maps surface points to the complete bilateral volume',
        omitted='interior tetrahedra and interior nodes; boundary triangle connectivity extracted once without decimation',
        bytes=total_bytes, size_fraction_of_volume=total_bytes / sum(r['bytes'] for r in source_rows),
        workers=selected, conversion_elapsed_s=elapsed, total_elapsed_s=time.perf_counter() - started,
        all_frame_readbacks_exact=True, paraview_reader_passed=True, source_files_unchanged=True,
        script_sha256=sha(Path(__file__)), volume_metadata_sha256=sha(_ROOT / 'metadata.json')))
    script = Path(__file__).read_bytes()
    with (_OUTPUT / 'export_r0_surface_view.py').open('xb') as stream:
        stream.write(script)
    text = f'''# 再生向けの軽量表面版（推奨）

ParaViewで **`{_OUTPUT / 'R0_surface_0p12_0p17.pvd'}`** を開いてApply → SurfaceまたはSurface With Edges → Play。
0.12000〜0.17000 s、0.00050 s間隔の101枚。`fold_id`は0=left、1=right。全15,358点・30,708三角形／枚、合計{total_bytes/2**20:.2f} MiB。変形済み絶対座標[mm]、倍率1倍でWarp By Vectorは不要。

検証済みの全体積VTUから**外部境界だけ**を抽出した表示用VTP。内部四面体・内部節点は省略するが、境界を間引き・平滑化せず、各時刻の同じ点の座標・変位をそのまま使用。点順序は表面部分集合へ変わるため`source_combined_point_id`で元の左右結合節点に対応付ける。体積セルの断面観察・内部場には親フォルダの全体積VTU/PVDを使う。最大／最小キー2枚は親の`R0_keyframes.pvd`で開く。

全101枚で元volumeの対応点と厳密一致・接続不変・有限性・左右ID・物理時刻を確認し、ParaViewヘッドレスPVDReaderにも合格。連成計算は0回。書込プロセス1/2/4の8枚実測から{selected}を採用、変換{elapsed:.3f} s。`frame_index.csv`に全ファイル・時刻・SHAを記録。`metadata.json`に省略内容と検証を記録。比較試行は別のbenchmarkフォルダに保持しPVD再生に含めない。
'''
    with (_OUTPUT / 'README.md').open('x') as stream:
        stream.write(text)
    # This README belongs to this new export, not to any original simulation.
    with (_ROOT / 'README.md').open('a') as stream:
        stream.write(f'\n## 軽量な表面版（アニメーション再生にはこちらを推奨）\n\n[表面版の説明](surface/README.md) と `surface/R0_surface_0p12_0p17.pvd` を追加。全101枚で{total_bytes/2**20:.2f} MiB。全体積VTUと同じ境界変位を保持し、内部四面体だけを省略した。\n')
    print(f'PASS lightweight surface series: {total_bytes/2**20:.2f} MiB, {len(rows)} frames', flush=True)


if __name__ == '__main__':
    main()
