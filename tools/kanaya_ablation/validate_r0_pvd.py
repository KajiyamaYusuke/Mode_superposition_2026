#!/usr/bin/env pvpython
"""Headless ParaView reader check; no rendering or solver invocation."""
import json
import sys
from pathlib import Path

from paraview import servermanager
from paraview.simple import PVDReader, UpdatePipeline


root = Path(sys.argv[1]).resolve()
results = {}
for name, count, samples in [("R0_0p12_0p17.pvd", 101, [0, 1, 50, 100]),
                              ("R0_keyframes.pvd", 2, [0, 1])]:
    reader = PVDReader(FileName=str(root / name))
    times = list(reader.TimestepValues)
    assert len(times) == count and all(a < b for a, b in zip(times, times[1:])), times
    frames = []
    for index in samples:
        UpdatePipeline(time=times[index], proxy=reader)
        grid = servermanager.Fetch(reader)
        assert grid.GetNumberOfPoints() > 0 and grid.GetNumberOfCells() > 0
        data = grid.GetPointData()
        folds = data.GetArray("fold_id")
        timestamp = data.GetArray("time_s")
        assert folds and timestamp and folds.GetRange() == (0., 1.)
        assert all(abs(t - times[index]) < 1e-14 for t in timestamp.GetRange())
        assert data.GetArray("displacement_mm").GetNumberOfComponents() == 3
        frames.append(dict(time_s=times[index], points=grid.GetNumberOfPoints(), cells=grid.GetNumberOfCells()))
    results[name] = dict(passed=True, time_count=len(times), start_s=times[0], end_s=times[-1], read_frames=frames)
results["passed"] = True
with (root / "paraview_reader_validation.json").open("x") as stream:
    json.dump(results, stream, indent=2)
print("PASS ParaView PVDReader: 101 physical times and two independent keyframes")
