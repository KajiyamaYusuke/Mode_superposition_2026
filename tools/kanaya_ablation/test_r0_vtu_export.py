"""Synthetic reconstruction/VTU/PVD tests, without a coupled model."""
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import export_r0_vtu as export


class ExportTests(unittest.TestCase):
    def test_relocated_archive_identity_is_hash_checked(self):
        entry = dict(source="/old/case/modal_native.bin", archive="saved/modal_native.bin", sha256="valid")
        records = {entry["source"]: entry}
        self.assertEqual(export.archived_reference_record(records, Path("/new/case/modal_native.bin"), "valid"), entry)
        with self.assertRaises(AssertionError):
            export.archived_reference_record(records, Path("/new/case/modal_native.bin"), "wrong")

    def setUp(self):
        export._CACHE.clear()
        export._COUNTS = [4, 4]
        export._DT = 1e-5
        export._REFERENCE = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                                      [0, 3, 0], [1, 3, 0], [0, 2, 0], [0, 3, 1]], dtype=float)
        export._MODES = [np.zeros((20, 4, 3)), np.zeros((20, 4, 3))]
        export._MODES[0][0, :, 1] = 2
        export._MODES[1][0, :, 1] = -2
        export._NATIVE = np.zeros((3, 2, 20, 7))
        export._NATIVE[1, :, 0, 4] = .001
        export._CONNECTIVITY = np.arange(8, dtype=np.int64)
        export._OFFSETS = np.array([0, 4, 8], dtype=np.int64)
        export._TYPES = np.array([10, 10], dtype=np.uint8)
        export._FOLDS = np.repeat(np.arange(2, dtype=np.uint8), 4)
        export._CELL_FOLDS = np.arange(2, dtype=np.uint8)
        export._POINT_IDS = np.tile(np.arange(4, dtype=np.int64), 2)

    def test_units_and_right_mirror(self):
        points, displacement = export.reconstruct(1)
        np.testing.assert_array_equal(displacement[:4, 1], np.full(4, 2.))
        np.testing.assert_array_equal(displacement[4:, 1], np.full(4, -2.))
        np.testing.assert_array_equal(points, export._REFERENCE + displacement)

    def test_compressed_vtu_roundtrip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="r0-vtu-unit-") as directory:
            export._OUTPUT = Path(directory)
            task = (1, str(Path(directory) / "frame.vtu"), "unit")
            row = export.write_frame(task)
            export.verify_frame(row)
            with self.assertRaises(FileExistsError):
                export.write_frame(task)

    def test_pvd_uses_physical_times(self):
        with tempfile.TemporaryDirectory(prefix="r0-pvd-unit-") as directory:
            path = Path(directory) / "series.pvd"
            export.write_pvd(path, [dict(step=50, time_s=.0005, file="frame50.vtu"),
                                    dict(step=0, time_s=0., file="frame0.vtu")])
            entries = ET.parse(path).findall("./Collection/DataSet")
            self.assertEqual([float(e.attrib["timestep"]) for e in entries], [0., .0005])

    def test_large_compressed_binary_roundtrip(self):
        # Small files missed the installed reader's appended-data failure.
        export._COUNTS = [500, 500]
        export._REFERENCE = np.random.default_rng(11).random((1000, 3))
        export._MODES = [np.zeros((20, 500, 3)), np.zeros((20, 500, 3))]
        export._CONNECTIVITY = np.arange(1000, dtype=np.int64)
        export._OFFSETS = np.arange(0, 1001, 4, dtype=np.int64)
        export._TYPES = np.full(250, 10, dtype=np.uint8)
        export._FOLDS = np.repeat(np.arange(2, dtype=np.uint8), 500)
        export._CELL_FOLDS = np.repeat(np.arange(2, dtype=np.uint8), 125)
        export._POINT_IDS = np.tile(np.arange(500, dtype=np.int64), 2)
        with tempfile.TemporaryDirectory(prefix="r0-vtu-large-unit-") as directory:
            export._OUTPUT = Path(directory)
            row = export.write_frame((1, str(Path(directory) / "large.vtu"), "unit"))
            self.assertGreater(row["bytes"], 8192)
            export.verify_frame(row)


if __name__ == "__main__":
    unittest.main()
