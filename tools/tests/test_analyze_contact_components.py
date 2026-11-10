import importlib.util
from pathlib import Path
import unittest

import numpy as np


MODULE_PATH = Path(__file__).parents[1] / "amplitude_asymmetry" / "analyze_contact_components.py"
SPEC = importlib.util.spec_from_file_location("analyze_contact_components", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class AnalyzeContactComponentsTest(unittest.TestCase):
    def test_interval_mask_uses_left_end_and_excludes_right_end(self):
        actual = np.flatnonzero(MODULE.interval_mask(6, 0.1, 0.2, 0.5))
        np.testing.assert_array_equal(actual, [2, 3, 4])

    def test_work_decomposition_and_sign_sums(self):
        count, nmode = 6, 2
        modal = np.zeros((count, 2, nmode, 7))
        component = np.zeros((count, 2, nmode, 2))
        modal[..., 6] = 1.0
        component[..., 0] = 2.0
        component[..., 1] = -0.5
        modal[..., 3] = component.sum(axis=-1)
        frequencies = [np.array([100.0, 200.0]), np.array([90.0, 180.0])]
        old_windows = MODULE.WINDOWS
        MODULE.WINDOWS = ((0.0, 0.5),)
        try:
            rows, checks = MODULE.summarize_work(modal, component, 0.1, frequencies)
        finally:
            MODULE.WINDOWS = old_windows
        left = next(row for row in rows if row["side"] == "L" and row["mode"] == "ALL")
        self.assertEqual(left["interval_count"], 5)
        self.assertEqual(left["W_spring_J"], 20.0)
        self.assertEqual(left["W_damp_J"], -5.0)
        self.assertEqual(left["W_total_J"], 15.0)
        self.assertEqual(left["W_reference_J"], 15.0)
        self.assertEqual(left["W_spring_positive_J"], 20.0)
        self.assertEqual(left["W_damp_negative_J"], -5.0)
        self.assertEqual(checks["max_abs_interval_mode_decomposition_residual_J"], 0.0)


if __name__ == "__main__":
    unittest.main()
