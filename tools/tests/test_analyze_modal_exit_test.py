#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analyze_modal_exit_test import (  # noqa: E402
    applied_fluid_force,
    classify,
    integrate_decimated,
)


class ModalExitTestTests(unittest.TestCase):
    def test_applied_fluid_is_separated_from_computed_contact(self):
        frame = pd.DataFrame({
            "applied_force": [13.0, 14.0],
            "total_force": [12.0, 15.0],
            "fluid_force": [10.0, 11.0],
        })
        np.testing.assert_allclose(applied_fluid_force(frame), [11.0, 10.0])

    def test_decimated_work_preserves_raw_dynamic_static_identity(self):
        time = np.linspace(0.0, 1.0, 101)
        q = time * time
        velocity = 2.0 * time
        reference_force = 3.0
        varying_force = reference_force + time
        frame = pd.DataFrame({
            "state_time_s": time,
            "q": q,
            "qdot": velocity,
            "Q_fluid_applied": varying_force,
        })
        result = integrate_decimated(frame, 0.0, 1.0, reference_force, 0.5)
        self.assertTrue(np.isclose(result["W_raw_decimated_trapezoid"],
                                   3.0 + 2.0 / 3.0, atol=1.0e-4))
        self.assertTrue(np.isclose(result["W_static"], 3.0))
        self.assertTrue(np.isclose(
            result["W_raw_decimated_trapezoid"],
            result["W_dynamic_decimated_trapezoid"] + result["W_static"]))
        self.assertGreaterEqual(result["D_decimated_trapezoid"], 0.0)

    def test_classification_requires_error_separation(self):
        self.assertEqual(classify(0.4, 1.0, 0.05, 0.05, 0.35, 0.45),
                         "EXCITING_INSUFFICIENT")
        self.assertEqual(classify(-0.4, 1.0, 0.05, 0.05, -0.45, -0.35),
                         "FLUID_DAMPING")
        self.assertEqual(classify(1.4, 1.0, 0.05, 0.05, 1.35, 1.45),
                         "INPUT_AT_LEAST_DAMPING")
        self.assertEqual(classify(0.02, 1.0, 0.05, 0.05, -0.03, 0.07),
                         "NEAR_ZERO_OR_INCONCLUSIVE")
        self.assertEqual(classify(0.4, 1.0, 0.05, 0.05, -0.2, 0.6),
                         "NEAR_ZERO_OR_INCONCLUSIVE")


if __name__ == "__main__":
    unittest.main()
