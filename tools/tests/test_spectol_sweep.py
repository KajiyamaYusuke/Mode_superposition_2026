#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spectol_sweep import flow_m3s_to_ls  # noqa: E402


class SpectolSweepUnitTests(unittest.TestCase):
    def test_cubic_metres_per_second_to_litres_per_second(self):
        converted = flow_m3s_to_ls(np.array([0.0, 1.0e-3, 2.5e-4]))
        np.testing.assert_allclose(converted, [0.0, 1.0, 0.25])


if __name__ == "__main__":
    unittest.main()
