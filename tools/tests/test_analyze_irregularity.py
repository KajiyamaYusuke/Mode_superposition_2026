#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analyze_irregularity import peak_table, phase_data, repeat_errors, spectrum  # noqa: E402


class IrregularityAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(0.0, 0.5, 1.0e-4)

    def test_periodic_peak_period_and_amplitude(self):
        values = 2.0 * np.sin(2.0 * np.pi * 100.0 * self.time)
        table = peak_table(self.time, values, 0.05)
        self.assertAlmostEqual(float(table.period_s.mean()), 0.01, places=6)
        self.assertAlmostEqual(float(table.peak_to_trough_amplitude.mean()), 4.0, places=4)

    def test_constant_phase_offset_has_no_drift_or_slips(self):
        left = np.sin(2.0 * np.pi * 100.0 * self.time)
        right = np.sin(2.0 * np.pi * 100.0 * self.time + 0.4)
        _, metrics = phase_data(self.time, left, right)
        self.assertAlmostEqual(metrics["mean_phase_diff_rad"], 0.4, places=3)
        self.assertAlmostEqual(metrics["phase_drift_rate_rad_s"], 0.0, places=3)
        self.assertEqual(metrics["phase_slip_count"], 0)

    def test_frequency_mismatch_produces_phase_drift(self):
        left = np.sin(2.0 * np.pi * 100.0 * self.time)
        right = np.sin(2.0 * np.pi * 102.0 * self.time)
        _, metrics = phase_data(self.time, left, right)
        self.assertAlmostEqual(metrics["phase_drift_rate_rad_s"], 4.0 * np.pi, delta=0.2)
        self.assertGreaterEqual(metrics["phase_slip_count"], 1)

    def test_repeat_error_identifies_period_two(self):
        rows = repeat_errors(np.array([1.0, 2.0] * 10), 3)
        self.assertGreater(rows[0]["mean_absolute_error"], 0.9)
        self.assertEqual(rows[1]["mean_absolute_error"], 0.0)

    def test_spectrum_peak(self):
        values = np.sin(2.0 * np.pi * 120.0 * self.time)
        table = spectrum(self.time, values)
        peak = table.loc[table.amplitude.idxmax(), "frequency_hz"]
        self.assertAlmostEqual(float(peak), 120.0, places=6)


if __name__ == "__main__":
    unittest.main()
