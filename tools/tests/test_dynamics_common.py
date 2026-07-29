#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dynamics_common import (  # noqa: E402
    cycle_table,
    detect_events,
    replace_pressure_parameter,
    spectrum_metrics,
)


class DynamicsCommonTests(unittest.TestCase):
    def setUp(self):
        self.dt = 1.0e-4
        self.time = np.arange(0.0, 1.0, self.dt)

    def cycles(self, values):
        events, _ = detect_events(
            self.time, values, "rising_threshold", section_level=0.5
        )
        return cycle_table(self.time, values, events)

    def test_period_one_has_stable_period_and_amplitude(self):
        values = 2.0 + np.sin(2 * np.pi * 40 * self.time)
        cycles = self.cycles(values)
        self.assertGreater(len(cycles), 30)
        self.assertLess(cycles["period_s"].std() / cycles["period_s"].mean(), 1e-8)
        self.assertLess(cycles["amplitude_mm2"].std(), 1e-8)

    def test_amplitude_modulation_changes_amplitude_more_than_period(self):
        carrier = np.sin(2 * np.pi * 40 * self.time)
        values = 2.0 + (1.0 + 0.3 * np.sin(2 * np.pi * 3 * self.time)) * carrier
        cycles = self.cycles(values)
        period_cv = cycles["period_s"].std() / cycles["period_s"].mean()
        amplitude_cv = cycles["amplitude_mm2"].std() / cycles["amplitude_mm2"].mean()
        self.assertGreater(amplitude_cv, 5 * period_cv)

    def test_alternating_cycles_form_two_amplitude_branches(self):
        phase = 2 * np.pi * 40 * self.time
        cycle_number = np.floor(40 * self.time).astype(int)
        amplitude = np.where(cycle_number % 2 == 0, 1.0, 0.55)
        cycles = self.cycles(2.0 + amplitude * np.sin(phase))
        even = cycles["amplitude_mm2"].to_numpy()[::2]
        odd = cycles["amplitude_mm2"].to_numpy()[1::2]
        self.assertGreater(abs(even.mean() - odd.mean()), 0.5)

    def test_quasiperiodic_signal_has_two_nonharmonic_spectral_peaks(self):
        values = np.sin(2 * np.pi * 40 * self.time) + 0.6 * np.sin(
            2 * np.pi * (40 * np.sqrt(2)) * self.time
        )
        _, metrics = spectrum_metrics(self.time, values)
        peaks = sorted(
            [metrics["dominant_frequency_1_hz"], metrics["dominant_frequency_2_hz"]]
        )
        self.assertAlmostEqual(peaks[0], 40, delta=2)
        self.assertAlmostEqual(peaks[1], 40 * np.sqrt(2), delta=2)

    def test_pressure_replacement_changes_only_target_value(self):
        original = "# another value\n1600\n# subglottal pressure Ps (Pa)\n1300\n# end\n"
        replaced = replace_pressure_parameter(original, 1450)
        self.assertIn("# another value\n1600", replaced)
        self.assertIn("# subglottal pressure Ps (Pa)\n1450", replaced)


if __name__ == "__main__":
    unittest.main()
