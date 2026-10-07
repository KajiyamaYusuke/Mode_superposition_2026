#!/usr/bin/env python3

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spectol_sweep import (  # noqa: E402
    condition_tag,
    flow_m3s_to_ls,
    load_complete_rows,
    minimum_glottal_area,
    update_param_file,
)
from steinecke_classification import (  # noqa: E402
    classify_locking,
    detect_peak_period,
    next_maximum_map,
)


class SpectolSweepUnitTests(unittest.TestCase):
    def test_cubic_metres_per_second_to_litres_per_second(self):
        converted = flow_m3s_to_ls(np.array([0.0, 1.0e-3, 2.5e-4]))
        np.testing.assert_allclose(converted, [0.0, 1.0, 0.25])

    def test_minimum_glottal_area_ignores_step_column(self):
        data = np.array([[0.0, 4.0, 2.0, 3.0], [5.0, 1.5, 2.5, 4.0]])
        np.testing.assert_allclose(minimum_glottal_area(data), [2.0, 1.5])

    def test_incomplete_trailing_row_is_ignored(self):
        with tempfile.NamedTemporaryFile(mode='w+', encoding='utf-8') as stream:
            stream.write('# step value\n0 1.0\n5 2.0\n10\n')
            stream.flush()
            result = load_complete_rows(stream.name, 2)
        np.testing.assert_allclose(result, [[0.0, 1.0], [5.0, 2.0]])

    def test_pressure_and_both_damping_values_are_updated(self):
        text = (
            "# damping coefficient zetaL\n0.05\n"
            "# damping coefficient zetaR\n0.06\n"
            "# subglottal pressure Ps (Pa)\n600\n"
        )
        with tempfile.NamedTemporaryFile(
            mode='w+', encoding='utf-8', delete=False
        ) as stream:
            stream.write(text)
            path = Path(stream.name)
        try:
            update_param_file(path, 750, 0.025)
            updated = path.read_text(encoding='utf-8')
        finally:
            path.unlink()
        self.assertIn("# damping coefficient zetaL\n0.025\n", updated)
        self.assertIn("# damping coefficient zetaR\n0.025\n", updated)
        self.assertIn("# subglottal pressure Ps (Pa)\n750\n", updated)

    def test_condition_tag_contains_pressure_and_damping(self):
        self.assertEqual(condition_tag(600, 0.025), "Ps_600Pa_zeta_0p025")

    def test_peak_period_detects_alternating_period_two(self):
        times = np.arange(12, dtype=float) * 0.01
        values = np.tile([1.0, 0.7], 6)
        period, errors = detect_peak_period(
            times, values, signal_scale=1.0, amplitude_rtol=1e-10,
            interval_rtol=1e-10,
        )
        self.assertEqual(period, 2)
        self.assertIn(2, errors)

    def test_next_maximum_map_preserves_order(self):
        x_value, next_value = next_maximum_map(np.array([1.0, 0.7, 1.1]))
        np.testing.assert_allclose(x_value, [1.0, 0.7])
        np.testing.assert_allclose(next_value, [0.7, 1.1])

    def test_classification_does_not_reduce_period_two_label(self):
        time = np.arange(0.0, 2.0, 0.0005)
        carrier = np.sin(2.0 * np.pi * 20.0 * time)
        envelope = 1.0 + 0.25 * np.cos(2.0 * np.pi * 10.0 * time)
        displacement = envelope * carrier
        result = classify_locking(
            time, displacement, displacement, transient_end=0.2,
            min_peaks=8, amplitude_rtol=0.03, interval_rtol=0.03,
        )
        self.assertEqual(result.get("base_label"), "2:2")
        self.assertNotEqual(result.get("base_label"), "1:1")


if __name__ == "__main__":
    unittest.main()
