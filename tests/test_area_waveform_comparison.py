"""Pure postprocessing tests; no coupled solver invocation."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/kanaya_ablation"))
from compare_area_waveforms import minimum_area, opening_peaks, spectrum


class AreaWaveformTests(unittest.TestCase):
    def test_invalid_sections_are_not_zero(self):
        rows = np.array([[0, np.nan, 3, 4], [5, 5, np.nan, 2]])
        np.testing.assert_array_equal(minimum_area(rows), [3, 2])

    def test_bad_area_rejected(self):
        for row in [[0, np.nan, np.nan], [0, np.inf, 1], [0, -1, 2]]:
            with self.assertRaises(ValueError):
                minimum_area(np.array([row]))

    def test_known_frequency_with_mean_offset(self):
        time = np.arange(2000) * 5e-5
        area = 5 + 2 * np.cos(2 * np.pi * 120 * time + .4)
        frequencies, _, index = spectrum(area, 5e-5)
        self.assertAlmostEqual(frequencies[index], 120)
        peaks, _, periods, frequency = opening_peaks(area, time, .20)
        # The final incomplete peak lacks the required right-hand prominence.
        self.assertEqual(len(peaks), 11)
        self.assertAlmostEqual(frequency, 120, delta=.1)
        self.assertAlmostEqual(periods.mean(), 1 / 120, delta=5e-5)


if __name__ == "__main__":
    unittest.main()
