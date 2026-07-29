#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from left_right_ratios import (  # noqa: E402
    RatioConfig,
    calculate_unsigned_ratio,
    detect_cycles,
    estimate_dominant_frequency,
)


class LeftRightRatioTests(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(0.0, 1.0, 1.0e-4)
        self.config = RatioConfig(smooth="none")

    def summarize(self, values):
        centered = values - np.mean(values)
        frequency = estimate_dominant_frequency(self.time, centered, 40.0, 500.0)
        peaks, cycles = detect_cycles(
            self.time, centered, centered, frequency, self.config
        )
        periods = cycles.loc[cycles["period_valid"], "period_raw_s"].to_numpy()
        amplitudes = cycles["amplitude_mm"].to_numpy()
        self.assertGreaterEqual(len(peaks), self.config.min_peak_count)
        return np.median(periods), np.median(amplitudes)

    def assert_ratios(self, left, right, expected_period, expected_amplitude):
        left_period, left_amplitude = self.summarize(left)
        right_period, right_amplitude = self.summarize(right)
        _, period_ratio = calculate_unsigned_ratio(left_period, right_period)
        amplitude_signed, amplitude_ratio = calculate_unsigned_ratio(
            left_amplitude, right_amplitude
        )
        self.assertAlmostEqual(period_ratio, expected_period, delta=0.03)
        self.assertAlmostEqual(amplitude_ratio, expected_amplitude, delta=0.03)
        return amplitude_signed

    def test_same_period_and_amplitude(self):
        wave = np.sin(2 * np.pi * 100 * self.time)
        self.assert_ratios(wave, wave, 1.0, 1.0)

    def test_two_to_one_period(self):
        left = np.sin(2 * np.pi * 100 * self.time)
        right = np.sin(2 * np.pi * 50 * self.time)
        self.assert_ratios(left, right, 2.0, 1.0)

    def test_two_to_one_amplitude(self):
        left = 2 * np.sin(2 * np.pi * 100 * self.time)
        right = np.sin(2 * np.pi * 100 * self.time)
        signed = self.assert_ratios(left, right, 1.0, 2.0)
        self.assertAlmostEqual(signed, 2.0, delta=0.03)

    def test_period_and_amplitude_two_to_one_with_noise(self):
        generator = np.random.default_rng(20260727)
        left = 2 * np.sin(2 * np.pi * 100 * self.time)
        right = np.sin(2 * np.pi * 50 * self.time)
        left += generator.normal(0, 0.01, len(self.time))
        right += generator.normal(0, 0.01, len(self.time))
        self.assert_ratios(left, right, 2.0, 2.0)


if __name__ == "__main__":
    unittest.main()
