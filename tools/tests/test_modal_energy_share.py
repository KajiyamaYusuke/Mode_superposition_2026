#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from plot_modal_energy_share import calculate_energy_share, summarize_modes  # noqa: E402


class ModalEnergyShareTests(unittest.TestCase):
    def test_energy_share_and_ranking(self):
        rows = []
        for time in (0.0, 0.1):
            rows.extend(
                [
                    {"time_s": time, "side": "L", "mode_index": 0,
                     "frequency_hz": 1.0, "q": 1.0, "qdot": 0.0},
                    {"time_s": time, "side": "L", "mode_index": 1,
                     "frequency_hz": 2.0, "q": 1.0, "qdot": 0.0},
                ]
            )
        energy = calculate_energy_share(pd.DataFrame(rows), None, None)
        grouped = energy.groupby("time_s")["energy_share"].sum().to_numpy()
        np.testing.assert_allclose(grouped, 1.0)
        summary = summarize_modes(energy)
        leader = summary.iloc[0]
        self.assertEqual(int(leader["mode_index"]), 1)
        self.assertEqual(int(leader["mean_energy_rank"]), 1)
        self.assertAlmostEqual(float(leader["mean_energy_share_percent"]), 80.0)

    def test_time_interval_is_applied(self):
        frame = pd.DataFrame(
            {
                "time_s": [0.0, 0.1, 0.2], "side": ["L"] * 3,
                "mode_index": [0] * 3, "frequency_hz": [100.0] * 3,
                "q": [1.0] * 3, "qdot": [0.0] * 3,
            }
        )
        energy = calculate_energy_share(frame, 0.1, 0.1)
        self.assertEqual(energy["time_s"].tolist(), [0.1])
        self.assertEqual(energy["energy_share_percent"].tolist(), [100.0])


if __name__ == "__main__":
    unittest.main()
