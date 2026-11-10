#!/usr/bin/env python3
"""Focused regression tests for the symmetric/asymmetric postprocessor."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/amplitude_asymmetry"))
from compare_symmetric_asymmetric import (  # noqa: E402
    covariance_metrics, mac_matrix, match_modes, signal_metrics, work_components,
)


class CompareDiagnosticsTests(unittest.TestCase):
    def test_mirrored_synthetic_signals_are_symmetric(self) -> None:
        left_solver_y = np.array([-1.0, 0.0, 1.0, 0.0])
        right_solver_y = -left_solver_y
        left = signal_metrics(-left_solver_y)
        right = signal_metrics(right_solver_y)
        self.assertEqual(left, right)
        self.assertEqual(right["half_ptp_mm"] / left["half_ptp_mm"], 1.0)

    def test_covariance_identity_and_phase_sign(self) -> None:
        x = np.linspace(-1, 1, 1001)
        reinforcing = covariance_metrics(np.column_stack([x, x]))
        canceling = covariance_metrics(np.column_stack([x, -x]))
        self.assertGreater(reinforcing["cross"], 0)
        self.assertLess(canceling["cross"], 0)
        self.assertTrue(np.isclose(reinforcing["total"],
                                   reinforcing["individual_sum"] + reinforcing["cross"]))
        self.assertTrue(np.isclose(canceling["total"],
                                   canceling["individual_sum"] + canceling["cross"]))
        self.assertTrue(np.isclose(canceling["total"], 0.0, atol=1e-15))

    def test_mac_accepts_sign_flip(self) -> None:
        reference = np.eye(2).reshape(2, 2, 1)
        candidate = -reference.copy()
        self.assertTrue(np.allclose(mac_matrix(reference, candidate), np.eye(2)))
        matches, _ = match_modes(reference, candidate, np.array([10., 20.]),
                                 np.array([10., 20.]))
        self.assertEqual(matches.match_quality.tolist(), ["high", "high"])
        self.assertEqual(matches.matched_mode.tolist(), [1, 2])

    def test_ambiguous_mode_is_not_forced(self) -> None:
        reference = np.eye(2).reshape(2, 2, 1)
        candidate = np.array([[[1.0], [1.0]], [[1.0], [-1.0]]])
        matches, _ = match_modes(reference, candidate, np.array([10., 20.]),
                                 np.array([15., 15.]))
        self.assertTrue(matches.ambiguous.all())
        self.assertTrue(matches.matched_mode.isna().all())

    def test_q_dq_split_qref_and_contact_zero(self) -> None:
        data = np.zeros((4, 2, 7))
        data[:, :, 0] = 2.0
        data[:, :, 1] = -0.5
        data[:, :, 2] = np.array([1., 2., 3., 4.])[:, None]
        data[:, :, 5] = 0.25
        data[:, :, 6] = 0.1
        rows = work_components(data, np.array([10., 20.]), 1e-3)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(np.isclose(row["W_pos_J"], 0.8) for row in rows))
        self.assertTrue(all(np.isclose(row["W_neg_J"], -0.2) for row in rows))
        self.assertTrue(all(np.isclose(row["W_raw_J"], 1.0) for row in rows))
        self.assertTrue(all(np.isclose(row["W_fluct_Qmean_J"], 0.0) for row in rows))
        self.assertTrue(all(row["W_contact_ref_J"] == 0 for row in rows))
        self.assertTrue(np.isclose(sum(row["W_raw_J"] for row in rows), 2.0))

    def test_existing_output_directory_is_refused(self) -> None:
        script = ROOT / "tools/amplitude_asymmetry/compare_symmetric_asymmetric.py"
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([
                sys.executable, str(script), "--asym-analysis", "missing",
                "--asym-production", "missing", "--sym-source", "missing",
                "--sym-production", "missing", "--output-dir", directory,
            ], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("refusing to overwrite", result.stderr)


if __name__ == "__main__":
    unittest.main()

