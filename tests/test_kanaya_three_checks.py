"""Synthetic checks only. No simulation executable or production data required."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

SPEC = importlib.util.spec_from_file_location('three_checks', Path(__file__).resolve().parents[1]/'tools/kanaya_diagnostics/three_checks.py')
POST = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(POST)


def surfaces(gap=2.):
    left = np.array([[[x, 0., z] for z in np.linspace(0, 17, 5)] for x in (0., 1., 2.)])
    right = left.copy(); right[..., 1] = gap
    return left, right


class GeometryChecks(unittest.TestCase):
    def test_small_width_loss_can_reverse_minimum_ranking(self):
        import pandas as pd
        spans = np.array([17., 16.3, 17.])
        gaps = np.array([.384, .398, .5])
        table = pd.DataFrame({'step': [10000]*3, 'plane_index': [41, 42, 43],
            'x': [10., 10.25, 10.5], 'valid': [True]*3,
            'common_span_mm': spans, 'mean_gap_common_mm': gaps,
            'A_reconstructed': spans*gaps})
        row = POST.ap_selection_sensitivity(table)[0]
        self.assertTrue(row['different_plane'])
        self.assertEqual(row['area_min_plane'], 42)
        self.assertEqual(row['gap_min_plane'], 41)
        self.assertAlmostEqual(row['x_shift_mm'], -.25)

    def test_constant_width_constant_gap(self):
        left, right = surfaces()
        x, cuts = POST.evaluate_sections(left, right, count=5)
        np.testing.assert_allclose(x, np.linspace(0, 2, 5))
        for cut in cuts:
            self.assertAlmostEqual(cut['A'], 34.)
            self.assertAlmostEqual(cut['G3'], 136.)
            self.assertAlmostEqual(cut['S_common'], 17.)
            self.assertAlmostEqual(cut['S_open'], 17.)

    def test_missing_one_side_ap_row_does_not_bridge(self):
        left, right = surfaces()
        right[:, 2, 0] += 10
        cut = POST.integrate_segments(POST.section_segments(left, 1., True),
                                      POST.section_segments(right, 1., False))
        self.assertAlmostEqual(cut['S_common'], 8.5)
        self.assertAlmostEqual(cut['A'], 17.)
        self.assertAlmostEqual(cut['A']/cut['S_common'], 2.)

    def test_true_constriction_constant_ap_width(self):
        left, right = surfaces()
        right[1, :, 1] = .5
        _, cuts = POST.evaluate_sections(left, right, count=3)
        self.assertAlmostEqual(cuts[1]['S_common'], 17.)
        self.assertAlmostEqual(cuts[1]['A'], 8.5)
        self.assertLess(cuts[1]['A'], cuts[0]['A'])

    def test_partial_contact_signed_linear_gap(self):
        left, right = surfaces()
        right[..., 1] = np.linspace(-1, 1, 5)
        cut = POST.integrate_segments(POST.section_segments(left, 1, True), POST.section_segments(right, 1, False))
        self.assertAlmostEqual(cut['S_common'], 17.)
        self.assertAlmostEqual(cut['S_open'], 8.5)
        self.assertAlmostEqual(cut['A'], 4.25)
        self.assertAlmostEqual(cut['G3'], 2.125)

    def test_no_intersection_is_nan_not_closed(self):
        left, right = surfaces()
        cuts = POST.integrate_segments(POST.section_segments(left, 3, True), POST.section_segments(right, 3, False))
        self.assertTrue(np.isnan(cuts['A']))
        self.assertEqual(cuts['S_common'], 0)
        self.assertFalse(cuts['A'] <= .01)

    def test_isolated_intersection_has_zero_width(self):
        left, right = surfaces()
        right[:, (0, 1, 3, 4), 0] += 10
        cut = POST.integrate_segments(POST.section_segments(left, 1, True), POST.section_segments(right, 1, False))
        self.assertTrue(np.isnan(cut['A']))
        self.assertEqual(cut['S_common'], 0)

    def test_nan_geometry_not_zero_area(self):
        left, right = surfaces()
        right[:] = np.nan
        cut = POST.integrate_segments(POST.section_segments(left, 1, True), POST.section_segments(right, 1, False))
        self.assertTrue(np.isnan(cut['A']))
        with self.assertRaises(ValueError): POST.positive_integrals(np.nan, 1, 2)

    def test_inner_endpoint_at_plane_aligned_edge(self):
        left, right = surfaces()
        left[0, :, 0] = 1.
        left[0, :, 1] = .5
        cut = POST.integrate_segments(POST.section_segments(left, 1, True), POST.section_segments(right, 1, False))
        self.assertAlmostEqual(cut['A'], 25.5)


class WorkChecks(unittest.TestCase):
    def test_constant_force_monotonic_displacement(self):
        force = np.full(100, 3.)
        dq = np.linspace(.01, .03, 100)
        total, mean, fluctuating = POST.work_decomposition(force, dq)
        self.assertAlmostEqual(total, 6.)
        self.assertAlmostEqual(mean, total)
        self.assertEqual(fluctuating, 0.)

    def test_in_phase_and_antiphase_sinusoidal_work(self):
        t = np.arange(10000)/10000
        velocity = np.sin(2*np.pi*100*t)
        for sign in (-1, 1):
            total, mean, fluctuating = POST.work_decomposition(sign*velocity, velocity/10000)
            self.assertAlmostEqual(total, sign*.5)
            self.assertAlmostEqual(mean, 0.)
            self.assertAlmostEqual(fluctuating, sign*.5)

    def test_constant_opening_no_vibration(self):
        q = np.full(100, .1)
        total, mean, fluctuating = POST.work_decomposition(np.ones(99), np.diff(q))
        self.assertEqual(total, 0.); self.assertEqual(fluctuating, 0.)
        self.assertAlmostEqual(np.std(q), 0., places=15)

    def test_vtu_one_step_offset_fails(self):
        points = np.array([[0., 0., 0.], [1., 2., 3.]])
        modes = np.ones((1, 2, 3))
        q0, q1 = np.array([0.]), np.array([.001])
        vtu_after = POST.displace(points, modes, q1)
        POST.close_gate(vtu_after, POST.displace(points, modes, q1))
        with self.assertRaises(AssertionError):
            POST.close_gate(vtu_after, POST.displace(points, modes, q0))

    def test_static_frequency_radians(self):
        omega = np.array([2*np.pi*50.])
        q = np.full((100, 1), .02)
        force = omega**2*q
        _, _, _, predicted, residual, _ = POST.static_prediction(q, np.zeros_like(q), force, omega, .06, .001)
        np.testing.assert_allclose(predicted, [.02])
        np.testing.assert_allclose(residual, [0.], atol=1e-10)

    def test_nonstationary_inertia_damping_correction(self):
        dt = .0001
        t = np.arange(1000)*dt
        omega, zeta = np.array([2*np.pi*50.]), .06
        q = (.003+.1*t+.5*2*t*t)[:, None]
        v = (.1+2*t)[:, None]
        force = 2+2*zeta*omega*v+omega**2*q
        _, _, acc, _, residual, corrected = POST.static_prediction(q, v, force, omega, zeta, dt)
        self.assertAlmostEqual(acc[0], 2.)
        self.assertGreater(abs(residual[0]), 1.)
        np.testing.assert_allclose(corrected, [0], atol=1e-10)

    def test_sign_split_and_window_additivity(self):
        t = np.arange(500)/500
        pos, neg = 2+np.sin(2*np.pi*t), -1+.5*np.cos(2*np.pi*t)
        dq = np.cos(6*np.pi*t)/500+.01
        allw = POST.work_decomposition(pos+neg, dq)[0]
        pieces = sum(POST.work_decomposition((pos+neg)[sl], dq[sl])[0] for sl in (slice(0, 250), slice(250, 500)))
        self.assertAlmostEqual(allw, pieces)
        self.assertAlmostEqual(allw, np.sum(pos*dq)+np.sum(neg*dq))


if __name__ == '__main__':
    unittest.main()
