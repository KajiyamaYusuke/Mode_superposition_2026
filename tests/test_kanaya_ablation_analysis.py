"""Guard against reporting an open equilibrium as periodic closure."""
import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/kanaya_ablation'))
from compare_runs import closure_periods, numeric_gate, phase_at_frequency


class AnalysisTests(unittest.TestCase):
    def test_open_equilibrium_has_no_cq_or_f0(self):
        result=closure_periods(np.ones(4999),1e-5)
        self.assertFalse(result['periodic_closure'])
        self.assertTrue(np.isnan(result['CQ']))
        self.assertTrue(np.isnan(result['f0_from_closures_Hz']))

    def test_permanently_closed_is_not_periodic(self):
        self.assertFalse(closure_periods(np.zeros(4999),1e-5)['periodic_closure'])

    def test_repeated_closure_has_measured_period(self):
        area=np.tile(np.r_[np.zeros(30),np.ones(70)],50)
        result=closure_periods(area,1e-5)
        self.assertTrue(result['periodic_closure'])
        self.assertAlmostEqual(result['CQ'],.3)
        self.assertAlmostEqual(result['f0_from_closures_Hz'],1000.)

    def test_numerical_gate_rejects_asymmetry_and_nan(self):
        for a,b in [([1.],[1.1]),([np.nan],[np.nan])]:
            with self.assertRaises(AssertionError):numeric_gate(a,b)

    def test_work_mean_and_fluctuation_add_to_total(self):
        rng=np.random.default_rng(17)
        force=rng.normal(size=(500,20))+10
        dq=rng.normal(size=(500,20))*1e-7
        total=np.sum(force*dq,axis=0)
        mean=force.mean(axis=0)*dq.sum(axis=0)
        fluct=np.sum((force-force.mean(axis=0))*dq,axis=0)
        numeric_gate(total,mean+fluct,atol=1e-18,rtol=1e-10)

    def test_downstream_nr_pressure_telescope(self):
        rng=np.random.default_rng(42)
        flow=rng.normal(size=1000)*1e-3
        downstream_previous=rng.normal(size=1000)*1e-3
        previous=np.r_[0.,flow[:-1]]
        factor=1.3e5
        accumulator=np.cumsum(previous-downstream_previous)
        pd=np.cumsum(factor*(flow-downstream_previous))
        numeric_gate(factor*accumulator,pd-factor*flow,atol=1e-9,rtol=1e-12)

    def test_force_velocity_phase_sign(self):
        dt=1e-4;t=np.arange(10000)*dt
        velocity=np.sin(2*np.pi*120*t)
        force=np.sin(2*np.pi*120*t+np.pi/6)
        self.assertAlmostEqual(phase_at_frequency(force,velocity,dt,120.),30.,places=3)


if __name__=='__main__':unittest.main()
