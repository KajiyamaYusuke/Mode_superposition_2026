import importlib.util
import unittest
from pathlib import Path
import numpy as np

path=Path(__file__).resolve().parents[1]/'body_ratio_jitter_261103/analyze_flow_jitter.py'
spec=importlib.util.spec_from_file_location('flow_jitter',path)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

class PeriodTests(unittest.TestCase):
    def test_constant_and_alternating(self):
        self.assertAlmostEqual(m.period_metrics(np.ones(20))['local'],0)
        self.assertAlmostEqual(m.period_metrics([1.,2.,1.,2.,1.])['local'],100/1.4)
        self.assertGreater(m.period_metrics([1.,2.,1.,2.,1.])['rap'],0)
    def test_ppq_and_window(self):
        self.assertAlmostEqual(m.period_metrics([1,1,2,1,1])['ppq5'],100*(2-1.2)/1.2)
        np.testing.assert_allclose(m.window_periods([.5,.6,.7,.75],.5,.75),[.1,.1])
    def test_sine_and_no_oscillation(self):
        t=np.arange(20000)*5e-5
        p=m.window_periods(m.events(np.sin(2*np.pi*125*t)),.5,1.)
        self.assertGreater(len(p),10); self.assertLess(m.period_metrics(p)['local'],1e-6)
        self.assertEqual(len(m.events(np.ones(len(t)))),0)
    def test_validation(self):
        with self.assertRaises(ValueError): m.events([0,np.nan,1])
        with self.assertRaises(ValueError): m.period_metrics([1,0,1])
        self.assertTrue(np.isnan(m.period_metrics([1,1])['ppq5']))

class QualityTests(unittest.TestCase):
    def test_stationary_and_silent(self):
        t=np.arange(20000)*5e-5
        self.assertEqual(m.analyze(np.sin(2*np.pi*125*t))['status'],'VALID_STATIONARY')
        self.assertEqual(m.analyze(np.zeros(len(t)))['status'],'NO_OSCILLATION')
    def test_nonstationary_frequency(self):
        t=np.arange(20000)*5e-5
        phase=2*np.pi*(100*t+15*t*t)
        self.assertEqual(m.analyze(np.sin(phase))['status'],'VALID_WINDOW_DEPENDENT')
    def test_invalid_input_and_short(self):
        self.assertEqual(m.analyze([0,np.nan,1])['status'],'NUMERIC_INVALID')
        self.assertEqual(m.analyze([0,1,0])['status'],'INSUFFICIENT_CYCLES')
    def test_scale_invariance(self):
        t=np.arange(20000)*5e-5; x=np.sin(2*np.pi*123.4*t)
        self.assertAlmostEqual(m.analyze(x)['metrics']['local'],m.analyze(x*1e-5)['metrics']['local'])

class SyntheticFailureTests(unittest.TestCase):
    def test_missed_and_double_peaks_and_noise(self):
        t=np.arange(20000)*5e-5; phase=2*np.pi*125*t
        missing=np.where((t>.63)&(t<.64),0,np.sin(phase))
        double=np.sin(phase)+.9*np.sin(2*phase)
        noisy=np.sin(phase)+np.random.default_rng(123).normal(0,1,len(t))
        for x in (missing,double,noisy):
            self.assertEqual(m.analyze(x)['status'],'AMBIGUOUS_PEAKS')
    def test_amplitude_modulation(self):
        t=np.arange(20000)*5e-5
        x=(1+.4*np.sin(2*np.pi*2*t))*np.sin(2*np.pi*125*t)
        self.assertEqual(m.analyze(x)['status'],'VALID_WINDOW_DEPENDENT')
    def test_flow_step_axis(self):
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w') as f:
            f.write('# step airflow[m^3/s]\n0 0.1\n5 0.2\n10 0.3\n'); f.flush()
            t,x=m.load_flow(f.name,nstep=15)
            np.testing.assert_allclose(t,[0,5e-5,1e-4]); np.testing.assert_allclose(x,[.1,.2,.3])
            with self.assertRaises(ValueError): m.load_flow(f.name,nstep=20)

class SamplingTests(unittest.TestCase):
    def test_interpolation_reduces_quantization_jitter(self):
        t=np.arange(20000)*5e-5; x=np.sin(2*np.pi*137.3*t)
        p1=m.window_periods(m.events(x,interpolate=True),.5,1.)
        p0=m.window_periods(m.events(x,interpolate=False),.5,1.)
        self.assertLess(m.period_metrics(p1)['local'],m.period_metrics(p0)['local'])
    def test_methods_separate_phase_but_share_period(self):
        t=np.arange(20000)*5e-5; x=np.sin(2*np.pi*125*t)
        a=m.events(x,method='A');b=m.events(x,method='B')
        self.assertGreater(abs(a[0]-b[0]),.001)
        self.assertAlmostEqual(np.median(np.diff(a)),np.median(np.diff(b)),places=8)
    def test_sample_boundary(self):
        self.assertEqual(len(m.window_periods([.49,.5,.6,.75,.76],.5,.75)),1)

class IndependentEvidenceTests(unittest.TestCase):
    def test_periodicity_sine(self):
        t=np.arange(10000)*5e-5
        d=m.periodicity(np.sin(2*np.pi*125*t),5e-5,1/125)
        self.assertGreater(d['acf_at_period'],.9)
        self.assertGreater(d['spectral_ratio'],.9)
        self.assertLess(m.periodicity(np.sin(2*np.pi*125*t),5e-5,1/90)['spectral_ratio'],.05)
    def test_event_correspondence(self):
        self.assertEqual(m.event_correspondence([0,1,2,3],[.5,1.5,2.5],0,4),1)
        self.assertEqual(m.event_correspondence([0,1,2,3],[.5,2.5],0,4),2/3)

class AlternatingPulseTests(unittest.TestCase):
    def test_alternating_periods_are_not_called_random(self):
        t=np.arange(20000)*5e-5; x=np.zeros_like(t); event=.001; i=0
        while event<1:
            x+=np.exp(-.5*((t-event)/.00045)**2)
            event += .006 if i%2==0 else .010
            i+=1
        result=m.analyze(x)
        self.assertEqual(result['status'],'AMBIGUOUS_PEAKS')
