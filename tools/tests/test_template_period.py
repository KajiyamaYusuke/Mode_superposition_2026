import importlib.util
from pathlib import Path
import unittest
import numpy as np
p=Path(__file__).resolve().parents[1]/'body_ratio_jitter_261103/template_period.py'
s=importlib.util.spec_from_file_location('template_period',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

class TemplateTests(unittest.TestCase):
    def test_regular_wave(self):
        dt=5e-5;t=np.arange(10000)*dt;x=np.sin(2*np.pi*125*t)
        r=m.estimate_period(x,dt);self.assertAlmostEqual(r['period_s'],.008,delta=dt)
        z=m.track(x,dt,r['period_s']);self.assertGreater(len(z['times']),50)
        self.assertAlmostEqual(np.median(np.diff(z['times'])),.008,delta=dt)
        self.assertGreater(np.min(z['similarity']),.99)
    def test_jitter_not_suppressed(self):
        dt=5e-5;t=np.arange(20000)*dt;times=[.01]
        while times[-1]<.99:
            times.append(times[-1]+(.0075 if len(times)%2 else .0085))
        x=np.zeros_like(t)
        for peak in times:x+=np.exp(-.5*((t-peak)/.0006)**2)
        r=m.estimate_period(x,dt)
        # ACF prefers 2T here; treating it as T suppresses alternating jitter.
        self.assertEqual(r['status'],'AMBIGUOUS_PERIOD')
        z=m.track(x,dt,r['period_s']/2); periods=np.diff(z['times'])
        self.assertGreater(np.mean(np.abs(np.diff(periods))),.0007)
    def test_silent_invalid(self):
        with self.assertRaises(ValueError):m.estimate_period(np.zeros(1000),.00005)

    def test_known_irregular_pulses_with_secondary_peak(self):
        dt=5e-5; t=np.arange(20000)*dt; rng=np.random.default_rng(42); truth=[.01]
        while truth[-1]<.99: truth.append(truth[-1]+.008*(1+rng.normal(0,.015)))
        x=np.zeros_like(t)
        for peak in truth:
            x+=np.exp(-.5*((t-peak)/.0006)**2)+.3*np.exp(-.5*((t-peak-.0015)/.00025)**2)
        r=m.estimate_period(x,dt); self.assertEqual(r['status'],'CANDIDATE')
        found=m.track(x,dt,r['period_s'])['times']; truth=np.array(truth)
        self.assertLess(np.max(np.min(abs(found[:,None]-truth[None,:]),axis=1)),dt)
        def jitter(z):
            periods=np.diff(z)
            return np.mean(abs(np.diff(periods)))/np.mean(periods)*100
        self.assertAlmostEqual(jitter(found),jitter(truth),delta=.05)
