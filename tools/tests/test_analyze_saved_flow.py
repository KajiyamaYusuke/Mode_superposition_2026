import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

directory=Path(__file__).resolve().parents[1]/'body_ratio_jitter_261103'
sys.path.insert(0,str(directory))
s=importlib.util.spec_from_file_location('analyze_saved',directory/'analyze_saved.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

class SavedTests(unittest.TestCase):
    def test_complete_sine_and_preserve(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);d=p/'output/dat';d.mkdir(parents=True)
            t=np.arange(20000)*5e-5; vals=np.sin(2*np.pi*125*t)*1e-4; step=np.arange(20000)*5
            np.savetxt(d/'airflow_vt.dat',np.column_stack((step,vals)),header='step airflow[m^3/s]')
            (p/'status.json').write_text(json.dumps({'simulation_status':'SIMULATION_COMPLETED','validation':{'flow_sha256':m.run_sweep.sha(d/'airflow_vt.dat')}}))
            np.savetxt(d/'pressure_vt.dat',np.column_stack((step,vals)),header='step pressure')
            m.analyze_case(p); result=json.loads((p/'jitter_analysis/metrics.json').read_text())
            self.assertEqual(result['jitter_status'],'VALID_STATIONARY')
            self.assertLess(result['flow_jitter']['local'],.01)
            with self.assertRaises(FileExistsError):m.analyze_case(p)

class TemplateAssessmentTests(unittest.TestCase):
    def evaluate(self, x):
        estimate=m.template.estimate_period(x,5e-5);tracks={}
        for multiplier in (.5,1,2):
            z=m.template.track(x,5e-5,estimate['period_s']*multiplier);times=z['times']+.5;p=np.diff(times)
            tracks[str(multiplier)]={'times':times,'similarity':z['similarity'],'periods':p,'metrics':m.flow.period_metrics(p)}
        evidence=m.flow.periodicity(x,5e-5,estimate['period_s'])
        return m.assess_template(estimate,tracks,evidence,x)
    def test_known_random_jitter(self):
        dt=5e-5;t=np.arange(10000)*dt;rng=np.random.default_rng(42);times=[.01]
        while times[-1]<.49:times.append(times[-1]+.008*(1+rng.normal(0,.015)))
        x=np.zeros_like(t)
        for v in times:x+=np.exp(-.5*((t-v)/.0006)**2)+.3*np.exp(-.5*((t-v-.0015)/.00025)**2)
        result=self.evaluate(x); self.assertIn(result['status'],('VALID_STATIONARY','VALID_WINDOW_DEPENDENT'))
        p=np.diff(times);truth=np.mean(abs(np.diff(p)))/np.mean(p)*100
        self.assertAlmostEqual(result['metrics']['local'],truth,delta=.15)
    def test_missing_cycles_rejected(self):
        dt=5e-5;t=np.arange(10000)*dt;x=np.zeros_like(t)
        for v in np.arange(.01,.49,.008):
            if .25<v<.27:continue
            x+=np.exp(-.5*((t-v)/.0006)**2)
        result=self.evaluate(x);self.assertEqual(result['status'],'UNDEFINED_JITTER')
    def test_strong_noise_rejected(self):
        t=np.arange(10000)*5e-5;x=np.sin(2*np.pi*125*t)+np.random.default_rng(21).normal(0,1,len(t))
        try:r=self.evaluate(x)
        except ValueError:return
        self.assertEqual(r['status'],'UNDEFINED_JITTER')
