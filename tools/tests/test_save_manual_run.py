import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

dirpath = Path(__file__).resolve().parents[1] / 'body_ratio_jitter_261103'
sys.path.insert(0, str(dirpath))
spec = importlib.util.spec_from_file_location('save_manual_run', dirpath / 'save_manual_run.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def fixture(root, duration_steps=10000, dt=1e-4, left_scale=2.0):
    run = root / 'output'; (run/'dat').mkdir(parents=True); (run/'txt').mkdir()
    param = '\n'.join(['20','18',str(duration_steps),'20',str(dt),'.06','.06',
                       '5 6000 5000','100','freq.txt','mode.vtu','surf.txt',
                       'input','result','0','900','1.13','1.83e-5','.001','343',
                       '0','100','-0.00001','leftModeFile = /tmp/b12c3.vtu',
                       'rightModeFile = /tmp/b6c3.vtu'])+'\n'
    (run/'txt/params_used.txt').write_text(param)
    steps=np.arange(0,duration_steps,5);t=steps*dt; x=np.sin(2*np.pi*100*t)
    np.savetxt(run/'dat/airflow_vt.dat',np.column_stack((steps,1e-4*x)),header='step airflow[m^3/s]')
    np.savetxt(run/'dat/displace.dat',np.column_stack((t,left_scale*x,.5*x)),header='time_s uyL_mm uyR_mm')
    np.savetxt(run/'dat/pressure_vt.dat',np.column_stack((steps,50*x)),header='step outlet_pressure[Pa]')
    return run


class ArchiveTests(unittest.TestCase):
    def test_finished_run_repeated_and_ratio(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root);dest=root/'archive'
            first,report=m.archive_run(run,dest,label='P900')
            second,_=m.archive_run(root,dest,label='P900')
            self.assertNotEqual(first,second)
            self.assertAlmostEqual(report['amplitude']['right_over_left'],.25,places=10)
            self.assertEqual(report['parameter_settings']['pressure_pa'],900)
            self.assertEqual(report['parameter_provenance'],'solver_recorded_params_used')
            self.assertEqual(report['window_samples'],1000)
            self.assertAlmostEqual(report['input_metadata']['body_token_inference']['right_over_left'],.5)
            self.assertTrue((first/'snapshot/dat/airflow_vt.dat').is_file())
            self.assertTrue((first/'waveforms.png').is_file())
            with (dest/'summary.csv').open() as f: self.assertEqual(len(list(csv.DictReader(f))),2)
            with (first/'waveforms.csv').open() as f: self.assertEqual(len(list(csv.DictReader(f))),2000)
    def test_mismatch_and_partial_are_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root); dest=root/'archive';other=root/'wrong.txt';other.write_text('different')
            with self.assertRaisesRegex(ValueError,'differs'):m.archive_run(run,dest,param_path=other)
            p=run/'dat/airflow_vt.dat';p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n')
            with self.assertRaisesRegex(ValueError,'incomplete'):m.archive_run(run,dest)
            self.assertFalse((dest/'summary.csv').exists())
            self.assertEqual(list(dest.glob('.building_*')),[])
    def test_user_param_and_custom_window(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root,duration_steps=3000);param=root/'param.txt';(run/'txt/params_used.txt').rename(param)
            out,report=m.archive_run(run,root/'archive',param_path=param,start=0,end=.25)
            self.assertEqual(report['parameter_provenance'],'user_supplied_not_independently_verified')
            self.assertEqual(report['window_s'],[0,.25]);self.assertEqual(report['window_samples'],500)
            self.assertTrue((out/'flow_cycles.csv').exists())
    def test_zero_amplitude_is_undefined(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root,left_scale=0)
            _,report=m.archive_run(run,root/'archive')
            self.assertIsNone(report['amplitude']['right_over_left'])
    def test_reject_incomplete_status(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root);(root/'status.json').write_text(json.dumps({'simulation_status':'INTERRUPTED'}))
            with self.assertRaisesRegex(ValueError,'INTERRUPTED'):m.archive_run(run,root/'archive')

    def test_reject_stale_completed_status(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run=fixture(root)
            (root/'status.json').write_text(json.dumps({'simulation_status':'SIMULATION_COMPLETED',
                                                        'validation':{'flow_sha256':'wrong'}}))
            with self.assertRaisesRegex(ValueError,'flow hash'):m.archive_run(run,root/'archive')
