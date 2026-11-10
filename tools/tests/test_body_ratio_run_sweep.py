import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
p=Path(__file__).resolve().parents[1]/'body_ratio_jitter_261103/run_sweep.py';spec=importlib.util.spec_from_file_location('run_sweep',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class ValidateTests(unittest.TestCase):
    def test_complete_and_missing_samples(self):
        with tempfile.TemporaryDirectory() as t:
            case=Path(t);root=case/'output';root.mkdir()
            (case/'stdout.log').write_text('  ps      = 900 [Pa]\n')
            for d in ('dat','txt','wav'):(root/d).mkdir()
            steps=np.arange(0,100000,5); flow=np.sin(steps/100)*.0001
            np.savetxt(root/'dat/airflow_vt.dat',np.column_stack((steps,flow)),header='step airflow')
            for name in ('displace.dat','displace_xy.dat','pressure_vt.dat'):
                axis=steps if name=='pressure_vt.dat' else steps*1e-5
                np.savetxt(root/'dat'/name,np.column_stack((axis,np.zeros(len(steps)))),header='step value')
            (root/'dat/area.dat').write_text('header\n');(root/'wav/test_sound.wav').write_bytes(b'notempty');(root/'txt/params_used.txt').write_text('param')
            (root/'txt/manifest.txt').write_text('nstep = 100000\ndt_s = 1e-05\nnmode = 20\npressure_schedule_final_pa = 900\nwrite_deformation_vtu = 0\n')
            self.assertEqual(m.validate(root,'900')['flow_samples'],20000)
            np.savetxt(root/'dat/airflow_vt.dat',np.column_stack((steps[:-1],flow[:-1])),header='step airflow')
            with self.assertRaises(ValueError):m.validate(root,'900')
    def test_no_overwrite_existing_case(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);(root/'input').mkdir(); param=root/'input/param.txt';param.write_text('')
            (root/'status.json').write_text('{}')
            with self.assertRaises(FileExistsError):m.run_case({'parameter_file':str(param)},Path('/bin/true'),4)
