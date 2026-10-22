import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/kanaya_ablation'))
from analyze_r0_residual_aperture import interval_metrics,sample_profile


def frame(rows):
    return pd.DataFrame([dict(z0_mm=a,z1_mm=b,h0_mm=h0,h1_mm=h1,
        yL0_mm=0.,yL1_mm=0.,yR0_mm=h0,yR1_mm=h1) for a,b,h0,h1 in rows])


class ApertureTests(unittest.TestCase):
    def test_exact_zero_crossing(self):
        result=interval_metrics(frame([(0.,1.,-1.,1.)]),-np.inf,np.inf)
        self.assertAlmostEqual(result['area_mm2'],.25)
        self.assertAlmostEqual(result['open_z_width_mm'],.5)
        self.assertAlmostEqual(result['negative_gap_z_width_mm'],.5)

    def test_holes_are_not_bridged(self):
        f=frame([(0.,1.,1.,1.),(2.,3.,2.,2.)])
        self.assertTrue(np.isnan(sample_profile(f,np.array([1.5]))).all())
        stats=interval_metrics(f,0.,3.)
        self.assertEqual(stats['area_mm2'],3.)
        self.assertEqual(stats['valid_z_width_mm'],2.)

    def test_region_partition_is_conservative(self):
        f=frame([(0.,1.,-1.,1.),(1.,3.,1.,2.)])
        total=interval_metrics(f,0.,3.)['area_mm2']
        parts=[interval_metrics(f,a,b)['area_mm2']for a,b in [(0.,.5),(.5,2.5),(2.5,3.)]]
        self.assertAlmostEqual(total,sum(parts))


if __name__=='__main__':unittest.main()
