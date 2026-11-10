import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd


PATH = Path(__file__).parents[1] / "amplitude_asymmetry" / "postprocess_contact_work_bitdelta.py"
SPEC = importlib.util.spec_from_file_location("postprocess_contact_work_bitdelta", PATH)
M = importlib.util.module_from_spec(SPEC); assert SPEC.loader; SPEC.loader.exec_module(M)


class PostprocessContactWorkBitdeltaTest(unittest.TestCase):
    def test_dat_exact_has_no_difference(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"same.dat"; path.write_text("# step value\n0 1\n")
            rows,summary=M.compare_dat_file(path,path,1e-5)
            self.assertEqual(rows,[]); self.assertEqual(summary["different_cells"],0)

    def test_dat_exact_text_only_and_numeric_delta(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); old=root/"old.dat"; new=root/"new.dat"
            old.write_text("# step value\n0 1.0\n5 2.0\n")
            new.write_text("# step value\n0 1.00\n5 2.1\n")
            rows,summary=M.compare_dat_file(old,new,1e-5)
            self.assertEqual([row["difference_kind"] for row in rows],["TEXT_ONLY","NUMERIC_DELTA"])
            self.assertEqual(summary["first_step"],0); self.assertEqual(summary["numeric_delta_cells"],1)

    def test_modal_one_ulp_reports_location_and_handles_signed_zero(self):
        old=np.zeros((2,2,2,7)); new=old.copy()
        new[1,1,0,4]=np.nextafter(0.0,1.0); new[0,0,0,0]=-0.0
        rows,summary,_=M.modal_delta_rows(old,new,.1)
        self.assertEqual(summary["first_different_step"],1)
        hit=next(row for row in rows if row["field"]=="q_n" and row["side"]=="R" and row["mode"]==1)
        self.assertEqual(hit["first_different_step"],1); self.assertEqual(hit["first_ulp_distance"],1)

    def test_work_split_windows_and_region_availability(self):
        old_windows=M.WINDOWS; M.WINDOWS=((0,.2),)
        try:
            modal=np.zeros((2,2,2,7)); component=np.zeros((2,2,2,2)); modal[...,6]=1
            component[:,0,:,0]=[[1,-.5],[2,-1]]; component[:,0,:,1]=.25
            component[:,1,:,0]=-.5; component[:,1,:,1]=-.25
            modal[...,3]=component.sum(axis=-1)
            rows=[]
            for step in range(2):
                for side in M.SIDES:
                    for ap in ("end","interior"):
                        rows.append({"step":step,"time_s":step*.1,"side":side,"ap_region":ap,"pairs":1,"max_penetration_m":1e-4,"spring_work_nodal_J":0.0,"damp_work_nodal_J":0.0})
            _,windows,region_mode,checks=M.work_products(modal,component,.1,pd.DataFrame(rows))
            self.assertEqual(checks["max_abs_component_sum_minus_reference_J"],0)
            self.assertTrue(any(row["availability"].endswith("NOT_IDENTIFIABLE") for row in region_mode))
            self.assertAlmostEqual(next(row for row in windows if row["side"]=="BOTH")["W_total_J"],-0.5)
        finally: M.WINDOWS=old_windows

    def test_window_left_closed_right_open(self):
        np.testing.assert_array_equal(np.flatnonzero(M.window_mask(4,.1,.1,.3)),[1,2])

    def test_ulp_nan_and_negative(self):
        self.assertIsNone(M.ulp_distance(float("nan"),1.0))
        self.assertEqual(M.ulp_distance(-1.0,np.nextafter(-1.0,-2.0)),1)

    def test_file_identity_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/"x"; p.write_bytes(b"abc"); self.assertEqual(M.file_identity(p),M.file_identity(p))

    def test_existing_output_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileExistsError): M.require_fresh_output(Path(directory))

    def test_missing_ap_rows_do_not_become_zero_work(self):
        old_windows=M.WINDOWS; M.WINDOWS=((0,.1),)
        try:
            modal=np.zeros((1,2,1,7)); component=np.zeros((1,2,1,2))
            empty=pd.DataFrame(columns=["step","time_s","side","ap_region","pairs","max_penetration_m","spring_work_nodal_J","damp_work_nodal_J"])
            _,_,rows,_=M.work_products(modal,component,.1,empty)
            self.assertFalse(any(row["classification"]=="AP_REGION" for row in rows))
        finally: M.WINDOWS=old_windows


if __name__ == "__main__": unittest.main()
