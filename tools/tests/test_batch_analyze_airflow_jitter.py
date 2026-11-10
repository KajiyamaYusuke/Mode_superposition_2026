"""Batch wrapper tests; opt in to the two preserved native-data regressions."""
import contextlib
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import batch_analyze_airflow_jitter as batch

np = batch.np


def read_csv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_manifest(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)
    return path


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='airflow_batch_test_')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def flow(self, name='run', missing=False):
        folder = self.root / name / 'snapshot/dat'
        folder.mkdir(parents=True)
        t = np.arange(2020) / 2000
        phase = np.mod(t, .01)
        x = np.maximum(1 - abs(phase - .004) / .002, 0)
        x[abs(x) < 1e-12] = 0
        if missing:
            x[(t >= .25) & (t < .35)] = np.maximum(x[(t >= .25) & (t < .35)], .08)
        np.savetxt(folder / 'airflow_vt.dat', np.c_[np.arange(len(t))*50, x], header='step airflow[m^3/s]')
        return folder / 'airflow_vt.dat'

    def manifest(self, rows):
        return write_manifest(self.root / 'cases.csv', rows)

    def execute(self, *arguments):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return batch.main(list(map(str, arguments)))

    def test_multiple_paths_and_nonascii_spaces_no_metadata_inference(self):
        a = self.flow('非ASCII 空白 r0p25_P900')
        b = self.flow('second')
        out = self.root / 'results'
        self.assertEqual(self.execute('--paths', a.parent.parent.parent, b, '--out-dir', out), 0)
        rows = read_csv(out / 'summary.csv')
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r['stiffness_ratio'] == r['pressure_pa'] == 'NA' for r in rows))
        self.assertNotEqual(rows[0]['case_id'], rows[1]['case_id'])
        self.assertEqual(read_csv(out / 'inputs_resolved.csv')[0]['flow_path'], str(a))
        self.assertFalse((out / 'figures').exists())
        self.assertEqual(read_csv(out / 'failures.csv'), [])
        self.assertFalse(list(out.rglob('*.json')))

    def test_manifest_single_multiple_bad_ids_and_both_source_columns(self):
        flow = self.flow()
        for identities in (['dup', 'dup'], ['../escape'], ['has space'], [''], ['a..b']):
            manifest = self.manifest([dict(case_id=i, flow_path=str(flow)) for i in identities])
            out = self.root / 'invalid'
            self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 2)
            self.assertFalse(out.exists())
        manifest = self.manifest([dict(case_id='both', flow_path=str(flow), run_path=str(flow.parent))])
        out = self.root / 'both'
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 1)
        self.assertEqual(read_csv(out / 'summary.csv')[0]['status'], 'FAILED')
        out = self.root / 'single'
        manifest = self.manifest([dict(case_id='single', flow_path=str(flow))])
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 0)
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 2)

    def test_dry_run_different_cwd_and_repo_relative_paths(self):
        flow = self.flow('日本語 space')
        manifest = self.manifest([dict(case_id='c', flow_path=os.path.relpath(flow, batch.REPO_ROOT))])
        out = self.root / 'not_written'
        with patch.object(batch.core, 'acf_guided_extract', side_effect=AssertionError('detector must not run')):
            self.assertEqual(self.execute('--manifest', manifest, '--dry-run', '--out-dir', out), 0)
        self.assertFalse(out.exists())
        call = subprocess.run([sys.executable, '-B', str(batch.REPO_ROOT/'tools/batch_analyze_airflow_jitter.py'),
                               '--manifest', str(manifest), '--dry-run'], cwd=self.root, capture_output=True, text=True)
        self.assertEqual(call.returncode, 0, call.stderr)
        self.assertIn(str(flow), call.stdout)

    def test_order_missing_file_and_fail_fast(self):
        flow = self.flow()
        manifest = self.manifest([dict(case_id='first', flow_path=str(flow)),
                                  dict(case_id='bad', flow_path=str(self.root/'missing/airflow_vt.dat')),
                                  dict(case_id='last', flow_path=str(flow))])
        out = self.root / 'continue'
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 2)
        rows = read_csv(out / 'summary.csv')
        self.assertEqual([r['case_id'] for r in rows], ['first', 'bad', 'last'])
        self.assertEqual([r['status'] for r in rows], ['COMPLETED', 'FAILED', 'COMPLETED'])
        self.assertEqual(read_csv(out/'failures.csv')[0]['phase'], 'PREFLIGHT')
        out = self.root / 'stop'
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out, '--fail-fast'), 2)
        self.assertEqual([r['status'] for r in read_csv(out/'summary.csv')], ['COMPLETED', 'FAILED', 'NOT_RUN_FAIL_FAST'])
        self.assertFalse((out/'cases/last').exists())

    def test_runtime_failure_and_interrupt_leave_no_partial_case(self):
        flow = self.flow()
        manifest = self.manifest([dict(case_id='bad', flow_path=str(flow)), dict(case_id='good', flow_path=str(flow))])
        real = batch.analyze_case
        def run(case, settings):
            if case['case_id'] == 'bad':
                raise IndexError('intentional unexpected case failure')
            return real(case, settings)
        out = self.root/'runtime'
        with patch.object(batch, 'analyze_case', side_effect=run):
            self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 2)
        self.assertEqual([r['status'] for r in read_csv(out/'summary.csv')], ['FAILED', 'COMPLETED'])
        self.assertFalse((out/'cases/bad').exists())
        self.assertEqual(list((out/'cases').glob('.*')), [])
        out = self.root/'interrupt'
        with patch.object(batch, 'analyze_case', side_effect=KeyboardInterrupt):
            self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 130)
        self.assertEqual([r['status'] for r in read_csv(out/'summary.csv')], ['INTERRUPTED', 'NOT_RUN_INTERRUPTED'])
        self.assertEqual(list((out/'cases').iterdir()), [])

    def test_format_steps_window_and_ambiguous_run_failures(self):
        flow = self.flow()
        a = self.flow('ambiguous/a')
        self.flow('ambiguous/b')
        manifest = self.manifest([dict(case_id='ambiguous', run_path=str(a.parents[3]))])
        self.assertEqual(self.execute('--manifest', manifest, '--dry-run'), 1)
        np.savetxt(flow, [[0, 0], [100, 1], [100, 0]])
        self.assertEqual(self.execute('--paths', flow, '--dry-run'), 1)
        np.savetxt(flow, [[0, 0], [100, 1], [200, 0]])
        self.assertEqual(self.execute('--paths', flow, '--dry-run'), 1)
        np.savetxt(flow, [[0, 0], [100, 1], [300, 0]])
        self.assertEqual(self.execute('--paths', flow, '--dry-run', '--start', 0, '--end', .001), 1)
        with self.assertRaises(SystemExit):
            batch.parser().parse_args(['--manifest', 'x', '--paths', 'y'])

    def test_all_intervals_and_terms_independently_recomputed(self):
        flow = self.flow(missing=True)
        out = self.root/'all_intervals'
        self.assertEqual(self.execute('--paths', flow, '--out-dir', out, '--start', .05, '--end', .95), 0)
        row = read_csv(out/'summary.csv')[0]
        times = [float(e['event_time_s']) for e in read_csv(out/'events_3pct.csv')]
        v = np.diff(times)
        self.assertGreater(max(v), .08)
        self.assertGreater(int(row['n_without_3pct_event']), 0)
        pp = read_csv(out/'periods_3pct.csv')
        np.testing.assert_array_equal(v, [float(p['interval_s']) for p in pp])
        self.assertEqual(sum(bool(p['jitter_contribution_pct']) for p in pp), len(v)-1)
        self.assertEqual(pp[-1]['next_T_s'], '')
        self.assertAlmostEqual(float(row['jitter_local_pct']), 100*np.mean(abs(np.diff(v)))/np.mean(v))
        terms = read_csv(out/'metric_terms_3pct.csv')
        for width, key, summary_key in ((2,'jitter_local_percent','jitter_local_pct'), (3,'rap_percent','rap_pct'), (5,'ppq5_percent','ppq5_pct')):
            nums = [abs(v[j+1]-v[j]) if width == 2 else abs(v[j+width//2]-sum(v[j:j+width])/width) for j in range(len(v)-width+1)]
            selected = [t for t in terms if t['metric'] == key]
            self.assertEqual(len(selected), len(v)-width+1)
            np.testing.assert_allclose([float(t['numerator_s']) for t in selected], nums, rtol=1e-14, atol=1e-16)
            self.assertAlmostEqual(float(row[summary_key]), 100*np.mean(nums)/np.mean(v))

    @staticmethod
    def selector_fixture():
        x=np.array([0,1,.02,0,.6,.02,0,1.3,.02,0,0,1,.1,0,0],float)
        data=dict(x=x,t=np.arange(len(x))*.001,fs=1000.,steps=np.arange(len(x))*100,dt_s=1e-5,stride=100)
        regions=[dict(region_id=0,start_s=0.,end_s=.01,search_scale_s=.01,signal='Ug'),
                 dict(region_id=1,start_s=.01,end_s=.015,search_scale_s=.005,signal='Ug')]
        pp=[]
        for i,prom,passed in ((1,.2,True),(4,.5,True),(7,.05,False),(11,1.,True)):
            pp.append(dict(event_id='Ug_max_'+str(i),sample_id=i,raw_time_s=float(data['t'][i]),event_time_s=float(data['t'][i]),
                sample_value=float(x[i]),prominence=prom,search_region_id=0 if i<10 else 1,method='Ug_max',inside_window=True,
                clarity_accepted=passed,detector_rejection_reasons=[] if passed else ['BELOW_LOCAL_PROMINENCE']))
        series=batch.core.all_event_extract(data,None,0,.015)[0];series['Ug_max']=[p for p in pp if p['clarity_accepted']]
        return data,dict(series=series,regions=regions,peak_candidates=pp)

    def test_peak_selectors_ties_raw_rejections_and_observed_override(self):
        data,base=self.selector_fixture()
        runs={s:batch.core.peak_associated_extract(data,0,.015,base,flow_threshold_pct=5,peak_selector=s) for s in batch.core.PEAK_SELECTORS}
        self.assertEqual([next(p for p in runs[s]['peaks'] if p['peak_role']=='MAIN')['sample_id'] for s in batch.core.PEAK_SELECTORS],[4,1,7])
        self.assertEqual([runs[s]['series']['main_zero'][0]['sample_id'] for s in batch.core.PEAK_SELECTORS],[6,3,9])
        raw_main=next(p for p in runs['raw-flow']['peaks'] if p['peak_role']=='MAIN')
        self.assertFalse(raw_main['clarity_accepted']);self.assertEqual(raw_main['detector_rejection_reasons'],['BELOW_LOCAL_PROMINENCE'])
        solo=batch.core.peak_associated_extract(data,0,.015,base,flow_threshold_pct=5,main_peak_overrides={0:'Ug_max_7'})
        self.assertEqual([e['event_time_s'] for e in solo['series']['main_3pct']],[e['event_time_s'] for e in runs['raw-flow']['series']['main_3pct']])
        for mode in ('clear-flow','raw-flow'):
            tied=[dict(p,sample_value=1.) for p in base['series']['Ug_max']]
            chosen,_=batch.core.peak_associated_select(tied,base['regions'],mode)
            self.assertEqual(next(p for p in chosen if p['peak_role']=='MAIN')['sample_id'],1)
        with self.assertRaises(ValueError):batch.core.peak_associated_extract(data,0,.015,base,main_peak_overrides={0:'NONEXISTENT'})
        with self.assertRaises(ValueError):batch.core.peak_associated_select(base['series']['Ug_max'],base['regions'],'invalid')

    def test_selector_empty_region_and_previous_branch_recomputed(self):
        import copy
        data,base=self.selector_fixture();empty=copy.deepcopy(base)
        empty['series']['Ug_max']=[p for p in empty['series']['Ug_max'] if p['search_region_id']==1]
        a=batch.core.peak_associated_extract(data,0,.015,empty)
        c=batch.core.peak_associated_extract(data,0,.015,empty,peak_selector='raw-flow')
        self.assertIsNone(a['cycles'][0]['main_peak_id']);self.assertEqual(c['cycles'][0]['main_peak_id'],'Ug_max_7')
        altered=copy.deepcopy(base);altered['regions']=[dict(region_id=i,start_s=l,end_s=r,search_scale_s=r-l,signal='Ug') for i,(l,r) in enumerate(((0,.004),(.004,.010),(.010,.015)))]
        for p in altered['peak_candidates']:p['search_region_id']=batch.core.acf_region_owner(p['event_time_s'],altered['regions'])
        altered['series']['Ug_max']=altered['peak_candidates']
        a=batch.core.peak_associated_extract(data,0,.015,altered)
        b=batch.core.peak_associated_extract(data,0,.015,altered,peak_selector='clear-flow')
        self.assertEqual(a['branches']['Ug_max_1']['branch_end_s'],.004)
        self.assertEqual(b['branches']['Ug_max_1']['branch_end_s'],.007)
        self.assertNotEqual(a['branches']['Ug_max_4']['selected']['3pct']['event_time_s'],b['branches']['Ug_max_7']['selected']['3pct']['event_time_s'])

    def test_time_first_event_comparison_records_reassignment(self):
        def e(at,rid,peak):return dict(event_time_s=at,search_region_id=rid,peak_id=peak,event_id=peak+'_event')
        old=[e(.51,0,'a'),e(.52,1,'b'),e(.54,3,'c')]
        new=[e(.51,1,'d'),e(.525,2,'e'),e(.541,3,'c')]
        rows=batch.selection_event_differences(old,new)
        self.assertEqual({r['change'] for r in rows},{'SAME_TIME_REASSIGNED','MISSING','ADDED','TIME_CHANGED_IN_REGION'})
        self.assertEqual(sum(r['A_event_id'] is not None for r in rows),len(old))
        self.assertEqual(sum(r['variant_event_id'] is not None for r in rows),len(new))

    def test_sensitivity_cli_atomic_outputs_baseline_and_full_terms(self):
        flow=self.flow();manifest=self.manifest([dict(case_id='synthetic',flow_path=str(flow),stiffness_ratio=.45,pressure_pa=900)])
        refs=self.root/'baseline'
        for pct,name in ((3,'allRatio_001'),(5,'allRatio_5pct_001'),(10,'allRatio_10pct_001')):
            self.assertEqual(self.execute('--manifest',manifest,'--flow-threshold-pct',pct,'--out-dir',refs/name),0)
        before=batch.core.sha256(refs/'allRatio_001/events_3pct.csv')
        out=self.root/'comparison'
        self.assertEqual(self.execute('--manifest',manifest,'--compare-peak-selectors','--comparison-thresholds',3,5,10,'--baseline-root',refs,'--out-dir',out),0)
        rows=read_csv(out/'comparison.csv');self.assertEqual(len(rows),9)
        self.assertEqual({r['selector'] for r in rows},{'A','B','C'})
        self.assertTrue(all(r['ids_parent_ids_times_intervals_metrics']=='EXACT' for r in read_csv(out/'regression_checks.csv')))
        events=read_csv(out/'events.csv');terms=read_csv(out/'metric_terms.csv')
        for r in rows:
            es=[e for e in events if e['selector']==r['selector'] and e['flow_threshold_pct']==r['flow_threshold_pct']]
            v=np.diff([float(e['event_time_s']) for e in es]);self.assertEqual(len(v),int(r['n_intervals']))
            tt=[t for t in terms if t['selector']==r['selector'] and t['flow_threshold_pct']==r['flow_threshold_pct']]
            self.assertEqual([sum(t['metric']==k for t in tt) for k in batch.core.METRIC_NAMES],[max(0,len(v)-1),max(0,len(v)-2),max(0,len(v)-4)])
        self.assertEqual(batch.core.sha256(refs/'allRatio_001/events_3pct.csv'),before)
        self.assertEqual(len(list((out/'figures').glob('*.png'))),3)
        self.assertEqual(self.execute('--manifest',manifest,'--compare-peak-selectors','--baseline-root',refs,'--out-dir',out),2)

    def test_flow_threshold_default_explicit_fractional_and_metadata(self):
        flow = self.flow(missing=True)
        baseline = self.root/'default_pct'
        explicit = self.root/'explicit_pct'
        self.assertEqual(self.execute('--paths', flow, '--start', .05, '--end', .95, '--out-dir', baseline), 0)
        self.assertEqual(self.execute('--paths', flow, '--start', .05, '--end', .95, '--flow-threshold-pct', 3, '--out-dir', explicit), 0)
        self.assertEqual(read_csv(baseline/'events_3pct.csv'), read_csv(explicit/'events_3pct.csv'))
        self.assertEqual(read_csv(baseline/'periods_3pct.csv'), read_csv(explicit/'periods_3pct.csv'))
        reference_peaks = read_csv(baseline/'peaks.csv')
        for pct, suffix in ((5.,'5pct'), (2.5,'2p5pct'), (10.,'10pct')):
            out = self.root/suffix
            self.assertEqual(self.execute('--paths', flow, '--start', .05, '--end', .95, '--flow-threshold-pct', pct, '--out-dir', out), 0)
            row = read_csv(out/'summary.csv')[0]
            self.assertEqual(float(row['flow_threshold_pct']), pct)
            self.assertEqual(row['event_method'], 'peak-associated-closure/main_'+suffix)
            self.assertIn('n_without_'+suffix+'_event', row)
            self.assertNotIn('n_without_3pct_event', row)
            self.assertEqual(read_csv(out/'peaks.csv'), reference_peaks)
            self.assertFalse((out/'events_3pct.csv').exists())
            self.assertFalse((out/'periods_3pct.csv').exists())
            self.assertEqual(float(read_csv(out/'inputs_resolved.csv')[0]['flow_threshold_pct']), pct)
            events = read_csv(out/('events_'+suffix+'.csv'))
            for e in events:
                self.assertEqual(float(e['threshold_pct']), pct)
                level = float(e['Upeak_m3_s'])*pct/100
                self.assertEqual(float(e['threshold_m3_s']), level)
                self.assertGreater(float(e['value_left']), level)
                self.assertGreaterEqual(level, float(e['value_right']))
                self.assertTrue(e['event_id'].startswith('main_'+suffix+'_'))
            times = np.array([float(e['event_time_s']) for e in events])
            periods = np.diff(times)
            self.assertAlmostEqual(float(row['jitter_local_pct']),100*np.mean(abs(np.diff(periods)))/np.mean(periods))
            self.assertIn('events_'+suffix+'.csv', (out/'report.md').read_text())
            self.assertIn(batch.threshold_label(pct)+'% airflow', (out/'report.md').read_text())
            area = read_csv(out/'area_alignment.csv')
            self.assertIn('has_'+suffix+'_event', area[0])
            self.assertIn('Amin_at_'+suffix+'_mm2', area[0])
            self.assertNotIn('has_3pct_event', area[0])
            if pct == 10:
                self.assertGreater(int(row['n_events']), len(read_csv(baseline/'events_3pct.csv')))

    def test_flow_threshold_invalid_and_comparison_group_separation(self):
        for value in ('0','100','-1','nan','inf','-inf','not-a-number'):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
                batch.parser().parse_args(['--paths','unused','--flow-threshold-pct',value])
            self.assertEqual(exc.exception.code,2)
        base = dict(settings_id='same',start_s=.5,end_s=1.,pressure_pa=900.,stiffness_ratio=.3)
        groups=batch.comparison_groups([dict(base,case_id='3',flow_threshold_pct=3.),dict(base,case_id='5',flow_threshold_pct=5.)])
        self.assertEqual(len(groups),2)
        self.assertEqual(batch.core.flow_threshold_suffix(3.0),'3pct')
        self.assertEqual(batch.core.flow_threshold_suffix(.5),'0p5pct')

    def test_fractional_threshold_area_remains_independent(self):
        flow=self.flow()
        args=batch.parser().parse_args(['--paths',str(flow),'--flow-threshold-pct','2.5'])
        def analyze():
            cases=batch.load_cases(args)
            with contextlib.redirect_stdout(io.StringIO()):
                batch.preflight(cases)
                return batch.analyze_case(cases[0],'test')
        missing=analyze()
        steps=np.arange(0,2020,2)*50
        sections=np.full((len(steps),50),10.)
        sections[:,1]=2+np.sin(np.arange(len(steps))*.2)
        np.savetxt(flow.with_name('area.dat'),np.c_[steps,sections],header='step area[mm^2]')
        aligned=analyze()
        self.assertEqual(missing['events_2p5pct'],aligned['events_2p5pct'])
        self.assertEqual(missing['summary'][0]['jitter_local_pct'],aligned['summary'][0]['jitter_local_pct'])
        for row in aligned['area_alignment']:
            if row['has_2p5pct_event']:
                self.assertAlmostEqual(row['lag_ms'],1000*(row['t_Amin_min_s']-row['t_Ug_2p5pct_s']))
                self.assertIsNotNone(row['Amin_at_2p5pct_mm2'])

    def test_area_missing_valid_invalid_does_not_change_events(self):
        flow = self.flow()
        args = batch.parser().parse_args(['--paths', str(flow)])
        def analyze():
            cases = batch.load_cases(args)
            with contextlib.redirect_stdout(io.StringIO()):
                batch.preflight(cases)
                return batch.analyze_case(cases[0], 'test')
        missing = analyze()
        self.assertEqual(missing['summary'][0]['area_status'], 'AREA_NOT_AVAILABLE')
        steps = np.arange(0, 2020, 2)*50
        sections = np.full((len(steps),50), 10.)
        sections[:, 0] = sections[:, -1] = 0
        sections[:,1] = 2 + np.sin(np.arange(len(steps))*.2)
        sections[:,2] = 2 + np.cos(np.arange(len(steps))*.2)
        sections[:,10] = np.nan
        np.savetxt(flow.with_name('area.dat'), np.c_[steps, sections], header='step area[mm^2]')
        valid = analyze()
        self.assertEqual(valid['summary'][0]['area_status'], 'AVAILABLE')
        self.assertGreater(valid['summary'][0]['n_area_aligned_events'], 0)
        self.assertTrue(any(r['min_section_switch_count'] for r in valid['area_alignment']))
        np.savetxt(flow.with_name('area.dat'), np.c_[steps, sections])
        invalid = analyze()
        self.assertEqual(invalid['summary'][0]['area_status'], 'AREA_FAILED')
        for other in (valid, invalid):
            self.assertEqual(missing['events_3pct'], other['events_3pct'])
            self.assertEqual(missing['summary'][0]['jitter_local_pct'], other['summary'][0]['jitter_local_pct'])

    def test_area_lag_native_stride_nan_section_plateau_and_missing_event(self):
        flowpath = self.flow()
        np.savetxt(flowpath, np.c_[np.arange(11)*100, np.ones(11)], header='step airflow')
        sections = np.full((6,50), 10.)
        sections[:,0] = sections[:,-1] = 0
        sections[:,10] = np.nan
        sections[1,1] = 2
        sections[2,2] = 1
        sections[3:5,1] = 0
        np.savetxt(flowpath.with_name('area.dat'), np.c_[np.arange(6)*200,sections], header='step area[mm^2]')
        case = dict(case_id='a', resolved_flow=flowpath, resolved_area=flowpath.with_name('area.dat'), dt_s=1e-5,
                    start_s=0., end_s=.01, area_status='AVAILABLE', area_reason='', area_sha256=batch.core.sha256(flowpath.with_name('area.dat')))
        flow = batch.read_flow(case)
        result = dict(peaks=[dict(peak_role='MAIN', event_id='p1', raw_time_s=.001, search_region_id=0),
                            dict(peak_role='MAIN', event_id='p2', raw_time_s=.006, search_region_id=1)],
                      cycles=[dict(main_peak_id='p1', right_censored=False),dict(main_peak_id='p2',right_censored=True)],
                      branches={'p1':dict(branch_end_s=.006),'p2':dict(branch_end_s=.01)},
                      series={'main_3pct':[dict(peak_id='p1', event_time_s=.003)]})
        rows = batch.area_alignment(case, flow, result)
        self.assertEqual(rows[0]['min_section_index'], 2)
        self.assertEqual(rows[0]['min_section_column_zero_based'], 3)
        self.assertAlmostEqual(rows[0]['t_Amin_min_s'], .004)
        self.assertAlmostEqual(rows[0]['lag_ms'], 1.)
        self.assertAlmostEqual(rows[0]['Amin_at_3pct_mm2'], 1.5)
        self.assertEqual(rows[0]['amin_at_event_rule'], 'LINEAR_NATIVE_AMIN_BRACKET')
        self.assertIsNone(rows[1]['t_Ug_3pct_s'])
        self.assertNotIn('lag_ms', rows[1])
        self.assertEqual(rows[1]['Amin_min_mm2'], 0)
        self.assertEqual(rows[1]['min_plateau_samples'], 2)
        self.assertTrue(rows[1]['right_censored'])
        self.assertGreater(rows[1]['rows_with_undefined_section'], 0)
        result['series']['main_3pct'][0]['event_time_s'] = .005
        negative = batch.area_alignment(case, flow, result)
        self.assertAlmostEqual(negative[0]['lag_ms'], -1.)
        # Declared m² gives the same mm² observations; no unit guessing.
        np.savetxt(flowpath.with_name('area.dat'), np.c_[np.arange(6)*200,sections*1e-6], header='step area[m^2]')
        case['area_sha256'] = batch.core.sha256(flowpath.with_name('area.dat'))
        converted = batch.area_alignment(case, flow, result)
        self.assertAlmostEqual(converted[0]['Amin_min_mm2'], 1.)
        # A final event beyond the last saved area sample is never extrapolated.
        np.savetxt(flowpath.with_name('area.dat'), np.c_[np.arange(5)*200,sections[:5]], header='step area[mm^2]')
        case['area_sha256'] = batch.core.sha256(flowpath.with_name('area.dat'))
        result['series']['main_3pct'].append(dict(peak_id='p2',event_time_s=.009))
        truncated = batch.area_alignment(case, flow, result)
        self.assertIsNone(truncated[1]['Amin_at_3pct_mm2'])
        self.assertEqual(truncated[1]['amin_at_event_rule'], 'OUTSIDE_SAMPLED_RANGE_NO_EXTRAPOLATION')

    def test_area_other_run_rejected_without_flow_failure(self):
        flow = self.flow()
        wrong = self.root/'other/area.dat'
        manifest = self.manifest([dict(case_id='a', flow_path=str(flow), area_path=str(wrong))])
        out = self.root/'area_wrong'
        self.assertEqual(self.execute('--manifest', manifest, '--out-dir', out), 0)
        row = read_csv(out/'summary.csv')[0]
        self.assertEqual(row['status'], 'COMPLETED')
        self.assertEqual(row['area_status'], 'AREA_FAILED')
        self.assertIn('SAME dat directory', row['area_reason'])

    def test_insufficient_events_blank_not_zero_and_all_failure_exit(self):
        flow = self.flow()
        data = np.loadtxt(flow)
        data[:,1] = 1
        np.savetxt(flow, data)
        out = self.root/'no_events'
        self.assertEqual(self.execute('--paths', flow, '--out-dir', out), 0)
        row = read_csv(out/'summary.csv')[0]
        self.assertEqual(row['status'], 'INSUFFICIENT_EVENTS')
        self.assertEqual([row[k] for k in ('jitter_local_pct','rap_pct','ppq5_pct')], ['','',''])
        self.assertEqual(read_csv(out/'events_3pct.csv'), [])
        self.assertEqual(read_csv(out/'periods_3pct.csv'), [])
        self.assertEqual(self.execute('--paths', self.root/'missing', '--out-dir', self.root/'all_failed'), 1)

    def test_comparison_groups_separate_pressure_window_unknowns_and_output_guard(self):
        base = dict(settings_id='fixed', status='COMPLETED', start_s=.5, end_s=1., pressure_pa=900., stiffness_ratio=.3)
        rows = [dict(base,case_id='a'), dict(base,case_id='b',stiffness_ratio=.4),
                dict(base,case_id='c',pressure_pa=1400), dict(base,case_id='d',start_s=.6),
                dict(base,case_id='unknown1',pressure_pa=None), dict(base,case_id='unknown2',pressure_pa=None),
                dict(base,case_id='no_ratio',stiffness_ratio=None)]
        groups = batch.comparison_groups(rows)
        self.assertEqual(sorted(map(len, groups.values())), [1,1,1,1,2])
        flow = self.flow()
        self.assertEqual(self.execute('--paths', flow, '--out-dir', flow.parent/'result'), 2)
        archive_out = batch.REPO_ROOT/'analysis_runs/manual_run_archive/should_never_be_created_by_batch'
        self.assertEqual(self.execute('--paths', flow, '--out-dir', archive_out), 2)
        self.assertFalse(archive_out.exists())


@unittest.skipUnless(os.environ.get('JITTER_BATCH_REAL_REGRESSION') == '1', 'set JITTER_BATCH_REAL_REGRESSION=1 for preserved 0.25/0.30 inputs')
class NativeRegression(unittest.TestCase):
    def test_batch_single_and_preserved_full_precision_events_metrics(self):
        with tempfile.TemporaryDirectory(prefix='airflow_batch_native_') as directory:
            out = Path(os.environ['JITTER_BATCH_TEST_OUTPUT']) if os.environ.get('JITTER_BATCH_TEST_OUTPUT') else Path(directory)/'batch'
            with contextlib.redirect_stdout(io.StringIO()):
                code = batch.main(['--manifest', 'configs/jitter_batch_cases.example.csv', '--out-dir', str(out)])
            self.assertEqual(code, 0)
            rows = read_csv(out/'summary.csv')
            events = read_csv(out/'events_3pct.csv')
            periods = read_csv(out/'periods_3pct.csv')
            reference = batch.REPO_ROOT/'analysis_runs/jitter_all_events_20261107/peak_associated_closure_20261107/results'
            expected = {'r0p25':(61,60,21.92513631206122), 'r0p30':(58,57,16.98616957590101)}
            for row in rows:
                case = row['case_id']
                prior = next(r for r in read_csv(reference/case/'comparison.csv') if r['method']=='main_3pct')
                self.assertEqual((int(row['n_events']),int(row['n_intervals'])), expected[case][:2])
                self.assertEqual(float(row['jitter_local_pct']), expected[case][2])
                for field,key in (('jitter_local_pct','jitter_local_percent'),('rap_pct','rap_percent'),('ppq5_pct','ppq5_percent')):
                    self.assertEqual(float(row[field]), float(prior[key]))
                current_events = [e for e in events if e['case_id']==case]
                prior_events = [e for e in read_csv(reference/case/'events.csv') if e['method']=='main_3pct']
                self.assertEqual([(e['event_id'],e['parent_peak_id'],float(e['event_time_s'])) for e in current_events],
                                 [(e['event_id'],e['peak_id'],float(e['event_time_s'])) for e in prior_events])
                self.assertEqual([float(p['interval_s']) for p in periods if p['case_id']==case],
                                 [float(p['interval_s']) for p in read_csv(reference/case/'periods.csv') if p['method']=='main_3pct'])
                # One-case invocation uses the same exact vector, not only the displayed metric.
                one = Path(directory)/case
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(batch.main(['--paths',read_csv(out/'inputs_resolved.csv')[rows.index(row)]['flow_path'], '--out-dir',str(one)]),0)
                self.assertEqual([(e['event_id'],float(e['event_time_s'])) for e in read_csv(one/'events_3pct.csv')],
                                 [(e['event_id'],float(e['event_time_s'])) for e in current_events])
                self.assertEqual(float(read_csv(one/'summary.csv')[0]['jitter_local_pct']),float(row['jitter_local_pct']))
                single_args = batch.core.parser().parse_args([
                    '--input', read_csv(out/'inputs_resolved.csv')[rows.index(row)]['flow_path'],
                    '--dt', '1e-5', '--start', '.5', '--end', '1', '--case-id', case,
                    '--tracking-method', 'peak-associated-closure', '--out-dir', str(Path(directory)/('single_'+case))])
                with patch.object(batch.core, 'peak_associated_figures'), contextlib.redirect_stdout(io.StringIO()):
                    single = batch.core.peak_associated_cli(single_args)
                single_events = single['result']['series']['main_3pct']
                self.assertEqual([(e['event_id'],e['peak_id'],e['event_time_s']) for e in single_events],
                                 [(e['event_id'],e['parent_peak_id'],float(e['event_time_s'])) for e in current_events])
                single_metric = next(r for r in single['rows'] if r['method']=='main_3pct')
                for field,key in (('jitter_local_pct','jitter_local_percent'),('rap_pct','rap_percent'),('ppq5_pct','ppq5_percent')):
                    self.assertEqual(float(row[field]), single_metric[key])
            self.assertLessEqual(len(list(out.rglob('*.png'))),1)
            self.assertFalse(list(out.rglob('*.json')))


if __name__ == '__main__':
    unittest.main()
