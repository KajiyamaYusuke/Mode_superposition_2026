#!/usr/bin/env python3
"""Analyze one completed flow file; preserve uncertainty and original output."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import analyze_flow_jitter as flow
import template_period as template
import run_sweep


def assess_template(estimate, tracks, evidence, selection):
    """Conservative B validation against a second template scale and spectrum."""
    one=tracks.get('1',{});half=tracks.get('0.5',{})
    if estimate.get('status')!='CANDIDATE' or 'times' not in one or 'times' not in half:
        return {'status':'UNDEFINED_JITTER','reason':'template or period ambiguous'}
    t=np.asarray(one['times']); q=np.asarray(one['similarity']);p=np.diff(t); other=np.asarray(half['times']); po=np.diff(other)
    if len(p)<20 or len(po)<20:return {'status':'UNDEFINED_JITTER','reason':'insufficient template cycles'}
    med=float(np.median(p)); expected=estimate['period_s']; match=float(np.mean(np.min(np.abs(t[:,None]-other[None,:]),axis=1)<.2*expected))
    m=flow.period_metrics(p); mh=flow.period_metrics(po)
    threshold=max(.2,.3*max(m['local'],mh['local']))
    if (abs(med/expected-1)>.05 or np.any((p<.5*med)|(p>1.5*med)) or
        np.quantile(q,.05)<.8 or match<.9 or abs(m['local']-mh['local'])>threshold or
        not np.isfinite(evidence['acf_at_period']) or evidence['acf_at_period']<.5 or
        not np.isfinite(evidence['spectral_ratio']) or evidence['spectral_ratio']<.05):
        return {'status':'UNDEFINED_JITTER','reason':'template quality, period or scale sensitivity'}
    halves=[]
    for start,end in ((.5,.75),(.75,1.)):
        periods=flow.window_periods(t,start,end)
        if len(periods)<10:return {'status':'UNDEFINED_JITTER','reason':'insufficient half-window cycles'}
        halves.append((float(np.mean(periods)),flow.period_metrics(periods)))
    amp=[]
    for part in np.array_split(selection,2):amp.append(float(np.percentile(part,95)-np.percentile(part,5)))
    j1,j2=halves[0][1]['local'],halves[1][1]['local']
    nonstationary=(abs(halves[0][0]-halves[1][0])/med>.05 or abs(amp[0]-amp[1])/max(amp)>.2 or abs(j1-j2)>max(.2,.5*max(j1,j2)))
    return {'status':'VALID_WINDOW_DEPENDENT' if nonstationary else 'VALID_STATIONARY', 'reason':'', 'metrics':m, 'event_match_fraction':match, 'similarity_p05':float(np.quantile(q,.05)), 'half_metrics':[halves[0][1],halves[1][1]]}


def plain(value):
    if isinstance(value,dict): return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [plain(x) for x in value]
    if isinstance(value,np.ndarray): return plain(value.tolist())
    if isinstance(value,(float,np.floating)) and not np.isfinite(value): return None
    if isinstance(value,(np.integer,)): return int(value)
    if isinstance(value,(np.floating,)): return float(value)
    return value


def analyze_case(case):
    status=json.loads((case/'status.json').read_text())
    if status['simulation_status']!='SIMULATION_COMPLETED': raise ValueError('case has not completed')
    if run_sweep.sha(case/'output/dat/airflow_vt.dat') != status.get('validation',{}).get('flow_sha256'):
        raise ValueError('flow hash does not match completed run')
    t,x=flow.load_flow(case/'output/dat/airflow_vt.dat');dt=5e-5; window=(t>=.5)&(t<1); selection=x[window]
    result=flow.analyze(x); report={'signal':'Ug','window_s':[.5,1.0],'run_status':'SIMULATION_COMPLETED','jitter_status':'UNDEFINED_JITTER','reason':'not validated'}
    report['analysis_code_sha256']={name:run_sweep.sha(Path(__file__).resolve().parent/name) for name in ('analyze_saved.py','analyze_flow_jitter.py','template_period.py')}
    config=run_sweep.STUDY/'analysis_config.json'
    if config.is_file(): report['analysis_config_sha256']=run_sweep.sha(config)
    candidates=result.pop('candidates',{});report['maxima_check']=result
    try:
        estimate=template.estimate_period(selection,dt);report['template_period_estimate']=estimate
        report['flow_periodicity']=flow.periodicity(selection,dt,estimate['period_s'])
        tracks={}
        for multiplier in (.5,1,2):
            try:
                z=template.track(selection,dt,estimate['period_s']*multiplier); times=z['times']+.5; periods=np.diff(times)
                tracks[str(multiplier)]={'times':times,'similarity':z['similarity'],'periods':periods,'metrics':flow.period_metrics(periods)}
            except ValueError as exc: tracks[str(multiplier)]={'error':str(exc)}
        report['template_candidates']=tracks
        b=tracks['1'];a=candidates.get('A_0.10_1');
        if estimate['status']=='CANDIDATE' and a and 'periods' in b and len(b['periods'])>=10 and result.get('status') in ('VALID_STATIONARY','VALID_WINDOW_DEPENDENT'):
            pa=a['periods'];pb=b['periods'];
            if abs(np.mean(pa)/np.mean(pb)-1)<.05 and abs(a['metrics']['local']-b['metrics']['local'])<max(.2,.25*max(a['metrics']['local'],b['metrics']['local'])):
                report['jitter_status']=result['status'];report['reason']='';report['flow_jitter']=a['metrics'];report['primary_method']='A_local_maximum'
            else: report['reason']='A/template method difference'
        else: report['reason']='event extraction or fundamental-period candidate ambiguous'
        bcheck=assess_template(estimate,tracks,report['flow_periodicity'],selection)
        report['template_validation']=bcheck
        if report['jitter_status']=='UNDEFINED_JITTER' and bcheck['status'] in ('VALID_STATIONARY','VALID_WINDOW_DEPENDENT'):
            report['jitter_status']=bcheck['status'];report['reason']='';report['flow_jitter']=bcheck['metrics'];report['primary_method']='B_template'
    except ValueError as exc: report['reason']=str(exc); tracks={}
    # Supplementary pressure evidence, never treated as Ug jitter.
    try:
        p=np.loadtxt(case/'output/dat/pressure_vt.dat',comments='#',ndmin=2)
        pressure=p[:,1][window]
        report['supplementary_acoustic_periodicity']=template.estimate_period(pressure,dt)
        if 'template_period_estimate' in report:
            ratio=report['supplementary_acoustic_periodicity']['period_s']/report['template_period_estimate']['period_s']
            report['supplementary_acoustic_periodicity']['period_ratio_to_Ug']=ratio
            report['supplementary_acoustic_periodicity']['period_candidate_agrees']=bool(abs(ratio-1)<.1)
            if abs(ratio-1)<.1 and 'times' in tracks.get('1',{}):
                aux=template.track(pressure,dt,report['template_period_estimate']['period_s'])['times']+.5
                report['supplementary_acoustic_periodicity']['one_aux_event_per_Ug_cycle_fraction']=flow.event_correspondence(tracks['1']['times'],aux,.5,1.)
    except (OSError,ValueError,IndexError) as exc: report['supplementary_acoustic_periodicity']={'status':'unavailable','reason':str(exc)}
    try:
        area=np.loadtxt(case/'output/dat/area.dat',comments='#',ndmin=2)
        if len(area)!=len(t) or area.shape[1]<4: raise ValueError('area samples unavailable')
        # The endpoints are boundary cuts; use only finite interior cuts.
        interior=area[:,2:-1]; valid=np.isfinite(interior); minimum=np.min(np.where(valid,interior,np.inf),axis=1)
        if not np.isfinite(minimum[window]).all(): raise ValueError('some area sections have no finite cut')
        report['supplementary_area_periodicity']=template.estimate_period(minimum[window],dt)
        if 'template_period_estimate' in report:
            ratio=report['supplementary_area_periodicity']['period_s']/report['template_period_estimate']['period_s']
            report['supplementary_area_periodicity']['period_ratio_to_Ug']=ratio
            report['supplementary_area_periodicity']['period_candidate_agrees']=bool(abs(ratio-1)<.1)
            if abs(ratio-1)<.1 and 'times' in tracks.get('1',{}):
                aux=template.track(minimum[window],dt,report['template_period_estimate']['period_s'])['times']+.5
                report['supplementary_area_periodicity']['one_aux_event_per_Ug_cycle_fraction']=flow.event_correspondence(tracks['1']['times'],aux,.5,1.)
    except (OSError,ValueError,IndexError) as exc: report['supplementary_area_periodicity']={'status':'unavailable','reason':str(exc)}
    out=case/'jitter_analysis'
    out.mkdir(exist_ok=False)
    with (out/'flow_peaks.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['method','time_s','flow_m3_s','similarity'])
        for key,val in candidates.items():
            for at in val['times']:
                if .5<=at<1:w.writerow([key,at,float(np.interp(at,t,x)),''])
        for key,val in tracks.items():
            if 'times' not in val:continue
            for at,similarity in zip(val['times'],val['similarity']):
                w.writerow(['template_'+key,at,float(np.interp(at,t,x)),similarity])
    with (out/'jitter_sensitivity.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['method','cycles','local_percent','rap_percent','ppq5_percent','cv_percent','validity'])
        for key,val in candidates.items():
            q=val['metrics'];w.writerow([key,len(val['periods']),q['local'],q['rap'],q['ppq5'],q['cv'],report['jitter_status']])
    with (out/'window_comparison.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['signal','method','start_s','end_s','cycles','mean_period_s','local_percent','rap_percent','ppq5_percent','cv_percent','validity'])
        for name,start,end in [('reference',.3,.5),('main',.5,1.),('first_half',.5,.75),('second_half',.75,1.)]:
            entry=candidates.get('A_0.10_1')
            if entry is None: continue
            periods=flow.window_periods(entry['times'],start,end);metrics=flow.period_metrics(periods)
            # Candidate values are diagnostics; a window is not certified separately.
            w.writerow(['Ug',name,start,end,len(periods),np.mean(periods) if len(periods) else '',metrics['local'],metrics['rap'],metrics['ppq5'],metrics['cv'],'CANDIDATE_ONLY'])
    with (out/'flow_cycles.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['method','t_start_s','t_end_s','period_s','window'])
        for key,val in candidates.items():
            ts=val['times'];ts=ts[(ts>=.5)&(ts<1)]
            for a,b in zip(ts,ts[1:]):w.writerow([key,a,b,b-a,'[0.50,1.00)'])
        for key,val in tracks.items():
            if 'times' not in val:continue
            for a,b in zip(val['times'],val['times'][1:]):w.writerow(['template_'+key,a,b,b-a,'[0.50,1.00)'])
    report['maxima_candidates']={k:{'cycles':len(v['periods']),'metrics':v['metrics']} for k,v in candidates.items()}
    (out/'metrics.json').write_text(json.dumps(plain(report),indent=2,allow_nan=False)+'\n')
    print(case.name,report['jitter_status'],report['reason'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);a=p.parse_args();analyze_case(a.case.resolve())
