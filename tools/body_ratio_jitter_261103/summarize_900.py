#!/usr/bin/env python3
"""Create the 900 Pa report from complete, separately validated cases."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT=Path(__file__).resolve().parents[2]/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'


def main():
    rows=[r for r in csv.DictReader((ROOT/'prepared_cases.csv').open()) if r['pressure_pa']=='900']
    rows.sort(key=lambda r:float(r['ratio']))
    if len(rows)!=18:raise ValueError('expected 18 900 Pa cases')
    reports=[]
    for row in rows:
        case=Path(row['parameter_file']).parent.parent
        status=json.loads((case/'status.json').read_text()); analysis=json.loads((case/'jitter_analysis/metrics.json').read_text())
        if status['simulation_status']!='SIMULATION_COMPLETED':raise ValueError('case not complete '+row['case_id'])
        reports.append((row,case,status,analysis))
    out=ROOT/'analysis_900_complete';out.mkdir(exist_ok=False)
    fields=['ratio','case_id','simulation_status','jitter_status','reason','primary_method','local_percent','rap_percent','ppq5_percent','cv_percent','flow_sha256']
    with (out/'jitter_summary_900.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row,case,status,a in reports:
            metric=a.get('flow_jitter',{});w.writerow(dict(ratio=row['ratio'],case_id=row['case_id'],simulation_status=status['simulation_status'],jitter_status=a['jitter_status'],reason=a.get('reason',''),primary_method=a.get('primary_method',''),local_percent=metric.get('local',''),rap_percent=metric.get('rap',''),ppq5_percent=metric.get('ppq5',''),cv_percent=metric.get('cv',''),flow_sha256=status['validation']['flow_sha256']))
    x=np.array([float(r['ratio']) for r,_,_,_ in reports]);y=np.array([a.get('flow_jitter',{}).get('local',np.nan) for _,_,_,a in reports],float)
    fig,(ax,qc)=plt.subplots(2,1,figsize=(10,6),sharex=True,gridspec_kw={'height_ratios':[5,1]});ax.plot(x,y,color='#155fad',lw=1.4)
    for ratio,_,_,a in reports:
        val=a.get('flow_jitter',{}).get('local')
        if val is None:continue
        mark='D' if a['jitter_status']=='VALID_WINDOW_DEPENDENT' else ('x' if a.get('primary_method')=='B_template' else 'o')
        ax.scatter(float(ratio['ratio']),val,marker=mark,color='#155fad',s=50)
    ax.set(ylabel='Flow jitter local [%]',title='900 Pa: flow jitter versus body modulus ratio');ax.grid(alpha=.2)
    ax.legend(handles=[Line2D([0],[0],marker='o',color='#155fad',linestyle='',label='A: flow maximum'),Line2D([0],[0],marker='x',color='#155fad',linestyle='',label='B: template'),Line2D([0],[0],marker='D',color='#155fad',linestyle='',label='Window dependent')],fontsize=8)
    for row,_,_,a in reports:
        marker='o' if a['jitter_status']=='VALID_STATIONARY' else ('D' if a['jitter_status']=='VALID_WINDOW_DEPENDENT' else 'x')
        qc.scatter(float(row['ratio']),0,marker=marker,color='#155fad' if marker!='x' else '#888888')
    qc.set(yticks=[0],yticklabels=['900 Pa'],xlabel="Body Young's modulus ratio $E_{B,R}/E_{B,L}$",ylim=(-1,1));qc.set_xticks(x);qc.tick_params(axis='x',rotation=45,labelsize=8)
    fig.text(.12,.01,'Ug, [0.50, 1.00) s; × in lower strip = undefined. Template thresholds are provisional.',fontsize=9);fig.subplots_adjust(bottom=.2,hspace=.08)
    for ext in ('png','svg'):fig.savefig(out/('jitter_vs_ratio_900.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for key,label in [('local','Local'),('rap','RAP'),('ppq5','PPQ5'),('cv','Period CV')]:
        values=[a.get('flow_jitter',{}).get(key,np.nan) for _,_,_,a in reports]
        ax.plot(x,values,marker='o',label=label)
    ax.set(xlabel='Body modulus ratio',ylabel='Percent',title='900 Pa: flow-period metrics');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(out/'period_metrics_900.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
    acf=[a.get('flow_periodicity',{}).get('acf_at_period',np.nan) for _,_,_,a in reports]
    similarity=[]
    for _,_,_,a in reports:
        vals=a.get('template_candidates',{}).get('1',{}).get('similarity',[])
        similarity.append(float(np.median(vals)) if vals else np.nan)
    axes[0].plot(x,acf,'o-');axes[0].set_ylabel('ACF at candidate period');axes[1].plot(x,similarity,'o-');axes[1].set(ylabel='Median template similarity',xlabel='Body modulus ratio')
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'periodicity_900.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for method,label in [('A_0.10_1','A: maxima'),('template','B: template')]:
        values=[]
        for _,_,_,a in reports:
            source=a.get('maxima_candidates',{}).get(method,{}) if method!='template' else a.get('template_candidates',{}).get('1',{})
            candidate=source.get('metrics',{}).get('local');values.append(np.nan if candidate is None else candidate)
        ax.scatter(x,values,label=label)
    ax.set(xlabel='Body modulus ratio',ylabel='Candidate local jitter [%]',title='900 Pa: event methods (includes unvalidated candidates)');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(out/'event_methods_900.png',dpi=180);plt.close(fig)
    eligible=[item for item in reports if item[3].get('flow_jitter',{}).get('local') is not None]
    eligible.sort(key=lambda item:item[3]['flow_jitter']['local'])
    examples=([eligible[0],eligible[-1]] if len(eligible)>1 else eligible[:])
    uncertain=next((item for item in reports if item[3]['jitter_status']=='UNDEFINED_JITTER'),None)
    if uncertain and uncertain not in examples:examples.append(uncertain)
    if not examples:examples=[reports[0]]
    fig,axes=plt.subplots(len(examples),1,figsize=(10,2.5*len(examples)),squeeze=False)
    for ax,(row,case,_,a) in zip(axes[:,0],examples):
        data=np.loadtxt(case/'output/dat/airflow_vt.dat',comments='#',ndmin=2);t=data[:,0]*1e-5;mask=(t>=.9)&(t<.95);ax.plot(t[mask],data[mask,1]*1e6,color='black',lw=.8)
        events=list(csv.DictReader((case/'jitter_analysis/flow_peaks.csv').open()))
        for method,marker in [('A_0.10_1','o'),('template_1','x')]:
            times=np.array([float(e['time_s']) for e in events if e['method']==method and .9<=float(e['time_s'])<.95])
            if len(times):ax.scatter(times,np.interp(times,t,data[:,1])*1e6,marker=marker,s=20,label=method)
        ax.set(ylabel='Ug [mL/s]',title=row['case_id']+' — '+a['jitter_status'])
        if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    axes[-1,0].set_xlabel('Time [s]');fig.tight_layout();fig.savefig(out/'flow_examples_900.png',dpi=180);plt.close(fig)
    valid=[(r,a) for r,_,_,a in reports if 'local' in a.get('flow_jitter',{})]
    lines=['# 900 Pa sweep','',f'- Simulations complete: {len(reports)}/18.',f'- Jitter numerically available: {len(valid)}/18; other results are undefined and were not replaced with zero.','- Flow Ug; main analysis window [0.50,1.00) s; original T3 d2 surface; VTU output disabled.','- Jitter may use method A or B; method and QC are recorded in CSV. Template thresholds were set after examining the first 900 Pa case and need sensitivity and independent validation.']
    if valid:
        valid.sort(key=lambda z:z[1]['flow_jitter']['local']);lo,hi=valid[0],valid[-1]
        lines.append(f"- Smallest available local jitter: {lo[1]['flow_jitter']['local']:.4g}% at ratio {lo[0]['ratio']}; largest: {hi[1]['flow_jitter']['local']:.4g}% at ratio {hi[0]['ratio']}.")
    lines+=['','High jitter alone does not prove chaos or a particular physical cause. Check candidate event times, waveforms and per-window comparisons before interpreting trends.','','Each case retains input hashes, raw flow, other solver outputs, logs, and analysis files.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')
    print('Saved 900 Pa results to',out)

if __name__=='__main__':main()
