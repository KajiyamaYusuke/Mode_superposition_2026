#!/usr/bin/env python3
"""Aggregate saved per-case analyses, leaving missing and ambiguous results blank."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'

def load():
    records=[]
    for r in csv.DictReader((ROOT/'prepared_cases.csv').open()):
        case=Path(r['parameter_file']).parent.parent; status=case/'status.json'; metric=case/'jitter_analysis/metrics.json'
        r['simulation_status']=json.loads(status.read_text()).get('simulation_status','UNKNOWN') if status.exists() else 'NOT_STARTED'
        if metric.exists():
            result=json.loads(metric.read_text());r['jitter_status']=result['jitter_status'];r['reason']=result['reason'];r['primary_method']=result.get('primary_method','');r['metrics']=result.get('flow_jitter',{});r['template']=result.get('template_candidates',{}).get('1',{});r['check']=result.get('maxima_check',{});r['periodicity']=result.get('flow_periodicity',{});r['candidates']=result.get('maxima_candidates',{})
        else:r.update(jitter_status='NOT_ANALYZED',reason='',primary_method='',metrics={},template={},check={},periodicity={},candidates={})
        records.append(r)
    return records

def main():
    records=load();summary=ROOT/'summary';summary.mkdir(exist_ok=False)
    fields=['case_id','ratio','pressure_pa','simulation_status','jitter_status','reason','primary_method','local','rap','ppq5','cv','cycles','acf_at_period','template_similarity_median']
    with (summary/'jitter_summary.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in records:
            sim=r['template'].get('similarity',[]);check=r['check'];w.writerow(dict(case_id=r['case_id'],ratio=r['ratio'],pressure_pa=r['pressure_pa'],simulation_status=r['simulation_status'],jitter_status=r['jitter_status'],reason=r['reason'],primary_method=r['primary_method'],cycles=check.get('cycles',''),acf_at_period=r['periodicity'].get('acf_at_period',''),template_similarity_median=np.median(sim) if sim else '',**{k:r['metrics'].get(k,'') for k in ('local','rap','ppq5','cv')}))
    with (summary/'jitter_method_comparison.csv').open('x',newline='') as f:
        fields=['case_id','ratio','pressure_pa','jitter_status','candidate','cycles','local','rap','ppq5','cv'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in records:
            sources=dict(r['candidates']);sources['template_T0']=r['template']
            for name,entry in sources.items():
                if not entry:continue
                metrics=entry.get('metrics',{});w.writerow(dict(case_id=r['case_id'],ratio=r['ratio'],pressure_pa=r['pressure_pa'],jitter_status=r['jitter_status'],candidate=name,cycles=entry.get('cycles',len(entry.get('periods',[]))),**{k:metrics.get(k,'') for k in ('local','rap','ppq5','cv')}))
    fig,ax=plt.subplots(figsize=(8,5));styles=[('900','tab:blue','o','-'),('1400','tab:red','s','--')]
    for pressure,color,marker,line in styles:
        rs=sorted((r for r in records if r['pressure_pa']==pressure),key=lambda r:float(r['ratio']))
        x=np.array([float(r['ratio']) for r in rs]);y=np.array([r['metrics'].get('local',np.nan) if r['jitter_status']=='VALID_STATIONARY' else np.nan for r in rs],float)
        ax.plot(x,y,color=color,marker=marker,linestyle=line,label=pressure+' Pa')
        for r in rs:
            if r['jitter_status']=='VALID_WINDOW_DEPENDENT' and 'local' in r['metrics']:ax.scatter(float(r['ratio']),r['metrics']['local'],marker='D',facecolors='none',edgecolors=color)
            elif r['jitter_status']=='VALID_STATIONARY' and r['primary_method']=='B_template' and 'local' in r['metrics']:ax.scatter(float(r['ratio']),r['metrics']['local'],marker='x',color=color,zorder=5)
    ax.set(xlabel="Body Young's modulus ratio E_B,R/E_B,L",ylabel='Flow Jitter local [%]',xlim=(.075,1.025));ax.legend();ax.grid(alpha=.25);fig.tight_layout();fig.savefig(summary/'jitter_vs_ratio_900_1400.png',dpi=180);fig.savefig(summary/'jitter_vs_ratio_900_1400.svg');plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for ax,(pressure,_,_,_) in zip(axs,styles):
        rs=sorted((r for r in records if r['pressure_pa']==pressure),key=lambda r:float(r['ratio']));x=[float(r['ratio']) for r in rs]
        for key in ('local','rap','ppq5','cv'):
            ax.plot(x,[r['metrics'].get(key,np.nan) for r in rs],marker='o',label=key)
        ax.set(ylabel='Percent',title=pressure+' Pa');ax.grid(alpha=.25);ax.legend(ncol=4)
    axs[-1].set_xlabel('Body modulus ratio');fig.tight_layout();fig.savefig(summary/'jitter_method_comparison.png',dpi=160);plt.close(fig)
    diff=[]
    for ratio in sorted(set(r['ratio'] for r in records),key=float):
        both=[r for r in records if r['ratio']==ratio]; pair={r['pressure_pa']:r for r in both}
        if len(pair)==2 and all('local' in pair[p]['metrics'] for p in ('900','1400')):diff.append((float(ratio),pair['1400']['metrics']['local']-pair['900']['metrics']['local']))
    fig,ax=plt.subplots();ax.axhline(0,color='gray',lw=.8)
    if diff:ax.scatter(*zip(*diff))
    ax.set(xlabel='Body modulus ratio',ylabel='J1400 − J900 [%]');fig.tight_layout();fig.savefig(summary/'pressure_difference.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for pressure,color,marker,line in styles:
        rs=sorted((r for r in records if r['pressure_pa']==pressure),key=lambda r:float(r['ratio']));x=[float(r['ratio']) for r in rs]
        acf=[r['periodicity'].get('acf_at_period',np.nan) for r in rs]
        similarity=[np.median(r['template']['similarity']) if r['template'].get('similarity') else np.nan for r in rs]
        axs[0].plot(x,acf,color=color,marker=marker,linestyle=line,label=pressure+' Pa')
        axs[1].plot(x,similarity,color=color,marker=marker,linestyle=line,label=pressure+' Pa')
    axs[0].set_ylabel('ACF at candidate period');axs[1].set_ylabel('Template similarity median');axs[1].set_xlabel('Body modulus ratio')
    for ax in axs:ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(summary/'periodicity_and_similarity.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for ax,(pressure,_,_,_) in zip(axs,styles):
        rs=sorted((r for r in records if r['pressure_pa']==pressure),key=lambda r:float(r['ratio']));x=[float(r['ratio']) for r in rs]
        for label,key in [('A maxima','A_0.10_1'),('B template','template')]:
            values=[]
            for r in rs:
                candidate=r['template'] if key=='template' else r['candidates'].get(key,{})
                value=candidate.get('metrics',{}).get('local')
                values.append(np.nan if value is None else value)
            ax.scatter(x,values,label=label)
        ax.set(ylabel='Candidate local [%]',title=pressure+' Pa — candidates include unvalidated estimates');ax.legend();ax.grid(alpha=.25)
    axs[-1].set_xlabel('Body modulus ratio');fig.tight_layout();fig.savefig(summary/'event_method_comparison.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(14,2.6));codes={'VALID_STATIONARY':'V','VALID_WINDOW_DEPENDENT':'W','UNDEFINED_JITTER':'?','NOT_ANALYZED':'—'}
    for iy,p in enumerate(('900','1400')):
        rs=sorted((r for r in records if r['pressure_pa']==p),key=lambda r:float(r['ratio']))
        for i,r in enumerate(rs):
            text=codes.get(r['jitter_status'],'!') if r['simulation_status']=='SIMULATION_COMPLETED' else ('…' if r['simulation_status'] in ('RUNNING','NOT_STARTED') else '×')
            ax.text(i,iy,text,ha='center',va='center',fontsize=13)
    ax.set(xticks=range(18),xticklabels=[f'{.1+i*.05:.2f}' for i in range(18)],yticks=[0,1],yticklabels=['900 Pa','1400 Pa'],title='QC: V valid, W window-dependent, ? undefined, … pending, × run failed');ax.set_xlim(-.5,17.5);ax.set_ylim(-.5,1.5);fig.tight_layout();fig.savefig(summary/'jitter_qc.png',dpi=160);plt.close(fig)
    examples=[]
    for pressure in ('900','1400'):
        subset=[r for r in records if r['pressure_pa']==pressure and r['simulation_status']=='SIMULATION_COMPLETED']
        valid=sorted((r for r in subset if r['metrics'].get('local') is not None),key=lambda r:r['metrics']['local'])
        chosen=([valid[0],valid[-1]] if len(valid)>1 else valid[:])
        uncertain=next((r for r in subset if r['jitter_status']=='UNDEFINED_JITTER'),None)
        if uncertain:chosen.append(uncertain)
        if not chosen and subset:chosen=[subset[0]]
        for r in chosen:
            if r not in examples:examples.append(r)
    if examples:
        fig,axs=plt.subplots(len(examples),1,figsize=(11,2.3*len(examples)),squeeze=False)
        for ax,r in zip(axs[:,0],examples):
            case=Path(r['parameter_file']).parent.parent;values=np.loadtxt(case/'output/dat/airflow_vt.dat',comments='#',ndmin=2);t=values[:,0]*1e-5;mask=(t>=.90)&(t<.95);ax.plot(t[mask],values[mask,1]*1e6,color='black',lw=.8)
            cycle_path=case/'jitter_analysis/flow_cycles.csv'
            if cycle_path.exists():
                cycles=list(csv.DictReader(cycle_path.open()))
                for name,marker in [('A_0.10_1','o'),('template_1','x')]:
                    ev=np.array([float(c['t_start_s']) for c in cycles if c['method']==name and .90<=float(c['t_start_s'])<.95])
                    if len(ev):ax.scatter(ev,np.interp(ev,t,values[:,1])*1e6,marker=marker,s=18,label=name)
                if ax.get_legend_handles_labels()[0]:ax.legend(loc='upper right')
            ax.set(ylabel='Ug [mL/s]',title=r['case_id']+' — '+r['jitter_status'])
        axs[-1,0].set_xlabel('Time [s]');fig.tight_layout();fig.savefig(summary/'flow_peak_examples.png',dpi=150);plt.close(fig)
    print('Saved',len(records),'case status records to',summary)
    if len(records)==36:
        import write_report
        write_report.main()

if __name__=='__main__':main()
