#!/usr/bin/env python3
"""Inspect completed representative runs and make waveform QC plots."""
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import analyze_flow_jitter as flow
import run_sweep

ROOT=Path(__file__).resolve().parents[2]/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'

def main():
    out=ROOT/'representative_qc';out.mkdir(exist_ok=False)
    records=[]
    for pressure in ('900','1400'):
        case=ROOT/'cases'/('r0p50_P'+pressure);status=json.loads((case/'status.json').read_text())
        if status['simulation_status']!='SIMULATION_COMPLETED': raise ValueError('representative did not complete')
        validation=run_sweep.validate(case/'output',pressure)
        if validation['flow_sha256']!=status['validation']['flow_sha256']:raise ValueError('flow hash changed')
        t,x=flow.load_flow(case/'output/dat/airflow_vt.dat'); sample=x[(t>=.5)&(t<1.)]; jump=np.diff(x); idx=np.argsort(np.abs(jump))[-5:][::-1]
        info={'pressure_pa':int(pressure),'validation':validation,'mean_flow_m3_s':float(np.mean(sample)),'flow_p05_m3_s':float(np.quantile(sample,.05)),'flow_p95_m3_s':float(np.quantile(sample,.95)),'maximum_jumps':[{'time_s':float(t[i]),'jump_m3_s':float(jump[i])} for i in idx]}
        analysis=flow.analyze(x);info['preliminary_peak_status']=analysis['status'];info['preliminary_peak_reason']=analysis.get('reason','');records.append(info)
        fig,ax=plt.subplots(2,1,figsize=(12,7));ax[0].plot(t,x*1e6,lw=.5);ax[0].axvspan(.5,1,color='orange',alpha=.12);ax[0].set(ylabel='Ug [mL/s]',title='r=0.50, '+pressure+' Pa')
        mask=(t>=.90)&(t<.95);ax[1].plot(t[mask],x[mask]*1e6,lw=.8)
        ref=x[(t>=.5)&(t<1)]
        for method,marker in [('A','o'),('B','v')]:
            e=flow.events(x,dt=5e-5,method=method,reference=ref);e=e[(e>=.90)&(e<.95)];ax[1].scatter(e,np.interp(e,t,x)*1e6,s=18,marker=marker,label=method)
        ax[1].set(xlabel='Time [s]',ylabel='Ug [mL/s]');ax[1].legend();fig.tight_layout();fig.savefig(out/('flow_P'+pressure+'.png'),dpi=150);plt.close(fig)
    (out/'summary.json').write_text(json.dumps(records,indent=2)+'\n');print('representative QC',out)

if __name__=='__main__':main()
