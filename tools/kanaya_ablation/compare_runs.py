#!/usr/bin/env python3
"""Comparable metrics, symmetry/work gates and plots for completed ablation runs."""
import argparse
import json
import os
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
from scipy import signal

from run_matrix import ROOT, OUT, dump, sha


def numeric_gate(a,b,atol=1e-12,rtol=1e-8):
    a,b=np.asarray(a),np.asarray(b)
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    delta=np.abs(a-b)
    scaled=delta/(atol+rtol*np.maximum(np.abs(a),np.abs(b)))
    assert scaled.max(initial=0)<=1, f'numerical gate failed: {scaled.max()}'
    return {'max_abs':float(delta.max(initial=0)),'max_scaled':float(scaled.max(initial=0))}


def closure_periods(area,dt):
    closed=area<=.01
    starts=np.flatnonzero(np.diff(closed.astype(int))==1)+1
    periods=np.diff(starts)*dt
    stable=len(periods)>=3 and np.std(periods)/np.mean(periods)<=.15 and closed.any() and (~closed).any()
    return {'closure_episodes':int(len(starts)), 'closure_period_CV':float(np.std(periods)/np.mean(periods)) if len(periods)>0 else np.nan,
            'periodic_closure':bool(stable), 'CQ':float(np.mean(closed)) if stable else np.nan,
            'f0_from_closures_Hz':float(1/np.mean(periods)) if stable else np.nan}


def spectrum(values,dt):
    values=np.asarray(values)
    ft=np.fft.rfft(signal.detrend(values)*np.hanning(len(values)))
    power=np.abs(ft)**2; power[0]=0
    if power.max()==0: return np.nan,0.
    peak=int(np.argmax(power))
    return float(np.fft.rfftfreq(len(values),dt)[peak]),float(power[max(1,peak-1):peak+2].sum()/power.sum())


def phase_at_frequency(force,velocity,dt,frequency):
    if not np.isfinite(frequency) or np.std(force)==0 or np.std(velocity)==0:
        return np.nan
    k=int(np.argmin(np.abs(np.fft.rfftfreq(len(force),dt)-frequency)))
    taper=np.hanning(len(force))
    a=np.fft.rfft(signal.detrend(force)*taper)[k]
    b=np.fft.rfft(signal.detrend(velocity)*taper)[k]
    return float(np.angle(a*np.conj(b),deg=True)) if abs(a)*abs(b)>0 else np.nan


def analyze(case):
    directory=OUT/case; base=directory/'output/kanaya_diagnostic'
    record=json.loads((directory/'manifest.json').read_text())
    meta=json.loads((base/'run_metadata.json').read_text())
    assert record['exit_code']==0 and meta['complete'] and meta['nmode']==[20,20]
    assert meta['nstep']==20000 and meta['integrated_intervals']==19999 and meta['excluded_intervals']==1
    n,dt=meta['integrated_intervals'],meta['dt_s']
    assert (base/'modal_native.bin').stat().st_size==n*2*20*7*8
    raw=np.memmap(base/'modal_native.bin',mode='r',dtype='<f8',shape=(n,2,20,7))
    flow=pd.read_csv(base/'flow_geometry.csv',float_precision='round_trip')
    summary=pd.read_csv(base/'modal_energy_summary.csv',float_precision='round_trip')
    consistency=pd.read_csv(base/'consistency.csv')
    included=consistency[consistency.included==1]
    assert len(included)==2*n and included.valid.eq(1).all() and included.scaled_max.max()<=1
    assert np.array_equal(flow.step,np.arange(20000))
    assert flow.iloc[:n].interval_included.eq(1).all() and flow.iloc[n].interval_included==0
    gates={'consistency_scaled_max':float(included.scaled_max.max())}
    # Reconstruct the ACTUAL NR equation from independently logged states.
    # The auxiliary pressure accumulator sums (U_(n-1)-Ud0_(n-1)); Pd0 sums
    # dt/Ca*(U_n-Ud0_(n-1)). With zero initial states their difference
    # telescopes: dt/Ca*Pu_aux_n = Pd0_n - dt/Ca*U_n. This also covers a
    # vocal tract without confusing the post-update Pd0 with the NR input.
    p=record['effective_parameters'];rho=float(p['rho_kg_m3'])
    lu=rho*float(p['subglottal_length_m'])/int(p['subglottal_sections'])/(2*np.pi*float(p['subglottal_radius_m'])**2)
    # New S3 manifests provide the actual initialized Lu (legacy has no /2).
    # Historical manifests retain the exact original fallback calculation.
    if 'acoustic_Lu' in p:lu=float(p['acoustic_Lu'])
    length=float(p['vocal_tract_length_m'])
    if length==0.:
        la,ca=.1,1e30
    else:
        av=np.pi*float(p['vocal_tract_radius_m'])**2
        dx=length/int(p['vocal_tract_sections'])
        la=rho*dx/av;ca=dx*av/(rho*float(p['c_sound_m_s'])**2)
    u=flow.Ug_m3_s.to_numpy();previous=np.r_[0.,u[:-1]]
    residual=flow.Rk*np.abs(u)*u+flow.Rv*u+(flow.Lg+la+lu)*(u-previous)/dt+flow.Pd0_Pa-dt/ca*u-flow.Pg_Pa
    active=flow.coefficients_active==1
    max_nr=float(residual[active].abs().max())
    assert max_nr<1e-8, f'coupled NR residual {max_nr} Pa'
    assert (u[~active]==0.).all()
    gates['coupled_NR_residual_max_Pa']=max_nr
    gates['NR_reconstruction_tolerance_Pa']=1e-8
    for k,field in [(0,'Qpos'),(1,'Qneg'),(2,'Qfluid'),(3,'Qcontact'),(4,'q'),(5,'qdot'),(6,'dq')]:
        gates['mirror_'+field]=numeric_gate(raw[:,0,:,k],raw[:,1,:,k])
    gates['mirror_uy']=numeric_gate(flow.uyL_mm.iloc[:n],-flow.uyR_mm.iloc[:n],atol=1e-9)
    gates['mirror_ux']=numeric_gate(flow.uxL_mm.iloc[:n],flow.uxR_mm.iloc[:n],atol=1e-9)
    gates['Q_sign_split']=numeric_gate(raw[...,0]+raw[...,1],raw[...,2])
    forces=pd.read_csv(base/'force_regions.csv')
    contacts=pd.read_csv(base/'contact_regions.csv')
    for frame,label in [(forces,'fluid'),(contacts,'contact')]:
        table=frame.groupby(['step','side'])[['Fx_N','Fy_N','opening_Fy_N']].sum().unstack('side')
        gates[label+'_mirror_Fx']=numeric_gate(table['Fx_N']['L'],table['Fx_N']['R'])
        gates[label+'_mirror_Fy']=numeric_gate(table['Fy_N']['L'],-table['Fy_N']['R'])
        gates[label+'_mirror_opening']=numeric_gate(table['opening_Fy_N']['L'],table['opening_Fy_N']['R'])
    assert contacts.step.max()==n-1 and len(contacts)==4*n
    contact_steps=contacts[contacts.side=='L'].groupby('step').agg(
        pairs=('pairs','sum'),max_penetration_m=('max_penetration_m','max'),
        opening_force_N=('opening_Fy_N','sum'))
    frequency=summary[(summary.side=='L')&(summary.window=='all')].sort_values('mode').f_Hz.to_numpy()
    omega=2*np.pi*frequency
    zeta=np.array([float(record['effective_parameters']['zetaL']),float(record['effective_parameters']['zetaR'])])
    work=[]
    for s,side in enumerate(['L','R']):
        for window,a,b in [('all',0,n),('ramp',0,5000),('transient',5000,10000),
                            ('early',10000,14999),('late',15000,n),('evaluation',10000,n)]:
            data=raw[a:b,s]
            force,q,v,dq=[data[...,k] for k in [2,4,5,6]]
            total=np.sum(force*dq,axis=0)
            mean=force.mean(axis=0)*dq.sum(axis=0)
            fluct=np.sum((force-force.mean(axis=0))*dq,axis=0)
            numeric_gate(total,mean+fluct,atol=1e-18,rtol=1e-10)
            wp=np.sum(data[...,0]*dq,axis=0); wn=np.sum(data[...,1]*dq,axis=0)
            contact=np.sum(data[...,3]*dq,axis=0)
            damping=np.sum(2*zeta[s]*omega*v*v*dt,axis=0)
            dac=np.sum(2*zeta[s]*omega*(v-v.mean(axis=0))**2*dt,axis=0)
            numeric_gate(total,wp+wn,atol=1e-18,rtol=1e-10)
            if window in ['all','ramp','transient','evaluation']:
                target=summary[(summary.side==side)&(summary.window==window)].sort_values('mode')
                for measured,key in [(total,'W_fluid'),(wp,'W_pos'),(wn,'W_neg'),(contact,'W_contact_reference'),(damping,'D_struct')]:
                    numeric_gate(measured,target[key],atol=1e-18,rtol=1e-9)
            for m in range(20):
                corr=float(np.corrcoef(force[:,m],v[:,m])[0,1]) if np.std(force[:,m])*np.std(v[:,m])>0 else np.nan
                peak,share=spectrum(q[:,m],dt)
                phase=phase_at_frequency(force[:,m],v[:,m],dt,peak)
                work.append({'case':case,'side':side,'mode':m+1,'window':window,'duration_s':(b-a)*dt,
                    'Wtotal':total[m],'Wmean':mean[m],'Wfluct':fluct[m],'Wpos':wp[m],'Wneg':wn[m],
                    'Wcontact_reference':contact[m],'D_total':damping[m],'D_ac_reference':dac[m],
                    'mean_q':q[:,m].mean(),'mean_Q':force[:,m].mean(),'q_AC_RMS':np.std(q[:,m]),
                    'v_AC_RMS':np.std(v[:,m]),'Q_v_correlation':corr,'q_peak_Hz':peak,'q_peak_power_share':share,
                    'phase_Q_ahead_of_v_deg_at_q_peak':phase})
    work=pd.DataFrame(work)
    rows=[]
    for label,a,b in [('early',10000,14999),('late',15000,n),('evaluation',10000,n)]:
        f=flow.iloc[a:b]
        c=contact_steps.loc[a:b-1]
        area=f.Amin_mm2.to_numpy(); assert np.isfinite(area).all()
        w=work[(work.side=='L')&(work.window==label)]
        force_stats=[r for r in meta['force_statistics'] if r['side']=='L' and r['window']==3]
        pos=next(r for r in force_stats if r['pressure_sign']==1)
        neg=next(r for r in force_stats if r['pressure_sign']==-1)
        peak,share=spectrum(f.Ug_m3_s,dt)
        row={'case':case,'window':label,'count':b-a,'duration_s':(b-a)*dt,
            'Amin_mean_mm2':area.mean(),'Amin_AC_RMS_mm2':np.std(area),'Amin_min_mm2':area.min(),
            'closure_fraction':float(np.mean(area<=.01)), **closure_periods(area,dt),
            'uyL_mean_mm':f.uyL_mm.mean(),'uyL_AC_RMS_mm':np.std(f.uyL_mm),
            'uxL_mean_mm':f.uxL_mm.mean(),'uxL_AC_RMS_mm':np.std(f.uxL_mm),
            'Ug_mean_ml_s':f.Ug_m3_s.mean()*1e6,'Ug_AC_RMS_ml_s':np.std(f.Ug_m3_s)*1e6,
            'Ug_peak_Hz_not_f0':peak,'Ug_peak_power_share':share,
            'Pg_mean_Pa':f.Pg_Pa.mean(),'Psep_mean_Pa':f.Psep_Pa.mean(),
            'xSep_mean_mm':f.xSep_mm.mean(),'xSep_min_mm':f.xSep_mm.min(),'xSep_max_mm':f.xSep_mm.max(),
            'Lg_mean':f.Lg.mean(),'Rv_mean':f.Rv.mean(),'Rk_mean':f.Rk.mean(),
            'solver_contact_flag_fraction':f.contact_flag.mean(),
            'adopted_contact_interval_fraction':float((c.pairs>0).mean()),
            'contact_max_penetration_m':c.max_penetration_m.max(),
            'contact_mean_opening_force_left_N':c.opening_force_N.mean(),
            'contact_iterations_mean':f.contact_iterations.mean(),
            'Wtotal_left':w.Wtotal.sum(),'Wmean_left':w.Wmean.sum(),'Wfluct_left':w.Wfluct.sum(),
            'Wfluct_rate_left_W':w.Wfluct.sum()/((b-a)*dt), 'D_total_left':w.D_total.sum(),
            'Wcontact_reference_left':w.Wcontact_reference.sum(),
            'modal_v_AC_RMS_norm':float(np.sqrt(np.sum(w.v_AC_RMS**2))),
            'elapsed_s':record['elapsed_s'],'threads':record['threads'],
            'profile':record['effective_parameters']['glottal_method_profile']}
        # Every-step force averages are available for the full evaluation
        # window. Do not pretend decimated CSV is every-step data in subwindows.
        if label=='evaluation':
            row.update(positive_opening_force_left_N=pos['mean_opening_N'],negative_opening_force_left_N=neg['mean_opening_N'],
                positive_pressure_area_left_mm2=pos['mean_area_m2']*1e6,
                negative_pressure_area_left_mm2=neg['mean_area_m2']*1e6)
        rows.append(row)
    result=pd.DataFrame(rows)
    early=result[result.window=='early'].iloc[0]; late=result[result.window=='late'].iloc[0]
    for key in ['Amin_AC_RMS_mm2','uyL_AC_RMS_mm','Ug_AC_RMS_ml_s','modal_v_AC_RMS_norm']:
        ratio=float(late[key]/early[key]) if early[key]>0 else np.nan
        result[key+'_late_over_early']=ratio
    ratio=float(late.uyL_AC_RMS_mm/early.uyL_AC_RMS_mm) if early.uyL_AC_RMS_mm>0 else np.nan
    state='growing' if ratio>1.1 else 'decaying' if ratio<.9 else 'approximately_stationary'
    result['fixed_probe_amplitude_trend']=state
    result['area_closure_threshold_mm2']=.01
    result['amplitude_log_rate_per_s']=np.log(ratio)/.05 if ratio>0 else np.nan
    gates.update(passed=True,work_summary_recomputed=True,final_interval_excluded=True,
                  exact_equal_windows=[4999,4999],model_mode_count=[20,20])
    return result,work,gates,flow


def figures(summary,works,flows,directory):
    os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='kanaya-ablation-mpl-'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names=list(flows)
    fig,axs=plt.subplots(4,1,figsize=(13,12),sharex=True)
    for name,flow in flows.items():
        f=flow.iloc[10000:19999]
        for ax,key in zip(axs,['Amin_mm2','Ug_m3_s','uyL_mm','xSep_mm']):
            ax.plot(f.time_s,f[key]*(1e6 if key=='Ug_m3_s' else 1),label=name,lw=.8)
    for ax,label in zip(axs,['Amin [mm^2]','Ug [ml/s]','fixed uyL [mm]','wall end / sep x [mm]']):
        ax.set_ylabel(label);ax.grid(alpha=.2)
    axs[0].axhline(.01,color='k',ls='--',lw=.7)
    axs[0].legend(fontsize=7,ncol=2);axs[-1].set_xlabel('time [s]')
    fig.tight_layout();fig.savefig(directory/'comparisons.png',dpi=160);plt.close(fig)
    total=works[(works.side=='L')&(works.window=='evaluation')].groupby('case')[['Wtotal','Wmean','Wfluct','D_total']].sum().loc[names]
    fig,ax=plt.subplots(figsize=(12,6));x=np.arange(len(total))
    for shift,key in [(-.3,'Wtotal'),(-.1,'Wmean'),(.1,'Wfluct'),(.3,'D_total')]:
        ax.bar(x+shift,total[key]*1e6,.2,label=key)
    ax.set_xticks(x,total.index,rotation=25,ha='right');ax.axhline(0,color='k',lw=.5)
    ax.set_ylabel('left normalized work [micro-J equivalent]');ax.legend();fig.tight_layout()
    fig.savefig(directory/'work_comparison.png',dpi=160);plt.close(fig)
    ev=summary[summary.window=='evaluation'].set_index('case').loc[names]
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for ax,key in zip(axs,['positive_opening_force_left_N','negative_opening_force_left_N']):
        ax.bar(names,ev[key]*1e3);ax.set_ylabel(key+' [mN]');ax.tick_params(axis='x',rotation=35)
    fig.tight_layout();fig.savefig(directory/'forces.png',dpi=160);plt.close(fig)
    bench=OUT/'benchmark_final_results.csv'
    if bench.exists():
        b=pd.read_csv(bench);b=b[b.name.str.match(r'benchmark/modern_t\d+$')].sort_values('threads')
        fig,axs=plt.subplots(1,2,figsize=(10,4))
        axs[0].plot(b.threads,b.elapsed_s,'o-');axs[0].set(xlabel='OMP threads',ylabel='whole run wall [s]')
        axs[1].plot(b.threads,b.core_loop_reported_s,'o-');axs[1].set(xlabel='OMP threads',ylabel='reported core loop [s]')
        fig.tight_layout();fig.savefig(directory/'cpu_benchmark.png',dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+');p.add_argument('--tag',default='comparison')
    args=p.parse_args()
    target=OUT/(args.tag+'_summary.csv')
    if target.exists():raise FileExistsError(target)
    summaries=[];works=[];gates={};flows={}
    for case in args.cases:
        s,w,g,f=analyze(case);summaries.append(s);works.append(w);gates[case]=g;flows[case]=f
    summary=pd.concat(summaries,ignore_index=True);work=pd.concat(works,ignore_index=True)
    summary.to_csv(target,index=False);work.to_csv(OUT/(args.tag+'_modal_work.csv'),index=False)
    dump(OUT/(args.tag+'_validation.json'),gates)
    plotdir=OUT/(args.tag+'_figures');plotdir.mkdir()
    figures(summary,work,flows,plotdir)
    print(summary[summary.window=='evaluation'][['case','Amin_mean_mm2','Amin_min_mm2','closure_fraction',
        'periodic_closure','uyL_AC_RMS_mm','fixed_probe_amplitude_trend','Wfluct_rate_left_W']].to_string(index=False))


if __name__=='__main__':main()
