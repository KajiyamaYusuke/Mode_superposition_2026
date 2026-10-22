#!/usr/bin/env python3
"""Same-window S2/S3 comparison and algebraic replay of unlogged S2 airway states."""
import json
import math
import os
import tempfile
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from next_common import ROOT,OUT,PREVIOUS,S2,sha,dump,physical_parameters
from compare_runs import analyze,numeric_gate


def coefficients(path):
    frame=pd.read_csv(path,float_precision='round_trip')
    return dict(zip(frame.coefficient,frame.value))


def replay_saved_s2(flow,c,params):
    # Reconstruct hidden circuit states from STORED Ug; no NR solve, modal
    # update, contact update or coupled simulation is performed here.
    dt=float(params['dt_s']);ng=int(params['subglottal_sections']);nv=int(params['vocal_tract_sections'])
    ps=float(params['ps_Pa']);ramp=float(params['pressure_ramp_time_sec'])
    uu=[0.]*(ng+1);pu=[0.]*(ng+2);ud=[0.]*nv;pd_=[0.]*nv
    rows=[]
    for r in flow.itertuples():
        n=r.step;u=r.Ug_m3_s;previous=0. if n==0 else flow.iloc[n-1].Ug_m3_s
        lung=ps*(.5*(1.-math.cos(math.pi*r.time_s/ramp)) if r.time_s<ramp else 1.)
        for j in range(ng):pu[j]+=uu[j]-uu[j+1]
        pu[ng]+=uu[ng]-previous;pu[ng+1]+=previous-ud[0]
        uu[0]-=dt/c['Lui']*(dt/c['Cui']*pu[0]-lung)
        uu[1]-=dt/(c['Lui']+c['Lu'])*(dt/c['Cu']*pu[1]-dt/c['Cui']*pu[0]+c['R2']*uu[1])
        for j in range(2,ng+1):uu[j]-=dt/(2.*c['Lu'])*(dt/c['Cu']*pu[j]-dt/c['Cu']*pu[j-1]+c['R2']*uu[j])
        pg=dt/c['Cu']*pu[ng]
        if r.coefficients_active:
            terms=[r.Rk*abs(u)*u,r.Rv*u,(r.Lg+c['La']+c['Lu'])*(u-previous)/dt,
                   dt/c['Ca']*pu[ng+1],-pg]
        else:terms=[0.]*5
        pd_[0]+=dt/c['Ca']*(u-ud[0])
        for j in range(1,nv):pd_[j]+=dt/c['Ca']*(ud[j-1]-ud[j])
        for j in range(nv-1):ud[j]+=dt/c['La']*(pd_[j]-pd_[j+1])
        ud[nv-1]=(pd_[nv-1]+(c['La']+c['Lr'])/dt*ud[nv-1])/((c['La']+c['Lr'])/dt+c['Rr'])
        row=dict(step=n,time_s=r.time_s,Ug_m3_s=u,Pg_Pa=pg,Pd0_after_Pa=pd_[0],Pout_Pa=pd_[-1],
                 lung_Pa=lung,NR_active=r.coefficients_active,NR_residual_Pa=sum(terms))
        row.update({f'Uu_{j}_m3_s':v for j,v in enumerate(uu)})
        row.update({f'Pu_accum_{j}':v for j,v in enumerate(pu)})
        row.update({f'Pu_physical_{j}_Pa':dt/(c['Cui'] if j==0 else c['Cu'])*pu[j] for j in range(ng+1)})
        row.update(dict(zip(['NR_quadratic_Pa','NR_viscous_Pa','NR_inertia_Pa','NR_downstream_Pa','NR_upstream_Pa'],terms)))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    target=OUT/'s2_vs_s3';target.mkdir(exist_ok=False)
    s3=OUT/'S3_legacy_upstream_circuit'
    first=json.loads((S2/'manifest.json').read_text());second=json.loads((s3/'manifest.json').read_text())
    def tokens(path):return [s.strip() for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    assert tokens(s3/'params.txt')==tokens(S2/'params.txt')+['upstreamAcousticProfile = legacy_main']
    assert sha(S2/'params.txt')==first['parameters_sha256']
    allowed={'inlet_radius_m','inlet_length_m','subglottal_length_m','subglottal_sections'}
    changes={k:{'S2':v,'S3':second['effective_parameters'][k]}for k,v in first['effective_parameters'].items()
             if k!='parameter_file' and second['effective_parameters'].get(k)!=v}
    assert set(changes)==allowed,changes
    effective=second['effective_parameters']
    assert effective['upstream_acoustic_profile']=='legacy_main' and effective['subglottal_Lu_definition']=='rho*dx/A'
    assert float(effective['vocal_tract_length_m'])==.175 and int(effective['subglottal_sections'])==3
    c2=coefficients(OUT/'validation/s2_omitted_v2_300/output/csv/acoustic_coefficients.csv')
    c3=coefficients(s3/'output/csv/acoustic_coefficients.csv')
    for key in ['R2','La','Ca','Lr','Rr']:assert c2[key]==c3[key]
    pd.DataFrame([dict(coefficient=k,S2=c2[k],S3=c3[k],ratio_S3_over_S2=c3[k]/c2[k])for k in c2]).to_csv(target/'acoustic_coefficient_comparison.csv',index=False)
    dump(target/'case_parameter_diff.json',dict(passed=True,only_change='complete upstream acoustic block',
        positional_and_existing_named_tokens_identical=True,appended_option='upstreamAcousticProfile = legacy_main',
        effective_changes=changes,S2_parameters_sha256=sha(S2/'params.txt'),S3_parameters_sha256=sha(s3/'params.txt'),
        S2_coefficient_basis='300-step same-S2 initialization; physical output bitwise-verified against archived previous binary',
        coefficients_S2=c2,coefficients_S3=c3,legacy_source_commit=second['legacy_source_commit']))
    summaries=[];works=[];validation={};flows={}
    for label,directory in [('S2',S2),('S3',s3)]:
        summary,work,gates,flow=analyze(str(directory))
        summary['case']=label;work['case']=label
        # Add explicit right fixed-material-probe statistics and local cycle minima.
        for window,a,b in [('early',10000,14999),('late',15000,19999),('evaluation',10000,19999)]:
            subset=flow.iloc[a:b];mask=summary.window==window
            summary.loc[mask,'uyR_mean_mm']=subset.uyR_mm.mean()
            summary.loc[mask,'uyR_AC_RMS_mm']=subset.uyR_mm.std(ddof=0)
            peaks,_=find_peaks(-subset.Amin_mm2.to_numpy(),distance=600,prominence=.1)
            pd.DataFrame(subset.iloc[peaks][['step','time_s','Amin_mm2','xMin_mm']]).to_csv(target/f'{label}_{window}_local_minima.csv',index=False)
        summaries.append(summary);works.append(work);validation[label]=gates;flows[label]=flow
    summary=pd.concat(summaries,ignore_index=True);work=pd.concat(works,ignore_index=True)
    summary.to_csv(target/'comparison_summary.csv',index=False);work.to_csv(target/'comparison_modal_work.csv',index=False)
    replay_parameters=dict(first['effective_parameters']);replay_parameters.update(physical_parameters(S2/'params.txt'))
    acoustic2=replay_saved_s2(flows['S2'],c2,replay_parameters)
    validation['S2_airway_reconstruction']={
        'Pg':numeric_gate(acoustic2.Pg_Pa,flows['S2'].Pg_Pa,atol=1e-7,rtol=1e-10),
        'Pd0':numeric_gate(acoustic2.Pd0_after_Pa,flows['S2'].Pd0_Pa,atol=1e-7,rtol=1e-10),
        'NR_max_Pa':float(acoustic2.loc[acoustic2.NR_active==1,'NR_residual_Pa'].abs().max()),
        'method':'algebraic state propagation driven by STORED Ug; no NR/modal/contact solve',
        'is_new_coupled_calculation':False}
    assert validation['S2_airway_reconstruction']['NR_max_Pa']<1e-7
    acoustic2.to_csv(target/'S2_acoustic_states_reconstructed.csv',index=False)
    acoustic3=pd.read_csv(s3/'output/csv/acoustic_trace.csv',float_precision='round_trip')
    assert len(acoustic3)==20000 and np.isfinite(acoustic3.select_dtypes('number')).all().all()
    assert len([c for c in acoustic3 if c.startswith('Uu_')])==4
    assert len([c for c in acoustic3 if c.startswith('Pu_accum_')])==5
    numeric_gate(acoustic3.Ug_m3_s,flows['S3'].Ug_m3_s)
    numeric_gate(acoustic3.Pg_Pa,flows['S3'].Pg_Pa)
    active=acoustic3.NR_active==1
    assert acoustic3.loc[active,'NR_converged'].eq(1).all()
    assert acoustic3.loc[active,'NR_residual_Pa'].abs().max()<1e-8
    validation['S3_actual_acoustic_trace']=dict(passed=True,rows=20000,
        maximum_NR_postcheck_Pa=float(acoustic3.loc[active,'NR_residual_Pa'].abs().max()),
        Uu_size=4,Pu_size=5,NR_iterations='-1 means not instrumented; original solver target1e-9, postcheck1e-8')
    old_summary=pd.read_csv(PREVIOUS/'comparison_summary.csv')
    old=old_summary[old_summary.case=='S2_legacy_vt0175'].set_index('window')
    new=summary[summary.case=='S2'].set_index('window')
    columns=[k for k in old.select_dtypes('number') if k in new]
    for column in columns:
        a=old.loc[new.index,column].to_numpy();b=new[column].to_numpy()
        assert np.array_equal(np.isnan(a),np.isnan(b)),column
        mask=np.isfinite(a)
        numeric_gate(a[mask],b[mask],atol=1e-12,rtol=1e-9)
    validation['passed']=True
    dump(target/'validation.json',validation)
    plots(target,summary,work,flows,acoustic2,acoustic3)
    print(summary[summary.window=='evaluation'][['case','Amin_mean_mm2','Amin_min_mm2','closure_fraction',
        'periodic_closure','uyL_AC_RMS_mm','Wfluct_left','D_total_left']].to_string(index=False))


def plots(target,summary,work,flows,a2,a3):
    os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='kanaya-s3-mpl-'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    path=target/'figures';path.mkdir()
    fig,axs=plt.subplots(3,1,figsize=(12,10),sharex=True)
    for label,flow in flows.items():
        f=flow.iloc[10000:19999]
        for ax,key in zip(axs,['Amin_mm2','uyL_mm','Ug_m3_s']):ax.plot(f.time_s,f[key]*(1e6 if key=='Ug_m3_s'else 1),label=label,lw=.9)
    axs[0].axhline(.01,color='k',ls='--');axs[0].legend()
    for ax,label in zip(axs,['Amin [mm^2]','fixed-material uyL [mm]','Ug [ml/s]']):ax.set_ylabel(label);ax.grid(alpha=.2)
    axs[-1].set_xlabel('time [s]');fig.tight_layout();fig.savefig(path/'area_and_displacement.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for label,a,ng in [('S2',a2,10),('S3',a3,3)]:
        a=a.iloc[10000:19999]
        axs[0].plot(a.time_s,a.Pu_physical_0_Pa,label=label+' inlet')
        axs[1].plot(a.time_s,a.Pg_Pa,label=label+' glottal');axs[1].plot(a.time_s,a.Pd0_after_Pa,ls='--',label=label+' downstream')
        axs[2].plot(a.time_s,a[f'Uu_{ng}_m3_s']*1e6,label=label+' last upstream Uu')
    for ax,label in zip(axs,['inlet pressure [Pa]','pressure [Pa]','upstream flow [ml/s]']):ax.set_ylabel(label);ax.legend(fontsize=8);ax.grid(alpha=.2)
    axs[-1].set_xlabel('time [s]');fig.tight_layout();fig.savefig(path/'upstream_pressure_flow.png',dpi=160);plt.close(fig)
    totals=work[(work.side=='L')&(work.window=='evaluation')].groupby('case')[['Wtotal','Wmean','Wfluct','Wpos','Wneg','D_total']].sum()
    fig,axs=plt.subplots(1,2,figsize=(12,5));totals.mul(1e3).plot.bar(ax=axs[0]);axs[0].set_ylabel('left modal work [mJ; D reference]')
    s=summary[summary.window=='evaluation'].set_index('case')
    s[['adopted_contact_interval_fraction']].plot.bar(ax=axs[1]);axs[1].set_ylabel('adopted local-contact fraction, not CQ')
    fig.tight_layout();fig.savefig(path/'work_and_contact.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
