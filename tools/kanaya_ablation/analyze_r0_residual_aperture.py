#!/usr/bin/env python3
"""All-mode offline geometry replay using unchanged production C++ intersections."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from next_common import ROOT,OUT,PREVIOUS,BUILD,environment,sha,dump
from compare_runs import numeric_gate,spectrum


def interval_metrics(frame,low,high):
    area=coverage=opened=negative=maxgap=0.
    for row in frame.itertuples():
        a=max(low,row.z0_mm);b=min(high,row.z1_mm)
        if b<=a:continue
        length=b-a;coverage+=length
        h0=row.h0_mm+(row.h1_mm-row.h0_mm)*(a-row.z0_mm)/(row.z1_mm-row.z0_mm)
        h1=row.h0_mm+(row.h1_mm-row.h0_mm)*(b-row.z0_mm)/(row.z1_mm-row.z0_mm)
        maxgap=max(maxgap,h0,h1)
        if h0<=0 and h1<=0:
            if min(h0,h1)<0:negative+=length
            continue
        if h0<0:
            positive=length*h1/(h1-h0);negative+=length-positive;h0=0.
        elif h1<0:
            positive=length*h0/(h0-h1);negative+=length-positive;h1=0.
        else:positive=length
        opened+=positive;area+=.5*(h0+h1)*positive
    return dict(area_mm2=area,valid_z_width_mm=coverage,open_z_width_mm=opened,
                negative_gap_z_width_mm=negative,max_positive_gap_mm=maxgap,
                mean_positive_gap_over_open_z_mm=area/opened if opened>0 else 0.)


def sample_profile(frame,z):
    values=np.full((len(z),4),np.nan)
    for row in frame.itertuples():
        mask=(z>=row.z0_mm)&(z<=row.z1_mm)
        t=(z[mask]-row.z0_mm)/(row.z1_mm-row.z0_mm)
        yl=row.yL0_mm+t*(row.yL1_mm-row.yL0_mm)
        yr=row.yR0_mm+t*(row.yR1_mm-row.yR0_mm)
        h=row.h0_mm+t*(row.h1_mm-row.h0_mm)
        values[mask]=np.stack([yl,yr,h,np.maximum(h,0.)],axis=-1)
    return values


def main():
    target=OUT/'r0_aperture';target.mkdir(parents=True,exist_ok=False)
    case=PREVIOUS/'R0_modern';base=case/'output/kanaya_diagnostic'
    meta=json.loads((base/'run_metadata.json').read_text())
    assert meta['nmode']==[20,20] and meta['integrated_intervals']==19999 and meta['excluded_step']==19999
    raw_path=base/'modal_native.bin';assert raw_path.stat().st_size==19999*2*20*7*8
    raw=np.memmap(raw_path,dtype='<f8',mode='r',shape=(19999,2,20,7))
    numeric_gate(raw[1:,:,:,4]-raw[:-1,:,:,4],raw[:-1,:,:,6],atol=1e-22,rtol=1e-12)
    flow=pd.read_csv(base/'flow_geometry.csv',float_precision='round_trip')
    evaluation=flow.iloc[10000:19999]
    peak_hz,_=spectrum(evaluation.Ug_m3_s.to_numpy(),meta['dt_s'])
    minimum_spacing=int(.7/(peak_hz*meta['dt_s']))
    prominence=.25*(np.quantile(evaluation.Amin_mm2,.95)-np.quantile(evaluation.Amin_mm2,.05))
    peaks,_=find_peaks(-evaluation.Amin_mm2.to_numpy(),distance=minimum_spacing,prominence=prominence)
    minima=np.array(peaks+10000,dtype=int)
    assert len(minima)>=5,'fewer than five distinct minima; inspect without inventing events'
    global_min=int(evaluation.Amin_mm2.idxmin())
    selected=[global_min]
    for step in minima[np.linspace(0,len(minima)-1,5,dtype=int)]:
        if int(step) not in selected:selected.append(int(step))
    requests=[(i,n) for i,n in enumerate(selected)]+[(99,15001)]
    (target/'requests.txt').write_text(''.join(f'{event} {step}\n' for event,step in requests))
    replay_cpus=os.environ.get('KANAYA_REPLAY_CPUS','0,2,4,6')
    env=environment(len(replay_cpus.split(',')));env['SIMULATION_RUN_DIR']=str(target/'initialization_output')
    env.pop('SIMULATION_ACOUSTIC_TRACE',None)
    (target/'work').mkdir()
    command=['taskset','-c',replay_cpus,str(BUILD/'kanaya_section_replay'),str(case/'params.txt'),str(raw_path),
             str(target/'requests.txt'),str(target)]
    with (target/'replay.log').open('x') as log:
        subprocess.run(command,cwd=target/'work',env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    sections=pd.read_csv(target/'replay_sections.csv',float_precision='round_trip')
    pieces=pd.read_csv(target/'replay_intervals.csv',float_precision='round_trip')
    checks=pd.read_csv(target/'replay_checks.csv',float_precision='round_trip')
    validation={};event_rows=[];regions=[];profiles=[];contact_rows=[]
    contact=pd.read_csv(base/'contact_regions.csv',float_precision='round_trip')
    for row in checks.itertuples():
        original=flow.iloc[row.step]
        validation[str(row.event)]={
            'area':numeric_gate([row.Amin_mm2],[original.Amin_mm2],atol=1e-10,rtol=1e-8),
            'x':numeric_gate([row.x_min_mm],[original.xMin_mm],atol=1e-10,rtol=1e-10),
            'fixed_probe':numeric_gate([row.probe_uxL_mm,row.probe_uyL_mm,row.probe_uxR_mm,row.probe_uyR_mm],
                                       [original.uxL_mm,original.uyL_mm,original.uxR_mm,original.uyR_mm],atol=1e-10)}
        event_rows.append(dict(event=row.event,step=row.step,time_s=original.time_s,
            purpose='VTU time-alignment control' if row.event==99 else 'distinct-cycle area minimum',
            min_index=row.min_index,x_min_mm=row.x_min_mm,A_solver_mm2=original.Amin_mm2,A_replay_mm2=row.Amin_mm2,
            absolute_area_error_mm2=abs(row.Amin_mm2-original.Amin_mm2),reference_z_min_mm=row.reference_z_min_mm,
            reference_z_max_mm=row.reference_z_max_mm))
        eligible=set(range(max(0,row.min_index-2),min(50,row.min_index+3)))
        adopted=contact[(contact.step==row.step)&(contact.side=='L')]
        event_sections=sections[sections.event==row.event]
        for c in adopted.itertuples():
            k=int(event_sections.iloc[np.argmin(np.abs(event_sections.x_mm-c.x_mean_mm))].section) if np.isfinite(c.x_mean_mm) else -1
            if k>=0:eligible.add(k)
            contact_rows.append(dict(event=row.event,geometry_step=row.step,geometry_time_s=original.time_s,
                contact_predicted_time_s=(row.step+1)*1e-5,ap_region=c.ap_region,pairs=c.pairs,
                max_penetration_m=c.max_penetration_m,contact_x_mean_mm=c.x_mean_mm,nearest_section=k,
                note='aggregate last-adopted t_(n+1) predicted pairs; not same-time signed-gap samples'))
        for k in sorted(eligible):
            frame=pieces[(pieces.event==row.event)&(pieces.section==k)]
            if frame.empty:continue
            section=event_sections[event_sections.section==k].iloc[0]
            total=interval_metrics(frame,-np.inf,np.inf)
            numeric_gate([total['area_mm2']],[section.area_mm2],atol=1e-10,rtol=1e-8)
            for width in [.5,1.,2.]:
                values=[]
                bounds=[('low_z_end',-np.inf,row.reference_z_min_mm+width),
                        ('interior',row.reference_z_min_mm+width,row.reference_z_max_mm-width),
                        ('high_z_end',row.reference_z_max_mm-width,np.inf)]
                for name,low,high in bounds:
                    stats=interval_metrics(frame,low,high);values.append(stats['area_mm2'])
                    regions.append(dict(event=row.event,step=row.step,section=k,x_mm=section.x_mm,
                        is_minimum=k==row.min_index,end_width_mm=width,region=name,
                        area_share_pct=100*stats['area_mm2']/total['area_mm2'] if total['area_mm2']>0 else np.nan,
                        support_z_min_mm=frame.z0_mm.min(),support_z_max_mm=frame.z1_mm.max(),**stats))
                numeric_gate([sum(values)],[total['area_mm2']],atol=1e-12)
            z=np.unique(np.r_[np.linspace(row.reference_z_min_mm,row.reference_z_max_mm,681),frame.z0_mm,frame.z1_mm])
            sampled=sample_profile(frame,z)
            for zi,values in zip(z,sampled):
                profiles.append(dict(event=row.event,section=k,x_mm=section.x_mm,z_mm=zi,valid=int(np.isfinite(values[0])),
                                     yL_mm=values[0],yR_mm=values[1],gap_signed_mm=values[2],gap_positive_mm=values[3]))
    events=pd.DataFrame(event_rows);region_frame=pd.DataFrame(regions);profile_frame=pd.DataFrame(profiles)
    events.to_csv(target/'selected_events.csv',index=False)
    region_frame.to_csv(target/'residual_aperture_regions.csv',index=False)
    profile_frame.to_csv(target/'residual_aperture_by_z.csv',index=False)
    pd.DataFrame(contact_rows).to_csv(target/'contact_x_comparisons.csv',index=False)
    pd.DataFrame([dict(threshold_mm2=t,official=t==.01,count=int((evaluation.Amin_mm2<=t).sum()),
        total=len(evaluation),fraction=float((evaluation.Amin_mm2<=t).mean()))for t in [.01,.05,.1]]).to_csv(target/'threshold_sensitivity.csv',index=False)
    flow.iloc[minima][['step','time_s','Amin_mm2','xMin_mm']].to_csv(target/'all_cycle_minima.csv',index=False)
    # Independent existing deformation output is AFTER n=15000, i.e. q_15001.
    vtu=PREVIOUS/'validation/pristine_modern_full/output/result/deform_combined0000.vtu'
    xml=ET.parse(vtu);points=np.fromstring(xml.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)
    reconstructed=pd.read_csv(target/'control_coordinates.csv',float_precision='round_trip')[['x_mm','y_mm','z_mm']].to_numpy()
    vtu_gate=numeric_gate(points,reconstructed,atol=1e-11,rtol=1e-13)
    dump(target/'validation.json',dict(passed=True,events=validation,native_shape=[19999,2,20,7],
        q_time='q_n from native field 4; production State::mode2uf followed by geometry-only displacement assignment',
        q_continuity_checked=True,vtu_time_control={'existing_output_step':15000,'native_q_step':15001,**vtu_gate},
        input_native_sha256=sha(raw_path),section_source_sha256=sha(ROOT/'src/fluid/channelSections.cpp'),
        section_reference_revision='bc62cf8f3a2c01fb1ed0a33f6d4a6864f6c1d2d5 (unchanged current file)',
        production_coupled_runs=0,anatomical_AP_orientation='unknown; only low-z/high-z labels used',
        cycle_selection={'Ug_spectrum_peak_Hz_not_f0':peak_hz,'minimum_peak_spacing_steps':minimum_spacing,
                         'minimum_prominence_mm2':prominence},
        no_missing_segment_bridging=True,official_closure_threshold_mm2=.01))
    plots(target,events,region_frame,profile_frame,sections,pieces,global_min)
    print('PASS R0 geometry replay; event minima:',events[['event','step','A_replay_mm2','absolute_area_error_mm2']].to_string(index=False))
    print(region_frame[(region_frame.is_minimum)&(region_frame.end_width_mm==1)&(region_frame.event!=99)][
        ['event','region','area_mm2','area_share_pct','open_z_width_mm','negative_gap_z_width_mm','max_positive_gap_mm']].to_string(index=False))


def plots(target,events,regions,profiles,sections,pieces,global_step):
    os.environ.setdefault('MPLCONFIGDIR',tempfile.mkdtemp(prefix='kanaya-aperture-mpl-'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figures=target/'figures';figures.mkdir()
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True)
    for row in events[events.event!=99].itertuples():
        p=profiles[(profiles.event==row.event)&(profiles.section==row.min_index)]
        axs[0].plot(p.z_mm,p.gap_signed_mm,label=f'n={row.step}, x={row.x_min_mm:.4f} mm')
        axs[1].plot(p.z_mm,p.gap_positive_mm)
    axs[0].axhline(0,color='k',lw=.7);axs[0].set_ylabel('signed gap [mm]');axs[0].legend(fontsize=8)
    axs[1].set(xlabel='physical z [mm] (anatomical direction unknown)',ylabel='positive gap [mm]')
    for ax in axs:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(figures/'r0_ap_gap_profiles.png',dpi=160);plt.close(fig)
    subset=regions[regions.is_minimum&(regions.end_width_mm==1)&(regions.event!=99)]
    pivot=subset.pivot(index='event',columns='region',values='area_mm2')[['low_z_end','interior','high_z_end']]
    ax=pivot.plot.bar(stacked=True,figsize=(10,5));ax.set(ylabel='open area contribution [mm^2]',xlabel='event; end bands = 1 mm')
    ax.figure.tight_layout();ax.figure.savefig(figures/'r0_ap_area_contributions.png',dpi=160);plt.close(ax.figure)
    event=int(events[events.step==global_step].iloc[0].event)
    s=sections[sections.event==event].sort_values('section');z=np.linspace(0,17,681)
    grid=np.stack([sample_profile(pieces[(pieces.event==event)&(pieces.section==k)],z)[:,2] for k in range(50)])
    pd.DataFrame([dict(section=k,x_mm=s.iloc[k].x_mm,z_mm=zi,gap_signed_mm=g)
                  for k in range(50) for zi,g in zip(z,grid[k])]).to_csv(target/'xz_gap_map.csv',index=False)
    fig,ax=plt.subplots(figsize=(12,6));v=min(float(np.nanmax(np.abs(grid))),1.)
    pcm=ax.pcolormesh(z,s.x_mm,np.ma.masked_invalid(grid),cmap='coolwarm',vmin=-v,vmax=v,shading='nearest')
    ax.axhline(float(events[events.event==event].x_min_mm.iloc[0]),color='k',ls='--',lw=.8)
    ax.set(xlabel='physical z [mm]',ylabel='physical x [mm]',title=f'R0 signed gap at n={global_step}; invalid cuts masked')
    fig.colorbar(pcm,ax=ax,label='signed gap [mm] (color clipped if >1 mm)')
    fig.tight_layout();fig.savefig(figures/'r0_xz_gap_map.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
