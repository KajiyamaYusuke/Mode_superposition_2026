#!/usr/bin/env python3
"""Single-file Ug event jitter and flow CPP, specification section 11 (2026-11-05).

No solver execution, automatic input discovery, or input modification. Jitter is
measured from actual candidate peaks; ACF is only a period proposal. CPP is a
research-specific 10log10(real-cepstrum**2) residual, not CPPS or clinical CPP.
"""
import argparse
import atexit
import contextlib
import copy
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import unittest
from datetime import datetime, timezone

sys.dont_write_bytecode = True
if 'MPLCONFIGDIR' not in os.environ:
    _mpl_cache = tempfile.TemporaryDirectory(prefix='airflow_metrics_mpl_')
    os.environ['MPLCONFIGDIR'] = _mpl_cache.name
    atexit.register(_mpl_cache.cleanup)
import numpy as np
from scipy import signal, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TOOL_ID = 'analyze_airflow_metrics_section11_v1'
TRACKER_VERSION = 'reanchor_context_recovery_v1'
# Separate, fixed prior for the new method; the legacy TRACKING_CONFIG is unchanged.
GLOBAL_CYCLE_CONFIG = {
    'version':'global_cycle_dag_v1', 'profile_name':'balanced-v1',
    'pitch_min_hz':60., 'pitch_max_hz':400., 'acf_window_sec':.100, 'acf_hop_sec':.020,
    'prominence_percent':10., 'excursion_percent':10., 'noise_sigma_multiplier':2.,
    'search_min_period_ratio':.55, 'search_max_period_ratio':1.65,
    'effective_shape_correlation':.65, 'incompatible_shape_correlation':-.25,
    'node_prominence_weight':1., 'node_excursion_weight':.75, 'node_acf_weight':.5,
    'node_basin_weight':.25, 'secondary_lobe_penalty':.75,
    'edge_shape_weight':.75, 'edge_phase_weight':.25, 'edge_acf_weight':.25,
    'period_regularization_weight':0., 'root_penalty':2., 'restart_penalty':2.,
    'intermediate_pulse_strength_ratio':.5,
    'phase_shift_warning_delta':.35,
    'ambiguity_score_margin':.20, 'high_acf_confidence':.50,
    'max_dp_candidates':2000, 'top_paths':3, 'max_alternative_paths':12,
}
CPP_CONFIG = {
    'target_fs_hz': 10000, 'max_frequency_hz': 5000,
    'frame_sec': .060, 'hop_sec': .005, 'window': 'hann', 'nfft': 2048,
    'pitch_min_hz': 60, 'pitch_max_hz': 400,
    'spectrum_epsilon_ratio': 1e-12, 'cepstrum_epsilon_ratio': 1e-12,
    'trend_start_sec': .001, 'trend_end_sec': .020, 'trend_fit': 'theil_sen',
    'peak_interpolation': 'quadratic', 'time_smoothing': False,
    'quefrency_smoothing': False, 'preemphasis': False,
}
TRACKING_CONFIG = {
    'pitch_min_hz': 60, 'pitch_max_hz': 400,
    'acf_window_sec': .100, 'acf_hop_sec': .020,
    'primary_prominence_percent': 10,
    'prominence_sweep_percent': [3, 5, 8, 10, 12, 15, 20],
    'candidate_min_distance_ratio': .20,
    'search_min_period_ratio': .65, 'search_max_period_ratio': 1.45,
    'template_halfwidth_ratio': .40, 'correlation_weight': .85,
    'prominence_weight': .15, 'min_correlation': .65,
    'min_score_margin': .10, 'high_confidence_correlation': .85,
    'reference_period_median_count': 5,
}
QC_CONFIG = {
    'min_cycles': 10, 'coverage_min': .90, 'median_correlation_min': .85,
    'ambiguous_fraction_max': .05, 'acf_relative_difference_max': .15,
    'threshold_event_match_min': .90, 'event_tolerance_raw_samples': 1,
    'threshold_jitter_absolute_pp': .2, 'threshold_jitter_relative': .20,
    'acf_minimum': .20, 'acf_near_best_tolerance': .05,
    'stronger_harmonic_margin': .02, 'anchor_center_radius_periods': 2,
    'padding_periods': 3, 'diagnostic_search_ratios': [.4, 2.5],
    'cpp_f0_mismatch_relative': .20,
}
CPP_DEFINITION = {
    'label': 'flow_CPP_dB', 'spectrum': 'P=abs(rfft((frame-mean(frame))*symmetric_Hann,nfft))**2',
    'real_cepstrum': 'c=irfft(ln(max(P,1e-12*max(P))),nfft)',
    'power_cepstrum_db': 'C_dB=10*log10(max(c**2,1e-12*max(c[q>0]**2)))',
    'background': 'Theil-Sen C_dB versus seconds, all positive q in [0.001,0.020]',
    'frame_cpp': 'maximum local peak of D=C_dB-B in pitch range; quadratic D refinement',
    'aggregate': 'arithmetic mean of valid frame CPP; no clipping or smoothing',
    'resampling': 'scipy.signal.resample_poly, Kaiser beta=5, constant padding, no upsampling below 10kHz',
    'not_equivalent_to': 'Praat CPPS, ADSV, or clinical voice CPP',
    'references': ['https://praat.org/manual/Spectrum__To_PowerCepstrum.html',
                   'https://praat.org/manual/PowerCepstrum__Get_peak_prominence___.html'],
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()


def clean_json(value):
    if isinstance(value, dict): return {str(k): clean_json(v) for k, v in value.items() if not str(k).startswith('_')}
    if isinstance(value, (list, tuple)): return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray): return clean_json(value.tolist())
    if isinstance(value, (bool, np.bool_)): return bool(value)
    if isinstance(value, (int, np.integer)): return int(value)
    if isinstance(value, (float, np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value, Path): return str(value)
    return value


def dump_json(path, value):
    Path(path).write_text(json.dumps(clean_json(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def write_csv(path, rows, columns=None):
    rows = list(rows)
    columns = columns or list(dict.fromkeys(k for r in rows for k in r if not k.startswith('_')))
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            cooked = {}
            for key in columns:
                v = clean_json(row.get(key))
                cooked[key] = json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v
            writer.writerow(cooked)


def resolve_dt(dt, manifest):
    other = None
    if manifest:
        values = []
        for line in Path(manifest).read_text().splitlines():
            if '=' in line:
                key, value = line.split('=', 1)
                if key.strip() == 'dt_s': values.append(float(value.strip()))
        if len(values) != 1: raise ValueError('manifest must contain exactly one dt_s')
        other = values[0]
    if dt is None: dt = other
    if dt is None or not np.isfinite(dt) or dt <= 0: raise ValueError('positive --dt or a manifest dt_s is required')
    if other is not None and not math.isclose(dt, other, rel_tol=1e-12, abs_tol=0):
        raise ValueError('--dt and manifest dt_s disagree')
    return float(dt)


def read_flow(path, dt):
    path = Path(path).expanduser().resolve()
    before = sha256(path)
    data = np.loadtxt(path, comments='#', ndmin=2)
    if data.ndim != 2 or data.shape[1] != 2 or len(data) < 3: raise ValueError('expected at least three ASCII two-column samples')
    if not np.isfinite(data).all(): raise ValueError('nonfinite step or airflow')
    steps = data[:, 0]
    if np.any(steps < 0) or np.any(steps != np.rint(steps)): raise ValueError('steps must be nonnegative integers')
    steps = steps.astype(np.int64); dif = np.diff(steps)
    if np.any(dif <= 0): raise ValueError('duplicate or decreasing step')
    if np.any(dif != dif[0]): raise ValueError('unequal step stride (including last interval); no interpolation performed: '+str(np.unique(dif).tolist()))
    if sha256(path) != before: raise ValueError('input changed while reading')
    spacing = int(dif[0])*dt; t = steps*dt; x = data[:, 1]
    return {'path': str(path), 'sha256': before, 'steps': steps, 't': t, 'x': x,
            'dt_s': dt, 'stride': int(dif[0]), 'sampling_interval_s': spacing, 'fs': 1/spacing}


def period_metrics(periods):
    p = np.asarray(periods, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any(p <= 0): raise ValueError('periods must be positive and finite')
    n = len(p); mean = float(np.mean(p)) if n else None
    result = {'n_cycles': n, 'mean_period_s': mean,
              'median_period_s': float(np.median(p)) if n else None,
              'median_f0_hz': float(1/np.median(p)) if n else None,
              'period_cv_percent': float(100*np.std(p)/mean) if n else None,
              'jitter_local_percent': None, 'rap_percent': None, 'ppq5_percent': None}
    if n >= 2: result['jitter_local_percent'] = float(100*np.mean(np.abs(np.diff(p)))/mean)
    if n >= 3: result['rap_percent'] = float(100*np.mean(np.abs(p[1:-1]-(p[:-2]+p[1:-1]+p[2:])/3))/mean)
    if n >= 5: result['ppq5_percent'] = float(100*np.mean(np.abs(p[2:-2]-np.convolve(p,np.ones(5)/5,'valid')))/mean)
    return result


def quadratic(values, index):
    if index <= 0 or index >= len(values)-1: return 0., float(values[index]), 'BOUNDARY_FALLBACK'
    a,b,c = values[index-1:index+2]; den = a-2*b+c
    scale = max(abs(a),abs(b),abs(c),np.finfo(float).tiny)
    if not np.isfinite(den) or den >= -32*np.finfo(float).eps*scale: return 0., float(b), 'FLAT_OR_NONCONCAVE_FALLBACK'
    offset = .5*(a-c)/den
    if not np.isfinite(offset) or abs(offset) > .5: return 0., float(b), 'UNSTABLE_FALLBACK'
    return float(offset), float(b-.25*(a-c)*offset), 'QUADRATIC'


def overlap_acf(x):
    x = np.asarray(x, dtype=float); x = x-np.mean(x); n = len(x)
    raw = signal.correlate(x,x,mode='full',method='fft')[n-1:]
    cs = np.r_[0.,np.cumsum(x)]; sq = np.r_[0.,np.cumsum(x*x)]
    lag = np.arange(n); count = n-lag
    s1 = cs[n-lag]; s2 = cs[n]-cs[lag]
    e1 = sq[n-lag]-s1*s1/count; e2 = sq[n]-sq[lag]-s2*s2/count
    cov = raw-s1*s2/count; denom = np.sqrt(np.maximum(e1,0)*np.maximum(e2,0))
    return np.clip(np.divide(cov,denom,out=np.zeros_like(raw),where=denom>np.finfo(float).tiny),-1,1)


def pitch_candidates(x, fs, config):
    x = np.asarray(x,dtype=float); amp = float(np.percentile(x,95)-np.percentile(x,5))
    if np.ptp(x) <= 32*np.finfo(float).eps*max(np.max(np.abs(x)),np.finfo(float).tiny):
        return {'status':'NO_OSCILLATION','period_s':None,'candidates':[],'alternatives':[], 'flags':['NO_OSCILLATION'], '_acf':np.zeros(len(x))}
    acf = overlap_acf(x)
    low = max(1,int(math.ceil(fs/config['pitch_max_hz']))); high = min(len(x)-3,int(math.floor(fs/config['pitch_min_hz'])))
    if high <= low: return {'status':'INSUFFICIENT_SAMPLES','period_s':None,'candidates':[],'alternatives':[], 'flags':['INSUFFICIENT_SAMPLES'],'_acf':acf}
    peaks,_ = signal.find_peaks(acf)
    peaks = peaks[(peaks>=low)&(peaks<=high)]
    power = np.abs(np.fft.rfft((x-x.mean())*np.hanning(len(x))))**2
    frequency = np.fft.rfftfreq(len(x),1/fs)
    band = (frequency>=config['pitch_min_hz'])&(frequency<=config['pitch_max_hz'])
    rows=[]
    def support(hz):
        mask = abs(frequency-hz)<=max(1/(len(x)/fs),.03*hz)
        return float(np.max(power[mask])/np.max(power[band])) if mask.any() and band.any() and np.max(power[band])>0 else 0.
    for p in peaks:
        rows.append({'lag_samples':int(p),'period_s':p/fs,'frequency_hz':fs/p,'acf':float(acf[p]),'spectrum_support':support(fs/p)})
    good = [r for r in rows if r['acf']>=QC_CONFIG['acf_minimum']]
    if not good: return {'status':'PERIOD_AMBIGUOUS','period_s':None,'candidates':rows,'alternatives':[],'flags':['WEAK_ACF'],'_acf':acf}
    strongest=max(r['acf'] for r in good)
    near=[r for r in good if r['acf']>=strongest-QC_CONFIG['acf_near_best_tolerance']]
    chosen=min(near,key=lambda r:r['period_s']); period=chosen['period_s']; flags=[];alternatives=[]
    for factor in (.5,1.,2.):
        q=period*factor; inside=1/config['pitch_max_hz']<=q<=1/config['pitch_min_hz']
        nearest=min(rows,key=lambda r:abs(r['period_s']-q)) if rows and inside else None
        if nearest and abs(nearest['period_s']/q-1)>.1:nearest=None
        row={'factor':factor,'target_period_s':q,'inside_search_range':inside,'nearest_candidate':nearest};alternatives.append(row)
        if factor!=1 and nearest and nearest['acf']>chosen['acf']+QC_CONFIG['stronger_harmonic_margin']:
            flags.append('PITCH_AMBIGUOUS')
        if factor==.5 and nearest and nearest['acf']>=QC_CONFIG['acf_minimum']:
            flags.append('PITCH_AMBIGUOUS')
    # Large T/2 and T alternatives are retained, even for an exactly periodic pulse train.
    if len(near)>1:flags.append('MULTIPLE_STRONG_ACF_CANDIDATES')
    dominant=float(frequency[np.where(band)[0][np.argmax(power[band])]]) if band.any() else None
    if dominant and (abs(dominant*period-2)<.2 or abs(dominant*period-.5)<.05):flags.append('PITCH_AMBIGUOUS')
    return {'status':'CANDIDATE','period_s':period,'chosen':chosen,'selection_rule':'earliest local maximum within 0.05 of strongest overlap-normalized ACF; no truth claim',
            'candidates':rows,'alternatives':alternatives,'flags':sorted(set(flags)),'amplitude_p95_p5':amp,
            'dominant_spectral_frequency_hz':dominant,
            '_acf':acf,'_frequencies':frequency,'_power':power}


def extract_track_legacy(t, x, fs, start, end, period, prominence_percent, config, method='A', diagnostics=False):
    """Bidirectional fixed-candidate tracking; gaps create new segments, never bridges."""
    polarity = 1 if method=='A' else -1
    y = polarity*x; inside=(t>=start)&(t<end)
    amplitude=float(np.percentile(x[inside],95)-np.percentile(x[inside],5))
    if period is None or amplitude<=0:
        empty={'method':method,'prominence_percent':prominence_percent,'events':[],'segments':[], 'attempts':[], 'breaks':[]}
        if diagnostics:empty['_diagnostic']={'breaks':[],'candidates':[],'reanchor_removals':[]}
        return empty
    pad=QC_CONFIG['padding_periods']*period;lo=max(0,int(np.searchsorted(t,start-pad)));hi=min(len(t),int(np.searchsorted(t,end+pad)))
    view=y[lo:hi]
    raw,props=signal.find_peaks(view,prominence=amplitude*prominence_percent/100)
    separated,_=signal.find_peaks(view,prominence=amplitude*prominence_percent/100,
                                 distance=max(1,int(math.ceil(config['candidate_min_distance_ratio']*period*fs))))
    selectable=set((separated+lo).tolist());half=max(1,int(round(config['template_halfwidth_ratio']*period*fs)))
    events=[];windows={}
    for j,(local,prom) in enumerate(zip(raw,props['prominences'])):
        k=int(local+lo);offset,_,interp=quadratic(y,k)
        row={'event_id':f'{method}_{k}','method':method,'sample_id':k,'raw_time_s':float(t[k]),
             'subsample_time_s':float(t[k]+offset/fs),'offset_us':offset/fs*1e6,'interpolation_status':interp,
             'prominence':float(prom),'normalized_prominence':float(prom/amplitude),
             'selectable':k in selectable,'accepted':False,'segment_id':None,
             'reason':'UNUSED_CANDIDATE' if k in selectable else 'MIN_DISTANCE_CANDIDATE',
             'match_score':None,'correlation':None,'anchor':False,'in_window':bool(start<=t[k]+offset/fs<end)}
        if k-half>=0 and k+half<len(x):
            w=x[k-half:k+half+1].copy();w-=w.mean();norm=np.linalg.norm(w)
            if norm>32*np.finfo(float).eps*max(np.linalg.norm(x[k-half:k+half+1]),np.finfo(float).tiny):windows[j]=w/norm
            else:row['reason']='FLAT_MATCH_WINDOW'
        else:row['reason']='INCOMPLETE_MATCH_WINDOW'
        events.append(row)
    available={j for j,e in enumerate(events) if e['selectable'] and j in windows}
    segments=[];attempts=[];breaks=[];segment_id=0
    diagnostic_breaks=[];diagnostic_candidates=[];reanchor_removals=[];blocked_sources={}
    while available:
        near=[j for j in available if abs(events[j]['raw_time_s']-(start+end)/2)<=QC_CONFIG['anchor_center_radius_periods']*period]
        if not near:
            nearest=min(available,key=lambda j:abs(events[j]['raw_time_s']-(start+end)/2))
            near=[j for j in available if abs(events[j]['raw_time_s']-events[nearest]['raw_time_s'])<=period]
        anchor=max(near,key=lambda j:events[j]['prominence']);segment_id+=1
        events[anchor].update(accepted=True,anchor=True,segment_id=segment_id,reason='ANCHOR',match_score=None,correlation=None)
        accepted={anchor};stop_limits={}
        for direction in (-1,1):
            current=anchor;template_id=anchor;recent=[];ref=period;previous_high=True
            while True:
                prev=events[current];choices=[]
                for j in available-accepted:
                    delta=direction*(events[j]['raw_time_s']-prev['raw_time_s'])
                    if config['search_min_period_ratio']*ref<=delta<=config['search_max_period_ratio']*ref:
                        rho=float(np.dot(windows[template_id],windows[j]));p=events[j]['normalized_prominence'];pp=prev['normalized_prominence']
                        sim=float(np.exp(-abs(np.log((p+1e-12)/(pp+1e-12)))))
                        score=config['correlation_weight']*rho+config['prominence_weight']*sim
                        choices.append((score,rho,j,delta))
                choices.sort(reverse=True)
                margin=choices[0][0]-choices[1][0] if len(choices)>1 else None
                accepted_next=bool(choices and choices[0][1]>=config['min_correlation'] and (margin is None or margin>=config['min_score_margin']))
                for rank,(score,rho,j,delta) in enumerate(choices):
                    attempts.append({'segment_id':segment_id,'direction':direction,'previous_event':prev['event_id'],
                                     'candidate_event':events[j]['event_id'],'template_event':events[template_id]['event_id'],
                                     'period_reference_s':ref,'delta_s':delta,'score':score,'rho':rho,'rank':rank+1,'margin':margin,
                                     'accepted':accepted_next and rank==0})
                    if events[j]['match_score'] is None or score>events[j]['match_score']:
                        events[j]['match_score']=score;events[j]['correlation']=rho
                if not accepted_next:
                    reason='AMBIGUOUS_MATCH' if choices else 'MISSING_EVENT_CANDIDATE'
                    expanded=[{'event_id':events[j]['event_id'],
                               'delta_s':direction*(events[j]['raw_time_s']-prev['raw_time_s']),
                               'delta_over_reference':direction*(events[j]['raw_time_s']-prev['raw_time_s'])/ref} for j in available-accepted if
                              QC_CONFIG['diagnostic_search_ratios'][0]*ref<=direction*(events[j]['raw_time_s']-prev['raw_time_s'])<=QC_CONFIG['diagnostic_search_ratios'][1]*ref]
                    # A failure beyond the padded data is an edge, not evidence of a missing pulse.
                    boundary=t[0] if direction<0 else t[-1]
                    edge=direction*(boundary-prev['raw_time_s']) < config['search_max_period_ratio']*ref
                    if diagnostics:
                        break_id=f'{method}_break_{len(diagnostic_breaks)+1:03d}'
                        best=choices[0] if choices else None;second=choices[1] if len(choices)>1 else None
                        failures=[]
                        if not choices:failures.append('NO_AVAILABLE_CANDIDATE_IN_SEARCH_RANGE')
                        if best and best[1]<config['min_correlation']:failures.append('CORRELATION_BELOW_MIN')
                        if margin is not None and margin<config['min_score_margin']:failures.append('SCORE_MARGIN_BELOW_MIN')
                        row={'break_id':break_id,'method':method,'segment_id':segment_id,'direction':direction,
                             'in_window':bool(start<=prev['subsample_time_s']<end),'break_time_s':prev['subsample_time_s'],
                             'previous_event':prev['event_id'],'previous_raw_time_s':prev['raw_time_s'],
                             'previous_time_s':prev['subsample_time_s'],'template_event':events[template_id]['event_id'],
                             'period_reference_s':ref,'acf_period_s':period,'predicted_next_time_s':prev['raw_time_s']+direction*ref,
                             'search_start_s':min(prev['raw_time_s']+direction*config['search_min_period_ratio']*ref,prev['raw_time_s']+direction*config['search_max_period_ratio']*ref),
                             'search_end_s':max(prev['raw_time_s']+direction*config['search_min_period_ratio']*ref,prev['raw_time_s']+direction*config['search_max_period_ratio']*ref),
                             'candidate_count':len(choices),'best_event':events[best[2]]['event_id'] if best else None,
                             'best_time_s':events[best[2]]['subsample_time_s'] if best else None,
                             'best_delta_s':best[3] if best else None,'best_rho':best[1] if best else None,
                             'best_score':best[0] if best else None,'score_margin':margin,
                             'runner_up_event':events[second[2]]['event_id'] if second else None,
                             'runner_up_score':second[0] if second else None,
                             'best_prominence':events[best[2]]['prominence'] if best else None,
                             'best_normalized_prominence':events[best[2]]['normalized_prominence'] if best else None,
                             'minimum_rho':config['min_correlation'],'minimum_score_margin':config['min_score_margin'],
                             'direct_reasons':failures,'tracking_reason':'DATA_EDGE' if edge else reason}
                        diagnostic_breaks.append(row)
                        score_by_index={v[2]:(rank+1,v) for rank,v in enumerate(choices)}
                        for j,event in enumerate(events):
                            delta=direction*(event['raw_time_s']-prev['raw_time_s'])
                            if not QC_CONFIG['diagnostic_search_ratios'][0]*ref<=delta<=QC_CONFIG['diagnostic_search_ratios'][1]*ref:continue
                            in_search=config['search_min_period_ratio']*ref<=delta<=config['search_max_period_ratio']*ref
                            why=('CURRENT_SEGMENT_ALREADY_ACCEPTED' if j in accepted else 'OTHER_SEGMENT_ALREADY_ACCEPTED' if event['accepted'] else
                                 'REANCHOR_EXCLUSION_REGION' if j in blocked_sources else 'MIN_DISTANCE_SUPPRESSED' if not event['selectable'] else
                                 'INVALID_SHAPE_WINDOW' if j not in windows else 'AVAILABLE')
                            rho=float(np.dot(windows[template_id],windows[j])) if j in windows else None
                            similarity=float(np.exp(-abs(np.log((event['normalized_prominence']+1e-12)/(prev['normalized_prominence']+1e-12)))))
                            score=config['correlation_weight']*rho+config['prominence_weight']*similarity if rho is not None else None
                            diagnostic_candidates.append({'break_id':break_id,'method':method,'event_id':event['event_id'],
                                'candidate_time_s':event['subsample_time_s'],'candidate_raw_time_s':event['raw_time_s'],
                                'delta_s':delta,'delta_over_reference':delta/ref,'in_search_range':in_search,
                                'available_at_break':j in available-accepted,'availability_reason':why,
                                'excluded_by_segment':blocked_sources.get(j),'rho':rho,'score':score,
                                'rank_in_actual_choices':score_by_index[j][0] if j in score_by_index else None,
                                'prominence':event['prominence'],'normalized_prominence':event['normalized_prominence']})
                    breaks.append({'segment_id':segment_id,'direction':direction,'previous_event':prev['event_id'],
                                   'time_s':prev['raw_time_s'],'reason':'DATA_EDGE' if edge else reason,
                                   'expanded_diagnostic_candidates':expanded,'period_reference_s':ref})
                    if choices:
                        for _,_,j,_ in choices:
                            if not events[j]['accepted']:events[j]['reason']='AMBIGUOUS_MATCH'
                    stop_limits[direction]=prev['raw_time_s']+direction*config['search_max_period_ratio']*ref
                    break
                score,rho,j,delta=choices[0];old=events[current]
                events[j].update(accepted=True,anchor=False,segment_id=segment_id,
                                 reason='HIGH_CONFIDENCE' if rho>=config['high_confidence_correlation'] else 'ACCEPTED_LOW_CONFIDENCE',
                                 match_score=score,correlation=rho)
                accepted.add(j);current=j
                if rho>=config['high_confidence_correlation']:
                    if previous_high:
                        recent.append(delta);recent=recent[-config['reference_period_median_count']:]
                        ref=float(np.median(recent))
                    else:recent=[]
                    template_id=j;previous_high=True
                else:recent=[];previous_high=False
                # Temporal changes are retained; ref only affects searching, not timestamps.
        ordered=sorted(accepted,key=lambda j:events[j]['sample_id']);segments.append({'segment_id':segment_id,'event_indices':ordered})
        lower=stop_limits.get(-1,events[ordered[0]]['raw_time_s']);upper=stop_limits.get(1,events[ordered[-1]]['raw_time_s'])
        removed={j for j in available if lower<=events[j]['raw_time_s']<=upper}
        for j in removed-accepted:
            if events[j]['reason']=='UNUSED_CANDIDATE':events[j]['reason']='MULTIPLE_PEAK_UNUSED'
            if diagnostics:
                blocked_sources[j]=segment_id
                reanchor_removals.append({'method':method,'removed_by_segment':segment_id,'event_id':events[j]['event_id'],
                    'time_s':events[j]['subsample_time_s'],'reason':events[j]['reason'],
                    'excluded_range_start_s':lower,'excluded_range_end_s':upper,
                    'meaning':'candidate removed from future anchor search; not proof of a secondary physical peak'})
        available-=removed|accepted
    result={'method':method,'prominence_percent':prominence_percent,'events':events,'segments':segments,
            'attempts':attempts,'breaks':breaks,'candidate_min_distance_samples':max(1,int(math.ceil(config['candidate_min_distance_ratio']*period*fs))),
            'template_halfwidth_samples':half,'amplitude_p95_p5':amplitude}
    if diagnostics:result['_diagnostic']={'breaks':diagnostic_breaks,'candidates':diagnostic_candidates,'reanchor_removals':reanchor_removals}
    return result


def phase_features(t,x,fs,track,period,config):
    """Measured pulse basins; no new fitted score weights or event timestamps."""
    y=x if track['method']=='A' else -x;amplitude=track.get('amplitude_p95_p5',0.)
    opposite,_=signal.find_peaks(-y,prominence=amplitude*track['prominence_percent']/100,
                                distance=max(1,int(math.ceil(config['candidate_min_distance_ratio']*period*fs))))
    rows={};cells={}
    for j,e in enumerate(track['events']):
        k=e['sample_id'];pos=int(np.searchsorted(opposite,k));left=int(opposite[pos-1]) if pos else None
        right=int(opposite[pos]) if pos<len(opposite) else None;cell=(left,right)
        rows[j]={'event_id':e['event_id'],'pulse_cell':cell,'opposite_before_s':float(t[left]) if left is not None else None,
                 'opposite_after_s':float(t[right]) if right is not None else None,
                 'relative_phase_in_basin':float((k-left)/(right-left)) if left is not None and right is not None else None,
                 'rise_duration_s':float((k-left)/fs) if left is not None else None,
                 'fall_duration_s':float((right-k)/fs) if right is not None else None,
                 'rise_height':float(y[k]-y[left]) if left is not None else None,
                 'fall_height':float(y[k]-y[right]) if right is not None else None}
        if e['selectable']:cells.setdefault(cell,[]).append(j)
    for j,f in rows.items():
        peers=cells.get(f['pulse_cell'],[]);f['n_selectable_in_basin']=len(peers)
        f['dominant_in_basin']=bool(peers and track['events'][j]['prominence']>=max(track['events'][k]['prominence'] for k in peers))
    return rows


def extract_track(t,x,fs,start,end,period,prominence_percent,config,method='A',diagnostics=False):
    policy=config.get('reanchor_policy','legacy')
    if policy not in ('legacy','recover-rejected'):raise ValueError('unknown reanchor policy')
    if policy=='legacy':return extract_track_legacy(t,x,fs,start,end,period,prominence_percent,config,method,diagnostics)
    track=extract_track_legacy(t,x,fs,start,end,period,prominence_percent,config,method,True)
    track.update(reanchor_policy=policy,tracker_version=TRACKER_VERSION,recovery_attempts=[],recovered_connections=[],search_limit=False)
    if period is None or not track['events']:return track
    events=track['events'];features=phase_features(t,x,fs,track,period,config);half=track['template_halfwidth_samples'];windows={}
    for j,e in enumerate(events):
        k=e['sample_id']
        if k-half>=0 and k+half<len(x):
            w=x[k-half:k+half+1].copy();w-=w.mean();norm=np.linalg.norm(w)
            if norm>32*np.finfo(float).eps*max(np.linalg.norm(x[k-half:k+half+1]),np.finfo(float).tiny):windows[j]=w/norm
    edges={tuple(sorted((a,b))) for s in track['segments'] for a,b in zip(s['event_indices'],s['event_indices'][1:])}
    used={j for j,e in enumerate(events) if e['accepted']};original_used=set(used);anchor_tried=set();visited=set();attempt_count=0
    cap=int(config.get('reanchor_attempt_cap',max(1000,20*len(events)**2)))
    recovery_breaks=[];recovery_candidates=[];state_reasons={};original_segments=copy.deepcopy(track['segments'])
    seed_breaks=track['_diagnostic']['breaks']
    for row in seed_breaks:row['stage']='LEGACY_SEED'
    def neighbors(j):return sorted(b if a==j else a for a,b in edges if a==j or b==j)
    def conflict(j):
        for a,b in edges:
            if a<j<b:return 'TRACK_CONFLICT'
        cell=features[j]['pulse_cell']
        if cell[0] is not None and cell[1] is not None and any(k!=j and features[k]['pulse_cell']==cell for k in used):return 'CONFLICTED_SAME_PULSE'
        if features[j]['n_selectable_in_basin']>1 and not features[j]['dominant_in_basin']:return 'PHASE_AMBIGUOUS_LOCAL_EXTREMUM'
        return None
    selectable={j for j,e in enumerate(events) if e['selectable'] and j in windows}
    next_context=max((s['segment_id'] for s in original_segments),default=0)
    while True:
        anchors=[j for j in selectable-used-anchor_tried if conflict(j) is None]
        if not anchors:break
        anchor=min(anchors,key=lambda j:(abs(events[j]['raw_time_s']-(start+end)/2),-events[j]['prominence'],j))
        anchor_tried.add(anchor);next_context+=1;context=next_context;connected=False
        for direction in (-1,1):
            current=anchor;template_id=anchor;ref=period;recent=[];previous_high=True
            while not track['search_limit']:
                if any(direction*(k-current)>0 for k in neighbors(current)):break
                prev=events[current];barriers=[j for j in used if direction*(j-current)>0]
                barrier=min(barriers,key=lambda j:abs(j-current)) if barriers else None
                choices=[];observed=[]
                for j in sorted(selectable):
                    delta=direction*(events[j]['raw_time_s']-prev['raw_time_s'])
                    if not config['search_min_period_ratio']*ref<=delta<=config['search_max_period_ratio']*ref:continue
                    why=None
                    if barrier is not None and direction*(j-barrier)>0:why='TRACK_CONFLICT'
                    elif j not in used:why=conflict(j)
                    elif any(direction*(k-j)<0 for k in neighbors(j)):why='TRACK_CONFLICT'
                    if j==current:continue
                    rho=float(np.dot(windows[template_id],windows[j]));p=events[j]['normalized_prominence'];pp=prev['normalized_prominence']
                    sim=float(np.exp(-abs(np.log((p+1e-12)/(pp+1e-12)))));score=config['correlation_weight']*rho+config['prominence_weight']*sim
                    safety=[]
                    if features[j]['n_selectable_in_basin']>1 and not features[j]['dominant_in_basin']:safety.append('PHASE_AMBIGUOUS_LOCAL_EXTREMUM')
                    if delta<config['search_min_period_ratio']*period or delta>config['search_max_period_ratio']*period:safety.append('ACF_PERIOD_CONFLICT')
                    if features[j]['pulse_cell']==features[current]['pulse_cell'] and all(v is not None for v in features[j]['pulse_cell']):safety.append('SAME_PULSE_PHASE_CONFLICT')
                    omitted=[k for k in selectable if min(current,j)<k<max(current,j) and features[k]['dominant_in_basin'] and
                             features[k]['pulse_cell'] not in (features[current]['pulse_cell'],features[j]['pulse_cell'])]
                    if omitted:safety.append('SKIPPED_PULSE_CELL')
                    observed.append((score,rho,j,delta,why,safety))
                    if why is None:choices.append((score,rho,j,delta,safety))
                choices.sort(reverse=True);margin=choices[0][0]-choices[1][0] if len(choices)>1 else None
                best=choices[0] if choices else None;failures=[]
                if not choices:failures.append('TRACK_CONFLICT' if observed else 'NO_AVAILABLE_CANDIDATE_IN_SEARCH_RANGE')
                if best and best[1]<config['min_correlation']:failures.append('LOW_SHAPE_CORR')
                if margin is not None and margin<config['min_score_margin']:failures.append('SCORE_MARGIN_BELOW_MIN')
                if best:failures.extend(best[4])
                attempts_start=len(track['recovery_attempts'])
                for rank,choice in enumerate(choices):
                    score,rho,j,delta,safety=choice;key=(current,j,direction,template_id,round(ref*fs,8))
                    if key in visited:continue
                    visited.add(key);attempt_count+=1
                    row={'segment_id':context,'stage':'RECOVERY','direction':direction,'previous_event':prev['event_id'],
                         'candidate_event':events[j]['event_id'],'template_event':events[template_id]['event_id'],
                         'period_reference_s':ref,'delta_s':delta,'score':score,'rho':rho,'rank':rank+1,'margin':margin,
                         'accepted':not failures and rank==0,'safety_reasons':safety,'context_key':list(key)}
                    track['recovery_attempts'].append(row);track['attempts'].append(row)
                    if attempt_count>=cap:track['search_limit']=True;failures.append('SEARCH_LIMIT');row['accepted']=False;break
                if failures:
                    for attempted in track['recovery_attempts'][attempts_start:]:attempted.update(accepted=False,rejection_reasons=list(failures))
                    bid=f'{method}_recovery_break_{len(recovery_breaks)+1:03d}'
                    row={'break_id':bid,'stage':'RECOVERY','method':method,'segment_id':context,'direction':direction,
                         'in_window':bool(start<=prev['subsample_time_s']<end),'break_time_s':prev['subsample_time_s'],
                         'previous_event':prev['event_id'],'previous_raw_time_s':prev['raw_time_s'],'previous_time_s':prev['subsample_time_s'],
                         'template_event':events[template_id]['event_id'],'period_reference_s':ref,'acf_period_s':period,
                         'predicted_next_time_s':prev['raw_time_s']+direction*ref,
                         'search_start_s':min(prev['raw_time_s']+direction*config['search_min_period_ratio']*ref,prev['raw_time_s']+direction*config['search_max_period_ratio']*ref),
                         'search_end_s':max(prev['raw_time_s']+direction*config['search_min_period_ratio']*ref,prev['raw_time_s']+direction*config['search_max_period_ratio']*ref),
                         'candidate_count':len(choices),'best_event':events[best[2]]['event_id'] if best else None,
                         'best_time_s':events[best[2]]['subsample_time_s'] if best else None,'best_delta_s':best[3] if best else None,
                         'best_rho':best[1] if best else None,'best_score':best[0] if best else None,'score_margin':margin,
                         'runner_up_event':events[choices[1][2]]['event_id'] if len(choices)>1 else None,'runner_up_score':choices[1][0] if len(choices)>1 else None,
                         'best_prominence':events[best[2]]['prominence'] if best else None,
                         'best_normalized_prominence':events[best[2]]['normalized_prominence'] if best else None,
                         'minimum_rho':config['min_correlation'],'minimum_score_margin':config['min_score_margin'],
                         'direct_reasons':failures,'tracking_reason':'LOW_SHAPE_CORR' if 'LOW_SHAPE_CORR' in failures else failures[0]}
                    recovery_breaks.append(row)
                    for score,rho,j,delta,why,safety in observed:
                        recovery_candidates.append({'break_id':bid,'method':method,'event_id':events[j]['event_id'],
                            'candidate_time_s':events[j]['subsample_time_s'],'candidate_raw_time_s':events[j]['raw_time_s'],
                            'delta_s':delta,'delta_over_reference':delta/ref,'in_search_range':True,'available_at_break':why is None,
                            'availability_reason':why or 'AVAILABLE','excluded_by_segment':None,'rho':rho,'score':score,
                            'rank_in_actual_choices':next((k+1 for k,c in enumerate(choices) if c[2]==j),None),
                            'prominence':events[j]['prominence'],'normalized_prominence':events[j]['normalized_prominence'],'safety_reasons':safety})
                    state_reasons.setdefault(current,set()).update(failures);break
                score,rho,j,delta,safety=best;was_used=j in used
                edges.add(tuple(sorted((current,j))));used.update((current,j));connected=True
                track['recovered_connections'].append({'method':method,'event_start':events[min(current,j)]['event_id'],
                    'event_end':events[max(current,j)]['event_id'],'start_s':events[min(current,j)]['subsample_time_s'],
                    'end_s':events[max(current,j)]['subsample_time_s'],'rho':rho,'score':score,'score_margin':margin,
                    'candidate_rank':1,'template_event':events[template_id]['event_id'],'period_reference_s':ref,
                    'interval_over_acf_period':delta/period,'stage':'RECOVERY','safety_reasons':safety,
                    'review_status':'PENDING_MANUAL_VALIDATION'})
                events[current].update(accepted=True,reason='RECOVERED_ANCHOR' if current==anchor and current not in original_used else events[current]['reason'])
                if j not in original_used:events[j].update(accepted=True,reason='RECOVERED_HIGH_CONFIDENCE' if rho>=config['high_confidence_correlation'] else 'RECOVERED_LOW_CONFIDENCE',correlation=rho,match_score=score,anchor=False)
                if anchor not in original_used:events[anchor]['anchor']=True
                current=j
                if was_used:break
                if rho>=config['high_confidence_correlation']:
                    if previous_high:recent=(recent+[delta])[-config['reference_period_median_count']:];ref=float(np.median(recent))
                    else:recent=[]
                    template_id=j;previous_high=True
                else:recent=[];previous_high=False
            if track['search_limit']:break
        if not connected:state_reasons.setdefault(anchor,set()).add('REJECTED_FOR_CONTEXT')
        if track['search_limit']:break
    # Connected components of the actual accepted edges, with unique event IDs.
    components=[];remaining=set(used)
    while remaining:
        stack=[min(remaining)];component=set()
        while stack:
            j=stack.pop()
            if j in component:continue
            component.add(j);stack.extend(k for k in neighbors(j) if k not in component)
        remaining-=component;components.append(sorted(component))
    components.sort(key=lambda c:c[0]);track['segments']=[]
    for sid,ids in enumerate(components,1):
        track['segments'].append({'segment_id':sid,'event_indices':ids})
        for j in ids:events[j]['segment_id']=sid
    for row in seed_breaks:
        j=next(k for k,e in enumerate(events) if e['event_id']==row['previous_event']);row['source_segment_id']=row['segment_id'];row['segment_id']=events[j]['segment_id']
        row['resolved_by_recovery']=any(row['direction']*(k-j)>0 for k in neighbors(j))
    for row in recovery_breaks:
        j=next(k for k,e in enumerate(events) if e['event_id']==row['previous_event']);e=events[j]
        row['context_segment_id']=row['segment_id'];row['segment_id']=e['segment_id'];row['resolved_by_recovery']=any(row['direction']*(k-j)>0 for k in neighbors(j))
    for row in track['breaks']:
        j=next(k for k,e in enumerate(events) if e['event_id']==row['previous_event'])
        row['stage']='LEGACY_SEED';row['resolved_by_recovery']=any(row['direction']*(k-j)>0 for k in neighbors(j))
    track['_diagnostic']['breaks']=seed_breaks+recovery_breaks;track['_diagnostic']['candidates'].extend(recovery_candidates)
    track['candidate_states']=[]
    for j,e in enumerate(events):
        why=conflict(j) if j not in used else None
        state='USED_BY_ACCEPTED_TRACK' if j in used else 'INSUFFICIENT_TEMPLATE_EDGE' if j not in windows else 'MIN_DISTANCE_SUPPRESSED' if not e['selectable'] else why or ('REJECTED_FOR_CONTEXT' if j in anchor_tried else 'UNATTEMPTED')
        track['candidate_states'].append({'event_id':e['event_id'],'method':method,'time_s':e['subsample_time_s'],'state':state,
            'was_legacy_adopted':j in original_used,'anchor_retried':j in anchor_tried,'context_failures':sorted(state_reasons.get(j,set())),**features[j]})
    for row in recovery_breaks:track['breaks'].append({'segment_id':row['segment_id'],'direction':row['direction'],'previous_event':row['previous_event'],
        'time_s':row['previous_raw_time_s'],'reason':row['tracking_reason'],'period_reference_s':row['period_reference_s'],'stage':'RECOVERY','resolved_by_recovery':row['resolved_by_recovery']})
    track['anchor_contexts_tried']=len(anchor_tried);track['context_attempt_count']=attempt_count
    return track


def pooled_metrics_v2(track,membership,start,end,interpolated=True):
    key='subsample_time_s' if interpolated else 'raw_time_s';details={};all_terms=[]
    for width,name in [(2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')]:
        terms=[];periods=[];identities=[];spans=[]
        for segment in track['segments']:
            sid=segment['segment_id'];ids=membership.get(sid,[]);times=[track['events'][j][key] for j in ids];v=np.diff(times)
            if len(v)<width:continue
            periods.extend(v);spans.append((times[0],times[-1]));identities.extend(f'{sid}:{i}' for i in range(len(v)))
            values=np.abs(np.diff(v)) if width==2 else np.abs(v[width//2:len(v)-width//2]-np.convolve(v,np.ones(width)/width,'valid'))
            terms.extend(values)
            for i,value in enumerate(values):all_terms.append({'method':track['method'],'metric':name,'version':'pooled_within_segments_v2','segment_id':sid,
                'first_period_id':i,'last_period_id':i+width-1,'event_start':track['events'][ids[i]]['event_id'],
                'event_end':track['events'][ids[i+width]]['event_id'],'start_s':times[i],'end_s':times[i+width],'numerator_term_s':float(value)})
        mean=float(np.mean(periods)) if periods else None
        details[name]={'value':float(100*np.mean(terms)/mean) if terms and mean else None,'mean_period_s':mean,'n_cycles':len(periods),
                       'period_cv_percent':period_metrics(periods)['period_cv_percent'],
                       'n_terms':len(terms),'sum_numerator_terms_s':float(sum(terms)),'covered_duration_s':float(sum(b-a for a,b in spans)),
                       'covered_time_fraction':float(sum(b-a for a,b in spans)/(end-start)),'period_ids':identities,
                       'eligible_segment_min_periods':width,'spans':spans}
    return {'version':'pooled_within_segments_v2','scope':'tracked portions only; each metric normalizes by periods in segments contributing its numerator terms',
            **{name:value['value'] for name,value in details.items()},'metric_details':details,'_terms':all_terms}


def summarize_track(track,start,end,period,min_cycles,fs,interpolated=True,frozen=None,force_main_segment=None):
    events=track['events'];key='subsample_time_s' if interpolated else 'raw_time_s';segment_rows=[];period_rows=[];membership={}
    all_periods=[];covered=0.;longest=None
    for segment in track['segments']:
        sid=segment['segment_id'];ids=(frozen.get(sid,[]) if frozen is not None else
             [i for i in segment['event_indices'] if start<=events[i][key]<end])
        membership[sid]=ids;times=np.array([events[i][key] for i in ids]);p=np.diff(times);metrics=period_metrics(p)
        length=float(times[-1]-times[0]) if len(times)>1 else 0.;covered+=length;all_periods.extend(p.tolist())
        row={'segment_id':sid,'start_s':float(times[0]) if len(times) else None,'end_s':float(times[-1]) if len(times) else None,
             'duration_s':length,'event_ids':[events[i]['event_id'] for i in ids],**metrics}
        segment_rows.append(row)
        if len(p)>0 and (longest is None or length>longest['duration_s']):longest=row
        median=float(np.median(p)) if len(p) else None
        for k,value in enumerate(p):
            anomaly=bool(period and (value<.65*period or value>1.45*period))
            period_rows.append({'method':track['method'],'segment_id':sid,'period_id':k,'event_start':events[ids[k]]['event_id'],
                                'event_end':events[ids[k+1]]['event_id'],'start_s':times[k],'end_s':times[k+1],
                                'period_s':value,'delta_previous_period_s':value-p[k-1] if k else None,
                                'flags':['PERIOD_JUMP'] if anomaly else [],'missing_event_possible':bool(period and value>1.45*period),
                                'selected_main':False})
    if force_main_segment is not None:
        longest=next((r for r in segment_rows if r['segment_id']==force_main_segment),None)
    main=longest or {'segment_id':None,'duration_s':0.,'event_ids':[],**period_metrics([])}
    for row in period_rows:row['selected_main']=row['segment_id']==main['segment_id']
    accepted=[e for e in events if e['accepted'] and start<=e[key]<end]
    correlations=[e['correlation'] for e in accepted if not e['anchor'] and e['correlation'] is not None]
    selectable=[e for e in events if e['selectable'] and start<=e[key]<end]
    ambiguous=[e for e in selectable if e['reason']=='AMBIGUOUS_MATCH']
    flags=[];coverage=covered/(end-start);medcorr=float(np.median(correlations)) if correlations else None
    consistency=abs(main['median_period_s']/period-1) if period and main['median_period_s'] else None
    if main['n_cycles']<min_cycles:flags.append('INSUFFICIENT_RESEARCH_CYCLES')
    if coverage<QC_CONFIG['coverage_min']:flags.append('LOW_COVERED_TIME')
    if medcorr is None or medcorr<QC_CONFIG['median_correlation_min']:flags.append('LOW_MATCH_CORRELATION')
    ambiguity=len(ambiguous)/len(selectable) if selectable else 0.
    if ambiguity>QC_CONFIG['ambiguous_fraction_max']:flags.append('AMBIGUOUS_MATCH_FRACTION')
    if consistency is None or consistency>QC_CONFIG['acf_relative_difference_max']:flags.append('ACF_EVENT_MISMATCH')
    if any(row['flags'] for row in period_rows):flags.append('PERIOD_JUMP')
    if any(e['reason']=='MULTIPLE_PEAK_UNUSED' and e['in_window'] for e in events):flags.append('DOUBLE_PEAKS')
    # Missing-candidate breaks are gaps, not automatically a physically missing pulse.
    if any(start<=b['time_s']<end and b['reason']=='MISSING_EVENT_CANDIDATE' and not b.get('resolved_by_recovery') for b in track['breaks']):flags.append('MISSING_EVENTS_POSSIBLE')
    if track.get('search_limit'):flags.append('SEARCH_LIMIT')
    qc_fail=[f for f in flags if f!='DOUBLE_PEAKS']
    expected=(end-start)/period if period else None
    pooled=period_metrics(all_periods)
    # Recompute differences exclusively within each segment; never bridge breaks.
    for metric,width in [('jitter_local_percent',2),('rap_percent',3),('ppq5_percent',5)]:
        terms=[]
        for row in segment_rows:
            ids=membership[row['segment_id']];v=np.diff([events[i][key] for i in ids]);n=len(v)
            if width==2 and n>=2:terms.extend(np.abs(np.diff(v)).tolist())
            if width==3 and n>=3:terms.extend(np.abs(v[1:-1]-np.convolve(v,np.ones(3)/3,'valid')).tolist())
            if width==5 and n>=5:terms.extend(np.abs(v[2:-2]-np.convolve(v,np.ones(5)/5,'valid')).tolist())
        pooled[metric]=float(100*np.mean(terms)/pooled['mean_period_s']) if terms and pooled['mean_period_s'] else None
    v2=pooled_metrics_v2(track,membership,start,end,interpolated)
    v2['status']='PENDING_MANUAL_VALIDATION' if v2['jitter_local_percent'] is not None else 'UNDEFINED'
    for row in period_rows:
        n=next(r['n_cycles'] for r in segment_rows if r['segment_id']==row['segment_id'])
        row.update(used_for_jitter_v2=n>=2,used_for_RAP_v2=n>=3,used_for_PPQ5_v2=n>=5)
    return {'main':main,'main_scope':'longest accepted contiguous segment; not whole-window jitter',
            'longest_segment':main,'pooled_within_segments_v2':v2,'research_primary_candidate':'pooled_within_segments_v2',
            'manual_validation_status':'PENDING_MANUAL_ANNOTATION',
            'segments':segment_rows,'pooled_within_segment_metrics':pooled,
            'n_complete_periods_all_segments':len(all_periods),'expected_periods_proxy':expected,
            'coverage_proxy':len(all_periods)/expected if expected else None,'covered_time_fraction':coverage,
            'coverage_definition':'sum of disjoint completed accepted segment spans / requested window duration',
            'median_match_correlation':medcorr,'ambiguous_event_count':len(ambiguous),'ambiguity_denominator_candidates':len(selectable),
            'ambiguous_event_fraction':ambiguity,'acf_consistency_relative_difference':consistency,
            'timestamp_resolution_us':1e6/fs,'median_samples_per_period':main['median_period_s']*fs if main['median_period_s'] else None,
            'n_unused_candidates':sum(not e['accepted'] for e in selectable),
            'n_short_periods':sum(bool(period and r['period_s']<.65*period) for r in period_rows),
            'n_long_periods':sum(bool(period and r['period_s']>1.45*period) for r in period_rows),
            'n_duplicate_adopted_samples':len(accepted)-len({e['sample_id'] for e in accepted}),
            'status':'QC_PASS' if not qc_fail else ('UNDEFINED' if main['jitter_local_percent'] is None else 'PROVISIONAL_QC_FAIL'),
            'flags':flags,'_membership':membership,'_periods':period_rows}


def match_events(a,b,tolerance):
    """One-to-one correspondence, no ground truth implied."""
    a=sorted(a);b=sorted(b);i=j=matched=0;errors=[]
    while i<len(a) and j<len(b):
        d=b[j]-a[i]
        if abs(d)<=tolerance:matched+=1;errors.append(d);i+=1;j+=1
        elif d<0:j+=1
        else:i+=1
    return {'matched':matched,'n_primary':len(a),'n_comparison':len(b),
            'correspondence_fraction':2*matched/(len(a)+len(b)) if len(a)+len(b) else None,
            'tolerance_s':tolerance,'errors_s':errors}


def local_acf_rows(t,x,start,end,fs,config):
    size=max(3,int(round(config['acf_window_sec']*fs)));hop=max(1,int(round(config['acf_hop_sec']*fs)))
    ids=np.flatnonzero((t>=start)&(t<end));rows=[]
    if len(ids)<size:return rows
    for a in range(ids[0],ids[-1]-size+2,hop):
        p=pitch_candidates(x[a:a+size],fs,config)
        rows.append({'time_s':float(t[a]+size/(2*fs)),'period_s':p['period_s'],'f0_hz':1/p['period_s'] if p['period_s'] else None,'flags':p['flags']})
    return rows


def cpp_analysis(x,fs,start,config,acf_period=None):
    target=config['target_fs_hz'];rows=[];examples={};resampling={}
    if fs<target*(1-1e-12):return {'status':'INVALID_SAMPLING_RATE','flow_CPP_dB':None,'n_frames':0,'n_valid_frames':0,'rows':rows,'_examples':examples}
    if np.ptp(x)<=32*np.finfo(float).eps*max(np.max(np.abs(x)),np.finfo(float).tiny):
        return {'status':'NO_SIGNAL','flow_CPP_dB':None,'n_frames':0,'n_valid_frames':0,'rows':rows,'_examples':examples,
                'flat_signal_check':'raw window is constant before resampling; FIR edge transients are not signal'}
    ratio=Fraction(target/fs).limit_denominator(1000000)
    actual=fs*ratio.numerator/ratio.denominator
    if not math.isclose(actual,target,rel_tol=1e-10):raise ValueError('cannot represent resampling ratio accurately')
    y=signal.resample_poly(x,ratio.numerator,ratio.denominator,window=('kaiser',5.0),padtype='constant') if ratio!=1 else np.asarray(x).copy()
    resampling={'raw_fs_hz':fs,'target_fs_hz':target,'actual_fs_hz':actual,'up':ratio.numerator,'down':ratio.denominator,
                'filter':'resample_poly default FIR, Kaiser beta 5; no additional smoothing','padtype':'constant'}
    length=int(round(config['frame_sec']*target));hop=int(round(config['hop_sec']*target));nfft=config['nfft']
    if length<3 or hop<1 or length>nfft:raise ValueError('CPP frame/hop invalid or frame longer than fixed nfft=2048')
    if len(y)<length:return {'status':'INSUFFICIENT_CPP_FRAMES','flow_CPP_dB':None,'n_frames':0,'n_valid_frames':0,'resampling':resampling,'rows':rows,'_examples':examples}
    q=np.arange(nfft//2+1)/target;positive=q>0
    trend=(q>=config['trend_start_sec'])&(q<=config['trend_end_sec'])
    search=(q>=1/config['pitch_max_hz'])&(q<=1/config['pitch_min_hz'])
    search_ids=np.flatnonzero(search)
    for frame_id,a in enumerate(range(0,len(y)-length+1,hop)):
        raw=y[a:a+length];flags=[]
        row={'frame_id':frame_id,'start_s':start+a/target,'center_s':start+(a+length/2)/target,
             'end_s':start+(a+length)/target,'CPP_frame_dB':None,'peak_q_s':None,'peak_frequency_hz':None,
             'second_peak_q_s':None,'second_peak_dB':None,'interpolation':None,'status':None,'flags':[]}
        if np.ptp(raw)<=32*np.finfo(float).eps*max(np.max(np.abs(raw)),np.finfo(float).tiny):
            row['status']='NO_SIGNAL';rows.append(row);continue
        frame=(raw-raw.mean())*np.hanning(length)
        power=np.abs(np.fft.rfft(frame,n=nfft))**2
        frame_freq=np.fft.rfftfreq(nfft,1/target)
        eps_p=config['spectrum_epsilon_ratio']*np.max(power)
        if eps_p<=0 or not np.isfinite(power).all():row['status']='NO_SIGNAL';rows.append(row);continue
        cep=np.fft.irfft(np.log(np.maximum(power,eps_p)),n=nfft)[:nfft//2+1]
        c2=cep*cep;eps_c=config['cepstrum_epsilon_ratio']*np.max(c2[positive])
        if eps_c<=0 or not np.isfinite(eps_c):row['status']='INVALID_CEPSTRUM';rows.append(row);continue
        db=10*np.log10(np.maximum(c2,eps_c))
        try:
            fit=stats.theilslopes(db[trend],q[trend]);slope=float(fit.slope);intercept=float(fit.intercept)
            if not np.isfinite([slope,intercept]).all():raise ValueError('nonfinite regression')
        except (ValueError,RuntimeError,FloatingPointError) as exc:
            row.update(status='TREND_FIT_FAILED',flags=[str(exc)]);rows.append(row);continue
        background=slope*q+intercept;residual=db-background
        peaks,_=signal.find_peaks(residual);peaks=peaks[search[peaks]]
        if not len(peaks):
            chosen=int(search_ids[np.argmax(residual[search])]);peaks=np.array([chosen]);flags.append('CEPSTRAL_PEAK_EDGE')
        else:chosen=int(peaks[np.argmax(residual[peaks])])
        if chosen in (int(search_ids[0]),int(search_ids[-1])) and 'CEPSTRAL_PEAK_EDGE' not in flags:flags.append('CEPSTRAL_PEAK_EDGE')
        ranked=sorted(peaks,key=lambda k:residual[k],reverse=True)
        offset,value,interp=quadratic(residual,chosen)
        peakq=(chosen+offset)/target
        if not (1/config['pitch_max_hz']<=peakq<=1/config['pitch_min_hz']):
            offset=0.;peakq=chosen/target;value=float(residual[chosen]);interp='SEARCH_EDGE_FALLBACK';flags.append('CEPSTRAL_PEAK_EDGE')
        hz=1/peakq
        if acf_period and abs(hz*acf_period-1)>QC_CONFIG['cpp_f0_mismatch_relative']:flags.append('CPP_F0_MISMATCH')
        dominant_bin_hz=float(frame_freq[np.argmax(power)])
        if dominant_bin_hz<config['pitch_min_hz']:flags.append('LOW_BAND_DOMINANT')
        if config['frame_sec']<3/config['pitch_min_hz']:flags.append('SHORT_FRAME_FOR_LOWEST_F0')
        if interp!='QUADRATIC':flags.append('PEAK_INTERPOLATION_FALLBACK')
        second=ranked[1] if len(ranked)>1 else None
        row.update(CPP_frame_dB=float(value),peak_q_s=peakq,peak_frequency_hz=hz,
                   second_peak_q_s=second/target if second is not None else None,
                   second_peak_dB=float(residual[second]) if second is not None else None,
                   peak_rank=1,n_peak_candidates=len(peaks),all_peak_candidates=[{'q_s':int(k)/target,'residual_dB':float(residual[k])} for k in ranked],
                   background_slope_db_per_s=slope,background_intercept_db=intercept,
                   background_at_peak_dB=slope*peakq+intercept,cepstrum_at_peak_dB=value+slope*peakq+intercept,
                   epsilon_P=float(eps_p),epsilon_C=float(eps_c),interpolation=interp,
                   dominant_spectrum_bin_hz=dominant_bin_hz,
                   below_pitch_floor_power_fraction=float(np.sum(power[frame_freq<config['pitch_min_hz']])/np.sum(power)),
                   status='QC_WARN' if flags else 'OK',flags=flags)
        rows.append(row)
        examples[frame_id]={'q':q,'C_dB':db,'background':background,'D':residual,'peaks':peaks,'row':row,'search':search}
    valid=[r['CPP_frame_dB'] for r in rows if r['CPP_frame_dB'] is not None]
    flags=sorted(set(f for r in rows for f in r['flags']))
    return {'status':('QC_WARN' if flags else 'OK') if valid else 'NO_VALID_CPP_FRAMES',
            'flow_CPP_dB':float(np.mean(valid)) if valid else None,'median_CPP_dB':float(np.median(valid)) if valid else None,
            'std_CPP_dB':float(np.std(valid)) if valid else None,'n_frames':len(rows),'n_valid_frames':len(valid),
            'n_invalid_frames':len(rows)-len(valid),'flags':flags,'resampling':resampling,
            'frame_samples':length,'hop_samples':hop,'actual_frame_sec':length/target,'actual_hop_sec':hop/target,
            'rows':rows,'_examples':examples}


def analyze_arrays(t,x,fs,start,end,tracking,cppconfig,min_cycles):
    inside=(t>=start)&(t<end)
    if np.count_nonzero(inside)<3:raise ValueError('analysis window has fewer than three samples')
    pitch=pitch_candidates(x[inside],fs,tracking)
    primary=float(tracking['primary_prominence_percent']);thresholds=sorted(set(map(float,tracking['prominence_sweep_percent']))|{primary})
    tracks={};summaries={};sensitivity=[]
    for method in ('A','B'):
        for threshold in thresholds:
            key=f'{method}:{threshold:g}';track=extract_track(t,x,fs,start,end,pitch['period_s'],threshold,tracking,method)
            result=summarize_track(track,start,end,pitch['period_s'],min_cycles,fs)
            tracks[key]=track;summaries[key]=result
    mainkey=f'A:{primary:g}';main=summaries[mainkey];basetrack=tracks[mainkey]
    for key,result in summaries.items():
        track=tracks[key];times=[e['raw_time_s'] for e in track['events'] if e['accepted'] and start<=e['subsample_time_s']<end]
        same_method_primary=tracks[f'{track["method"]}:{primary:g}']
        same_method_times=[e['raw_time_s'] for e in same_method_primary['events'] if e['accepted'] and start<=e['subsample_time_s']<end]
        correspondence=match_events(same_method_times,times,QC_CONFIG['event_tolerance_raw_samples']/fs)
        metric=result['main']
        sensitivity.append({'method':track['method'],'prominence_percent':track['prominence_percent'],
                            **{k:metric.get(k) for k in ('jitter_local_percent','rap_percent','ppq5_percent','n_cycles','mean_period_s')},
                            'n_events':len(times),'main_segment_id':metric['segment_id'],'coverage_proxy':result['coverage_proxy'],
                            'n_selectable_candidates':sum(e['selectable'] and e['in_window'] for e in track['events']),
                            'n_unused_candidates':result['n_unused_candidates'],
                            'event_acceptance_fraction':len(times)/sum(e['selectable'] and e['in_window'] for e in track['events']) if any(e['selectable'] and e['in_window'] for e in track['events']) else None,
                            'covered_time_fraction':result['covered_time_fraction'],'acf_consistency':result['acf_consistency_relative_difference'],
                            'status':result['status'],'flags':list(result['flags']),'event_correspondence_to_primary':correspondence['correspondence_fraction'],
                            'event_match_tolerance_s':1/fs})
        sensitivity[-1]['status_basis']='structural QC before cross-threshold sensitivity'
        sensitivity[-1]['correspondence_basis']='same method at fixed primary prominence; A/B phase difference not treated as threshold mismatch'
    valid=[r for r in sensitivity if r['method']=='A' and r['status']=='QC_PASS' and r['jitter_local_percent'] is not None]
    values=np.array([r['jitter_local_percent'] for r in valid]);value=main['main']['jitter_local_percent']
    sensitivity_summary={'n_thresholds':len(thresholds),'n_valid_thresholds':len(values),'n_excluded_thresholds':len(thresholds)-len(values),
                         'validity_rule':'pre-sensitivity structural QC_PASS; all excluded values remain in CSV',
                         'median':float(np.median(values)) if len(values) else None,'min':float(np.min(values)) if len(values) else None,
                         'max':float(np.max(values)) if len(values) else None,'IQR':float(np.percentile(values,75)-np.percentile(values,25)) if len(values) else None,
                         'CV_percent':float(100*np.std(values)/np.mean(values)) if len(values) and np.mean(values)!=0 else None,
                         'interpretation':'method sensitivity, not a statistical confidence interval'}
    if len(valid)<2:main['flags'].append('INSUFFICIENT_VALID_THRESHOLDS')
    if any(r['event_correspondence_to_primary'] is None or r['event_correspondence_to_primary']<.9 for r in valid):main['flags'].append('THRESHOLD_EVENT_MISMATCH')
    if value is not None and any(abs(r['jitter_local_percent']-value)>max(.2,.2*abs(value)) for r in valid):main['flags'].append('THRESHOLD_JITTER_SENSITIVE')
    if 'PITCH_AMBIGUOUS' in pitch['flags']:main['flags'].append('PITCH_AMBIGUOUS')
    if not pitch['period_s']:main['flags'].extend(pitch['flags'])
    if main['flags'] and any(f!='DOUBLE_PEAKS' for f in main['flags']):
        main['status']='UNDEFINED' if value is None else 'PROVISIONAL_QC_FAIL'
    comparison=[]
    for interpolated in (False,True):
        summary=summarize_track(basetrack,start,end,pitch['period_s'],min_cycles,fs,interpolated=interpolated,frozen=main['_membership'],force_main_segment=main['main']['segment_id'])
        selected=[e for e in basetrack['events'] if e['event_id'] in main['main']['event_ids']]
        offsets=[e['offset_us'] for e in selected]
        comparison.append({'subsample':'quadratic' if interpolated else 'off',**{k:summary['main'].get(k) for k in ('jitter_local_percent','rap_percent','ppq5_percent','mean_period_s','n_cycles')},
                           'pooled_v2':summary['pooled_within_segments_v2'],
                           'event_identity_changed':False,'event_ids':main['main']['event_ids'],'offset_mean_us':float(np.mean(offsets)) if offsets and interpolated else 0.,
                           'offset_max_abs_us':float(np.max(np.abs(offsets))) if offsets and interpolated else 0.,
                           'n_interpolation_fallback':sum(e['interpolation_status']!='QUADRATIC' for e in selected),
                           'boundary_membership_changed_without_freezing':any((start<=e['raw_time_s']<end)!=(start<=e['subsample_time_s']<end) for e in selected),
                           'timestamp_resolution_us':1e6/fs,
                           'boundary_identity_note':'same ON-selected event IDs for OFF; no re-detection or period filtering'})
    # Control B uses its own events; never fills missing A events.
    control=summaries[f'B:{primary:g}'];alignment=[]
    tb=np.array([e['subsample_time_s'] for e in tracks[f'B:{primary:g}']['events'] if e['accepted'] and e['in_window']])
    for i,cycle in enumerate(main['_periods']):
        a0,a1=cycle['start_s'],cycle['end_s'];b=tb[(tb>=a0)&(tb<a1)]
        alignment.append({'A_cycle_id':i,'segment_id':cycle['segment_id'],'A_start_s':a0,'A_end_s':a1,'n_B_events':len(b),'B_times_s':b.tolist(),
                          'status':'ONE_CONTROL_EVENT' if len(b)==1 else 'PHASE_EVENT_DISAGREEMENT'})
    cpp=cpp_analysis(x[inside],fs,float(t[inside][0]),cppconfig,pitch['period_s'])
    return {'pitch':pitch,'primary':main,'control':control,'sensitivity':sensitivity,'sensitivity_summary':sensitivity_summary,
            'subsample_comparison':comparison,'event_alignment':alignment,'cpp':cpp,
            'local_acf':local_acf_rows(t,x,start,end,fs,tracking),'_tracks':tracks,'_summaries':summaries,
            '_primary_key':mainkey,'_control_key':f'B:{primary:g}'}


def sampling_comparison(data,reference,start,end,tracking,cppconfig,min_cycles,primary):
    if reference is None:return [{'status':'NOT_AVAILABLE','reason':'no actual stride-1 reference supplied; stride-5 samples cannot reconstruct it'}]
    if reference['stride']!=1:raise ValueError('--reference-fullrate must have actual step stride 1')
    if not math.isclose(data['dt_s'],reference['dt_s'],rel_tol=1e-12):raise ValueError('reference solver dt differs')
    offset=data['steps']-reference['steps'][0]
    if np.any(offset<0) or np.any(offset>=len(reference['x'])):raise ValueError('full-rate reference does not cover input steps')
    expected=reference['x'][offset]
    scale=max(np.ptp(data['x']),np.max(np.abs(data['x'])),np.finfo(float).tiny)
    if not np.allclose(expected,data['x'],rtol=1e-9,atol=scale*1e-10):raise ValueError('reference and input differ at matching steps; not a sampling-only comparison')
    rows=[]
    original_ids=primary['primary']['main']['event_ids']
    original=[e['raw_time_s'] for e in primary['_tracks'][primary['_primary_key']]['events'] if e['event_id'] in original_ids]
    for phase in (-1,0,1,2,3,4):
        use=np.ones(len(reference['x']),dtype=bool) if phase==-1 else reference['steps']%5==phase
        t=reference['t'][use];x=reference['x'][use];fs=reference['fs'] if phase==-1 else reference['fs']/5
        result=analyze_arrays(t,x,fs,start,end,tracking,cppconfig,min_cycles)
        for sub in result['subsample_comparison']:
            adopted=result['primary']['main']['event_ids']
            other=[e['raw_time_s'] for e in result['_tracks'][result['_primary_key']]['events'] if e['event_id'] in adopted]
            match=match_events(original,other,max(1/data['fs'],1/fs))
            rows.append({'status':'ANALYZED','source':'actual fullrate' if phase==-1 else 'decimated actual fullrate',
                         'step_stride':1 if phase==-1 else 5,'phase':None if phase==-1 else phase,'fs_hz':fs,
                         **sub,'flow_CPP_dB':result['cpp']['flow_CPP_dB'],'jitter_qc':result['primary']['status'],
                         'coverage_proxy':result['primary']['coverage_proxy'],'f0_hz':result['primary']['main']['median_f0_hz'],
                         'event_correspondence_to_input':match['correspondence_fraction'],
                         'reference_match_note':'same step values verified; no separate model provenance inferred from DAT alone'})
    return rows


def save_plots(out,data,start,end,result,zoom_start,zoom_duration,tracking,sampling):
    t=data['t'];x=data['x'];fs=data['fs'];primary=result['primary'];pitch=result['pitch'];cpp=result['cpp']
    events=result['_tracks'][result['_primary_key']]['events'];control=result['_tracks'][result['_control_key']]['events']
    accepted=[e for e in events if e['accepted']];rejected=[e for e in events if not e['accepted']]
    fig,ax=plt.subplots(figsize=(11,4));ax.plot(t,x*1e6,lw=.6,color='black');ax.axvspan(start,end,color='orange',alpha=.13,label='analysis window')
    if accepted:ax.scatter([e['subsample_time_s'] for e in accepted],np.interp([e['subsample_time_s'] for e in accepted],t,x)*1e6,s=9,label='A accepted')
    ax.set(xlabel='Time [s]',ylabel='Ug [mL/s]',title='Whole saved flow; initial flat/ramp section is not removed automatically');ax.legend();fig.tight_layout();fig.savefig(out/'flow_overview.png',dpi=160);plt.close(fig)
    mainids=set(primary['main']['event_ids']);main=[e for e in events if e['event_id'] in mainids];main.sort(key=lambda e:e['subsample_time_s'])
    if zoom_start is None:zoom_start=max(start,((main[0]['subsample_time_s']+main[-1]['subsample_time_s'])/2-zoom_duration/2) if main else (start+end-zoom_duration)/2)
    zoom_end=min(end,zoom_start+zoom_duration);mask=(t>=zoom_start)&(t<zoom_end)
    fig,ax=plt.subplots(figsize=(12,5));ax.plot(t[mask],x[mask]*1e6,color='black',lw=1)
    for es,marker,label,color in [(accepted,'o','A accepted','tab:blue'),(rejected,'x','unused / rejected candidate','gray'),([e for e in control if e['accepted']],'^','B control minima','tab:green')]:
        use=[e for e in es if zoom_start<=e['subsample_time_s']<zoom_end]
        if use:at=[e['subsample_time_s'] for e in use];ax.scatter(at,np.interp(at,t,x)*1e6,s=32,marker=marker,label=label,color=color,zorder=4)
    visible=[(i,e) for i,e in enumerate(main) if zoom_start<=e['subsample_time_s']<zoom_end]
    for i,e in visible:ax.annotate(str(i+1),(e['subsample_time_s'],np.interp(e['subsample_time_s'],t,x)*1e6),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
    for (i,a),(j,b) in zip(visible,visible[1:]):
        if j==i+1:
            xpos=(a['subsample_time_s']+b['subsample_time_s'])/2
            ax.text(xpos,.92,f'{1e3*(b["subsample_time_s"]-a["subsample_time_s"]):.2f} ms',transform=ax.get_xaxis_transform(),ha='center',fontsize=7,rotation=20)
    if len(visible)>=6:
        first,last=visible[0][1]['subsample_time_s'],visible[5][1]['subsample_time_s']
        ax.axvspan(first,last,color='tab:blue',alpha=.08,label='5 consecutive main periods (PPQ5 example)')
    if pitch['period_s']:
        q=pitch['period_s'];ax.annotate('',xy=(zoom_start+q,.07),xytext=(zoom_start,.07),xycoords=ax.get_xaxis_transform(),arrowprops={'arrowstyle':'<->','color':'purple'});ax.text(zoom_start+q/2,.02,'ACF T0',transform=ax.get_xaxis_transform(),ha='center',color='purple',fontsize=8)
    ax.set(xlabel='Time [s]',ylabel='Ug [mL/s]',title=f'{Path(data["path"]).name}\nCandidate peaks and tracked periods — {primary["status"]}')
    if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(0,-.18),ncol=2)
    ax.grid(alpha=.15);fig.tight_layout();fig.savefig(out/'flow_events_zoom.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4));acf=pitch['_acf'];lags=np.arange(len(acf))/fs
    limit=min((end-start)/2,2/tracking['pitch_min_hz']);use=lags<=limit;ax.plot(lags[use]*1000,acf[use])
    for row in pitch['candidates']:ax.scatter(row['period_s']*1000,row['acf'],s=20)
    colors=['tab:orange','tab:red','tab:purple']
    for alt,color in zip(pitch['alternatives'],colors):
        q=alt['target_period_s'];ax.axvline(q*1000,color=color,ls='--',label=f'{alt["factor"]:g} T0: {1/q:.1f} Hz'+(' (outside search)' if not alt['inside_search_range'] else ''))
    ax.axvspan(1000/tracking['pitch_max_hz'],1000/tracking['pitch_min_hz'],alpha=.08,color='green');ax.set(xlabel='Lag [ms]',ylabel='Overlap-normalized Pearson ACF',title='Period proposals, not event times')
    if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    fig.tight_layout()
    for name in ('acf_period_candidates.png','autocorrelation_period_candidates.png'):fig.savefig(out/name,dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5));examples=cpp['_examples']
    if examples:
        chosen=sorted(examples)[len(examples)//2];ex=examples[chosen];row=ex['row'];q=ex['q'];use=(q>0)&(q<=.022)
        ax.plot(q[use]*1000,ex['C_dB'][use],label='10log10(c²)');ax.plot(q[use]*1000,ex['background'][use],label='Theil–Sen background B(q)',color='tab:orange')
        ax.axvspan(2.5,1000/60,color='green',alpha=.08,label='CPP search (60–400 Hz)')
        pk=ex['peaks'];ax.scatter(q[pk]*1000,ex['C_dB'][pk],marker='x',color='gray',label='local residual-peak candidates')
        qp=row['peak_q_s'];bp=row['background_at_peak_dB'];cp=row['cepstrum_at_peak_dB']
        ax.scatter([qp*1000],[cp],color='red',zorder=5);ax.annotate('',xy=(qp*1000,cp),xytext=(qp*1000,bp),arrowprops={'arrowstyle':'<->','color':'red'})
        ax.annotate(f'CPP={row["CPP_frame_dB"]:.2f} dB\nq={qp*1000:.3f} ms, 1/q={1/qp:.1f} Hz',
                    (qp*1000,cp),xytext=(.48,.78),textcoords='axes fraction',fontsize=9,
                    bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'},arrowprops={'arrowstyle':'-','color':'red'})
        ax.margins(y=.12)
        ax.set_title(f'{Path(data["path"]).name}\nFrame {chosen}, t={row["center_s"]:.4f} s; CPP is vertical difference at same q')
        ax.legend(fontsize=8)
    else:ax.text(.5,.5,cpp['status'],transform=ax.transAxes,ha='center');ax.set_title('No valid CPP frame; no fabricated cepstrum')
    ax.set(xlabel='Quefrency [ms]',ylabel='Power cepstrum [dB]');fig.tight_layout();fig.savefig(out/'cepstrum_example.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True);rows=cpp['rows']
    if rows:
        axes[0].plot([r['center_s'] for r in rows],[r['CPP_frame_dB'] for r in rows],'.-');axes[1].plot([r['center_s'] for r in rows],[r['peak_frequency_hz'] for r in rows],'.-')
    else:axes[0].text(.5,.5,cpp['status'],transform=axes[0].transAxes,ha='center')
    axes[0].set_ylabel('flow frame CPP [dB]');axes[1].set(ylabel='CPP peak 1/q [Hz]',xlabel='Frame time [s]');fig.tight_layout();fig.savefig(out/'cpp_frames.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
    for method,style in [('A','-'),('B','--')]:
        rows=[r for r in result['sensitivity'] if r['method']==method];p=[r['prominence_percent'] for r in rows]
        for metric,marker in [('jitter_local_percent','o'),('rap_percent','s'),('ppq5_percent','^')]:axes[0].plot(p,[r[metric] for r in rows],style,marker=marker,label=method+' '+metric)
        axes[1].plot(p,[r['n_cycles'] for r in rows],style,marker='o',label=method)
        axes[2].plot(p,[r['covered_time_fraction'] for r in rows],style,marker='o',label=method)
        for r in rows:axes[1].annotate('P' if r['status']=='QC_PASS' else '!',(r['prominence_percent'],r['n_cycles']),fontsize=8)
    axes[0].set_ylabel('Candidate metric [%]');axes[0].legend(fontsize=7,ncol=2);axes[1].set_ylabel('Main segment cycles');axes[1].legend();axes[2].set(xlabel='Prominence [% of P95−P5]',ylabel='Covered time fraction (proxy)');axes[2].legend();fig.tight_layout();fig.savefig(out/'threshold_sensitivity.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4));compar=result['subsample_comparison']
    for row in compar:
        axes[0].plot(range(3),[row[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')],marker='o',label=row['subsample'])
    axes[0].set(xticks=range(3),xticklabels=['Jitter','RAP','PPQ5'],ylabel='Percent',title='Same event IDs, no re-detection');axes[0].legend()
    offsets=[e['offset_us'] for e in events if e['event_id'] in mainids]
    if offsets:axes[1].hist(offsets,bins=20)
    axes[1].set(xlabel='Quadratic − integer time [µs]',ylabel='Main events');fig.tight_layout();fig.savefig(out/'subsample_comparison.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4));analyzed=[r for r in sampling if r['status']=='ANALYZED']
    if analyzed:
        for sub in ('off','quadratic'):
            rows=[r for r in analyzed if r['subsample']==sub];ax.plot(range(len(rows)),[r['jitter_local_percent'] for r in rows],'o-',label=sub)
            ax.set_xticks(range(len(rows)));ax.set_xticklabels(['1 step' if r['step_stride']==1 else f'5 step phase {r["phase"]}' for r in rows],rotation=25,fontsize=8)
        ax.set_ylabel('Jitter local [%]');ax.legend()
    else:ax.text(.5,.5,'NOT_AVAILABLE\nNo actual stride-1 reference. No reconstruction from stride-5.',ha='center',transform=ax.transAxes)
    ax.set_title('Sampling comparison');fig.tight_layout();fig.savefig(out/'sampling_comparison.png',dpi=160);plt.close(fig)



def diagnostic_block_metrics(blocks):
    """Pool numerators inside blocks; never form a difference across blocks."""
    periods=[v for block in blocks for v in block]
    metric=period_metrics(periods);terms={2:[],3:[],5:[]}
    for block in blocks:
        v=np.asarray(block,dtype=float)
        if len(v)>=2:terms[2].extend(np.abs(np.diff(v)).tolist())
        if len(v)>=3:terms[3].extend(np.abs(v[1:-1]-np.convolve(v,np.ones(3)/3,'valid')).tolist())
        if len(v)>=5:terms[5].extend(np.abs(v[2:-2]-np.convolve(v,np.ones(5)/5,'valid')).tolist())
    for width,name in ((2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')):
        metric[name]=float(100*np.mean(terms[width])/metric['mean_period_s']) if terms[width] and metric['mean_period_s'] else None
    metric.update(n_adjacent_period_differences=len(terms[2]),sum_adjacent_abs_difference_s=float(sum(terms[2])),
                  n_RAP_terms=len(terms[3]),n_PPQ5_terms=len(terms[5]),sum_period_s=float(sum(periods)),
                  numerator_pooling='within each contiguous segment only; normalization mean includes all adopted periods')
    return metric


def diagnostic_clip_blocks(track,left,right):
    blocks=[];details=[]
    for segment in track['segments']:
        ev=[track['events'][j] for j in segment['event_indices']]
        chosen=[]
        for a,b in zip(ev,ev[1:]):
            t0,t1=a['subsample_time_s'],b['subsample_time_s']
            if t0>=left-1e-12 and t1<=right+1e-12:
                chosen.append(t1-t0)
                details.append({'method':track['method'],'segment_id':segment['segment_id'],
                    'event_start':a['event_id'],'event_end':b['event_id'],'start_s':t0,'end_s':t1,'period_s':t1-t0})
        if chosen:blocks.append(chosen)
    return blocks,details


def diagnostic_gaps(summary,start,end,method):
    intervals=sorted((r['start_s'],r['end_s'],r['segment_id']) for r in summary['segments'] if r['n_cycles']>0)
    rows=[];cursor=start;previous=None
    for left,right,sid in intervals:
        if left>cursor+1e-12:
            rows.append({'method':method,'start_s':cursor,'end_s':left,'duration_s':left-cursor,
                         'window_fraction':(left-cursor)/(end-start),'previous_segment':previous,'next_segment':sid,
                         'kind':'LEADING_EDGE' if previous is None else 'INTERNAL_TRACKING_GAP'})
        cursor=max(cursor,right);previous=sid
    if cursor<end-1e-12:
        rows.append({'method':method,'start_s':cursor,'end_s':end,'duration_s':end-cursor,
                     'window_fraction':(end-cursor)/(end-start),'previous_segment':previous,'next_segment':None,'kind':'TRAILING_EDGE'})
    return rows


def diagnostic_same_windows(tracks,summaries,start,end):
    windows=[]
    for method in ('A','B'):
        main=summaries[method]['main']
        if main['n_cycles']>0:windows.append((method+'_longest',main['start_s'],main['end_s']))
    for a in summaries['A']['segments']:
        for b in summaries['B']['segments']:
            if not a['n_cycles'] or not b['n_cycles']:continue
            left=max(a['start_s'],b['start_s'],start);right=min(a['end_s'],b['end_s'],end)
            if right>left:windows.append((f'common_A{a["segment_id"]}_B{b["segment_id"]}',left,right))
    rows=[];periods=[]
    for name,left,right in windows:
        for method in ('A','B'):
            blocks,detail=diagnostic_clip_blocks(tracks[method],left,right)
            m=diagnostic_block_metrics(blocks);duration=sum(v['period_s'] for v in detail)
            rows.append({'window_name':name,'method':method,'window_start_s':left,'window_end_s':right,
                         'window_duration_s':right-left,'covered_duration_s':duration,'covered_time_fraction':duration/(right-left),
                         'n_contiguous_blocks':len(blocks),'boundary_rule':'closed comparison interval [left,right]; complete cycles only; no segment bridge',**m})
            periods.extend(dict(r,window_name=name,window_start_s=left,window_end_s=right) for r in detail)
    return rows,periods


def diagnostic_break_details(track,summary,start,end):
    info=track.get('_diagnostic',{});events={e['event_id']:e for e in track['events']};out=[]
    adopted=[e for e in track['events'] if e['accepted']]
    details=info.get('candidates',[])
    for original in info.get('breaks',[]):
        row=dict(original);direction=row['direction'];prev=row['previous_raw_time_s'];sid=row['segment_id']
        candidates=[d for d in details if d['break_id']==row['break_id']]
        next_actual=sorted((d for d in candidates if d['delta_s']>0),key=lambda d:d['delta_s'])
        row['next_candidate_time_s']=row['best_time_s'] if row['best_time_s'] is not None else (next_actual[0]['candidate_time_s'] if next_actual else None)
        row['next_candidate_delta_s']=row['best_delta_s'] if row['best_delta_s'] is not None else (next_actual[0]['delta_s'] if next_actual else None)
        row['n_raw_threshold_candidates_in_search']=sum(d['in_search_range'] for d in candidates)
        row['n_reanchor_blocked_candidates_in_search']=sum(d['in_search_range'] and d['availability_reason']=='REANCHOR_EXCLUSION_REGION' for d in candidates)
        resumed=[e for e in adopted if e['segment_id']!=sid and direction*(e['raw_time_s']-prev)>0]
        resumed.sort(key=lambda e:direction*(e['raw_time_s']-prev))
        row['tracking_resumes_in_direction']=bool(resumed)
        row['resume_event']=resumed[0]['event_id'] if resumed else None
        row['resume_time_s']=resumed[0]['subsample_time_s'] if resumed else None
        row['resume_segment_id']=resumed[0]['segment_id'] if resumed else None
        row['resume_gap_s']=direction*(resumed[0]['subsample_time_s']-row['previous_time_s']) if resumed else None
        for side in ('before','after'):
            es=[e for e in adopted if (e['raw_time_s']<prev if side=='before' else e['raw_time_s']>prev) and e['segment_id']!=sid]
            target=min(es,key=lambda e:abs(e['raw_time_s']-prev)) if es else None
            row['other_segment_'+side+'_time_s']=target['subsample_time_s'] if target else None
        out.append(row)
    return out


def diagnostic_segment_aggregation(track,summary,start,end):
    rows=[];blocks=[]
    for segment in summary['segments']:
        ids=summary['_membership'][segment['segment_id']];times=[track['events'][i]['subsample_time_s'] for i in ids]
        p=np.diff(times);blocks.append(p)
        m=diagnostic_block_metrics([p]);duration=float(sum(p))
        rows.append({'method':track['method'],'scope':'segment','segment_id':segment['segment_id'],
                     'start_s':times[0] if times else None,'end_s':times[-1] if times else None,
                     'covered_duration_s':duration,'window_fraction':duration/(end-start),**m})
    main=next((r for r in rows if r['segment_id']==summary['main']['segment_id']),None)
    if main:rows.append(dict(main,scope='longest_segment'))
    rows.append({'method':track['method'],'scope':'pooled_within_segments','segment_id':None,'start_s':start,'end_s':end,
                 'covered_duration_s':sum(sum(p) for p in blocks),'window_fraction':sum(sum(p) for p in blocks)/(end-start),
                 **diagnostic_block_metrics(blocks)})
    return rows


def diagnostic_plot_break(out,data,start,end,pitch,tracks,summaries,row,context_periods):
    t=data['t'];x=data['x'];period=row['acf_period_s'];center=row['break_time_s']
    left=max(t[0],center-context_periods*period);right=min(t[-1],center+context_periods*period)
    fig,ax=plt.subplots(figsize=(13,5.5));mask=(t>=left)&(t<=right);ax.plot(t[mask],x[mask]*1e6,color='black',lw=.9,label='raw Ug')
    for method,color,marker in [('A','tab:blue','o'),('B','tab:green','v')]:
        tr=tracks[method]
        for accepted in (True,False):
            ev=[e for e in tr['events'] if e['accepted']==accepted and left<=e['subsample_time_s']<=right]
            at=[e['subsample_time_s'] for e in ev]
            if at:ax.scatter(at,np.interp(at,t,x)*1e6,color=color if accepted else 'gray',marker=marker if accepted else 'x',s=30 if accepted else 22,label=method+(' adopted' if accepted else ' unused candidates'),alpha=1 if accepted else .65)
    ax.axvspan(max(left,row['search_start_s']),min(right,row['search_end_s']),color='orange',alpha=.15,label='allowed search interval')
    ax.axvline(row['predicted_next_time_s'],color='red',ls='--',label='predicted next (no snapping)')
    ax.axvline(row['previous_time_s'],color='purple',ls=':',label='previous/reference event')
    rho='none' if row['best_rho'] is None else f'{row["best_rho"]:.4f}'
    margin='none (single/no candidate)' if row['score_margin'] is None else f'{row["score_margin"]:.4f}'
    ax.set_title(f'{row["break_id"]} [{row.get("stage","TRACKING")}, resolved={row.get("resolved_by_recovery",False)}], direction {row["direction"]:+d}, t={center:.6f} s — {" / ".join(row["direct_reasons"])}\n'
                 f'rho={rho} (min {row["minimum_rho"]:.2f}), margin={margin} (min {row["minimum_score_margin"]:.2f}); '
                 f'candidates={row["candidate_count"]}; ACF T0={period*1e3:.3f} ms, Tref={row["period_reference_s"]*1e3:.3f} ms',fontsize=10)
    if row['best_time_s'] is not None:
        ax.scatter([row['best_time_s']],[np.interp(row['best_time_s'],t,x)*1e6],s=65,facecolors='none',edgecolors='red',marker='s',label='rejected best candidate',zorder=5)
        ax.annotate('best candidate',(row['best_time_s'],np.interp(row['best_time_s'],t,x)*1e6),xytext=(12,14),textcoords='offset points',arrowprops={'arrowstyle':'->','color':'red'},fontsize=8)
    if row['resume_time_s'] is not None and left<=row['resume_time_s']<=right:ax.axvline(row['resume_time_s'],color='tab:brown',ls='-.',label=f'resumed segment {row["resume_segment_id"]}')
    ax.set(xlabel='Time [s]',ylabel='Ug [mL/s]',xlim=(left,right));ax.grid(alpha=.15);ax.legend(fontsize=8,ncol=3,loc='upper left',bbox_to_anchor=(0,-.17));fig.tight_layout();fig.savefig(out/(row['break_id']+'.png'),dpi=150);plt.close(fig)


def diagnostic_overview(out,data,start,end,tracks,summaries,gaps):
    fig,axes=plt.subplots(2,1,figsize=(14,8),sharex=True);t=data['t'];mask=(t>=start)&(t<end);x=data['x']
    for ax,method in zip(axes,('A','B')):
        ax.plot(t[mask],x[mask]*1e6,color='black',lw=.6)
        summary=summaries[method];mainid=summary['main']['segment_id']
        for s in summary['segments']:
            if s['n_cycles']>0:ax.axvspan(s['start_s'],s['end_s'],color='tab:blue' if s['segment_id']==mainid else 'tab:green',alpha=.16 if s['segment_id']==mainid else .10)
        for gap in gaps:
            if gap['method']==method:ax.axvspan(gap['start_s'],gap['end_s'],color='tab:red',alpha=.12)
        adopted=[e for e in tracks[method]['events'] if e['accepted'] and e['in_window']]
        ax.scatter([e['subsample_time_s'] for e in adopted],np.interp([e['subsample_time_s'] for e in adopted],t,x)*1e6,s=8,color='tab:blue')
        ax.set(ylabel='Ug [mL/s]',title=f'{method}: blue = longest/main; green = other tracked segments; red = not covered; covered {summary["covered_time_fraction"]:.1%}')
        ax.grid(alpha=.12)
    axes[-1].set_xlabel('Time [s]');fig.tight_layout();fig.savefig(out/'tracking_coverage_overview.png',dpi=180);plt.close(fig)


def diagnostic_new_connection_figures(out,data,period,tracks,variant_tracks,connections):
    """Show every distinct newly connected interval; keep each scenario's scores in CSV."""
    out.mkdir();groups={};index=[];t=data['t'];x=data['x']
    for row in connections:groups.setdefault((row['method'],row['event_start'],row['event_end']),[]).append(row)
    for (method,a,b),rows in groups.items():
        row=rows[0];tr=variant_tracks[(method,row['variant_id'])];center=(row['start_s']+row['end_s'])/2
        left=max(t[0],center-4*period);right=min(t[-1],center+4*period);mask=(t>=left)&(t<=right)
        fig,ax=plt.subplots(figsize=(12,5));ax.plot(t[mask],x[mask]*1e6,color='black',lw=.9,label='raw Ug')
        for original,label,color,marker in [(tracks[method],'strict adopted','tab:blue','o'),(tr,'diagnostic adopted','tab:orange','x'),
                (tracks['B' if method=='A' else 'A'],'strict opposite phase','tab:green','v')]:
            at=[e['subsample_time_s'] for e in original['events'] if e['accepted'] and left<=e['subsample_time_s']<=right]
            if at:ax.scatter(at,np.interp(at,t,x)*1e6,color=color,marker=marker,s=35,label=label)
        ax.axvspan(row['start_s'],row['end_s'],color='orange',alpha=.15,label='newly connected interval')
        for at in (row['start_s'],row['end_s']):ax.axvline(at,color='red',ls='--',lw=.8)
        rho='none' if row['rho'] is None else f'{row["rho"]:.4f}';margin='none (one candidate)' if row['score_margin'] is None else f'{row["score_margin"]:.4f}'
        ax.set_title(f'DIAGNOSTIC ONLY: {method} {a} to {b}; representative {row["variant_id"]}\n'
                     f'rho={rho}; margin={margin}; interval/T0={row["interval_over_acf_period"]:.3f}; '
                     f'opposite-phase events={row["n_strict_opposite_phase_events"]} (not ground truth)',fontsize=10)
        ax.set(xlabel='Time [s]',ylabel='Ug [mL/s]',xlim=(left,right));ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.15);fig.tight_layout()
        name=f'{method}_{a}_{b}.png';fig.savefig(out/name,dpi=130);plt.close(fig)
        index.append({'method':method,'event_start':a,'event_end':b,'figure':name,'representative_variant':row['variant_id'],
                      'all_variants':[r['variant_id'] for r in rows],'review_note':'visual review required; coverage is not validation'})
    write_csv(out.parent/'new_connection_figure_index.csv',index)


def diagnostic_analysis(out,data,start,end,result,tracking,min_cycles,baseline_path=None,context_periods=5):
    """Observation + one-factor diagnostic reruns; main outputs/configuration are preserved."""
    out.mkdir();period=result['pitch']['period_s'];tracks={};summaries={}
    for method in ('A','B'):
        tr=extract_track(data['t'],data['x'],data['fs'],start,end,period,tracking['primary_prominence_percent'],tracking,method,diagnostics=True)
        sm=summarize_track(tr,start,end,period,min_cycles,data['fs'])
        original=result['_tracks'][result['_primary_key'] if method=='A' else result['_control_key']]
        if [(e['event_id'],e['accepted'],e['segment_id']) for e in tr['events']] != [(e['event_id'],e['accepted'],e['segment_id']) for e in original['events']]:raise RuntimeError('diagnostic instrumentation changed strict event selection')
        tracks[method]=tr;summaries[method]=sm
    baseline_check={'instrumented_strict_event_selection_identical':True,'supplied_baseline':str(baseline_path) if baseline_path else None}
    if baseline_path:
        prior=json.loads(baseline_path.read_text());meta=prior['metadata']
        if meta['input_sha256']!=data['sha256'] or meta['start_s']!=start or meta['end_s']!=end or meta['dt_s']!=data['dt_s']:raise ValueError('diagnostic baseline input/window/dt mismatch')
        if prior['tracking_config']!=tracking:raise ValueError('diagnostic baseline tracking settings mismatch')
        for method,key in (('A','primary'),('B','control')):
            old=prior[key];new=result[key]
            if old['main']['event_ids']!=new['main']['event_ids']:raise RuntimeError('baseline main event IDs differ')
            for metric in ('jitter_local_percent','rap_percent','ppq5_percent'):
                a,b=old['main'][metric],new['main'][metric]
                if (a is None)!=(b is None) or (a is not None and not np.isclose(a,b,rtol=1e-10,atol=1e-10)):raise RuntimeError('baseline metric differs: '+method+' '+metric)
        if not np.isclose(prior['cpp']['flow_CPP_dB'],result['cpp']['flow_CPP_dB'],rtol=1e-9,atol=1e-9):raise RuntimeError('CPP differs from existing baseline')
        baseline_check['prior_metrics_and_event_ids_reproduced']=True
        baseline_check['baseline_summary_sha256']=sha256(baseline_path)
    dump_json(out/'strict_baseline_check.json',baseline_check)
    dump_json(out/'strict_results.json',{'tracking_config':tracking,'qc_config':dict(QC_CONFIG,min_cycles=min_cycles),
        'primary':result['primary'],'control':result['control'],'sensitivity':result['sensitivity'],
        'cpp_unchanged_reference':{k:result['cpp'].get(k) for k in ('flow_CPP_dB','status','n_frames','n_valid_frames')},
        'main_definition':'longest continuous segment; QC flags do not remove individual periods',
        'warning_semantics':{'DOUBLE_PEAKS':'at least one unused candidate was removed from future anchor search inside an exclusion region; not proof of true double peaks',
                             'MISSING_EVENTS_POSSIBLE':'no AVAILABLE candidate within the forward/backward search range; raw candidate may exist outside the range, be blocked, or be below prominence'}})
    breaks=[];candidates=[];removals=[];gaps=[];aggregates=[]
    for method in ('A','B'):
        tr=tracks[method];sm=summaries[method]
        breaks.extend(diagnostic_break_details(tr,sm,start,end));candidates.extend(tr['_diagnostic']['candidates']);removals.extend(tr['_diagnostic']['reanchor_removals'])
        gaps.extend(diagnostic_gaps(sm,start,end,method));aggregates.extend(diagnostic_segment_aggregation(tr,sm,start,end))
    raw_details=[]
    for method in ('A','B'):
        y=data['x'] if method=='A' else -data['x'];ix,prop=signal.find_peaks(y,prominence=0)
        amp=tracks[method].get('amplitude_p95_p5',1.);threshold=amp*tracking['primary_prominence_percent']/100
        candidate_samples={e['sample_id'] for e in tracks[method]['events']}
        for row in breaks:
            if row['method']!=method:continue
            ids=[k for k,i in enumerate(ix) if row['search_start_s']<=data['t'][i]<=row['search_end_s']]
            row['n_all_raw_peaks_in_search']=len(ids);row['n_raw_peaks_passing_prominence']=sum(prop['prominences'][k]>=threshold for k in ids)
            row['absence_explanation']='NOT_A_MISSING_CANDIDATE_DECISION'
            if row['candidate_count']==0:
                row['absence_explanation']=('NO_RAW_PEAK_IN_SEARCH_INTERVAL' if not ids else
                    'ALL_RAW_PEAKS_BELOW_PROMINENCE' if row['n_raw_peaks_passing_prominence']==0 else
                    'QUALIFIED_RAW_PEAK_EXISTS_BUT_NOT_AVAILABLE_OR_INVALID_WINDOW')
            if ids:
                best=max(ids,key=lambda k:prop['prominences'][k]);row['strongest_raw_peak_in_search_time_s']=float(data['t'][ix[best]])
                row['strongest_raw_peak_in_search_prominence']=float(prop['prominences'][best])
            else:row['strongest_raw_peak_in_search_time_s']=None;row['strongest_raw_peak_in_search_prominence']=None
            for k in ids:
                raw_details.append({'break_id':row['break_id'],'method':method,'sample_id':int(ix[k]),
                    'raw_time_s':float(data['t'][ix[k]]),'prominence':float(prop['prominences'][k]),
                    'normalized_prominence':float(prop['prominences'][k]/amp) if amp else None,
                    'passes_prominence':bool(prop['prominences'][k]>=threshold),'present_in_tracking_candidates':int(ix[k]) in candidate_samples})
    write_csv(out/'break_diagnostics.csv',breaks);write_csv(out/'break_candidates.csv',candidates);write_csv(out/'reanchor_exclusions.csv',removals)
    write_csv(out/'raw_peaks_at_breaks.csv',raw_details)
    write_csv(out/'uncovered_intervals.csv',gaps);write_csv(out/'segment_jitter_comparison.csv',aggregates)
    same,used=diagnostic_same_windows(tracks,summaries,start,end);write_csv(out/'same_window_metrics.csv',same);write_csv(out/'same_window_periods.csv',used)
    diagnostic_overview(out,data,start,end,tracks,summaries,gaps)
    plots=out/'break_figures';plots.mkdir()
    for row in breaks:
        if row['in_window']:diagnostic_plot_break(plots,data,start,end,result['pitch'],tracks,summaries,row,context_periods)
    # Sweep one factor at a time; unchanged downstream QC is retained in all runs.
    variants=[('strict','baseline',None)]
    variants += [('rho','min_correlation',v) for v in (.50,.55,.60,.65,.70)]
    variants += [('margin','min_score_margin',v) for v in (.00,.05,.10,.15)]
    variants += [('prominence','primary_prominence_percent',v) for v in tracking['prominence_sweep_percent']]
    sensitivity=[];connections=[];variant_tracks={};variant_periods=[];variant_segments=[];variant_events=[]
    for axis,key,value in variants:
        cfg=dict(tracking)
        if value is not None:cfg[key]=value
        variant_id=axis if value is None else f'{axis}_{value:g}'
        for method in ('A','B'):
            tr=tracks[method] if axis=='strict' else extract_track(data['t'],data['x'],data['fs'],start,end,period,cfg['primary_prominence_percent'],cfg,method)
            sm=summaries[method] if axis=='strict' else summarize_track(tr,start,end,period,min_cycles,data['fs'])
            variant_tracks[(method,variant_id)]=tr
            variant_periods.extend(dict(r,method=method,variant_id=variant_id,diagnostic_only=axis!='strict') for r in sm['_periods'])
            variant_segments.extend(dict(r,method=method,variant_id=variant_id,diagnostic_only=axis!='strict') for r in sm['segments'])
            variant_events.extend(dict(r,variant_id=variant_id,diagnostic_only=axis!='strict') for r in tr['events'])
            original=tracks[method];original_pairs={(r['event_start'],r['event_end']) for r in summaries[method]['_periods']}
            accepted_ids={e['event_id'] for e in tr['events'] if e['accepted'] and e['in_window']};strict_ids={e['event_id'] for e in original['events'] if e['accepted'] and e['in_window']}
            newly=[r for r in sm['_periods'] if (r['event_start'],r['event_end']) not in original_pairs]
            attempt={(r['previous_event'],r['candidate_event']):r for r in tr['attempts'] if r['accepted']}
            lookup={e['event_id']:e for e in tr['events']};bad_rho=bad_margin=0
            for r in newly:
                a,b=r['event_start'],r['event_end'];at=attempt.get((a,b)) or attempt.get((b,a));rho=at['rho'] if at else None;margin=at['margin'] if at else None
                bad_rho+=rho is not None and rho<tracking['min_correlation'];bad_margin+=margin is not None and margin<tracking['min_score_margin']
                opposite=tracks['B' if method=='A' else 'A']['events']
                control_count=sum(e['accepted'] and r['start_s']<=e['subsample_time_s']<r['end_s'] for e in opposite)
                connections.append(dict(r,variant_id=variant_id,diagnostic_only=axis!='strict',rho=rho,score_margin=margin,
                    template_event=at['template_event'] if at else None,reference_period_s=at['period_reference_s'] if at else None,
                    score=at['score'] if at else None,
                    below_strict_rho=rho is not None and rho<tracking['min_correlation'],below_strict_margin=margin is not None and margin<tracking['min_score_margin'],
                    interval_over_acf_period=r['period_s']/period if period else None,n_strict_opposite_phase_events=control_count,
                    review_note='new connection versus strict baseline; not automatically an improvement'))
            sensitivity.append({'variant_id':variant_id,'method':method,'axis':axis,'diagnostic_only':axis!='strict',
                'reanchor_policy':cfg.get('reanchor_policy','legacy'),
                'event_correspondence_to_selected_strict':match_events([e['subsample_time_s'] for e in original['events'] if e['accepted'] and e['in_window']],
                    [e['subsample_time_s'] for e in tr['events'] if e['accepted'] and e['in_window']],100e-6)['correspondence_fraction'],
                **{'pooled_v2_'+k:v for k,v in sm['pooled_within_segments_v2'].items() if k in ('jitter_local_percent','rap_percent','ppq5_percent','metric_details')},
                'min_correlation':cfg['min_correlation'],'min_score_margin':cfg['min_score_margin'],'prominence_percent':cfg['primary_prominence_percent'],
                'n_tracked_cycles':sm['n_complete_periods_all_segments'],'main_duration_s':sm['main']['duration_s'],
                'main_window_fraction':sm['main']['duration_s']/(end-start),'covered_time_fraction':sm['covered_time_fraction'],
                'ambiguous_events':sm['ambiguous_event_count'],'n_new_event_ids':len(accepted_ids-strict_ids),'n_removed_event_ids':len(strict_ids-accepted_ids),
                'n_new_connections':len(newly),'n_new_connections_below_strict_rho':int(bad_rho),'n_new_connections_below_strict_margin':int(bad_margin),
                'n_short_periods':sm['n_short_periods'],'n_long_periods':sm['n_long_periods'],
                'jitter_local_percent':sm['main']['jitter_local_percent'],'rap_percent':sm['main']['rap_percent'],'ppq5_percent':sm['main']['ppq5_percent'],
                'n_main_cycles':sm['main']['n_cycles'],'pooled_jitter_percent':sm['pooled_within_segment_metrics']['jitter_local_percent'],
                'status_structural_QC':sm['status'],'flags':sm['flags'],'strict_QC_thresholds_changed':False})
    write_csv(out/'diagnostic_sensitivity.csv',sensitivity);write_csv(out/'new_connections_vs_strict.csv',connections)
    write_csv(out/'sensitivity_periods.csv',variant_periods);write_csv(out/'sensitivity_segments.csv',variant_segments)
    write_csv(out/'sensitivity_events.csv',variant_events)
    diagnostic_new_connection_figures(out/'new_connection_figures',data,period,tracks,variant_tracks,connections)
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for col,axis in enumerate(('rho','margin','prominence')):
        for method,style in [('A','o-'),('B','s--')]:
            rs=[r for r in sensitivity if r['axis']==axis and r['method']==method];x=[r[{'rho':'min_correlation','margin':'min_score_margin','prominence':'prominence_percent'}[axis]] for r in rs]
            axes[0,col].plot(x,[r['covered_time_fraction'] for r in rs],style,label=method)
            axes[1,col].plot(x,[r['jitter_local_percent'] for r in rs],style,label=method)
        axes[0,col].set(title=axis+' (one factor only)',ylabel='Covered time fraction');axes[1,col].set(xlabel=axis,ylabel='Longest-segment Jitter [%]')
        for ax in axes[:,col]:ax.grid(alpha=.15);ax.legend()
    fig.suptitle('Diagnostic tracking thresholds; improved coverage is not proof of correct matching');fig.tight_layout();fig.savefig(out/'diagnostic_threshold_effects.png',dpi=160);plt.close(fig)
    diagnostic_report(out,start,end,result,tracks,summaries,breaks,gaps,aggregates,same,sensitivity)
    return {'strict_baseline_check':baseline_check,'n_breaks_in_window':sum(r['in_window'] for r in breaks),
            'diagnostic_only':True,'diagnostic_sensitivity_runs':len(sensitivity),'cpp_modified':False}


def diagnostic_report(out,start,end,result,tracks,summaries,breaks,gaps,aggregates,same,sensitivity):
    def f(v):return 'undefined' if v is None else f'{v:.8g}'
    lines=['# Tracking/QC diagnostic report','', 'Diagnostics retain the explicitly selected policy '+str(result['_tracks'][result['_primary_key']].get('reanchor_policy','legacy'))+'. Variants are diagnostic-only; no synthetic events, time snapping, or bridging of segment gaps. CPP is not recomputed under altered rules.', '', '## Strict segmentation and selection','']
    for method in ('A','B'):
        s=summaries[method];m=s['main'];pool=s['pooled_within_segment_metrics'];rows=[r for r in breaks if r['method']==method and r['in_window']]
        reasons={why:sum(why in r['direct_reasons'] for r in rows) for why in sorted({v for r in rows for v in r['direct_reasons']})}
        lines += [f"{method}: {s['n_complete_periods_all_segments']} tracked cycles, longest {m['n_cycles']} cycles [{m.get('start_s')}, {m.get('end_s')}], main Jitter={f(m['jitter_local_percent'])}%, pooled-within-segment Jitter={f(pool['jitter_local_percent'])}%.",
                  f"Covered={s['covered_time_fraction']:.6%}; not covered={1-s['covered_time_fraction']:.6%}; longest fraction={m['duration_s']/(end-start):.6%}; interruption decisions by direct reason: {reasons}. Forward/backward observations of a gap are not independent physical failures.", '']
    lines+=['QC flags do not discard tracked periods. The main-value reduction is the deliberate longest-contiguous-segment selection; structural QC labels reliability. DOUBLE_PEAKS is set for unused candidates excluded from future reanchoring, even when no physical double peak has been proven. MISSING_EVENTS_POSSIBLE means no available candidate in the search interval, not necessarily a missing visible pulse. break_candidates.csv separates range, availability, and shape/score failures.','', '## Exact uncovered intervals','']
    if tracks['A'].get('reanchor_policy')=='recover-rejected':lines.insert(4,'Recovery first freezes the legacy seed. LEGACY_SEED interruption/exclusion records are historical, and resolved_by_recovery identifies later connections. Rejected raw candidates remain present and are revisited unless an explicit accepted-track/pulse-basin conflict applies. These diagnostic counts include unsuccessful contexts, not independent missing physical pulses.')
    for r in gaps:lines.append(f"- {r['method']} [{r['start_s']:.9f}, {r['end_s']:.9f}] s, {r['duration_s']*1000:.3f} ms, {r['kind']}.")
    lines+=['','## A/B on identical time windows','','Each row uses only complete event intervals inside the closed comparison window, retaining original segment boundaries. Different-phase endpoint cycles and counts can differ; exact periods are in same_window_periods.csv.']
    for r in same:
        if r['window_name'] in ('A_longest','B_longest'):lines.append(f"- {r['window_name']} [{r['window_start_s']:.9f}, {r['window_end_s']:.9f}]: {r['method']} {r['n_cycles']} periods, {r['n_adjacent_period_differences']} adjacent differences; Jitter={f(r['jitter_local_percent'])}%, RAP={f(r['rap_percent'])}%, PPQ5={f(r['ppq5_percent'])}%.")
    lines+=['','## Diagnostic sensitivity','','One parameter is changed at a time. Status columns still use the original structural QC. New connections, scores below original cutoffs, event changes, and possible skipped periods are saved for inspection. More coverage alone is not an improvement criterion.']
    for r in sensitivity:
        if r['method']=='A':lines.append(f"- {r['variant_id']}: tracked={r['n_tracked_cycles']}, longest={r['n_main_cycles']}, covered={r['covered_time_fraction']:.4%}, Jitter={f(r['jitter_local_percent'])}%, pooled={f(r['pooled_jitter_percent'])}%, new connections={r['n_new_connections']}, below strict rho={r['n_new_connections_below_strict_rho']}, below strict margin={r['n_new_connections_below_strict_margin']}.")
    lines+=['','## Interpretation limits','','The strict main value samples a selected interval and may differ greatly from other segments. Pooled-within-segment differences avoid bridging but still omit untracked intervals and do not recover whole-window truth. Human-readable break figures and new-connection checks are required before recommending threshold changes. No permanent threshold change is made by this diagnostic option.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')

def optimal_event_assignment(reference,automatic,tolerance):
    """Maximum-cardinality, minimum-total-absolute-error one-to-one assignment."""
    from scipy.optimize import linear_sum_assignment
    a=np.asarray(reference,float);b=np.asarray(automatic,float);n=len(a);k=len(b)
    if not n or not k:return []
    penalty=n+k+1.;cost=np.full((n+k,n+k),penalty*10*(n+k+1));delta=np.abs(a[:,None]-b[None,:])
    cost[:n,:k]=np.where(delta<=tolerance+1e-12,delta/tolerance,cost[:n,:k])
    for i in range(n):cost[i,k+i]=penalty
    for j in range(k):cost[n+j,j]=penalty
    cost[n:,k:]=0.;rr,cc=linear_sum_assignment(cost)
    return sorted((int(i),int(j),float(b[j]-a[i])) for i,j in zip(rr,cc) if i<n and j<k and delta[i,j]<=tolerance+1e-12)


def read_manual_reference(path,case_id,data,start,end):
    with path.open(newline='') as f:
        reader=csv.DictReader(f);fields=reader.fieldnames or [];rows=list(reader)
    legacy='time_s' in fields and 'event_time_s' not in fields
    required={'case_id','phase','event_time_s','confidence','reviewer','notes'}
    if not legacy and not required.issubset(fields):raise ValueError('manual CSV requires '+','.join(sorted(required)))
    result=[];previous={};polarity_peaks={}
    for phase in ('A','B'):
        ix,p=signal.find_peaks(data['x'] if phase=='A' else -data['x'],plateau_size=True)
        polarity_peaks[phase]=(data['t'][p['left_edges']],data['t'][p['right_edges']])
    for source in rows:
        if not legacy and source['case_id']!=case_id:continue
        phase='A' if legacy else source['phase'];confidence='high' if legacy else source['confidence']
        if phase not in ('A','B') or confidence not in ('high','uncertain'):raise ValueError('manual phase must be A/B and confidence high/uncertain')
        at=float(source['time_s' if legacy else 'event_time_s'])
        if legacy and np.isfinite(at) and not start<=at<end:continue
        if not np.isfinite(at) or not start<=at<end:raise ValueError('manual event time must be finite seconds inside [start,end)')
        if phase in previous and at<=previous[phase]:raise ValueError('manual times must be strictly increasing per phase, without duplicates')
        previous[phase]=at;left,right=polarity_peaks[phase]
        distance=np.maximum(np.maximum(left-at,at-right),0.)
        error=float(np.min(distance)) if len(distance) else None
        if error is None or error>150e-6+1e-12:raise ValueError('manual event does not agree with the requested flow maximum/minimum phase within 150 us')
        result.append({'case_id':case_id,'phase':phase,'event_time_s':at,'confidence':confidence,'reviewer':source.get('reviewer','legacy-manual'),
                       'notes':source.get('notes',''),'phase_extremum_error_s':error,'input_format':'legacy_time_s' if legacy else 'phase_reference_v1'})
    if rows and not result:raise ValueError('manual CSV contains no events for --case-id '+case_id)
    return result


def evaluate_manual_reference(rows,tracks,start,end):
    statistics=[];matches=[];period_checks=[];metric_periods=[]
    for phase in ('A','B'):
        all_ref=[r for r in rows if r['phase']==phase]
        if not all_ref:continue
        left=min(r['event_time_s'] for r in all_ref);right=max(r['event_time_s'] for r in all_ref)
        for scope in ('high_only','all_annotations'):
            refs=[r for r in all_ref if scope=='all_annotations' or r['confidence']=='high']
            if not refs:continue
            for tol in (50e-6,100e-6,150e-6):
                uncertain=[r['event_time_s'] for r in all_ref if r['confidence']=='uncertain'] if scope=='high_only' else []
                auto=[e for e in tracks[phase]['events'] if e['accepted'] and left-tol<=e['subsample_time_s']<=right+tol and
                      not any(abs(e['subsample_time_s']-u)<=tol for u in uncertain)]
                a=[r['event_time_s'] for r in refs];b=[r['subsample_time_s'] for r in auto];pairs=optimal_event_assignment(a,b,tol)
                by_ref={i:(j,error) for i,j,error in pairs};used_auto={j for _,j,_ in pairs};errors=[abs(v[2]) for v in pairs]
                for i,row in enumerate(refs):
                    j,error=by_ref.get(i,(None,None));matches.append({'phase':phase,'scope':scope,'tolerance_us':tol*1e6,
                        'classification':'TP' if j is not None else 'FN','manual_time_s':row['event_time_s'],
                        'automatic_event_id':auto[j]['event_id'] if j is not None else None,'automatic_time_s':b[j] if j is not None else None,
                        'signed_error_s':error,'confidence':row['confidence']})
                for j,event in enumerate(auto):
                    if j not in used_auto:matches.append({'phase':phase,'scope':scope,'tolerance_us':tol*1e6,'classification':'FP',
                        'automatic_event_id':event['event_id'],'automatic_time_s':b[j],'manual_time_s':None,
                        'note':'unmatched event; phase shift is a hypothesis, not an automatic diagnosis'})
                # Uncertain annotations split high-confidence manual period sequences.
                blocks=[];time_blocks=[];block=[]
                for row in all_ref:
                    if scope=='high_only' and row['confidence']!='high':
                        if block:blocks.append(np.diff(block));time_blocks.append(block);block=[]
                    else:block.append(row['event_time_s'])
                if block:blocks.append(np.diff(block));time_blocks.append(block)
                manual_metrics=metric_blocks_v2(blocks,right-left if right>left else 1.)
                auto_blocks=[];auto_detail_chunks=[]
                for seg in tracks[phase]['segments']:
                    chunk=[];ev=[tracks[phase]['events'][j] for j in seg['event_indices']]
                    for x0,x1 in zip(ev,ev[1:]):
                        t0,t1=x0['subsample_time_s'],x1['subsample_time_s']
                        if left<=t0 and t1<=right and not any(t0-tol<=u<=t1+tol for u in uncertain):
                            chunk.append({'start_s':t0,'end_s':t1,'period_s':t1-t0,'event_start':x0['event_id'],'event_end':x1['event_id'],'source_segment_id':seg['segment_id']})
                        elif chunk:auto_detail_chunks.append(chunk);chunk=[]
                    if chunk:auto_detail_chunks.append(chunk)
                auto_blocks=[[r['period_s'] for r in chunk] for chunk in auto_detail_chunks]
                for kind,chunks in [('automatic',auto_detail_chunks),('manual',[[{'start_s':t0,'end_s':t1,'period_s':t1-t0} for t0,t1 in zip(ts,ts[1:])] for ts in time_blocks])]:
                    for number,chunk in enumerate(chunks):
                        for detail in chunk:metric_periods.append(dict(detail,phase=phase,scope=scope,tolerance_us=tol*1e6,reference=kind,contiguous_block_id=number,
                            used_for_jitter_v2=len(chunk)>=2,used_for_RAP_v2=len(chunk)>=3,used_for_PPQ5_v2=len(chunk)>=5))
                auto_metrics=metric_blocks_v2(auto_blocks,right-left if right>left else 1.)
                tp=len(pairs);fp=len(auto)-tp;fn=len(refs)-tp;precision=tp/len(auto) if auto else None;recall=tp/len(refs) if refs else None
                missing_periods=doubled=0
                for seg in tracks[phase]['segments']:
                    ev=[tracks[phase]['events'][j] for j in seg['event_indices']]
                    for x0,x1 in zip(ev,ev[1:]):
                        t0,t1=x0['subsample_time_s'],x1['subsample_time_s']
                        if t0<left or t1>right:continue
                        interior=[at for at in a if t0+tol<at<t1-tol]
                        missing_periods+=bool(interior)
                        period_checks.append({'phase':phase,'scope':scope,'tolerance_us':tol*1e6,'event_start':x0['event_id'],
                            'event_end':x1['event_id'],'start_s':t0,'end_s':t1,'n_interior_manual_events':len(interior),
                            'status':'POSSIBLE_SKIPPED_MANUAL_CYCLE' if interior else 'NO_INTERIOR_MANUAL_EVENT'})
                for t0,t1 in zip(a,a[1:]):doubled+=sum(t0+tol<at<t1-tol for at in b)>0
                statistics.append({'phase':phase,'scope':scope,'tolerance_us':tol*1e6,'status':'MANUAL_REFERENCE_EVALUATED',
                    'annotation_extent_start_s':left,'annotation_extent_end_s':right,'requested_window_fraction':(right-left)/(end-start),
                    'TP':tp,'FP':fp,'FN':fn,'precision':precision,'recall':recall,
                    'F1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
                    'median_absolute_time_error_s':float(np.median(errors)) if errors else None,
                    'p95_absolute_time_error_s':float(np.percentile(errors,95)) if errors else None,
                    'n_automatic_periods_spanning_manual_event':missing_periods,'n_manual_intervals_with_extra_auto_event':doubled,
                    'high_reference_events':sum(r['confidence']=='high' for r in all_ref),
                    'reference_quantity_warning':'FEWER_THAN_31_HIGH_EVENTS' if sum(r['confidence']=='high' for r in all_ref)<31 else None,
                    **{'manual_'+k:v for k,v in manual_metrics.items()},**{'auto_'+k:v for k,v in auto_metrics.items()},
                    'jitter_difference_pp':auto_metrics['jitter_local_percent']-manual_metrics['jitter_local_percent'] if auto_metrics['jitter_local_percent'] is not None and manual_metrics['jitter_local_percent'] is not None else None,
                    'scope_note':'same annotated time extent and phase; complete periods only, original automatic breaks retained; high-only excludes uncertain neighborhoods'})
    return statistics,matches,period_checks,metric_periods


def metric_blocks_v2(blocks,window_duration):
    out={}
    for width,name in [(2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')]:
        eligible=[np.asarray(b) for b in blocks if len(b)>=width];periods=[v for b in eligible for v in b];terms=[]
        for v in eligible:terms.extend(np.abs(np.diff(v)) if width==2 else np.abs(v[width//2:len(v)-width//2]-np.convolve(v,np.ones(width)/width,'valid')))
        mean=float(np.mean(periods)) if periods else None
        out[name]=float(100*np.mean(terms)/mean) if terms and mean else None
        out[name+'_mean_period_s']=mean;out[name+'_n_cycles']=len(periods);out[name+'_n_terms']=len(terms)
        out[name+'_covered_time_fraction']=float(sum(periods)/window_duration)
    return out


def export_manual_review(out,data,start,end,result,case_id,review_start=None,review_end=None):
    left=start if review_start is None else review_start;right=end if review_end is None else review_end
    if not start<=left<right<=end:raise ValueError('manual review extent must be inside analysis window')
    directory=out/'manual_review';directory.mkdir();tracks={k:result['_tracks'][result['_primary_key'] if k=='A' else result['_control_key']] for k in ('A','B')}
    period=result['pitch']['period_s'] or .008;panels=[];cursor=left;k=0
    while cursor<right-1e-12:
        k+=1;stop=min(cursor+.060,right);panels.append((f'panel_{k:03d}',cursor,stop,'complete-window panel'))
        if stop>=right:break
        cursor=stop-min(period,.015)
    important=[(.5497,.5644),(.6815,.6965),(.7491,.7664),(.84,.89),(.93,.98)] if case_id=='r0p30' else [(.604,.624),(.785,.812),(.5857,.600428),(.645,.69),(.945,.975)] if case_id=='r0p25' else []
    for i,(a,b) in enumerate(important,1):
        a=max(left,a-period);b=min(right,b+period)
        if a<b:panels.append((f'focus_{i:03d}',a,b,'specified failure/regression or control'))
    index=[];t=data['t'];x=data['x']
    for panel_id,a,b,note in panels:
        breaks=[dict(row,method=method) for method,tr in tracks.items() for row in tr.get('_diagnostic',{}).get('breaks',[]) if a<=row['break_time_s']<=b]
        for blind in (True,False):
            fig,ax=plt.subplots(figsize=(14,5.5));mask=(t>=a)&(t<=b);ax.plot(t[mask],x[mask]*1e6,color='black',lw=.9,label='raw Ug (20 kHz; no smoothing)')
            if not blind:
                for method,color,marker in [('A','tab:blue','o'),('B','tab:green','v')]:
                    tr=tracks[method]
                    for adopted in (False,True):
                        ev=[e for e in tr['events'] if e['accepted']==adopted and a<=e['subsample_time_s']<=b];at=[e['subsample_time_s'] for e in ev]
                        ax.scatter(at,np.interp(at,t,x)*1e6,color=color if adopted else 'gray',s=26 if adopted else 17,marker=marker if adopted else 'x',label=f'{method} '+('adopted' if adopted else 'candidate'))
                    for segment in tr['segments']:
                        at=[tr['events'][j]['subsample_time_s'] for j in segment['event_indices']]
                        if len(at)>1 and max(at)>a and min(at)<b:ax.axvspan(max(a,min(at)),min(b,max(at)),color=color,alpha=.045)
                        for at0 in (min(at),max(at)):
                            if a<=at0<=b:ax.axvline(at0,color=color,alpha=.4,lw=.7,ls=':')
                text=[]
                for row in breaks:
                    if row.get('resolved_by_recovery'):continue
                    at=row['predicted_next_time_s']
                    if a<=at<=b:ax.axvline(at,color='red',ls='--',lw=.7,alpha=.4)
                    rho='none' if row['best_rho'] is None else f'{row["best_rho"]:.3f}';score='none' if row['best_score'] is None else f'{row["best_score"]:.3f}'
                    text.append(f'{row["break_id"]}: rho={rho}, score={score}, {"/".join(row["direct_reasons"])}')
                ax.text(0,-.19,'\n'.join(text[:7])+ ('\n(more contexts in tracking_breaks_detailed.csv)' if len(text)>7 else ''),transform=ax.transAxes,fontsize=7,va='top')
                ax.plot([],[],color='red',ls='--',label=f'predicted next; ACF T0={period*1e3:.3f} ms (not a detected event)')
            ax.set(title=f'{case_id} {panel_id} [{a:.6f}, {b:.6f}] s — '+('BLIND: mark consistent A maxima / B minima' if blind else 'DIAGNOSTIC: events, segment spans/boundaries, failed contexts'),xlabel='Time [s]',ylabel='Ug [mL/s]',xlim=(a,b));ax.grid(alpha=.15)
            if not blind:ax.legend(fontsize=7,ncol=3,loc='upper left')
            fig.tight_layout(rect=(0,.15 if not blind else 0,1,1));fig.savefig(directory/f'{panel_id}_{"blind" if blind else "diagnostic"}.png',dpi=150);plt.close(fig)
        index.append({'case_id':case_id,'panel_id':panel_id,'t_start_s':a,'t_end_s':b,'figure_blind':f'manual_review/{panel_id}_blind.png',
                      'figure_diagnostic':f'manual_review/{panel_id}_diagnostic.png','break_ids':[r['break_id'] for r in breaks],'notes':note})
    write_csv(out/'manual_review_index.csv',index)
    write_csv(out/'manual_events_template.csv',[],['case_id','phase','event_time_s','confidence','reviewer','notes'])
    (out/'manual_annotation_instructions.md').write_text(f'# Manual reference: {case_id}\n\nStart with the blind panels in manual_review_index.csv. Record actual corresponding maxima (A) or minima (B) in seconds, strictly increasing within each phase. Use confidence high or uncertain and preserve a consistent phase. Do not copy diagnostic event markers as truth. The template contains a header only: no fabricated events. Aim for at least 31 high-confidence events (30 periods) per case, preferably the complete window, covering difficult and control intervals. Overlapping panels must not duplicate annotations. Uncertain events split high-only metric sequences. A long flat minimum may not identify a unique time: mark uncertain.\n\nColumns: case_id,phase,event_time_s,confidence,reviewer,notes. Use --case-id {case_id} --manual-events <filled.csv> and a NEW --out-dir with the same input/dt/window/policy. Precision/recall use same-phase one-to-one maximum-cardinality/minimum-error assignments at 50/100/150 us; 100 us is primary. Unannotated time outside first/last manual event is excluded; annotate all visible events inside that range or label uncertain. CSV validation does not certify physical event identity.\n')
    return index


def reanchor_validation_outputs(out,data,start,end,result,tracking,min_cycles,case_id,baseline_path=None):
    """Freeze legacy results, expose both aggregation definitions and every recovery edge."""
    legacy_cfg={k:v for k,v in tracking.items() if k not in ('reanchor_policy','reanchor_attempt_cap')}
    legacy=analyze_arrays(data['t'],data['x'],data['fs'],start,end,legacy_cfg,CPP_CONFIG if result.get('_cpp_config') is None else result['_cpp_config'],min_cycles)
    if clean_json(legacy['cpp'])!=clean_json(result['cpp']):raise RuntimeError('CPP changed between legacy and recovery; refusing to publish')
    baseline=out/'baseline';baseline.mkdir();dump_json(baseline/'strict_results.json',legacy)
    baseline_check={'CPP_exactly_identical_current_runtime':True,'baseline_path':str(baseline_path) if baseline_path else None}
    if baseline_path:
        old=json.loads(baseline_path.read_text());meta=old['metadata']
        if meta['input_sha256']!=data['sha256'] or (meta['dt_s'],meta['start_s'],meta['end_s'])!=(data['dt_s'],start,end):raise ValueError('baseline waveform/dt/window mismatch')
        for method,key in [('A','primary'),('B','control')]:
            if old[key]['main']['event_ids']!=legacy[key]['main']['event_ids']:raise RuntimeError('legacy baseline main event IDs changed')
            for metric in ('jitter_local_percent','rap_percent','ppq5_percent'):
                a,b=old[key]['main'][metric],legacy[key]['main'][metric]
                if a is not None and not np.isclose(a,b,rtol=1e-10,atol=1e-10):raise RuntimeError('legacy baseline metric changed')
        a,b=old['cpp']['flow_CPP_dB'],legacy['cpp']['flow_CPP_dB']
        if (a is None)!=(b is None) or a is not None and not np.isclose(a,b,rtol=1e-9,atol=1e-9):raise RuntimeError('CPP differs from stored baseline')
        baseline_check['stored_baseline_values_and_main_IDs_reproduced']=True;baseline_check['baseline_sha256']=sha256(baseline_path)
    comparisons=[];jitter=[];states=[];connections=[];phase_rows=[];breaks=[];terms=[];halves=[];same=[];same_periods=[];frozen_events=[];frozen_periods=[];frozen_segments=[]
    tracks={};summaries={};baseline_tracks={};baseline_summaries={}
    for method,key,private in [('A','primary','_primary_key'),('B','control','_control_key')]:
        tr=result['_tracks'][result[private]];sm=result[key];oldtr=legacy['_tracks'][legacy[private]];oldsm=legacy[key]
        tracks[method]=tr;summaries[method]=sm;baseline_tracks[method]=oldtr;baseline_summaries[method]=oldsm
        oldids={e['event_id'] for e in oldtr['events'] if e['accepted'] and e['in_window']};newids={e['event_id'] for e in tr['events'] if e['accepted'] and e['in_window']}
        oldpairs={(r['event_start'],r['event_end']) for r in oldsm['_periods']};newpairs={(r['event_start'],r['event_end']) for r in sm['_periods']}
        comparisons.append({'case_id':case_id,'method':method,'baseline_policy':'legacy','current_policy':tracking.get('reanchor_policy','legacy'),
            'old_adopted_ids':sorted(oldids),'new_adopted_ids':sorted(newids),'new_event_ids':sorted(newids-oldids),'removed_event_ids':sorted(oldids-newids),
            'old_candidate_count':sum(e['in_window'] for e in oldtr['events']),'new_candidate_count':sum(e['in_window'] for e in tr['events']),
            'old_tracked_cycles':oldsm['n_complete_periods_all_segments'],'new_tracked_cycles':sm['n_complete_periods_all_segments'],
            'old_covered_time_fraction':oldsm['covered_time_fraction'],'new_covered_time_fraction':sm['covered_time_fraction'],
            'n_new_connections':len(newpairs-oldpairs),'n_removed_connections':len(oldpairs-newpairs),
            'new_short_periods':sm['n_short_periods'],'new_long_periods':sm['n_long_periods'],'duplicate_adopted_events':sm['n_duplicate_adopted_samples'],
            'search_limit':tr.get('search_limit',False),'status':'PENDING_MANUAL_VALIDATION'})
        for policy,tt,ss in [('legacy',oldtr,oldsm),('current',tr,sm)]:
            for version in ('longest_legacy','pooled_v1_baseline' if policy=='legacy' else 'pooled_v1_current','pooled_v2_baseline' if policy=='legacy' else 'pooled_v2_current'):
                row={'case_id':case_id,'method':method,'policy':policy,'version':version,'all_tracked_coverage':ss['covered_time_fraction'],'flow_CPP_dB':result['cpp']['flow_CPP_dB']}
                if version.startswith('longest'):
                    row.update({k:ss['main'][k] for k in ('jitter_local_percent','rap_percent','ppq5_percent','mean_period_s','n_cycles','period_cv_percent')})
                    row.update(covered_time_fraction=ss['main']['duration_s']/(end-start),start_s=ss['main'].get('start_s'),end_s=ss['main'].get('end_s'))
                    for width,name in [(2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')]:
                        row[name+'_n_terms']=max(0,ss['main']['n_cycles']-width+1);row[name+'_n_cycles']=ss['main']['n_cycles'] if ss['main']['n_cycles']>=width else 0
                        row[name+'_mean_period_s']=ss['main']['mean_period_s'] if ss['main']['n_cycles']>=width else None;row[name+'_covered_time_fraction']=row['covered_time_fraction'] if ss['main']['n_cycles']>=width else 0.
                elif version.startswith('pooled_v1'):
                    blocks=[np.diff([tt['events'][j]['subsample_time_s'] for j in ss['_membership'][s['segment_id']]]) for s in tt['segments']]
                    row.update(diagnostic_block_metrics(blocks));row['covered_time_fraction']=ss['covered_time_fraction']
                else:
                    v2=ss['pooled_within_segments_v2'];row.update({k:v2[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')})
                    for name,detail in v2['metric_details'].items():
                        for k in ('n_terms','n_cycles','mean_period_s','covered_time_fraction'):row[name+'_'+k]=detail[k]
                    row['period_cv_percent']=v2['metric_details']['jitter_local_percent']['period_cv_percent']
                    row['covered_time_fraction']=v2['metric_details']['jitter_local_percent']['covered_time_fraction']
                jitter.append(row)
            frozen_events.extend(dict(e,policy=policy) for e in tt['events']);frozen_periods.extend(dict(r,policy=policy) for r in ss['_periods']);frozen_segments.extend(dict(r,method=method,policy=policy) for r in ss['segments'])
        states.extend(tr.get('candidate_states',[]));terms.extend(sm['pooled_within_segments_v2']['_terms'])
        features=phase_features(data['t'],data['x'],data['fs'],tr,result['pitch']['period_s'],tracking) if result['pitch']['period_s'] else {}
        phase_rows.extend(dict(f,method=method) for f in features.values());by_id={e['event_id']:i for i,e in enumerate(tr['events'])}
        recovered={(r['event_start'],r['event_end']):r for r in tr.get('recovered_connections',[])}
        for p in sm['_periods']:
            pair=(p['event_start'],p['event_end'])
            if pair in oldpairs:continue
            info=recovered.get(pair,{});a,b=(by_id[k] for k in pair);fa,fb=features[a],features[b]
            other=legacy['_tracks'][legacy['_control_key'] if method=='A' else legacy['_primary_key']]
            row=dict(p,**{k:v for k,v in info.items() if k not in p});row.update(variant_id='recover-rejected',diagnostic_only=False,
                interval_over_acf_period=p['period_s']/result['pitch']['period_s'],n_strict_opposite_phase_events=sum(e['accepted'] and p['start_s']<=e['subsample_time_s']<p['end_s'] for e in other['events']),
                phase_before=fa['relative_phase_in_basin'],phase_after=fb['relative_phase_in_basin'],
                basin_start_before=fa['opposite_before_s'],basin_end_before=fa['opposite_after_s'],basin_start_after=fb['opposite_before_s'],basin_end_after=fb['opposite_after_s'],
                phase_gates='non-overlap; separate opposite-extremum basins; dominant candidate; no intervening distinct pulse basin; ACF interval bounds',
                unimplemented_gates='no trained phase classifier; no guarantee against subtle within-lobe shifts',review_status='AUTO_GATES_PASS_PENDING_MANUAL')
            if p['period_s']/result['pitch']['period_s']<tracking['search_min_period_ratio']-1e-2 or p['period_s']/result['pitch']['period_s']>tracking['search_max_period_ratio']+1e-2:raise RuntimeError('new connection violates ACF safety bounds')
            connections.append(row)
        detailed=diagnostic_break_details(tr,sm,start,end)
        y=data['x'] if method=='A' else -data['x'];ix,prop=signal.find_peaks(y,prominence=0);threshold=tr.get('amplitude_p95_p5',0.)*tr['prominence_percent']/100
        for row in detailed:
            indexes=[k for k,v in enumerate(ix) if row['search_start_s']<=data['t'][v]<=row['search_end_s']];qualified=[k for k in indexes if prop['prominences'][k]>=threshold]
            ev=next(e for e in tr['events'] if e['event_id']==row['previous_event']);row['previous_reference_adopted']=ev['accepted']
            row['n_all_raw_extrema_in_search']=len(indexes);row['n_qualified_raw_extrema_in_search']=len(qualified)
            reasons=row['direct_reasons']
            if 'LOW_SHAPE_CORR' in reasons or 'CORRELATION_BELOW_MIN' in reasons:category='LOW_SHAPE_CORR'
            elif row['candidate_count']==0:
                category='NO_RAW_PEAK' if not indexes else 'BELOW_PROMINENCE' if not qualified else 'REANCHOR_EXCLUSION_CONFLICT' if row.get('n_reanchor_blocked_candidates_in_search',0) else 'INSUFFICIENT_TEMPLATE_EDGE' if row['tracking_reason']=='DATA_EDGE' else 'TRACK_CONFLICT'
            else:category=' / '.join(reasons)
            row['reason_category']=category;breaks.append(row)
        for label,left,right in [('early_half',start,(start+end)/2),('late_half',(start+end)/2,end)]:
            for policy,tt in [('legacy',oldtr),('current',tr)]:
                blocks,details=diagnostic_clip_blocks(tt,left,right)
                halves.append({'case_id':case_id,'method':method,'policy':policy,'window_name':label,'start_s':left,'end_s':right,
                               'tracked_covered_time_fraction':sum(r['period_s'] for r in details)/(right-left),**metric_blocks_v2(blocks,right-left)})
    for policy,tt,ss in [('legacy',baseline_tracks,baseline_summaries),('current',tracks,summaries)]:
        rows,ps=diagnostic_same_windows(tt,ss,start,end)
        for row in rows:
            blocks,_=diagnostic_clip_blocks(tt[row['method']],row['window_start_s'],row['window_end_s'])
            v2=metric_blocks_v2(blocks,row['window_end_s']-row['window_start_s'])
            same.append(dict(row,policy=policy,aggregation_version='pooled_within_segments_v2',
                **{'pooled_v1_'+k:row[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent','mean_period_s','n_cycles')},**v2,
                mean_period_s=v2['jitter_local_percent_mean_period_s'],n_cycles=v2['jitter_local_percent_n_cycles'],
                covered_time_fraction=v2['jitter_local_percent_covered_time_fraction']))
        same_periods.extend(dict(p,policy=policy) for p in ps)
        for left,right,name in [(start,end,'full_requested'),(start,(start+end)/2,'early_half'),((start+end)/2,end,'late_half')]:
            for method in ('A','B'):
                blocks,details=diagnostic_clip_blocks(tt[method],left,right);same.append({'policy':policy,'window_name':name,'method':method,'window_start_s':left,'window_end_s':right,'aggregation_version':'pooled_within_segments_v2',**metric_blocks_v2(blocks,right-left)})
                same_periods.extend(dict(p,policy=policy,window_name=name) for p in details)
    period_counts={}
    for row in same_periods:
        key=(row['policy'],row['window_name'],row['method'],row['segment_id']);period_counts[key]=period_counts.get(key,0)+1
    for row in same_periods:
        count=period_counts[(row['policy'],row['window_name'],row['method'],row['segment_id'])]
        row.update(used_for_jitter_v2=count>=2,used_for_RAP_v2=count>=3,used_for_PPQ5_v2=count>=5)
    write_csv(baseline/'events_and_comparison.csv',frozen_events);write_csv(baseline/'periods_and_comparison.csv',frozen_periods);write_csv(baseline/'segments_and_comparison.csv',frozen_segments)
    write_csv(baseline/'events.csv',[r for r in frozen_events if r['policy']=='legacy']);write_csv(baseline/'periods.csv',[r for r in frozen_periods if r['policy']=='legacy']);write_csv(baseline/'segments.csv',[r for r in frozen_segments if r['policy']=='legacy'])
    write_csv(out/'reanchor_comparison.csv',comparisons);write_csv(out/'jitter_comparison.csv',jitter);write_csv(out/'reanchor_candidate_states.csv',states)
    write_csv(out/'tracking_breaks_detailed.csv',breaks);write_csv(out/'recovery_context_attempts.csv',[dict(r,method=method) for method,tr in tracks.items() for r in tr.get('recovery_attempts',[])])
    audits=[]
    for row in connections:
        tr=tracks[row['method']];lookup={e['event_id']:e for e in tr['events']};a,b=lookup[row['event_start']],lookup[row['event_end']]
        opposite=tracks['B' if row['method']=='A' else 'A']
        row['n_current_opposite_phase_events']=sum(e['accepted'] and row['start_s']<=e['subsample_time_s']<row['end_s'] for e in opposite['events'])
        row['absolute_basin_phase_change']=abs(row['phase_after']-row['phase_before']) if row['phase_before'] is not None and row['phase_after'] is not None else None
        raw_ratio=(b['raw_time_s']-a['raw_time_s'])/result['pitch']['period_s']
        checks={'correlation_pass':row['rho']>=tracking['min_correlation'],
                'margin_pass':row['score_margin'] is None or row['score_margin']>=tracking['min_score_margin'],
                'candidate_rank_one':row['candidate_rank']==1,'no_safety_rejection':not row['safety_reasons'],
                'raw_ACF_interval_pass':tracking['search_min_period_ratio']-1e-12<=raw_ratio<=tracking['search_max_period_ratio']+1e-12,
                'real_unique_ordered_events':a['sample_id']<b['sample_id'] and a['accepted'] and b['accepted']}
        if not all(checks.values()):raise RuntimeError('new connection failed independent automatic audit: '+row['event_start']+' '+row['event_end'])
        audits.append(dict(row,**checks,raw_interval_over_ACF_period=raw_ratio,
            manual_status='PENDING_MANUAL_ANNOTATION',audit_interpretation='automatic checks only; different-phase agreement and basin phase do not prove correct event identity'))
    write_csv(out/'phase_features.csv',phase_rows);write_csv(out/'new_connections.csv',connections);write_csv(out/'new_connection_audit.csv',audits);write_csv(out/'metric_terms_v2.csv',terms)
    write_csv(out/'half_window_metrics.csv',halves);write_csv(out/'same_absolute_window_metrics.csv',same);write_csv(out/'same_absolute_window_periods.csv',same_periods)
    regressions=[]
    points=[('local_peak_7975','A',.7974,.7976),('B_2T_connection','B',.5857,.6004278977)] if case_id=='r0p25' else [('local_peak_55845','A',.55835,.55855)] if case_id=='r0p30' else []
    for label,method,left,right in points:
        tr=tracks[method];sm=summaries[method]
        if label=='B_2T_connection':
            matches=[r for r in sm['_periods'] if abs(r['start_s']-left)<100e-6 and abs(r['end_s']-right)<100e-6]
            regressions.append({'case_id':case_id,'regression':label,'method':method,'start_s':left,'end_s':right,'accepted':bool(matches),'result':'RISK_CONNECTION_PRESENT' if matches else 'RISK_CONNECTION_NOT_ADOPTED'})
        else:
            ev=[e for e in tr['events'] if left<=e['raw_time_s']<=right]
            regressions.append({'case_id':case_id,'regression':label,'method':method,'start_s':left,'end_s':right,'accepted':any(e['accepted'] for e in ev),
                'event_ids':[e['event_id'] for e in ev],'candidate_reasons':[e['reason'] for e in ev],'result':'RISK_CANDIDATE_ADOPTED' if any(e['accepted'] for e in ev) else 'RISK_CANDIDATE_NOT_ADOPTED' if ev else 'NOT_A_CURRENT_PROMINENCE_CANDIDATE'})
    write_csv(out/'mandatory_regression_cases.csv',regressions)
    diagnostic_new_connection_figures(out/'recovery_connection_figures',data,result['pitch']['period_s'],baseline_tracks,{(m,'recover-rejected'):tr for m,tr in tracks.items()},connections) if result['pitch']['period_s'] else None
    for method in ('A','B'):pooled_coverage_figure(out,data,start,end,tracks[method],summaries[method])
    dump_json(out/'baseline_check.json',baseline_check)
    return {'baseline_check':baseline_check,'comparisons':comparisons,'jitter_comparison':jitter,'manual_validation_status':'PENDING_MANUAL_ANNOTATION',
            'missing_reason_categories':[{'method':row['method'],'stage':row.get('stage'),'resolved_by_recovery':row.get('resolved_by_recovery',False),'reason':row['reason_category']} for row in breaks if row['in_window']],
            'new_connection_status':'AUTO_GATES_PASS_PENDING_MANUAL','shape_gate_limits':'basin diagnostics do not establish human phase identity',
            'research_primary_candidate':'pooled_within_segments_v2'}


def pooled_coverage_figure(out,data,start,end,track,summary):
    from matplotlib.patches import Patch
    fig,axes=plt.subplots(3,1,figsize=(15,9),sharex=True);mask=(data['t']>=start)&(data['t']<end)
    for ax,(width,name) in zip(axes,[(2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')]):
        ax.plot(data['t'][mask],data['x'][mask]*1e6,color='black',lw=.6)
        for gap in diagnostic_gaps(summary,start,end,track['method']):ax.axvspan(gap['start_s'],gap['end_s'],color='tab:red',alpha=.12)
        for seg in summary['segments']:
            if not seg['n_cycles']:continue
            ax.axvspan(seg['start_s'],seg['end_s'],color='tab:green' if seg['n_cycles']>=width else 'gray',alpha=.12)
        detail=summary['pooled_within_segments_v2']['metric_details'][name]
        value='undefined' if detail['value'] is None else f'{detail["value"]:.4f}%'
        ax.set(title=f'{track["method"]} {name}: pooled_v2={value}; {detail["n_cycles"]} periods, {detail["n_terms"]} terms; used time {detail["covered_time_fraction"]:.2%}',ylabel='Ug [mL/s]',xlim=(start,end));ax.grid(alpha=.1)
    axes[0].legend(handles=[Patch(color='tab:green',alpha=.2,label='used for this metric'),Patch(color='gray',alpha=.2,label='tracked segment too short'),Patch(color='tab:red',alpha=.2,label='no complete tracked period')],fontsize=8,ncol=3)
    axes[-1].set_xlabel('Time [s]');fig.suptitle('Metric-specific coverage; gaps are never bridged; manual reference status is recorded in summary.json');fig.tight_layout();fig.savefig(out/f'pooled_v2_coverage_{track["method"]}.png',dpi=170);plt.close(fig)


def global_local_periods(data,start,end,cfg):
    centers=np.unique(np.r_[np.arange(start,end,cfg['acf_hop_sec']),end]);frames=[];candidates=[]
    t=data['t'];x=data['x'];fs=data['fs']
    for number,center in enumerate(centers):
        left=max(start,center-cfg['acf_window_sec']/2);right=min(end,center+cfg['acf_window_sec']/2)
        mask=(t>=left)&(t<right);v=x[mask];p=pitch_candidates(v,fs,cfg)
        amplitude=float(np.percentile(v,95)-np.percentile(v,5)) if len(v) else 0.
        dif=np.diff(v);noise=float(1.4826*np.median(abs(dif-np.median(dif)))/np.sqrt(2)) if len(dif) else 0.
        confidence=p.get('chosen',{}).get('acf',max((r['acf'] for r in p['candidates']),default=0.))
        frames.append({'frame_id':number,'center_s':float(center),'start_s':left,'end_s':right,'n_samples':int(mask.sum()),
            'period_s':p['period_s'],'acf_confidence':confidence,'local_amplitude':amplitude,'difference_MAD_noise_sigma':noise,
            'status':p['status'],'flags':p['flags'],'dominant_spectral_frequency_hz':p.get('dominant_spectral_frequency_hz'),
            'alternatives':p['alternatives']})
        for rank,row in enumerate(sorted(p['candidates'],key=lambda r:(-r['acf'],r['period_s'])),1):
            candidates.append(dict(row,frame_id=number,center_s=float(center),window_start_s=left,window_end_s=right,
                rank=rank,selected=row.get('period_s')==p['period_s'],status=p['status'],flags=p['flags'],alternatives=p['alternatives']))
    return frames,candidates


def global_period_at(at,frames):
    centers=np.array([r['center_s'] for r in frames]);valid=[r for r in frames if r['period_s'] is not None]
    nearest=min(frames,key=lambda r:abs(r['center_s']-at));confidence=float(np.interp(at,centers,[r['acf_confidence'] for r in frames]))
    period=float(np.interp(at,[r['center_s'] for r in valid],[r['period_s'] for r in valid])) if valid else None
    supported=nearest['period_s'] is not None and confidence>=QC_CONFIG['acf_minimum']
    return period,confidence,supported,nearest['flags']


def global_raw_candidates(data,start,end,frames,cfg):
    t=data['t'];x=data['x'];fs=data['fs'];lo=max(0,np.searchsorted(t,start-.05));hi=min(len(t),np.searchsorted(t,end+.05));out={}
    for phase,polarity in [('A',1),('B',-1)]:
        y=polarity*x;ix,props=signal.find_peaks(y[lo:hi],prominence=0,plateau_size=True);rows=[]
        for k,(local,prom) in enumerate(zip(ix,props['prominences'])):
            i=int(local+lo);at=float(t[i]);period,confidence,supported,flags=global_period_at(at,frames)
            radius=1.5*(period or 1/cfg['pitch_min_hz']);a=max(0,np.searchsorted(t,at-radius));b=min(len(t),np.searchsorted(t,at+radius))
            v=x[a:b];amp=float(np.percentile(v,95)-np.percentile(v,5));dif=np.diff(v)
            noise=float(1.4826*np.median(abs(dif-np.median(dif)))/np.sqrt(2)) if len(dif) else 0.
            offset,_,status=quadratic(y,i);pleft=int(props['left_edges'][k]+lo);pright=int(props['right_edges'][k]+lo)
            rows.append({'event_id':f'{phase}_{i}','method':phase,'sample_id':i,'raw_time_s':at,'subsample_time_s':at+offset/fs,
                'offset_us':offset/fs*1e6,'interpolation_status':status,'prominence':float(prom),'local_amplitude':amp,
                'local_window_start_s':float(t[a]),'local_window_end_s':float(t[b-1]),'local_window_samples':b-a,
                'normalized_prominence':float(prom/amp) if amp>0 else None,'noise_sigma':noise,
                'plateau_start_s':float(t[pleft]),'plateau_end_s':float(t[pright]),'plateau_width_s':float((pright-pleft)/fs),
                'timing_uncertainty_s':float((pright-pleft+1)/(2*fs)),
                'local_period_s':period,'local_acf_confidence':confidence,'period_supported':supported,'local_period_flags':flags,
                'in_window':bool(start<=at+offset/fs<end)})
        out[phase]=rows
    return out


def global_candidate_features(data,raw,phase,cfg):
    rows=copy.deepcopy(raw[phase]);opposite=raw['B' if phase=='A' else 'A'];x=data['x'];t=data['t'];polarity=1 if phase=='A' else -1
    def preliminary(e):
        return e['normalized_prominence'] is not None and e['normalized_prominence']>=cfg['prominence_percent']/100 and e['prominence']>=cfg['noise_sigma_multiplier']*e['noise_sigma']
    opposites=[e for e in opposite if preliminary(e)];opp_times=np.array([e['raw_time_s'] for e in opposites]);cells={};windows={}
    for j,e in enumerate(rows):
        pos=int(np.searchsorted(opp_times,e['raw_time_s']));left=opposites[pos-1] if pos else None;right=opposites[pos] if pos<len(opposites) else None
        heights=[polarity*(x[e['sample_id']]-x[o['sample_id']]) for o in (left,right) if o is not None]
        delta=max(0.,min(heights)) if len(heights)==2 else max([e['prominence'],*heights])
        cell=(left['sample_id'] if left else None,right['sample_id'] if right else None);e['pulse_cell']=cell
        e['opposite_before_s']=left['raw_time_s'] if left else None;e['opposite_after_s']=right['raw_time_s'] if right else None
        e['peak_trough_amplitude']=float(delta);e['normalized_excursion']=float(delta/e['local_amplitude']) if e['local_amplitude']>0 else None
        e['relative_phase_in_basin']=(e['raw_time_s']-left['raw_time_s'])/(right['raw_time_s']-left['raw_time_s']) if left and right else None
        e['rise_duration_s']=e['raw_time_s']-left['raw_time_s'] if left else None;e['fall_duration_s']=right['raw_time_s']-e['raw_time_s'] if right else None
        reasons=[]
        if e['normalized_prominence'] is None:reasons.append('UNDEFINED_LOCAL_AMPLITUDE')
        elif e['normalized_prominence']<cfg['prominence_percent']/100:reasons.append('BELOW_LOCAL_PROMINENCE')
        if e['normalized_excursion'] is None or e['normalized_excursion']<cfg['excursion_percent']/100:reasons.append('BELOW_LOCAL_EXCURSION')
        if e['prominence']<cfg['noise_sigma_multiplier']*e['noise_sigma']:reasons.append('BELOW_MAD_NOISE_FLOOR')
        if not e['period_supported']:reasons.append('PERIOD_UNRESOLVED')
        e['selectable']=not reasons;e['rejection_reasons']=reasons;e['accepted']=False;e['anchor']=False;e['segment_id']=None;e['match_score']=None;e['correlation']=None
        e.update(threshold_prominence_percent=cfg['prominence_percent'],threshold_excursion_percent=cfg['excursion_percent'],
                 local_prominence_threshold=e['local_amplitude']*cfg['prominence_percent']/100,
                 local_excursion_threshold=e['local_amplitude']*cfg['excursion_percent']/100,
                 MAD_prominence_floor=cfg['noise_sigma_multiplier']*e['noise_sigma'])
        e['effective_prominence_threshold']=max(e['local_prominence_threshold'],e['MAD_prominence_floor'])
        e['reason']='ELIGIBLE_GLOBAL_CANDIDATE' if not reasons else ' / '.join(reasons)
        if preliminary(e):cells.setdefault(cell,[]).append(j)
        if e['local_period_s']:
            positions=e['raw_time_s']+np.linspace(-.4,.4,129)*e['local_period_s'];w=np.interp(positions,t,x);w-=w.mean();norm=np.linalg.norm(w)
            if norm>np.finfo(float).tiny:windows[j]=w/norm
            e['template_edge']=bool(positions[0]<t[0] or positions[-1]>t[-1])
        else:e['template_edge']=True
    for j,e in enumerate(rows):
        peers=cells.get(e['pulse_cell'],[]);rank=1+sum(rows[k]['prominence']>e['prominence'] for k in peers)
        secondary=bool(rank>1);e['basin_candidate_rank']=rank;e['n_candidates_in_basin']=len(peers);e['secondary_lobe']=secondary
        parts={'node_prominence':cfg['node_prominence_weight']*min(e['normalized_prominence'] or 0.,1.),
               'node_excursion':cfg['node_excursion_weight']*min(e['normalized_excursion'] or 0.,1.),
               'node_acf':cfg['node_acf_weight']*max(0.,e['local_acf_confidence']),
               'node_basin':cfg['node_basin_weight']*(not secondary),'secondary_penalty':-cfg['secondary_lobe_penalty']*secondary}
        e.update(parts);e['node_score']=float(sum(parts.values()));e['shape_flags']=(['SECONDARY_LOBE'] if secondary else [])+(['UNBOUNDED_BASIN'] if None in e['pulse_cell'] else [])+(['TEMPLATE_EDGE'] if e['template_edge'] else [])
    return rows,windows


def global_graph(rows,windows,frames,cfg):
    ids=[i for i,e in enumerate(rows) if e['selectable'] and e['in_window']];n=len(ids)
    if n>cfg['max_dp_candidates']:raise ValueError('SEARCH_LIMIT: too many eligible global candidates; raw data not discarded')
    transitions=np.full((n,n),-np.inf);connected=np.zeros((n,n),dtype=bool);audits=[];edge_lookup={}
    for b,j in enumerate(ids):
        ej=rows[j]
        for a in range(b):
            i=ids[a];ei=rows[i];delta=ej['raw_time_s']-ei['raw_time_s'];period,confidence,supported,flags=global_period_at((ei['raw_time_s']+ej['raw_time_s'])/2,frames)
            same=ei['pulse_cell']==ej['pulse_cell'] and None not in ei['pulse_cell']
            if delta>=cfg['search_min_period_ratio']*max(ej['local_period_s'] or 0.,ei['local_period_s'] or 0.) and not same:transitions[a,b]=-cfg['restart_penalty']
            if not period or delta>2.25*period or delta<.35*period:continue
            ratio=delta/period;rho=float(np.dot(windows[i],windows[j])) if i in windows and j in windows else None
            reasons=[]
            if not supported:reasons.append('PERIOD_UNRESOLVED')
            if not cfg['search_min_period_ratio']<=ratio<=cfg['search_max_period_ratio']:reasons.append('OUTSIDE_LOCAL_PERIOD_RANGE')
            if 1.8<=ratio<=2.2:reasons.append('DOUBLE_PERIOD_POSSIBLE')
            if same:reasons.append('SAME_PULSE_CELL')
            if rho is None or rho<cfg['incompatible_shape_correlation']:reasons.append('INCOMPATIBLE_SHAPE')
            intermediate=[k for k in ids[a+1:b] if not rows[k]['secondary_lobe'] and rows[k]['pulse_cell'] not in (ei['pulse_cell'],ej['pulse_cell'])
                and rows[k]['prominence']>=cfg['intermediate_pulse_strength_ratio']*min(ei['prominence'],ej['prominence'])
                and rows[k]['peak_trough_amplitude']>=cfg['intermediate_pulse_strength_ratio']*min(ei['peak_trough_amplitude'],ej['peak_trough_amplitude'])]
            if intermediate:reasons.append('SKIPPED_STRONG_INTERMEDIATE_PULSE')
            phase_difference=abs(ei['relative_phase_in_basin']-ej['relative_phase_in_basin']) if ei['relative_phase_in_basin'] is not None and ej['relative_phase_in_basin'] is not None else None
            shape=cfg['edge_shape_weight']*(rho if rho is not None else -1.)
            phase=cfg['edge_phase_weight']*(1-2*phase_difference) if phase_difference is not None else 0.
            acf=cfg['edge_acf_weight']*max(confidence,0.);score=shape+phase+acf
            row={'event_start':ei['event_id'],'event_end':ej['event_id'],'start_s':ei['subsample_time_s'],'end_s':ej['subsample_time_s'],
                 'raw_delta_s':delta,'local_period_s':period,'interval_over_local_period':ratio,'shape_correlation':rho,
                 'phase_difference':phase_difference,'edge_shape_score':shape,'edge_phase_score':phase,'edge_acf_score':acf,
                 'period_regularization_score':0.,'edge_score':score,'legal':not reasons,'rejection_reasons':reasons,'local_period_flags':flags}
            row['intermediate_strong_events']=[rows[k]['event_id'] for k in intermediate]
            row['warning_flags']=(['PHASE_SHIFT_POSSIBLE'] if phase_difference is not None and phase_difference>cfg['phase_shift_warning_delta'] else [])+(['WEAK_SHAPE_CORR'] if rho is not None and rho<cfg['effective_shape_correlation'] else [])
            audits.append(row)
            if not reasons:transitions[a,b]=score;connected[a,b]=True;edge_lookup[(a,b)]=row
    return ids,transitions,connected,audits,edge_lookup


def global_dp(node_scores,transitions,connected,cfg,forbidden=()):
    n=len(node_scores);K=cfg['top_paths'];scores=np.full((n,K),-np.inf);states=[[] for _ in range(n)]
    for j in range(n):
        if j in forbidden:continue
        proposals=[(float(node_scores[j]-cfg['root_penalty']),(),(),-1,-1)]
        if j:
            values=(scores[:j]+transitions[:j,j,None]+node_scores[j]).ravel();order=np.argsort(-values,kind='stable')
            for flat in order[:K*3]:
                if not np.isfinite(values[flat]):continue
                i,rank=divmod(int(flat),K)
                if rank>=len(states[i]):continue
                old=states[i][rank];proposals.append((float(values[flat]),old['nodes'],old['kinds'],i,rank))
        proposals.sort(key=lambda r:(-r[0],r[1],r[3],r[4]));seen=set()
        for score,prefix,kinds,i,rank in proposals:
            nodes=prefix+(j,)
            if nodes in seen:continue
            seen.add(nodes);state={'score':score,'nodes':nodes,'kinds':kinds+(('CONNECT' if connected[i,j] else 'RESTART'),) if i>=0 else ()}
            scores[j,len(states[j])]=score;states[j].append(state)
            if len(states[j])==K:break
    endings=[{'score':0.,'nodes':(),'kinds':()}]+[r for ss in states for r in ss]
    endings.sort(key=lambda r:(-r['score'],r['nodes'],r['kinds']));best=endings[:K]
    suffix=np.zeros(n)
    for i in range(n-1,-1,-1):
        if i+1<n:suffix[i]=max(0.,float(np.max(transitions[i,i+1:]+np.asarray(node_scores[i+1:])+suffix[i+1:])))
    return best,scores[:,0]+suffix


def global_path_track(rows,ids,path,edge_lookup,cfg):
    events=copy.deepcopy(rows);segments=[];current=[];attempts=[]
    for pos,node in enumerate(path['nodes']):
        j=ids[node];event=events[j]
        if pos and path['kinds'][pos-1]=='RESTART':
            if current:segments.append({'segment_id':len(segments)+1,'event_indices':current});current=[]
        current.append(j);event.update(accepted=True,reason='GLOBAL_ESTIMATED',anchor=len(current)==1)
        if pos and path['kinds'][pos-1]=='CONNECT':
            edge=edge_lookup[(path['nodes'][pos-1],node)];event['correlation']=edge['shape_correlation'];event['match_score']=edge['edge_score']
            attempts.append({'previous_event':events[ids[path['nodes'][pos-1]]]['event_id'],'candidate_event':event['event_id'],
                             'rho':edge['shape_correlation'],'score':edge['edge_score'],'accepted':True,'rank':1,'margin':None})
    if current:segments.append({'segment_id':len(segments)+1,'event_indices':current})
    for segment in segments:
        for j in segment['event_indices']:events[j]['segment_id']=segment['segment_id']
    return {'method':rows[0]['method'] if rows else 'A','prominence_percent':cfg['prominence_percent'],'events':events,'segments':segments,
            'attempts':attempts,'breaks':[],'tracking_method':'global-cycle','tracker_version':cfg['version']}


def global_time_partition(track,summary,start,end,data,frames):
    periods=summary['_periods'];intervals=[];cursor=start
    for r in sorted(periods,key=lambda r:r['start_s']):
        if r['start_s']>cursor:
            period=global_period_at((cursor+r['start_s'])/2,frames)[0]
            kind='WINDOW_EDGE' if cursor==start and period and r['start_s']-cursor<=1.65*period else 'UNRESOLVED'
            intervals.append({'start_s':cursor,'end_s':r['start_s'],'classification':kind,'reason':'before first completed period' if kind=='WINDOW_EDGE' else 'no legal coherent event edge'})
        ev={e['event_id']:e for e in track['events']};a,b=ev[r['event_start']],ev[r['event_end']]
        intervals.append({'start_s':r['start_s'],'end_s':r['end_s'],'classification':'ADOPTED_PERIOD',
            'quality':'HIGH_CONFIDENCE' if a['confidence_class']==b['confidence_class']=='HIGH_CONFIDENCE' else 'ESTIMATED',
            'ambiguous':a['ambiguous_assignment'] or b['ambiguous_assignment'] or 'PITCH_AMBIGUOUS' in a['local_period_flags'] or 'PITCH_AMBIGUOUS' in b['local_period_flags'],
            'event_start':r['event_start'],'event_end':r['event_end'],'segment_id':r['segment_id']})
        cursor=r['end_s']
    if cursor<end:
        period=global_period_at((cursor+end)/2,frames)[0];kind='WINDOW_EDGE' if periods and period and end-cursor<=1.65*period else 'UNRESOLVED'
        intervals.append({'start_s':cursor,'end_s':end,'classification':kind,'reason':'after last completed period' if kind=='WINDOW_EDGE' else 'no supported coherent event edge'})
    samples=data['t'][(data['t']>=start)&(data['t']<end)];membership=np.zeros(len(samples),int)
    for row in intervals:
        row.update(method=track['method'],duration_s=row['end_s']-row['start_s'],window_fraction=(row['end_s']-row['start_s'])/(end-start))
        mask=(samples>=row['start_s'])&(samples<row['end_s']);membership+=mask
        center=(row['start_s']+row['end_s'])/2;near=min(frames,key=lambda r:abs(r['center_s']-center))
        row['amplitude_present']=near['local_amplitude']>3*near['difference_MAD_noise_sigma'] and near['local_amplitude']>0
    if not np.all(membership==1) or not np.isclose(sum(r['duration_s'] for r in intervals),end-start):raise RuntimeError('global time partition has a gap or overlap')
    return intervals


def global_cycle_track(data,start,end,cfg,prepared=None,alternatives=True):
    began=time.perf_counter()
    if prepared is None:
        frames,period_candidates=global_local_periods(data,start,end,cfg);raw=global_raw_candidates(data,start,end,frames,cfg)
    else:frames,period_candidates,raw=prepared
    supported_frames=[f for f in frames if f['period_s'] is not None]
    period_status='PERIOD_AMBIGUOUS' if supported_frames and all('PITCH_AMBIGUOUS' in f['flags'] for f in supported_frames) else 'LOCAL_PERIOD_CANDIDATES' if supported_frames else 'NO_OSCILLATION' if all(f['status']=='NO_OSCILLATION' for f in frames) else 'APERIODIC_OR_UNRESOLVED'
    output={};all_edges=[]
    for phase in ('A','B'):
        rows,windows=global_candidate_features(data,raw,phase,cfg);ids,T,C,edge_rows,lookup=global_graph(rows,windows,frames,cfg)
        paths,marginal=global_dp([rows[j]['node_score'] for j in ids],T,C,cfg);path=paths[0];best_score=path['score'];selected=set(path['nodes'])
        margin=best_score-paths[1]['score'] if len(paths)>1 else None
        track=global_path_track(rows,ids,path,lookup,cfg);track['method']=phase
        chosen_edges=[lookup[(path['nodes'][k],path['nodes'][k+1])] for k,kind in enumerate(path['kinds']) if kind=='CONNECT']
        phase_warning_ids={identity for edge in chosen_edges if 'PHASE_SHIFT_POSSIBLE' in edge['warning_flags'] for identity in (edge['event_start'],edge['event_end'])}
        for node,j in enumerate(ids):
            e=track['events'][j];near=[k for k in range(len(ids)) if k not in selected and abs(rows[ids[k]]['raw_time_s']-e['raw_time_s'])<.55*(e['local_period_s'] or 0.) and marginal[k]>=best_score-cfg['ambiguity_score_margin']]
            e.update(path_score_margin=margin,alternative_event_times=[rows[ids[k]]['subsample_time_s'] for k in near],
                     ambiguous_assignment=bool(near),max_marginal_score=float(marginal[node]) if np.isfinite(marginal[node]) else None)
        for e in track['events']:
            e.setdefault('ambiguous_assignment',False);e.setdefault('alternative_event_times',[])
            linked=[r['rho'] for r in track['attempts'] if r['previous_event']==e['event_id'] or r['candidate_event']==e['event_id']]
            correlation=min(linked) if linked else None;e['minimum_adjacent_shape_correlation']=correlation
            high=bool(e['accepted'] and linked and min(linked)>=cfg['effective_shape_correlation'] and e['local_acf_confidence']>=cfg['high_acf_confidence'] and not e['ambiguous_assignment'] and not e['shape_flags'] and e['plateau_width_s']<=2/data['fs'] and 'PITCH_AMBIGUOUS' not in e['local_period_flags'] and e['event_id'] not in phase_warning_ids)
            e['confidence_class']='HIGH_CONFIDENCE' if high else 'ESTIMATED' if e['accepted'] else 'NOT_ADOPTED'
            e['selection_flags']=(['AMBIGUOUS_ASSIGNMENT'] if e['ambiguous_assignment'] else [])+(['WEAK_SHAPE_CORR'] if correlation is not None and correlation<cfg['effective_shape_correlation'] else [])+(['PLATEAU_TIMING_UNCERTAINTY'] if e['plateau_width_s']>2/data['fs'] else [])+(['PERIOD_AMBIGUOUS'] if 'PITCH_AMBIGUOUS' in e['local_period_flags'] else [])
            if e['event_id'] in phase_warning_ids:e['selection_flags'].append('PHASE_SHIFT_POSSIBLE')
            if e['accepted']:e['reason']='GLOBAL_'+e['confidence_class']
        summary=summarize_track(track,start,end,None,10,data['fs']);by_event={e['event_id']:e for e in track['events']}
        selected_lookup={(edge['event_start'],edge['event_end']):edge for edge in chosen_edges}
        for p in summary['_periods']:
            edge=selected_lookup[(p['event_start'],p['event_end'])];a,b=by_event[p['event_start']],by_event[p['event_end']]
            p.update(local_period_s=edge['local_period_s'],interval_over_local_period=edge['interval_over_local_period'],
                     shape_correlation=edge['shape_correlation'],phase_difference=edge['phase_difference'],
                     quality='HIGH_CONFIDENCE' if a['confidence_class']==b['confidence_class']=='HIGH_CONFIDENCE' else 'ESTIMATED',
                     period_identification_status=period_status,flags=sorted(set(a['selection_flags']+b['selection_flags'])))
        intervals=global_time_partition(track,summary,start,end,data,frames)
        metrics=copy.deepcopy(summary['pooled_within_segments_v2']);metrics['status']='PROVISIONAL' if metrics['jitter_local_percent'] is not None else 'UNDEFINED_JITTER'
        if period_status=='PERIOD_AMBIGUOUS':
            metrics['conditional_metrics']={k:metrics[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')};metrics['status']='PERIOD_AMBIGUOUS'
            for name in ('jitter_local_percent','rap_percent','ppq5_percent'):
                metrics[name]=None;metrics['metric_details'][name]['conditional_value']=metrics['metric_details'][name]['value'];metrics['metric_details'][name]['value']=None
        alternative_rows=[];alternative_periods=[];seen=set()
        def add_alternative(p,label):
            if p['nodes'] in seen:return
            seen.add(p['nodes']);tr=global_path_track(rows,ids,p,lookup,cfg);tr['method']=phase
            sm=summarize_track(tr,start,end,None,10,data['fs']);mm=sm['pooled_within_segments_v2']
            alternative_rows.append({'phase':phase,'path_id':label,'path_score':p['score'],'score_loss':best_score-p['score'],
                'event_times':[tr['events'][j]['subsample_time_s'] for seg in tr['segments'] for j in seg['event_indices']],
                'n_cycles':sm['n_complete_periods_all_segments'],'time_coverage':sm['covered_time_fraction'],
                **{k:mm[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')},'status':'PROVISIONAL_ALTERNATIVE_NOT_A_CONFIDENCE_INTERVAL'})
            alternative_periods.extend(dict(r,path_id=label) for r in sm['_periods'])
        if alternatives:
            for rank,p in enumerate(paths):add_alternative(p,f'top_{rank+1}')
            ambiguous=[node for node in path['nodes'] if track['events'][ids[node]]['ambiguous_assignment']]
            for node in ambiguous[:cfg['max_alternative_paths']]:
                alt,_=global_dp([rows[j]['node_score'] for j in ids],T,C,cfg,(node,));add_alternative(alt[0],'exclude_'+rows[ids[node]]['event_id'])
        for edge in edge_rows:edge.update(method=phase,adopted=any(edge['event_start']==p['event_start'] and edge['event_end']==p['event_end'] for p in summary['_periods']))
        all_edges.extend(edge_rows)
        adopted=[e for e in track['events'] if e['accepted']];duration=end-start
        output[phase]={'track':track,'summary':summary,'metrics':metrics,'intervals':intervals,'path_score':best_score,'path_score_margin':margin,
            'n_raw_candidates':len(rows),'n_eligible_candidates':len(ids),'dp_states':len(ids)*cfg['top_paths'],'matrix_bytes':T.nbytes+C.nbytes,'period_identification_status':period_status,
            'n_adopted_events':len(adopted),'high_confidence_event_fraction':sum(e['confidence_class']=='HIGH_CONFIDENCE' for e in adopted)/len(adopted) if adopted else None,
            'estimated_event_fraction':sum(e['confidence_class']=='ESTIMATED' for e in adopted)/len(adopted) if adopted else None,
            'time_coverage':summary['covered_time_fraction'],'unresolved_time_fraction':sum(r['duration_s'] for r in intervals if r['classification']=='UNRESOLVED')/duration,
            'ambiguous_time_fraction':sum(r['duration_s'] for r in intervals if r.get('ambiguous'))/duration,
            'amplitude_present_unresolved_fraction':sum(r['duration_s'] for r in intervals if r['classification']=='UNRESOLVED' and r['amplitude_present'])/duration,
            'classified_time_fraction':sum(r['duration_s'] for r in intervals)/duration,'alternative_paths':alternative_rows,'alternative_periods':alternative_periods,
            'alternative_exclusion_cap':cfg['max_alternative_paths'],'n_ambiguous_events':sum(e['ambiguous_assignment'] for e in adopted)}
        ratios=[r['interval_over_local_period'] for r in chosen_edges]
        output[phase].update(phase_shift_warning_event_fraction=sum('PHASE_SHIFT_POSSIBLE' in e['selection_flags'] for e in adopted)/len(adopted) if adopted else None,
            n_double_period_connections=sum(1.8<=r<=2.2 for r in ratios),n_duplicate_adopted_samples=len(adopted)-len({e['sample_id'] for e in adopted}),
            median_period_over_local_acf=float(np.median(ratios)) if ratios else None,p95_absolute_period_ratio_deviation=float(np.percentile(abs(np.array(ratios)-1),95)) if ratios else None)
        output[phase]['alternative_exclusions_attempted']=min(output[phase]['n_ambiguous_events'],cfg['max_alternative_paths']) if alternatives else 0
        output[phase]['alternative_exclusions_not_evaluated']=max(0,output[phase]['n_ambiguous_events']-output[phase]['alternative_exclusions_attempted'])
    return {'methods':output,'frames':frames,'period_candidates':period_candidates,'edge_audit':all_edges,'elapsed_s':time.perf_counter()-began,
            'profile':cfg,'status':'CALIBRATION_PENDING_MANUAL','_prepared':(frames,period_candidates,raw)}


def global_threshold_profiles(base,full=True):
    profiles=[];keys={}
    def add(label,**changes):
        cfg=dict(base,**changes);key=json.dumps({k:v for k,v in cfg.items() if k!='profile_name'},sort_keys=True)
        if key in keys:keys[key]['aliases'].append(label);return
        cfg['profile_name']=label;row={'profile_id':label,'aliases':[],'config':cfg};profiles.append(row);keys[key]=row
    add(base['profile_name'])
    if full:
        for p in (3,5,8,10,12,15,20):
            for amplitude in (3,5,10,15,20):add(f'prom{p}_amp{amplitude}',prominence_percent=float(p),excursion_percent=float(amplitude))
        for window in (.06,.1,.14):add(f'acf_{window:g}',acf_window_sec=window)
        add('strict_edge_range',search_min_period_ratio=.65,search_max_period_ratio=1.45)
        for rho in (.50,.55,.60,.65,.70):add(f'confidence_rho_{rho:g}',effective_shape_correlation=rho)
        add('MAD_floor_off',noise_sigma_multiplier=0.)
    return profiles


def resolve_global_profile(specification,args):
    cfg=dict(GLOBAL_CYCLE_CONFIG,pitch_min_hz=args.fmin,pitch_max_hz=args.fmax)
    if specification!='balanced-v1':
        path=Path(specification)
        if not path.is_file():raise ValueError('threshold profile must be balanced-v1 or an explicit JSON profile file')
        document=json.loads(path.read_text())
        if not isinstance(document,dict):raise ValueError('global profile must be a JSON object')
        supplied=document.get('config',document)
        if not isinstance(supplied,dict):raise ValueError('global config must be a JSON object')
        if set(supplied)-set(cfg):raise ValueError('unknown global profile keys')
        cfg.update(supplied)
    for key,default in GLOBAL_CYCLE_CONFIG.items():
        if isinstance(default,str):
            if not isinstance(cfg[key],str) or not cfg[key]:raise ValueError('global string settings must be nonempty')
        elif not isinstance(cfg[key],(int,float)) or isinstance(cfg[key],bool):raise ValueError('global numeric settings must be numbers')
    numeric=[v for v in cfg.values() if isinstance(v,(float,int))]
    if not np.isfinite(numeric).all():raise ValueError('global profile numbers must be finite')
    if not 0<cfg['pitch_min_hz']<cfg['pitch_max_hz'] or min(cfg['acf_window_sec'],cfg['acf_hop_sec'])<=0:raise ValueError('invalid global ACF range/window')
    if not all(0<=cfg[k]<=100 for k in ('prominence_percent','excursion_percent')):raise ValueError('percent thresholds must be in [0,100]')
    if not 0<cfg['search_min_period_ratio']<cfg['search_max_period_ratio']<1.8:raise ValueError('global edge interval must be positive and exclude ordinary 2T connections')
    if cfg['period_regularization_weight']!=0:raise ValueError('period regularization is disabled to avoid suppressing true jitter')
    if not 1<=cfg['max_dp_candidates']<=2000 or not 1<=cfg['top_paths']<=5:raise ValueError('invalid global search cap')
    if not all(isinstance(cfg[k],int) for k in ('max_dp_candidates','top_paths','max_alternative_paths')) or not 0<=cfg['max_alternative_paths']<=50:raise ValueError('global search caps must be bounded integers')
    if not -1<=cfg['incompatible_shape_correlation']<=cfg['effective_shape_correlation']<=1:raise ValueError('invalid global shape thresholds')
    if not 0<=cfg['phase_shift_warning_delta']<=1 or not 0<cfg['intermediate_pulse_strength_ratio']<=1:raise ValueError('invalid relative phase/strength setting')
    if any(cfg[k]<0 for k in cfg if k.endswith('_weight') or k.endswith('_penalty') or k=='noise_sigma_multiplier'):raise ValueError('score weights and noise multiplier must be nonnegative')
    return cfg


def global_profile_row(profile,phase,entry,reference):
    tr=entry['track'];events=[e for e in tr['events'] if e['accepted']];baseline=[e for e in reference['track']['events'] if e['accepted']];m=entry['metrics']
    row={'profile_id':profile['profile_id'],'method':phase,**profile['config'],
         'n_events':len(events),'n_cycles':entry['summary']['n_complete_periods_all_segments'],'time_coverage':entry['time_coverage'],
         'unresolved_time_fraction':entry['unresolved_time_fraction'],'ambiguous_time_fraction':entry['ambiguous_time_fraction'],
         'amplitude_present_unresolved_fraction':entry['amplitude_present_unresolved_fraction'],'classified_time_fraction':entry['classified_time_fraction'],
         'high_confidence_event_fraction':entry['high_confidence_event_fraction'],'estimated_event_fraction':entry['estimated_event_fraction'],
         'secondary_lobe_event_fraction':sum(e['secondary_lobe'] for e in events)/len(events) if events else None,
         'event_uniqueness_fraction':1-entry['n_ambiguous_events']/len(events) if events else None,
         'path_score_margin':entry['path_score_margin'],'status':'CALIBRATION_PENDING_MANUAL',
         'metric_status':m['status'],
         **{k:entry[k] for k in ('phase_shift_warning_event_fraction','n_double_period_connections','n_duplicate_adopted_samples','median_period_over_local_acf','p95_absolute_period_ratio_deviation')},
         **{k:m[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')}}
    for name,detail in m['metric_details'].items():
        for key in ('n_cycles','n_terms','mean_period_s','covered_time_fraction'):row[name+'_'+key]=detail[key]
    for tolerance in (50,100,150):
        match=optimal_event_assignment([e['subsample_time_s'] for e in baseline],[e['subsample_time_s'] for e in events],tolerance*1e-6)
        row[f'event_correspondence_{tolerance}us']=2*len(match)/(len(baseline)+len(events)) if baseline or events else None
    return row


def manual_common_period_comparison(manual,tracks,matches,start,end):
    """Compare the same matched consecutive periods, not merely the same time extent."""
    statistics=[];period_rows=[]
    for phase in ('A','B'):
        refs=[r for r in manual if r['phase']==phase]
        if not refs:continue
        edges={}
        for segment in tracks[phase]['segments']:
            ev=[tracks[phase]['events'][j] for j in segment['event_indices']]
            for a,b in zip(ev,ev[1:]):edges[(a['event_id'],b['event_id'])]=(a,b,segment['segment_id'])
        for scope in ('high_only','all_annotations'):
            for tolerance in (50,100,150):
                mapping={r['manual_time_s']:r['automatic_event_id'] for r in matches if r['phase']==phase and r['scope']==scope and abs(r['tolerance_us']-tolerance)<1e-8 and r['classification']=='TP'}
                blocks=[];block=[];previous=None
                for a,b in zip(refs,refs[1:]):
                    eligible=(scope=='all_annotations' or a['confidence']==b['confidence']=='high') and a['event_time_s'] in mapping and b['event_time_s'] in mapping
                    pair=(mapping.get(a['event_time_s']),mapping.get(b['event_time_s']));edge=edges.get(pair)
                    if eligible and edge:
                        aa,bb,sid=edge
                        if not (start<=aa['subsample_time_s']<bb['subsample_time_s']<end):eligible=False
                    if not eligible or edge is None:
                        if block:blocks.append(block);block=[]
                        previous=None;continue
                    if previous is not None and (previous['manual_end_s']!=a['event_time_s'] or previous['event_end']!=pair[0] or previous['segment_id']!=sid):
                        if block:blocks.append(block);block=[]
                    row={'phase':phase,'scope':scope,'tolerance_us':tolerance,'manual_start_s':a['event_time_s'],'manual_end_s':b['event_time_s'],
                         'auto_start_s':aa['subsample_time_s'],'auto_end_s':bb['subsample_time_s'],'event_start':pair[0],'event_end':pair[1],'segment_id':sid,
                         'manual_period_s':b['event_time_s']-a['event_time_s'],'auto_period_s':bb['subsample_time_s']-aa['subsample_time_s']}
                    block.append(row);previous=row
                if block:blocks.append(block)
                m=metric_blocks_v2([[r['manual_period_s'] for r in b] for b in blocks],end-start);a=metric_blocks_v2([[r['auto_period_s'] for r in b] for b in blocks],end-start)
                statistics.append({'phase':phase,'scope':scope,'tolerance_us':tolerance,'comparison_basis':'same TP-matched consecutive period IDs; no uncertain/unmatched/segment-gap bridge',
                    **{'manual_common_'+k:v for k,v in m.items()},**{'auto_common_'+k:v for k,v in a.items()},
                    **{name+'_difference_pp':a[name]-m[name] if a[name] is not None and m[name] is not None else None for name in ('jitter_local_percent','rap_percent','ppq5_percent')}})
                for number,block in enumerate(blocks):
                    for row in block:period_rows.append(dict(row,common_block_id=number,used_for_jitter_v2=len(block)>=2,used_for_RAP_v2=len(block)>=3,used_for_PPQ5_v2=len(block)>=5))
    return statistics,period_rows


def global_review_figures(out,data,start,end,result,legacy,recover,case_id,review_start=None,review_end=None):
    t=data['t'];x=data['x'];mask=(t>=start)&(t<end);methods=result['methods'];frames=result['frames']
    def overlay(ax,left,right,blind=False,labels=False):
        take=(t>=left)&(t<right);ax.plot(t[take],x[take]*1e6,color='black',lw=.8,label='raw Ug; no smoothing')
        if not blind:
            for phase,color,marker in [('A','tab:blue','o'),('B','tab:green','v')]:
                tr=methods[phase]['track']
                for adopted in (False,True):
                    ev=[e for e in tr['events'] if e['accepted']==adopted and left<=e['subsample_time_s']<right];at=[e['subsample_time_s'] for e in ev]
                    ax.scatter(at,np.interp(at,t,x)*1e6,marker=marker if adopted else 'x',color=color if adopted else 'gray',s=23 if adopted else 9,alpha=1 if adopted else .45,label=phase+(' adopted' if adopted else ' raw/rejected'))
                    if labels and adopted:
                        for e in ev:
                            height=np.interp(e['subsample_time_s'],t,x)*1e6
                            ax.annotate(e['event_id'],(e['subsample_time_s'],height),xytext=(0,7 if phase=='A' else -12),textcoords='offset points',fontsize=6,rotation=30)
                            ax.vlines(e['subsample_time_s'],height-(1 if phase=='A' else -1)*e['effective_prominence_threshold']*1e6,height,color=color,alpha=.5,lw=2)
                for row in methods[phase]['intervals']:
                    if row['classification']=='UNRESOLVED' and row['end_s']>left and row['start_s']<right:ax.axvspan(max(left,row['start_s']),min(right,row['end_s']),color=color,alpha=.10)
                if labels:
                    for p in methods[phase]['summary']['_periods']:
                        if left<=p['start_s'] and p['end_s']<=right:
                            ax.annotate(f'{p["period_s"]*1e3:.2f}ms',((p['start_s']+p['end_s'])/2,ax.get_ylim()[1]*(.95 if phase=='A' else .85)),ha='center',fontsize=6,color=color)
            ax.legend(fontsize=7,ncol=3,loc='upper left')
        ax.set(xlim=(left,right),xlabel='Time [s]',ylabel='Ug [mL/s]');ax.grid(alpha=.12)
    fig,axes=plt.subplots(2,1,figsize=(16,8),sharex=True,gridspec_kw={'height_ratios':[3,1]});overlay(axes[0],start,end)
    for phase,color in [('A','tab:orange'),('B','tab:red')]:
        tr=legacy['_tracks'][legacy['_primary_key'] if phase=='A' else legacy['_control_key']];ev=[e['subsample_time_s'] for e in tr['events'] if e['accepted'] and e['in_window']]
        axes[0].scatter(ev,np.interp(ev,t,x)*1e6,s=40,facecolors='none',edgecolors=color,linewidths=.6,label='legacy '+phase)
    axes[0].legend(fontsize=7,ncol=4);axes[0].set_title(case_id+' full-window raw candidates and global events; shaded = unresolved')
    axes[1].plot([r['center_s'] for r in frames],[r['period_s']*1e3 if r['period_s'] else np.nan for r in frames],'o-',ms=3);axes[1].set(ylabel='Local T0 [ms]',xlabel='Time [s]');axes[1].grid(alpha=.2);fig.tight_layout();fig.savefig(out/'flow_full_window_events.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(13,6),sharex=True)
    for row in result['period_candidates']:axes[0].scatter(row['center_s'],row['period_s']*1e3,c='tab:blue' if row['selected'] else 'gray',s=16 if row['selected'] else 5,alpha=.8 if row['selected'] else .3)
    for factor,style in [(.5,'--'),(1.,'-'),(2.,':')]:axes[0].plot([r['center_s'] for r in frames],[r['period_s']*factor*1e3 if r['period_s'] else np.nan for r in frames],style,label=f'{factor:g}T0')
    axes[0].set_ylabel('Period candidates [ms]');axes[0].legend();axes[1].plot([r['center_s'] for r in frames],[r['acf_confidence'] for r in frames],'o-');axes[1].set(xlabel='Time [s]',ylabel='Normalized ACF',ylim=(-.05,1.05))
    for ax in axes:ax.grid(alpha=.15);ax.set_xlim(start,end)
    fig.suptitle(case_id+' local ACF: search prediction only; alternatives are not event times');fig.tight_layout();fig.savefig(out/'acf_timeline.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(15,4.5),sharex=True)
    for ax,phase in zip(axes,('A','B')):
        for row in methods[phase]['intervals']:
            color='tab:green' if row.get('quality')=='HIGH_CONFIDENCE' else 'tab:orange' if row['classification']=='ADOPTED_PERIOD' else 'tab:red' if row['classification']=='UNRESOLVED' else 'gray'
            ax.axvspan(row['start_s'],row['end_s'],color=color,alpha=.65)
        ax.set(yticks=[],ylabel=phase,xlim=(start,end));ax.set_title(f'{phase}: green=high confidence, orange=estimated, red=unresolved, gray=edge; classified=100%')
    axes[-1].set_xlabel('Time [s]');fig.tight_layout();fig.savefig(out/'event_timeline.png',dpi=160);plt.close(fig)
    panels=out/'event_panels';panels.mkdir();focus=out/'unresolved_focus';focus.mkdir();index=[]
    left=start if review_start is None else review_start;right=end if review_end is None else review_end
    if not start<=left<right<=end:raise ValueError('manual review window outside analysis extent')
    windows=[];cursor=left;number=0;typical=np.median([r['period_s'] for r in frames if r['period_s']]) if any(r['period_s'] for r in frames) else .008
    while cursor<right-1e-12:
        stop=min(cursor+.060,right)
        if stop-cursor<.05 and number:cursor=max(left,right-.060)
        number+=1;windows.append((f'panel_{number:03d}',cursor,stop,panels,'full-window'))
        if stop>=right:break
        cursor=stop-min(typical,.015)
    important=[(.574,.581),(.591,.600),(.607,.621),(.790,.806)] if case_id=='r0p25' else [(.549,.565),(.681,.697),(.749,.767)] if case_id=='r0p30' else []
    for phase in ('A','B'):
        for row in methods[phase]['intervals']:
            if row['classification']=='UNRESOLVED':important.append((row['start_s'],row['end_s']))
    for k,(a,b) in enumerate(important,1):windows.append((f'focus_{k:03d}',max(start,a-5*typical),min(end,b+5*typical),focus,'unresolved / required regression'))
    # Every event newly adopted versus recover-rejected receives a dedicated neighborhood.
    for phase in ('A','B'):
        old={e['event_id'] for e in recover['_tracks'][recover['_primary_key'] if phase=='A' else recover['_control_key']]['events'] if e['accepted']}
        added=[e for e in methods[phase]['track']['events'] if e['accepted'] and e['event_id'] not in old]
        for e in added:windows.append(('new_'+e['event_id'],max(start,e['subsample_time_s']-5*typical),min(end,e['subsample_time_s']+5*typical),focus,'new event vs recover-rejected; requires review'))
        old_pairs={(r['event_start'],r['event_end']) for r in recover['_summaries'][recover['_primary_key'] if phase=='A' else recover['_control_key']]['_periods']}
        for p in methods[phase]['summary']['_periods']:
            if (p['event_start'],p['event_end']) not in old_pairs and not any(a<=p['start_s'] and p['end_s']<=b for _,a,b,_,_ in windows):
                windows.append(('edge_'+p['event_start']+'_'+p['event_end'],max(start,p['start_s']-5*typical),min(end,p['end_s']+5*typical),focus,'new connection vs recover-rejected'))
    for name,a,b,directory,note in windows:
        for blind in (True,False):
            fig,ax=plt.subplots(figsize=(15,6));overlay(ax,a,b,blind,not blind);ax.set_title(f'{case_id} {name} [{a:.6f}, {b:.6f}] '+('BLIND' if blind else 'DIAGNOSTIC'))
            if not blind:
                visible=[e for phase in ('A','B') for e in methods[phase]['track']['events'] if e['accepted'] and a<=e['subsample_time_s']<b]
                lines=[f'{e["event_id"]}: Alocal={e["local_amplitude"]*1e6:.1f}mL/s, p/A={e["normalized_prominence"]:.3f}, dU/A={e["normalized_excursion"]:.3f}, T={e["local_period_s"]*1e3:.2f}ms, {e["confidence_class"]}' for e in visible]
                split=(len(lines)+1)//2;bottom=min(.60,max(.20,.018*split+.08));fig.subplots_adjust(left=.06,right=.99,top=.90,bottom=bottom)
                fig.text(.06,.025,'\n'.join(lines[:split]),fontsize=6,va='bottom');fig.text(.53,.025,'\n'.join(lines[split:]),fontsize=6,va='bottom')
            else:fig.tight_layout()
            fig.savefig(directory/f'{name}_{"blind" if blind else "diagnostic"}.png',dpi=140);plt.close(fig)
        index.append({'case_id':case_id,'panel_id':name,'t_start_s':a,'t_end_s':b,'figure_blind':str((directory/(name+'_blind.png')).relative_to(out)),
                      'figure_diagnostic':str((directory/(name+'_diagnostic.png')).relative_to(out)),'notes':note})
    write_csv(out/'manual_review_index.csv',index);write_csv(out/'manual_events_template.csv',[],['case_id','phase','event_time_s','confidence','reviewer','notes'])
    # A compact selection example contains actual raw events and measured intervals.
    fig,axes=plt.subplots(2,1,figsize=(14,8))
    examples=[(.59,.63),(.79,.83)] if start<=.59 and end>=.83 else [(start+.15*(end-start),start+.35*(end-start)),(start+.60*(end-start),start+.80*(end-start))]
    for ax,(a,b) in zip(axes,examples):overlay(ax,a,b,labels=True)
    fig.suptitle(case_id+' event selection examples: secondary candidates retained, actual intervals labeled');fig.tight_layout();fig.savefig(out/'event_selection_examples.png',dpi=160);plt.close(fig)
    return index


def global_threshold_figure(out,rows):
    fig,axes=plt.subplots(2,4,figsize=(17,8));ps=[3,5,8,10,12,15,20];amps=[3,5,10,15,20]
    for row,phase in enumerate(('A','B')):
        for col,key in enumerate(('time_coverage','event_uniqueness_fraction','jitter_local_percent','ppq5_percent')):
            values=np.full((len(amps),len(ps)),np.nan)
            for iy,a in enumerate(amps):
                for ix,p in enumerate(ps):
                    r=next((r for r in rows if r['method']==phase and r['prominence_percent']==p and r['excursion_percent']==a and r['acf_window_sec']==.1 and r['search_min_period_ratio']==.55 and r['effective_shape_correlation']==.65 and r['noise_sigma_multiplier']==2.),None)
                    if r and r[key] is not None:values[iy,ix]=r[key]
            fraction=key in ('time_coverage','event_uniqueness_fraction')
            im=axes[row,col].imshow(values,origin='lower',aspect='auto',vmin=0 if fraction else None,vmax=1 if fraction else None);axes[row,col].set(xticks=range(len(ps)),xticklabels=ps,yticks=range(len(amps)),yticklabels=amps,xlabel='Local prominence [%]',ylabel='Excursion [%]',title=phase+' '+key);fig.colorbar(im,ax=axes[row,col],shrink=.8)
    fig.suptitle('Shared threshold grid: coverage and correspondence are not accuracy; Jitter is not minimized');fig.tight_layout();fig.savefig(out/'threshold_grid.png',dpi=155);plt.close(fig)


def global_comparison_figure(out,rows):
    groups=sorted({(r['case_id'],r['phase']) for r in rows});fig,axes=plt.subplots(len(groups),4,figsize=(16,3.3*len(groups)),squeeze=False)
    labels={'jitter_local_percent':'Jitter local [%]','rap_percent':'RAP [%]','ppq5_percent':'PPQ5 [%]','flow_CPP_dB':'Flow CPP [dB]'}
    fig.subplots_adjust(left=.07,right=.98,bottom=.10,top=.92,wspace=.35,hspace=1.0)
    for line,(case,phase) in enumerate(groups):
        subset=[r for r in rows if r['case_id']==case and r['phase']==phase]
        for col,key in enumerate(('jitter_local_percent','rap_percent','ppq5_percent','flow_CPP_dB')):
            values=[r[key] if r[key] is not None else np.nan for r in subset];axes[line,col].bar(range(len(subset)),values,color=['gray','tab:blue','tab:orange'][:len(subset)])
            axes[line,col].set(xticks=range(len(subset)),xticklabels=['Legacy\nlongest','Recover\npooled','Global\nestimated'][:len(subset)],title=case+' '+phase+' '+labels[key])
            axes[line,col].grid(axis='y',alpha=.15)
        note='; '.join(f'{name}: tracked {r["tracked_time_fraction"]:.1%}, J-used {r["jitter_local_percent_covered_time_fraction"]:.1%}' for name,r in zip(('Legacy','Recover','Global'),subset))
        fig.text(.07,axes[line,0].get_position().y0-.072,note,fontsize=8)
    fig.suptitle('Event/aggregation comparison: global estimates provisional; CPP unchanged',y=.985);fig.savefig(out/'jitter_rap_ppq5_cpp_comparison.png',dpi=155);plt.close(fig)


def global_cycle_cli(args):
    if args.input is None:raise ValueError('--input is required')
    if args.min_cycles<1:raise ValueError('positive research min-cycles required')
    if args.primary_prominence!=10 or args.prominences!=TRACKING_CONFIG['prominence_sweep_percent']:raise ValueError('global-cycle thresholds are configured by --threshold-profile; legacy prominence flags must keep their defaults')
    began=time.perf_counter();data=read_flow(args.input,resolve_dt(args.dt,args.manifest));start=data['t'][0] if args.start is None else args.start;end=data['t'][-1]+data['sampling_interval_s'] if args.end is None else args.end
    if not np.isfinite([start,end]).all() or start>=end or start<data['t'][0] or end>data['t'][-1]+data['sampling_interval_s']+1e-12:raise ValueError('global window outside saved data support')
    if np.count_nonzero((data['t']>=start)&(data['t']<end))<3:raise ValueError('global window needs at least three saved samples')
    if args.cpp_frame_length!=CPP_CONFIG['frame_sec'] or args.cpp_hop!=CPP_CONFIG['hop_sec']:raise ValueError('global validation keeps CPP settings fixed')
    cfg=resolve_global_profile(args.threshold_profile,args);case_id=args.case_id or next((part.split('_P')[0] for part in Path(args.input).parts if part.startswith('r') and '_P' in part),Path(args.input).stem)
    out=(args.out_dir or Path('analysis_runs')/('global_cycle_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve()
    protected=[Path(args.input).resolve()]+[Path(p).resolve() for p in (args.manifest,args.manual_events,args.diagnostic_baseline,args.reference_fullrate) if p]
    if any(p==out or out in p.parents for p in protected):raise ValueError('global output contains protected input')
    if out.exists():raise ValueError('global-cycle requires a new --out-dir, even with --overwrite')
    tracking=dict(TRACKING_CONFIG,pitch_min_hz=cfg['pitch_min_hz'],pitch_max_hz=cfg['pitch_max_hz'])
    legacy=analyze_arrays(data['t'],data['x'],data['fs'],start,end,tracking,CPP_CONFIG,args.min_cycles)
    recover=analyze_arrays(data['t'],data['x'],data['fs'],start,end,dict(tracking,reanchor_policy='recover-rejected'),CPP_CONFIG,args.min_cycles)
    if clean_json(legacy['cpp'])!=clean_json(recover['cpp']):raise RuntimeError('CPP regression mismatch')
    if args.diagnostic_baseline:
        prior=json.loads(args.diagnostic_baseline.read_text());meta=prior['metadata']
        if meta['input_sha256']!=data['sha256'] or (meta['dt_s'],meta['start_s'],meta['end_s'])!=(data['dt_s'],start,end):raise ValueError('global baseline input/time mismatch')
        source=recover if meta.get('reanchor_policy')=='recover-rejected' else legacy
        for key in ('primary','control'):
            if source[key]['main']['event_ids']!=prior[key]['main']['event_ids']:raise RuntimeError('old baseline event IDs changed')
            for metric in ('jitter_local_percent','rap_percent','ppq5_percent'):
                a,b=source[key]['pooled_within_segments_v2'][metric],prior[key]['pooled_within_segments_v2'][metric]
                if (a is None)!=(b is None) or a is not None and not np.isclose(a,b,rtol=1e-10,atol=1e-10):raise RuntimeError('old baseline metric changed')
        if not np.isclose(legacy['cpp']['flow_CPP_dB'],prior['cpp']['flow_CPP_dB'],rtol=1e-9,atol=1e-9):raise RuntimeError('stored CPP changed')
    main=global_cycle_track(data,start,end,cfg);profiles=global_threshold_profiles(cfg,not args.skip_threshold_sweep);cache={cfg['acf_window_sec']:main['_prepared']};sweep=[];variant_events=[];variant_periods=[];manual_sweep=[]
    manual=read_manual_reference(args.manual_events,case_id,data,start,end) if args.manual_events else []
    for profile in profiles:
        pcfg=profile['config'];window=pcfg['acf_window_sec']
        if window not in cache:
            frames,candidates=global_local_periods(data,start,end,pcfg);cache[window]=(frames,candidates,global_raw_candidates(data,start,end,frames,pcfg))
        run=main if profile['profile_id']==cfg['profile_name'] else global_cycle_track(data,start,end,pcfg,cache[window],False)
        for phase in ('A','B'):
            v=run['methods'][phase];sweep.append(global_profile_row(profile,phase,v,main['methods'][phase]))
            variant_events.extend(dict(e,profile_id=profile['profile_id']) for e in v['track']['events'])
            variant_periods.extend(dict(p,profile_id=profile['profile_id']) for p in v['summary']['_periods'])
        if manual:
            track_map={p:run['methods'][p]['track'] for p in ('A','B')}
            stats,matches,_,_=evaluate_manual_reference(manual,track_map,start,end)
            common,_=manual_common_period_comparison(manual,track_map,matches,start,end)
            for row in stats:
                linked=next((r for r in common if r['phase']==row['phase'] and r['scope']==row['scope'] and abs(r['tolerance_us']-row['tolerance_us'])<1e-8),{})
                row.update({k:v for k,v in linked.items() if k not in ('phase','scope','tolerance_us')})
            manual_sweep.extend(dict(r,profile_id=profile['profile_id']) for r in stats)
    comparison=[];alignments=[];changes=[];new_connections=[];subsamples=[];plateau_sensitivity=[];same_windows=[];same_periods=[]
    for phase,key,private in [('A','primary','_primary_key'),('B','control','_control_key')]:
        entry=main['methods'][phase];tr=entry['track'];sm=entry['summary'];baseline=recover['_tracks'][recover[private]]
        for label,source,scope in [('longest_legacy',legacy,legacy[key]['main']),('recover_rejected_pooled_v2',recover,recover[key]['pooled_within_segments_v2']),('global_cycle_estimated',None,entry['metrics'])]:
            row={'case_id':case_id,'phase':phase,'method':label,'flow_CPP_dB':legacy['cpp']['flow_CPP_dB'],'CPP_frames':legacy['cpp']['n_valid_frames'],
                 'tracked_time_fraction':entry['time_coverage'] if source is None else source[key]['covered_time_fraction'],
                 **{k:scope[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')},
                 'status':entry['metrics']['status'] if source is None else source[key]['status']}
            for width,name in [(2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')]:
                detail=scope.get('metric_details',{}).get(name)
                if detail:
                    for attr in ('n_cycles','n_terms','mean_period_s','covered_time_fraction'):row[name+'_'+attr]=detail[attr]
                else:
                    row[name+'_n_cycles']=scope['n_cycles'] if scope['n_cycles']>=width else 0;row[name+'_n_terms']=max(0,scope['n_cycles']-width+1);row[name+'_mean_period_s']=scope['mean_period_s'] if scope['n_cycles']>=width else None
                    row[name+'_covered_time_fraction']=scope['duration_s']/(end-start) if scope['n_cycles']>=width else 0.
            row.update(unresolved_time_fraction=entry['unresolved_time_fraction'] if source is None else None,
                       ambiguous_time_fraction=entry['ambiguous_time_fraction'] if source is None else None,
                       classified_time_fraction=entry['classified_time_fraction'] if source is None else None,
                       high_confidence_event_fraction=entry['high_confidence_event_fraction'] if source is None else None,
                       estimated_event_fraction=entry['estimated_event_fraction'] if source is None else None)
            comparison.append(row)
        old_ids={e['event_id']:e for e in baseline['events'] if e['accepted'] and e['in_window']};new_ids={e['event_id']:e for e in tr['events'] if e['accepted']}
        for identity in sorted(set(old_ids)|set(new_ids)):
            e=new_ids.get(identity,old_ids.get(identity));changes.append({'phase':phase,'event_id':identity,'time_s':e['subsample_time_s'],'recover_adopted':identity in old_ids,'global_adopted':identity in new_ids,'change':'SHARED' if identity in old_ids and identity in new_ids else 'NEW' if identity in new_ids else 'REMOVED'})
        high=[e for e in baseline['events'] if e['accepted'] and e['in_window'] and e['correlation'] is not None and e['correlation']>=.85]
        for tol in (50,100,150):
            pairs=optimal_event_assignment([e['subsample_time_s'] for e in high],[e['subsample_time_s'] for e in new_ids.values()],tol*1e-6)
            alignments.append({'phase':phase,'tolerance_us':tol,'recover_high_confidence_events':len(high),'matched':len(pairs),'high_confidence_reference_alignment':len(pairs)/len(high) if high else None,'not_manual_accuracy':True})
        old_pairs={(p['event_start'],p['event_end']) for p in recover['_summaries'][recover[private]]['_periods']}
        for p in sm['_periods']:
            if (p['event_start'],p['event_end']) not in old_pairs:new_connections.append(dict(p,review_status='PENDING_MANUAL_ANNOTATION'))
        for mode in ('quadratic','off'):
            clone=copy.deepcopy(tr)
            if mode=='off':
                for e in clone['events']:e['subsample_time_s']=e['raw_time_s']
            blocks,periods=diagnostic_clip_blocks(clone,start,end);subsamples.append({'phase':phase,'subsample':mode,'same_event_ids':True,'period_identification_status':entry['period_identification_status'],**metric_blocks_v2(blocks,end-start)})
        if phase=='B':
            for representative in ('plateau_start_s','raw_time_s','plateau_end_s'):
                clone=copy.deepcopy(tr)
                for e in clone['events']:e['subsample_time_s']=e[representative]
                blocks,periods=diagnostic_clip_blocks(clone,start,end);plateau_sensitivity.append({'phase':'B','representative':representative,'actual_same_plateau_event_ids':True,'period_identification_status':entry['period_identification_status'],**metric_blocks_v2(blocks,end-start)})
        windows=[('full_requested',start,end),('early_half',start,(start+end)/2),('late_half',(start+end)/2,end)]
        for p in ('A','B'):
            m=recover['primary' if p=='A' else 'control']['main']
            if m['n_cycles']:windows.append(('recover_'+p+'_longest',m['start_s'],m['end_s']))
        for name,left,right in windows:
            for label,track in [('recover-rejected',baseline),('global-cycle',tr)]:
                blocks,periods=diagnostic_clip_blocks(track,left,right);same_windows.append({'phase':phase,'method':label,'window_name':name,'start_s':left,'end_s':right,**metric_blocks_v2(blocks,right-left)})
                same_periods.extend(dict(p,method_name=label,window_name=name) for p in periods)
    manual_status='MANUAL_REFERENCE_EVALUATED' if manual else 'PENDING_MANUAL_ANNOTATION'
    manual_stats,manual_matches,manual_checks,manual_periods=evaluate_manual_reference(manual,{p:main['methods'][p]['track'] for p in ('A','B')},start,end) if manual else ([],[],[],[])
    common_manual_stats,common_manual_periods=manual_common_period_comparison(manual,{p:main['methods'][p]['track'] for p in ('A','B')},manual_matches,start,end) if manual else ([],[])
    sampling=[{'status':'NOT_AVAILABLE','reason':'no actual stride-1 reference; no synthetic upsampling claimed'}]
    if args.reference_fullrate:
        ref=read_flow(args.reference_fullrate,data['dt_s']);common=np.searchsorted(ref['steps'],data['steps'])
        if ref['stride']!=1 or np.any(common>=len(ref['steps'])) or not np.array_equal(ref['steps'][common],data['steps']) or not np.allclose(ref['x'][common],data['x'],rtol=1e-10,atol=1e-14):raise ValueError('invalid actual full-rate reference')
        sampling=[]
        for shift in (-1,0,1,2,3,4):
            dd=ref if shift==-1 else dict(ref,t=ref['t'][shift::5],x=ref['x'][shift::5],fs=ref['fs']/5)
            rr=global_cycle_track(dd,start,end,cfg,alternatives=False)
            for phase,v in rr['methods'].items():sampling.append({'status':'ACTUAL_FULLRATE_REFERENCE','phase':phase,'decimation_phase':shift,'sample_rate_hz':dd['fs'],**{k:v['metrics'][k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')}})
    summary={'metadata':{'case_id':case_id,'input_path':data['path'],'input_sha256':data['sha256'],'dt_s':data['dt_s'],'step_stride':data['stride'],
        'sample_rate_hz':data['fs'],'start_s':start,'end_s':end,'window_semantics':'[start,end)','tracking_method':'global-cycle','script_sha256':sha256(__file__),
        'created_utc':datetime.now(timezone.utc).isoformat(),'settings_sha256':hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest(),
        'elapsed_s':time.perf_counter()-began,'sampling_reference_status':sampling[0]['status']},
        'tracking_config_global_cycle':cfg,'calibration_status':'CALIBRATION_PENDING_MANUAL','manual_validation_status':manual_status,
        'manual_reference':{'path':str(args.manual_events.resolve()),'sha256':sha256(args.manual_events),'n_validated_events':len(manual)} if args.manual_events else None,
        'manual_validation':manual_stats,
        'manual_common_period_validation':common_manual_stats,
        'global_cycle_estimated':{p:{k:v for k,v in entry.items() if k not in ('track','summary','alternative_periods')} for p,entry in main['methods'].items()},
        'cpp':legacy['cpp'],'cpp_config':CPP_CONFIG,'CPP_definition':CPP_DEFINITION,
        'baseline_legacy':{p:legacy['primary' if p=='A' else 'control'] for p in ('A','B')},
        'baseline_recover_rejected':{p:recover['primary' if p=='A' else 'control'] for p in ('A','B')},
        'comparison':comparison,'sampling_comparison':sampling,'search_notes':{'full_window_scanned':True,'period_prediction_replaces_event_times':False,
        'regularization_weight':0.,'scores_fixed_before_sweep':True,'high_confidence_is_not_calibrated_accuracy':True,'threshold_profiles':len(profiles),
        'sweep_design':'full 7x5 local prominence/excursion grid plus one-factor ACF/range/confidence/noise tests; no automatic best-Jitter profile selection'}}
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.airflow_global_',dir=out.parent) as directory:
        temp=Path(directory);dump_json(temp/'summary.json',summary);dump_json(temp/'analysis_metadata.json',summary['metadata']);dump_json(temp/'selected_profile.json',{'config':cfg});dump_json(temp/'threshold_profiles.json',profiles)
        write_csv(temp/'threshold_profiles.csv',[dict(profile_id=p['profile_id'],aliases=p['aliases'],**p['config']) for p in profiles]);write_csv(temp/'threshold_comparison.csv',sweep)
        write_csv(temp/'threshold_events.csv',variant_events);write_csv(temp/'threshold_periods.csv',variant_periods)
        write_csv(temp/'local_period_frames.csv',main['frames']);write_csv(temp/'local_period_candidates.csv',main['period_candidates']);write_csv(temp/'edge_score_audit.csv',main['edge_audit'])
        candidates=[dict(e,case_id=case_id) for phase in ('A','B') for e in main['methods'][phase]['track']['events']];write_csv(temp/'event_candidates_all.csv',candidates);write_csv(temp/'candidate_score_audit.csv',candidates)
        classifications=[];metric_terms=[];alternatives_rows=[];alternatives_periods=[]
        for phase,entry in main['methods'].items():
            event_columns=list(dict.fromkeys(k for e in candidates for k in e)) if candidates else ['event_id','subsample_time_s','case_id']
            write_csv(temp/f'global_events_{phase}.csv',[dict(e,case_id=case_id) for e in entry['track']['events'] if e['accepted']],event_columns)
            write_csv(temp/f'global_periods_{phase}.csv',entry['summary']['_periods'],None if entry['summary']['_periods'] else ['method','segment_id','event_start','event_end','start_s','end_s','period_s','used_for_jitter_v2','used_for_RAP_v2','used_for_PPQ5_v2']);classifications.extend(entry['intervals']);metric_terms.extend(entry['metrics']['_terms']);alternatives_rows.extend(entry['alternative_paths']);alternatives_periods.extend(entry['alternative_periods'])
        write_csv(temp/'time_classification.csv',classifications);write_csv(temp/'unresolved_intervals.csv',[r for r in classifications if r['classification']=='UNRESOLVED'],list(classifications[0]))
        write_csv(temp/'metric_terms_global.csv',metric_terms);write_csv(temp/'alternative_paths.csv',alternatives_rows);write_csv(temp/'alternative_periods.csv',alternatives_periods)
        write_csv(temp/'method_comparison.csv',comparison);write_csv(temp/'high_confidence_alignment.csv',alignments);write_csv(temp/'event_change_comparison.csv',changes)
        write_csv(temp/'global_subsample_comparison.csv',subsamples);write_csv(temp/'plateau_representative_sensitivity.csv',plateau_sensitivity);write_csv(temp/'sampling_comparison.csv',sampling)
        write_csv(temp/'same_window_metrics.csv',same_windows);write_csv(temp/'same_window_periods.csv',same_periods)
        write_csv(temp/'cpp_frames.csv',legacy['cpp']['rows']);dump_json(temp/'baseline_legacy.json',legacy);dump_json(temp/'baseline_recover_rejected.json',recover)
        for label,baseline_result in [('legacy',legacy),('recover-rejected',recover)]:
            baseline_events=[];baseline_periods=[]
            for private in ('_primary_key','_control_key'):
                baseline_events.extend(baseline_result['_tracks'][baseline_result[private]]['events']);baseline_periods.extend(baseline_result['_summaries'][baseline_result[private]]['_periods'])
            write_csv(temp/('baseline_'+label+'_events.csv'),baseline_events);write_csv(temp/('baseline_'+label+'_periods.csv'),baseline_periods)
        write_csv(temp/'manual_validation.csv',manual_stats or [{'status':manual_status,'case_id':case_id}]);write_csv(temp/'manual_event_matches.csv',manual_matches)
        write_csv(temp/'manual_period_checks.csv',manual_checks);write_csv(temp/'manual_metric_periods.csv',manual_periods);write_csv(temp/'manual_threshold_comparison.csv',manual_sweep or [{'status':'PENDING_MANUAL_ANNOTATION'}])
        write_csv(temp/'manual_common_period_metrics.csv',common_manual_stats or [{'status':'PENDING_MANUAL_ANNOTATION'}]);write_csv(temp/'manual_common_periods.csv',common_manual_periods)
        write_csv(temp/'manual_reference_validated.csv',manual,['case_id','phase','event_time_s','confidence','reviewer','notes','phase_extremum_error_s','input_format'])
        index=global_review_figures(temp,data,start,end,main,legacy,recover,case_id,args.manual_review_start,args.manual_review_end)
        for row in new_connections:
            matching=[r for r in index if r['t_start_s']<=row['start_s'] and row['end_s']<=r['t_end_s']]
            if not matching:raise RuntimeError('new connection has no complete review panel')
            panel=min(matching,key=lambda r:r['t_end_s']-r['t_start_s']);row['review_figure']=panel['figure_diagnostic']
        write_csv(temp/'new_connections_review.csv',new_connections);global_threshold_figure(temp,sweep);global_comparison_figure(temp,comparison)
        dump_json(temp/'regression_check.json',{'CPP_exact_current_legacy_vs_recover':clean_json(legacy['cpp'])==clean_json(recover['cpp']),
            'baseline_supplied':str(args.diagnostic_baseline) if args.diagnostic_baseline else None,'input_sha256_after':sha256(data['path']),
            'full_window_classified':all(np.isclose(v['classified_time_fraction'],1.) for v in main['methods'].values()),'manual_status':manual_status})
        (temp/'manual_annotation_instructions.md').write_text('# Manual reference\n\nFirst inspect blind panels. Fill manual_events_template.csv with case_id,phase,event_time_s,confidence,reviewer,notes; time in seconds; A maxima, B minima; high/uncertain; unique increasing times per phase. Do not copy automatic markers as truth. Aim for at least 31 high-confidence A events per case, including difficult and control intervals. Prefer complete-window annotation. Overlap panels share event IDs. Use --manual-events <filled.csv> with the same --tracking-method global-cycle, --threshold-profile selected_profile.json, input/dt/window and a NEW --out-dir; --skip-threshold-sweep can shorten a manual check. Global events remain provisional pending actual reference comparison. Plateau intervals do not define precise B times.\n')
        (temp/'report.md').write_text(global_report_text(summary,sweep))
        if sha256(data['path'])!=data['sha256']:raise RuntimeError('source flow changed during analysis')
        if out.exists():raise ValueError('global output appeared concurrently')
        summary['metadata']['elapsed_total_s']=time.perf_counter()-began
        summary['metadata']['elapsed_s_note']='elapsed_s measures analysis before export; elapsed_total_s includes plots and all exports'
        dump_json(temp/'summary.json',summary);dump_json(temp/'analysis_metadata.json',summary['metadata'])
        dump_json(temp/'run_manifest.json',{'files':{str(p.relative_to(temp)):sha256(p) for p in temp.rglob('*') if p.is_file()},'input_sha256':data['sha256'],'script_sha256':sha256(__file__)})
        temp.rename(out)
    print('global-cycle:',case_id,'CALIBRATION_PENDING_MANUAL;',manual_status)
    for phase,v in main['methods'].items():print(phase,'events',v['n_adopted_events'],'coverage',v['time_coverage'],'Jitter/RAP/PPQ5',[v['metrics'][k] for k in ('jitter_local_percent','rap_percent','ppq5_percent')])
    print('Saved:',out)
    return summary


def global_report_text(summary,sweep):
    lines=['# Full-window global-cycle inference','', 'Status: CALIBRATION_PENDING_MANUAL; '+summary['manual_validation_status']+'. All global metrics are PROVISIONAL. Local ACF predicts search regions only; event times are actual raw extrema with quadratic refinement. No period-time regularization or snapping.','',
        '## Method comparison','','| Phase | Method | Jitter % | RAP % | PPQ5 % | Tracked time | Jitter-used time |','|---|---|---:|---:|---:|---:|---:|']
    def f(v):return 'UNDEFINED' if v is None else f'{v:.6f}'
    for r in summary['comparison']:lines.append(f'| {r["phase"]} | {r["method"]} | {f(r["jitter_local_percent"])} | {f(r["rap_percent"])} | {f(r["ppq5_percent"])} | {r["tracked_time_fraction"]:.3%} | {r["jitter_local_percent_covered_time_fraction"]:.3%} |')
    lines+=['','Full-window classification is distinct from adopted-period coverage. time_classification.csv partitions every sample/time into adopted period, unresolved or explicit window edge; ambiguous flags can overlap provisional adopted periods. Metrics pool only within independent segments and use metric-specific eligible-segment means, identical to pooled_v2. All actual periods/terms are saved.','',
        '## Review and alternatives','','Inspect flow_full_window_events.png, event_panels paired blind/diagnostic images, unresolved_focus, new_connections_review.csv, event_selection_examples.png, acf_timeline.png and event_timeline.png. New edges are provisional, not manual truth. alternative_paths.csv and alternative_periods.csv show top paths and bounded forced-exclusion sensitivity. The recorded cap can leave additional exclusions unaudited; this is not a statistical confidence interval.','',
        '## Thresholds and sampling','','The same fixed score weights/profile are used across cases. threshold_profiles.json/csv and selected_profile.json preserve all settings. A full 7x5 local prominence/excursion grid plus one-factor ACF/range/confidence/noise sweeps is used unless --skip-threshold-sweep is explicit. Confidence-rho changes are confidence classification only, not hard selection gates. No profile is chosen by Jitter minimum. MAD of first differences can include physical slope: its floor is a sensitivity setting, not a measured instrument-noise calibration.','',
        'global_subsample_comparison.csv retains the same event IDs. plateau_representative_sensitivity.csv tests actual plateau endpoints/centers and may change complete boundary periods; B is not assumed precise. Without actual stride-1 data sampling comparison is NOT_AVAILABLE. CPP function/configuration are unchanged; frames and warnings are retained.','',
        'Manual annotations use the existing phase-specific one-to-one maximum-cardinality/minimum-time-error matcher at 50/100/150 us. manual_common_period_metrics.csv compares exactly the same matched consecutive periods, with missing/uncertain/automatic-segment boundaries preserved; manual_common_periods.csv records every pair. Full-annotation-extent metrics are retained separately. No annotation means PENDING_MANUAL_ANNOTATION, not precision/recall inferred from another detector. Score weights and amplitude gates may still bias selection; unresolved and ambiguous fractions must accompany each estimate. No physical stabilization mechanism is considered.','']
    return '\n'.join(lines)


CLOSURE_CONFIG = {
    'version':'closure_local_peak_downcross_v1', 'baseline_m3_s':0.,
    'stable_fraction_T0':.03, 'minimum_below_samples':2,
    'pulse_area_weight':.20, 'pulse_width_weight':.20, 'closing_evidence_weight':.25,
    'recovery_edge_penalty':.25, 'manual_tolerance_samples':[1,2],
    'closing_search_max_T0':.95, 'next_pulse_guard_min_T0':.55,
}
METRIC_NAMES = ('jitter_local_percent','rap_percent','ppq5_percent')


def read_minimum_area(path,dt,column='interior-min',time_column=0,unit='auto'):
    """Actual section areas: exclude boundary sections as selectSeparationByAreaRatio does."""
    path=Path(path).resolve();digest=sha256(path);a=np.loadtxt(path,comments='#',ndmin=2)
    header=' '.join(line.strip() for line in path.open() if line.lstrip().startswith('#'))
    if len(a)<3 or a.shape[1]<2 or not 0<=time_column<a.shape[1]:raise ValueError('invalid area matrix or step column')
    steps=a[:,time_column]
    if not np.isfinite(steps).all() or np.any(steps<0) or np.any(steps!=np.rint(steps)):raise ValueError('area steps must be finite nonnegative integers')
    diff=np.diff(steps)
    if np.any(diff<=0) or np.any(diff!=diff[0]):raise ValueError('duplicate, decreasing or variable-stride area steps; no interpolation')
    detected='mm2' if 'mm^2' in header or 'mm2' in header else 'm2' if 'm^2' in header or 'm2' in header else None
    if unit=='auto':
        if detected is None:raise ValueError('area unit unknown; supply --area-unit')
        unit=detected
    elif detected and unit!=detected:raise ValueError('explicit area unit disagrees with header')
    sections=[i for i in range(a.shape[1]) if i!=time_column]
    if column=='interior-min':
        if len(sections)<3:raise ValueError('interior-min requires at least three section columns; specify actual Amin column')
        selected=sections[1:-1]
    elif column=='all-min':selected=sections
    else:
        try:selected=[int(column)]
        except (ValueError,TypeError):raise ValueError('area-column must be interior-min, all-min or zero-based integer')
        if selected[0]==time_column or not 0<=selected[0]<a.shape[1]:raise ValueError('invalid area value column')
    v=a[:,selected]
    if np.isinf(v).any() or np.any(v[np.isfinite(v)]<0):raise ValueError('negative or infinite selected area; not clipped')
    finite=np.isfinite(v)
    if np.any(~finite.any(axis=1)):raise ValueError('area row has no finite selected section; not interpolated')
    x=np.min(np.where(finite,v,np.inf),axis=1)*(1e6 if unit=='m2' else 1.)
    fs=1/(diff[0]*dt)
    if sha256(path)!=digest:raise ValueError('area source changed while reading')
    return {'path':str(path),'sha256':digest,'steps':steps.astype(np.int64),'t':steps*dt,'x':x,'dt_s':dt,
            'stride':int(diff[0]),'sampling_interval_s':1/fs,'fs':fs,
            'area_metadata':{'header':header,'input_columns':a.shape[1],'time_column_zero_based':time_column,
                'selected_area_columns_zero_based':selected,'reduction':column,'input_unit':unit,'output_unit':'mm2',
                'undefined_selected_sections':int((~finite).sum()),'rows_with_undefined_selected_section':int((~finite).any(axis=1).sum()),
                'zero_samples':int((x==0).sum()),'negative_samples':0,'minimum_mm2':float(x.min()),'maximum_mm2':float(x.max()),
                'zero_interpretation':'stored lower limit; no inference of physical contact or unrecorded negative area'}}


def closure_archive_provenance(flow,area):
    provenance={'status':'EXTERNAL_INPUTS_NOT_ARCHIVE_VERIFIED','time_offset_s':None}
    p=Path(flow['path']);snapshot=p.parent.parent;record=snapshot.parent/'metrics.json'
    if snapshot.name=='snapshot' and record.is_file():
        saved=json.loads(record.read_text());hashes=saved.get('snapshot_sha256',{})
        for data in (flow,area):
            if data is None:continue
            path=Path(data['path'])
            if not path.is_relative_to(snapshot):raise ValueError('area/flow do not belong to the same archived snapshot')
            key=str(path.relative_to(snapshot))
            if hashes.get(key)!=data['sha256']:raise ValueError('snapshot hash mismatch for '+key)
        provenance.update(status='SAME_ARCHIVED_RUN_HASH_VERIFIED',archive_id=saved.get('archive_id'),archive_metrics_path=str(record),archive_metrics_sha256=sha256(record))
    if area:
        a,b=flow['t'],area['t'];common=np.intersect1d(flow['steps'],area['steps'])
        provenance.update(flow_stride=flow['stride'],area_stride=area['stride'],shared_step_count=len(common),
            first_sample_offset_s=float(b[0]-a[0]),last_sample_offset_s=float(b[-1]-a[-1]),time_offset_s=0.,
            time_note='both use saved step times dt; no resampling or lag correction')
    return provenance


def closure_crossing_candidates(data,peak_index,right_index,threshold_pct,period):
    """All real downward brackets, with stability checked on original samples."""
    t=data['t'];x=data['x'];i=int(peak_index);right=min(int(right_index),len(x)-1)
    peak=float(x[i]);threshold=peak*threshold_pct/100.;rows=[]
    if peak<=0 or right<=i:return rows
    hold=max(CLOSURE_CONFIG['minimum_below_samples'],int(np.ceil(CLOSURE_CONFIG['stable_fraction_T0']*(period or 0)*data['fs']))+1)
    for k in range(i,right):
        if not x[k]>threshold>=x[k+1]:continue
        at=float(t[k]+(x[k]-threshold)/(x[k]-x[k+1])*(t[k+1]-t[k]));stop=k+1+hold
        stability=stop<=len(x) and bool(np.all(x[k+1:stop]<=threshold))
        zero_stop=min(len(x),k+1+hold+2);zero=np.flatnonzero(x[k+1:zero_stop]==0)
        plateau=False
        if len(zero):
            z=k+1+int(zero[0]);plateau=z+1<len(x) and x[z+1]==0 and bool(np.all(x[k+1:z+2]<=threshold))
        rows.append({'sample_left':k,'sample_right':k+1,'step_left':int(data['steps'][k]),'step_right':int(data['steps'][k+1]),
            't_left_s':float(t[k]),'t_right_s':float(t[k+1]),'Ug_left_m3_s':float(x[k]),'Ug_right_m3_s':float(x[k+1]),
            'threshold_pct':threshold_pct,'threshold_m3_s':threshold,'cross_time_s':at,'slope_m3_s2':float((x[k+1]-x[k])/(t[k+1]-t[k])),
            'stability_required_samples':hold,'stability_required_s':CLOSURE_CONFIG['stable_fraction_T0']*(period or 0),
            'stable_below_threshold':stability,'immediate_zero_plateau':plateau,'stable':bool(stability or plateau),
            'stability_check_end_s':float(t[min(stop-1,len(t)-1)]),'selected':False})
    stable=[r for r in rows if r['stable']]
    if stable:stable[0]['selected']=True
    return rows


def closure_pulse_candidates(data,start,end,base,cfg,thresholds):
    """Independent common pulse DAG; global A events are candidates, never fixed seeds."""
    frames,_,raw=base['_prepared'];rows,windows=global_candidate_features(data,raw,'A',cfg);t=data['t'];x=data['x']
    all_crossings=[];by_peak={};strong=[r for r in rows if r['selectable']]
    for e in rows:
        i=e['sample_id'];T=e['local_period_s'];left=e['pulse_cell'][0];right=e['pulse_cell'][1]
        left=left if left is not None else max(0,int(np.searchsorted(t,e['raw_time_s']-1.1*(T or .01))))
        basin_right=right
        # A secondary hump can create an early valley without ending the physical pulse.
        # Continue inspecting real samples, but stop before the next distinct strong peak.
        next_peaks=[r['raw_time_s'] for r in strong if r['raw_time_s']-e['raw_time_s']>=CLOSURE_CONFIG['next_pulse_guard_min_T0']*(T or .01)]
        limit=min(e['raw_time_s']+CLOSURE_CONFIG['closing_search_max_T0']*(T or .01),min(next_peaks,default=float(t[-1])))
        right=max(i,min(len(x)-1,int(np.searchsorted(t,limit,side='left'))-1))
        valley=i+int(np.argmin(x[i:right+1]))
        e.update(pulse_id='P_'+e['event_id'],Upeak_m3_s=float(x[i]),Umin_m3_s=float(np.min(x[i:right+1])),
            pulse_start_s=float(t[left]),pulse_end_s=float(t[valley]),pulse_width_s=float(t[valley]-t[left]),
            closing_search_limit_s=float(t[right]),original_basin_end_s=float(t[basin_right]) if basin_right is not None else None,
            pulse_area_m3=float(np.trapz(np.maximum(x[left:valley+1],0),t[left:valley+1])),
            Upeak_definition='saved raw sample of selected main local maximum; no global amplitude normalization',
            in_window=bool(start-.03<=e['raw_time_s']<end+.03))
        peak_of_branch=i+int(np.argmax(x[i:valley+1]));e['largest_peak_before_valley_m3_s']=float(x[peak_of_branch]);e['largest_peak_before_valley_time_s']=float(t[peak_of_branch])
        branch={}
        for pct in thresholds:
            choices=closure_crossing_candidates(data,i,right,pct,T)
            for k,r in enumerate(choices):r.update(pulse_id=e['pulse_id'],peak_event_id=e['event_id'],crossing_candidate_id=e['pulse_id']+'_'+str(pct)+'_'+str(k))
            all_crossings.extend(choices);branch[pct]=choices
        by_peak[e['event_id']]=branch
        support=any(any(r['stable'] for r in choices) for choices in branch.values())
        e['any_stable_closure']=support
        parts={'node_pulse_area':CLOSURE_CONFIG['pulse_area_weight']*min(e['pulse_area_m3']/max(e['Upeak_m3_s']*(T or .01),np.finfo(float).tiny),1.),
               'node_pulse_width':CLOSURE_CONFIG['pulse_width_weight']*min(e['pulse_width_s']/(T or .01),1.),
               'node_closing_evidence':CLOSURE_CONFIG['closing_evidence_weight']*support}
        e.update(parts);e['node_score']+=sum(parts.values())
        e['alternative_peak_candidates']=[r['event_id'] for r in raw['A'] if e['pulse_start_s']<r['raw_time_s']<e['pulse_end_s'] and r['event_id']!=e['event_id']]
    ids,transitions,connected,audit,lookup=global_graph(rows,windows,frames,cfg);position={rows[j]['event_id']:n for n,j in enumerate(ids)}
    # Recover only a waveform-shape veto supported by real stable closures at both ends.
    # Time, same-basin, intermediate-pulse and double-period vetoes remain intact.
    for edge in audit:
        if edge['rejection_reasons']!=['INCOMPATIBLE_SHAPE']:continue
        a,b=position[edge['event_start']],position[edge['event_end']];ea,eb=rows[ids[a]],rows[ids[b]]
        jointly=[]
        for pct in thresholds:
            ca=next((r for r in by_peak[ea['event_id']][pct] if r['selected']),None);cb=next((r for r in by_peak[eb['event_id']][pct] if r['selected']),None)
            if ca and cb and cfg['search_min_period_ratio']<=(cb['cross_time_s']-ca['cross_time_s'])/edge['local_period_s']<=cfg['search_max_period_ratio']:jointly.append(pct)
        if not jointly or ea['secondary_lobe'] or eb['secondary_lobe']:continue
        edge.update(legal=True,recovery='DIRECT_STABLE_CLOSURE_EVIDENCE',recovery_supported_thresholds=jointly,
                    edge_score=edge['edge_phase_score']+edge['edge_acf_score']-CLOSURE_CONFIG['recovery_edge_penalty'])
        edge['warning_flags']=sorted(set(edge['warning_flags']+['DIRECT_CLOSURE_RECOVERY_PROVISIONAL']))
        transitions[a,b]=edge['edge_score'];connected[a,b]=True;lookup[(a,b)]=edge
    paths,marginal=global_dp([rows[j]['node_score'] for j in ids],transitions,connected,cfg);path=paths[0]
    track=global_path_track(rows,ids,path,lookup,cfg);chosen={e['event_id'] for e in track['events'] if e['accepted']}
    incident={}
    for k,kind in enumerate(path['kinds']):
        if kind=='CONNECT':
            edge=lookup[(path['nodes'][k],path['nodes'][k+1])]
            for identity in (edge['event_start'],edge['event_end']):incident.setdefault(identity,[]).append(edge)
    selected_nodes=set(path['nodes'])
    for node,j in enumerate(ids):
        e=track['events'][j];edges=incident.get(e['event_id'],[]);flags=list(e['shape_flags'])
        for edge in edges:flags.extend(edge['warning_flags'])
        alternatives=[k for k in range(len(ids)) if k not in selected_nodes and abs(rows[ids[k]]['raw_time_s']-e['raw_time_s'])<.55*(e['local_period_s'] or 0.) and np.isfinite(marginal[k])]
        margin=path['score']-max(marginal[k] for k in alternatives) if alternatives else None
        if margin is not None and margin<cfg['ambiguity_score_margin']:flags.append('AMBIGUOUS_PULSE_ASSIGNMENT')
        if e['largest_peak_before_valley_m3_s']>e['Upeak_m3_s']+32*np.finfo(float).eps*max(e['Upeak_m3_s'],np.finfo(float).tiny):flags.append('LARGER_ALTERNATIVE_PEAK_BEFORE_VALLEY')
        if 'PITCH_AMBIGUOUS' in e['local_period_flags']:flags.append('LOCAL_PERIOD_AMBIGUOUS')
        high=e['local_acf_confidence']>=cfg['high_acf_confidence'] and edges and not flags
        e.update(confidence_class='HIGH_CONFIDENCE' if high else 'ESTIMATED',selection_flags=sorted(set(flags)),
            path_score_margin=path['score']-paths[1]['score'] if len(paths)>1 else None,competing_peak_score_margin=margin,
            alternative_peak_times=[rows[ids[k]]['raw_time_s'] for k in alternatives],max_marginal_score=float(marginal[node]) if np.isfinite(marginal[node]) else None)
    for e in track['events']:
        e.setdefault('confidence_class','ESTIMATED');e.setdefault('selection_flags',e['shape_flags'])
    selected=sorted((e for e in track['events'] if e['accepted']),key=lambda e:e['raw_time_s'])
    for order,e in enumerate(selected):e['pulse_order']=order
    for r in all_crossings:r['pulse_selected']=r['peak_event_id'] in chosen
    return track,by_peak,all_crossings,audit


def closure_threshold_tracks(data,pulses,branches,thresholds,frames,start,end,cfg):
    selected=sorted((e for e in pulses['events'] if e['accepted']),key=lambda e:e['raw_time_s']);all_rows=[];tracks={}
    for pct in thresholds:
        name='Ug_downcross_'+format(pct,'g').replace('.','p')+'pct';events=[];segments=[];current=[];previous=None
        for p in selected:
            options=branches[p['event_id']][pct];chosen=next((r for r in options if r['selected']),None)
            flags=list(p['selection_flags']);stable_count=sum(r['stable'] for r in options)
            if len(options)>1:flags.append('ALTERNATIVE_CROSSINGS')
            if chosen is None:flags.append('NO_THRESHOLD_CROSSING' if not options else 'NO_STABLE_THRESHOLD_CROSSING')
            at=chosen['cross_time_s'] if chosen else None;inside=at is not None and start<=at<end
            e=dict(p,event_id=p['pulse_id']+'_'+str(pct)+'pct',method=name,threshold_pct=pct,threshold_m3_s=p['Upeak_m3_s']*pct/100,
                peak_event_id=p['event_id'],peak_time_s=p['raw_time_s'],subsample_time_s=at,raw_time_s=chosen['t_right_s'] if chosen else None,
                accepted=bool(inside),event_found=chosen is not None,n_crossings=len(options),n_stable_crossings=stable_count,
                confidence_class='HIGH_CONFIDENCE' if p['confidence_class']=='HIGH_CONFIDENCE' and not len(options)>1 else 'ESTIMATED',
                selection_flags=sorted(set(flags)),reason='STABLE_DESCENDING_CROSSING' if inside else 'OUTSIDE_WINDOW' if chosen else flags[-1],
                pulse_segment_id=p['segment_id'],segment_id=None,crossing=chosen)
            if chosen:e.update({k:v for k,v in chosen.items() if k not in ('selected','pulse_selected')})
            connected=False;reason='NO_PREVIOUS_EVENT'
            if inside and previous is not None:
                T,conf,supported,periodflags=global_period_at((previous['subsample_time_s']+at)/2,frames)
                ratio=(at-previous['subsample_time_s'])/T if T else None
                connected=(p['pulse_order']==previous['pulse_order']+1 and p['segment_id']==previous['pulse_segment_id'] and supported and ratio is not None and cfg['search_min_period_ratio']<=ratio<=cfg['search_max_period_ratio'])
                reason='CONNECTED_REAL_CROSSINGS' if connected else 'MISSING_PULSE_OR_PARENT_SEGMENT_OR_INTERVAL_RANGE'
                e['interval_over_local_period']=ratio
            e['previous_connection_reason']=reason
            if not inside or not connected:
                if current:segments.append({'segment_id':len(segments)+1,'event_indices':current});current=[]
            index=len(events);events.append(e)
            if inside:current.append(index);previous=e
            else:previous=None
            all_rows.append(dict(e,selected_pulse=True))
        if current:segments.append({'segment_id':len(segments)+1,'event_indices':current})
        for s in segments:
            for j in s['event_indices']:events[j]['segment_id']=s['segment_id']
        tracks[name]={'method':name,'events':events,'segments':segments}
    # Rejected real peak candidates are retained for all three thresholds, never silently discarded.
    for p in pulses['events']:
        if p['accepted'] or not p['in_window']:continue
        for pct in thresholds:
            choices=branches[p['event_id']][pct];hit=next((r for r in choices if r['selected']),None)
            all_rows.append(dict(p,method='Ug_downcross_'+format(pct,'g').replace('.','p')+'pct',threshold_pct=pct,
                threshold_m3_s=p['Upeak_m3_s']*pct/100,selected_pulse=False,accepted=False,event_found=hit is not None,
                cross_time_s=hit['cross_time_s'] if hit else None,n_crossings=len(choices),reason='PULSE_NOT_ON_COMMON_DAG_PATH',crossing=hit))
    return tracks,all_rows


def closure_filter_track(track,start,end,confirmed=False,representative=None):
    result=copy.deepcopy(track);segments=[]
    for s in result['segments']:
        block=[]
        for j in s['event_indices']:
            e=result['events'][j]
            if representative=='plateau_center_s':
                if e.get('plateau_width_s',0)>0:e['subsample_time_s']=.5*(e['plateau_start_s']+e['plateau_end_s'])
            elif representative in ('plateau_start_s','plateau_end_s'):
                if e.get('plateau_width_s',0)>0:e['subsample_time_s']=e[representative]
            elif representative is not None:e['subsample_time_s']=e[representative]
            valid=e['subsample_time_s'] is not None and start<=e['subsample_time_s']<end and (not confirmed or e.get('confidence_class')=='HIGH_CONFIDENCE')
            if valid:block.append(j)
            elif block:segments.append({'segment_id':len(segments)+1,'event_indices':block});block=[]
        if block:segments.append({'segment_id':len(segments)+1,'event_indices':block})
    result['segments']=segments;adopted={j for s in segments for j in s['event_indices']}
    for j,e in enumerate(result['events']):e['accepted']=j in adopted
    return result


def closure_statistics(track,start,end,confirmed=False,representative=None,period_status=None):
    tr=closure_filter_track(track,start,end,confirmed,representative);membership={s['segment_id']:s['event_indices'] for s in tr['segments']}
    metrics=pooled_metrics_v2(tr,membership,start,end);periods=[];events=[e for e in tr['events'] if e['accepted']];duration=end-start
    for s in tr['segments']:
        es=[tr['events'][j] for j in s['event_indices']]
        for k,(a,b) in enumerate(zip(es,es[1:])):
            periods.append({'method':tr['method'],'segment_id':s['segment_id'],'period_id':k,'event_start':a['event_id'],'event_end':b['event_id'],
                'start_s':a['subsample_time_s'],'end_s':b['subsample_time_s'],'period_s':b['subsample_time_s']-a['subsample_time_s'],
                'quality':'HIGH_CONFIDENCE' if a.get('confidence_class')==b.get('confidence_class')=='HIGH_CONFIDENCE' else 'ESTIMATED',
                'flags':sorted(set(a.get('selection_flags',[])+b.get('selection_flags',[])))})
    ambiguous=period_status=='PERIOD_AMBIGUOUS'
    if ambiguous:
        for name in METRIC_NAMES:metrics['metric_details'][name]['conditional_value']=metrics[name];metrics[name]=None;metrics['metric_details'][name]['value']=None
    row={'method':tr['method'],'scope':'confirmed_model_only' if confirmed else 'including_provisional','start_s':start,'end_s':end,
         'status':'PERIOD_AMBIGUOUS' if ambiguous else 'PROVISIONAL' if metrics['jitter_local_percent'] is not None else 'UNDEFINED_JITTER',
         'n_events':len(events),'n_periods':len(periods),'n_segments':len(tr['segments']),
         'n_provisional_events':sum(e.get('confidence_class')!='HIGH_CONFIDENCE' for e in events),
         'mean_period_s':float(np.mean([p['period_s'] for p in periods])) if periods else None,
         'tracked_time_fraction':sum(p['period_s'] for p in periods)/duration,'manual_status':'PENDING_MANUAL_ANNOTATION'}
    for name in METRIC_NAMES:
        detail=metrics['metric_details'][name];row[name]=metrics[name]
        for key in ('n_cycles','n_terms','mean_period_s','covered_time_fraction','sum_numerator_terms_s'):row[name+'_'+key]=detail[key]
    return {'track':tr,'metrics':metrics,'periods':periods,'row':row}


def closure_classify(stat,data,start,end,frames):
    periods=stat['periods'];out=[];cursor=start
    def gap(a,b,edge=False):
        T=global_period_at((a+b)/2,frames)[0];kind='WINDOW_EDGE' if edge and T and b-a<=1.65*T else 'UNRESOLVED'
        return {'start_s':a,'end_s':b,'classification':kind,'reason':'explicit boundary fragment' if kind=='WINDOW_EDGE' else 'no consecutive observed events; missing crossing, rejected candidate, parent gap or incompatible interval'}
    for p in periods:
        if p['start_s']<cursor-1e-12:raise RuntimeError('overlapping periods')
        if p['start_s']>cursor:out.append(gap(cursor,p['start_s'],cursor==start))
        out.append(dict(p,classification='VALID_EVENT_PERIOD' if p['quality']=='HIGH_CONFIDENCE' else 'PROVISIONAL_EVENT_PERIOD',reason='observed event pair'))
        cursor=p['end_s']
    if cursor<end:out.append(gap(cursor,end,bool(periods)))
    samples=data['t'][(data['t']>=start)&(data['t']<end)];count=np.zeros(len(samples),int)
    for r in out:
        r.update(method=stat['row']['method'],duration_s=r['end_s']-r['start_s'],window_fraction=(r['end_s']-r['start_s'])/(end-start))
        count+=(samples>=r['start_s'])&(samples<r['end_s'])
    if not np.all(count==1) or not np.isclose(sum(r['duration_s'] for r in out),end-start):raise RuntimeError('closure/area classification gap or overlap')
    stat['row'].update(classification_coverage=1.,unresolved_fraction=sum(r['window_fraction'] for r in out if r['classification']=='UNRESOLVED'),
        provisional_time_fraction=sum(r['window_fraction'] for r in out if r['classification']=='PROVISIONAL_EVENT_PERIOD'),
        ambiguous_time_fraction=sum(r['window_fraction'] for r in out if any('AMBIGUOUS' in f for f in r.get('flags',[]))),
        n_unresolved=sum(r['classification']=='UNRESOLVED' for r in out))
    return out


def closure_flow_area_matches(stats,frames):
    rows=[];differences=[]
    flow={k:v for k,v in stats.items() if k.startswith('Ug_downcross_')}
    for a,sa in flow.items():
        ea={e['pulse_id']:e for e in sa['track']['events'] if e['accepted']}
        for b,sb in flow.items():
            if a>=b:continue
            eb={e['pulse_id']:e for e in sb['track']['events'] if e['accepted']}
            for identity in sorted(ea.keys()|eb.keys()):
                x,y=ea.get(identity),eb.get(identity)
                differences.append({'pulse_id':identity,'method_a':a,'method_b':b,'found_a':x is not None,'found_b':y is not None,
                    'time_a_s':x['subsample_time_s'] if x else None,'time_b_s':y['subsample_time_s'] if y else None,
                    'delta_b_minus_a_s':y['subsample_time_s']-x['subsample_time_s'] if x and y else None,
                    'same_Upeak':bool(x and y and x['Upeak_m3_s']==y['Upeak_m3_s'])})
    if 'Amin_minima' not in stats:return rows,differences
    minima=[e for e in stats['Amin_minima']['track']['events'] if e['accepted']]
    T=float(np.median([f['period_s'] for f in frames if f['period_s']])) if any(f['period_s'] for f in frames) else None
    for name,s in flow.items():
        es=[e for e in s['track']['events'] if e['accepted']];pairs=optimal_event_assignment([e['subsample_time_s'] for e in es],[e['subsample_time_s'] for e in minima],.45*T) if T else []
        paired={i:(j,error) for i,j,error in pairs};used={j for _,j,_ in pairs}
        for i,e in enumerate(es):
            j,error=paired.get(i,(None,None));area=minima[j] if j is not None else None
            rows.append({'method':name,'pulse_id':e['pulse_id'],'flow_event_id':e['event_id'],'flow_time_s':e['subsample_time_s'],
                'area_event_id':area['event_id'] if area else None,'area_time_s':area['subsample_time_s'] if area else None,
                'area_minus_flow_s':area['subsample_time_s']-e['subsample_time_s'] if area else None,
                'area_plateau_width_s':area['plateau_width_s'] if area else None,
                'status':'PROVISIONAL_PHASE_CORRESPONDENCE' if area else 'UNMATCHED_FLOW',
                'matching_tolerance_s':.45*T if T else None,'not_manual_truth':True})
        for j,e in enumerate(minima):
            if j not in used:rows.append({'method':name,'pulse_id':None,'flow_event_id':None,'flow_time_s':None,'area_event_id':e['event_id'],'area_time_s':e['subsample_time_s'],'area_minus_flow_s':None,'status':'UNMATCHED_AREA','not_manual_truth':True})
    return rows,differences


def closure_manual_reference(path,case_id,data,area,stats,start,end):
    if path is None:return [],[],[],[],[]
    with Path(path).open(newline='') as f:reader=csv.DictReader(f);fields=set(reader.fieldnames or []);original=list(reader)
    rows=[]
    if 'signal' not in fields:
        old=read_manual_reference(path,case_id,data,start,end)
        for r in old:rows.append(dict(r,method='Ug_global_maxima_A' if r['phase']=='A' else 'Ug_global_minima_B'))
    else:
        required={'case_id','signal','event_type','threshold_pct','event_time_s','confidence','reviewer','notes'}
        if not required.issubset(fields):raise ValueError('manual CSV missing required closure/area fields')
        for r in original:
            if r['case_id']!=case_id:continue
            sig=r['signal'];kind=r['event_type'];pct=r['threshold_pct'].strip()
            if sig=='Ug' and kind=='downcross':method='Ug_downcross_'+format(float(pct),'g').replace('.','p')+'pct'
            elif sig=='Ug' and kind in ('maxima','minima'):method='Ug_global_maxima_A' if kind=='maxima' else 'Ug_global_minima_B'
            elif sig=='Amin' and kind in ('maxima','minima'):method='Amin_'+kind
            else:raise ValueError('unknown manual signal/event_type')
            if method not in stats:raise ValueError('manual method unavailable in current analysis: '+method)
            at=float(r['event_time_s'])
            if not np.isfinite(at) or not start<=at<end or r['confidence'] not in ('high','uncertain'):raise ValueError('invalid manual time or confidence')
            rows.append(dict(r,event_time_s=at,method=method))
    validation=[];matches=[];common=[];common_periods=[]
    for method in sorted({r['method'] for r in rows}):
        refs=[r for r in rows if r['method']==method];times=[r['event_time_s'] for r in refs]
        if any(b<=a for a,b in zip(times,times[1:])):raise ValueError('manual times must be unique and increasing within each event definition')
        refs=[dict(r,phase='A') for r in refs];tracks={'A':stats[method]['track'],'B':{'events':[],'segments':[]}}
        sr,mt,_,_=evaluate_manual_reference(refs,tracks,start,end)
        cs,cp=manual_common_period_comparison(refs,tracks,mt,start,end)
        # Existing phase matcher retains 50/100/150 us; actual source-stride 1/2-sample comparison below.
        fs=(area if method.startswith('Amin') else data)['fs']
        for scope in ('high_only','all_annotations'):
            rr=[r for r in refs if scope=='all_annotations' or r['confidence']=='high']
            if not rr:continue
            lo,hi=rr[0]['event_time_s'],rr[-1]['event_time_s']
            for samples in CLOSURE_CONFIG['manual_tolerance_samples']:
                tol=samples/fs;uncertain=[r['event_time_s'] for r in refs if scope=='high_only' and r['confidence']=='uncertain']
                es=[e for e in tracks['A']['events'] if e['accepted'] and lo-tol<=e['subsample_time_s']<=hi+tol and not any(abs(e['subsample_time_s']-at)<=tol for at in uncertain)]
                pairs=optimal_event_assignment([r['event_time_s'] for r in rr],[e['subsample_time_s'] for e in es],tol);tp=len(pairs);fp=len(es)-tp;fn=len(rr)-tp
                validation.append({'method':method,'scope':scope,'tolerance_samples':samples,'tolerance_us':tol*1e6,'TP':tp,'FP':fp,'FN':fn,
                    'precision':tp/len(es) if es else None,'recall':tp/len(rr),'F1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
                    'status':'MANUAL_REFERENCE_EVALUATED','annotation_extent_start_s':lo,'annotation_extent_end_s':hi,
                    'median_absolute_error_s':float(np.median([abs(p[2]) for p in pairs])) if pairs else None})
        matches.extend(dict(r,method=method) for r in mt);common.extend(dict(r,method=method) for r in cs);common_periods.extend(dict(r,method=method) for r in cp)
    return rows,validation,matches,common,common_periods


def closure_comparison_figure(out,rows):
    groups=sorted({r['case_id'] for r in rows});names=list(dict.fromkeys(r['method'] for r in rows if r['scope']=='including_provisional' and r.get('window_name','full')=='full'))
    labels=[n.replace('Ug_downcross_','Ug ').replace('pct','%').replace('0p5','0.5').replace('Ug_global_','Ug ').replace('_','\n') for n in names]
    fig,axes=plt.subplots(len(groups),4,figsize=(18,4.5*len(groups)),squeeze=False)
    for row,case in enumerate(groups):
        by={r['method']:r for r in rows if r['case_id']==case and r['scope']=='including_provisional' and r.get('window_name','full')=='full'}
        for col,key in enumerate((*METRIC_NAMES,'jitter_local_percent_covered_time_fraction')):
            values=[by.get(n,{}).get(key) for n in names];ax=axes[row,col]
            ax.bar(range(len(names)),[v if v is not None else np.nan for v in values],color=['tab:red','tab:orange','tab:purple','tab:blue','tab:green','gray','silver'][:len(names)])
            ax.set(xticks=range(len(names)),xticklabels=labels,title=case+' '+{'jitter_local_percent':'Jitter [%]','rap_percent':'RAP [%]','ppq5_percent':'PPQ5 [%]','jitter_local_percent_covered_time_fraction':'Jitter-used time fraction'}[key]);ax.tick_params(axis='x',labelsize=7)
            if col==3:ax.set_ylim(0,1.05)
            ax.grid(axis='y',alpha=.15)
    fig.suptitle('Observed events, all eligible segments; estimates provisional; no gap bridging');fig.tight_layout();fig.savefig(out/'jitter_comparison.png',dpi=160);plt.close(fig)


def closure_figures(out,data,area,stats,intervals,pulses,crossings,matches,frames,case,start,end):
    colors=['tab:red','tab:orange','tab:purple'];closure_names=[n for n in stats if n.startswith('Ug_downcross_')];t=data['t'];x=data['x']
    def flow_plot(ax,left,right,blind=False,diagnostic=False):
        mask=(t>=left)&(t<right);ax.plot(t[mask],x[mask]*1e6,c='black',lw=.85,label='Ug raw')
        if not blind:
            ps=[e for e in pulses['events'] if e['accepted'] and left<=e['raw_time_s']<right]
            ax.scatter([e['raw_time_s'] for e in ps],[e['Upeak_m3_s']*1e6 for e in ps],s=19,c='tab:blue',label='common pulse Upeak')
            for n,color,marker in zip(closure_names,colors,('o','x','v')):
                es=[e for e in stats[n]['track']['events'] if e['accepted'] and left<=e['subsample_time_s']<right]
                ax.scatter([e['subsample_time_s'] for e in es],[e['threshold_m3_s']*1e6 for e in es],c=color,marker=marker,s=25,label=n)
                for r in intervals[n]:
                    if r['end_s']>left and r['start_s']<right and r['classification'] in ('UNRESOLVED','PROVISIONAL_EVENT_PERIOD'):ax.axvspan(max(left,r['start_s']),min(right,r['end_s']),color=color if r['classification']=='UNRESOLVED' else 'tab:orange',alpha=.05 if r['classification']=='UNRESOLVED' else .02)
                if diagnostic:
                    for e in es:ax.annotate(f'{e["pulse_id"]} {e["threshold_pct"]:g}%',(e['subsample_time_s'],e['threshold_m3_s']*1e6),xytext=(3,14+closure_names.index(n)*14),textcoords='offset points',fontsize=6,rotation=35,arrowprops={'arrowstyle':'->','lw':.5,'color':color})
            ax.legend(fontsize=7,ncol=2,loc='upper right')
        ax.set(xlim=(left,right),ylabel='Ug [mL/s]',xlabel='Time [s]');ax.grid(alpha=.1)
    def area_plot(ax,left,right,blind=False,diagnostic=False):
        if area is None:ax.text(.5,.5,'MISSING_INPUT: Amin',transform=ax.transAxes,ha='center');ax.set_xlim(left,right);return
        tt,aa=area['t'],area['x'];mask=(tt>=left)&(tt<right);ax.plot(tt[mask],aa[mask],c='black',lw=.85,label='Amin raw')
        if not blind:
            for name,color,marker in [('Amin_maxima','tab:blue','o'),('Amin_minima','tab:green','v')]:
                es=[e for e in stats[name]['track']['events'] if e['accepted'] and left<=e['subsample_time_s']<right]
                ax.scatter([e['subsample_time_s'] for e in es],np.interp([e['subsample_time_s'] for e in es],tt,aa),c=color,marker=marker,s=20,label=name)
                for e in es:
                    if e['plateau_width_s']>0:
                        for at in (e['plateau_start_s'],e['plateau_end_s']):ax.plot(at,np.interp(at,tt,aa),marker='|',c=color,ms=10)
                    if diagnostic:ax.annotate(e['event_id']+' '+e['confidence_class'],(e['subsample_time_s'],np.interp(e['subsample_time_s'],tt,aa)),xytext=(2,8),textcoords='offset points',fontsize=6,rotation=30)
                for r in intervals[name]:
                    if r['end_s']>left and r['start_s']<right and r['classification']=='UNRESOLVED':ax.axvspan(max(left,r['start_s']),min(right,r['end_s']),color=color,alpha=.08)
            ax.legend(fontsize=7)
        ax.set(xlim=(left,right),ylabel='Amin [mm²]',xlabel='Time [s]');ax.grid(alpha=.1)
    for filename,fn in [('flow_threshold_crossings_full.png',flow_plot),('min_area_extrema_full.png',area_plot)]:
        fig,ax=plt.subplots(figsize=(16,4.5));fn(ax,start,end);fig.suptitle(case+' '+filename[:-4]);fig.tight_layout();fig.savefig(out/filename,dpi=170);plt.close(fig)
    fig,axes=plt.subplots(len(closure_names),1,figsize=(15,2*len(closure_names)),squeeze=False)
    for ax,n in zip(axes[:,0],closure_names):
        for r in intervals[n]:ax.axvspan(r['start_s'],r['end_s'],color={'VALID_EVENT_PERIOD':'tab:green','PROVISIONAL_EVENT_PERIOD':'tab:orange','UNRESOLVED':'tab:red','WINDOW_EDGE':'gray'}[r['classification']],alpha=.65)
        ax.set(xlim=(start,end),yticks=[],ylabel=n)
    axes[-1,0].set_xlabel('Time [s]');fig.suptitle(case+' classified whole window: green=valid model, orange=provisional, red=unresolved, gray=edge');fig.tight_layout();fig.savefig(out/'threshold_event_timeline.png',dpi=160);plt.close(fig)
    focus=[(.574,.600),(.607,.621),(.790,.806)] if case=='r0p25' else [(.549,.565),(.681,.697),(.749,.767)] if case=='r0p30' else []
    T=np.median([f['period_s'] for f in frames if f['period_s']]) if any(f['period_s'] for f in frames) else .008
    for folder in ('flow_threshold_crossings_zoom','min_area_extrema_zoom','flow_area_aligned_panels','manual_review','crossing_interpolation_zoom'):(out/folder).mkdir()
    index=[];windows=[];cursor=start;i=0
    while cursor<end-1e-12:
        stop=min(cursor+.060,end)
        if stop-cursor<.05 and i:cursor=max(start,end-.060)
        i+=1;windows.append((f'panel_{i:03d}',cursor,stop))
        if stop>=end:break
        cursor=stop-min(T,.015)
    for k,(a,b) in enumerate(focus,1):
        left,right=max(start,a-5*T),min(end,b+5*T)
        for folder,fn in [('flow_threshold_crossings_zoom',flow_plot),('min_area_extrema_zoom',area_plot)]:
            fig,ax=plt.subplots(figsize=(16,5));fn(ax,left,right,diagnostic=True);fig.suptitle(case+f' focus {k}: actual events, no T0 snapping');fig.tight_layout();fig.savefig(out/folder/f'focus_{k:03d}.png',dpi=170);plt.close(fig)
    from matplotlib.patches import ConnectionPatch
    for name,left,right in windows:
        for blind in (True,False):
            fig,axes=plt.subplots(2,1,figsize=(15,8),sharex=True);flow_plot(axes[0],left,right,blind,not blind);area_plot(axes[1],left,right,blind,not blind)
            fig.suptitle(case+' '+name+(' BLIND: annotate actual downcrossings and area extrema' if blind else ' DIAGNOSTIC: estimates require manual validation'))
            fig.tight_layout();file='manual_review/'+name+('_blind.png' if blind else '_diagnostic.png');fig.savefig(out/file,dpi=155);plt.close(fig)
            if not blind:index.append({'case_id':case,'panel_id':name,'start_s':left,'end_s':right,'blind':'manual_review/'+name+'_blind.png','diagnostic':file})
        fig,axes=plt.subplots(2,1,figsize=(15,7),sharex=True);flow_plot(axes[0],left,right);area_plot(axes[1],left,right)
        if area:
            for r in matches:
                if r['method']==closure_names[0] and r.get('area_time_s') is not None and r.get('flow_time_s') is not None and left<=r['flow_time_s']<right and left<=r['area_time_s']<right:
                    fig.add_artist(ConnectionPatch((r['flow_time_s'],np.interp(r['flow_time_s'],t,x)*1e6),(r['area_time_s'],np.interp(r['area_time_s'],area['t'],area['x'])),coordsA='data',coordsB='data',axesA=axes[0],axesB=axes[1],color='tab:purple',alpha=.35,lw=.6))
        fig.suptitle(case+' provisional same-cycle correspondence; no lag correction');fig.tight_layout();fig.savefig(out/'flow_area_aligned_panels'/f'{name}.png',dpi=155);plt.close(fig)
    chosen=[r for r in crossings if r['selected'] and r['pulse_selected'] and start<=r['cross_time_s']<end]
    # One bracket per threshold per difficult interval, plus a deterministic control.
    audit_choices=[]
    for pct in sorted({r['threshold_pct'] for r in chosen}):
        options=[r for r in chosen if r['threshold_pct']==pct]
        for a,b in focus+[(start,end)]:
            r=min(options,key=lambda r:abs(r['cross_time_s']-(a+b)/2)) if options else None
            if r and r['crossing_candidate_id'] not in {z['crossing_candidate_id'] for z in audit_choices}:audit_choices.append(r)
    for k,r in enumerate(audit_choices,1):
        i=r['sample_left'];a=max(0,i-4);b=min(len(t),i+6);fig,ax=plt.subplots(figsize=(9,4))
        ax.plot(t[a:b],x[a:b]*1e6,'ko-',ms=4,label='saved samples');ax.axhline(r['threshold_m3_s']*1e6,c='tab:red',ls='--',label=f'{r["threshold_pct"]:g}% local Upeak')
        ax.scatter([r['cross_time_s']],[r['threshold_m3_s']*1e6],marker='x',s=65,c='tab:purple',label='linear downcross');ax.scatter(t[i:i+2],x[i:i+2]*1e6,s=70,facecolors='none',edgecolors='tab:blue',label='bracket pair')
        ax.set(xlabel='Time [s]',ylabel='Ug [mL/s]',title=case+' '+r['crossing_candidate_id']);ax.ticklabel_format(axis='x',useOffset=False);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'crossing_interpolation_zoom'/f'crossing_{k:03d}.png',dpi=170);plt.close(fig)
    write_csv(out/'crossing_interpolation_review.csv',audit_choices);write_csv(out/'manual_review_index.csv',index)
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    pcts=[stats[n]['track']['events'][0]['threshold_pct'] if stats[n]['track']['events'] else float(n.split('_')[-1][:-3].replace('p','.')) for n in closure_names]
    axes[0].plot(pcts,[stats[n]['row']['n_events'] for n in closure_names],'o-');axes[0].set_ylabel('Observed events')
    for metric in METRIC_NAMES:axes[1].plot(pcts,[stats[n]['row'][metric] for n in closure_names],'o-',label=metric)
    axes[1].legend(fontsize=7);axes[1].set_ylabel('Metric [%]')
    for n in closure_names:
        es=[e for e in stats[n]['track']['events'] if e['accepted']];axes[2].scatter([e['subsample_time_s'] for e in es],[(e['subsample_time_s']-e['peak_time_s'])*1e3 for e in es],s=7,label=n)
    axes[2].set(xlabel='Time [s]',ylabel='Closure time minus Upeak [ms]');axes[2].legend(fontsize=7)
    for ax in axes[:2]:ax.set_xlabel('Local peak threshold [%]')
    fig.suptitle(case+' closure threshold sensitivity; no profile chosen by minimum Jitter');fig.tight_layout();fig.savefig(out/'closure_threshold_sensitivity.png',dpi=160);plt.close(fig)


def closure_validate_statistics(stats,data,area,threshold_rows,intervals,start,end):
    """Independent scalar-loop audit of all real brackets, periods and numerator terms."""
    checked=0
    for r in threshold_rows:
        if not r.get('accepted') or not r.get('selected_pulse'):continue
        k=r['sample_left'];threshold=r['Upeak_m3_s']*r['threshold_pct']/100.;x=data['x'];t=data['t']
        if not x[k]>threshold>=x[k+1]:raise RuntimeError('invalid downcross bracket')
        expected=t[k]+(threshold-x[k])*(t[k+1]-t[k])/(x[k+1]-x[k])
        if not math.isclose(expected,r['subsample_time_s'],abs_tol=1e-13):raise RuntimeError('linear interpolation mismatch')
        if r['sample_right']!=k+1 or r['slope_m3_s2']>=0:raise RuntimeError('non-descending crossing')
        checked+=1
    for name,stat in stats.items():
        tr=stat['track'];by={e['event_id']:e for e in tr['events']};period_by={}
        for p in stat['periods']:
            a,b=by[p['event_start']],by[p['event_end']]
            if p['period_s']<=0 or not math.isclose(b['subsample_time_s']-a['subsample_time_s'],p['period_s'],abs_tol=1e-13):raise RuntimeError('invalid exported period')
            if name.startswith('Ug_downcross_') and (b['pulse_order']!=a['pulse_order']+1 or b['pulse_segment_id']!=a['pulse_segment_id']):raise RuntimeError('missing pulse or parent boundary bridged')
            period_by.setdefault(p['segment_id'],[]).append(p)
        for width,metric in ((2,METRIC_NAMES[0]),(3,METRIC_NAMES[1]),(5,METRIC_NAMES[2])):
            terms=[];values=[];valid_pairs=set()
            for sid,periods in period_by.items():
                if len(periods)<width:continue
                pp=[r['period_s'] for r in periods];values.extend(pp)
                for j in range(len(pp)-width+1):
                    value=abs(pp[j+1]-pp[j]) if width==2 else abs(pp[j+width//2]-sum(pp[j:j+width])/width)
                    terms.append(value);valid_pairs.add((sid,periods[j]['event_start'],periods[j+width-1]['event_end']))
            detail=stat['metrics']['metric_details'][metric];mean=sum(values)/len(values) if values else None
            if detail['n_terms']!=len(terms) or detail['n_cycles']!=len(values):raise RuntimeError('metric count mismatch')
            expected=100*sum(terms)/len(terms)/mean if terms else None
            actual=stat['metrics'][metric]
            if actual is not None and (expected is None or not math.isclose(expected,actual,rel_tol=1e-10,abs_tol=1e-10)):raise RuntimeError('metric formula mismatch')
            if mean and not math.isclose(mean,detail['mean_period_s'],rel_tol=1e-12):raise RuntimeError('metric denominator mismatch')
            for term in stat['metrics']['_terms']:
                if term['metric']==metric and (term['segment_id'],term['event_start'],term['event_end']) not in valid_pairs:raise RuntimeError('term crossed a segment boundary')
            if not math.isclose(sum(values)/(end-start),detail['covered_time_fraction'],abs_tol=1e-12):raise RuntimeError('metric coverage mismatch')
        ts=(area if name.startswith('Amin') else data)['t'];ts=ts[(ts>=start)&(ts<end)];membership=np.zeros(len(ts),int)
        for row in intervals[name]:membership+=(ts>=row['start_s'])&(ts<row['end_s'])
        if not np.all(membership==1):raise RuntimeError('incomplete whole-window classification')
    return {'status':'PASS','real_crossing_brackets_checked':checked,'method_count':len(stats),'all_metric_terms_and_denominators_recomputed':True,
        'no_parent_gap_or_missing_pulse_bridges':True,'all_samples_classified_exactly_once':True,'not_manual_ground_truth':True}


def closure_area_cli(args):
    if args.input is None:raise ValueError('--input is required')
    dt=resolve_dt(args.dt,args.manifest);data=read_flow(args.input,dt)
    start=data['t'][0] if args.start is None else args.start;end=data['t'][-1]+1/data['fs'] if args.end is None else args.end
    if not np.isfinite([start,end]).all() or not data['t'][0]<=start<end<=data['t'][-1]+1/data['fs']+1e-12:raise ValueError('invalid closure analysis window')
    try:thresholds=sorted(float(v.strip()) for v in args.closure_thresholds.split(','))
    except ValueError:raise ValueError('--closure-thresholds requires comma-separated percentages')
    if not thresholds or len(set(thresholds))!=len(thresholds) or any(not np.isfinite(v) or not 0<v<100 for v in thresholds):raise ValueError('unique closure thresholds must lie between 0 and 100 percent')
    if args.cpp_frame_length!=CPP_CONFIG['frame_sec'] or args.cpp_hop!=CPP_CONFIG['hop_sec']:raise ValueError('closure analysis preserves CPP settings')
    if args.primary_prominence!=10 or args.prominences!=TRACKING_CONFIG['prominence_sweep_percent']:raise ValueError('use --threshold-profile for global settings; legacy prominence overrides not used')
    cfg=resolve_global_profile(args.threshold_profile,args);began=time.perf_counter()
    case=args.case_id or next((s.split('_P')[0] for s in Path(data['path']).parts if s.startswith('r') and '_P' in s),Path(data['path']).stem)
    out=(args.out_dir or Path('analysis_runs')/('closure_area_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve()
    area_path=args.area_input.resolve() if args.area_input else Path(data['path']).with_name('area.dat')
    protected=[Path(data['path']),area_path]+[p.resolve() for p in (args.manual_events,args.manifest,args.diagnostic_baseline) if p]
    if out.exists():raise ValueError('closure output already exists; choose a new --out-dir, including with --overwrite')
    if any(p==out or p.is_relative_to(out) for p in protected):raise ValueError('output contains protected source input')
    area=read_minimum_area(area_path,dt,args.area_column,args.area_time_column,args.area_unit) if area_path.is_file() else None
    if area and (start<area['t'][0] or end>area['t'][-1]+1/area['fs']+1e-12):raise ValueError('area does not cover requested full window')
    provenance=closure_archive_provenance(data,area)
    global_flow=global_cycle_track(data,start,end,cfg)
    # Existing baseline implementations are called without modification.
    track_cfg=dict(TRACKING_CONFIG,pitch_min_hz=args.fmin,pitch_max_hz=args.fmax)
    old=analyze_arrays(data['t'],data['x'],data['fs'],start,end,track_cfg,CPP_CONFIG,args.min_cycles)
    recover=analyze_arrays(data['t'],data['x'],data['fs'],start,end,dict(track_cfg,reanchor_policy='recover-rejected'),CPP_CONFIG,args.min_cycles)
    pulses,branches,crossings,pulse_edges=closure_pulse_candidates(data,start,end,global_flow,cfg,thresholds)
    tracks,threshold_rows=closure_threshold_tracks(data,pulses,branches,thresholds,global_flow['frames'],start,end,cfg)
    period_status=global_flow['methods']['A']['period_identification_status']
    statuses={name:period_status for name in tracks};global_area=global_cycle_track(area,start,end,cfg) if area else None
    if global_area:
        for phase,name in (('A','Amin_maxima'),('B','Amin_minima')):
            tracks[name]=copy.deepcopy(global_area['methods'][phase]['track']);tracks[name]['method']=name;statuses[name]=global_area['methods'][phase]['period_identification_status']
    for phase,name in (('A','Ug_global_maxima_A'),('B','Ug_global_minima_B')):
        tracks[name]=copy.deepcopy(global_flow['methods'][phase]['track']);tracks[name]['method']=name;statuses[name]=global_flow['methods'][phase]['period_identification_status']
    stats={name:closure_statistics(track,start,end,period_status=statuses[name]) for name,track in tracks.items()}
    intervals={name:closure_classify(stat,area if name.startswith('Amin') else data,start,end,global_area['frames'] if name.startswith('Amin') else global_flow['frames']) for name,stat in stats.items()}
    for name,stat in stats.items():
        for event in stat['track']['events']:event['method']=name
    rows=[];all_periods=[];all_terms=[]
    for name,stat in stats.items():
        for confirmed in (False,True):
            for window_name,left,right in [('full',start,end),('early_half',start,(start+end)/2),('late_half',(start+end)/2,end)]:
                item=stat if not confirmed and window_name=='full' else closure_statistics(tracks[name],left,right,confirmed,period_status=statuses[name])
                if item is not stat:closure_classify(item,area if name.startswith('Amin') else data,left,right,global_area['frames'] if name.startswith('Amin') else global_flow['frames'])
                row=dict(item['row'],case_id=case,window_name=window_name,flow_CPP_dB=old['cpp']['flow_CPP_dB']);rows.append(row)
                all_periods.extend(dict(p,case_id=case,scope=row['scope'],window_name=window_name) for p in item['periods'])
                all_terms.extend(dict(p,method=name,case_id=case,scope=row['scope'],window_name=window_name) for p in item['metrics']['_terms'])
    matching,differences=closure_flow_area_matches(stats,global_flow['frames'])
    annotation,manual_stats,manual_matches,manual_common,manual_common_periods=closure_manual_reference(args.manual_events,case,data,area,stats,start,end)
    manual_status='MANUAL_REFERENCE_EVALUATED' if annotation else 'PENDING_MANUAL_ANNOTATION'
    for row in rows:row['manual_status']='MANUAL_REFERENCE_EVALUATED' if any(r['method']==row['method'] for r in annotation) else 'PENDING_MANUAL_ANNOTATION'
    sensitivities=[];sensitivity_events=[];representatives=[]
    if global_area:
        cached={cfg['acf_window_sec']:global_area['_prepared']}
        for profile in global_threshold_profiles(cfg,not args.skip_threshold_sweep):
            pcfg=profile['config'];window=pcfg['acf_window_sec']
            if window not in cached:
                ff,cc=global_local_periods(area,start,end,pcfg);cached[window]=(ff,cc,global_raw_candidates(area,start,end,ff,pcfg))
            result=global_area if pcfg==cfg else global_cycle_track(area,start,end,pcfg,cached[window],False)
            for phase,name in (('A','Amin_maxima'),('B','Amin_minima')):
                entry=result['methods'][phase];s=closure_statistics(dict(entry['track'],method=name),start,end,period_status=entry['period_identification_status'])
                sensitivities.append(dict(s['row'],profile_id=profile['profile_id'],config=pcfg))
                sensitivity_events.extend({'profile_id':profile['profile_id'],'method':name,'event_id':e['event_id'],'time_s':e['subsample_time_s']} for e in s['track']['events'] if e['accepted'])
        for name in ('Amin_maxima','Amin_minima'):
            keys=('subsample_time_s','raw_time_s') if name=='Amin_maxima' else ('subsample_time_s','raw_time_s','plateau_start_s','plateau_center_s','plateau_end_s')
            for key in keys:
                s=closure_statistics(tracks[name],start,end,representative=key,period_status=statuses[name]);representatives.append(dict(s['row'],representative=key,same_event_ids=True))
    baseline_check={'cpp_legacy_vs_recover_exact':clean_json(old['cpp'])==clean_json(recover['cpp']),'saved_baseline_supplied':bool(args.diagnostic_baseline)}
    if args.diagnostic_baseline:
        saved=json.loads(args.diagnostic_baseline.read_text())
        if saved['metadata']['input_sha256']!=data['sha256'] or saved['metadata']['start_s']!=start or saved['metadata']['end_s']!=end:raise ValueError('saved global baseline input/window mismatch')
        for phase in ('A','B'):
            before=saved['global_cycle_estimated'][phase];after=global_flow['methods'][phase]
            for key in ('metrics','n_adopted_events','time_coverage','unresolved_time_fraction','period_identification_status','path_score','path_score_margin'):
                if clean_json(before[key])!=clean_json(after[key]):raise RuntimeError('existing global baseline changed: '+phase+' '+key)
        if clean_json(saved['cpp'])!=clean_json(old['cpp']):raise RuntimeError('saved CPP frames/settings changed')
        baseline_check['saved_global_metrics_counts_scores_and_CPP_exact']=True
    validation=closure_validate_statistics(stats,data,area,threshold_rows,intervals,start,end)
    summary={'metadata':{'case_id':case,'start_s':start,'end_s':end,'window_semantics':'[start,end)','dt_s':dt,'flow_input_path':data['path'],
        'flow_input_sha256':data['sha256'],'flow_stride':data['stride'],'flow_fs_hz':data['fs'],'area_input_path':area['path'] if area else str(area_path),
        'area_input_sha256':area['sha256'] if area else None,'area_stride':area['stride'] if area else None,'area_fs_hz':area['fs'] if area else None,
        'area_status':'AVAILABLE' if area else 'MISSING_INPUT','area_definition':area['area_metadata'] if area else None,
        'script_sha256':sha256(__file__),'provenance':provenance,'created_utc':datetime.now(timezone.utc).isoformat()},
        'closure_config':CLOSURE_CONFIG,'thresholds_pct':thresholds,'global_config':cfg,'manual_status':manual_status,
        'manual_validation':manual_stats,'global_calibration':'CALIBRATION_PENDING_MANUAL','comparison':rows,
        'methods':{name:{'statistics':stat['row'],'metrics':stat['metrics']} for name,stat in stats.items()},
        'flow_CPP_reference':old['cpp'],'CPP_config':CPP_CONFIG,'CPP_definition':CPP_DEFINITION,
        'baseline_check':baseline_check,'validation':validation,'area_sensitivity_profiles':len(sensitivities)//2,
        'limitations':['closure-like flow marker is not physical contact','model-confirmed is not human ground truth',
            'shared pulse selection may remain ambiguous; thresholds do not define pulse truth','missing crossings split all metric terms; no fabricated events',
            'area min-plateau representative sensitivity is separate from closure thresholds']}
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.closure_area_',dir=out.parent) as tmp:
        temp=Path(tmp)
        dump_json(temp/'summary.json',summary);dump_json(temp/'analysis_metadata.json',summary['metadata']);dump_json(temp/'selected_profile.json',{'config':cfg})
        write_csv(temp/'flow_threshold_events.csv',threshold_rows);write_csv(temp/'flow_crossing_candidates.csv',crossings)
        write_csv(temp/'flow_pulse_candidates.csv',pulses['events']);write_csv(temp/'flow_pulse_edge_audit.csv',pulse_edges)
        write_csv(temp/'jitter_metrics_by_event_definition.csv',rows);write_csv(temp/'periods_by_definition.csv',all_periods);write_csv(temp/'metric_terms_by_definition.csv',all_terms)
        write_csv(temp/'qc_intervals.csv',[dict(r,case_id=case) for group in intervals.values() for r in group])
        write_csv(temp/'event_matching_flow_vs_area.csv',matching);write_csv(temp/'threshold_crossing_differences.csv',differences)
        write_csv(temp/'area_threshold_sensitivity.csv',sensitivities);write_csv(temp/'area_sensitivity_events.csv',sensitivity_events);write_csv(temp/'area_representative_sensitivity.csv',representatives)
        write_csv(temp/'flow_local_period_frames.csv',global_flow['frames']);write_csv(temp/'flow_local_period_candidates.csv',global_flow['period_candidates'])
        if global_area:
            es=[dict(e,signal='Amin',phase=phase) for phase,entry in global_area['methods'].items() for e in entry['track']['events']]
            write_csv(temp/'min_area_candidates.csv',es);write_csv(temp/'min_area_events_A_B.csv',[e for e in es if e['accepted']]);write_csv(temp/'min_area_plateaus.csv',[e for e in es if e['method']=='B' and e['plateau_width_s']>0],list(es[0]) if es else None)
            write_csv(temp/'area_local_period_frames.csv',global_area['frames']);write_csv(temp/'area_local_period_candidates.csv',global_area['period_candidates']);write_csv(temp/'area_edge_score_audit.csv',global_area['edge_audit'])
            mask=(area['t']>=start)&(area['t']<end);write_csv(temp/'minimum_area_waveform.csv',({'step':int(step),'time_s':at,'Amin_mm2':value} for step,at,value in zip(area['steps'][mask],area['t'][mask],area['x'][mask])))
        write_csv(temp/'events_by_definition.csv',[dict(e,definition=name) for name,s in stats.items() for e in s['track']['events'] if e['accepted']])
        write_csv(temp/'manual_events_template.csv',[],['case_id','signal','event_type','threshold_pct','event_time_s','confidence','reviewer','notes'])
        write_csv(temp/'manual_reference_validated.csv',annotation);write_csv(temp/'manual_validation.csv',manual_stats or [{'status':manual_status}]);write_csv(temp/'manual_matches.csv',manual_matches)
        write_csv(temp/'manual_common_period_metrics.csv',manual_common);write_csv(temp/'manual_common_periods.csv',manual_common_periods)
        write_csv(temp/'cpp_frames_reference.csv',old['cpp']['rows'])
        dump_json(temp/'baseline_legacy.json',old);dump_json(temp/'baseline_recover_rejected.json',recover)
        closure_figures(temp,data,area,stats,intervals,pulses,crossings,matching,global_flow['frames'],case,start,end);closure_comparison_figure(temp,rows)
        (temp/'manual_annotation_instructions.md').write_text('# Manual annotation\n\nInspect the overlapping blind Ug/Amin panels first. Fill manual_events_template.csv with case_id,signal,event_type,threshold_pct,event_time_s,confidence,reviewer,notes. signal=Ug/Amin; event_type=downcross/maxima/minima; threshold_pct=0.5/2/5 only for Ug downcross; otherwise blank. Time in seconds, confidence=high/uncertain, increasing unique times within each definition. Peaks use the same local Upeak and zero baseline for all three thresholds. Do not annotate a trough as an absent downcross. Plateau minima have an interval uncertainty. Re-run the same command with --manual-events filled.csv and a NEW --out-dir; --skip-threshold-sweep is available. Existing phase=A/B manual files remain supported for baseline Ug extrema. Precision/recall at one/two native samples and common TP-matched period metrics are separate from automatic QC. No annotations means PENDING_MANUAL_ANNOTATION.\n')
        (temp/'report.md').write_text(closure_report_text(summary,representatives,matching))
        for source in (data,area):
            if source and sha256(source['path'])!=source['sha256']:raise RuntimeError('source changed during closure/area analysis')
        validation['exported_period_rows']=len(all_periods);validation['exported_metric_terms']=len(all_terms)
        dump_json(temp/'validation.json',validation)
        summary['metadata']['elapsed_total_s']=time.perf_counter()-began;dump_json(temp/'summary.json',summary);dump_json(temp/'analysis_metadata.json',summary['metadata'])
        dump_json(temp/'run_manifest.json',{'files':{str(p.relative_to(temp)):sha256(p) for p in temp.rglob('*') if p.is_file()},'script_sha256':sha256(__file__),'flow_input_sha256':data['sha256'],'area_input_sha256':area['sha256'] if area else None})
        if out.exists():raise ValueError('output appeared concurrently')
        temp.rename(out)
    print('closure-threshold-and-area:',case,manual_status)
    for name,s in stats.items():print(name,'J/R/P',[s['row'][k] for k in METRIC_NAMES],'periods',s['row']['n_periods'],'J-covered',s['row']['jitter_local_percent_covered_time_fraction'])
    print('Saved:',out)
    return summary


def closure_report_text(summary,representatives,matches):
    lines=['# Flow downcrossing and minimum-area event comparison','',
        'Status: PROVISIONAL; '+summary['manual_status']+'. Whole window searched, all eligible segments pooled; no fabricated crossing or gap-bridging term. Model-confirmed events are not manual truth.','',
        '## Full-window metrics','','| Definition | Jitter % | RAP % | PPQ5 % | Events | Periods | Jitter-used time | Unresolved |','|---|---:|---:|---:|---:|---:|---:|---:|']
    def f(v):return 'UNDEFINED' if v is None else f'{v:.5f}'
    for r in summary['comparison']:
        if r['scope']=='including_provisional' and r['window_name']=='full':lines.append(f'| {r["method"]} | {f(r[METRIC_NAMES[0]])} | {f(r[METRIC_NAMES[1]])} | {f(r[METRIC_NAMES[2]])} | {r["n_events"]} | {r["n_periods"]} | {r["jitter_local_percent_covered_time_fraction"]:.2%} | {r["unresolved_fraction"]:.2%} |')
    lines+=['','![Comparison](jitter_comparison.png)','','The statistics table includes early/late halves and confirmed_model_only separately. Numerators, eligible period means, actual period IDs and metric-specific coverage are exported. Every native sample is classified once as VALID_EVENT_PERIOD / PROVISIONAL_EVENT_PERIOD / UNRESOLVED / WINDOW_EDGE; classification coverage is not Jitter coverage.','',
        '## Event definition and limitations','','The common pulse DAG starts from all real raw maxima and valley basins, not fixed global-A seeds. Fixed pulse area/width and observed closing evidence supplement existing balanced-v1 node scores. Shape-incompatible edges can be recovered only with actual stable same-threshold crossings at both ends, no secondary lobes, and all existing time/basin/intermediate-pulse vetoes intact. Recovery is provisional. Pulse choices and alternatives, scores and recovery edges are exported. Three thresholds use the same selected Upeak sample and pulse ID, zero flow baseline, linear interpolation on real adjacent downward samples, and a below-threshold hold of at least 0.03T0 / two sample points or an immediate true zero plateau. Multiple real crossings are retained; the first stable crossing is selected and flagged. No crossing means missing, not minimum/zero/ACF replacement. The physically correct pulse phase is still subject to manual review.','',
        '## Area input and representation','',json.dumps(summary['metadata']['area_definition'],ensure_ascii=False),
        '', 'area.dat contains section areas; interior-min uses finite interior section columns, excludes two boundary cuts, and follows solver separation selection / saved archive convention. Undefined cuts are not zero. Original area steps and unit are retained (mm² output); no resampling. Area maxima and minima independently use unchanged global-cycle balanced-v1; area sensitivity uses the same shared profile grid. Min-plateau start/center/end and maxima interpolation ON/OFF keep event identities fixed. A minimum is not automatically physical contact.','',
        '## Review files','','Inspect flow_threshold_crossings_full.png, min_area_extrema_full.png, threshold_event_timeline.png, closure_threshold_sensitivity.png, difficult zooms, aligned panels, and crossing_interpolation_zoom. manual_review contains complete overlapping blind/diagnostic panels. Annotation format and reanalysis instructions are in manual_annotation_instructions.md. Flow-versus-area correspondence is a one-to-one time-proximity diagnostic within 0.45 median local T0, with unmatched events retained; it is neither enforced simultaneity nor a causal delay estimate.','',
        'CPP reference: '+f(summary['flow_CPP_reference']['flow_CPP_dB'])+' dB. CPP code/configuration and old event trackers are unchanged. Area CPP is not computed.','',
        'Area representative sensitivity:','']
    for r in representatives:lines.append(f'- {r["method"]} / {r["representative"]}: J={f(r[METRIC_NAMES[0]])}%, RAP={f(r[METRIC_NAMES[1]])}%, PPQ5={f(r[METRIC_NAMES[2]])}%, coverage={r["jitter_local_percent_covered_time_fraction"]:.2%}')
    return '\n'.join(lines)+'\n'


ALL_EVENT_METHODS = ('Ug_max','Ug_min','Ug_zero','Ug_3pct','Ug_10pct','Amin_max','Amin_min')
ALL_EVENT_RULE = 'raw local reversals; prominence=0; no distance, ACF, shape, pulse or interval gate'


def all_event_extrema(data,prefix):
    """Every saved local reversal; one flat extremum is one center-time event."""
    output={}
    for kind,polarity in (('max',1),('min',-1)):
        ids,props=signal.find_peaks(polarity*data['x'],prominence=0,plateau_size=True);rows=[]
        for k,i in enumerate(ids):
            left,right=int(props['left_edges'][k]),int(props['right_edges'][k]);offset,height,status=quadratic(polarity*data['x'],int(i))
            at=.5*(data['t'][left]+data['t'][right]) if right>left else data['t'][i]+offset/data['fs']
            rows.append({'method':prefix+'_'+kind,'event_id':prefix+'_'+kind+'_'+str(i),'event_time_s':float(at),
                'raw_time_s':float(data['t'][i]),'sample_id':int(i),'step':int(data['steps'][i]),'sample_value':float(data['x'][i]),
                'prominence':float(props['prominences'][k]),'plateau_left_index':left,'plateau_right_index':right,
                'plateau_start_s':float(data['t'][left]),'plateau_end_s':float(data['t'][right]),
                'plateau_width_s':float(data['t'][right]-data['t'][left]),'timing_uncertainty_s':float((right-left+1)/(2*data['fs'])),
                'interpolation':'PLATEAU_CENTER' if right>left else status,'all_event_detector_rule':ALL_EVENT_RULE})
        output[kind]=rows
    return output


def all_event_extract(flow,area,start,end):
    extrema=all_event_extrema(flow,'Ug');all_series={'Ug_max':extrema['max'],'Ug_min':extrema['min']}
    x,t=flow['x'],flow['t'];nonpositive=x<=0;starts=np.flatnonzero(nonpositive & np.r_[False,x[:-1]>0]);zero=[]
    for j in starts:
        right=int(j)
        while right+1<len(x) and x[right+1]<=0:right+=1
        k=int(j)-1;alpha=float(x[k]/(x[k]-x[j]));at=float(t[k]+alpha*(t[j]-t[k]))
        zero.append({'method':'Ug_zero','event_id':'Ug_zero_'+str(j),'event_time_s':at,'raw_time_s':float(t[j]),
            'sample_id':int(j),'step':int(flow['steps'][j]),'sample_value':0.,'sample_left':k,'sample_right':int(j),
            'value_left':float(x[k]),'value_right':float(x[j]),'linear_fraction':alpha,
            'plateau_start_s':float(t[j]),'plateau_end_s':float(t[right]),'plateau_width_s':float(t[right]-t[j]),
            'region_last_sample':right,'interpolation':'LINEAR_ZERO_CROSSING' if x[j]<0 else 'EXACT_FIRST_ZERO_SAMPLE',
            'zero_rule':'positive to nonpositive; one event per continuous nonpositive run; epsilon=0'})
    all_series['Ug_zero']=zero;branch_rows=[];candidate_crossings=[]
    minima=extrema['min'];min_ids=np.array([e['sample_id'] for e in minima],int)
    for pct in (3.,10.):
        name='Ug_'+format(pct,'g')+'pct';found=[];used_brackets=set()
        for peak in extrema['max']:
            i=peak['sample_id'];position=np.searchsorted(min_ids,peak['plateau_right_index'],side='right')
            valley=minima[position] if position<len(minima) else None;right=valley['plateau_left_index'] if valley else len(x)-1
            if t[i]>=end or t[right]<start:continue
            # Reuse the real-bracket/linear-interpolation implementation, ignoring all hold gates.
            # Its legacy stable/selected fields never decide adoption in all-events mode.
            options=closure_crossing_candidates(flow,i,right,pct,None);hit=options[-1] if options else None
            if hit:
                bracket=(hit['sample_left'],hit['sample_right'])
                if bracket in used_brackets:raise RuntimeError('same threshold bracket assigned to two local peaks')
                used_brackets.add(bracket)
                alpha=(hit['cross_time_s']-hit['t_left_s'])/(hit['t_right_s']-hit['t_left_s'])
                e={'method':name,'event_id':name+'_peak_'+str(i),'event_time_s':hit['cross_time_s'],
                    'raw_time_s':hit['t_right_s'],'sample_id':hit['sample_right'],'step':hit['step_right'],'sample_value':hit['threshold_m3_s'],
                    'sample_left':hit['sample_left'],'sample_right':hit['sample_right'],'value_left':hit['Ug_left_m3_s'],
                    'value_right':hit['Ug_right_m3_s'],'linear_fraction':alpha,'peak_id':peak['event_id'],'peak_time_s':float(t[i]),
                    'peak_sample_id':i,'Upeak_m3_s':float(x[i]),'threshold_pct':pct,'threshold_m3_s':hit['threshold_m3_s'],
                    'branch_end_sample':right,'branch_end_s':float(t[right]),'interpolation':'LINEAR_LAST_DOWNCROSS',
                    'crossing_count':len(options),'selection_rule':'last real downward bracket before next local valley; no stability gate'}
                found.append(e)
            for rank,option in enumerate(options):
                candidate_crossings.append({'method':name,'peak_id':peak['event_id'],'peak_sample_id':i,'peak_time_s':float(t[i]),
                    'Upeak_m3_s':float(x[i]),'threshold_pct':pct,'threshold_m3_s':option['threshold_m3_s'],
                    'branch_end_sample':right,'branch_end_s':float(t[right]),'sample_left':option['sample_left'],'sample_right':option['sample_right'],
                    't_left_s':option['t_left_s'],'t_right_s':option['t_right_s'],'value_left':option['Ug_left_m3_s'],'value_right':option['Ug_right_m3_s'],
                    'event_time_s':option['cross_time_s'],'crossing_rank':rank+1,'selected_last_crossing':rank==len(options)-1,
                    'inside_window':bool(start<=option['cross_time_s']<end),'stability_gate_applied':False})
            branch_rows.append({'method':name,'peak_id':peak['event_id'],'peak_time_s':float(t[i]),'Upeak_m3_s':float(x[i]),
                'threshold_m3_s':float(x[i])*pct/100,'branch_end_sample':right,'branch_end_s':float(t[right]),
                'valley_id':valley['event_id'] if valley else None,'right_censored_branch':valley is None,'crossing_count':len(options),
                'first_crossing_s':options[0]['cross_time_s'] if options else None,'last_crossing_s':hit['cross_time_s'] if hit else None,
                'event_found':bool(hit),'event_inside_window':bool(hit and start<=hit['cross_time_s']<end),
                'reason':'LAST_DESCENDING_CROSSING' if hit else 'NONPOSITIVE_PEAK' if x[i]<=0 else 'NO_THRESHOLD_CROSSING'})
        all_series[name]=found
    if area:
        aa=all_event_extrema(area,'Amin');all_series.update(Amin_max=aa['max'],Amin_min=aa['min'])
    series={};counts=[]
    for name in ALL_EVENT_METHODS:
        candidates=[e for e in all_series.get(name,[]) if start<=e['event_time_s']<end];candidates.sort(key=lambda e:(e['event_time_s'],e['event_id']))
        unique=[]
        for e in candidates:
            if not unique or e['event_time_s']!=unique[-1]['event_time_s']:unique.append(e)
        series[name]=unique;counts.append({'method':name,'detected_in_window':len(candidates),'unique_events_used':len(unique),
            'duplicate_time_count':len(candidates)-len(unique),'noise_excluded':0,'ACF_excluded':0,'shape_excluded':0,'length_excluded':0,
            'all_detected_events_used':'YES','status':'MISSING_AREA_INPUT' if name.startswith('Amin') and area is None else 'AUTOMATED_NOT_MANUALLY_VERIFIED'})
    return series,branch_rows,candidate_crossings,counts


def all_event_metrics(events,method,start,end):
    """One sorted event vector, all N intervals, one mean denominator; never segmented."""
    times=np.array([e['event_time_s'] for e in events]);v=np.diff(times)
    if np.any(v<=0):raise ValueError('all-event times must be strictly increasing')
    m=period_metrics(v);N=len(v);mean=m['mean_period_s'];median=float(np.median(v)) if N else None
    row=dict(m,method=method,estimator='J_all_detected',n_events=len(events),n_intervals=N,n_jitter_terms=max(0,N-1),n_RAP_terms=max(0,N-2),n_PPQ5_terms=max(0,N-4),
        first_event_s=float(times[0]) if len(times) else None,last_event_s=float(times[-1]) if len(times) else None,
        min_interval_s=float(v.min()) if N else None,max_interval_s=float(v.max()) if N else None,std_interval_s=float(np.std(v)) if N else None,
        mean_frequency_hz=1/mean if mean else None,max_over_median_interval=float(v.max()/median) if N else None,
        event_span_fraction=float((times[-1]-times[0])/(end-start)) if N else 0.,scan_fraction=1.,
        start_s=start,end_s=end,all_detected_events_used='YES',intervals_dropped=0,intervals_split=0,segments_used=1 if len(events) else 0,
        manual_status='MANUAL_VALIDATION_PENDING',status='AUTOMATED_NOT_MANUALLY_VERIFIED' if N>=2 else 'NOT_DEFINED_INSUFFICIENT_EVENTS')
    periods=[];terms=[]
    for i,value in enumerate(v):periods.append({'method':method,'interval_id':i,'event_start':events[i]['event_id'],'event_end':events[i+1]['event_id'],
        'start_s':float(times[i]),'end_s':float(times[i+1]),'interval_s':float(value),'over_median_interval':float(value/median),
        'long_interval_diagnostic_only':bool(value>2*median),'included_in_all_event_sequence':True})
    for width,key in ((2,'jitter_local_percent'),(3,'rap_percent'),(5,'ppq5_percent')):
        for i in range(max(0,N-width+1)):
            numerator=abs(v[i+1]-v[i]) if width==2 else abs(v[i+width//2]-sum(v[i:i+width])/width)
            terms.append({'method':method,'metric':key,'term_id':i,'first_interval_id':i,'last_interval_id':i+width-1,
                'event_start':events[i]['event_id'],'event_end':events[i+width]['event_id'],'start_s':float(times[i]),'end_s':float(times[i+width]),
                'time_s':float(times[i+1]) if width==2 else float(times[i+width//2]),'numerator_s':float(numerator),
                'mean_all_intervals_s':mean,'period_pair_jitter_contribution_percent':float(100*numerator/mean) if width==2 else None,
                'normalized_term_percent':float(100*numerator/mean),'included':True})
    js=[r['normalized_term_percent'] for r in terms if r['metric']==METRIC_NAMES[0]]
    if js and not math.isclose(sum(js)/len(js),m['jitter_local_percent'],rel_tol=1e-12,abs_tol=1e-12):raise RuntimeError('J_i average differs from all-event Jitter')
    return row,periods,terms


def all_event_figures(out,data,area,series,metrics,contributions,start,end,case):
    t,x=data['t'],data['x'];colors={'Ug_max':'tab:blue','Ug_min':'tab:green','Ug_zero':'black','Ug_3pct':'tab:orange','Ug_10pct':'tab:purple','Amin_max':'tab:blue','Amin_min':'tab:green'}
    def overlay(ax,names,signal_data,left,right):
        if signal_data is None:ax.text(.5,.5,'MISSING_AREA_INPUT',transform=ax.transAxes,ha='center');return
        tt,xx=signal_data['t'],signal_data['x'];scale=1 if names[0].startswith('Amin') else 1e6;mask=(tt>=left)&(tt<right)
        ax.plot(tt[mask],xx[mask]*scale,c='black',lw=.7,label='raw signal')
        for k,name in enumerate(names):
            es=[e for e in series[name] if left<=e['event_time_s']<right];et=[e['event_time_s'] for e in es]
            ax.scatter(et,np.interp(et,tt,xx)*scale,s=10 if name.endswith(('max','min')) else 23,marker=('o','v','x')[k%3],color=colors[name],label=name,zorder=3)
        ax.set(xlim=(left,right),xlabel='Time [s]',ylabel='Amin [mm²]' if scale==1 else 'Ug [mL/s]');ax.legend(fontsize=7,ncol=3);ax.grid(alpha=.12)
    for stem,names,dd in [('flow_extrema_full',('Ug_max','Ug_min'),data),('flow_closure_full',('Ug_zero','Ug_3pct','Ug_10pct'),data),('area_extrema_full',('Amin_max','Amin_min'),area)]:
        fig,ax=plt.subplots(figsize=(16,4));overlay(ax,names,dd,start,end);ax.set_title(case+' all independently detected events; not manually verified');fig.tight_layout();fig.savefig(out/(stem+'.png'),dpi=155);plt.close(fig)
    fig,axes=plt.subplots(7,1,figsize=(15,14),sharex=True)
    for ax,name in zip(axes,ALL_EVENT_METHODS):
        rr=[r for r in contributions if r['method']==name];ax.plot([r['time_s'] for r in rr],[r['period_pair_jitter_contribution_percent'] for r in rr],'.-',ms=2,lw=.5,color=colors[name]);ax.set_ylabel(name+'\nJ_i [%]',fontsize=8);ax.grid(alpha=.12);ax.set_xlim(start,end)
    axes[-1].set_xlabel('Time [s]');fig.suptitle(case+' J_i contributions, mean over all adjacent interval pairs = J_all_detected');fig.tight_layout();fig.savefig(out/'jitter_contributions.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(7,1,figsize=(15,17));examples=[]
    for ax,name in zip(axes,ALL_EVENT_METHODS):
        rr=[r for r in contributions if r['method']==name]
        if not rr:ax.text(.5,.5,name+': insufficient events',transform=ax.transAxes,ha='center');continue
        top=max(rr,key=lambda r:r['period_pair_jitter_contribution_percent']);pad=max(.001,.15*(top['end_s']-top['start_s']));left=max(start,top['start_s']-pad);right=min(end,top['end_s']+pad)
        overlay(ax,(name,),area if name.startswith('Amin') else data,left,right);ax.set_title(f'{name}: largest J_i={top["period_pair_jitter_contribution_percent"]:.3f}% at {top["time_s"]:.6f}s; INCLUDED',fontsize=9);examples.append(top)
    fig.suptitle(case+' largest-contribution examples, retained in main results');fig.tight_layout();fig.savefig(out/'largest_jitter_examples.png',dpi=145);plt.close(fig)
    write_csv(out/'largest_jitter_examples.csv',examples)
    # Three important contexts, not a new stack of overlapping review panels.
    fig,axes=plt.subplots(3,1,figsize=(15,10));zero=series['Ug_zero']
    if len(zero)>1:
        k=int(np.argmax(np.diff([e['event_time_s'] for e in zero])));a,b=zero[k]['event_time_s'],zero[k+1]['event_time_s'];overlay(axes[0],('Ug_zero','Ug_3pct','Ug_10pct'),data,max(start,a-.003),min(end,b+.003));axes[0].set_title(f'Longest zero-to-zero gap: {(b-a)*1000:.3f} ms; used as ONE interval')
    lo,hi=(.574,.600) if case=='r0p25' else (.675,.700)
    if not start<=lo<hi<=end:lo,hi=start+.25*(end-start),start+.40*(end-start)
    overlay(axes[1],('Ug_max','Ug_min'),data,lo,hi);axes[1].set_title('Multiple observed reversals are not merged to a common cycle')
    lo,hi=(.627,.644) if case=='r0p25' else (.837,.854)
    if not start<=lo<hi<=end:lo,hi=start+.6*(end-start),start+.75*(end-start)
    overlay(axes[2],('Ug_zero','Ug_3pct','Ug_10pct'),data,lo,hi);axes[2].set_title('Independent thresholds / zero entry: no artificial events at missing crossings')
    fig.suptitle(case+' selected diagnostic examples; automated, not manually validated');fig.tight_layout();fig.savefig(out/'important_event_examples.png',dpi=155);plt.close(fig)


def all_events_cli(args):
    if args.input is None:raise ValueError('--input is required')
    data=read_flow(args.input,resolve_dt(args.dt,args.manifest));start=.5 if args.start is None else args.start;end=1. if args.end is None else args.end
    if not np.isfinite([start,end]).all() or not data['t'][0]<=start<end<=data['t'][-1]+1/data['fs']+1e-12:raise ValueError('all-events window outside input support')
    if sorted(float(v) for v in args.closure_thresholds.split(',')) not in ([.5,2.,5.],[3.,10.]):raise ValueError('all-events uses independent fixed 3%,10% thresholds')
    if args.manual_events:raise ValueError('all-events does not gate on manual references; use exported events for independent review')
    case=args.case_id or next((s.split('_P')[0] for s in Path(data['path']).parts if s.startswith('r') and '_P' in s),Path(data['path']).stem)
    out=(args.out_dir or Path('analysis_runs')/('all_events_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve();area_path=args.area_input.resolve() if args.area_input else Path(data['path']).with_name('area.dat')
    if out.exists():raise ValueError('all-events requires a NEW output directory, including with --overwrite')
    if any(p==out or p.is_relative_to(out) for p in [Path(data['path']),area_path]+([args.cpp_reference.resolve()] if args.cpp_reference else [])):raise ValueError('output contains source input')
    area=read_minimum_area(area_path,data['dt_s'],args.area_column,args.area_time_column,args.area_unit) if area_path.is_file() else None
    if area and not area['t'][0]<=start<end<=area['t'][-1]+1/area['fs']+1e-12:raise ValueError('area lacks full requested window')
    provenance=closure_archive_provenance(data,area);series,branches,crossings,counts=all_event_extract(data,area,start,end)
    rows=[];periods=[];terms=[];events=[];half=[];sensitivity=[]
    cpp=None
    if args.cpp_reference:
        with args.cpp_reference.open(newline='') as f:references=list(csv.DictReader(f))
        matched=[r for r in references if r['case_id']==case]
        if len(matched)!=1:raise ValueError('CPP reference must contain exactly one matching case')
        ref=matched[0]
        if ref.get('input_sha256') and ref['input_sha256']!=data['sha256']:raise ValueError('CPP reference input SHA mismatch')
        if ref.get('start_s') and (float(ref['start_s'])!=start or float(ref['end_s'])!=end):raise ValueError('CPP reference window mismatch')
        cpp=float(ref['flow_CPP_dB'])
    for name in ALL_EVENT_METHODS:
        es=series[name];row,pp,tt=all_event_metrics(es,name,start,end);row.update(case_id=case,flow_CPP_reference_dB=cpp,CPP_recomputed=False,
            n_threshold_not_reached=sum(b['method']==name and b['reason']=='NO_THRESHOLD_CROSSING' for b in branches),
            n_multiple_crossing_branches=sum(b['method']==name and b['crossing_count']>1 for b in branches))
        if name.startswith('Amin') and area is None:row['status']='MISSING_AREA_INPUT';row['scan_fraction']=0.
        rows.append(row);periods.extend(dict(p,case_id=case) for p in pp);terms.extend(dict(t,case_id=case) for t in tt);events.extend(dict(e,case_id=case) for e in es)
        for label,left,right in (('early_half',start,(start+end)/2),('late_half',(start+end)/2,end)):
            hh,_,_=all_event_metrics([e for e in es if left<=e['event_time_s']<right],name,left,right);half.append(dict(hh,case_id=case,window_name=label,reference_only=True))
        if name.endswith(('max','min')):
            for rule in ('raw_time_s','plateau_start_s','plateau_center_s','plateau_end_s'):
                replacement=[]
                for e in es:
                    at=e['raw_time_s'] if rule=='raw_time_s' else .5*(e['plateau_start_s']+e['plateau_end_s']) if rule=='plateau_center_s' and e['plateau_width_s']>0 else e[rule] if rule!='plateau_center_s' and e['plateau_width_s']>0 else e['event_time_s']
                    replacement.append(dict(e,event_time_s=at))
                replacement.sort(key=lambda e:e['event_time_s']);ss,_,_=all_event_metrics(replacement,name,start,end);sensitivity.append(dict(ss,time_rule=rule,same_event_IDs=True,reference_only=True))
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.all_events_',dir=out.parent) as directory:
        temp=Path(directory);figures=temp/'figures';figures.mkdir()
        write_csv(temp/'comparison.csv',rows);write_csv(temp/'events.csv',events);write_csv(temp/'periods.csv',periods);write_csv(temp/'metric_terms.csv',terms)
        contributions=[r for r in terms if r['metric']=='jitter_local_percent'];write_csv(temp/'jitter_contributions.csv',contributions)
        write_csv(temp/'candidate_branches.csv',branches);write_csv(temp/'candidate_crossings.csv',crossings);write_csv(temp/'detector_counts.csv',counts)
        write_csv(temp/'half_window_reference.csv',half);write_csv(temp/'extremum_time_sensitivity.csv',sensitivity)
        metadata={'case_id':case,'input_path':data['path'],'input_sha256':data['sha256'],'area_path':area['path'] if area else str(area_path),'area_sha256':area['sha256'] if area else None,
            'dt_s':data['dt_s'],'flow_stride':data['stride'],'area_stride':area['stride'] if area else None,'start_s':start,'end_s':end,'window':'[start,end)',
            'detector_rule':ALL_EVENT_RULE,'plateau_rule':'one extremum at center; zero entry at region start / actual zero crossing',
            'threshold_rule':'3%,10% of EACH raw local peak; last downward crossing before next raw local valley; no hold gate',
            'area_definition':area['area_metadata'] if area else 'MISSING_AREA_INPUT','run_provenance':provenance,
            'nonpositive_flow_samples':int((data['x']<=0).sum()),'negative_flow_samples':int((data['x']<0).sum()),
            'CPP_reference':str(args.cpp_reference) if args.cpp_reference else 'NOT_PROVIDED','CPP_recomputed':False,'script_sha256':sha256(__file__),
            'manual_status':'MANUAL_VALIDATION_PENDING'}
        write_csv(temp/'analysis_metadata.csv',[{'key':k,'value':json.dumps(clean_json(v),ensure_ascii=False)} for k,v in metadata.items()])
        all_event_figures(figures,data,area,series,rows,contributions,start,end,case)
        lines=['# Independent all-detected-event intervals','',case+'; '+str(start)+'–'+str(end)+' s; AUTOMATED_NOT_MANUALLY_VERIFIED.','',
            'One sorted event vector per method. All M−1 intervals, including long gaps, use one mean denominator. No ACF, shape, DP, segmentation, hold-time gate or interval trimming. Event-span fraction is (last−first)/(end−start), not physical event occupancy. CPP is a saved reference only.','',
            '| Method | Events | Intervals | Jitter % | RAP % | PPQ5 % |','|---|---:|---:|---:|---:|---:|']
        for r in rows:lines.append('| '+r['method']+' | '+str(r['n_events'])+' | '+str(r['n_intervals'])+' | '+' | '.join('NOT_DEFINED_INSUFFICIENT_EVENTS' if r[k] is None else f'{r[k]:.6f}' for k in METRIC_NAMES)+' |')
        (temp/'report.md').write_text('\n'.join(lines)+'\n')
        for source in (data,area):
            if source and sha256(source['path'])!=source['sha256']:raise RuntimeError('source changed during all-events analysis')
        write_csv(temp/'output_manifest.csv',[{'path':str(p.relative_to(temp)),'sha256':sha256(p),'bytes':p.stat().st_size} for p in temp.rglob('*') if p.is_file()])
        if out.exists():raise ValueError('output appeared concurrently')
        temp.rename(out)
    print('all-events',case,'AUTOMATED_NOT_MANUALLY_VERIFIED; CPP not recomputed')
    for row in rows:print(row['method'],'M/N',row['n_events'],row['n_intervals'],'J/R/P',[row[k] for k in METRIC_NAMES])
    print('Saved:',out)
    return {'rows':rows,'metadata':metadata,'series':series}


ACF_GUIDED_ORDER = ('Ug_zero','Ug_max','Ug_min','Ug_3pct','Ug_10pct','Amin_max','Amin_min')


def acf_search_regions(start,end,frames,cfg):
    """Disjoint ownership regions; zero/many detections allowed, no phase assignment."""
    regions=[];cursor=float(start)
    while cursor<end:
        T,confidence,supported,flags=global_period_at(cursor,frames)
        scale=T if T is not None and T>0 else 1/cfg['pitch_min_hz']
        stop=min(end,cursor+scale)
        if not stop>cursor:raise RuntimeError('invalid search scale')
        regions.append({'region_id':len(regions),'start_s':cursor,'end_s':stop,'predicted_T0_s':T,'search_scale_s':scale,
            'scale_source':'LOCAL_ACF' if T else 'FALLBACK_1_FMIN_SEARCH_ONLY','acf_confidence':confidence,
            'acf_supported_diagnostic_only':supported,'acf_flags':flags})
        cursor=stop
    return regions


def acf_region_owner(at,regions):
    if not regions or not regions[0]['start_s']<=at<regions[-1]['end_s']:return None
    return int(np.searchsorted([r['end_s'] for r in regions],at,side='right'))


def acf_clear_extrema(data,start,end,frames,cfg,prefix,regions):
    # Reuse exactly the old local percentile / difference-MAD and valley/peak features.
    # Extract with the complete saved context, avoiding a find_peaks cut at region borders.
    raw=global_raw_candidates(data,float(data['t'][0]),float(data['t'][-1]+1/data['fs']),frames,cfg)
    accepted={};audit=[]
    for phase,kind in (('A','max'),('B','min')):
        features,_=global_candidate_features(data,raw,phase,cfg);clear=[]
        for e in features:
            reasons=[reason for reason in e['rejection_reasons'] if reason!='PERIOD_UNRESOLVED']
            at=.5*(e['plateau_start_s']+e['plateau_end_s']) if e['plateau_width_s']>0 else e['subsample_time_s']
            row={k:e[k] for k in ('sample_id','raw_time_s','prominence','local_amplitude','local_window_start_s','local_window_end_s',
                'local_window_samples','normalized_prominence','noise_sigma','plateau_start_s','plateau_end_s','plateau_width_s',
                'timing_uncertainty_s','local_period_s','local_acf_confidence','period_supported','local_period_flags',
                'peak_trough_amplitude','normalized_excursion','opposite_before_s','opposite_after_s',
                'local_prominence_threshold','local_excursion_threshold','MAD_prominence_floor','effective_prominence_threshold')}
            row.update(method=prefix+'_'+kind,event_id=prefix+'_'+kind+'_'+str(e['sample_id']),event_time_s=float(at),
                sample_value=float(data['x'][e['sample_id']]),step=int(data['steps'][e['sample_id']]),
                interpolation='PLATEAU_CENTER' if e['plateau_width_s']>0 else e['interpolation_status'],
                clarity_accepted=not reasons,detector_rejection_reasons=reasons,inside_window=bool(start<=at<end),
                search_region_id=acf_region_owner(at,regions),period_gate_applied=False,shape_gate_applied=False,
                peak_distance_gate_applied=False,dp_applied=False,
                excursion_context='TWO_SIDED_PROMINENT_OPPOSITE_EXTREMA' if e['opposite_before_s'] is not None and e['opposite_after_s'] is not None else 'BOUNDARY_ONE_SIDED_AND_PROMINENCE_FALLBACK')
            if not reasons:clear.append(row)
            if start-.05<=at<end+.05:audit.append(row)
        accepted[kind]=clear
    return accepted,audit


def acf_guided_extract(flow,area,start,end,cfg=None):
    cfg=dict(GLOBAL_CYCLE_CONFIG) if cfg is None else dict(cfg)
    raw_series,_,_,_=all_event_extract(flow,area,start,end)
    series={'Ug_zero':raw_series['Ug_zero']};audit=[];region_rows=[];frames_rows=[];period_candidates=[];signal_details={}
    for prefix,data in (('Ug',flow),('Amin',area)):
        if data is None:
            series.update(Amin_max=[],Amin_min=[]);continue
        frames,candidates=global_local_periods(data,start,end,cfg);regions=acf_search_regions(start,end,frames,cfg)
        clear,local_audit=acf_clear_extrema(data,start,end,frames,cfg,prefix,regions);audit.extend(local_audit)
        for kind in ('max','min'):series[prefix+'_'+kind]=[e for e in clear[kind] if start<=e['event_time_s']<end]
        for region in regions:
            group=[e for e in local_audit if e['search_region_id']==region['region_id']]
            region.update(signal=prefix,n_max_candidates=sum(e['method']==prefix+'_max' for e in group),
                n_min_candidates=sum(e['method']==prefix+'_min' for e in group),
                n_max_accepted=sum(e['method']==prefix+'_max' and e['clarity_accepted'] for e in group),
                n_min_accepted=sum(e['method']==prefix+'_min' and e['clarity_accepted'] for e in group))
        native=data['t'][(data['t']>=start)&(data['t']<end)];membership=np.zeros(len(native),int)
        for region in regions:membership+=(native>=region['start_s'])&(native<region['end_s'])
        if not np.all(membership==1):raise RuntimeError('ACF search regions lost or duplicated native samples')
        region_rows.extend(regions);frames_rows.extend(dict(f,signal=prefix) for f in frames);period_candidates.extend(dict(c,signal=prefix) for c in candidates)
        signal_details[prefix]={'clear':clear,'frames':frames,'regions':regions}
    # Each clear maximum owns its real descent, independently of T0 and other event series.
    # A following clear maximum is only a real-wave branch boundary, never a common pulse seed.
    branches=[];crossings=[];maxima=signal_details['Ug']['clear']['max'];minima=signal_details['Ug']['clear']['min'];t,x=flow['t'],flow['x']
    for pct in (3.,10.):
        name='Ug_'+format(pct,'g')+'pct';events=[];used=set()
        for j,peak in enumerate(maxima):
            i=peak['sample_id'];next_i=maxima[j+1]['sample_id'] if j+1<len(maxima) else len(x)-1
            opposite=next((e for e in minima if i<e['sample_id']<next_i),None)
            if opposite:
                right=int(np.searchsorted(t,opposite['plateau_start_s']));boundary='NEXT_CLEAR_LOCAL_MINIMUM'
            else:
                limit=max(i+1,next_i);right=i+int(np.argmin(x[i:limit+1]));boundary='OBSERVED_TROUGH_BEFORE_NEXT_CLEAR_MAX_OR_DATA_EDGE'
            if t[i]>=end or t[right]<start:continue
            options=closure_crossing_candidates(flow,i,right,pct,None);hit=options[-1] if options else None
            for rank,c in enumerate(options):crossings.append({'method':name,'peak_id':peak['event_id'],'peak_time_s':float(t[i]),
                'Upeak_m3_s':float(x[i]),'threshold_pct':pct,'threshold_m3_s':c['threshold_m3_s'],
                'sample_left':c['sample_left'],'sample_right':c['sample_right'],'value_left':c['Ug_left_m3_s'],'value_right':c['Ug_right_m3_s'],
                't_left_s':c['t_left_s'],'t_right_s':c['t_right_s'],'event_time_s':c['cross_time_s'],
                'crossing_rank':rank+1,'selected_last':rank==len(options)-1,'stability_gate_applied':False})
            branch={'method':name,'peak_id':peak['event_id'],'peak_sample_id':i,'peak_time_s':float(t[i]),'Upeak_m3_s':float(x[i]),
                'branch_end_sample':right,'branch_end_s':float(t[right]),'branch_boundary_rule':boundary,
                'valley_id':opposite['event_id'] if opposite else None,'next_clear_peak_id':maxima[j+1]['event_id'] if j+1<len(maxima) else None,
                'n_crossings':len(options),'first_crossing_s':options[0]['cross_time_s'] if options else None,'last_crossing_s':hit['cross_time_s'] if hit else None,
                'reason':'LAST_ACTUAL_DOWNCROSSING' if hit else 'NO_THRESHOLD_CROSSING'};branches.append(branch)
            if not hit or not start<=hit['cross_time_s']<end:continue
            bracket=(hit['sample_left'],hit['sample_right'])
            if bracket in used:raise RuntimeError('same threshold bracket assigned to two clear peaks')
            used.add(bracket)
            events.append({'method':name,'event_id':name+'_peak_'+str(i),'event_time_s':hit['cross_time_s'],'raw_time_s':hit['t_right_s'],
                'sample_id':hit['sample_right'],'step':hit['step_right'],'sample_value':hit['threshold_m3_s'],
                'sample_left':hit['sample_left'],'sample_right':hit['sample_right'],'value_left':hit['Ug_left_m3_s'],'value_right':hit['Ug_right_m3_s'],
                'linear_fraction':(hit['cross_time_s']-hit['t_left_s'])/(hit['t_right_s']-hit['t_left_s']),
                'peak_id':peak['event_id'],'peak_sample_id':i,'Upeak_m3_s':float(x[i]),'threshold_pct':pct,'threshold_m3_s':hit['threshold_m3_s'],
                'branch_end_sample':right,'branch_boundary_rule':boundary,'n_crossings':len(options),'interpolation':'LINEAR_LAST_DOWNCROSS',
                'search_region_id':acf_region_owner(hit['cross_time_s'],signal_details['Ug']['regions'])})
        series[name]=events
    for name,events in series.items():
        events.sort(key=lambda e:e['event_time_s'])
        if len({e['event_time_s'] for e in events})!=len(events):raise RuntimeError('duplicate ACF-guided event time')
    return {'series':series,'peak_candidates':audit,'regions':region_rows,'frames':frames_rows,
        'period_candidates':period_candidates,'branches':branches,'crossings':crossings,'signal_details':signal_details,'raw_series':raw_series}


def acf_guided_figures(out,flow,area,result,periods,contrib,start,end,case):
    series=result['series'];colors={'Ug_max':'tab:blue','Ug_min':'tab:green','Ug_zero':'black','Ug_3pct':'tab:orange','Ug_10pct':'tab:purple','Amin_max':'tab:blue','Amin_min':'tab:green'}
    def overlay(ax,names,data,left,right,raw_points=False):
        if data is None:ax.text(.5,.5,'MISSING_AREA_INPUT',transform=ax.transAxes,ha='center');return
        scale=1 if names[0].startswith('Amin') else 1e6;t,x=data['t'],data['x'];take=(t>=left)&(t<right)
        ax.plot(t[take],x[take]*scale,'k.-' if raw_points else 'k-',ms=3,lw=.8,label='saved raw signal')
        for name in names:
            es=[e for e in series[name] if left<=e['event_time_s']<right];times=[e['event_time_s'] for e in es]
            ax.scatter(times,np.interp(times,t,x)*scale,c=colors[name],s=25,marker='v' if name.endswith('min') else 'x' if name.endswith('pct') else 'o',label=name,zorder=4)
        ax.set(xlim=(left,right),xlabel='Time [s]',ylabel='Amin [mm²]' if scale==1 else 'Ug [mL/s]');ax.grid(alpha=.12);ax.legend(fontsize=7,ncol=3)
    for name,methods,data in [('flow_extrema_full',('Ug_max','Ug_min'),flow),('flow_closure_full',('Ug_zero','Ug_3pct','Ug_10pct'),flow),('area_extrema_full',('Amin_max','Amin_min'),area)]:
        fig,ax=plt.subplots(figsize=(16,4));overlay(ax,methods,data,start,end);fig.suptitle(case+' fixed local prominence/excursion 10%; all accepted intervals retained');fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=160);plt.close(fig)
    zz=[p for p in periods if p['method']=='Ug_zero'];jj=[p for p in contrib if p['method']=='Ug_zero']
    fig,axes=plt.subplots(2,1,figsize=(15,8),sharex=True)
    axes[0].plot([p['end_s'] for p in zz],[p['interval_s']*1000 for p in zz],'ko-',ms=3,lw=.7,label='ALL zero-entry intervals')
    frames=result['signal_details']['Ug']['frames'];axes[0].plot([f['center_s'] for f in frames],[f['period_s']*1000 if f['period_s'] else np.nan for f in frames],ls='--',c='gray',label='local ACF T0, diagnostic only')
    axes[0].set_ylabel('Zero-entry interval [ms]');axes[0].legend(fontsize=8)
    axes[1].plot([p['time_s'] for p in jj],[p['period_pair_jitter_contribution_percent'] for p in jj],'o-',ms=3,lw=.7,c='tab:red');axes[1].set(xlabel='Time [s]',ylabel='Zero-entry J_i [%]')
    for ax in axes:ax.set_xlim(start,end);ax.grid(alpha=.12)
    fig.suptitle(case+' primary Ug_zero: every interval and every J_i; no interval removal');fig.tight_layout();fig.savefig(out/'zero_intervals_and_contributions.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(4,1,figsize=(15,13))
    for ax,wanted in zip(axes[:2],('multiple','empty')):
        regions=[r for r in result['regions'] if r['signal']=='Ug' and (r['n_max_accepted']>=2 if wanted=='multiple' else r['n_max_accepted']==0)]
        data=flow;prefix='Ug'
        if not regions and area is not None:regions=[r for r in result['regions'] if r['signal']=='Amin' and (r['n_max_accepted']>=2 if wanted=='multiple' else r['n_max_accepted']==0)];data=area;prefix='Amin'
        if regions:
            reg=regions[0];pad=reg['search_scale_s']*.5;overlay(ax,(prefix+'_max',prefix+'_min'),data,max(start,reg['start_s']-pad),min(end,reg['end_s']+pad))
            for at in (reg['start_s'],reg['end_s']):ax.axvline(at,c='gray',ls='--',lw=.7)
            ax.set_title(f'{prefix} search region {reg["region_id"]}: {reg["n_max_accepted"]} clear maxima; no one-per-cycle rule',fontsize=10)
        else:ax.text(.5,.5,'No observed '+wanted+'-maximum search region in this dataset',transform=ax.transAxes,ha='center')
    for ax,which in zip(axes[2:],('shortest','longest')):
        if not zz:ax.text(.5,.5,'Insufficient zero-entry events',transform=ax.transAxes,ha='center');continue
        interval=(min if which=='shortest' else max)(zz,key=lambda p:p['interval_s']);pad=max(6/flow['fs'],.1*interval['interval_s'])
        overlay(ax,('Ug_zero',),flow,max(start,interval['start_s']-pad),min(end,interval['end_s']+pad),raw_points=which=='shortest')
        ax.set_title(f'{which} zero-entry interval: {interval["interval_s"]*1000:.6f} ms, ratio to T0={interval["interval_over_T0"]}; INCLUDED',fontsize=10)
    fig.suptitle(case+' representative cases; automated, manual validation pending');fig.tight_layout();fig.savefig(out/'representative_examples.png',dpi=155);plt.close(fig)


PEAK_SELECTORS = ('prominence','clear-flow','raw-flow')


def peak_associated_select(clear,regions,peak_selector='prominence',main_peak_overrides=None):
    """One observed main peak per region; nondefault selectors are sensitivity analyses."""
    if peak_selector not in PEAK_SELECTORS:raise ValueError('unknown peak selector')
    overrides=main_peak_overrides or {}
    if set(overrides)-{r['region_id'] for r in regions}:raise ValueError('override outside search regions')
    peaks=[];cycles=[]
    for region in regions:
        candidates=sorted((e for e in clear if e['search_region_id']==region['region_id']),
            key=lambda e:(-e['prominence'] if peak_selector=='prominence' else -e['sample_value'],e['event_time_s'],e['event_id']))
        main=candidates[0] if candidates else None
        overridden=region['region_id'] in overrides
        if overridden:
            main=next((e for e in candidates if e['event_id']==overrides[region['region_id']]),None)
            if main is None:raise ValueError('override is not an observed peak in this region')
        for rank,e in enumerate(candidates):
            row=dict(e,peak_role='MAIN' if e['event_id']==main['event_id'] else 'SECONDARY',prominence_rank=rank+1,
                parent_main_id=main['event_id'],main_selection_rule='MAX_PROMINENCE_THEN_EARLIEST_THEN_ID')
            if peak_selector!='prominence' or overridden:
                prominence_order=sorted(candidates,key=lambda c:(-c['prominence'],c['event_time_s'],c['event_id']))
                row.update(selection_rank=rank+1,prominence_rank=next(k+1 for k,c in enumerate(prominence_order) if c['event_id']==e['event_id']),
                    main_selection_rule='SINGLE_REGION_OBSERVED_OVERRIDE' if overridden else 'MAX_FLOW_THEN_EARLIEST_THEN_ID_'+peak_selector,
                    peak_selector=peak_selector,diagnostic_only=True)
            peaks.append(row)
        cycles.append(dict(region,main_peak_id=main['event_id'] if main else None,
            main_peak_time_s=main['event_time_s'] if main else None,Upeak_m3_s=main['sample_value'] if main else None,
            secondary_peak_count=max(0,len(candidates)-1),peak_status='OBSERVED' if main else 'NO_CLEAR_PEAK_NO_FILL'))
    return sorted(peaks,key=lambda e:e['event_time_s']),cycles


def peak_associated_descent_ids(x):
    """Physical non-increasing runs, separated only by real positive increments."""
    return np.r_[0,np.cumsum(np.diff(x)>0)].astype(int)


def flow_threshold_percent(value):
    """A positive downcross percentage; exact zero remains a separate event rule."""
    pct=float(value)
    if not np.isfinite(pct) or not 0<pct<100:
        raise ValueError('flow threshold percentage must be finite and strictly between 0 and 100')
    return pct


def flow_threshold_suffix(value):
    pct=flow_threshold_percent(value)
    token=str(int(pct)) if pct.is_integer() else repr(pct)
    return token.replace('.','p')+'pct'


def peak_associated_thresholds(value=3.):
    pct=flow_threshold_percent(value)
    values=[3,10]+([pct] if pct not in (3,10) else [])
    return [(p,flow_threshold_suffix(p)) for p in values]


def peak_associated_branch(flow,peak,boundary_s,raw_zero,descent_ids,flow_threshold_pct=3.):
    """All real reaches and a fixed representative; never apply hold/T0 gates."""
    t,x=flow['t'],flow['x'];i=peak['sample_id']
    right=min(len(t)-1,int(np.searchsorted(t,boundary_s,side='left')))
    right=max(i,right);observed_right=max(i,min(len(t)-1,int(np.searchsorted(t,boundary_s,side='left'))-1))
    minimum=i+int(np.argmin(x[i:observed_right+1]))
    zero=[dict(e) for e in raw_zero if t[i]<e['event_time_s']<boundary_s]
    anchor=zero[0]['sample_id'] if zero else minimum
    anchor_time=zero[0]['event_time_s'] if zero else float(t[anchor])
    options={'zero':zero};selected={};details=[];candidates=[]
    thresholds=peak_associated_thresholds(flow_threshold_pct)
    percentages=dict((suffix,pct) for pct,suffix in thresholds)
    for pct,suffix in thresholds:
        options[suffix]=[]
        for c in closure_crossing_candidates(flow,i,right,float(pct),None):
            if not c['cross_time_s']<boundary_s:continue
            options[suffix].append({'event_time_s':c['cross_time_s'],'sample_id':c['sample_right'],
                'sample_left':c['sample_left'],'sample_right':c['sample_right'],'value_left':c['Ug_left_m3_s'],
                'value_right':c['Ug_right_m3_s'],'threshold_m3_s':c['threshold_m3_s'],
                'raw_time_s':c['t_right_s'],'interpolation':'LINEAR_REAL_DOWNCROSS_NO_HOLD',
                'event_id':'cross_'+str(pct)+'_'+str(c['sample_left'])+'_'+peak['event_id']})
    for suffix,es in options.items():
        eligible=es if suffix=='zero' else [e for e in es if e['event_time_s']<=anchor_time]
        chosen=(eligible[0] if suffix=='zero' else eligible[-1]) if eligible else None
        rule='FIRST_EXACT_ZERO_ENTRY' if suffix=='zero' else 'LAST_DOWNCROSS_BEFORE_FIRST_ZERO' if zero else 'LAST_DOWNCROSS_TO_OBSERVED_GLOBAL_TROUGH_NO_ZERO'
        selected[suffix]=None
        for rank,c in enumerate(es):
            c.update(peak_id=peak['event_id'],peak_role=peak['peak_role'],search_region_id=peak['search_region_id'],
                peak_time_s=peak['event_time_s'],peak_raw_time_s=float(t[i]),Upeak_m3_s=float(x[i]),
                threshold_name=suffix,threshold_pct=0 if suffix=='zero' else percentages[suffix],
                threshold_m3_s=0. if suffix=='zero' else c['threshold_m3_s'],
                candidate_rank=rank+1,descent_run_id=int(descent_ids[c['sample_left']]),
                branch_end_s=boundary_s,terminal_anchor_s=anchor_time,terminal_anchor_kind='FIRST_ZERO' if zero else 'OBSERVED_GLOBAL_MINIMUM',
                selected_by_fixed_rule=c is chosen,selection_rule=rule,
                candidate_id=peak['event_id']+'_'+suffix+'_'+str(rank+1))
            candidates.append(c)
            if c is chosen:selected[suffix]=c
        # With a positive peak, a second downcross implies an intervening real upcross.
        details.append({'peak_id':peak['event_id'],'peak_role':peak['peak_role'],'search_region_id':peak['search_region_id'],
            'threshold_name':suffix,'threshold_pct':0 if suffix=='zero' else percentages[suffix],
            'branch_start_s':float(t[i]),'branch_end_s':boundary_s,'first_reach_s':es[0]['event_time_s'] if es else None,
            'last_reach_s':es[-1]['event_time_s'] if es else None,'reach_count':len(es),
            'rises_above_after_first':bool(es and np.any(x[es[0]['sample_right']+1:observed_right+1]>es[0]['threshold_m3_s'])),
            'repeated_downcross_after_rise':len(es)>1,'selected_reach_s':chosen['event_time_s'] if chosen else None,
            'first_vs_selected_delta_s':chosen['event_time_s']-es[0]['event_time_s'] if chosen else None,
            'selection_rule':rule,'status':'REACHED' if chosen else 'NOT_REACHED_TO_TERMINAL_ANCHOR',
            'terminal_anchor_s':anchor_time,'terminal_anchor_kind':'FIRST_ZERO' if zero else 'OBSERVED_GLOBAL_MINIMUM'})
    return {'selected':selected,'options':options,'details':details,'candidates':candidates,
        'post_peak_min_m3_s':float(x[i:observed_right+1].min()),'post_peak_max_m3_s':float(x[i:observed_right+1].max()),
        'post_peak_min_time_s':float(t[minimum]),'branch_end_s':boundary_s}


def peak_associated_extract(flow,start,end,base,flow_threshold_pct=3.,peak_selector='prominence',main_peak_overrides=None):
    """Main branches are peak-to-next-main. Secondary branches use next clear peak.

    Candidate history is not ownership. Main representative descents have priority;
    a physical falling run can belong to only one peak in the auxiliary union.
    """
    thresholds=peak_associated_thresholds(flow_threshold_pct)
    suffixes=('zero',*[suffix for pct,suffix in thresholds])
    clear=base['series']['Ug_max'];regions=[r for r in base['regions'] if r['signal']=='Ug']
    if peak_selector=='raw-flow' or main_peak_overrides:
        raw=[e for e in base['peak_candidates'] if e['method']=='Ug_max' and e['inside_window']]
        if peak_selector=='raw-flow':clear=raw
        elif main_peak_overrides:
            lookup={e['event_id']:e for e in raw};present={e['event_id'] for e in clear};clear=list(clear)
            for rid,identity in main_peak_overrides.items():
                if identity not in lookup or lookup[identity]['search_region_id']!=rid:raise ValueError('override is not a recorded raw peak in its region')
                if identity not in present:clear.append(lookup[identity]);present.add(identity)
    peaks,cycles=peak_associated_select(clear,regions,peak_selector,main_peak_overrides);mains=[p for p in peaks if p['peak_role']=='MAIN']
    raw_zero=base['series']['Ug_zero'];descent_ids=peak_associated_descent_ids(flow['x'])
    branches={};history=[];candidates=[];owned={};series={};order=[]
    def event(c,name):
        return dict(c,method=name,event_id=name+'_'+c['candidate_id'])
    for suffix in suffixes:
        for scope in ('main','including_secondary'):series[scope+'_'+suffix]=[]
    for pct,suffix in thresholds:series['first_crossing_'+suffix]=[]
    series['last_zero_sensitivity']=[]
    for j,p in enumerate(mains):
        boundary=min(end,mains[j+1]['raw_time_s']) if j+1<len(mains) else end
        b=peak_associated_branch(flow,p,boundary,raw_zero,descent_ids,flow_threshold_pct);branches[p['event_id']]=b
        history.extend(b['details']);candidates.extend(b['candidates'])
        cycle=cycles[p['search_region_id']]
        cycle.update(next_main_peak_id=mains[j+1]['event_id'] if j+1<len(mains) else None,
            next_main_peak_time_s=mains[j+1]['event_time_s'] if j+1<len(mains) else None,
            next_main_peak_interval_s=mains[j+1]['event_time_s']-p['event_time_s'] if j+1<len(mains) else None,
            branch_end_s=boundary,right_censored=j+1==len(mains),
            post_peak_min_m3_s=b['post_peak_min_m3_s'],post_peak_max_m3_s=b['post_peak_max_m3_s'],
            minimum_to_Upeak_ratio=b['post_peak_min_m3_s']/p['sample_value'],
            minimum_to_post_peak_maximum_ratio=b['post_peak_min_m3_s']/b['post_peak_max_m3_s'])
        for suffix in suffixes:
            c=b['selected'][suffix];d=next(v for v in b['details'] if v['threshold_name']==suffix)
            cycle.update({suffix+'_reached':c is not None,suffix+'_reach_count':d['reach_count'],suffix+'_representative_s':c['event_time_s'] if c else None})
            if c:
                series['main_'+suffix].append(event(c,'main_'+suffix))
                run=c['descent_run_id']
                if run in owned and owned[run]!=p['event_id']:raise RuntimeError('two main peaks own the same physical descent')
                owned[run]=p['event_id']
            if suffix!='zero' and b['options'][suffix]:series['first_crossing_'+suffix].append(event(b['options'][suffix][0],'first_crossing_'+suffix))
        if b['options']['zero']:series['last_zero_sensitivity'].append(event(b['options']['zero'][-1],'last_zero_sensitivity'))
        at={s:b['selected'][s]['event_time_s'] if b['selected'][s] else None for s in ('zero','3pct','10pct')}
        violations=[]
        for earlier,later in (('10pct','3pct'),('3pct','zero'),('10pct','zero')):
            if at[earlier] is not None and at[later] is not None and at[earlier]>at[later]+1e-12:violations.append(earlier+'_AFTER_'+later)
        order.append({'peak_id':p['event_id'],'peak_time_s':p['event_time_s'],'search_region_id':p['search_region_id'],
            'zero_s':at['zero'],'pct3_s':at['3pct'],'pct10_s':at['10pct'],'order_violations':violations,
            'status':'ORDER_CONTRADICTION' if violations else 'ORDER_OK' if at['zero'] is not None else 'NO_ZERO_OBSERVED_TROUGH_REFERENCE',
            'selected_descent_runs':{s:b['selected'][s]['descent_run_id'] if b['selected'][s] else None for s in at},
            'all_three_on_one_monotone_run':all(at[s] is not None for s in at) and len({b['selected'][s]['descent_run_id'] for s in at if b['selected'][s]})==1})
    for suffix in suffixes:series['including_secondary_'+suffix]=[dict(e,method='including_secondary_'+suffix,event_id='including_secondary_'+e['event_id']) for e in series['main_'+suffix]]
    for j,p in enumerate(peaks):
        if p['peak_role']=='MAIN':continue
        boundary=min(end,peaks[j+1]['raw_time_s']) if j+1<len(peaks) else end
        b=peak_associated_branch(flow,p,boundary,raw_zero,descent_ids,flow_threshold_pct);branches[p['event_id']]=b
        history.extend(b['details']);candidates.extend(b['candidates'])
        # Decide ownership for the entire physical run before emitting any threshold.
        selected=b['selected'];z=selected['zero']
        same_closure_already_owned=bool(z and z['descent_run_id'] in owned and owned[z['descent_run_id']]!=p['event_id'])
        free={c['descent_run_id'] for c in selected.values() if c and c['descent_run_id'] not in owned and not same_closure_already_owned}
        for run in free:owned[run]=p['event_id']
        for suffix,c in selected.items():
            if c and not same_closure_already_owned and owned.get(c['descent_run_id'])==p['event_id']:series['including_secondary_'+suffix].append(event(c,'including_secondary_'+suffix))
        for d in b['details']:
            c=selected[d['threshold_name']]
            d['auxiliary_adoption']='ADOPTED' if c and not same_closure_already_owned and owned.get(c['descent_run_id'])==p['event_id'] else 'CLOSURE_EPISODE_ALREADY_OWNED' if c and same_closure_already_owned else 'DESCENT_ALREADY_OWNED_BY_MAIN_OR_OTHER_PEAK' if c else 'NOT_REACHED'
        for c in b['candidates']:c['auxiliary_episode_conflict']=same_closure_already_owned
    for c in candidates:
        c['assigned_owner_peak_id']=owned.get(c['descent_run_id'])
        c['representative_owned_by_this_peak']=c['selected_by_fixed_rule'] and c['assigned_owner_peak_id']==c['peak_id'] and not c.get('auxiliary_episode_conflict',False)
    zero_match=[]
    for z in raw_zero:
        parents=[p for p in mains if float(flow['t'][p['sample_id']])<z['event_time_s']<branches[p['event_id']]['branch_end_s']]
        parent=parents[0] if parents else None;run=int(descent_ids[z['sample_left']]);owner=owned.get(run)
        primary=parent and branches[parent['event_id']]['selected']['zero']
        selected=bool(primary and primary['sample_id']==z['sample_id'])
        zero_match.append(dict(z,parent_main_id=parent['event_id'] if parent else None,parent_main_time_s=parent['event_time_s'] if parent else None,
            assigned_owner_peak_id=owner,main_representative=selected,
            reason='FIRST_ZERO_OF_MAIN_PULSE' if selected else 'REPEATED_ZERO_IN_MAIN_BRANCH' if parent else 'LEFT_WINDOW_EDGE_BEFORE_FIRST_MAIN_NO_FILL'))
    for name,es in series.items():
        es.sort(key=lambda e:e['event_time_s'])
        if len({e['event_time_s'] for e in es})!=len(es):raise RuntimeError('duplicate closure representative')
    return {'series':series,'peaks':peaks,'cycles':cycles,'reach_history':history,'crossing_candidates':candidates,
        'order_diagnostics':order,'zero_matching':zero_match,'branches':branches,'base':base}


def peak_associated_figures(out,flow,result,periods,terms,rows,start,end,case):
    t,x=flow['t'],flow['x'];series=result['series'];peaks=result['peaks']
    def overlay(ax,left,right,regions=False):
        take=(t>=left)&(t<right);ax.plot(t[take],x[take]*1e6,'k-',lw=.8,label='saved Ug')
        for role,marker,color in (('MAIN','^','tab:blue'),('SECONDARY','v','tab:gray')):
            pp=[p for p in peaks if p['peak_role']==role and left<=p['event_time_s']<right]
            ax.scatter([p['event_time_s'] for p in pp],[p['sample_value']*1e6 for p in pp],s=25,marker=marker,c=color,label=role,zorder=4)
        for suffix,color,marker in (('zero','black','o'),('3pct','tab:orange','x'),('10pct','tab:green','+')):
            es=[e for e in series['main_'+suffix] if left<=e['event_time_s']<right];at=[e['event_time_s'] for e in es]
            ax.scatter(at,[e['threshold_m3_s']*1e6 for e in es],s=35,marker=marker,c=color,label='main '+suffix,zorder=5)
        if regions:
            for r in result['cycles']:
                if left<=r['start_s']<right:ax.axvline(r['start_s'],ls='--',c='gray',lw=.6)
        ax.set(xlim=(left,right),xlabel='Time [s]',ylabel='Ug [mL/s]');ax.grid(alpha=.12);ax.legend(fontsize=7,ncol=3)
    fig,ax=plt.subplots(figsize=(16,4));overlay(ax,start,end);ax.set_title(case+' main-peak-associated closure; fixed rules; manual validation pending');fig.tight_layout();fig.savefig(out/'full_window.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(14,11))
    multi=next((r for r in result['cycles'] if r['secondary_peak_count']),None)
    r=multi or result['cycles'][0];pad=2*r['search_scale_s'];overlay(axes[0],max(start,r['start_s']-pad),min(end,r['end_s']+pad),True);axes[0].set_title('ACF search regions: maximum prominence selects MAIN; all SECONDARY peaks retained')
    complete=[r for r in result['order_diagnostics'] if r['zero_s'] is not None and r['pct3_s'] is not None and r['pct10_s'] is not None]
    differing={d['peak_id'] for d in result['reach_history'] if d['peak_role']=='MAIN' and d['threshold_name']!='zero' and d['first_vs_selected_delta_s'] and d['first_vs_selected_delta_s']>0}
    representative=next((r for r in complete if r['peak_id'] in differing),complete[0] if complete else result['order_diagnostics'][0] if result['order_diagnostics'] else None)
    if representative:
        p=next(p for p in peaks if p['event_id']==representative['peak_id']);b=result['branches'][p['event_id']]
        overlay(axes[1],max(start,p['raw_time_s']-.001),min(end,b['branch_end_s']+.001))
        for suffix,color in (('3pct','tab:orange'),('10pct','tab:green')):
            es=b['options'][suffix]
            if es:axes[1].scatter([es[0]['event_time_s']],[es[0]['threshold_m3_s']*1e6],facecolors='none',edgecolors=color,s=75,label='first '+suffix)
        axes[1].legend(fontsize=7,ncol=3);axes[1].set_title('Same main pulse: first crossings (open circles) versus terminal-to-first-zero selections')
    repeated=[d for d in result['reach_history'] if d['peak_role']=='MAIN' and d['threshold_name']=='zero' and d['reach_count']>1]
    special=next((d for d in repeated if d['first_reach_s']<=.95955+1e-9 and d['last_reach_s']>=.95965-1e-9),repeated[0] if repeated else None)
    if special:
        b=result['branches'][special['peak_id']];es=b['options']['zero'];span=es[-1]['event_time_s']-es[0]['event_time_s'];pad=max(.0003,.3*span)
        left=max(start,es[0]['event_time_s']-pad);right=min(end,es[-1]['event_time_s']+pad);overlay(axes[2],left,right)
        take=(t>=left)&(t<right);axes[2].plot(t[take],x[take]*1e6,'k.',ms=3)
        axes[2].scatter([e['event_time_s'] for e in es],np.zeros(len(es)),facecolors='none',edgecolors='tab:red',s=75,label='all exact zero entries')
        axes[2].legend(fontsize=7,ncol=3);axes[2].set_title('Repeated zero entries in one main branch: first selected, complete history retained')
    else:axes[2].text(.5,.5,'No repeated zero entry in a main branch',ha='center',transform=axes[2].transAxes)
    fig.tight_layout();fig.savefig(out/'association_and_repeated_zero.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(14,7),sharex=True)
    pp=[p for p in periods if p['method']=='main_zero'];jj=[r for r in terms if r['method']=='main_zero' and r['metric']==METRIC_NAMES[0]]
    axes[0].plot([p['end_s'] for p in pp],[p['interval_s']*1000 for p in pp],'ko-',ms=3,lw=.8);axes[0].set_ylabel('Main zero interval [ms]')
    axes[1].plot([p['time_s'] for p in jj],[p['normalized_term_percent'] for p in jj],'o-',c='tab:red',ms=3,lw=.8);axes[1].set(xlabel='Time [s]',ylabel='J_i [%]')
    for ax in axes:ax.set_xlim(start,end);ax.grid(alpha=.15)
    fig.suptitle(case+' all actual main-zero intervals and all Jitter terms; no gap splitting');fig.tight_layout();fig.savefig(out/'main_zero_intervals_and_contributions.png',dpi=160);plt.close(fig)
    names=['reference_Ug_zero','main_zero','main_3pct','main_10pct'];lookup={r['method']:r for r in rows}
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,key in zip(axes,METRIC_NAMES):
        ax.bar(range(4),[lookup[n][key] if lookup[n][key] is not None else np.nan for n in names]);ax.set_xticks(range(4),['all zero','main zero','main 3%','main 10%'],rotation=20);ax.set(ylabel='[%]',title=key);ax.grid(axis='y',alpha=.12)
    fig.suptitle(case+' event definitions differ; auxiliary sensitivity is CSV only');fig.tight_layout();fig.savefig(out/'comparison.png',dpi=160);plt.close(fig)
    contradictions=[d for d in result['order_diagnostics'] if d['order_violations']]
    if contradictions:
        fig,axes=plt.subplots(len(contradictions),1,figsize=(14,3*len(contradictions)),squeeze=False)
        for ax,d in zip(axes[:,0],contradictions):
            p=next(p for p in peaks if p['event_id']==d['peak_id']);overlay(ax,max(start,p['raw_time_s']-.001),min(end,result['branches'][p['event_id']]['branch_end_s']+.001));ax.set_title(str(d['order_violations']))
        fig.tight_layout();fig.savefig(out/'order_contradictions.png',dpi=150);plt.close(fig)


def peak_associated_cli(args):
    if args.input is None:raise ValueError('--input is required')
    if args.primary_prominence!=10 or args.threshold_profile!='balanced-v1' or args.fmin!=60 or args.fmax!=400:raise ValueError('peak-associated-closure uses unchanged fixed 10% balanced-v1 ACF/clarity')
    if args.manual_events:raise ValueError('manual annotations cannot gate this event definition')
    flow=read_flow(args.input,resolve_dt(args.dt,args.manifest));start=.5 if args.start is None else args.start;end=1. if args.end is None else args.end
    if not np.isfinite([start,end]).all() or not flow['t'][0]<=start<end<=flow['t'][-1]+1/flow['fs']+1e-12:raise ValueError('invalid peak-associated window')
    if np.any(flow['x']<0):raise ValueError('exact zero association requires nonnegative Ug; signed-flow semantics must be specified separately')
    out=(args.out_dir or Path('analysis_runs/jitter_all_events_20261107')/('peak_associated_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve()
    area_path=args.area_input.resolve() if args.area_input else Path(flow['path']).with_name('area.dat')
    if out.exists():raise ValueError('peak-associated-closure requires a NEW output directory')
    archive_root=next((p for p in Path(flow['path']).parents if p.name=='manual_run_archive'),None)
    if archive_root and out.is_relative_to(archive_root):raise ValueError('manual_run_archive is read-only')
    if any(p==out or p.is_relative_to(out) for p in (Path(flow['path']),area_path)):raise ValueError('output contains source input')
    area=read_minimum_area(area_path,flow['dt_s'],args.area_column,args.area_time_column,args.area_unit) if area_path.is_file() else None
    provenance=closure_archive_provenance(flow,area);case=args.case_id or next((s.split('_P')[0] for s in Path(flow['path']).parts if s.startswith('r') and '_P' in s),Path(flow['path']).stem)
    base=acf_guided_extract(flow,area,start,end);result=peak_associated_extract(flow,start,end,base)
    series=dict(result['series'])
    series.update({'reference_'+name:es for name,es in base['series'].items()})
    rows=[];periods=[];terms=[];events=[]
    for name,es in series.items():
        # Every series uses the unchanged all-event metric kernel and its full interval vector.
        row,pp,tt=all_event_metrics(es,name,start,end);scope='MAIN' if name.startswith('main_') else 'REFERENCE_ACF_GUIDED' if name.startswith('reference_') else 'AUXILIARY_NUMERIC_ONLY'
        row.update(case_id=case,scope=scope,estimator='ALL_ACTUAL_ASSOCIATED_EVENTS_NO_INTERVAL_GATE');rows.append(row)
        periods.extend(dict(p,case_id=case,scope=scope) for p in pp);terms.extend(dict(r,case_id=case,scope=scope) for r in tt)
        events.extend(dict(e,method=name,case_id=case,scope=scope) for e in es)
    main_by_id={p['event_id']:p for p in result['peaks'] if p['peak_role']=='MAIN'}
    event_lookup={e['event_id']:e for e in events if e['method'].startswith('main_')}
    for p in periods:
        T,confidence,supported,flags=global_period_at((p['start_s']+p['end_s'])/2,base['signal_details']['Ug']['frames'])
        p.update(predicted_T0_s=T,interval_over_T0=p['interval_s']/T if T else None,interval_gate_applied=False)
        if p['scope']=='MAIN':
            ea,eb=event_lookup[p['event_start']],event_lookup[p['event_end']]
            p['search_regions_between_events']=eb['search_region_id']-ea['search_region_id']-1
            p['main_peaks_without_same_event_between']=sum(ea['peak_raw_time_s']<v['raw_time_s']<eb['peak_raw_time_s'] for v in main_by_id.values())
    reference_checks=[]
    if args.comparison_reference:
        ref=args.comparison_reference.resolve();meta={r['key']:json.loads(r['value']) for r in csv.DictReader((ref.parent/'analysis_metadata.csv').open())}
        if meta['input_sha256']!=flow['sha256'] or meta['start_s']!=start or meta['end_s']!=end:raise ValueError('comparison reference input/window mismatch')
        prior=list(csv.DictReader(ref.open()))
        for old in prior:
            new=next(r for r in rows if r['method']=='reference_'+old['method'])
            for key in ('n_events','n_intervals',*METRIC_NAMES,'mean_period_s','min_interval_s','max_interval_s'):
                expected=float(old[key]) if old[key] not in ('','None') else None
                if new[key]!=expected:raise RuntimeError('existing ACF result changed: '+old['method']+' '+key)
            reference_checks.append({'method':old['method'],'all_counts_and_metrics_exact':'PASS'})
    metadata={'case_id':case,'input_path':flow['path'],'input_sha256':flow['sha256'],'area_path':area['path'] if area else None,'area_sha256':area['sha256'] if area else None,
        'start_s':start,'end_s':end,'dt_s':flow['dt_s'],'stride':flow['stride'],'script_sha256':sha256(__file__),'provenance':provenance,
        'main_selection':'maximum prominence per unchanged half-open ACF search region; ties earliest then ID',
        'zero_rule':'first exact positive-to-zero entry in main-peak-to-next-main-peak branch; no epsilon/hold/reopening merge in candidate history',
        'threshold_rule':'last real downcross before first exact zero; if no zero, before first observed global branch minimum; no stability gate',
        'branch_boundary':'next main RAW peak sample exclusive; last branch end exclusive, right censored; left edge before first main unassigned',
        'secondary_rule':'next clear peak bounds secondary branch; main-selected physical falling runs own priority; no replacement on ownership conflict',
        'comparison_reference':str(args.comparison_reference) if args.comparison_reference else None,'CPP_recomputed':False,'manual_status':'MANUAL_VALIDATION_PENDING'}
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.peak_associated_',dir=out.parent) as directory:
        temp=Path(directory);(temp/'figures').mkdir()
        for name,data in [('comparison',[r for r in rows if r['scope']!='AUXILIARY_NUMERIC_ONLY']),('auxiliary_comparison',[r for r in rows if r['scope']=='AUXILIARY_NUMERIC_ONLY']),
            ('events',events),('periods',periods),('metric_terms',terms),('jitter_contributions',[r for r in terms if r['metric']==METRIC_NAMES[0]]),
            ('peaks',result['peaks']),('cycles',result['cycles']),('reach_history',result['reach_history']),('crossing_candidates',result['crossing_candidates']),
            ('order_diagnostics',result['order_diagnostics']),('zero_matching',result['zero_matching']),('reference_checks',reference_checks),
            ('peak_candidates',base['peak_candidates']),('local_period_frames',base['frames'])]:write_csv(temp/(name+'.csv'),data)
        write_csv(temp/'analysis_metadata.csv',[{'key':k,'value':json.dumps(clean_json(v),ensure_ascii=False)} for k,v in metadata.items()])
        peak_associated_figures(temp/'figures',flow,result,periods,terms,rows,start,end,case)
        for d in (flow,area):
            if d and sha256(d['path'])!=d['sha256']:raise RuntimeError('input changed during association')
        write_csv(temp/'output_manifest.csv',[{'path':str(p.relative_to(temp)),'sha256':sha256(p),'bytes':p.stat().st_size} for p in temp.rglob('*') if p.is_file()])
        if out.exists():raise ValueError('output appeared concurrently')
        temp.rename(out)
    print(case,'peak-associated-closure; MANUAL_VALIDATION_PENDING')
    for r in rows:
        if r['scope']=='MAIN':print(r['method'],'M/N',r['n_events'],r['n_intervals'],'J/R/P',[r[k] for k in METRIC_NAMES])
    print('Saved:',out)
    return {'rows':rows,'metadata':metadata,'result':result}


def acf_guided_cli(args):
    if args.input is None:raise ValueError('--input is required')
    if args.primary_prominence!=10 or args.threshold_profile!='balanced-v1' or args.fmin!=60 or args.fmax!=400:raise ValueError('acf-guided-events uses fixed 10% clarity and existing balanced-v1 60–400 Hz scale; no optimization')
    if args.manual_events:raise ValueError('manual annotations must not gate this analysis; inspect exported events independently')
    flow=read_flow(args.input,resolve_dt(args.dt,args.manifest));start=.5 if args.start is None else args.start;end=1. if args.end is None else args.end
    if not np.isfinite([start,end]).all() or not flow['t'][0]<=start<end<=flow['t'][-1]+1/flow['fs']+1e-12:raise ValueError('invalid ACF-guided window')
    out=(args.out_dir or Path('analysis_runs/jitter_all_events_20261107')/('acf_guided_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve()
    area_path=args.area_input.resolve() if args.area_input else Path(flow['path']).with_name('area.dat')
    if out.exists():raise ValueError('acf-guided-events requires a NEW output directory')
    if any(p==out or p.is_relative_to(out) for p in [Path(flow['path']),area_path]+([args.cpp_reference.resolve()] if args.cpp_reference else [])):raise ValueError('output contains protected source')
    area=read_minimum_area(area_path,flow['dt_s'],args.area_column,args.area_time_column,args.area_unit) if area_path.is_file() else None
    if area and not area['t'][0]<=start<end<=area['t'][-1]+1/area['fs']+1e-12:raise ValueError('area lacks full window')
    provenance=closure_archive_provenance(flow,area);case=args.case_id or next((s.split('_P')[0] for s in Path(flow['path']).parts if s.startswith('r') and '_P' in s),Path(flow['path']).stem)
    cfg=dict(GLOBAL_CYCLE_CONFIG);result=acf_guided_extract(flow,area,start,end,cfg);rows=[];periods=[];terms=[];events=[];counts=[];zero_details=[]
    cpp=None
    if args.cpp_reference:
        with args.cpp_reference.open(newline='') as f:refs=[r for r in csv.DictReader(f) if r['case_id']==case]
        if len(refs)!=1 or refs[0].get('input_sha256') not in (None,'',flow['sha256']):raise ValueError('invalid saved CPP reference')
        if refs[0].get('start_s') and (float(refs[0]['start_s'])!=start or float(refs[0]['end_s'])!=end):raise ValueError('CPP reference window mismatch')
        cpp=float(refs[0]['flow_CPP_dB'])
    for name in ACF_GUIDED_ORDER:
        es=result['series'][name];row,pp,tt=all_event_metrics(es,name,start,end);row.update(case_id=case,estimator='J_all_clear_detected_zero_independent',flow_CPP_reference_dB=cpp,CPP_recomputed=False)
        if name.startswith('Amin') and area is None:row['status']='MISSING_AREA_INPUT';row['scan_fraction']=0.
        rows.append(row);events.extend(dict(e,case_id=case) for e in es);terms.extend(dict(t,case_id=case) for t in tt)
        frames=result['signal_details']['Amin' if name.startswith('Amin') and area is not None else 'Ug']['frames']
        for p in pp:
            T,confidence,supported,flags=global_period_at((p['start_s']+p['end_s'])/2,frames)
            p.update(case_id=case,predicted_T0_s=T,interval_over_T0=p['interval_s']/T if T else None,acf_confidence_diagnostic_only=confidence,
                period_supported_diagnostic_only=supported,interval_gate_applied=False)
            periods.append(p)
            if name=='Ug_zero':zero_details.append(p)
        cc=[e for e in result['peak_candidates'] if e['method']==name and e['inside_window']]
        counts.append({'case_id':case,'method':name,'raw_candidates':len(cc) if name.endswith(('max','min')) else len(es),
            'accepted_events':len(es),'rejected_by_fixed_detector':sum(not e['clarity_accepted'] for e in cc),'interval_excluded':0,
            'zero_independent_of_clarity':name=='Ug_zero','accepted_all_intervals_used':True})
    if result['series']['Ug_zero']!=result['raw_series']['Ug_zero']:raise RuntimeError('zero-entry events changed')
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.acf_guided_',dir=out.parent) as directory:
        temp=Path(directory);figures=temp/'figures';figures.mkdir()
        write_csv(temp/'comparison.csv',rows);write_csv(temp/'events.csv',events);write_csv(temp/'periods.csv',periods);write_csv(temp/'metric_terms.csv',terms)
        contributions=[t for t in terms if t['metric']=='jitter_local_percent'];write_csv(temp/'jitter_contributions.csv',contributions)
        write_csv(temp/'peak_candidates.csv',result['peak_candidates']);write_csv(temp/'detector_counts.csv',counts);write_csv(temp/'search_regions.csv',result['regions'])
        write_csv(temp/'local_period_frames.csv',result['frames']);write_csv(temp/'local_period_candidates.csv',result['period_candidates']);write_csv(temp/'zero_interval_details.csv',zero_details)
        write_csv(temp/'candidate_branches.csv',result['branches']);write_csv(temp/'candidate_crossings.csv',result['crossings'])
        metadata={'case_id':case,'input_path':flow['path'],'input_sha256':flow['sha256'],'area_path':area['path'] if area else str(area_path),'area_sha256':area['sha256'] if area else None,
            'dt_s':flow['dt_s'],'flow_stride':flow['stride'],'area_stride':area['stride'] if area else None,'start_s':start,'end_s':end,
            'prominence_fraction':.1,'excursion_fraction':.1,'local_amplitude':'P95-P5 over +/-1.5T0; original data outside search regions available',
            'noise_floor':'2 * difference MAD sigma, sigma=1.4826*median(abs(diff-median(diff)))/sqrt(2); variable absolute per candidate',
            'opposite_extrema':'nearest preceding/following opposite extrema meeting preliminary prominence 10% and MAD floor, same as existing global_candidate_features',
            'T0_missing_fallback':'1/fmin=16.6667 ms for search/local-amplitude scale only; no ACF-support adoption gate',
            'area_definition':area['area_metadata'] if area else 'MISSING_AREA_INPUT','provenance':provenance,
            'shape_or_DP_or_interval_gate':False,'zero_rule':'unchanged all_event_extract exact positive-to-nonpositive; no epsilon/hold/merging',
            'negative_flow_samples':int((flow['x']<0).sum()),'script_sha256':sha256(__file__),'CPP_recomputed':False,'CPP_reference':str(args.cpp_reference) if args.cpp_reference else None,
            'manual_status':'MANUAL_VALIDATION_PENDING'}
        write_csv(temp/'analysis_metadata.csv',[{'key':k,'value':json.dumps(clean_json(v),ensure_ascii=False)} for k,v in metadata.items()])
        acf_guided_figures(figures,flow,area,result,periods,contributions,start,end,case)
        lines=['# ACF-scale, fixed 10% clear events; independent zero-entry Jitter','',case+'; primary=Ug_zero; MANUAL_VALIDATION_PENDING.','',
            'Local ACF guides search scale only. Every region permits zero/one/multiple clear events; every accepted interval is retained without segmentation or length/shape/ACF gating. Zero entries are unchanged.','',
            '| Method | Events | Intervals | Jitter % | RAP % | PPQ5 % |','|---|---:|---:|---:|---:|---:|']
        for row in rows:lines.append('| '+row['method']+' | '+str(row['n_events'])+' | '+str(row['n_intervals'])+' | '+' | '.join('UNDEFINED' if row[k] is None else f'{row[k]:.6f}' for k in METRIC_NAMES)+' |')
        (temp/'report.md').write_text('\n'.join(lines)+'\n')
        for data in (flow,area):
            if data and sha256(data['path'])!=data['sha256']:raise RuntimeError('source changed during ACF-guided analysis')
        write_csv(temp/'output_manifest.csv',[{'path':str(p.relative_to(temp)),'sha256':sha256(p),'bytes':p.stat().st_size} for p in temp.rglob('*') if p.is_file()])
        if out.exists():raise ValueError('output appeared concurrently')
        temp.rename(out)
    print('acf-guided-events',case,'PRIMARY Ug_zero; MANUAL_VALIDATION_PENDING')
    for row in rows:print(row['method'],'M/N',row['n_events'],row['n_intervals'],'J/R/P',[row[k] for k in METRIC_NAMES])
    print('Saved:',out)
    return {'rows':rows,'metadata':metadata,'result':result}


def analysis_tree_snapshot(path):
    """Read-only content identity; symlinks are recorded and never followed."""
    path=Path(path);records=[];links=[];total=allocated=0;empty=[]
    for base,dirs,names in os.walk(path,followlinks=False):
        if not dirs and not names:empty.append(str(Path(base)))
        for name in sorted(dirs+names):
            p=Path(base)/name
            if p.is_symlink():links.append(str(p));records.append((str(p.relative_to(path)),'SYMLINK',os.readlink(p)))
            elif p.is_file():
                st=p.stat();total+=st.st_size;allocated+=st.st_blocks*512;records.append((str(p.relative_to(path)),st.st_size,sha256(p)))
    records.sort();digest=hashlib.sha256(json.dumps(records,ensure_ascii=False).encode()).hexdigest()
    return {'size_bytes':total,'allocated_file_bytes':allocated,'n_files':sum(r[1]!='SYMLINK' for r in records),
            'tree_sha256':digest,'symlinks':links,'empty_leaf_directories':empty,'files':records}


def cleanup_analysis_allowlist(repo,out,allowlist,protected_paths,tracked_paths):
    """Delete only validated explicit targets, after a saved dry run; never delete a run root."""
    repo=Path(repo).resolve();base=repo/'analysis_runs';out=Path(out);out.mkdir(parents=True,exist_ok=True)
    protected=[Path(p).resolve() for p in protected_paths];tracked=[(repo/p).resolve() for p in tracked_paths];dry=[]
    for item in allowlist:
        relative=Path(item['relative_path']);candidate=repo/relative;reason=None
        if relative.is_absolute() or '..' in relative.parts or candidate==base or base not in candidate.resolve().parents:reason='NOT_STRICTLY_INSIDE_ANALYSIS_RUNS'
        elif not candidate.is_dir():reason='NOT_AN_EXISTING_DIRECTORY'
        elif any(p.is_symlink() for p in [candidate,*list(candidate.parents)[:len(candidate.parts)-len(base.parts)]]):reason='SYMLINK_PATH'
        elif any(candidate.resolve()==p or candidate.resolve() in p.parents or p.is_dir() and p in candidate.resolve().parents for p in protected):reason='PROTECTED_PATH_OVERLAP'
        elif any(candidate.resolve()==p or candidate.resolve() in p.parents for p in tracked):reason='TRACKED_CONTENT'
        snapshot=analysis_tree_snapshot(candidate) if reason is None else None
        if snapshot:
            if snapshot['symlinks']:reason='SYMLINK_CONTENT'
            elif snapshot['tree_sha256']!=item['expected_tree_sha256']:reason='CONTENT_HASH_MISMATCH'
            elif any(Path(r[0]).name=='airflow_vt.dat' or Path(r[0]).suffix.lower() in ('.dat','.vtu','.wav','.msh','.bin') for r in snapshot['files']):reason='PRIMARY_DATA_PRESENT'
            elif item['kind']=='EMPTY_BUILD_CACHE':
                if snapshot['n_files'] or any(candidate.iterdir()) or 'CMakeFiles' not in candidate.parts or candidate.name not in ('CMakeScratch','pkgRedirects','tmp','Temporary'):reason='NOT_A_DISPOSABLE_EMPTY_BUILD_CACHE'
            elif item['kind']=='OWNED_DERIVED_DEVELOPMENT':
                if candidate.name!='development' or not candidate.parent.name.startswith('airflow_global_cycle_validation_') or not item.get('replacement_verified'):reason='DEVELOPMENT_REPLACEMENT_NOT_PROVEN'
                else:
                    for file in candidate.rglob('*.csv'):
                        with file.open(newline='') as f:
                            reader=csv.DictReader(f);fields=set(reader.fieldnames or []);rows=list(reader)
                        if {'case_id','phase','event_time_s','confidence'}.issubset(fields) and rows or file.name=='manual_validation.csv' and any(r.get('status') not in ('PENDING_MANUAL_ANNOTATION',None) for r in rows):reason='HUMAN_ANNOTATION_PRESENT';break
            else:reason='UNKNOWN_DELETION_KIND'
        dry.append(dict(item,status='KEEP' if reason else 'SAFE_TO_DELETE',safety_reason=reason,
                        realpath=str(candidate.resolve()),size_bytes=snapshot['size_bytes'] if snapshot else None,
                        n_files=snapshot['n_files'] if snapshot else None,tree_sha256=snapshot['tree_sha256'] if snapshot else None))
    write_csv(out/'cleanup_allowlist_dry_run.csv',dry)
    (out/'cleanup_dry_run.md').write_text('# Explicit cleanup dry run\n\nOnly SAFE_TO_DELETE entries in cleanup_allowlist_dry_run.csv are eligible. Primary data, tracked/protected content, symlinks, unknown purposes and changed hashes are blocked. Empty build caches use rmdir; owned derived development outputs require replacement proof and no human annotations. No wildcard or run-root deletion is used.\n\n'+'\n'.join(f'- {r["status"]}: {r["relative_path"]}: {r.get("safety_reason") or r.get("reason")}' for r in dry)+'\n')
    deleted=[];kept=[];before=shutil.disk_usage(base)
    for row in dry:
        if row['status']!='SAFE_TO_DELETE':kept.append(row);continue
        p=repo/row['relative_path'];now=analysis_tree_snapshot(p)
        if now['tree_sha256']!=row['tree_sha256'] or now['symlinks']:kept.append(dict(row,safety_reason='CHANGED_AFTER_DRY_RUN'));continue
        allocated=now['allocated_file_bytes']+p.stat().st_blocks*512
        if row['kind']=='EMPTY_BUILD_CACHE':p.rmdir()
        else:shutil.rmtree(p)
        deleted.append(dict(row,deleted_utc=datetime.now(timezone.utc).isoformat(),released_allocated_bytes_estimate=allocated,verified_absent=not p.exists()))
    after=shutil.disk_usage(base);write_csv(out/'cleanup_deleted.csv',deleted,['relative_path','size_bytes','tree_sha256','reason','kind','deleted_utc','released_allocated_bytes_estimate','verified_absent']);write_csv(out/'cleanup_kept_allowlist.csv',kept)
    dump_json(out/'cleanup_space.json',{'before_free_bytes':before.free,'after_free_bytes':after.free,'observed_free_change_bytes':after.free-before.free,
        'deleted_logical_bytes':sum(r['size_bytes'] for r in deleted),'released_allocated_bytes_estimate':sum(r['released_allocated_bytes_estimate'] for r in deleted),
        'note':'filesystem free-space difference can include concurrent activity; estimates are separate from logical file bytes','deleted_directories':len(deleted)})
    return deleted,kept


def parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--input',type=Path,help='explicit two-column step/Ug DAT; required except --self-test')
    p.add_argument('--dt',type=float,help='solver timestep [s]; no inferred default')
    p.add_argument('--manifest',type=Path,help='explicit manifest containing dt_s; must agree with --dt')
    p.add_argument('--start',type=float,help='analysis start [s], inclusive; default first sample')
    p.add_argument('--end',type=float,help='analysis end [s], exclusive; default last sample plus sample interval')
    p.add_argument('--out-dir',type=Path,help='results directory; default analysis_runs/airflow_metrics/<stem>_<UTC>')
    p.add_argument('--prominences',type=float,nargs='+',default=TRACKING_CONFIG['prominence_sweep_percent'],help='threshold sweep [%% of P95−P5]')
    p.add_argument('--primary-prominence',type=float,default=10,help='fixed primary prominence [%%]')
    p.add_argument('--fmin',type=float,default=60,help='ACF/event pitch floor [Hz]; CPP remains section-11 60–400 Hz')
    p.add_argument('--fmax',type=float,default=400,help='ACF/event pitch ceiling [Hz]')
    p.add_argument('--min-cycles',type=int,default=10,help='research QC minimum periods, separate from formula minima 2/3/5')
    p.add_argument('--zoom-start',type=float,help='QC plot start [s]; default center of main segment')
    p.add_argument('--zoom-duration',type=float,default=.050,help='QC zoom duration [s]')
    p.add_argument('--cpp-frame-length',type=float,default=.060,help='CPP frame duration [s], must fit fixed nfft2048 at 10kHz')
    p.add_argument('--cpp-hop',type=float,default=.005,help='CPP hop duration [s]')
    p.add_argument('--reference-fullrate',type=Path,help='explicit actual stride-1 DAT of same waveform; no synthetic upsampling')
    p.add_argument('--manual-events',type=Path,help='manual reference CSV: case_id,phase,event_time_s,confidence,reviewer,notes; legacy time_s is also supported')
    p.add_argument('--overwrite',action='store_true',help='replace only this tool\'s previously generated files; protect source inputs')
    p.add_argument('--self-test',action='store_true',help='run embedded fixed-seed synthetic and CLI tests without existing runs')
    p.add_argument('--diagnose-tracking',action='store_true',help='add separate interruption, common-window, segment-pooling and one-factor threshold diagnostics; strict main outputs unchanged')
    p.add_argument('--diagnostic-baseline',type=Path,help='optional existing summary.json whose input, window, settings, event IDs and metrics must be reproduced')
    p.add_argument('--diagnostic-context-periods',type=float,default=5,help='periods before AND after each interruption shown in diagnostic plots')
    p.add_argument('--reanchor-policy',choices=['legacy','recover-rejected'],default='legacy',help='legacy reproduces previous tracking; recovery revisits rejected candidates with explicit conflict/safety checks')
    p.add_argument('--case-id',help='manual reference case ID; otherwise derived from archived run name or input stem')
    p.add_argument('--export-manual-review',action='store_true',help='export full-window overlapping blind/diagnostic panels and empty manual CSV template')
    p.add_argument('--manual-review-start',type=float,help='optional review panel start [s], inside analysis window')
    p.add_argument('--manual-review-end',type=float,help='optional review panel end [s], inside analysis window')
    p.add_argument('--tracking-method',choices=['existing','global-cycle','closure-threshold-and-area','all-events','acf-guided-events','peak-associated-closure'],default='existing',help='legacy, independent all-events, or main-peak-associated closure events')
    p.add_argument('--comparison-reference',type=Path,help='peak-associated-closure: saved per-case acf-guided comparison.csv; require exact old metrics and counts')
    p.add_argument('--cpp-reference',type=Path,help='all-events: saved case_id/flow_CPP_dB reference CSV; CPP is never recomputed')
    p.add_argument('--closure-thresholds',default='0.5,2,5',help='comma-separated local Upeak percentages, with zero baseline')
    p.add_argument('--area-input',type=Path,help='explicit area DAT; default verified sibling area.dat when present')
    p.add_argument('--area-column',default='interior-min',help='interior-min of section columns (solver convention), all-min, or explicit zero-based Amin column')
    p.add_argument('--area-time-column',type=int,default=0,help='zero-based step column; time is step times dt')
    p.add_argument('--area-unit',choices=['auto','mm2','m2'],default='auto',help='input area unit; auto requires a recognized header')
    p.add_argument('--area-event-method',choices=['global-cycle-extrema'],default='global-cycle-extrema')
    p.add_argument('--threshold-profile',default='balanced-v1',help='global-cycle fixed prior: balanced-v1 or explicit selected_profile.json')
    p.add_argument('--skip-threshold-sweep',action='store_true',help='global-cycle: run the chosen profile only; full shared grid is the default')
    return p


def report_text(summary):
    m=summary['primary']['main'];p=summary['primary'];c=summary['cpp'];meta=summary['metadata'];sens=summary['sensitivity_summary']
    def fmt(v):return 'undefined' if v is None else f'{v:.8g}'
    lines=['# Airflow metrics and reliability','',f"Input: `{Path(meta['input_path']).name}`; SHA256 `{meta['input_sha256']}`.",
           f"Solver dt={meta['dt_s']:.8g} s; step stride={meta['step_stride']}; sample rate={meta['sample_rate_hz']:.8g} Hz.",
           f"Window [{meta['start_s']:.8g}, {meta['end_s']:.8g}) s; {meta['window_samples']} samples; Ug range [{meta['window_min_m3_s']:.8g}, {meta['window_max_m3_s']:.8g}] m³/s.",
           '', '## Four metrics', '',
           f"Jitter local={fmt(m['jitter_local_percent'])}%; RAP={fmt(m['rap_percent'])}%; PPQ5={fmt(m['ppq5_percent'])}%.",
           f"Event QC: **{p['status']}**; flags: {', '.join(p['flags']) or 'none'}.",
           f"These values use the same events in the longest accepted contiguous segment (ID {m['segment_id']}, {m['n_cycles']} periods, {m['duration_s']:.8g} s); they are not whole-window jitter. No period trimming or gap bridging is performed.",
           f"Mean period denominator={fmt(m['mean_period_s'])} s; period CV={fmt(m['period_cv_percent'])}%; median F0={fmt(m['median_f0_hz'])} Hz.",
           f"Control B (different minima phase), same requested window: Jitter={fmt(summary['control']['main']['jitter_local_percent'])}%, RAP={fmt(summary['control']['main']['rap_percent'])}%, PPQ5={fmt(summary['control']['main']['ppq5_percent'])}%; {summary['control']['main']['n_cycles']} periods. Phase-dependent differences are not corrected by averaging A and B; inspect event_alignment.csv.",
           f"A-to-A cycles with exactly one B event: {sum(r['n_B_events']==1 for r in summary['event_alignment'])}/{len(summary['event_alignment'])}; this is phase correspondence, not ground-truth recall.",
           f"flow_CPP_dB={fmt(c['flow_CPP_dB'])} dB; status **{c['status']}**; valid frames {c['n_valid_frames']}/{c['n_frames']}.",
           f"CPP median={fmt(c.get('median_CPP_dB'))} dB; std={fmt(c.get('std_CPP_dB'))} dB. CPP flags: {', '.join(c.get('flags',[])) or 'none'}.",
           *[f"- CPP warning frame {r['frame_id']}, center {r['center_s']:.6g} s: peak={fmt(r['peak_frequency_hz'])} Hz, CPP={fmt(r['CPP_frame_dB'])} dB, {r['flags']}. Value retained per fixed specification." for r in c['rows'] if r['flags']],
           '', '## Period candidates and tracking', '',
           f"ACF proposed period={fmt(summary['pitch']['period_s'])} s; selection: {summary['pitch'].get('selection_rule','no candidate')}. Alternatives and all evidence are in acf_candidates.csv and summary.json.",
           f"Expected periods proxy=(window duration)/T_ACF={fmt(p['expected_periods_proxy'])}; detected completed periods across separate segments={p['n_complete_periods_all_segments']}.",
           f"coverage_proxy={fmt(p['coverage_proxy'])}; covered_time_fraction={fmt(p['covered_time_fraction'])} (sum of completed segment spans / requested duration). These are not true detection rates.",
           f"Median accepted non-anchor shape correlation={fmt(p['median_match_correlation'])}; ambiguous fraction={fmt(p['ambiguous_event_fraction'])} ({p['ambiguous_event_count']}/{p['ambiguity_denominator_candidates']} selectable in-window candidates). Missing-candidate breaks are separately logged.",
           f"ACF/event relative median-period difference={fmt(p['acf_consistency_relative_difference'])}. Pitch flags: {summary['pitch']['flags']}.",
           '', '## Threshold sensitivity', '',
           f"Primary prominence fixed before analysis: {summary['tracking_config']['primary_prominence_percent']}%. All sweep values: {summary['tracking_config']['prominence_sweep_percent']}%.",
           f"Eligible structural-QC thresholds={sens['n_valid_thresholds']}; excluded={sens['n_excluded_thresholds']}. Jitter median/min/max/IQR/CV: {fmt(sens['median'])}, {fmt(sens['min'])}, {fmt(sens['max'])}, {fmt(sens['IQR'])}, {fmt(sens['CV_percent'])}.",
           'This is method sensitivity, not a confidence interval. All provisional values, event counts, correspondence (one raw sample tolerance), coverage and QC are retained in threshold_sensitivity.csv.', '',
           '## Subsample timing', '',
           f"Raw timestamp resolution={p['timestamp_resolution_us']:.8g} µs; median samples per period={fmt(p['median_samples_per_period'])}.",
           'OFF and quadratic ON use the same event IDs and the same segment; peaks are not reselected. Interpolation does not prove closeness to true event times. Fallbacks and offsets are saved.',
           *[f"- {r['subsample']}: Jitter={fmt(r['jitter_local_percent'])}%, RAP={fmt(r['rap_percent'])}%, PPQ5={fmt(r['ppq5_percent'])}%; max offset={fmt(r['offset_max_abs_us'])} µs." for r in summary['subsample_comparison']],
           '', '## Sampling interval', '']
    if summary['sampling_comparison'][0]['status']=='NOT_AVAILABLE':lines.append('NOT_AVAILABLE: no actual stride-1 reference. Stride-1 data cannot be recovered from stride-5 data; no fake comparison was made.')
    else:
        lines.append('Actual stride-1 and five stride-5 decimation phases were compared, with identical CPP frame seconds and 10 kHz CPP bandwidth. Values at common saved steps were verified before comparison.')
        for row in summary['sampling_comparison']:lines.append(f"- stride={row['step_stride']}, phase={row['phase']}, {row['subsample']}: Jitter={fmt(row['jitter_local_percent'])}%, CPP={fmt(row['flow_CPP_dB'])} dB, QC={row['jitter_qc']}.")
    lines+=['', '## Visual checks and limitations', '',
            'Inspect flow_events_zoom.png: raw candidates, accepted/rejected peaks, event numbers, actual periods and a five-period PPQ5 example. Markers are at real candidate times, never snapped to T0. Unselected peaks may be real secondary lobes, not necessarily noise.',
            'Inspect cepstrum_example.png: the measured 10log10(c²), one fixed Theil–Sen background, residual-peak candidates, adopted q and vertical CPP difference. CPP is independent of event tracking and is not CPPS.',
            'CPP definition: '+CPP_DEFINITION['frame_cpp']+'. '+CPP_DEFINITION['aggregate']+'. Reference algorithm settings are fully recorded in summary.json.',
            'Harmonic alternatives do not prove period doubling. Symmetry does not guarantee zero jitter. Detector ambiguity, true temporal variation and sample quantization must be distinguished. QC_PASS is not proof of truth; clinical acoustic thresholds do not apply to simulated airflow.',
            '', '## Warnings', '', *['- '+w for w in meta['warnings']]]
    return '\n'.join(lines)+'\n'


def analyze_cli(args):
    if args.tracking_method=='peak-associated-closure':return peak_associated_cli(args)
    if args.tracking_method=='acf-guided-events':return acf_guided_cli(args)
    if args.tracking_method=='all-events':return all_events_cli(args)
    if args.tracking_method=='closure-threshold-and-area':return closure_area_cli(args)
    if args.tracking_method=='global-cycle':return global_cycle_cli(args)
    if args.input is None:raise ValueError('--input is required except with --self-test')
    dt=resolve_dt(args.dt,args.manifest);data=read_flow(args.input,dt)
    if not np.isfinite([args.fmin,args.fmax]).all() or not 0<args.fmin<args.fmax:raise ValueError('require 0 < fmin < fmax')
    if args.min_cycles<1 or not args.prominences or any(not np.isfinite(v) or v<=0 for v in args.prominences+[args.primary_prominence]):raise ValueError('positive finite prominences and min-cycles are required')
    if not np.isfinite([args.cpp_frame_length,args.cpp_hop,args.zoom_duration]).all() or min(args.cpp_frame_length,args.cpp_hop,args.zoom_duration)<=0:raise ValueError('positive finite frame/hop/zoom durations are required')
    start=data['t'][0] if args.start is None else args.start;end=data['t'][-1]+data['sampling_interval_s'] if args.end is None else args.end
    if not np.isfinite([start,end]).all() or start<0 or start>=end:raise ValueError('invalid half-open time window')
    tol=data['sampling_interval_s']*1e-6
    if start<data['t'][0]-tol or end>data['t'][-1]+data['sampling_interval_s']+tol:raise ValueError('window outside saved data support')
    if args.zoom_start is not None and (not np.isfinite(args.zoom_start) or args.zoom_start<start or args.zoom_start>=end):raise ValueError('--zoom-start must be inside analysis window')
    if any(v is not None and not np.isfinite(v) for v in (args.manual_review_start,args.manual_review_end)):raise ValueError('manual review bounds must be finite')
    if (args.manual_review_start is not None or args.manual_review_end is not None) and not args.export_manual_review:raise ValueError('manual review bounds require --export-manual-review')
    case_id=args.case_id or next((part.split('_P')[0] for part in Path(args.input).parts if part.startswith('r') and '_P' in part),Path(args.input).stem)
    reference=read_flow(args.reference_fullrate,dt) if args.reference_fullrate else None
    out=(args.out_dir or Path('analysis_runs/airflow_metrics')/(Path(args.input).stem+'_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))).resolve()
    if args.diagnostic_baseline and not (args.diagnose_tracking or args.reanchor_policy=='recover-rejected'):raise ValueError('--diagnostic-baseline requires --diagnose-tracking or recovery policy')
    if not np.isfinite(args.diagnostic_context_periods) or args.diagnostic_context_periods<=0:raise ValueError('positive finite diagnostic context required')
    protected=[Path(args.input).resolve()]+[Path(p).resolve() for p in (args.manifest,args.reference_fullrate,args.manual_events,args.diagnostic_baseline) if p]
    if any(p==out or out in p.parents for p in protected):raise ValueError('output directory contains a protected input; choose a separate output directory')
    if out.exists():
        if args.reanchor_policy=='recover-rejected' or args.export_manual_review:raise ValueError('recovery/review requires a new --out-dir; existing results are protected')
        if args.diagnose_tracking:raise ValueError('diagnostic mode requires a new --out-dir; existing analysis results are protected even with --overwrite')
        if not args.overwrite:raise ValueError('output directory already exists; use --overwrite only for generated results')
        if (out/'tracking_diagnostics').exists():raise ValueError('existing diagnostic results are append-protected; choose a new --out-dir')
        marker=out/'analysis_metadata.json'
        if not marker.is_file() or json.loads(marker.read_text()).get('tool_id')!=TOOL_ID:raise ValueError('refusing to overwrite an unrecognized output directory')
    tracking=dict(TRACKING_CONFIG,pitch_min_hz=args.fmin,pitch_max_hz=args.fmax,
                  primary_prominence_percent=args.primary_prominence,prominence_sweep_percent=list(args.prominences))
    if args.reanchor_policy!='legacy':tracking['reanchor_policy']=args.reanchor_policy
    cppconfig=dict(CPP_CONFIG,frame_sec=args.cpp_frame_length,hop_sec=args.cpp_hop)
    result=analyze_arrays(data['t'],data['x'],data['fs'],start,end,tracking,cppconfig,args.min_cycles)
    result['_cpp_config']=cppconfig
    sampling=sampling_comparison(data,reference,start,end,tracking,cppconfig,args.min_cycles,result)
    inside=(data['t']>=start)&(data['t']<end);selected=data['x'][inside]
    warnings=[]
    if args.start is None or args.end is None:warnings.append('Window not fully explicit; initial static interval or pressure ramp may be included. Use a stationary window for primary evaluation.')
    if args.start is None or args.end is None:print('WARNING: '+warnings[0],file=sys.stderr)
    if data['fs']/args.fmax<10:warnings.append('Few raw samples per highest requested period; sub-sample timing is weakly constrained.')
    static=bool(np.any(data['x']==0));warnings.extend(result['pitch']['flags'])
    meta={'tool_id':TOOL_ID,'created_utc':datetime.now(timezone.utc).isoformat(),'input_path':data['path'],'input_sha256':data['sha256'],
          'case_id':case_id,'tracker_version':TRACKER_VERSION if args.reanchor_policy!='legacy' else 'section11_legacy',
          'reanchor_policy':args.reanchor_policy,'settings_sha256':hashlib.sha256(json.dumps({'tracking':tracking,'cpp':cppconfig,'qc':dict(QC_CONFIG,min_cycles=args.min_cycles)},sort_keys=True).encode()).hexdigest(),
          'manifest_path':str(args.manifest.resolve()) if args.manifest else None,'manifest_sha256':sha256(args.manifest) if args.manifest else None,
          'script_sha256':sha256(__file__),'dt_s':dt,'step_stride':data['stride'],'step_stride_validation':'all intervals identical including last; no interpolation of missing samples',
          'sampling_interval_s':data['sampling_interval_s'],'sample_rate_hz':data['fs'],'first_step':int(data['steps'][0]),'last_step':int(data['steps'][-1]),
          'saved_time_span_s':float(data['t'][-1]-data['t'][0]),'saved_support_end_s':float(data['t'][-1]+data['sampling_interval_s']),
          'n_samples':len(data['x']),'start_s':start,'end_s':end,'window_semantics':'[start,end)', 'window_samples':int(np.count_nonzero(inside)),
          'window_min_m3_s':float(selected.min()),'window_max_m3_s':float(selected.max()),'input_min_m3_s':float(data['x'].min()),'input_max_m3_s':float(data['x'].max()),
          'has_exact_zero_samples':static,'initial_exact_zero_samples':int(np.argmax(data['x']!=0)) if np.any(data['x']!=0) else len(data['x']),
          'warnings':warnings,'reference_fullrate_sha256':reference['sha256'] if reference else None,
          'dependencies':{'numpy':np.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__}}
    summary={'tool_id':TOOL_ID,'metadata':meta,'cpp_config':cppconfig,'tracking_config':tracking,'qc_config':dict(QC_CONFIG,min_cycles=args.min_cycles),
             'CPP_definition':CPP_DEFINITION,'period_formulas':{'jitter':'100*mean(abs(diff(T)))/mean(T), N>=2','RAP':'100*mean(abs(T[1:-1]-3-cycle local mean))/mean(T), N>=3','PPQ5':'100*mean(abs(T[2:-2]-5-cycle local mean))/mean(T), N>=5',
             'primary_scope':'longest contiguous accepted segment','pooled_scope':'diagnostic only: within-segment numerators pooled; denominator mean of all within-segment periods; never bridge segment boundaries'},
             **{k:v for k,v in result.items() if not k.startswith('_')},'sampling_comparison':sampling}
    manual_rows=read_manual_reference(args.manual_events,case_id,data,start,end) if args.manual_events else []
    manual_stats=[];manual_matches=[];manual_period_checks=[];manual_metric_periods=[]
    current_tracks={method:result['_tracks'][result[private]] for method,private in [('A','_primary_key'),('B','_control_key')]}
    if manual_rows:
        manual_stats,manual_matches,manual_period_checks,manual_metric_periods=evaluate_manual_reference(manual_rows,current_tracks,start,end)
        summary['manual_validation_status']='MANUAL_REFERENCE_EVALUATED'
        for phase,key in [('A','primary'),('B','control')]:
            if any(row['phase']==phase for row in manual_rows):
                summary[key]['manual_validation_status']='MANUAL_REFERENCE_EVALUATED'
                summary[key]['manual_validation_scope']='statistics within supplied annotated phase/time extent; not certification of the entire requested window'
        summary['manual_validation']={'status':'MANUAL_REFERENCE_EVALUATED','reference_sha256':sha256(args.manual_events),'rows':manual_stats}
        main_stat=next((r for r in manual_stats if r['phase']=='A' and r['scope']=='high_only' and np.isclose(r['tolerance_us'],100)),manual_stats[0])
        summary['manual_event_validation']={'precision':main_stat['precision'],'recall':main_stat['recall'],'F1':main_stat['F1'],'matched':main_stat['TP'],'scope':'phase-specific annotation extent, 100 us primary'}
    else:summary['manual_validation_status']='PENDING_MANUAL_ANNOTATION'
    out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.airflow_metrics_',dir=out.parent) as name:
        temp=Path(name)
        dump_json(temp/'analysis_metadata.json',meta);dump_json(temp/'summary.json',summary)
        m=result['primary']['main'];c=result['cpp']
        row={k:m[k] for k in ('jitter_local_percent','rap_percent','ppq5_percent','n_cycles','mean_period_s','period_cv_percent')}
        row.update(flow_CPP_dB=c['flow_CPP_dB'],jitter_status=result['primary']['status'],cpp_status=c['status'],start_s=start,end_s=end,
                   main_segment_id=m['segment_id'],main_segment_duration_s=m['duration_s'],scope=result['primary']['main_scope'],cpp_valid_frames=c['n_valid_frames'],coverage_proxy=result['primary']['coverage_proxy'],covered_time_fraction=result['primary']['covered_time_fraction'])
        row.update(research_primary_candidate='pooled_within_segments_v2',manual_validation_status=summary['manual_validation_status'],
                   **{'pooled_v2_'+name:detail['value'] for name,detail in result['primary']['pooled_within_segments_v2']['metric_details'].items()},
                   **{'pooled_v2_'+name+'_'+key:detail[key] for name,detail in result['primary']['pooled_within_segments_v2']['metric_details'].items() for key in ('mean_period_s','n_cycles','n_terms','covered_time_fraction')})
        write_csv(temp/'summary.csv',[row]);write_csv(temp/'cpp_frames.csv',c['rows'],list(c['rows'][0]) if c['rows'] else ['frame_id','CPP_frame_dB','status'])
        write_csv(temp/'acf_candidates.csv',result['pitch']['candidates'],['lag_samples','period_s','frequency_hz','acf','spectrum_support'])
        write_csv(temp/'local_acf.csv',result['local_acf'],['time_s','period_s','f0_hz','flags'])
        write_csv(temp/'threshold_sensitivity.csv',result['sensitivity']);write_csv(temp/'subsample_comparison.csv',result['subsample_comparison']);write_csv(temp/'sampling_comparison.csv',sampling)
        write_csv(temp/'event_alignment.csv',result['event_alignment'],['A_cycle_id','segment_id','A_start_s','A_end_s','n_B_events','B_times_s','status'])
        records=[];periods=[];segments=[];matches=[];breaks=[]
        for key in (result['_primary_key'],result['_control_key']):
            tr=result['_tracks'][key];summ=result['_summaries'][key]
            for e in tr['events']:records.append(dict(e,selected_main=e['event_id'] in summ['main']['event_ids']))
            periods.extend(summ['_periods']);segments.extend(dict(r,method=tr['method']) for r in summ['segments']);matches.extend(dict(r,method=tr['method']) for r in tr['attempts']);breaks.extend(dict(r,method=tr['method']) for r in tr['breaks'])
        write_csv(temp/'events.csv',records,['event_id','method','sample_id','raw_time_s','subsample_time_s','offset_us','interpolation_status','prominence','normalized_prominence','selectable','accepted','segment_id','reason','match_score','correlation','anchor','in_window','selected_main'])
        write_csv(temp/'candidate_events.csv',records,['event_id','method','sample_id','raw_time_s','prominence','selectable','accepted','reason','segment_id'])
        write_csv(temp/'periods.csv',periods,['method','segment_id','period_id','event_start','event_end','start_s','end_s','period_s','delta_previous_period_s','flags','missing_event_possible','selected_main','used_for_jitter_v2','used_for_RAP_v2','used_for_PPQ5_v2'])
        write_csv(temp/'segments.csv',segments,['method','segment_id','start_s','end_s','duration_s','n_cycles','mean_period_s','jitter_local_percent','rap_percent','ppq5_percent','event_ids'])
        write_csv(temp/'candidate_matches.csv',matches,['method','segment_id','direction','previous_event','candidate_event','template_event','period_reference_s','delta_s','score','rho','rank','margin','accepted'])
        write_csv(temp/'tracking_breaks.csv',breaks,['method','segment_id','direction','previous_event','time_s','reason','expanded_diagnostic_candidates','period_reference_s'])
        if 'manual_event_validation' in summary:dump_json(temp/'manual_event_validation.json',summary['manual_event_validation'])
        if args.reanchor_policy=='recover-rejected' or args.export_manual_review or args.manual_events:
            write_csv(temp/'manual_validation.csv',manual_stats or [{'status':'PENDING_MANUAL_ANNOTATION','case_id':case_id}])
            write_csv(temp/'manual_event_matches.csv',manual_matches);write_csv(temp/'manual_period_checks.csv',manual_period_checks)
            write_csv(temp/'manual_metric_periods.csv',manual_metric_periods)
            write_csv(temp/'manual_reference_validated.csv',manual_rows,['case_id','phase','event_time_s','confidence','reviewer','notes','phase_extremum_error_s','input_format'])
            (temp/'manual_validation_report.md').write_text('# Manual event validation\n\n'+summary['manual_validation_status']+'\n\n'+('No human event annotations have been supplied. No real-data precision/recall or phase correctness has been established. Use the blind review panels and fill the empty CSV template.' if not manual_rows else 'See manual_validation.csv for high-only/all annotations at 50/100/150 us; manual_event_matches.csv for TP/FP/FN; manual_period_checks.csv for skipped-period candidates. Metrics are limited to the annotated time extent. A partial annotation is not whole-window validation.')+'\n')
        save_plots(temp,data,start,end,result,args.zoom_start,args.zoom_duration,tracking,sampling)
        (temp/'report.md').write_text(report_text(summary))
        if args.diagnose_tracking:
            diagnostic_summary=diagnostic_analysis(temp/'tracking_diagnostics',data,start,end,result,tracking,args.min_cycles,
                                                  args.diagnostic_baseline if args.reanchor_policy=='legacy' else None,args.diagnostic_context_periods)
            dump_json(temp/'tracking_diagnostics'/'diagnostic_summary.json',diagnostic_summary)
        if args.reanchor_policy=='recover-rejected':
            validation=reanchor_validation_outputs(temp,data,start,end,result,tracking,args.min_cycles,case_id,args.diagnostic_baseline)
            validation['manual_validation_status']=summary['manual_validation_status'];summary['reanchor_validation']=validation
            dump_json(temp/'summary.json',summary)
            with (temp/'report.md').open('a') as f:f.write('\n## Reanchor validation and research aggregation\n\nLegacy fields main/summary.csv jitter retain the longest-segment definition. The new research candidate pooled_within_segments_v2 uses only segments that contribute terms to each metric, with metric-specific normalization and coverage. See jitter_comparison.csv, metric_terms_v2.csv and periods.csv. New connections satisfy conservative automatic gates but require manual phase validation. Manual status: '+summary['manual_validation_status']+'.\n')
        if args.export_manual_review:export_manual_review(temp,data,start,end,result,case_id,args.manual_review_start,args.manual_review_end)
        if sha256(data['path'])!=data['sha256']:raise ValueError('input changed during analysis; refusing to publish results')
        if not out.exists():temp.rename(out)
        else:
            if args.reanchor_policy=='recover-rejected' or args.export_manual_review:raise ValueError('recovery/review output appeared concurrently; choose a new directory')
            if args.diagnose_tracking:raise ValueError('diagnostic results require a new --out-dir; refusing to overwrite concurrent output')
            if not (temp/'manual_event_validation.json').exists() and (out/'manual_event_validation.json').is_file():
                (out/'manual_event_validation.json').unlink()
            for file in temp.iterdir():os.replace(file,out/file.name)
    print(f"Jitter local [%]: {m['jitter_local_percent']}; RAP [%]: {m['rap_percent']}; PPQ5 [%]: {m['ppq5_percent']}; {result['primary']['status']}")
    print(f"flow_CPP_dB: {c['flow_CPP_dB']}; {c['status']}; valid frames {c['n_valid_frames']}/{c['n_frames']}")
    print(f"Main segment periods: {m['n_cycles']}; covered time fraction (proxy): {result['primary']['covered_time_fraction']:.6g}")
    if args.reanchor_policy=='recover-rejected':
        v2=result['primary']['pooled_within_segments_v2'];d=v2['metric_details']['jitter_local_percent']
        print(f"Research candidate pooled_v2: Jitter={v2['jitter_local_percent']}; RAP={v2['rap_percent']}; PPQ5={v2['ppq5_percent']}; Jitter cycles={d['n_cycles']}, terms={d['n_terms']}, time fraction={d['covered_time_fraction']:.6g}; {summary['manual_validation_status']}")
    print('Saved:',out)
    return summary


def self_test():
    """Fixed-seed tests live here to keep implementation to one source file."""
    class Tests(unittest.TestCase):
        @staticmethod
        def pulses(fs=20000,duration=.6,jitter=False,double=False,missing=False):
            t=np.arange(round(fs*duration))/fs;truth=[.012];i=0
            while truth[-1]<duration-.015:
                truth.append(truth[-1]+(.008*(1+.025*np.sin(i*.91)) if jitter else .008));i+=1
            x=np.zeros_like(t)
            for i,at in enumerate(truth):
                if missing and i==len(truth)//2:continue
                x+=np.exp(-.5*((t-at)/.00065)**2)
                if double:x+=.35*np.exp(-.5*((t-at-.0025)/.0004)**2)
            return t,x,np.array(truth)

        @staticmethod
        def closure_signal(jitter=False,amplitude=False,missing=False,fs=20000,duration=.6):
            t=np.arange(round(duration*fs))/fs;x=np.zeros_like(t);peaks=[];sample=int(.012*fs);i=0
            while sample/ fs<duration-.012:
                at=sample/fs;amp=1+.3*np.sin(i*.7) if amplitude else 1.;rise=.0015;fall=.0025
                pulse=np.where(t<at,1+(t-at)/rise,1-(t-at)/fall);x=np.maximum(x,amp*np.maximum(pulse,0));peaks.append(sample)
                sample+=round(.008*fs)+(round(5*np.sin(i*.9)) if jitter else 0);i+=1
            x[np.abs(x)<1e-12]=0.  # Exact synthetic plateaus, not floating-point ramp residue.
            if missing:
                middle=peaks[len(peaks)//2];x[max(0,middle-round(.0035*fs)):middle+round(.008*fs)]=np.maximum(x[max(0,middle-round(.0035*fs)):middle+round(.008*fs)],.10)
            dt=1e-5;stride=round(1/fs/dt)
            return {'t':t,'x':x,'fs':fs,'steps':np.arange(len(t))*stride,'stride':stride,'dt_s':dt,'sampling_interval_s':1/fs},peaks

        @staticmethod
        def closure_pipeline(data,start=.1,end=.5):
            cfg=dict(GLOBAL_CYCLE_CONFIG);base=global_cycle_track(data,start,end,cfg,alternatives=False)
            pulses,branches,candidates,edges=closure_pulse_candidates(data,start,end,base,cfg,[.5,2.,5.])
            tracks,rows=closure_threshold_tracks(data,pulses,branches,[.5,2.,5.],base['frames'],start,end,cfg)
            return base,pulses,tracks,rows,candidates

        def test_closure_regular_zero_plateau_is_not_event(self):
            data,peaks=self.closure_signal();base,pulses,tracks,rows,candidates=self.closure_pipeline(data)
            for name,tr in tracks.items():
                stat=closure_statistics(tr,.1,.5);self.assertGreater(stat['row']['n_periods'],40);self.assertLess(stat['row']['jitter_local_percent'],1e-7)
                pct=next(e['threshold_pct'] for e in tr['events']);expected=np.array(peaks)/data['fs']+.0025*(1-pct/100)
                actual=np.array([e['subsample_time_s'] for e in tr['events'] if e['accepted']])
                self.assertLess(np.max(np.min(abs(actual[:,None]-expected),axis=1)),1e-12)
                self.assertTrue(all(e['Ug_left_m3_s']>e['threshold_m3_s']>=e['Ug_right_m3_s'] for e in tr['events'] if e['accepted']))
                self.assertGreater(sum(e['confidence_class']=='HIGH_CONFIDENCE' for e in tr['events'] if e['accepted']),30)

        def test_closure_known_jitter_and_local_amplitude(self):
            data,peaks=self.closure_signal(jitter=True,amplitude=True);_,_,tracks,_,_=self.closure_pipeline(data)
            for tr in tracks.values():
                es=[e for e in tr['events'] if e['accepted']];pct=es[0]['threshold_pct'];truth=np.array(peaks)/data['fs']+.0025*(1-pct/100);truth=truth[(truth>=.1)&(truth<.5)]
                self.assertEqual(len(es),len(truth));np.testing.assert_allclose([e['subsample_time_s'] for e in es],truth,atol=1e-12)
                known=period_metrics(np.diff(truth));actual=closure_statistics(tr,.1,.5)['row']
                for key in METRIC_NAMES:self.assertAlmostEqual(actual[key],known[key],places=8)

        def test_closure_missing_crossing_splits_then_resumes(self):
            data,peaks=self.closure_signal(missing=True);base,pulses,tracks,rows,_=self.closure_pipeline(data)
            stats={n:closure_statistics(tr,.1,.5) for n,tr in tracks.items()};intervals={n:closure_classify(st,data,.1,.5,base['frames']) for n,st in stats.items()}
            for tr in tracks.values():
                missing=[e for e in tr['events'] if e['reason']=='NO_THRESHOLD_CROSSING'];self.assertTrue(missing)
                self.assertTrue(any(e['accepted'] and e['subsample_time_s']>missing[0]['peak_time_s']+.012 for e in tr['events']))
                self.assertGreater(len(tr['segments']),1)
            self.assertEqual(closure_validate_statistics(stats,data,None,rows,intervals,.1,.5)['status'],'PASS')

        def test_closure_secondary_hump_is_not_pulse_end(self):
            data,peaks=self.closure_signal()
            for index in peaks:data['x']+=.32*np.exp(-.5*((data['t']-(index/data['fs']+.0012))/.00018)**2)
            _,_,tracks,_,_=self.closure_pipeline(data)
            for tr in tracks.values():
                es=[e for e in tr['events'] if e['accepted']];self.assertGreater(len(es),45)
                self.assertLess(closure_statistics(tr,.1,.5)['row']['jitter_local_percent'],.05)
                self.assertTrue(all(e['Upeak_m3_s']>.99 for e in es))

        def test_closure_bounce_first_stable_and_no_fabrication(self):
            data={'t':np.arange(12)/10000,'x':np.array([1,.01,.04,.01,.03,0,0,0,0,0,0,0]),'fs':10000,'steps':np.arange(12)}
            rows=closure_crossing_candidates(data,0,8,2.,.008);chosen=[r for r in rows if r['selected']]
            self.assertEqual(len(chosen),1);self.assertEqual(chosen[0]['sample_left'],4);self.assertGreater(len(rows),1)
            data['x']=np.ones(12)*.1;self.assertEqual(closure_crossing_candidates(data,0,8,.5,.008),[])
            data['x']=np.zeros(12);self.assertEqual(closure_crossing_candidates(data,0,8,5,.008),[])

        def test_closure_area_columns_units_and_step_qc(self):
            with tempfile.TemporaryDirectory() as directory:
                p=Path(directory)/'area.dat';p.write_text('# step area[mm^2]\n0 0 2 3 nan\n5 0 3 4 nan\n10 0 4 5 nan\n')
                d=read_minimum_area(p,1e-5);np.testing.assert_array_equal(d['x'],[2,3,4]);self.assertEqual(d['stride'],5)
                d=read_minimum_area(p,1e-5,'all-min');np.testing.assert_array_equal(d['x'],[0,0,0])
                p.write_text('# step area[m^2]\n0 0 0.000002 0.000003 nan\n10 0 0.000003 0.000004 nan\n20 0 0.000004 0.000005 nan\n')
                d=read_minimum_area(p,1e-5);np.testing.assert_allclose(d['x'],[2,3,4]);self.assertEqual(d['stride'],10)
                with self.assertRaises(ValueError):read_minimum_area(p,1e-5,unit='mm2')
                p.write_text('# step area[mm^2]\n0 0 nan nan 0\n5 0 2 3 0\n10 0 3 4 0\n')
                with self.assertRaises(ValueError):read_minimum_area(p,1e-5)
                p.write_text('# step area[mm^2]\n0 0 -1 3 0\n5 0 2 3 0\n10 0 3 4 0\n')
                with self.assertRaises(ValueError):read_minimum_area(p,1e-5)
                p.write_text('# step area[mm^2]\n0 0 1 3 0\n5 0 2 3 0\n11 0 3 4 0\n')
                with self.assertRaises(ValueError):read_minimum_area(p,1e-5)

        def test_closure_area_independent_extrema_plateau_and_stride(self):
            flow,_=self.closure_signal();area,_=self.closure_signal(fs=10000)
            rr=global_cycle_track(area,.1,.5,dict(GLOBAL_CYCLE_CONFIG),alternatives=False)
            for phase in ('A','B'):
                stat=closure_statistics(rr['methods'][phase]['track'],.1,.5);self.assertGreater(stat['row']['n_events'],40)
                self.assertLess(stat['row']['jitter_local_percent'],1e-6)
            minima=[e for e in rr['methods']['B']['track']['events'] if e['accepted']];self.assertTrue(all(e['plateau_width_s']>0 for e in minima))
            self.assertEqual(area['stride'],2*flow['stride'])
            times=[]
            for key in ('plateau_start_s','plateau_center_s','plateau_end_s'):
                s=closure_statistics(rr['methods']['B']['track'],.1,.5,representative=key);self.assertLess(s['row']['jitter_local_percent'],1e-6);times.append(s['track']['events'][s['track']['segments'][0]['event_indices'][0]]['subsample_time_s'])
            self.assertLess(times[0],times[1]);self.assertLess(times[1],times[2])
            for event in rr['methods']['A']['track']['events']:event['plateau_width_s']=0
            original=closure_statistics(rr['methods']['A']['track'],.1,.5)
            unchanged=closure_statistics(rr['methods']['A']['track'],.1,.5,representative='plateau_start_s')
            self.assertEqual(original['periods'],unchanged['periods'])

        def test_closure_confirmed_filter_never_reconnects_rejected_event(self):
            data,_=self.closure_signal();base,_,tracks,_,_=self.closure_pipeline(data);tr=next(iter(tracks.values()))
            inside=[e for e in tr['events'] if e['accepted']]
            for e in inside:e['confidence_class']='HIGH_CONFIDENCE'
            mid=inside[len(inside)//2];mid['confidence_class']='ESTIMATED'
            stat=closure_statistics(tr,.1,.5,confirmed=True)
            self.assertTrue(all(not p['start_s']<mid['subsample_time_s']<p['end_s'] for p in stat['periods']))
            self.assertGreater(stat['row']['n_segments'],1)

        def test_closure_cli_export_and_manual_compatibility(self):
            data,peaks=self.closure_signal(duration=.3)
            with tempfile.TemporaryDirectory() as directory,contextlib.redirect_stdout(io.StringIO()):
                root=Path(directory);flow=root/'airflow_vt.dat';ar=root/'area.dat'
                np.savetxt(flow,np.c_[data['steps'],data['x']],header='step airflow[m^3/s]')
                np.savetxt(ar,np.c_[data['steps'],np.zeros(len(data['x'])),data['x']*20,data['x']*30,np.zeros(len(data['x']))],header='step area[mm^2]')
                before={str(p):sha256(p) for p in (flow,ar)}
                argv=['--input',str(flow),'--dt','1e-5','--start','.1','--end','.25','--case-id','synthetic','--tracking-method','closure-threshold-and-area','--skip-threshold-sweep','--out-dir',str(root/'out')]
                args=parser().parse_args(argv);s=closure_area_cli(args)
                self.assertEqual(s['manual_status'],'PENDING_MANUAL_ANNOTATION');self.assertEqual(s['metadata']['area_status'],'AVAILABLE')
                for p in ('flow_threshold_events.csv','min_area_candidates.csv','min_area_events_A_B.csv','min_area_plateaus.csv','periods_by_definition.csv','metric_terms_by_definition.csv','qc_intervals.csv','flow_threshold_crossings_full.png','min_area_extrema_full.png','jitter_comparison.png','manual_events_template.csv'):
                    self.assertTrue((root/'out'/p).is_file(),p)
                manifest=json.loads((root/'out/run_manifest.json').read_text());self.assertTrue(all(sha256(root/'out'/p)==h for p,h in manifest['files'].items()))
                self.assertEqual(before,{str(p):sha256(p) for p in (flow,ar)})
                with self.assertRaises(ValueError):closure_area_cli(args)
                # Actual hand-defined linear-ramp crossing times, independent of detector outputs.
                refs=root/'manual.csv';rows=[{'case_id':'synthetic','signal':'Ug','event_type':'downcross','threshold_pct':.5,'event_time_s':i/data['fs']+.0025*.995,'confidence':'high','reviewer':'synthetic truth','notes':'analytic triangle'} for i in peaks if .1<=i/data['fs']+.0025*.995<.25]
                write_csv(refs,rows)
                fake_tracks={}
                for pct in (.5,2,5):
                    name='Ug_downcross_'+format(pct,'g').replace('.','p')+'pct';times=[i/data['fs']+.0025*(1-pct/100) for i in peaks if .1<=i/data['fs']+.0025*(1-pct/100)<.25]
                    es=[{'accepted':True,'event_id':str(i),'subsample_time_s':at} for i,at in enumerate(times)]
                    fake_tracks[name]={'track':{'events':es,'segments':[{'segment_id':1,'event_indices':list(range(len(es)))}]}}
                validated,ss,mm,cc,cp=closure_manual_reference(refs,'synthetic',data,None,fake_tracks,.1,.25)
                self.assertTrue(validated);self.assertTrue(all(r['F1']==1 for r in ss));self.assertTrue(cc)

        def test_closure_missing_area_and_unmatched_area_render(self):
            from unittest.mock import patch
            data,_=self.closure_signal(duration=.3)
            with tempfile.TemporaryDirectory() as directory,contextlib.redirect_stdout(io.StringIO()):
                root=Path(directory);flow=root/'airflow_vt.dat';np.savetxt(flow,np.c_[data['steps'],data['x']],header='step airflow[m^3/s]')
                args=parser().parse_args(['--input',str(flow),'--dt','1e-5','--start','.1','--end','.25','--tracking-method','closure-threshold-and-area','--skip-threshold-sweep','--out-dir',str(root/'out')])
                with patch.dict(closure_area_cli.__globals__,closure_figures=lambda *a:None):result=closure_area_cli(args)
                self.assertEqual(result['metadata']['area_status'],'MISSING_INPUT');self.assertEqual(len(result['methods']),5)
                self.assertTrue(all(result['methods'][n]['statistics']['classification_coverage']==1 for n in result['methods']))
                es=[{'accepted':True,'event_id':'B_'+str(i),'subsample_time_s':at,'plateau_width_s':0} for i,at in enumerate((.12,.128,.136))]
                area_stat={'track':{'events':es,'segments':[{'segment_id':1,'event_indices':[0,1,2]}]}}
                flow_stat={'track':{'events':[{'accepted':True,'event_id':'c','pulse_id':'p','subsample_time_s':.121,'Upeak_m3_s':1.}], 'segments':[{'segment_id':1,'event_indices':[0]}]}}
                rows,_=closure_flow_area_matches({'Ug_downcross_0p5pct':flow_stat,'Amin_minima':area_stat},[{'period_s':.008}])
                unmatched=[r for r in rows if r['status']=='UNMATCHED_AREA'];self.assertEqual(len(unmatched),2);self.assertTrue(all(r['flow_time_s'] is None for r in unmatched))

        @staticmethod
        def all_event_data(values,fs=1000):
            x=np.array(values,dtype=float);dt=1e-5;stride=round(1/fs/dt)
            return {'x':x,'t':np.arange(len(x))/fs,'fs':fs,'steps':np.arange(len(x))*stride,'dt_s':dt,'stride':stride,'sampling_interval_s':1/fs}

        def test_all_event_long_short_intervals_and_every_term(self):
            times=[.501,.502,.5021,.510,.511,.710,.7102,.999]
            events=[{'event_time_s':at,'event_id':str(i)} for i,at in enumerate(times)];m,p,t=all_event_metrics(events,'synthetic',.5,1.)
            self.assertEqual(m['n_intervals'],7);self.assertEqual(m['n_jitter_terms'],6);self.assertEqual(m['n_RAP_terms'],5);self.assertEqual(m['n_PPQ5_terms'],3)
            self.assertEqual(len(p),7);self.assertEqual(m['intervals_dropped'],0);self.assertGreater(max(e['interval_s'] for e in p),.2)
            mean=sum(np.diff(times))/7
            for width,name in ((2,METRIC_NAMES[0]),(3,METRIC_NAMES[1]),(5,METRIC_NAMES[2])):
                expected=[];v=np.diff(times)
                for i in range(len(v)-width+1):expected.append(abs(v[i+1]-v[i]) if width==2 else abs(v[i+width//2]-sum(v[i:i+width])/width))
                self.assertAlmostEqual(m[name],100*sum(expected)/len(expected)/mean)
                self.assertEqual(len([r for r in t if r['metric']==name]),len(expected))

        def test_all_event_multipeaks_are_not_merged(self):
            d=self.all_event_data([0,1,.3,.8,0,1,.3,.8,0,1,.3,.8,0]);series,_,_,counts=all_event_extract(d,None,0,.013)
            self.assertEqual(len(series['Ug_max']),6);self.assertEqual(len(series['Ug_min']),5)
            self.assertTrue(all(r['ACF_excluded']==r['shape_excluded']==r['length_excluded']==0 for r in counts))
            m,periods,_=all_event_metrics(series['Ug_max'],'Ug_max',0,.013);self.assertEqual(m['n_intervals'],5)

        def test_all_event_zero_plateau_and_missing_zero_gap(self):
            d=self.all_event_data([1,0,0,0,1,.2,1,.2,1,0,0,1,0,0]);series,_,_,_=all_event_extract(d,None,0,.014)
            self.assertEqual([e['sample_id'] for e in series['Ug_zero']],[1,9,12]);m,periods,terms=all_event_metrics(series['Ug_zero'],'Ug_zero',0,.014)
            self.assertAlmostEqual(periods[0]['interval_s'],.008);self.assertEqual(m['n_intervals'],2);self.assertEqual(m['n_jitter_terms'],1);self.assertIsNone(m['rap_percent'])
            self.assertEqual(len(terms),1);self.assertEqual(m['segments_used'],1)

        def test_all_event_thresholds_local_branches_and_no_hold_gate(self):
            d=self.all_event_data([0,1,.05,.8,.01,.7,.2,1,0,0,1,.05,1,0]);series,branches,crossings,_=all_event_extract(d,None,0,.014)
            self.assertGreater(len(series['Ug_10pct']),len(series['Ug_3pct']))
            for pct in (3.,10.):
                rr=series['Ug_'+format(pct,'g')+'pct'];self.assertEqual(len({(r['sample_left'],r['sample_right']) for r in rr}),len(rr))
                for e in rr:
                    k=e['sample_left'];level=e['Upeak_m3_s']*pct/100;self.assertGreater(d['x'][k],level);self.assertGreaterEqual(level,d['x'][k+1])
                    self.assertAlmostEqual(e['event_time_s'],d['t'][k]+(level-d['x'][k])/(d['x'][k+1]-d['x'][k])/d['fs'])
            self.assertTrue(any(b['reason']=='NO_THRESHOLD_CROSSING' for b in branches));self.assertTrue(all(not r['stability_gate_applied'] for r in crossings))
            self.assertTrue(any(e['sample_left']==3 for e in series['Ug_3pct']))  # Only one low sample before a rebound: still included.

        def test_all_event_windows_plateau_and_short_sequence(self):
            d=self.all_event_data([1,0,0,0,1,0,1,0,1]);series,_,_,_=all_event_extract(d,None,.001,.007)
            self.assertEqual([e['event_time_s'] for e in series['Ug_zero']],[.001,.005]);self.assertTrue(all(.001<=e['event_time_s']<.007 for es in series.values() for e in es))
            minima=all_event_extrema(d,'Ug')['min'];self.assertEqual(minima[0]['event_time_s'],.002)
            for n in (0,1,2):
                m,_,_=all_event_metrics([{'event_time_s':.1+i*.01,'event_id':str(i)} for i in range(n)],'s',0,1)
                self.assertIsNone(m['jitter_local_percent']);self.assertEqual(m['n_intervals'],max(0,n-1))

        def test_all_event_cli_reuses_inputs_but_no_tracker_or_CPP(self):
            from unittest.mock import patch
            data,_=self.closure_signal(duration=.3)
            with tempfile.TemporaryDirectory() as directory,contextlib.redirect_stdout(io.StringIO()):
                root=Path(directory);flow=root/'airflow_vt.dat';area=root/'area.dat'
                np.savetxt(flow,np.c_[data['steps'],data['x']],header='step airflow[m^3/s]');np.savetxt(area,np.c_[data['steps'],np.zeros(len(data['x'])),data['x'],data['x']*2,np.zeros(len(data['x']))],header='step area[mm^2]')
                original={str(p):sha256(p) for p in (flow,area)}
                args=parser().parse_args(['--input',str(flow),'--dt','1e-5','--start','.1','--end','.25','--tracking-method','all-events','--closure-thresholds','3,10','--out-dir',str(root/'out')])
                def forbidden(*a,**k):raise AssertionError('legacy tracker/ACF/CPP must not run')
                with patch.dict(all_events_cli.__globals__,global_cycle_track=forbidden,pooled_metrics_v2=forbidden,pitch_candidates=forbidden,cpp_analysis=forbidden):result=all_events_cli(args)
                self.assertEqual(len(result['rows']),7);self.assertTrue(all(r['all_detected_events_used']=='YES' for r in result['rows']))
                es=list(csv.DictReader((root/'out/events.csv').open()));pp=list(csv.DictReader((root/'out/periods.csv').open()));tt=list(csv.DictReader((root/'out/jitter_contributions.csv').open()))
                for row in result['rows']:
                    times=[float(e['event_time_s']) for e in es if e['method']==row['method']];v=np.diff(times)
                    self.assertEqual(len(v),row['n_intervals']);self.assertEqual(sum(e['method']==row['method'] for e in pp),len(v));self.assertEqual(sum(e['method']==row['method'] for e in tt),max(0,len(v)-1))
                    self.assertEqual(period_metrics(v)['jitter_local_percent'],row['jitter_local_percent'])
                self.assertEqual(original,{str(p):sha256(p) for p in (flow,area)})
                self.assertEqual(len(list((root/'out/figures').glob('*.png'))),6)
                with self.assertRaises(ValueError):all_events_cli(args)

        def test_acf_guided_zero_one_many_regions_and_boundaries(self):
            from unittest.mock import patch
            d=self.all_event_data(np.zeros(200),fs=1000)
            for i in (102,108,115):d['x'][i]=1.
            frames=[{'center_s':v,'period_s':.01,'acf_confidence':.9,'flags':[]} for v in (.1,.14)]
            with patch.dict(acf_guided_extract.__globals__,global_local_periods=lambda *a:(frames,[])):
                result=acf_guided_extract(d,None,.1,.14)
            regions=[r for r in result['regions'] if r['signal']=='Ug'];counts=[r['n_max_accepted'] for r in regions]
            self.assertIn(0,counts);self.assertIn(1,counts);self.assertIn(2,counts)
            self.assertEqual(sum(counts),len(result['series']['Ug_max']));self.assertEqual(len(result['series']['Ug_max']),3)
            self.assertEqual(acf_region_owner(.11,[{'start_s':.1,'end_s':.11},{'start_s':.11,'end_s':.12}]),1)
            self.assertIsNone(acf_region_owner(.12,[{'start_s':.1,'end_s':.11},{'start_s':.11,'end_s':.12}]))

        def test_acf_guided_small_bump_rejected_clear_double_kept(self):
            t,x,truth=self.pulses(duration=.6,double=True)
            x+=sum(.02*np.exp(-.5*((t-at-.0045)/.00015)**2) for at in truth)
            d=self.all_event_data(x,fs=20000);result=acf_guided_extract(d,None,.1,.5)
            es=result['series']['Ug_max'];primary=sum(np.min(abs(e['event_time_s']-truth))<.0002 for e in es);secondary=sum(np.min(abs(e['event_time_s']-(truth+.0025)))<.0002 for e in es)
            self.assertGreater(primary,40);self.assertGreater(secondary,40)
            self.assertFalse(any(np.min(abs(e['event_time_s']-(truth+.0045)))<.0002 for e in es))
            self.assertTrue(any(not e['clarity_accepted'] and 'BELOW_LOCAL_PROMINENCE' in e['detector_rejection_reasons'] for e in result['peak_candidates']))

        def test_acf_guided_low_ACF_confidence_never_rejects_clear_peak(self):
            from unittest.mock import patch
            d=self.all_event_data(np.zeros(200),fs=1000)
            for i in (102,115,128):d['x'][i-1:i+2]=[.5,1.,.5]
            frames=[{'center_s':v,'period_s':None,'acf_confidence':0.,'flags':['PITCH_AMBIGUOUS']} for v in (.1,.14)]
            with patch.dict(acf_guided_extract.__globals__,global_local_periods=lambda *a:(frames,[])):result=acf_guided_extract(d,None,.1,.14)
            self.assertEqual(len(result['series']['Ug_max']),3);self.assertTrue(all(not e['period_supported'] for e in result['series']['Ug_max']))
            self.assertTrue(all(r['scale_source']=='FALLBACK_1_FMIN_SEARCH_ONLY' for r in result['regions']))

        def test_acf_guided_zero_independent_and_long_intervals_preserved(self):
            data,_=self.closure_signal(missing=True);normal=acf_guided_extract(data,None,.1,.5)
            strict=acf_guided_extract(data,None,.1,.5,dict(GLOBAL_CYCLE_CONFIG,prominence_percent=95.,excursion_percent=95.))
            self.assertEqual(normal['series']['Ug_zero'],strict['series']['Ug_zero'])
            self.assertEqual(normal['series']['Ug_zero'],all_event_extract(data,None,.1,.5)[0]['Ug_zero'])
            m,periods,terms=all_event_metrics(normal['series']['Ug_zero'],'Ug_zero',.1,.5)
            self.assertGreater(max(p['interval_s'] for p in periods),.012);self.assertEqual(len(periods),m['n_events']-1)
            self.assertEqual(len([t for t in terms if t['metric']=='jitter_local_percent']),len(periods)-1)

        def test_acf_guided_crossings_use_clear_peaks_and_real_brackets(self):
            data,_=self.closure_signal();result=acf_guided_extract(data,None,.1,.5);clear={e['event_id'] for e in result['signal_details']['Ug']['clear']['max']}
            for name in ('Ug_3pct','Ug_10pct'):
                es=result['series'][name];self.assertGreater(len(es),40)
                for e in es:
                    self.assertIn(e['peak_id'],clear);k=e['sample_left'];level=e['threshold_m3_s'];self.assertGreater(data['x'][k],level);self.assertGreaterEqual(level,data['x'][k+1])
                    expected=data['t'][k]+(level-data['x'][k])/(data['x'][k+1]-data['x'][k])/data['fs'];self.assertAlmostEqual(expected,e['event_time_s'],places=12)

        def test_acf_guided_cli_no_DP_shape_CPP_or_segment_gates(self):
            from unittest.mock import patch
            data,_=self.closure_signal(duration=.3)
            with tempfile.TemporaryDirectory() as directory,contextlib.redirect_stdout(io.StringIO()):
                root=Path(directory);flow=root/'airflow_vt.dat';area=root/'area.dat'
                np.savetxt(flow,np.c_[data['steps'],data['x']],header='step airflow[m^3/s]');np.savetxt(area,np.c_[data['steps'],np.zeros(len(data['x'])),data['x'],data['x']*2,np.zeros(len(data['x']))],header='step area[mm^2]')
                before={str(p):sha256(p) for p in (flow,area)}
                args=parser().parse_args(['--input',str(flow),'--dt','1e-5','--start','.1','--end','.25','--case-id','synthetic','--tracking-method','acf-guided-events','--out-dir',str(root/'out')])
                def forbidden(*a,**k):raise AssertionError('DP/shape-path/CPP/segment gate must not be called')
                with patch.dict(acf_guided_cli.__globals__,global_graph=forbidden,global_dp=forbidden,global_cycle_track=forbidden,pooled_metrics_v2=forbidden,cpp_analysis=forbidden):result=acf_guided_cli(args)
                self.assertEqual(result['rows'][0]['method'],'Ug_zero');self.assertEqual(len(result['rows']),7)
                self.assertEqual(len(list((root/'out/figures').glob('*.png'))),5)
                es=list(csv.DictReader((root/'out/events.csv').open()));pp=list(csv.DictReader((root/'out/periods.csv').open()));jj=list(csv.DictReader((root/'out/jitter_contributions.csv').open()))
                for row in result['rows']:
                    times=[float(e['event_time_s']) for e in es if e['method']==row['method']];v=np.diff(times)
                    self.assertEqual(len(v),row['n_intervals']);self.assertEqual(sum(p['method']==row['method'] for p in pp),len(v));self.assertEqual(sum(p['method']==row['method'] for p in jj),max(0,len(v)-1))
                    self.assertEqual(period_metrics(v)['jitter_local_percent'],row['jitter_local_percent'])
                self.assertEqual(before,{str(p):sha256(p) for p in (flow,area)})
                with self.assertRaises(ValueError):acf_guided_cli(args)

        @staticmethod
        def associated_fixture(values,peak_ids,region_bounds,prominences=None):
            data=Tests.all_event_data(values);regions=[];clear=[]
            for j,(left,right) in enumerate(zip(region_bounds[:-1],region_bounds[1:])):
                regions.append({'region_id':j,'start_s':left,'end_s':right,'search_scale_s':right-left,'predicted_T0_s':right-left,'signal':'Ug'})
            for k,i in enumerate(peak_ids):
                at=float(data['t'][i]);clear.append({'event_id':'Ug_max_'+str(i),'sample_id':i,'event_time_s':at,'raw_time_s':at,
                    'sample_value':float(data['x'][i]),'prominence':prominences[k] if prominences is not None else float(data['x'][i]),
                    'search_region_id':acf_region_owner(at,regions)})
            series=all_event_extract(data,None,region_bounds[0],region_bounds[-1])[0];series['Ug_max']=clear
            base={'series':series,'regions':regions}
            return data,base

        def test_peak_associated_main_secondary_empty_and_tie(self):
            data,base=self.associated_fixture([0,1,0,.6,0,1,0,0,0,0],(1,3,5),(0,.005,.007,.01),(1.,.8,1.))
            peaks,cycles=peak_associated_select(base['series']['Ug_max'],base['regions'])
            self.assertEqual([p['sample_id'] for p in peaks if p['peak_role']=='MAIN'],[1,5])
            self.assertEqual([p['sample_id'] for p in peaks if p['peak_role']=='SECONDARY'],[3])
            self.assertEqual(cycles[2]['peak_status'],'NO_CLEAR_PEAK_NO_FILL')
            self.assertEqual(cycles[0]['secondary_peak_count'],1)
            tied=[dict(e,prominence=1.) for e in base['series']['Ug_max']]
            self.assertEqual(peak_associated_select(tied,base['regions'])[0][0]['peak_role'],'MAIN')

        def test_peak_associated_first_zero_terminal_crossings_and_history(self):
            data,base=self.associated_fixture([0,1,.08,.02,.6,.09,.02,0,.2,0,1,.5,0,0,0],(1,4,8,10),(0,.01,.015))
            result=peak_associated_extract(data,0,.015,base);series=result['series']
            self.assertEqual([e['sample_id'] for e in series['main_zero']],[7,12])
            self.assertEqual([e['sample_id'] for e in series['including_secondary_zero']],[7,9,12])
            details=[d for d in result['reach_history'] if d['peak_id']=='Ug_max_1' and d['threshold_name']=='zero'][0]
            self.assertEqual(details['reach_count'],2);self.assertTrue(details['rises_above_after_first']);self.assertTrue(details['repeated_downcross_after_rise'])
            for suffix in ('3pct','10pct'):
                e=series['main_'+suffix][0];first=series['first_crossing_'+suffix][0]
                self.assertGreater(e['event_time_s'],first['event_time_s']);self.assertLess(e['event_time_s'],.007)
                k=e['sample_left'];level=e['threshold_m3_s'];self.assertGreater(data['x'][k],level);self.assertGreaterEqual(level,data['x'][k+1])
                self.assertAlmostEqual(e['event_time_s'],data['t'][k]+(level-data['x'][k])/(data['x'][k+1]-data['x'][k])*.001)
            self.assertTrue(all(not d['order_violations'] for d in result['order_diagnostics']))
            secondary=[d for d in result['reach_history'] if d['peak_id']=='Ug_max_4']
            self.assertTrue(all(d['auxiliary_adoption']=='CLOSURE_EPISODE_ALREADY_OWNED' for d in secondary))
            owned={}
            for name,es in series.items():
                if not name.startswith('including_secondary_'):continue
                for e in es:
                    run=e['descent_run_id'];self.assertIn(owned.get(run,e['peak_id']),(e['peak_id'],));owned[run]=e['peak_id']

        def test_peak_associated_nonclosure_independent_thresholds_and_edges(self):
            data,base=self.associated_fixture([.2,1,.08,.05,.2,1,.7,.2,0,1],(1,5),(0,.005,.008))
            result=peak_associated_extract(data,0,.008,base)
            self.assertEqual(result['series']['main_zero'],[]);self.assertEqual(result['series']['main_3pct'],[])
            self.assertEqual(len(result['series']['main_10pct']),1)
            self.assertTrue(all(not r['zero_reached'] for r in result['cycles']))
            self.assertTrue(all(e['event_time_s']<.008 for es in result['series'].values() for e in es))

        def test_peak_associated_missing_cycles_all_intervals_and_terms(self):
            x=np.full(150,.2);ids=list(range(1,142,10));x[ids]=1.
            for i in (3,13,43,53,113,123,143):x[i]=0.
            data,base=self.associated_fixture(x,ids,np.arange(0,.151,.01))
            result=peak_associated_extract(data,0,.15,base);es=result['series']['main_zero']
            m,pp,tt=all_event_metrics(es,'main_zero',0,.15)
            self.assertEqual(m['n_events'],7);self.assertEqual(m['n_intervals'],6);self.assertAlmostEqual(m['max_interval_s'],.06)
            self.assertEqual([sum(t['metric']==k for t in tt) for k in METRIC_NAMES],[5,4,2])
            self.assertEqual(m['intervals_dropped'],0);self.assertEqual(m['intervals_split'],0)
            self.assertTrue(any(not c['zero_reached'] for c in result['cycles']))
            v=np.diff([e['event_time_s'] for e in es]);self.assertAlmostEqual(m[METRIC_NAMES[0]],100*np.mean(abs(np.diff(v)))/np.mean(v))

        def test_peak_associated_cli_exact_baseline_and_no_CPP_or_DP(self):
            from unittest.mock import patch
            data,_=self.closure_signal(duration=.3)
            with tempfile.TemporaryDirectory() as directory,contextlib.redirect_stdout(io.StringIO()):
                root=Path(directory);flow=root/'airflow_vt.dat';np.savetxt(flow,np.c_[data['steps'],data['x']],header='step airflow[m^3/s]')
                original=sha256(flow);common=['--input',str(flow),'--dt','1e-5','--start','.1','--end','.25','--case-id','synthetic']
                old=acf_guided_cli(parser().parse_args(common+['--tracking-method','acf-guided-events','--out-dir',str(root/'old')]))
                old_hash=sha256(root/'old/comparison.csv')
                args=parser().parse_args(common+['--tracking-method','peak-associated-closure','--comparison-reference',str(root/'old/comparison.csv'),'--out-dir',str(root/'new')])
                def forbidden(*a,**k):raise AssertionError('CPP/DP/segmentation must not run')
                with patch.dict(peak_associated_cli.__globals__,cpp_analysis=forbidden,global_graph=forbidden,global_dp=forbidden,global_cycle_track=forbidden,pooled_metrics_v2=forbidden):result=peak_associated_cli(args)
                self.assertEqual(sha256(flow),original);self.assertEqual(sha256(root/'old/comparison.csv'),old_hash)
                self.assertEqual(len(list((root/'new/figures').glob('*.png'))),4)
                es=list(csv.DictReader((root/'new/events.csv').open()));pp=list(csv.DictReader((root/'new/periods.csv').open()));tt=list(csv.DictReader((root/'new/metric_terms.csv').open()))
                for row in result['rows']:
                    v=np.diff([float(e['event_time_s']) for e in es if e['method']==row['method']])
                    self.assertEqual(len(v),row['n_intervals']);self.assertEqual(sum(p['method']==row['method'] for p in pp),len(v))
                    self.assertEqual(sum(p['method']==row['method'] and p['metric']==METRIC_NAMES[0] for p in tt),max(0,len(v)-1))
                    self.assertEqual(period_metrics(v)[METRIC_NAMES[0]],row[METRIC_NAMES[0]])
                self.assertEqual([r for r in result['rows'] if r['method']=='reference_Ug_zero'][0][METRIC_NAMES[0]],old['rows'][0][METRIC_NAMES[0]])
                with self.assertRaises(ValueError):peak_associated_cli(args)

        def test_period_formulas(self):
            periods=np.array([.008,.0081,.0079,.0082,.008,.0084,.0078]);m=period_metrics(periods)
            mean=sum(periods)/len(periods)
            local=100*sum(abs(periods[i+1]-periods[i]) for i in range(len(periods)-1))/(len(periods)-1)/mean
            rap=100*sum(abs(periods[i]-sum(periods[i-1:i+2])/3) for i in range(1,len(periods)-1))/(len(periods)-2)/mean
            ppq=100*sum(abs(periods[i]-sum(periods[i-2:i+3])/5) for i in range(2,len(periods)-2))/(len(periods)-4)/mean
            self.assertAlmostEqual(m['jitter_local_percent'],local);self.assertAlmostEqual(m['rap_percent'],rap);self.assertAlmostEqual(m['ppq5_percent'],ppq)
            self.assertIsNone(period_metrics([1])['jitter_local_percent']);self.assertIsNone(period_metrics([1,1])['rap_percent']);self.assertIsNone(period_metrics([1]*4)['ppq5_percent'])

        def test_regular_and_subsample_identity(self):
            t,x,truth=self.pulses();pitch=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            self.assertAlmostEqual(pitch['period_s'],.008,delta=1/20000)
            track=extract_track(t,x,20000,.1,.5,pitch['period_s'],10,TRACKING_CONFIG)
            on=summarize_track(track,.1,.5,pitch['period_s'],10,20000);off=summarize_track(track,.1,.5,pitch['period_s'],10,20000,False,on['_membership'])
            self.assertGreater(on['main']['n_cycles'],40);self.assertLess(on['main']['jitter_local_percent'],1e-6)
            self.assertEqual(on['main']['event_ids'],off['main']['event_ids'])
            times=np.array([e['subsample_time_s'] for e in track['events'] if e['event_id'] in on['main']['event_ids']])
            self.assertLess(np.max(np.min(abs(times[:,None]-truth),axis=1)),1/20000)

        def test_known_jitter_not_suppressed(self):
            t,x,truth=self.pulses(jitter=True,double=True);p=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            tr=extract_track(t,x,20000,.1,.5,p['period_s'],10,TRACKING_CONFIG);r=summarize_track(tr,.1,.5,p['period_s'],10,20000)
            ev={e['event_id']:e for e in tr['events']};times=np.array([ev[k]['subsample_time_s'] for k in r['main']['event_ids']]);indices=np.argmin(abs(times[:,None]-truth),axis=1)
            self.assertGreater(len(indices),30);self.assertTrue(np.all(np.diff(indices)==1))
            actual=period_metrics(np.diff(truth[indices]));estimated=r['main']
            for key in ('jitter_local_percent','rap_percent','ppq5_percent'):self.assertAlmostEqual(estimated[key],actual[key],delta=.03)

        def test_double_peaks_not_double_cycles(self):
            t,x,_=self.pulses(double=True);p=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            tr=extract_track(t,x,20000,.1,.5,p['period_s'],10,TRACKING_CONFIG);r=summarize_track(tr,.1,.5,p['period_s'],10,20000)
            self.assertAlmostEqual(p['period_s'],.008,delta=.00005);self.assertGreater(r['main']['n_cycles'],40);self.assertLess(r['main']['n_cycles'],55)
            self.assertIn('DOUBLE_PEAKS',r['flags']);self.assertLess(r['main']['jitter_local_percent'],.01)
            self.assertGreater(len(tr['events']),sum(e['accepted'] for e in tr['events']))

        def test_missing_no_segment_bridge(self):
            t,x,truth=self.pulses(missing=True);p=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            tr=extract_track(t,x,20000,.1,.5,p['period_s'],10,TRACKING_CONFIG);r=summarize_track(tr,.1,.5,p['period_s'],10,20000)
            self.assertGreaterEqual(sum(s['n_cycles']>0 for s in r['segments']),2)
            self.assertIn('MISSING_EVENTS_POSSIBLE',r['flags'])
            gap=truth[len(truth)//2]
            self.assertFalse(any(p['start_s']<gap<p['end_s'] for p in r['_periods']))

        def test_period_doubling_candidates_retained(self):
            t,x,truth=self.pulses();x*=1+.45*np.cos(2*np.pi*(t-.012)/.016)
            p=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            self.assertTrue(any(abs(r['period_s']-.016)<.00015 for r in p['candidates']))
            self.assertEqual(len(p['alternatives']),3)
            self.assertIn('PITCH_AMBIGUOUS',p['flags'])
            # If ACF selects 2T, an actual T candidate is retained; no fabricated half-period events.
            tr=extract_track(t,x,20000,.1,.5,p['period_s'],10,TRACKING_CONFIG)
            for e in tr['events']:
                self.assertLessEqual(abs(e['subsample_time_s']-e['raw_time_s']),.5/20000+1e-12)

        def test_drift_modulation_shape_and_spurious_peak(self):
            t,x,truth=self.pulses(jitter=True,double=True);x=(1+.2*np.sin(2*np.pi*2*t))*x+.2*t
            x+=.4*np.exp(-.5*((t-.333)/.0001)**2)
            p=pitch_candidates(x[(t>=.1)&(t<.5)],20000,TRACKING_CONFIG)
            tr=extract_track(t,x,20000,.1,.5,p['period_s'],10,TRACKING_CONFIG);r=summarize_track(tr,.1,.5,p['period_s'],10,20000)
            self.assertGreater(len(tr['events']),30);self.assertTrue(all(e['sample_id']>=0 for e in tr['events']))
            self.assertTrue(r['flags'] or r['covered_time_fraction']>=.9)

        def test_initial_static_ramp_and_changing_shape(self):
            fs=20000;t=np.arange(12000)/fs;truth=np.arange(.108,.585,.008);x=np.zeros_like(t)
            for i,at in enumerate(truth):
                width=.0006*(1+.12*np.sin(i*.2))
                x+=np.exp(-.5*((t-at)/width)**2)
            x*=np.clip((t-.08)/.04,0,1)*(1+.15*np.sin(2*np.pi*2*t))
            x[t<.08]=0
            r=analyze_arrays(t,x,fs,.2,.5,TRACKING_CONFIG,CPP_CONFIG,10)
            self.assertGreater(r['primary']['main']['n_cycles'],25)
            self.assertLess(r['primary']['main']['jitter_local_percent'],.1)
            self.assertGreater(r['cpp']['n_valid_frames'],0)

        def test_cpp_definition_scale_and_noise(self):
            t,x,_=self.pulses(duration=.3);cfg=dict(CPP_CONFIG)
            regular=cpp_analysis(x,20000,0,cfg,.008);scaled=cpp_analysis(x*1e-5,20000,0,cfg,.008)
            noise=cpp_analysis(np.random.default_rng(123).normal(0,1,len(x)),20000,0,cfg,.008)
            self.assertGreater(regular['n_valid_frames'],10)
            self.assertAlmostEqual(regular['flow_CPP_dB'],scaled['flow_CPP_dB'],delta=1e-5)
            self.assertGreater(regular['flow_CPP_dB'],noise['flow_CPP_dB'])
            periods=[r['peak_q_s'] for r in regular['rows'] if r['peak_q_s']]
            self.assertAlmostEqual(np.median(periods),.008,delta=.0006)
            for ex in regular['_examples'].values():
                row=ex['row'];self.assertAlmostEqual(row['CPP_frame_dB'],row['cepstrum_at_peak_dB']-row['background_at_peak_dB'],places=10)
            self.assertEqual(cpp_analysis(x,5000,0,cfg)['status'],'INVALID_SAMPLING_RATE')
            self.assertEqual(cpp_analysis(np.zeros(2000),20000,0,cfg)['n_valid_frames'],0)
            self.assertEqual(cpp_analysis(x[:100],20000,0,cfg)['status'],'INSUFFICIENT_CPP_FRAMES')
            self.assertEqual(cpp_analysis(np.ones(2000)*.0001,20000,0,cfg)['n_valid_frames'],0)

        def test_interpolation_quantization_and_boundaries(self):
            fs=20000;t=np.arange(12000)/fs;x=np.sin(2*np.pi*137.3*t);p=pitch_candidates(x[(t>=.1)&(t<.5)],fs,TRACKING_CONFIG)
            tr=extract_track(t,x,fs,.1,.5,p['period_s'],10,TRACKING_CONFIG);on=summarize_track(tr,.1,.5,p['period_s'],10,fs);off=summarize_track(tr,.1,.5,p['period_s'],10,fs,False,on['_membership'])
            self.assertLess(on['main']['jitter_local_percent'],off['main']['jitter_local_percent'])
            self.assertTrue(all(.1<=r['start_s']<r['end_s']<.5 for r in on['_periods']))
            self.assertNotEqual(quadratic(np.array([1.,1.,1.]),1)[2],'QUADRATIC')

        def test_cpp_harmonics_epsilon_and_frame_sensitivity(self):
            fs=20000;t=np.arange(6000)/fs
            x=sum(np.cos(2*np.pi*125*k*t)/k for k in range(1,25))
            epsilon_values=[]
            for epsilon in (1e-10,1e-12,1e-14):
                cfg=dict(CPP_CONFIG,spectrum_epsilon_ratio=epsilon)
                c=cpp_analysis(x,fs,0,cfg,.008);epsilon_values.append(c['flow_CPP_dB'])
                self.assertTrue(np.isfinite(c['flow_CPP_dB']))
                self.assertAlmostEqual(np.median([r['peak_q_s'] for r in c['rows'] if r['peak_q_s']]),.008,delta=.0002)
            frame_values=[]
            for duration in (.05,.06,.08):
                c=cpp_analysis(x,fs,0,dict(CPP_CONFIG,frame_sec=duration),.008)
                frame_values.append(c['flow_CPP_dB']);self.assertGreater(c['n_valid_frames'],0)
            # These settings change the absolute result; do not assert software/clinical interchangeability.
            self.assertGreater(np.ptp(frame_values),.01)
            print('CPP synthetic epsilon sensitivity [dB]:',epsilon_values)
            print('CPP synthetic frame-length sensitivity [dB]:',frame_values)

        def test_input_and_dt_errors(self):
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);p=root/'with spaces.dat';p.write_text('# step airflow\n\n0 0\n5 1e-4\n10 2e-4\n')
                data=read_flow(p,1e-5);self.assertEqual(data['stride'],5);self.assertAlmostEqual(data['fs'],20000)
                for text in ('0 0\n5 nan\n10 0\n','0 0\n5 0\n5 0\n','0 0\n5 0\n11 0\n','0 0 0\n5 0 0\n10 0 0\n'):
                    p.write_text(text)
                    with self.assertRaises(ValueError):read_flow(p,1e-5)
                with self.assertRaises(ValueError):resolve_dt(None,None)
                manifest=root/'manifest';manifest.write_text('dt_s = 1e-5\n')
                self.assertEqual(resolve_dt(None,manifest),1e-5)
                with self.assertRaises(ValueError):resolve_dt(2e-5,manifest)

        def test_zero_noise_and_short_qc(self):
            fs=20000;t=np.arange(3000)/fs
            r=analyze_arrays(t,np.zeros(len(t)),fs,0,.15,TRACKING_CONFIG,CPP_CONFIG,10)
            self.assertIsNone(r['primary']['main']['jitter_local_percent']);self.assertEqual(r['pitch']['status'],'NO_OSCILLATION')
            noise=np.random.default_rng(8).normal(0,1,len(t))
            n=analyze_arrays(t,noise,fs,0,.15,TRACKING_CONFIG,CPP_CONFIG,10)
            self.assertNotEqual(n['primary']['status'],'QC_PASS')
            short=analyze_arrays(t[:100],np.sin(2*np.pi*125*t[:100]),fs,0,.005,TRACKING_CONFIG,CPP_CONFIG,10)
            self.assertNotEqual(short['primary']['status'],'QC_PASS')

        def test_fullrate_phases_and_mismatched_reference(self):
            fs=100000;t,x,_=self.pulses(fs=fs,duration=.16)
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);full=root/'full.dat';sparse=root/'sparse.dat'
                steps=np.arange(len(x));np.savetxt(full,np.c_[steps,x],header='step airflow');np.savetxt(sparse,np.c_[steps[::5],x[::5]],header='step airflow')
                d=read_flow(sparse,1/fs);ref=read_flow(full,1/fs);r=analyze_arrays(d['t'],d['x'],d['fs'],.02,.14,TRACKING_CONFIG,CPP_CONFIG,10)
                rows=sampling_comparison(d,ref,.02,.14,TRACKING_CONFIG,CPP_CONFIG,10,r)
                self.assertEqual(len(rows),12);self.assertEqual({r['phase'] for r in rows if r['phase'] is not None},{0,1,2,3,4})
                for row in rows:
                    self.assertAlmostEqual(row['mean_period_s'],.008,delta=.0001)
                self.assertEqual(sampling_comparison(d,None,.02,.14,TRACKING_CONFIG,CPP_CONFIG,10,r)[0]['status'],'NOT_AVAILABLE')
                ref['x']=ref['x']+1
                with self.assertRaises(ValueError):sampling_comparison(d,ref,.02,.14,TRACKING_CONFIG,CPP_CONFIG,10,r)

        def test_cli_artifacts_stdout_overwrite_and_protection(self):
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);path=root/'synthetic flow.dat';t,x,truth=self.pulses(duration=.30)
                np.savetxt(path,np.c_[np.arange(len(t))*5,x],header='step airflow[m^3/s]');out=root/'result'
                args=parser().parse_args(['--input',str(path),'--dt','1e-5','--start','.06','--end','.26','--out-dir',str(out)])
                original=sha256(path);buffer=io.StringIO()
                with contextlib.redirect_stdout(buffer):r=analyze_cli(args)
                self.assertIn('flow_CPP_dB',buffer.getvalue());self.assertIn('Jitter local',buffer.getvalue())
                for name in ('summary.json','analysis_metadata.json','summary.csv','events.csv','periods.csv','threshold_sensitivity.csv','subsample_comparison.csv','sampling_comparison.csv','flow_events_zoom.png','cepstrum_example.png','cpp_frames.csv','report.md'):
                    self.assertTrue((out/name).is_file(),name)
                with self.assertRaises(ValueError):analyze_cli(args)
                args.overwrite=True
                human=root/'human_events.csv'
                with human.open('w',newline='') as f:
                    w=csv.writer(f);w.writerow(['time_s']);w.writerows([[v] for v in truth])
                args.manual_events=human
                with contextlib.redirect_stdout(io.StringIO()):analyze_cli(args)
                human_metrics=json.loads((out/'manual_event_validation.json').read_text())
                self.assertGreater(human_metrics['precision'],.9);self.assertGreater(human_metrics['recall'],.9)
                args.manual_events=None
                with contextlib.redirect_stdout(io.StringIO()):analyze_cli(args)
                self.assertFalse((out/'manual_event_validation.json').exists())
                self.assertEqual(sha256(path),original)
                args.out_dir=root
                with self.assertRaises(ValueError):analyze_cli(args)

        def test_diagnostic_instrumentation_preserves_tracking(self):
            t,x,_=self.pulses(jitter=True,double=True,missing=True);fs=20000
            pitch=pitch_candidates(x[(t>=.1)&(t<.5)],fs,TRACKING_CONFIG)
            a=extract_track(t,x,fs,.1,.5,pitch['period_s'],10,TRACKING_CONFIG)
            b=extract_track(t,x,fs,.1,.5,pitch['period_s'],10,TRACKING_CONFIG,diagnostics=True)
            self.assertEqual(clean_json(a),clean_json(b))
            self.assertTrue(b['_diagnostic']['breaks'])
            for row in b['_diagnostic']['breaks']:
                if row['candidate_count']==0:self.assertIn('NO_AVAILABLE_CANDIDATE_IN_SEARCH_RANGE',row['direct_reasons'])
            cfg=dict(TRACKING_CONFIG,min_correlation=1.00001)
            c=extract_track(t,x,fs,.1,.5,pitch['period_s'],10,cfg,diagnostics=True)
            self.assertTrue(any('CORRELATION_BELOW_MIN' in r['direct_reasons'] for r in c['_diagnostic']['breaks']))

        def test_diagnostic_pooling_never_bridges_and_exact_windows(self):
            result=diagnostic_block_metrics([np.ones(5),np.ones(5)*2])
            self.assertEqual(result['n_cycles'],10);self.assertEqual(result['n_adjacent_period_differences'],8)
            self.assertEqual(result['n_RAP_terms'],6);self.assertEqual(result['n_PPQ5_terms'],2)
            self.assertEqual(result['jitter_local_percent'],0);self.assertEqual(result['rap_percent'],0);self.assertEqual(result['ppq5_percent'],0)
            track={'method':'A','events':[{'event_id':str(i),'subsample_time_s':t} for i,t in enumerate([0,.1,.2,.5,.6,.7])],
                   'segments':[{'segment_id':1,'event_indices':[0,1,2]},{'segment_id':2,'event_indices':[3,4,5]}]}
            blocks,details=diagnostic_clip_blocks(track,.1,.6)
            self.assertEqual(len(blocks),2);self.assertEqual(len(details),2)
            self.assertTrue(all(r['period_s']<.11 for r in details))

        def test_diagnostic_cli_keeps_original_baseline(self):
            self.assertIn('--diagnose-tracking',parser().format_help())
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);path=root/'flow.dat';t,x,_=self.pulses(duration=.30,missing=True)
                np.savetxt(path,np.c_[np.arange(len(t))*5,x],header='step airflow')
                argv=['--input',str(path),'--dt','1e-5','--start','.06','--end','.26','--out-dir',str(root/'original')]
                with contextlib.redirect_stdout(io.StringIO()):before=analyze_cli(parser().parse_args(argv))
                protected={p:sha256(p) for p in (root/'original').rglob('*') if p.is_file()}
                argv[-1]=str(root/'diagnosis');argv+=['--diagnose-tracking','--diagnostic-baseline',str(root/'original/summary.json')]
                with contextlib.redirect_stdout(io.StringIO()):after=analyze_cli(parser().parse_args(argv))
                self.assertEqual(before['primary']['main']['event_ids'],after['primary']['main']['event_ids'])
                self.assertEqual(before['cpp']['flow_CPP_dB'],after['cpp']['flow_CPP_dB'])
                for p,digest in protected.items():self.assertEqual(sha256(p),digest)
                for filename in ('break_diagnostics.csv','break_candidates.csv','uncovered_intervals.csv','same_window_metrics.csv','same_window_periods.csv','segment_jitter_comparison.csv','diagnostic_sensitivity.csv','new_connections_vs_strict.csv','sensitivity_periods.csv','sensitivity_segments.csv','sensitivity_events.csv','new_connection_figure_index.csv'):
                    self.assertTrue((root/'diagnosis/tracking_diagnostics'/filename).is_file(),filename)
                args=parser().parse_args(argv+['--overwrite']);args.out_dir=root/'original';args.diagnostic_baseline=None
                with self.assertRaisesRegex(ValueError,'diagnostic mode requires a new'):analyze_cli(args)
                for p,digest in protected.items():self.assertEqual(sha256(p),digest)

        @staticmethod
        def changing_pulses():
            t=np.arange(12000)/20000;truth=np.arange(.012,.585,.008);x=np.zeros_like(t)
            for at in truth:
                width=.0002 if at<.31 else .0005 if at<.319 else .0012
                x+=(1.2 if .29<at<.305 else 1.)*np.exp(-.5*((t-at)/width)**2)
            return t,x,truth

        def test_reanchor_context_recovers_without_relaxing_or_moving_events(self):
            t,x,truth=self.changing_pulses();cfg=dict(TRACKING_CONFIG,reanchor_policy='recover-rejected')
            a=extract_track_legacy(t,x,20000,.1,.5,.008,10,TRACKING_CONFIG,'A',True)
            b=extract_track(t,x,20000,.1,.5,.008,10,cfg,'A',True)
            self.assertEqual([(e['event_id'],e['subsample_time_s']) for e in a['events']],[(e['event_id'],e['subsample_time_s']) for e in b['events']])
            self.assertTrue(b['recovered_connections']);self.assertLess(b['anchor_contexts_tried'],len(b['events'])+1)
            for r in b['recovered_connections']:
                self.assertGreaterEqual(r['rho'],.65);self.assertTrue(r['score_margin'] is None or r['score_margin']>=.1)
                self.assertTrue(.65<=r['interval_over_acf_period']<=1.45)
            ss=summarize_track(b,.1,.5,.008,10,20000)
            self.assertGreater(ss['n_complete_periods_all_segments'],summarize_track(a,.1,.5,.008,10,20000)['n_complete_periods_all_segments'])
            self.assertLess(ss['pooled_within_segments_v2']['jitter_local_percent'],1e-4)
            ids=[j for s in b['segments'] for j in s['event_indices']];self.assertEqual(len(ids),len(set(ids)))
            spans=sorted((s['event_indices'][0],s['event_indices'][-1]) for s in b['segments'])
            self.assertTrue(all(a1<b0 for (_,a1),(b0,_) in zip(spans,spans[1:])))
            for row in b['recovered_connections']:
                self.assertLess(min(abs(truth-row['start_s'])),30e-6)
                self.assertLess(min(abs(truth-row['end_s'])),30e-6)
            self.assertTrue(any(r.get('resolved_by_recovery') for r in b['_diagnostic']['breaks']))

        def test_reanchor_stress_missing_multiple_peaks_and_phase(self):
            cfg=dict(TRACKING_CONFIG,reanchor_policy='recover-rejected')
            for jitter,double,missing in [(False,False,False),(True,True,False),(True,False,True),(False,True,True)]:
                t,x,truth=self.pulses(jitter=jitter,double=double,missing=missing)
                tr=extract_track(t,x,20000,.1,.5,.008,10,cfg,'A',True)
                for seg in tr['segments']:
                    times=[tr['events'][j]['subsample_time_s'] for j in seg['event_indices']]
                    for a,b in zip(times,times[1:]):
                        self.assertTrue(.65*.008-1e-4<=b-a<=1.45*.008+1e-4)
                    for at in times:self.assertLess(min(abs(truth-at)),.0005)
                for row in tr['recovered_connections']:self.assertEqual(row['safety_reasons'],[])

        def test_reanchor_cap_and_strict_gate_boundary(self):
            t,x,_=self.changing_pulses();cfg=dict(TRACKING_CONFIG,reanchor_policy='recover-rejected',reanchor_attempt_cap=1)
            tr=extract_track(t,x,20000,.1,.5,.008,10,cfg,'A',True)
            self.assertTrue(tr['search_limit']);self.assertIn('SEARCH_LIMIT',summarize_track(tr,.1,.5,.008,10,20000)['flags'])
            # Tune a real observed failed context to either side of its measured rho.
            seed=extract_track_legacy(t,x,20000,.1,.5,.008,10,TRACKING_CONFIG,'A',True)
            rho=next(r['best_rho'] for r in seed['_diagnostic']['breaks'] if r['best_rho'] is not None)
            for threshold in (rho-1e-8,rho+1e-8,.65):
                cfg=dict(TRACKING_CONFIG,reanchor_policy='recover-rejected',min_correlation=threshold)
                tr=extract_track(t,x,20000,.1,.5,.008,10,cfg,'A',True)
                self.assertTrue(all(r['rho']>=threshold for r in tr['recovered_connections']))

        def test_reanchor_shape_correlation_near_point65(self):
            axis=np.arange(-64,65)/20000;small=np.exp(-.5*(axis/.0002)**2);small-=small.mean();small/=np.linalg.norm(small)
            lo,hi=.0005,.0012
            for _ in range(35):
                mid=(lo+hi)/2;w=np.exp(-.5*(axis/mid)**2);w-=w.mean();rho=float(np.dot(small,w/np.linalg.norm(w)))
                if rho>.65:lo=mid
                else:hi=mid
            t=np.arange(12000)/20000;x=np.zeros_like(t)
            for at in np.arange(.012,.585,.008):
                width=.0002 if at<.31 else .0005 if at<.319 else mid
                x+=(1.2 if .29<at<.305 else 1.)*np.exp(-.5*((t-at)/width)**2)
            cfg=dict(TRACKING_CONFIG,reanchor_policy='recover-rejected');tr=extract_track(t,x,20000,.1,.5,.008,10,cfg,'A',True)
            rhos=[r['rho'] for r in tr['attempts']]
            self.assertLess(min(abs(v-.65) for v in rhos),.001)
            self.assertTrue(all(r['rho']>=.65 for r in tr['attempts'] if r['accepted']))

        def test_reanchor_incomplete_template_edges_remain_unaccepted(self):
            t=np.arange(12000)/20000;x=np.zeros_like(t)
            for at in [.0005,*np.arange(.012,.585,.008),.599]:x+=np.exp(-.5*((t-at)/.00035)**2)
            tr=extract_track(t,x,20000,0,.6,.008,10,dict(TRACKING_CONFIG,reanchor_policy='recover-rejected'),'A',True)
            edges=[r for r in tr['candidate_states'] if r['state']=='INSUFFICIENT_TEMPLATE_EDGE']
            self.assertGreaterEqual(len(edges),2)
            edge_ids={r['event_id'] for r in edges};self.assertFalse(any(e['accepted'] and e['event_id'] in edge_ids for e in tr['events']))

        def test_reanchor_amplitude_shape_and_competing_directions(self):
            t,x,_=self.changing_pulses();x*=1+.35*np.sin(2*np.pi*3*t)
            x+=.25*np.exp(-.5*((t-.334)/.00016)**2)
            tr=extract_track(t,x,20000,.1,.5,.008,10,dict(TRACKING_CONFIG,reanchor_policy='recover-rejected'),'A',True)
            spans=sorted((s['event_indices'][0],s['event_indices'][-1]) for s in tr['segments'])
            self.assertTrue(all(a1<b0 for (_,a1),(b0,_) in zip(spans,spans[1:])))
            self.assertEqual(len(tr['events']),len({e['event_id'] for e in tr['events']}))
            self.assertTrue(any(r['state']=='INSUFFICIENT_TEMPLATE_EDGE' for r in tr['candidate_states']) or any(not e['in_window'] for e in tr['events']))
            keys=[tuple(r['context_key']) for r in tr['recovery_attempts']];self.assertEqual(len(keys),len(set(keys)))

        def test_metric_v2_denominators_exclude_noncontributing_segments(self):
            blocks=[[.010],[.008,.009],[.012,.011,.010],[.020,.021,.022,.023,.024]]
            v1=diagnostic_block_metrics(blocks);v2=metric_blocks_v2(blocks,1.)
            self.assertEqual(v1['n_cycles'],11);self.assertEqual(v2['jitter_local_percent_n_cycles'],10)
            self.assertEqual(v2['rap_percent_n_cycles'],8);self.assertEqual(v2['ppq5_percent_n_cycles'],5)
            self.assertEqual(v2['jitter_local_percent_n_terms'],7);self.assertEqual(v2['rap_percent_n_terms'],4);self.assertEqual(v2['ppq5_percent_n_terms'],1)
            self.assertAlmostEqual(v2['jitter_local_percent'],100*.001/np.mean(blocks[1]+blocks[2]+blocks[3]))
            self.assertAlmostEqual(v2['rap_percent'],0,places=10);self.assertAlmostEqual(v2['ppq5_percent'],0,places=10)

        def test_manual_assignment_global_nearest_and_tolerance(self):
            pairs=optimal_event_assignment([0.,150e-6],[90e-6],100e-6)
            self.assertEqual([(a,b) for a,b,_ in pairs],[(1,0)])
            self.assertEqual(len(optimal_event_assignment([0,.001],[.00008,.00108],100e-6)),2)
            self.assertEqual(len(optimal_event_assignment([0,.001],[.00008,.00108],50e-6)),0)
            self.assertEqual(len(optimal_event_assignment([0,.00015],[.00008,.00017],100e-6)),2)

        def test_manual_reference_validation_and_uncertain_period_breaks(self):
            with tempfile.TemporaryDirectory() as name:
                path=Path(name)/'manual.csv';t=np.arange(8000)/20000;x=np.cos(2*np.pi*125*t);data={'t':t,'x':x}
                rows=[{'case_id':'synthetic','phase':'A','event_time_s':.104+i*.008,'confidence':'uncertain' if i==3 else 'high','reviewer':'test','notes':''} for i in range(8)]
                write_csv(path,rows);parsed=read_manual_reference(path,'synthetic',data,.1,.2);self.assertEqual(len(parsed),8)
                events=[{'event_id':str(i),'subsample_time_s':r['event_time_s'],'accepted':True} for i,r in enumerate(rows)]
                track={'method':'A','events':events,'segments':[{'segment_id':1,'event_indices':list(range(8))}]}
                stats,matched,checks,periods=evaluate_manual_reference(parsed,{'A':track,'B':{'events':[],'segments':[]}},.1,.2)
                high=next(r for r in stats if r['scope']=='high_only' and np.isclose(r['tolerance_us'],100))
                self.assertEqual(high['FP'],0);self.assertEqual(high['FN'],0);self.assertAlmostEqual(high['precision'],1)
                self.assertEqual(high['manual_jitter_local_percent_n_cycles'],5)
                self.assertEqual(high['auto_jitter_local_percent_n_cycles'],5)
                self.assertTrue(all(not (r['start_s']<.128<r['end_s']) for r in periods if r['scope']=='high_only'))
                for mutation in ('duplicate','outside','phase','confidence','case'):
                    changed=copy.deepcopy(rows)
                    if mutation=='duplicate':changed[1]['event_time_s']=changed[0]['event_time_s']
                    if mutation=='outside':changed[0]['event_time_s']=104000
                    if mutation=='phase':changed[0]['phase']='B'
                    if mutation=='confidence':changed[0]['confidence']='certain'
                    if mutation=='case':
                        for r in changed:r['case_id']='different'
                    write_csv(path,changed)
                    with self.assertRaises(ValueError):read_manual_reference(path,'synthetic',data,.1,.2)
                write_csv(path,[],['case_id','phase','event_time_s','confidence','reviewer','notes']);self.assertEqual(read_manual_reference(path,'synthetic',data,.1,.2),[])

        def test_manual_missing_extra_phase_and_no_cross_gap_metrics(self):
            refs=[{'phase':'A','event_time_s':.104+i*.008,'confidence':'high'} for i in range(8)]
            times=[r['event_time_s'] for r in refs];automatic=times[:3]+times[4:]
            track={'method':'A','events':[{'event_id':str(i),'subsample_time_s':at,'accepted':True} for i,at in enumerate(automatic)],
                   'segments':[{'segment_id':1,'event_indices':[0,1,2]},{'segment_id':2,'event_indices':[3,4,5,6]}]}
            stats,_,checks,periods=evaluate_manual_reference(refs,{'A':track,'B':{'events':[],'segments':[]}},.1,.2)
            main=next(r for r in stats if r['scope']=='high_only' and np.isclose(r['tolerance_us'],100))
            self.assertEqual(main['TP'],7);self.assertEqual(main['FN'],1);self.assertEqual(main['FP'],0)
            self.assertAlmostEqual(main['auto_jitter_local_percent'],0,places=10)
            self.assertEqual(main['auto_jitter_local_percent_n_terms'],3)
            # A known merged edge spans the missing reference; a spurious event adds an FP.
            track['events'].append({'event_id':'spurious','subsample_time_s':.132,'accepted':True})
            track['segments']=[{'segment_id':1,'event_indices':[0,1,2,7,3,4,5,6]}]
            stats,_,checks,_=evaluate_manual_reference(refs,{'A':track,'B':{'events':[],'segments':[]}},.1,.2)
            main=next(r for r in stats if r['scope']=='high_only' and np.isclose(r['tolerance_us'],100))
            self.assertEqual(main['FP'],1);self.assertGreater(main['n_automatic_periods_spanning_manual_event'],0)
            self.assertTrue(all(r['phase']=='A' for r in stats))

        def test_recovery_cli_review_pending_baseline_and_cpp(self):
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);path=root/'flow.dat';t,x,_=self.changing_pulses();np.savetxt(path,np.c_[np.arange(len(t))*5,x])
                args=parser().parse_args(['--input',str(path),'--dt','1e-5','--start','.1','--end','.5','--out-dir',str(root/'new'),
                                          '--reanchor-policy','recover-rejected','--export-manual-review','--case-id','synthetic'])
                digest=sha256(path)
                with contextlib.redirect_stdout(io.StringIO()):s=analyze_cli(args)
                self.assertEqual(s['manual_validation_status'],'PENDING_MANUAL_ANNOTATION');self.assertEqual(sha256(path),digest)
                for fn in ['reanchor_comparison.csv','jitter_comparison.csv','metric_terms_v2.csv','manual_review_index.csv','manual_events_template.csv','manual_validation.csv','tracking_breaks_detailed.csv','phase_features.csv','new_connections.csv']:
                    self.assertTrue((root/'new'/fn).is_file(),fn)
                self.assertTrue(s['reanchor_validation']['baseline_check']['CPP_exactly_identical_current_runtime'])
                with (root/'new/manual_events_template.csv').open() as f:self.assertEqual(list(csv.DictReader(f)),[])
                with (root/'new/manual_review_index.csv').open() as f:panels=[r for r in csv.DictReader(f) if r['panel_id'].startswith('panel')]
                self.assertAlmostEqual(float(panels[0]['t_start_s']),.1);self.assertAlmostEqual(float(panels[-1]['t_end_s']),.5)
                for a,b in zip(panels,panels[1:]):self.assertLess(float(b['t_start_s']),float(a['t_end_s']))
                args.overwrite=True
                with self.assertRaises(ValueError):analyze_cli(args)

                reference=root/'filled_synthetic_reference.csv';truth=np.arange(.012,.585,.008)
                write_csv(reference,[{'case_id':'synthetic','phase':'A','event_time_s':float(at),'confidence':'high','reviewer':'synthetic_test','notes':'known synthetic center'} for at in truth if .1<=at<.5])
                args.out_dir=root/'manual_new';args.export_manual_review=False;args.manual_events=reference
                with contextlib.redirect_stdout(io.StringIO()):validated=analyze_cli(args)
                self.assertEqual(validated['manual_validation_status'],'MANUAL_REFERENCE_EVALUATED')
                self.assertEqual(validated['primary']['manual_validation_status'],'MANUAL_REFERENCE_EVALUATED')
                self.assertEqual(validated['control']['manual_validation_status'],'PENDING_MANUAL_ANNOTATION')
                self.assertGreater(validated['manual_event_validation']['precision'],.99)
                self.assertGreater(validated['manual_event_validation']['recall'],.99)
                self.assertTrue((root/'manual_new/manual_metric_periods.csv').is_file())

        def test_global_constant_period_entire_window_real_extrema(self):
            t,x,truth=self.pulses();data={'t':t,'x':x,'fs':20000};r=global_cycle_track(data,.1,.5,GLOBAL_CYCLE_CONFIG)
            for phase,v in r['methods'].items():
                self.assertAlmostEqual(v['classified_time_fraction'],1);self.assertLess(v['metrics']['jitter_local_percent'],1e-5)
                self.assertGreater(v['time_coverage'],.94)
                ids=[e['event_id'] for e in v['track']['events'] if e['accepted']];self.assertEqual(len(ids),len(set(ids)))
                for e in v['track']['events']:
                    self.assertLessEqual(abs(e['subsample_time_s']-e['raw_time_s']),.5/20000+1e-12)
                    if e['accepted'] and phase=='A':self.assertLess(min(abs(truth-e['subsample_time_s'])),1e-5)

        def test_global_known_jitter_not_regularized_by_DP(self):
            fs=20000;t=np.arange(12000)/fs;truth=[.012]
            for i in range(80):truth.append(truth[-1]+.008*(1+.08*np.sin(i*.7)))
            truth=np.array(truth);x=sum(np.exp(-.5*((t-at)/.00045)**2) for at in truth);r=global_cycle_track({'t':t,'x':x,'fs':fs},.1,.5,GLOBAL_CYCLE_CONFIG)
            expected=period_metrics(np.diff(truth[(truth>=.1)&(truth<.5)]))
            self.assertAlmostEqual(r['methods']['A']['metrics']['jitter_local_percent'],expected['jitter_local_percent'],delta=.002)
            self.assertEqual(GLOBAL_CYCLE_CONFIG['period_regularization_weight'],0)

        def test_global_double_peaks_are_not_independent_cycles(self):
            t,x,truth=self.pulses(double=True,jitter=True);r=global_cycle_track({'t':t,'x':x,'fs':20000},.1,.5,GLOBAL_CYCLE_CONFIG)
            tr=r['methods']['A']['track'];adopted=[e for e in tr['events'] if e['accepted']]
            self.assertGreater(len(adopted),40)
            for e in adopted:self.assertLess(min(abs(truth-e['subsample_time_s'])),.0005)
            for p in r['methods']['A']['summary']['_periods']:self.assertGreater(p['period_s'],.0044)

        def test_global_plateau_uncertainty_not_precision(self):
            t=np.arange(12000)/20000;x=np.maximum(np.sin(2*np.pi*125*t),0.)
            r=global_cycle_track({'t':t,'x':x,'fs':20000},.1,.5,GLOBAL_CYCLE_CONFIG)
            b=[e for e in r['methods']['B']['track']['events'] if e['accepted']]
            self.assertTrue(b);self.assertTrue(any(e['plateau_width_s']>.001 for e in b))
            self.assertTrue(all(e['confidence_class']=='ESTIMATED' for e in b))
            for e in b:self.assertTrue(e['plateau_start_s']<=e['raw_time_s']<=e['plateau_end_s'])

        def test_global_modulation_shape_change_low_amplitude(self):
            t,x,truth=self.changing_pulses();x*=1+.65*np.sin(2*np.pi*2*t)
            r=global_cycle_track({'t':t,'x':x,'fs':20000},.1,.5,GLOBAL_CYCLE_CONFIG)
            a=r['methods']['A'];self.assertGreater(a['time_coverage'],.90)
            for e in a['track']['events']:
                if e['accepted']:self.assertLess(min(abs(truth-e['subsample_time_s'])),.0002)
            self.assertAlmostEqual(a['classified_time_fraction'],1)

        def test_global_aperiodic_subharmonic_and_no_oscillation(self):
            t=np.arange(12000)/20000;data={'t':t,'x':np.zeros_like(t),'fs':20000}
            zero=global_cycle_track(data,.1,.5,GLOBAL_CYCLE_CONFIG)
            self.assertIsNone(zero['methods']['A']['metrics']['jitter_local_percent']);self.assertAlmostEqual(zero['methods']['A']['unresolved_time_fraction'],1)
            data['x']=np.random.default_rng(77).normal(size=len(t));noise=global_cycle_track(data,.1,.5,GLOBAL_CYCLE_CONFIG)
            self.assertIsNone(noise['methods']['A']['metrics']['jitter_local_percent'])
            data['x']=sum((1 if i%2==0 else .2)*np.exp(-.5*((t-at)/.00045)**2) for i,at in enumerate(np.arange(.012,.59,.008)))
            alternating=global_cycle_track(data,.1,.5,GLOBAL_CYCLE_CONFIG)
            self.assertEqual(alternating['methods']['A']['period_identification_status'],'PERIOD_AMBIGUOUS')
            self.assertIsNone(alternating['methods']['A']['metrics']['jitter_local_percent'])
            self.assertIn('conditional_metrics',alternating['methods']['A']['metrics'])

        def test_global_missing_pulse_restarts_not_twoT(self):
            t,x,_=self.pulses(missing=True);x[(t>.28)&(t<.34)]=0.
            r=global_cycle_track({'t':t,'x':x,'fs':20000},.1,.5,GLOBAL_CYCLE_CONFIG)
            a=r['methods']['A'];self.assertGreater(a['unresolved_time_fraction'],0)
            self.assertTrue(any(e['accepted'] and e['subsample_time_s']>.45 for e in a['track']['events']))
            self.assertTrue(all(not (p['start_s']<.28 and p['end_s']>.34) for p in a['summary']['_periods']))
            self.assertFalse(any(row['adopted'] and 'DOUBLE_PERIOD_POSSIBLE' in row['rejection_reasons'] for row in r['edge_audit']))

        def test_global_dp_uses_future_evidence_and_alternative_margin(self):
            cfg=dict(GLOBAL_CYCLE_CONFIG,root_penalty=0.,restart_penalty=2.);T=np.full((4,4),-np.inf);C=np.zeros((4,4),bool)
            # The locally stronger middle candidate has poor future fit.
            for i,j,value in [(0,1,.5),(0,2,.3),(1,3,-.8),(2,3,.8)]:T[i,j]=value;C[i,j]=True
            paths,marginal=global_dp([1.,1.1,1.,1.],T,C,cfg)
            self.assertEqual(paths[0]['nodes'],(0,2,3));self.assertGreater(paths[0]['score'],paths[1]['score'])
            self.assertGreater(marginal[2],marginal[1])

        def test_global_profiles_shared_and_confidence_only_sensitivity(self):
            t,x,_=self.pulses();data={'t':t,'x':x,'fs':20000};a=global_cycle_track(data,.1,.5,GLOBAL_CYCLE_CONFIG)
            b=global_cycle_track(data,.1,.5,dict(GLOBAL_CYCLE_CONFIG,effective_shape_correlation=.7),a['_prepared'],False)
            for phase in ('A','B'):self.assertEqual([e['event_id'] for e in a['methods'][phase]['track']['events'] if e['accepted']],[e['event_id'] for e in b['methods'][phase]['track']['events'] if e['accepted']])
            profiles=global_threshold_profiles(GLOBAL_CYCLE_CONFIG);grid={(p['config']['prominence_percent'],p['config']['excursion_percent']) for p in profiles}
            self.assertTrue(all((p,a) in grid for p in (3,5,8,10,12,15,20) for a in (3,5,10,15,20)))
            self.assertEqual(len(profiles),43)

        def test_global_search_cap_and_profile_validation(self):
            t,x,_=self.pulses()
            with self.assertRaisesRegex(ValueError,'SEARCH_LIMIT'):global_cycle_track({'t':t,'x':x,'fs':20000},.1,.5,dict(GLOBAL_CYCLE_CONFIG,max_dp_candidates=1))
            with tempfile.TemporaryDirectory() as name:
                path=Path(name)/'profile.json';dump_json(path,{'config':dict(GLOBAL_CYCLE_CONFIG,period_regularization_weight=1.)})
                args=parser().parse_args([])
                with self.assertRaises(ValueError):resolve_global_profile(str(path),args)

        def test_global_sampling_phases_same_real_events(self):
            t=np.arange(50000)/100000;truth=np.arange(.012,.49,.008);x=sum(np.exp(-.5*((t-at)/.00045)**2) for at in truth)
            full=global_cycle_track({'t':t,'x':x,'fs':100000},.1,.4,GLOBAL_CYCLE_CONFIG,alternatives=False)
            for phase in range(5):
                dec=global_cycle_track({'t':t[phase::5],'x':x[phase::5],'fs':20000},.1,.4,GLOBAL_CYCLE_CONFIG,alternatives=False)
                self.assertLess(dec['methods']['A']['metrics']['jitter_local_percent'],.001)
                self.assertAlmostEqual(dec['methods']['A']['metrics']['jitter_local_percent'],full['methods']['A']['metrics']['jitter_local_percent'],delta=.001)

        def test_global_cli_panels_metrics_manifest_and_source_protection(self):
            with tempfile.TemporaryDirectory() as name:
                root=Path(name);t,x,_=self.pulses(duration=.32,jitter=True);source=root/'flow.dat';np.savetxt(source,np.c_[np.arange(len(t))*5,x]);digest=sha256(source)
                args=parser().parse_args(['--input',str(source),'--dt','1e-5','--start','.08','--end','.28','--tracking-method','global-cycle','--skip-threshold-sweep','--case-id','synthetic','--out-dir',str(root/'out')])
                with contextlib.redirect_stdout(io.StringIO()):summary=analyze_cli(args)
                self.assertEqual(summary['manual_validation_status'],'PENDING_MANUAL_ANNOTATION');self.assertEqual(sha256(source),digest)
                for p in ('flow_full_window_events.png','threshold_grid.png','acf_timeline.png','event_timeline.png','event_selection_examples.png','jitter_rap_ppq5_cpp_comparison.png','run_manifest.json','event_candidates_all.csv','time_classification.csv','metric_terms_global.csv'):
                    self.assertTrue((root/'out'/p).is_file(),p)
                with (root/'out/manual_review_index.csv').open() as f:panels=[r for r in csv.DictReader(f) if r['panel_id'].startswith('panel')]
                self.assertEqual(float(panels[0]['t_start_s']),.08);self.assertEqual(float(panels[-1]['t_end_s']),.28)
                self.assertTrue(all(.05-1e-12<=float(r['t_end_s'])-float(r['t_start_s'])<=.07+1e-12 for r in panels))
                manifest=json.loads((root/'out/run_manifest.json').read_text());self.assertTrue(all(sha256(root/'out'/p)==h for p,h in manifest['files'].items()))
                args.overwrite=True
                with self.assertRaises(ValueError):analyze_cli(args)
                reference=root/'manual.csv';truth=self.pulses(duration=.32,jitter=True)[2]
                write_csv(reference,[{'case_id':'synthetic','phase':'A','event_time_s':float(at),'confidence':'high','reviewer':'synthetic_test','notes':'known center'} for at in truth if .08<=at<.28])
                args.out_dir=root/'manual_check';args.manual_events=reference
                with contextlib.redirect_stdout(io.StringIO()):validated=analyze_cli(args)
                self.assertEqual(validated['manual_validation_status'],'MANUAL_REFERENCE_EVALUATED')
                stats=next(r for r in validated['manual_validation'] if r['phase']=='A' and r['scope']=='high_only' and np.isclose(r['tolerance_us'],100))
                self.assertGreater(stats['precision'],.95);self.assertGreater(stats['recall'],.95)

        def test_cleanup_blocks_roots_sources_symlinks_tracked_and_hash_mismatch(self):
            with tempfile.TemporaryDirectory() as name:
                repo=Path(name);base=repo/'analysis_runs';base.mkdir();safe=base/'work/CMakeFiles/CMakeScratch';safe.mkdir(parents=True)
                primary=base/'primary';primary.mkdir();(primary/'airflow_vt.dat').write_text('0 1\n');link=base/'escape';link.symlink_to(repo,target_is_directory=True)
                tracked=base/'tracked';tracked.mkdir();(tracked/'x.txt').write_text('tracked')
                items=[]
                for p in [base,primary,link,tracked,safe]:items.append({'relative_path':str(p.relative_to(repo)),'kind':'EMPTY_BUILD_CACHE','expected_tree_sha256':analysis_tree_snapshot(p)['tree_sha256'],'reason':'test'})
                items[-1]['expected_tree_sha256']='wrong'
                deleted,kept=cleanup_analysis_allowlist(repo,repo/'audit',items,[],['analysis_runs/tracked/x.txt'])
                self.assertEqual(deleted,[]);self.assertEqual(len(kept),5);self.assertTrue((primary/'airflow_vt.dat').exists());self.assertTrue(safe.exists())

        def test_cleanup_only_allowlisted_empty_cache_deleted(self):
            with tempfile.TemporaryDirectory() as name:
                repo=Path(name);p=repo/'analysis_runs/work/CMakeFiles/CMakeScratch';p.mkdir(parents=True);keep=repo/'analysis_runs/keep';keep.mkdir();(keep/'data.dat').write_text('keep')
                row={'relative_path':str(p.relative_to(repo)),'kind':'EMPTY_BUILD_CACHE','expected_tree_sha256':analysis_tree_snapshot(p)['tree_sha256'],'reason':'empty regenerable CMake scratch'}
                deleted,kept=cleanup_analysis_allowlist(repo,repo/'audit',[row],[keep],[])
                self.assertEqual(len(deleted),1);self.assertFalse(p.exists());self.assertTrue((keep/'data.dat').exists())

        def test_cleanup_owned_development_requires_replacement_and_no_human_data(self):
            with tempfile.TemporaryDirectory() as name:
                repo=Path(name);p=repo/'analysis_runs/airflow_global_cycle_validation_test/development';p.mkdir(parents=True);(p/'summary.json').write_text('{}')
                row={'relative_path':str(p.relative_to(repo)),'kind':'OWNED_DERIVED_DEVELOPMENT','expected_tree_sha256':analysis_tree_snapshot(p)['tree_sha256'],'replacement_verified':False,'reason':'synthetic development output'}
                deleted,kept=cleanup_analysis_allowlist(repo,repo/'audit1',[row],[],[]);self.assertFalse(deleted);self.assertTrue(p.exists())
                row['replacement_verified']=True
                write_csv(p/'manual_events.csv',[{'case_id':'r','phase':'A','event_time_s':.1,'confidence':'high','reviewer':'human','notes':''}])
                row['expected_tree_sha256']=analysis_tree_snapshot(p)['tree_sha256']
                deleted,kept=cleanup_analysis_allowlist(repo,repo/'audit2',[row],[],[]);self.assertFalse(deleted);self.assertTrue((p/'manual_events.csv').exists())
                (p/'manual_events.csv').unlink();row['expected_tree_sha256']=analysis_tree_snapshot(p)['tree_sha256']
                deleted,kept=cleanup_analysis_allowlist(repo,repo/'audit3',[row],[],[]);self.assertEqual(len(deleted),1);self.assertFalse(p.exists())

        def test_global_manual_common_periods_do_not_bridge_missing_reference(self):
            refs=[{'phase':'A','event_time_s':.104+i*.008,'confidence':'high'} for i in range(8)]
            times=[r['event_time_s'] for i,r in enumerate(refs) if i!=3]
            track={'events':[{'event_id':str(i),'accepted':True,'subsample_time_s':at+80e-6} for i,at in enumerate(times)],
                   'segments':[{'segment_id':1,'event_indices':list(range(7))}]}
            track_map={'A':track,'B':{'events':[],'segments':[]}}
            _,matches,_,_=evaluate_manual_reference(refs,track_map,.1,.2)
            stats,periods=manual_common_period_comparison(refs,track_map,matches,.1,.2)
            row=next(r for r in stats if r['scope']=='high_only' and r['tolerance_us']==100)
            self.assertEqual(row['manual_common_jitter_local_percent_n_cycles'],5);self.assertEqual(row['auto_common_jitter_local_percent_n_cycles'],5)
            self.assertEqual(row['manual_common_jitter_local_percent_n_terms'],3)
            self.assertAlmostEqual(row['jitter_local_percent_difference_pp'],0,places=9)
            self.assertTrue(all(r['auto_end_s']-r['auto_start_s']<.009 for r in periods))

    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():return 1
    print('SELF_TEST_PASS: fixed-seed synthetic metrics, tracking, CPP, sampling phases and protected CLI output')
    return 0


def main(argv=None):
    p=parser();args=p.parse_args(argv)
    if args.self_test:return self_test()
    try:analyze_cli(args)
    except (ValueError,OSError,RuntimeError,json.JSONDecodeError) as exc:
        p.exit(2,'ERROR: '+str(exc)+'\n')
    return 0

if __name__=='__main__':raise SystemExit(main())
