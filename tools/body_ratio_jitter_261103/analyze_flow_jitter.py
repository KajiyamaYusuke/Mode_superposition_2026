#!/usr/bin/env python3
"""Conservative Ug period primitives; no automatic production QC approval."""
import numpy as np
from scipy.signal import find_peaks, correlate


def period_metrics(periods):
    p = np.asarray(periods, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any(p <= 0):
        raise ValueError('periods must be a finite positive vector')
    n = len(p)
    if n < 2:
        return {'local': np.nan, 'rap': np.nan, 'ppq5': np.nan, 'cv': np.nan}
    mean = p.mean()
    return {'local': 100*np.mean(np.abs(np.diff(p)))/mean,
            'rap': 100*np.mean(np.abs(p[1:-1] - (p[:-2]+p[1:-1]+p[2:])/3))/mean if n>=3 else np.nan,
            'ppq5': 100*np.mean(np.abs(p[2:-2] - np.convolve(p, np.ones(5)/5, mode='valid')))/mean if n>=5 else np.nan,
            'cv': 100*np.std(p)/mean}


def events(values, dt=5e-5, prominence_fraction=0.1, method='A', interpolate=True, reference=None):
    x=np.asarray(values,dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all() or not np.isfinite(dt) or dt <= 0:
        raise ValueError('invalid waveform')
    if len(x)<3: return np.array([])
    reference = x if reference is None else np.asarray(reference, dtype=float)
    if reference.size == 0 or not np.isfinite(reference).all(): raise ValueError('invalid reference window')
    amplitude=np.percentile(reference,95)-np.percentile(reference,5)
    if amplitude <= max(1e-15, np.max(np.abs(reference))*1e-12): return np.array([])
    if method not in ('A','B'): raise ValueError('unknown method')
    signal=x if method=='A' else -x
    indices,_=find_peaks(signal,prominence=prominence_fraction*amplitude,distance=max(1,int(1/(500*dt))))
    offset=np.zeros(len(indices))
    if interpolate and len(indices):
        left,center,right=signal[indices-1],signal[indices],signal[indices+1]
        denom=left-2*center+right
        valid=denom<0
        offset[valid]=0.5*(left[valid]-right[valid])/denom[valid]
        offset=np.clip(offset,-0.5,0.5)
    return (indices+offset)*dt


def window_periods(t, start, end):
    t=np.asarray(t,dtype=float)
    if not np.isfinite(t).all() or np.any(np.diff(t)<=0) or start>=end: raise ValueError('invalid event times or window')
    # Both events belong to the half-open analysis window.
    inside=t[(t>=start)&(t<end)]
    return np.diff(inside)


def load_flow(path, dt=1e-5, nstep=100000):
    with open(path) as f:
        header=f.readline()
    if not header.startswith('#') or 'airflow' not in header.lower():
        raise ValueError('missing airflow header')
    data=np.loadtxt(path,comments='#',ndmin=2)
    expected=np.arange(0,nstep,5)
    if data.shape != (len(expected),2) or not np.isfinite(data).all() or not np.array_equal(data[:,0],expected):
        raise ValueError('invalid step series or nonfinite flow')
    return data[:,0]*dt,data[:,1]


def periodicity(values, dt, period):
    """Compare an event period to the ACF and spectrum without using FFT as events."""
    x = np.asarray(values, dtype=float)
    if x.size < 3 or not np.isfinite(x).all() or period <= 0:
        raise ValueError('invalid periodicity input')
    x = x - x.mean()
    energy = float(np.dot(x, x))
    if energy <= 0:
        return {'acf_at_period': np.nan, 'spectral_ratio': np.nan}
    lag = int(round(period / dt))
    if lag < 1 or lag >= len(x):
        return {'acf_at_period': np.nan, 'spectral_ratio': np.nan}
    acf = correlate(x, x, mode='full', method='fft')[len(x)-1:]
    acf_value = acf[lag] / acf[0]
    power = np.abs(np.fft.rfft(x * np.hanning(len(x)))) ** 2
    freq = np.fft.rfftfreq(len(x), dt)
    band = (freq >= 30) & (freq <= 500)
    if not band.any() or power[band].max() <= 0:
        ratio = np.nan
    else:
        target = abs(freq - 1.0/period) <= max(5.0, 2.0/(dt*len(x)))
        ratio = float(power[band & target].max() / power[band].max()) if np.any(band & target) else 0.0
    return {'acf_at_period': float(acf_value), 'spectral_ratio': ratio}


def event_correspondence(a, b, start, end):
    """Fraction of A-to-A cycles containing exactly one B event."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[(a >= start) & (a < end)]
    if len(a) < 2:
        return np.nan
    counts = np.searchsorted(b, a[1:], side='left') - np.searchsorted(b, a[:-1], side='left')
    return float(np.mean(counts == 1))


def analyze(values, dt=5e-5, start=.5, end=1.):
    """Exploratory analysis: conservative QC; unresolved cases must be reviewed."""
    x=np.asarray(values,dtype=float)
    if not np.isfinite(x).all():
        return {'status':'NUMERIC_INVALID'}
    selection=x[int(round(start/dt)):int(round(end/dt))]
    if len(selection)<3: return {'status':'INSUFFICIENT_CYCLES'}
    amp=np.percentile(selection,95)-np.percentile(selection,5)
    if amp<=max(1e-15,np.max(np.abs(selection))*1e-12): return {'status':'NO_OSCILLATION'}
    candidates={}
    for method, frac, interpolation in [('A',.05,True),('A',.10,True),('A',.15,True),('A',.10,False),('B',.10,True)]:
        key=f'{method}_{frac:.2f}_{int(interpolation)}'
        times=events(x,dt,frac,method,interpolation,selection)
        periods=window_periods(times,start,end)
        candidates[key]={'periods':periods,'times':times,'metrics':period_metrics(periods)}
    main=candidates['A_0.10_1']; p=main['periods']
    if len(p)<10: return {'status':'INSUFFICIENT_CYCLES','cycles':len(p),'candidates':candidates}
    med=np.median(p)
    # A period outside a broad 2:1 range indicates missing or extra events.
    if np.any((p<.5*med)|(p>1.5*med)):
        return {'status':'AMBIGUOUS_PEAKS','cycles':len(p),'reason':'outlying period','candidates':candidates}
    for name in ['A_0.05_1','A_0.15_1','B_0.10_1']:
        q=candidates[name]['periods']
        if len(q)<10 or abs(len(q)-len(p))>1 or abs(np.median(q)/med-1)>.1:
            return {'status':'AMBIGUOUS_PEAKS','cycles':len(p),'reason':'method or threshold disagreement','candidates':candidates}
    correspondence = event_correspondence(main['times'], candidates['B_0.10_1']['times'], start, end)
    evidence = periodicity(selection, dt, med)
    if correspondence < .9 or not np.isfinite(evidence['acf_at_period']) or evidence['acf_at_period'] < .2 or not np.isfinite(evidence['spectral_ratio']) or evidence['spectral_ratio'] < .05:
        return {'status':'AMBIGUOUS_PEAKS','cycles':len(p),'reason':'event correspondence or periodicity check', 'event_correspondence':correspondence, 'periodicity':evidence,'candidates':candidates}
    t=events(x,dt,.1,'A',True,selection)
    first=window_periods(t,start,(start+end)/2); second=window_periods(t,(start+end)/2,end)
    if min(len(first),len(second))<10:
        return {'status':'INSUFFICIENT_CYCLES','cycles':len(p),'candidates':candidates}
    m1=period_metrics(first);m2=period_metrics(second)
    amp1=np.percentile(selection[:len(selection)//2],95)-np.percentile(selection[:len(selection)//2],5)
    amp2=np.percentile(selection[len(selection)//2:],95)-np.percentile(selection[len(selection)//2:],5)
    drift=abs(np.mean(first)-np.mean(second))/med>.05 or abs(amp1-amp2)/max(amp1,amp2)>.2
    jitter_change=abs(m1['local']-m2['local'])>max(.2,.5*max(m1['local'],m2['local']))
    # Recognize consistent alternating periods only; call it a candidate.
    alt=np.mean(p[::2])-np.mean(p[1::2])
    multiplied=abs(alt)>.05*med and np.std(p[::2])+np.std(p[1::2])<abs(alt)
    # Alternating intervals can also arise from two peaks per physical cycle.
    # Until waveform morphology and subharmonics are reviewed, do not label it.
    status='AMBIGUOUS_PEAKS' if multiplied else ('VALID_WINDOW_DEPENDENT' if drift or jitter_change else 'VALID_STATIONARY')
    return {'status':status,'reason':'alternating intervals require waveform review' if multiplied else '', 'cycles':len(p),'metrics':main['metrics'],'half_metrics':[m1,m2], 'event_correspondence':correspondence, 'periodicity':evidence, 'candidates':candidates}
