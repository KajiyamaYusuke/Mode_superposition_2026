#!/usr/bin/env python3
"""Candidate waveform-template tracking; a detected period is not automatically valid."""
import numpy as np
from scipy.signal import correlate, find_peaks, fftconvolve


def estimate_period(x, dt):
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or not np.isfinite(x).all() or len(x)<4 or dt<=0: raise ValueError('invalid waveform')
    x=x-x.mean(); energy=np.dot(x,x)
    if energy<=1e-30: raise ValueError('no oscillation')
    corr=correlate(x,x,mode='full',method='fft')[len(x)-1:]/energy
    low=max(1,int(1/(500*dt))); high=min(len(x)-2,int(1/(30*dt)))
    peaks,_=find_peaks(corr[low:high+1]); peaks=peaks+low
    if len(peaks)==0: raise ValueError('no autocorrelation maximum')
    # Select earliest near-best ACF peak, retain alternatives for manual review.
    best=np.max(corr[peaks]); near=peaks[corr[peaks]>=max(.2,best*.85)]
    if len(near)==0: raise ValueError('weak autocorrelation')
    lag=int(near[0]); f=np.fft.rfftfreq(len(x),dt);power=np.abs(np.fft.rfft(x*np.hanning(len(x))))**2
    spectral_hz = float(f[np.argmax(power[1:])+1])
    ratio = spectral_hz * lag * dt
    return {'period_s':lag*dt,'acf':float(corr[lag]),'candidate_lags_samples':peaks.tolist(),'spectrum_peak_hz':spectral_hz,'acf_near_best_lags':near.tolist(), 'status':'AMBIGUOUS_PERIOD' if ratio < .8 or ratio > 1.2 else 'CANDIDATE', 'spectral_to_acf_frequency_ratio':ratio}


def track(x,dt,period_s):
    """Matched filter plus weak half-period exclusion; returns candidates and similarity.

    Fixed-template correlation can miss evolving pulses. Caller must assess confidence,
    half/double periods, peak selection, edge effects and differences from method A.
    """
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or not np.isfinite(x).all() or dt<=0 or period_s<=0: raise ValueError('invalid inputs')
    n=int(round(period_s/dt));
    if n<4 or len(x)<3*n: raise ValueError('insufficient samples')
    # Template is one full raw waveform around a strong maximum, avoiding end effects.
    centers,_=find_peaks(x,distance=max(1,n//2),prominence=(np.percentile(x,95)-np.percentile(x,5))*.1)
    centers=centers[(centers>=n//2)&(centers+n-n//2<len(x))]
    if len(centers)==0: raise ValueError('no usable template')
    center=int(centers[np.argmax(x[centers])]); template=x[center-n//2:center-n//2+n].copy(); template-=template.mean()
    norm=np.linalg.norm(template)
    if norm<=1e-15: raise ValueError('flat template')
    template/=norm
    dots=fftconvolve(x,template[::-1],mode='valid'); ones=np.ones(n); sums=fftconvolve(x,ones,mode='valid'); sq=fftconvolve(x*x,ones,mode='valid'); var=np.maximum(sq-sums*sums/n,0); denom=np.sqrt(var); quality=np.divide(dots,denom,out=np.zeros_like(dots),where=denom>1e-15)
    indices,_=find_peaks(quality,height=.5,distance=max(1,n//2))
    offsets=np.zeros(len(indices))
    for j,k in enumerate(indices):
        d=quality[k-1]-2*quality[k]+quality[k+1]
        if d<0: offsets[j]=np.clip(.5*(quality[k-1]-quality[k+1])/d,-.5,.5)
    times=(indices+offsets+n//2)*dt
    return {'times':times,'similarity':quality[indices], 'template_center_s':center*dt,'template':template}
