#!/usr/bin/env python3
"""Fixed-grid airflow threshold jitter analysis; see README_simple_fixed_period_jitter.md."""
import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]


def read_samples(path, dt, area=False):
    """Read uniformly spaced solver steps and values in their declared units."""
    path = Path(path)
    with path.open(encoding='utf-8') as f:
        header = ''.join(line for line in f if line.lstrip().startswith('#')).lower()
    unit = ('mm^2' if 'mm^2' in header else 'm^2' if 'm^2' in header else '') if area else ('m^3/s' if 'm^3/s' in header else '')
    if not unit:
        raise ValueError(f'{path}: missing or unsupported unit in header')
    data = np.loadtxt(path, comments='#', ndmin=2)
    if len(data) < 3 or (data.shape[1] < 50 if area else data.shape[1] != 2):
        raise ValueError(f'{path}: expected at least 3 rows and {"50 area columns" if area else "2 flow columns"}')
    steps = data[:, 0]
    dif = np.diff(steps)
    if not np.isfinite(steps).all() or np.any(steps < 0) or np.any(steps != np.rint(steps)) or np.any(dif <= 0) or np.any(dif != dif[0]):
        raise ValueError(f'{path}: invalid or nonuniform integer steps')
    values = data[:, 2:50] if area else data[:, 1]
    if area:
        if np.isinf(values).any() or np.any(values[np.isfinite(values)] < 0) or np.any(~np.isfinite(values).any(axis=1)):
            raise ValueError(f'{path}: invalid area values')
        values = np.min(np.where(np.isfinite(values), values, np.inf), axis=1) * (1e6 if unit == 'm^2' else 1)
    elif not np.isfinite(values).all():
        raise ValueError(f'{path}: nonfinite flow')
    return steps * dt, values, int(dif[0])


def estimate_initial_period(time, flow, start, end, fmin, fmax):
    """Choose the strongest positive local ACF peak in the frequency range."""
    x = flow[(time >= start) & (time < min(end, start + .10))]
    if len(x) < 3:
        raise ValueError('insufficient initial waveform')
    x = x - x.mean()
    acf = signal.correlate(x, x, mode='full', method='fft')[len(x)-1:]
    if acf[0] <= np.finfo(float).tiny:
        raise ValueError('constant initial waveform')
    acf /= acf[0]
    dt = time[1] - time[0]
    lo, hi = max(1, math.ceil(1 / fmax / dt)), min(len(acf)-1, math.floor(1 / fmin / dt))
    peaks, _ = signal.find_peaks(acf)
    peaks = peaks[(peaks >= lo) & (peaks <= hi) & (acf[peaks] > 0)]
    if not len(peaks):
        raise ValueError('no positive ACF peak in frequency range')
    best = peaks[np.argmax(acf[peaks])]
    return best * dt, [(float(p * dt), float(acf[p])) for p in peaks], float(acf[best])


def find_initial_valleys(time, flow, start, period):
    """Choose the first prominent valley and the lowest valleys near the next two periods."""
    indices = np.flatnonzero((time >= start) & (time < start + 4 * period))
    if len(indices) < 3:
        raise ValueError('insufficient samples for initial valleys')
    x = flow[indices]
    amplitude = np.percentile(x, 95) - np.percentile(x, 5)
    peaks, _ = signal.find_peaks(-x, prominence=.1 * amplitude)
    candidates = indices[peaks]
    if amplitude <= 0 or not len(candidates):
        return [], candidates
    chosen = [candidates[0]]
    for k in (1, 2):
        near = candidates[(time[candidates] > time[chosen[-1]]) & (np.abs(time[candidates] - (time[chosen[0]] + k * period)) <= .35 * period)]
        if not len(near):
            return [], candidates
        chosen.append(near[np.argmin(flow[near])])
    return chosen, candidates


def make_fixed_bins(start, end, anchor, period):
    """Intersect the immutable anchor + k * period grid with the evaluation window."""
    first = math.floor((start - anchor) / period)
    last = math.ceil((end - anchor) / period)
    return [(k, max(start, anchor + k*period), min(end, anchor + (k+1)*period),
             anchor + k*period, anchor + (k+1)*period) for k in range(first, last)
            if min(end, anchor + (k+1)*period) > max(start, anchor + k*period)]


def peak_candidates(time, flow, period):
    """Find global local maxima with prominence and both base rises >= 10% of local P95-P5."""
    peaks, props = signal.find_peaks(flow, prominence=0)
    accepted = []
    for i, p in enumerate(peaks):
        left = np.searchsorted(time, time[p] - 1.5*period)
        right = np.searchsorted(time, time[p] + 1.5*period, side='right')
        local = flow[left:right]
        amp = np.percentile(local, 95) - np.percentile(local, 5)
        if amp > 0 and props['prominences'][i] >= .1*amp and flow[p]-flow[props['left_bases'][i]] >= .1*amp and flow[p]-flow[props['right_bases'][i]] >= .1*amp:
            accepted.append(p)
    return np.asarray(accepted, dtype=int)


def find_downcross(time, flow, peak, stop, fraction):
    """Last linearly interpolated descent before first zero or observed branch minimum."""
    stop_index = np.searchsorted(time, stop, side='left')
    if stop_index <= peak + 1 or flow[peak] <= 0:
        return None, 0
    branch = flow[peak:stop_index]
    zeros = np.flatnonzero(branch <= 0)
    finish = int(zeros[0]) if len(zeros) else int(np.argmin(branch))
    threshold = flow[peak] * fraction
    x = branch[:finish+1]
    crossings = np.flatnonzero((x[:-1] > threshold) & (x[1:] <= threshold))
    if not len(crossings):
        return None, 0
    i = peak + crossings[-1]
    t = time[i] + (time[i+1]-time[i]) * (flow[i]-threshold)/(flow[i]-flow[i+1])
    return float(t), len(crossings)


def period_metrics(events):
    """Calculate metrics using every adjacent event interval."""
    p = np.diff(np.asarray(events, dtype=float))
    n = len(p)
    mean = float(np.mean(p)) if n else None
    result = dict(intervals=n, jitter_terms=max(0,n-1), rap_terms=max(0,n-2), ppq5_terms=max(0,n-4),
                  jitter_pct=None, rap_pct=None, ppq5_pct=None, mean_interval_ms=mean*1000 if n else None,
                  min_interval_ms=float(np.min(p)*1000) if n else None, max_interval_ms=float(np.max(p)*1000) if n else None)
    if n >= 2:
        result['jitter_pct'] = float(100*np.mean(np.abs(np.diff(p)))/mean)
    if n >= 3:
        result['rap_pct'] = float(100*np.mean(np.abs(p[1:-1]-np.convolve(p, np.ones(3)/3, 'valid')))/mean)
    if n >= 5:
        result['ppq5_pct'] = float(100*np.mean(np.abs(p[2:-2]-np.convolve(p, np.ones(5)/5, 'valid')))/mean)
    return result, p


def plot_case(path, time, flow, bins, events, start, end):
    """Show the full evaluation waveform and all fixed boundaries and events."""
    fig, ax = plt.subplots(figsize=(18, 5))
    mask = (time >= start) & (time < end)
    ax.plot(time[mask], flow[mask], color='0.25', lw=.7, label='flow')
    for b in bins:
        ax.axvline(b['boundary_time_s'], color='tab:red' if 'BOUNDARY_NEAR_PEAK' in b['flags'] else '0.75', lw=.5)
    for key, color, marker in [('main_peak_time_s','tab:blue','o'), ('event_3pct_s','tab:orange','v'), ('event_5pct_s','tab:green','^')]:
        ts = [e[key] for e in events if e[key] is not None]
        ax.scatter(ts, np.interp(ts, time, flow), s=12, c=color, marker=marker, label=key, zorder=3)
    ax.set(xlim=(start,end), xlabel='Time [s]', ylabel='Airflow [m³/s]', title=path.stem + ' (red: boundary near peak)')
    ax.legend(loc='upper right', fontsize=8); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def analyze_case(case, opt, figures):
    """Analyze one case; area is a separate, non-selecting diagnostic."""
    cid = case['case_id']
    time, flow, stride = read_samples(case['flow_path'], opt.dt)
    if time[0] > opt.start or time[-1] + stride*opt.dt < opt.end - 1e-12:
        raise ValueError('flow does not cover evaluation window')
    acf, candidates, acf_value = estimate_initial_period(time, flow, opt.start, opt.end, opt.fmin, opt.fmax)
    valleys, valley_candidates = find_initial_valleys(time, flow, opt.start, acf)
    init = dict(T_acf_ms=1000*acf, acf_selected=acf_value,
                acf_candidates=';'.join(f'{t:.9g}:{v:.6g}' for t,v in candidates),
                valley_candidates_s=';'.join(f'{time[i]:.9g}' for i in valley_candidates), stride=stride)
    if not valleys:
        return dict(status='INIT_FAILED', error='three initial valleys not found', init=init, bins=[], events=[], area=[])
    anchor = float(time[valleys[0]])
    period = float((time[valleys[2]]-anchor)/2)
    init.update(T_ref_ms=1000*period, anchor_s=anchor, valley1_s=float(time[valleys[1]]), valley2_s=float(time[valleys[2]]), T_ref_minus_acf_ms=1000*(period-acf))
    peaks = peak_candidates(time, flow, period)
    grid = make_fixed_bins(opt.start, opt.end, anchor, period)
    bins, mains = [], []
    for k, left, right, boundary, next_boundary in grid:
        inside = peaks[(time[peaks] >= left) & (time[peaks] < right)]
        main = inside[np.argmax(flow[inside])] if len(inside) else None
        nearby = peaks[np.abs(time[peaks]-boundary) <= .1*period]
        flags = []
        if left > boundary + 1e-12 or right < next_boundary - 1e-12: flags.append('EDGE_PARTIAL')
        if main is None: flags.append('NO_MAIN_PEAK')
        if len(nearby): flags.append('BOUNDARY_NEAR_PEAK')
        row = dict(case_id=cid, bin_id=k, boundary_time_s=boundary, bin_start_s=left, bin_end_s=right,
                   next_boundary_s=next_boundary, main_peak_count=int(main is not None), secondary_count=max(0,len(inside)-1),
                   secondary_times_s=';'.join(f'{time[p]:.9g}' for p in inside if p != main),
                   main_peak_time_s=float(time[main]) if main is not None else None,
                   main_flow=float(flow[main]) if main is not None else None,
                   peak_to_boundary_distance_s=float(min(time[main]-boundary,next_boundary-time[main])) if main is not None else None,
                   Ug_at_boundary=float(np.interp(boundary,time,flow)) if time[0] <= boundary <= time[-1] else None, flags=';'.join(flags))
        bins.append(row)
        if main is not None: mains.append((main, row))
    events = []
    for i, (p, b) in enumerate(mains):
        stop = float(time[mains[i+1][0]]) if i+1 < len(mains) else opt.end
        row = dict(case_id=cid, main_peak_id=i+1, bin_id=b['bin_id'], main_peak_time_s=float(time[p]), main_flow=float(flow[p]), branch_end_s=stop)
        for pct in (3, 5):
            t, count = find_downcross(time, flow, p, stop, pct/100)
            row[f'event_{pct}pct_s'], row[f'crossings_{pct}pct'] = t, count
            if t is None: b['flags'] += (';' if b['flags'] else '') + f'NO_{pct}PCT'
        events.append(row)
    area_rows, area_error = [], ''
    if case.get('area_path'):
        try:
            at, av, _ = read_samples(case['area_path'], opt.dt, area=True)
            if at[0] > opt.start or at[-1] + (at[1]-at[0]) < opt.end - 1e-12:
                raise ValueError('area does not cover evaluation window')
            for e in events:
                indices = np.flatnonzero((at >= e['main_peak_time_s']) & (at < e['branch_end_s']))
                if len(indices):
                    j = indices[np.argmin(av[indices])]
                    row = dict(case_id=cid, main_peak_id=e['main_peak_id'], area_min_time_s=float(at[j]), area_min_mm2=float(av[j]))
                    for pct in (3, 5):
                        t = e[f'event_{pct}pct_s']
                        row[f'area_minus_{pct}pct_ms'] = (float(at[j])-t)*1000 if t is not None else None
                    area_rows.append(row)
        except (ValueError, OSError) as exc:
            area_error = str(exc)
    plot_case(figures / f'{cid}.png', time, flow, bins, events, opt.start, opt.end)
    return dict(status='OK', error='', init=init, bins=bins, events=events, area=area_rows, area_error=area_error)


def write_csv(path, rows, fields):
    """Write stable headers even when a table is empty."""
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def run_batch(cases, opt):
    """Run cases independently and write compact, auditable tables."""
    out = Path(opt.out_dir); out.mkdir(parents=True, exist_ok=False)
    figures = out / 'figures'; figures.mkdir()
    summaries, all_bins, all_events, all_periods, all_area = [], [], [], [], []
    for case in cases:
        cid = case['case_id']
        try:
            result = analyze_case(case, opt, figures)
        except (ValueError, OSError) as exc:
            result = dict(status='FAILED', error=str(exc), init={}, bins=[], events=[], area=[])
        all_bins.extend(result['bins']); all_events.extend(result['events']); all_area.extend(result['area'])
        for pct in (3, 5):
            times = sorted(e[f'event_{pct}pct_s'] for e in result['events'] if e[f'event_{pct}pct_s'] is not None)
            metrics, periods = period_metrics(times)
            flags = [b['flags'] for b in result['bins']]
            row = dict(case_id=cid, stiffness_ratio=case.get('stiffness_ratio',''), pressure_pa=case.get('pressure_pa',''), threshold_pct=pct,
                       status=result['status'], events=len(times), flagged_bins=sum('BOUNDARY_NEAR_PEAK' in f for f in flags),
                       edge_partial=sum('EDGE_PARTIAL' in f for f in flags), no_main_peak=sum('NO_MAIN_PEAK' in f for f in flags),
                       no_3pct=sum('NO_3PCT' in f for f in flags), no_5pct=sum('NO_5PCT' in f for f in flags),
                       error=result['error'], area_error=result.get('area_error',''), **result['init'], **metrics)
            row['mean_period_warning'] = bool(len(periods) and result['init'].get('T_ref_ms') and abs(metrics['mean_interval_ms']/result['init']['T_ref_ms']-1) > .1)
            summaries.append(row)
            for i, p in enumerate(periods):
                all_periods.append(dict(case_id=cid, threshold_pct=pct, interval_id=i+1, start_s=times[i], end_s=times[i+1], interval_s=p,
                                        jitter_contribution_pct=100*abs(periods[i+1]-p)/np.mean(periods) if i+1 < len(periods) else None))
    write_csv(out/'summary.csv', summaries, 'case_id stiffness_ratio pressure_pa threshold_pct status events intervals jitter_terms rap_terms ppq5_terms jitter_pct rap_pct ppq5_pct mean_interval_ms min_interval_ms max_interval_ms T_ref_ms anchor_s valley1_s valley2_s T_acf_ms T_ref_minus_acf_ms acf_selected acf_candidates valley_candidates_s stride flagged_bins edge_partial no_main_peak no_3pct no_5pct mean_period_warning error area_error'.split())
    write_csv(out/'events.csv', all_events, 'case_id main_peak_id bin_id main_peak_time_s main_flow branch_end_s event_3pct_s crossings_3pct event_5pct_s crossings_5pct'.split())
    write_csv(out/'bins.csv', all_bins, 'case_id bin_id boundary_time_s bin_start_s bin_end_s next_boundary_s main_peak_count secondary_count secondary_times_s main_peak_time_s main_flow peak_to_boundary_distance_s Ug_at_boundary flags'.split())
    write_csv(out/'periods.csv', all_periods, 'case_id threshold_pct interval_id start_s end_s interval_s jitter_contribution_pct'.split())
    if any(c.get('area_path') for c in cases):
        write_csv(out/'area_alignment.csv', all_area, 'case_id main_peak_id area_min_time_s area_min_mm2 area_minus_3pct_ms area_minus_5pct_ms'.split())
    plot_comparison(figures/'jitter_vs_ratio.png', summaries)
    (out/'README.md').write_text('Fixed boundaries and all observed events; see tools/README_simple_fixed_period_jitter.md for definitions and limitations.\n', encoding='utf-8')
    return summaries


def plot_comparison(path, rows):
    """Plot ratios grouped by pressure and threshold, when numeric ratios exist."""
    groups = {}
    for r in rows:
        try:
            ratio = float(r['stiffness_ratio'])
            if not math.isfinite(ratio) or r['jitter_pct'] is None: continue
        except (ValueError, TypeError): continue
        groups.setdefault((r['pressure_pa'] or 'unspecified', r['threshold_pct']), []).append((ratio, r['jitter_pct']))
    if not groups: return
    fig, ax = plt.subplots(figsize=(8,5))
    for (pressure, pct), points in sorted(groups.items()):
        points.sort(); ax.plot(*zip(*points), marker='o', label=f'{pressure} Pa, {pct}%')
    ax.set(xlabel='Stiffness ratio', ylabel='Jitter local [%]'); ax.legend(); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--manifest'); source.add_argument('--case', action='append', metavar='ID:PATH')
    parser.add_argument('--out-dir', required=True); parser.add_argument('--start', type=float, default=.5); parser.add_argument('--end', type=float, default=1.)
    parser.add_argument('--dt', type=float, default=1e-5); parser.add_argument('--fmin', type=float, default=60); parser.add_argument('--fmax', type=float, default=400)
    opt = parser.parse_args()
    if not all(math.isfinite(x) for x in (opt.start,opt.end,opt.dt,opt.fmin,opt.fmax)) or not (0 <= opt.start < opt.end and opt.dt > 0 and 0 < opt.fmin < opt.fmax):
        parser.error('invalid window, dt or frequency range')
    if opt.manifest:
        with open(opt.manifest, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or not {'case_id','flow_path'} <= set(reader.fieldnames): parser.error('manifest needs case_id and flow_path')
            cases = list(reader)
    else:
        try: cases = [dict(zip(('case_id','flow_path'), s.split(':',1))) for s in opt.case]
        except ValueError: parser.error('--case requires ID:PATH')
    if not cases or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', c.get('case_id','')) or c.get('case_id') in ('.','..') or not c.get('flow_path') for c in cases):
        parser.error('cases require a safe case_id and flow_path')
    if len({c['case_id'] for c in cases}) != len(cases): parser.error('case_id must be unique')
    for case in cases:
        for key in ('flow_path','area_path'):
            if case.get(key): case[key] = str((ROOT / Path(case[key]).expanduser()).resolve())
    try: summaries = run_batch(cases, opt)
    except FileExistsError: parser.error('out-dir already exists')
    for r in summaries: print(f"{r['case_id']} {r['threshold_pct']}% {r['status']}: events={r['events']} jitter={r['jitter_pct']} error={r['error']}")


if __name__ == '__main__':
    main()
