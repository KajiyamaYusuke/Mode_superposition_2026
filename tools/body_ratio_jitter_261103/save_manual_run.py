#!/usr/bin/env python3
"""Archive and analyze a finished manually run simulation without running it.

The output directory is resolved from --run-dir, which may be either the solver's
output directory or a case directory containing an output/ subdirectory.
"""
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import analyze_flow_jitter as flow
import template_period as template

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DESTINATION = ROOT / 'analysis_runs/manual_run_archive'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def plain(value):
    if isinstance(value, dict): return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [plain(v) for v in value]
    if isinstance(value, np.ndarray): return plain(value.tolist())
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer): return int(value)
    if isinstance(value, np.bool_): return bool(value)
    return value


def parse_param(path):
    lines = [s.strip() for s in path.read_text(encoding='utf-8-sig').splitlines()
             if s.strip() and not s.lstrip().startswith('#')]
    if len(lines) < 23: raise ValueError('parameter file has too few positional settings')
    nstep, dt, ps = int(lines[2]), float(lines[4]), float(lines[15])
    if nstep <= 0 or not np.isfinite(dt) or dt <= 0 or not np.isfinite(ps):
        raise ValueError('invalid nstep, dt or pressure')
    named = {}
    for line in lines[23:]:
        if '=' in line:
            key, value = line.split('=', 1)
            named[key.strip().lower()] = value.strip()
    return {'nmode': int(lines[0]), 'nsurfz': int(lines[1]), 'nstep': nstep,
            'nwrite': int(lines[3]), 'dt_s': dt, 'pressure_pa': ps,
            'zeta_left': float(lines[5]), 'zeta_right': float(lines[6]),
            'contact_coefficients': [float(v) for v in lines[7].split()],
            'ncont': int(lines[8]), 'named': named}


def body_token(path):
    hits = re.findall(r'(?<![A-Za-z0-9])b(\d+)(?:p(\d+))?c(\d+)(?:p(\d+))?(?![A-Za-z0-9])', Path(path).stem, flags=re.I)
    if len(hits) != 1: return None
    integer, fraction, cover, cover_fraction = hits[0]
    return {'body_kpa': str(Decimal(integer + ('.' + fraction if fraction else ''))),
            'cover_kpa': str(Decimal(cover + ('.' + cover_fraction if cover_fraction else '')))}


def input_metadata(param, manifest_path, original_param):
    mapping = {'left_frequency': ('leftfrequencyfile', 'left_frequency'),
               'right_frequency': ('rightfrequencyfile', 'right_frequency'),
               'left_mode': ('leftmodefile', 'left_mode_vtu'),
               'right_mode': ('rightmodefile', 'right_mode_vtu'),
               'left_surface': ('leftsurfacenasfile', 'left_surface_nas'),
               'right_surface': ('rightsurfacenasfile', 'right_surface_nas')}
    manifest = {}
    if manifest_path and manifest_path.exists():
        for line in manifest_path.read_text().splitlines():
            if '=' in line:
                k, v = line.split('=', 1); manifest[k.strip()] = v.strip().strip('"')
    result = {}
    for label, (param_key, manifest_key) in mapping.items():
        raw = manifest.get(manifest_key) or param['named'].get(param_key)
        if not raw: continue
        path = Path(raw).expanduser()
        # The copied params_used.txt changes the base for relative paths; only
        # use its original source location if explicitly supplied by the user.
        if not path.is_absolute() and original_param is not None:
            path = (original_param.parent / path).resolve()
        result[label] = {'recorded_path': raw, 'sha256_now': sha(path) if path.is_absolute() and path.is_file() else None,
                         'hash_meaning': 'current file, not independently proven at simulation time'}
    left = body_token(result.get('left_mode', {}).get('recorded_path', ''))
    right = body_token(result.get('right_mode', {}).get('recorded_path', ''))
    if left and right:
        result['body_token_inference'] = {'left': left, 'right': right,
            'right_over_left': float(Decimal(right['body_kpa']) / Decimal(left['body_kpa'])) if Decimal(left['body_kpa']) else None,
            'note': 'inferred from filenames, not measured material properties'}
    return result


def resolve_run(path):
    path = path.expanduser().resolve()
    if (path / 'dat').is_dir(): return path
    if (path / 'output' / 'dat').is_dir(): return path / 'output'
    raise ValueError(f'{path}: expected dat/ or output/dat/')


def stable_copy(src, dst):
    before = sha(src); shutil.copy2(src, dst); copied = sha(dst); after = sha(src)
    if before != copied or before != after:
        raise RuntimeError(f'file changed during copy: {src}; archive the finished run instead')
    return before


def load_signals(dat, nstep, dt):
    expected = np.arange(0, nstep, 5)
    def read(name, columns):
        path = dat / name
        data = np.loadtxt(path, comments='#', ndmin=2)
        if data.shape != (len(expected), columns) or not np.isfinite(data).all():
            raise ValueError(f'{path}: incomplete, malformed, or nonfinite')
        return data
    with (dat / 'airflow_vt.dat').open() as f: header = f.readline()
    if not header.startswith('#') or 'airflow' not in header.lower(): raise ValueError('invalid flow header')
    ug = read('airflow_vt.dat', 2)
    if not np.array_equal(ug[:, 0], expected): raise ValueError('flow step indices are not the complete 5-step sequence')
    displacement = read('displace.dat', 3)
    time = expected * dt
    if not np.allclose(displacement[:, 0], time, rtol=0, atol=max(1e-9, dt * 1e-4)):
        raise ValueError('displacement time axis does not match dt and flow samples')
    result = {'time_s': time, 'flow_m3_s': ug[:, 1], 'uy_left_mm': displacement[:, 1], 'uy_right_mm': displacement[:, 2]}
    pressure_path = dat / 'pressure_vt.dat'
    if pressure_path.exists():
        pressure = read('pressure_vt.dat', 2)
        if not np.array_equal(pressure[:, 0], expected): raise ValueError('acoustic pressure step indices differ')
        result['acoustic_pressure_pa'] = pressure[:, 1]
    area_path = dat / 'area.dat'
    if area_path.exists():
        area = np.loadtxt(area_path, comments='#', ndmin=2)
        if area.shape[0] != len(expected) or area.shape[1] < 4 or not np.array_equal(area[:, 0], expected):
            raise ValueError('area output is incomplete or has invalid steps')
        # Invalid cuts may be NaN; use finite interior sections only.
        interior = area[:, 2:-1]; finite = np.isfinite(interior)
        minimum = np.min(np.where(finite, interior, np.inf), axis=1)
        result['minimum_area_mm2'] = np.where(np.isfinite(minimum), minimum, np.nan)
    return result


def amplitude(values):
    p5, p95 = np.percentile(values, [5, 95])
    centered = values - np.mean(values)
    return {'half_peak_to_peak_mm': float(np.ptp(values) / 2),
            'half_p95_minus_p5_mm': float((p95-p5)/2),
            'ac_rms_mm': float(np.sqrt(np.mean(centered**2))),
            'mean_mm': float(np.mean(values))}


def analyze_jitter(values, dt, start, end):
    selected = values[int(round(start/dt)):int(round(end/dt))]
    result = flow.analyze(values, dt=dt, start=start, end=end)
    candidates = result.pop('candidates', {})
    report = {'status': 'UNDEFINED_JITTER', 'reason': result.get('reason', result['status']),
              'maxima_qc': result, 'method': None, 'metrics': None,
              'provisional': True, 'note': 'fixed template thresholds were selected after inspection of a previous representative case'}
    tracks = {}
    try:
        estimate = template.estimate_period(selected, dt)
        report['period_estimate'] = estimate
        report['periodicity'] = flow.periodicity(selected, dt, estimate['period_s'])
        for multiplier in (.5, 1, 2):
            try:
                track = template.track(selected, dt, estimate['period_s'] * multiplier)
                times = track['times'] + start; periods = np.diff(times)
                tracks[str(multiplier)] = {'times': times, 'similarity': track['similarity'],
                                           'periods': periods, 'metrics': flow.period_metrics(periods)}
            except ValueError as exc: tracks[str(multiplier)] = {'error': str(exc)}
        main, half = tracks['1'], tracks['0.5']; a = candidates.get('A_0.10_1')
        if estimate['status'] == 'CANDIDATE' and a and 'periods' in main and len(main['periods']) >= 10 and result.get('status') in ('VALID_STATIONARY', 'VALID_WINDOW_DEPENDENT'):
            ma, mb = a['metrics'], main['metrics']
            if abs(np.mean(a['periods']) / np.mean(main['periods']) - 1) < .05 and abs(ma['local'] - mb['local']) < max(.2, .25 * max(ma['local'], mb['local'])):
                report.update(status=result['status'], reason='', method='A_local_maximum', metrics=ma)
        if report['status'] == 'UNDEFINED_JITTER' and estimate['status'] == 'CANDIDATE' and 'times' in main and 'times' in half:
            t = main['times']; other = half['times']; periods = main['periods']; po = half['periods']
            if len(periods) >= 20 and len(po) >= 20:
                med = float(np.median(periods)); mt = main['metrics']; mh = half['metrics']; evidence = report['periodicity']
                match = float(np.mean(np.min(abs(t[:, None]-other[None, :]), axis=1) < .2*estimate['period_s']))
                quality = (abs(med/estimate['period_s']-1) <= .05 and np.all((periods >= .5*med) & (periods <= 1.5*med)) and
                           np.quantile(main['similarity'], .05) >= .8 and match >= .9 and
                           abs(mt['local']-mh['local']) <= max(.2, .3*max(mt['local'], mh['local'])) and
                           np.isfinite(evidence['acf_at_period']) and evidence['acf_at_period'] >= .5 and
                           np.isfinite(evidence['spectral_ratio']) and evidence['spectral_ratio'] >= .05)
                middle = (start+end)/2
                p1 = flow.window_periods(t, start, middle); p2 = flow.window_periods(t, middle, end)
                if quality and min(len(p1), len(p2)) >= 10:
                    j1, j2 = flow.period_metrics(p1)['local'], flow.period_metrics(p2)['local']
                    a1, a2 = [float(np.percentile(v, 95)-np.percentile(v, 5)) for v in np.array_split(selected, 2)]
                    changing = (abs(np.mean(p1)-np.mean(p2))/med > .05 or abs(a1-a2)/max(a1,a2) > .2 or
                                abs(j1-j2) > max(.2,.5*max(j1,j2)))
                    report.update(status='VALID_WINDOW_DEPENDENT' if changing else 'VALID_STATIONARY',
                                  reason='', method='B_template', metrics=mt,
                                  template_scale_event_match=match)
        report['template_candidates'] = {k: {'cycles': len(v.get('periods', [])), 'metrics': v.get('metrics'),
                                                   'error': v.get('error')} for k, v in tracks.items()}
    except ValueError as exc:
        report['reason'] = str(exc)
    report['maxima_candidates'] = {k: {'cycles': len(v['periods']), 'metrics': v['metrics']} for k, v in candidates.items()}
    return report, candidates, tracks


def save_tables(out, signal, start, end, candidates, tracks):
    keys = list(signal)
    with (out / 'waveforms.csv').open('x', newline='') as f:
        w = csv.writer(f); w.writerow(keys)
        w.writerows(zip(*(signal[k] for k in keys)))
    with (out / 'flow_events.csv').open('x', newline='') as f:
        w = csv.writer(f); w.writerow(['signal', 'method', 'time_s', 'flow_m3_s', 'template_similarity'])
        for key, entry in candidates.items():
            for t in entry['times']:
                if start <= t < end:
                    w.writerow(['Ug', key, t, np.interp(t, signal['time_s'], signal['flow_m3_s']), ''])
        for key, entry in tracks.items():
            if 'times' not in entry: continue
            for t, similarity in zip(entry['times'], entry['similarity']):
                if start <= t < end:
                    w.writerow(['Ug', 'template_' + key, t, np.interp(t, signal['time_s'], signal['flow_m3_s']), similarity])
    with (out / 'flow_cycles.csv').open('x', newline='') as f:
        w = csv.writer(f); w.writerow(['signal', 'method', 'start_s', 'end_s', 'period_s'])
        for key, entry in list(candidates.items()) + [('template_' + k, v) for k, v in tracks.items()]:
            if 'times' not in entry: continue
            times = entry['times'][(entry['times'] >= start) & (entry['times'] < end)]
            for a, b in zip(times, times[1:]): w.writerow(['Ug', key, a, b, b-a])


def save_figures(out, signal, start, end, candidates, tracks):
    t = signal['time_s']; ug = signal['flow_m3_s']; left = signal['uy_left_mm']; right = signal['uy_right_mm']
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    ax[0].plot(t, ug*1e6, lw=.65); ax[0].set_ylabel('Ug [mL/s]')
    ax[1].plot(t, left, lw=.65, label='Left uy'); ax[1].plot(t, right, lw=.65, label='Right uy'); ax[1].set_ylabel('Monitor uy [mm]'); ax[1].legend()
    if 'acoustic_pressure_pa' in signal: ax[2].plot(t, signal['acoustic_pressure_pa'], lw=.65); ax[2].set_ylabel('Outlet pressure [Pa]')
    elif 'minimum_area_mm2' in signal: ax[2].plot(t, signal['minimum_area_mm2'], lw=.65); ax[2].set_ylabel('Minimum area [mm²]')
    else: ax[2].text(.5, .5, 'No pressure or area output', transform=ax[2].transAxes, ha='center')
    for a in ax: a.axvspan(start, end, alpha=.10, color='orange'); a.grid(alpha=.15)
    ax[-1].set_xlabel('Time [s]'); fig.tight_layout(); fig.savefig(out / 'waveforms.png', dpi=150); plt.close(fig)
    zoom_start = max(start, end-.05); mask = (t >= zoom_start) & (t < end)
    fig, ax = plt.subplots(figsize=(11, 3.8)); ax.plot(t[mask], ug[mask]*1e6, color='black', lw=.8)
    for label, entry, marker in [('A: maximum', candidates.get('A_0.10_1', {}), 'o'), ('B: template', tracks.get('1', {}), 'x')]:
        if 'times' not in entry: continue
        points = entry['times']; points = points[(points >= zoom_start) & (points < end)]
        if len(points): ax.scatter(points, np.interp(points, t, ug)*1e6, marker=marker, s=22, label=label)
    ax.set(xlabel='Time [s]', ylabel='Ug [mL/s]', title='Candidate events (not necessarily valid periods)'); ax.grid(alpha=.15)
    if ax.get_legend_handles_labels()[0]: ax.legend()
    fig.tight_layout(); fig.savefig(out / 'flow_events.png', dpi=150); plt.close(fig)
    selection = ug[(t >= start) & (t < end)]; centered = selection-np.mean(selection)
    power = np.abs(np.fft.rfft(centered*np.hanning(len(centered))))**2
    freq = np.fft.rfftfreq(len(centered), np.median(np.diff(t)))
    fig, ax = plt.subplots(figsize=(9, 3.5)); ax.plot(freq, power, lw=.7); ax.set(xlim=(0, min(1000, freq[-1])), xlabel='Frequency [Hz]', ylabel='Power [arbitrary]', title='Ug spectrum: analysis window'); ax.grid(alpha=.15); fig.tight_layout(); fig.savefig(out / 'flow_spectrum.png', dpi=150); plt.close(fig)


def append_index(destination, row):
    fields = ['archive_id', 'created_utc', 'label', 'pressure_pa', 'body_ratio_from_filename', 'window_start_s', 'window_end_s',
              'jitter_status', 'jitter_method', 'jitter_local_percent', 'rap_percent', 'ppq5_percent', 'period_cv_percent',
              'amplitude_left_mm', 'amplitude_right_mm', 'amplitude_ratio_R_over_L', 'archive_path']
    with (destination / '.index.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = destination / 'summary.csv'; exists = path.exists()
        if exists:
            with path.open(newline='') as f:
                if next(csv.reader(f), None) != fields: raise ValueError('existing summary.csv has different columns')
        with path.open('a', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            if not exists: w.writeheader()
            w.writerow({k: row.get(k, '') for k in fields}); f.flush(); os.fsync(f.fileno())
        fcntl.flock(lock, fcntl.LOCK_UN)


def archive_run(run_dir, destination, param_path=None, label='', start=.5, end=1.):
    run = resolve_run(run_dir)
    recorded = run / 'txt' / 'params_used.txt'
    supplied = param_path.expanduser().resolve() if param_path else None
    if supplied and not supplied.is_file(): raise ValueError(f'parameter file not found: {supplied}')
    if recorded.is_file():
        if supplied and sha(supplied) != sha(recorded):
            raise ValueError('supplied --param differs from recorded params_used.txt')
        source_param = recorded; provenance = 'solver_recorded_params_used'
    elif supplied:
        source_param = supplied; provenance = 'user_supplied_not_independently_verified'
    else: raise ValueError('no txt/params_used.txt; supply --param for an explicitly unverified parameter file')
    parent_status = run.parent / 'status.json' if run.name == 'output' else run / 'status.json'
    if parent_status.is_file():
        saved_status = json.loads(parent_status.read_text())
        state = saved_status.get('simulation_status', saved_status.get('status'))
        if state in ('RUNNING', 'INTERRUPTED', 'RUN_FAILED', 'simulation_failed', 'running'):
            raise ValueError(f'source run status is {state}; incomplete runs are not archived as complete')
    else: saved_status = None
    if not np.isfinite(start) or not np.isfinite(end) or start < 0 or end <= start:
        raise ValueError('analysis window must satisfy 0 <= start < end')
    if label and not re.fullmatch(r'[A-Za-z0-9_.-]+', label): raise ValueError('label may contain only letters, digits, dot, dash or underscore')
    destination = destination.expanduser().resolve(); destination.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc); archive_id = stamp.strftime('%Y%m%dT%H%M%S') + '_' + uuid.uuid4().hex[:8]
    temp = Path(tempfile.mkdtemp(prefix='.building_', dir=destination)); final = destination / ((label + '_' if label else '') + archive_id)
    published = False
    try:
        snapshot = temp / 'snapshot'; (snapshot / 'dat').mkdir(parents=True); (snapshot / 'txt').mkdir()
        sources = {'txt/params_used.txt': source_param, 'dat/airflow_vt.dat': run/'dat/airflow_vt.dat', 'dat/displace.dat': run/'dat/displace.dat'}
        for name in ('pressure_vt.dat', 'area.dat', 'displace_xy.dat', 'pressure.dat'):
            path = run/'dat'/name
            if path.is_file(): sources['dat/'+name] = path
        manifest = run/'txt/manifest.txt'
        if manifest.is_file(): sources['txt/manifest.txt'] = manifest
        if parent_status.is_file(): sources['txt/source_status.json'] = parent_status
        hashes = {}
        for relative, source in sources.items():
            if not source.is_file(): raise ValueError(f'required output missing: {source}')
            hashes[relative] = stable_copy(source, snapshot/relative)
        for relative, source in sources.items():
            if sha(source) != hashes[relative]: raise RuntimeError(f'source changed while archiving: {source}')
        if saved_status:
            expected_param_hash = saved_status.get('parameter_sha256')
            expected_flow_hash = saved_status.get('validation', {}).get('flow_sha256')
            if expected_param_hash and expected_param_hash != hashes['txt/params_used.txt']:
                raise ValueError('recorded status parameter hash does not match params_used.txt')
            if expected_flow_hash and expected_flow_hash != hashes['dat/airflow_vt.dat']:
                raise ValueError('recorded status flow hash does not match output')
        params = parse_param(snapshot/'txt/params_used.txt'); dt = params['dt_s']; nstep = params['nstep']; sample_dt = 5*dt
        if end > nstep*dt + 1e-12: raise ValueError('analysis window extends beyond simulation duration')
        if any(abs(v/sample_dt-round(v/sample_dt)) > 1e-7 for v in (start,end)):
            raise ValueError('window boundaries must align with the 5-step flow sampling interval')
        signal = load_signals(snapshot/'dat', nstep, dt)
        indices = np.arange(len(signal['time_s']))
        mask = (indices >= round(start/sample_dt)) & (indices < round(end/sample_dt))
        if np.count_nonzero(mask) < 10: raise ValueError('analysis window has fewer than 10 samples')
        left = amplitude(signal['uy_left_mm'][mask]); right = amplitude(signal['uy_right_mm'][mask])
        left_amp = left['half_peak_to_peak_mm']; right_amp = right['half_peak_to_peak_mm']
        ratio = right_amp/left_amp if left_amp > 1e-12 else None
        jitter, candidates, tracks = analyze_jitter(signal['flow_m3_s'], sample_dt, start, end)
        metadata = input_metadata(params, snapshot/'txt/manifest.txt', supplied)
        save_tables(temp, signal, start, end, candidates, tracks); save_figures(temp, signal, start, end, candidates, tracks)
        created = stamp.isoformat(); metric = jitter.get('metrics') or {}
        report = {'archive_id': archive_id, 'label': label, 'created_utc': created, 'source_run': str(run),
                  'parameter_provenance': provenance, 'source_status': saved_status,
                  'completion_check': 'all expected step samples present; solver exit status not independently known' if not saved_status else 'all expected step samples present; source status recorded',
                  'window_s': [start, end], 'window_samples': int(np.count_nonzero(mask)), 'parameter_settings': params, 'input_metadata': metadata,
                  'snapshot_sha256': hashes,
                  'analysis_code_sha256': {p.name: sha(p) for p in [Path(__file__), Path(flow.__file__), Path(template.__file__)]},
                  'jitter': jitter, 'amplitude': {'definition': 'monitor displacement half peak-to-peak in same half-open window; R/L',
                                                 'left': left, 'right': right, 'right_over_left': ratio,
                                                 'warning': 'monitor-point amplitude, not full-surface amplitude'},
                  'notes': ['Jitter uses Ug; candidate event CSVs are not automatically validated.',
                            'Input-file hashes are current hashes and may not establish the files used by an older run.',
                            'Template thresholds are provisional; do not interpret high jitter as proof of chaos.']}
        (temp/'metrics.json').write_text(json.dumps(plain(report), ensure_ascii=False, indent=2, allow_nan=False)+'\n')
        metrics = {'archive_id': archive_id, 'created_utc': created, 'label': label, 'pressure_pa': params['pressure_pa'],
                   'body_ratio_from_filename': metadata.get('body_token_inference', {}).get('right_over_left'),
                   'window_start_s': start, 'window_end_s': end, 'jitter_status': jitter['status'], 'jitter_method': jitter['method'],
                   'jitter_local_percent': metric.get('local'), 'rap_percent': metric.get('rap'), 'ppq5_percent': metric.get('ppq5'),
                   'period_cv_percent': metric.get('cv'), 'amplitude_left_mm': left_amp, 'amplitude_right_mm': right_amp,
                   'amplitude_ratio_R_over_L': ratio, 'archive_path': str(final)}
        with (temp/'summary_row.csv').open('x', newline='') as f:
            w=csv.DictWriter(f, fieldnames=list(metrics));w.writeheader();w.writerow(metrics)
        temp.rename(final); published = True
        append_index(destination, metrics)
        return final, report
    finally:
        if not published and temp.exists(): shutil.rmtree(temp)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir', type=Path, default=ROOT/'output', help='solver output directory or case directory (default: repository output/)')
    p.add_argument('--param', type=Path, help='original parameter file; checked against recorded params_used.txt when present')
    p.add_argument('--save-dir', type=Path, default=DEFAULT_DESTINATION, help='archive directory')
    p.add_argument('--label', default='', help='optional identifier: ASCII letters, digits, dot, dash, underscore')
    p.add_argument('--window', type=float, nargs=2, metavar=('START_S','END_S'), default=(.5,1.), help='half-open analysis window in seconds (default: 0.50 1.00)')
    args=p.parse_args()
    try: final, report=archive_run(args.run_dir,args.save_dir,args.param,args.label,*args.window)
    except (ValueError, OSError, RuntimeError, json.JSONDecodeError) as exc: p.exit(2, f'ERROR: {exc}\n')
    print('Saved archive:', final)
    print('Jitter (provisional):', report['jitter']['status'], report['jitter']['reason'])
    print('Monitor displacement amplitude ratio R/L:', report['amplitude']['right_over_left'])
    print('Parameter provenance:', report['parameter_provenance'])

if __name__=='__main__': main()
