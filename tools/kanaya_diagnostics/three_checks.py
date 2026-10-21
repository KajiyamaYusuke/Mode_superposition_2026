#!/usr/bin/env python3
"""Read-only replay of an existing Kanaya run; never launches the solver.

Run: python3 -B tools/kanaya_diagnostics/three_checks.py
All source data are hashed before/after. An existing output directory is refused.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy import signal
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = Path('/home/kajiyama/code/simulation-main')
TOL = 1e-12
WINDOWS = {'all': (0, 19999), 'ramp': (0, 5000),
           'transient': (5000, 10000), 'early': (10000, 15000),
           'late': (15000, 19999), 'evaluation': (10000, 19999)}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def git_state(root):
    return {key: subprocess.check_output(['git', *args], cwd=root, text=True)
            for key, args in {'status': ['status', '--short'], 'branch': ['branch', '--show-current'],
                              'head': ['rev-parse', 'HEAD'], 'diff': ['diff', '--binary'],
                              'index_diff': ['diff', '--cached', '--binary']}.items()}


def metadata_tree(root):
    rows = []
    for base, dirs, names in os.walk(root):
        dirs.sort()
        for name in sorted(names):
            path = Path(base) / name
            st = path.stat()
            rows.append([str(path.relative_to(root)), st.st_size, st.st_mtime_ns])
    return {'files': len(rows), 'path_size_mtime_sha256':
            hashlib.sha256(json.dumps(rows).encode()).hexdigest()}


def original_snapshot():
    return {'git': git_state(ORIGINAL), 'trees': {name: metadata_tree(ORIGINAL / name)
            for name in ('input', 'output', 'result', 'build', 'tools/diagnostics')}}


def inventory(paths):
    result = {}
    for path in sorted(set(paths)):
        st = path.stat()
        result[str(path)] = {'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns,
                             'resolved_path': str(path.resolve()), 'sha256': sha(path)}
    return result


def close_gate(a, b, *, atol=1e-12, rtol=1e-8, label='numeric gate'):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.array_equal(np.isnan(a), np.isnan(b)):
        raise AssertionError(f'{label}: shape/NaN mismatch')
    mask = np.isfinite(a) & np.isfinite(b)
    if np.isinf(a).any() or np.isinf(b).any():
        raise AssertionError(f'{label}: infinity')
    delta = np.abs(a[mask] - b[mask])
    scaled = delta / (atol + rtol * np.maximum(np.abs(a[mask]), np.abs(b[mask])))
    if scaled.max(initial=0) > 1:
        raise AssertionError(f'{label}: max scaled residual {scaled.max():.6g}')
    return {'max_abs': float(delta.max(initial=0)),
            'max_scaled': float(scaled.max(initial=0)), 'atol': atol, 'rtol': rtol,
            'matching_nan': int(np.isnan(a).sum()), 'passed': True}


def rounded_gate(reconstructed, observed, label):
    """C++ defaultfloat precision 6: half of the last printed significant digit."""
    a, b = np.asarray(reconstructed), np.asarray(observed)
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        raise AssertionError(f'{label}: NaN mask mismatch')
    finite = np.isfinite(a)
    scale = np.maximum(np.abs(a[finite]), np.abs(b[finite]))
    quantum = 10. ** (np.floor(np.log10(np.maximum(scale, 1e-300))) - 5)
    tolerance = .50001 * quantum + 2e-10
    error = np.abs(a[finite] - b[finite])
    if np.any(error > tolerance):
        raise AssertionError(f'{label}: rounded mismatch {np.max(error/tolerance):.6g}')
    return {'passed': True, 'max_abs': float(error.max(initial=0)),
            'max_rounding_fraction': float((error/tolerance).max(initial=0)),
            'tolerance': '0.50001*10^(floor(log10(abs(value)))-5)+2e-10; six significant digits',
            'matching_nan': int(np.isnan(a).sum())}


def positive_integrals(h0, h1, length):
    if not np.isfinite([h0, h1, length]).all():
        raise ValueError('nonfinite gap, not a closed section')
    if h0 <= 0 and h1 <= 0:
        return 0., 0., 0.
    if h0 < 0:
        length *= h1 / (h1 - h0)
        h0 = 0.
    elif h1 < 0:
        length *= h0 / (h0 - h1)
        h1 = 0.
    return (.5*(h0+h1)*length,
            .25*(h0**3+h0*h0*h1+h0*h1*h1+h1**3)*length, length)


def section_segments(surface, plane, maximum_y):
    """surface[i,j,xyz]; do not bridge invalid i edges or missing AP rows."""
    intersections = []
    for j in range(surface.shape[1]):
        candidates = []
        for i in range(surface.shape[0]-1):
            a, b = surface[i, j], surface[i+1, j]
            if not np.isfinite([a, b]).all():
                continue
            if plane < min(a[0], b[0])-TOL or plane > max(a[0], b[0])+TOL:
                continue
            dx = b[0]-a[0]
            if abs(dx) <= TOL:
                new = [a, b] if abs(plane-a[0]) <= TOL and abs(plane-b[0]) <= TOL else []
            else:
                new = [a + np.clip((plane-a[0])/dx, 0., 1.)*(b-a)]
            for p in new:
                if not any(np.hypot(*(p[1:]-q[1:])) <= TOL for q in candidates):
                    candidates.append(p)
        intersections.append((max if maximum_y else min)(candidates, key=lambda p: p[1])
                             if candidates else None)
    last = None
    segments = []
    for j, p in enumerate(intersections):
        if p is not None:
            if last is not None and p[2] < last[2]-TOL:
                raise AssertionError('AP z ordering reversed')
            last = p
            if j and intersections[j-1] is not None and p[2] > intersections[j-1][2]:
                segments.append((intersections[j-1], p))
    return segments


def integrate_segments(left, right, gap_sign=1):
    a = g3 = common = opened = 0.
    il = ir = 0
    profile = []
    while il < len(left) and ir < len(right):
        la, lb = left[il]
        ra, rb = right[ir]
        low, high = max(la[2], ra[2]), min(lb[2], rb[2])
        if high > low:
            yl = la[1] + (np.array([low, high])-la[2])/(lb[2]-la[2])*(lb[1]-la[1])
            yr = ra[1] + (np.array([low, high])-ra[2])/(rb[2]-ra[2])*(rb[1]-ra[1])
            h = gap_sign*(yr-yl)
            da, dg, ds = positive_integrals(h[0], h[1], high-low)
            a += da; g3 += dg; opened += ds; common += high-low
            profile.append([low, high, *yl, *yr, *h])
        if lb[2] <= rb[2]:
            il += 1
        if rb[2] <= lb[2]:
            ir += 1
    return {'A': a if common > 0 else np.nan, 'G3': g3 if common > 0 else np.nan,
            'S_common': common, 'S_open': opened, 'profile': profile}


def evaluate_sections(left, right, count=50, gap_sign=1):
    low = max(np.nanmin(left[..., 0]), np.nanmin(right[..., 0]))
    high = min(np.nanmax(left[..., 0]), np.nanmax(right[..., 0]))
    if not high > low:
        raise AssertionError('no whole-surface common x range')
    x = low + (high-low)*np.arange(count)/(count-1)
    result = [integrate_segments(section_segments(left, v, gap_sign > 0),
                                 section_segments(right, v, gap_sign < 0), gap_sign) for v in x]
    return x, result


def load_modes(path, mass, nmode):
    """Stream ASCII VTU, matching frequency-key order and (raw Z,X,Y) mapping."""
    groups, active, points = {}, {}, None
    for _, element in ET.iterparse(path, events=('end',)):
        if element.tag != 'DataArray':
            continue
        name = element.get('Name', '')
        match = re.match(r'Displacement_field,_([XYZ])-?component_@_(.*)', name)
        if match:
            component, raw = match.groups()
            if component == 'X':
                key, index = raw, 2
                while key in groups:
                    key = raw + '_v' + str(index); index += 1
                active[raw] = key
                groups[key] = {}
            groups[active[raw]][component] = np.fromstring(element.text or '', sep=' ')
        elif name not in ('connectivity', 'offsets', 'types') and element.get('NumberOfComponents') == '3':
            points = np.fromstring(element.text or '', sep=' ').reshape(-1, 3)[:, [2, 0, 1]].copy()
        element.clear()
    if points is None:
        raise ValueError('VTU Points not found')
    keys = sorted(groups, key=lambda k: (float(k.split('_Hz')[0]), k))[:nmode]
    modes = np.stack([np.stack([groups[k][c] for c in ('Z', 'X', 'Y')], axis=-1) for k in keys])
    assert modes.shape == (nmode, len(points), 3)
    modes *= (1/np.sqrt(mass/len(points)*np.sum(modes*modes, axis=(1, 2))))[:, None, None]
    return points, modes, keys


def surface_from_nas(path, points):
    grid, quads = {}, []
    for line in path.read_text().splitlines():
        fields = line.strip().split(',')
        if fields[0] == 'GRID':
            grid[int(fields[1])] = [float(fields[5]), float(fields[3]), float(fields[4])]
        elif fields[0] == 'CQUAD4':
            quads.append([int(v) for v in fields[3:7]])
    tree = cKDTree(points)
    mapped = {}
    for key in sorted({v for q in quads for v in q}):
        candidates = sorted(tree.query_ball_point(grid[key], 1e-3))
        candidates = [c for c in candidates if np.sum((points[c]-grid[key])**2) < 1e-6]
        if not candidates:
            raise AssertionError('NAS/VTU mapping failed')
        mapped[key] = candidates[0]  # first qualifying VTU ID, as in Geometry.cpp
    quads = [[mapped[v] for v in q] for q in quads]
    unique = sorted({v for q in quads for v in q})
    clusters = []
    for z in sorted(points[unique, 2]):
        if clusters and z-clusters[-1][-1] <= 1e-3:
            clusters[-1].append(z)
        else:
            clusters.append([z])
    levels = np.array([np.mean(c) for c in clusters])
    row_of = {p: int(np.argmin(np.abs(levels-points[p, 2]))) for p in unique}
    assert all(abs(levels[row_of[p]]-points[p, 2]) <= 1e-2 for p in unique)
    adj = {p: set() for p in unique}
    for q in quads:
        for j in sorted({row_of[p] for p in q}):
            ids = [p for p in q if row_of[p] == j]
            assert len(ids) == 2, 'unexpected quad z topology'
            a, b = ids; adj[a].add(b); adj[b].add(a)
    rows = []
    for j in range(len(levels)):
        ids = [p for p in unique if row_of[p] == j]
        ends = [p for p in ids if len(adj[p]) <= 1]
        # This mesh is a single chain per AP row. Other topologies stop rather
        # than guessing the C++ unordered-component traversal.
        assert len(ends) == 2 and max(len(adj[p]) for p in ids) <= 2
        cur = min(ends, key=lambda p: points[p, 0])
        ordered, seen = [], set()
        while cur is not None:
            ordered.append(cur); seen.add(cur)
            remaining = sorted(adj[cur]-seen)
            cur = remaining[0] if remaining else None
        assert len(ordered) == len(ids), 'disconnected AP row requires exact traversal helper'
        if rows:
            a, b, c, d = points[[ordered[0], ordered[-1], rows[-1][0], rows[-1][-1]]]
            if np.sum((a-d)**2)+np.sum((b-c)**2) < np.sum((a-c)**2)+np.sum((b-d)**2):
                ordered.reverse()
        rows.append(ordered)
    surface = np.full((max(map(len, rows)), len(rows)), -1, dtype=int)
    for j, row in enumerate(rows):
        surface[:len(row), j] = row
    assert (surface >= 0).all(), 'missing nodes require masked geometry replay'
    maximum_i = np.unravel_index(np.argmax(points[surface, 0]), surface.shape)[0]
    return surface, maximum_i+1


def displace(points, modes, q):
    return points + np.einsum('m,mnc->nc', q, modes, optimize=True)*1e3


def work_decomposition(force, dq):
    total = np.sum(force*dq, axis=0)
    mean = np.mean(force, axis=0)*np.sum(dq, axis=0)
    fluctuating = np.sum((force-np.mean(force, axis=0))*dq, axis=0)
    close_gate(total, mean+fluctuating, atol=1e-18, rtol=1e-10, label='exact work decomposition')
    return total, mean, fluctuating


def static_prediction(q, v, force, omega, damping, dt):
    qb, vb, fb = q.mean(axis=0), v.mean(axis=0), force.mean(axis=0)
    acc = np.gradient(v, dt, axis=0, edge_order=2).mean(axis=0)
    residual = omega**2*qb-fb
    return qb, vb, acc, fb/omega**2, residual, residual+acc+2*damping*omega*vb


def phase_diagnostic(q, v, force, dt):
    def spectrum(values):
        ft = np.fft.rfft(signal.detrend(values)*np.hanning(len(values)))
        power = np.abs(ft)**2
        power[0] = 0
        peak = int(np.argmax(power))
        share = float(power[max(1, peak-1):peak+2].sum()/max(power.sum(), 1e-300))
        return ft, peak, share
    fq, kq, sq = spectrum(q)
    fv, kv, sv = spectrum(v)
    ff, kf, sf = spectrum(force)
    frequencies = np.fft.rfftfreq(len(q), dt)
    # Two-segment Welch coherence is only descriptive for this short record.
    _, coh = signal.coherence(v, force, fs=1/dt, nperseg=max(16, len(q)//3))
    ck = int(round(frequencies[kv]/(1/dt/max(16, len(q)//3))))
    coherence = float(coh[min(ck, len(coh)-1)])
    accepted = abs(kv-kf) <= 1 and min(sv, sf) >= .65 and coherence >= .7
    corr = float(np.corrcoef(v, force)[0, 1]) if np.std(v)*np.std(force) > 0 else np.nan
    return {'q_AC_RMS': float(np.std(q)), 'qdot_AC_RMS': float(np.std(v)),
            'Q_AC_RMS': float(np.std(force)), 'Q_qdot_correlation': corr,
            'q_peak_Hz': frequencies[kq], 'qdot_peak_Hz': frequencies[kv],
            'Q_peak_Hz': frequencies[kf], 'q_spectral_share': sq, 'qdot_spectral_share': sv,
            'Q_spectral_share': sf, 'coherence_descriptive': coherence,
            'phase_accepted': accepted,
            'phase_deg': float(np.angle(ff[kv]*fv[kv].conjugate(), deg=True)) if accepted else np.nan,
            'phase_reason': 'matched concentrated spectra; short-record descriptive' if accepted else
                            'multipeak/mismatched band/weak coherence: no scalar phase'}


def json_write(path, data):
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(w) for k, w in v.items()}
        if isinstance(v, (list, tuple)): return [clean(w) for w in v]
        if isinstance(v, np.ndarray): return clean(v.tolist())
        if isinstance(v, np.generic): return clean(v.item())
        if isinstance(v, float) and not np.isfinite(v): return None
        return v
    with path.open('x') as stream:
        json.dump(clean(data), stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def read_run(root):
    base = root/'output/kanaya_diagnostic'
    required = ['modal_native.bin', 'run_metadata.json', 'flow_geometry.csv',
                'modal_energy_summary.csv', 'modal_selected_timeseries.csv',
                'contact_regions.csv', 'force_regions.csv', 'consistency.csv']
    missing = [v for v in required if not (base/v).is_file()]
    if missing:
        raise FileNotFoundError(f'missing_inputs: {missing}')
    meta = json.loads((base/'run_metadata.json').read_text())
    assert meta['complete'] and meta['integrated_intervals'] == 19999
    assert meta['excluded_intervals'] == 1 and meta['excluded_step'] == 19999
    assert meta['nmode'] == [100, 100] and abs(meta['dt_s']-1e-5) < 1e-18
    assert abs(meta['actual_ramp_time_s']-.05) < 1e-15
    assert meta['modal_binary_fields'] == ['Q_pos', 'Q_neg', 'Q_fluid', 'Q_contact_adopted', 'q_n', 'qdot_n', 'delta_q']
    n, nm, dt = meta['integrated_intervals'], meta['nmode'][0], meta['dt_s']
    assert (base/'modal_native.bin').stat().st_size == n*2*nm*7*8
    raw = np.memmap(base/'modal_native.bin', dtype='<f8', mode='r', shape=(n, 2, nm, 7))
    for start in range(0, n, 1000):
        assert np.isfinite(raw[start:start+1000]).all(), 'nonfinite native data'
    consistency = pd.read_csv(base/'consistency.csv')
    valid = consistency[consistency.included == 1]
    assert len(valid) == 2*n and valid.valid.eq(1).all() and valid.scaled_max.max() <= 1
    assert consistency[consistency.included == 0].step.to_list() == [n, n]
    flow = pd.read_csv(base/'flow_geometry.csv')
    assert np.array_equal(flow.step, np.arange(n+1))
    assert np.array_equal(flow.interval_included, np.r_[np.ones(n), 0])
    close_gate(flow.time_s, flow.step.to_numpy()*dt, label='flow time')
    summary = pd.read_csv(base/'modal_energy_summary.csv')
    manifest = dict(line.split(' = ', 1) for line in (root/'output/txt/manifest.txt').read_text().splitlines() if ' = ' in line)
    frequencies = []
    for side in ('left', 'right'):
        table = np.loadtxt(manifest[side+'_frequency'], comments='%')[:nm]
        frequencies.append(table[:, 0])
    omega = 2*np.pi*np.array(frequencies)
    zeta = np.array([float(manifest['zetaL']), float(manifest['zetaR'])])
    validation = {'metadata_complete': True, 'missing_inputs': [],
                  'consistency_scaled_max': float(valid.scaled_max.max()), 'native_bytes': raw.nbytes,
                  'Q_split': close_gate(raw[..., 0]+raw[..., 1], raw[..., 2], label='Q split'),
                  'q_continuity': close_gate(raw[:-1, ..., 4]+raw[:-1, ..., 6], raw[1:, ..., 4],
                                             atol=1e-18, rtol=1e-10, label='q continuity')}
    summary_checks = []
    for s, side in enumerate(('L', 'R')):
        for window in ('all', 'ramp', 'transient', 'evaluation'):
            a, b = WINDOWS[window]
            target = summary[(summary.side == side) & (summary.window == window)].sort_values('mode')
            assert len(target) == nm and target.interval_count.eq(b-a).all()
            close_gate(target.f_Hz, frequencies[s], label='frequencies')
            for col, key in enumerate(('W_pos', 'W_neg', 'W_fluid', 'W_contact_reference')):
                measured = np.sum(raw[a:b, s, :, col]*raw[a:b, s, :, 6], axis=0)
                result = close_gate(measured, target[key], atol=2e-18, rtol=1e-9, label=f'{side}/{window}/{key}')
                summary_checks.append({'side': side, 'window': window, 'field': key, **result})
            damping = np.sum(raw[a:b, s, :, 5]**2*2*zeta[s]*omega[s]*dt, axis=0)
            summary_checks.append({'side': side, 'window': window, 'field': 'D_struct',
                **close_gate(damping, target.D_struct, atol=2e-18, rtol=1e-9, label='stored damping')})
    validation['stored_summary_checks'] = summary_checks
    dat = {name: np.loadtxt(root/f'output/dat/{name}.dat') for name in ('area', 'section_x', 'gap_cubed', 'pressure', 'separation')}
    steps = dat['area'][:, 0].astype(int)
    assert np.array_equal(steps, np.arange(0, n, 5))
    for key, values in dat.items():
        assert np.array_equal(values[:, 0], steps)
        assert values.shape[1] == (9 if key == 'separation' else 51)
    minimum = np.nanmin(dat['area'][:, 2:-1], axis=1)  # interior sections 1..48
    idx = np.nanargmin(dat['area'][:, 2:-1], axis=1)+1
    validation['dat_minimum_vs_flow'] = rounded_gate(flow.Amin_mm2.to_numpy()[steps], minimum, 'DAT minimum')
    validation['dat_minimum_x_vs_flow'] = rounded_gate(flow.xMin_mm.to_numpy()[steps], dat['section_x'][np.arange(len(steps)), idx+1], 'DAT minimum x')
    validation['dat_sep_vs_flow'] = rounded_gate(flow.xSep_mm.to_numpy()[steps], dat['separation'][:, 2], 'DAT separation')
    # Confirm exact source-stage dynamics without performing a time advance:
    # recover recorded acceleration from the velocity history, gamma=.55.
    accelerations = np.zeros((n, 2, nm))
    for k in range(n-1):
        accelerations[k+1] = ((raw[k+1, ..., 5]-raw[k, ..., 5])/dt-.45*accelerations[k])/.55
    residual = accelerations[1:]+2*zeta[:, None]*omega*raw[1:, ..., 5]+omega**2*raw[1:, ..., 4]
    validation['recorded_newmark_equation'] = close_gate(residual, raw[:-1, ..., 2]+raw[:-1, ..., 3],
                                                        atol=2e-10, rtol=1e-7, label='recorded Newmark equation')
    prior = json.loads((base/'manifest.json').read_text())
    source_hashes = {name: sha(root/name) for name in prior['source_sha256']}
    assert source_hashes == prior['source_sha256'], 'current model sources differ from recorded run'
    provenance = json.loads((base/'validation/production_provenance.json').read_text())
    for entry in provenance['inputs']:
        assert sha(Path(entry['path'])) == entry['sha256'], 'run input provenance mismatch'
    assert sha(root/'output/txt/params_used.txt') == sha(root/'input/kanaya_param.txt')
    validation['source_and_input_run_provenance'] = True
    contact = pd.read_csv(base/'contact_regions.csv')
    assert len(contact) == 4*n and contact.step.max() == n-1
    assert (contact.groupby('step').size() == 4).all()
    forces = pd.read_csv(base/'force_regions.csv')
    assert forces.wrong_normal_count.eq(0).all()
    evaluation_forces = forces[(forces.side == 'L') & (forces.time_s >= .1)]
    count = evaluation_forces.step.nunique()
    force_means = evaluation_forces.groupby(['pressure_sign', 'zone'])[['opening_Fy_N', 'area_m2']].sum()/count
    validation['decimated_force_regions_mean'] = [dict(pressure_sign=int(k[0]), zone=k[1],
        opening_force_N=float(v.opening_Fy_N), applied_area_mm2=float(v.area_m2*1e6)) for k, v in force_means.iterrows()]
    for entry in meta['force_statistics']:
        if entry['side'] == 'L' and entry['window'] == 3:
            measured = evaluation_forces[evaluation_forces.pressure_sign == entry['pressure_sign']].opening_Fy_N.sum()/count
            assert abs(measured-entry['mean_opening_N']) < .001*abs(entry['mean_opening_N'])
    validation['force_region_decimation_vs_full_mean_within_0p1percent'] = True
    selected_checks = []
    for chunk in pd.read_csv(base/'modal_selected_timeseries.csv', chunksize=40000):
        indices = chunk.step.to_numpy(int)
        sides = chunk.side.map({'L': 0, 'R': 1}).to_numpy(int)
        modes = chunk['mode'].to_numpy(int)-1
        for k, field in enumerate(['Q_pos', 'Q_neg', 'Q_fluid', 'Q_contact_adopted_predicted_stage', 'q_n', 'qdot_n', 'delta_q']):
            result = close_gate(raw[indices, sides, modes, k], chunk[field].to_numpy(),
                                atol=1e-18, rtol=1e-9, label='shared selected timeseries')
            selected_checks.append(result['max_abs'])
    validation['selected_timeseries_max_abs'] = max(selected_checks)
    return raw, meta, flow, summary, manifest, omega, zeta, dat, validation


def geometry_setup(root, manifest, raw, meta, validation):
    text = (root/'output/txt/params_used.txt').read_text().splitlines()
    mass = None
    for i, line in enumerate(text):
        if 'mass of vocal fold' in line:
            mass = float(text[i+1])
    assert mass is not None
    points, modes, surfaces, nxs = [], [], [], []
    cache = {}
    for side in ('left', 'right'):
        key = manifest[side+'_mode_vtu']
        if key not in cache:
            cache[key] = load_modes(Path(key), mass, raw.shape[2])
        p, phi, _ = cache[key]
        sf, nx = surface_from_nas(Path(manifest[side+'_surface_nas']), p)
        points.append(p.copy()); modes.append(phi.copy()); surfaces.append(sf); nxs.append(nx)
    ymid = np.max(points[1][:, 1])
    points[1][:, 1] = 2*ymid-points[1][:, 1]
    modes[1][:, :, 1] *= -1
    ni, nj = min(nxs), min(s.shape[1] for s in surfaces)
    il, ir = surfaces[0][:ni, :nj], surfaces[1][:ni, :nj]
    gaps = points[1][ir, 1]-points[0][il, 1]
    sign = -1 if gaps.sum() < 0 else 1
    correction = float(manifest['target_initial_minimum_gap_mm'])-np.min(sign*gaps)
    points[0][:, 1] -= .5*sign*correction
    points[1][:, 1] += .5*sign*correction
    assert sign == meta['gapSign']
    exported = pd.read_csv(root/'output/csv/surfp_output.csv', comment='#', header=None,
                           names=['i', 'j', 'id', 'x', 'y', 'z'], skipinitialspace=True)
    assert np.array_equal(surfaces[0][exported.i, exported.j], exported.id)
    gram = mass/len(points[0]) * (modes[0].reshape(raw.shape[2], -1) @ modes[0].reshape(raw.shape[2], -1).T)
    for side, m1, m2, value in meta['gram']:
        close_gate(np.array([gram[m1-1, m2-1]]), np.array([value]), atol=1e-13, rtol=1e-10, label='Gram')
    geometry = {'npoints': [len(p) for p in points], 'nsurfl': [s.shape[0] for s in surfaces],
                'nsurfz': [s.shape[1] for s in surfaces], 'nxsup': nxs, 'mass_kg': mass,
                'gap_sign': sign, 'initial_gap_translation_per_side_mm': .5*correction,
                'full_gram_max_off_diagonal': float(np.max(np.abs(gram-np.diag(np.diag(gram))))),
                'full_gram_max_diagonal_error': float(np.max(np.abs(np.diag(gram)-1))),
                'surface_NAS_export_ID_match': True}
    # Independent VTU replay gate, including an explicit wrong-time control.
    vtu_checks = []
    for k in (0, 125, 249):
        path = root/f'output/result/deform_combined{k:04d}.vtu'
        values = ET.parse(path).find('.//Points/DataArray')
        obs = np.fromstring(values.text, sep=' ').reshape(-1, 3)
        n = 15000+20*k
        correct = np.concatenate([displace(points[s], modes[s], raw[n+1, s, :, 4]) for s in range(2)])
        wrong = np.concatenate([displace(points[s], modes[s], raw[n, s, :, 4]) for s in range(2)])
        gate = close_gate(correct, obs, atol=2e-11, rtol=1e-12, label='VTU q_(n+1)')
        wrong_error = float(np.max(np.abs(wrong-obs)))
        assert wrong_error > 1e-7, 'wrong-time control did not discriminate'
        vtu_checks.append({'file': str(path), 'solver_step': n, 'coordinate_step': n+1,
                           'wrong_q_n_max_abs_mm': wrong_error, **gate})
    validation['geometry'] = geometry
    validation['VTU_time_offset'] = vtu_checks
    return points, modes, surfaces, nxs, sign


def ap_analysis(raw, flow, dat, geom, validation):
    points, modes, surfaces, nxs, sign = geom
    steps = dat['area'][:, 0].astype(int)
    ev = flow.iloc[10000:19999]
    desired = [0, 5, 5000, 10000, 12500, 15000, 17500, 19500]
    for key in ('Amin_mm2', 'xSep_mm'):
        desired += [int(ev[key].idxmin()), int(ev[key].idxmax())]
    desired += [11000, 14000, 18000, 19990]
    chosen = sorted({int(steps[np.argmin(np.abs(steps-v))]) for v in desired})
    rows, frames, gates = [], {}, []
    for step in chosen:
        coords = [displace(points[s], modes[s], raw[step, s, :, 4]) for s in range(2)]
        x, sections = evaluate_sections(coords[0][surfaces[0][:nxs[0]]],
                                        coords[1][surfaces[1][:nxs[1]]], gap_sign=sign)
        ar = np.array([r['A'] for r in sections]); gr = np.array([r['G3'] for r in sections])
        d = np.flatnonzero(steps == step)[0]
        gates.append({'step': step, 'area': rounded_gate(ar, dat['area'][d, 1:], 'area replay'),
                      'x': rounded_gate(x, dat['section_x'][d, 1:], 'x replay'),
                      'G3': rounded_gate(gr, dat['gap_cubed'][d, 1:], 'G3 replay'),
                      'Amin_flow': close_gate(np.array([np.nanmin(ar[1:-1])]),
                          np.array([flow.Amin_mm2.iloc[step]]), atol=1e-9, rtol=1e-9, label='Amin replay')})
        for k, r in enumerate(sections):
            valid = r['S_common'] > 0
            rows.append({'step': step, 'time_s': flow.time_s.iloc[step], 'plane_index': k, 'x': x[k],
                'A_observed': dat['area'][d, k+1], 'A_reconstructed': ar[k],
                'G3_observed': dat['gap_cubed'][d, k+1], 'G3_reconstructed': gr[k],
                'common_span_mm': r['S_common'], 'open_span_mm': r['S_open'],
                'mean_gap_common_mm': ar[k]/r['S_common'] if valid else np.nan,
                'mean_gap_open_mm': ar[k]/r['S_open'] if r['S_open'] > 0 else np.nan,
                'valid': valid, 'area_error': ar[k]-dat['area'][d, k+1],
                'reason': 'common AP intersection' if valid else 'no common AP interval; undefined, not closure'})
        frames[step] = (x, sections, coords)
    validation['AP_reconstruction'] = gates
    coverage = pd.DataFrame(rows)
    minimum_rows, comparisons = [], []
    for step, group in coverage.groupby('step'):
        interior = group[(group.plane_index > 0) & (group.plane_index < 49) & group.valid]
        minimum = interior.loc[interior.A_reconstructed.idxmin()]
        minimum_rows.append(minimum.to_dict())
        # Width-only comparison uses nearest interior plane with >=95% of
        # this frame's maximum common span, not an assumed uniform gap.
        maxspan = interior.common_span_mm.max()
        nearfull = interior[(interior.common_span_mm >= .995*maxspan)
                            & (interior.plane_index <= minimum.plane_index-2)]
        if nearfull.empty:
            nearfull = interior[(interior.common_span_mm >= .995*maxspan)
                                & (interior.plane_index != minimum.plane_index)]
        assert not nearfull.empty, 'no distinct near-full-width reference plane'
        ref = nearfull.loc[(nearfull.x-minimum.x).abs().idxmin()]
        comparisons.append({'step': int(step), 'minimum_plane': int(minimum.plane_index),
            'minimum_x_mm': minimum.x, 'minimum_area_mm2': minimum.A_reconstructed,
            'common_span_mm': minimum.common_span_mm, 'open_span_mm': minimum.open_span_mm,
            'common_span_fraction_of_frame_max': minimum.common_span_mm/maxspan,
            'reference_plane': int(ref.plane_index), 'reference_span_mm': ref.common_span_mm,
            'reference_area_mm2': ref.A_reconstructed,
            'span_ratio': minimum.common_span_mm/ref.common_span_mm,
            'open_span_ratio': minimum.open_span_mm/ref.open_span_mm if ref.open_span_mm > 0 else np.nan,
            'gap_ratio': minimum.mean_gap_common_mm/ref.mean_gap_common_mm if ref.mean_gap_common_mm > 0 else np.nan,
            'area_ratio': minimum.A_reconstructed/ref.A_reconstructed,
            'reference_to_min_log_width_fraction':
                np.log(minimum.common_span_mm/ref.common_span_mm)/np.log(minimum.A_reconstructed/ref.A_reconstructed)
                if minimum.A_reconstructed > 0 and minimum.A_reconstructed != ref.A_reconstructed else np.nan})
    validation['AP_minimum_comparisons'] = comparisons
    validation['AP_selected_frames'] = chosen
    validation['AP_all_DAT_undefined_counts'] = np.isnan(dat['area'][:, 1:]).sum(axis=0).tolist()
    validation['AP_width_selection_sensitivity'] = ap_selection_sensitivity(coverage)
    return coverage, frames


def ap_selection_sensitivity(coverage):
    """Compare rankings only; never replace measured areas or run the flow model."""
    result = []
    for step, group in coverage.groupby('step'):
        valid = group[(group.plane_index > 0) & (group.plane_index < 49)
                      & group.valid & (group.common_span_mm > 0)]
        a = valid.loc[valid.A_reconstructed.idxmin()]
        h = valid.loc[valid.mean_gap_common_mm.idxmin()]
        result.append({'step': int(step), 'area_min_plane': int(a.plane_index),
            'gap_min_plane': int(h.plane_index), 'different_plane': bool(a.plane_index != h.plane_index),
            'area_min_x_mm': float(a.x), 'gap_min_x_mm': float(h.x),
            'x_shift_mm': float(h.x-a.x), 'area_min_A_mm2': float(a.A_reconstructed),
            'gap_min_A_mm2': float(h.A_reconstructed), 'area_min_span_mm': float(a.common_span_mm),
            'gap_min_span_mm': float(h.common_span_mm),
            'area_min_mean_gap_mm': float(a.mean_gap_common_mm), 'minimum_mean_gap_mm': float(h.mean_gap_common_mm)})
    return result


def modal_analysis(raw, meta, flow, omega, zeta, validation):
    dt = meta['dt_s']
    static, dynamic = [], []
    for s, side in enumerate(('L', 'R')):
        for name, (a, b) in WINDOWS.items():
            data = raw[a:b, s]
            qp, qn, force, contact, q, v, dq = [data[..., k] for k in range(7)]
            qb, vb, acc, qs, residual, corrected = static_prediction(q, v, force+contact, omega[s], zeta[s], dt)
            wt, wm, wf = work_decomposition(force, dq)
            wpos, wneg, wc = [np.sum(x*dq, axis=0) for x in (qp, qn, contact)]
            close_gate(wt, wpos+wneg, atol=2e-18, rtol=1e-10, label='work sign split')
            dtotal = np.sum(2*zeta[s]*omega[s]*v*v*dt, axis=0)
            dac = np.sum(2*zeta[s]*omega[s]*(v-vb)**2*dt, axis=0)
            for m in range(raw.shape[2]):
                scale = max(abs(omega[s, m]**2*qb[m]), abs((force[:, m]+contact[:, m]).mean()))
                static.append({'side': side, 'mode': m+1, 'window': name, 'start_s': a*dt,
                    'end_s': b*dt, 'count': b-a, 'f_Hz': omega[s, m]/(2*np.pi),
                    'mean_Qfluid': force[:, m].mean(), 'mean_Qcontact': contact[:, m].mean(),
                    'mean_abs_Q': np.abs(force[:, m]+contact[:, m]).mean(),
                    'mean_q': qb[m], 'mean_abs_q': np.abs(q[:, m]).mean(),
                    'mean_qdot': vb[m], 'qdot_RMS_AC': np.std(v[:, m]),
                    'mean_qddot_finite_difference': acc[m], 'qstatic': qs[m],
                    'qstatic_fluid': force[:, m].mean()/omega[s, m]**2,
                    'qstatic_contact': contact[:, m].mean()/omega[s, m]**2,
                    'static_residual': residual[m], 'static_residual_scale': scale,
                    'static_residual_normalized': residual[m]/scale if scale > 1e-8 else np.nan,
                    'motion_corrected_residual': corrected[m],
                    'motion_corrected_normalized': corrected[m]/scale if scale > 1e-8 else np.nan})
                stats = phase_diagnostic(q[:, m], v[:, m], force[:, m], dt)
                relations = {}
                for column in ('Amin_mm2', 'Ug_m3_s', 'xSep_mm'):
                    values = flow[column].to_numpy()[a:b]
                    relations[f'{column}_q_correlation'] = float(np.corrcoef(values, q[:, m])[0, 1])
                    relations[f'{column}_qdot_correlation'] = float(np.corrcoef(values, v[:, m])[0, 1])
                dynamic.append({'side': side, 'mode': m+1, 'window': name, 'start_s': a*dt,
                    'end_s': b*dt, 'duration_s': (b-a)*dt, 'count': b-a,
                    'Wtotal': wt[m], 'Wmean_force': wm[m], 'Wfluctuating_force': wf[m],
                    'Wpos': wpos[m], 'Wneg': wneg[m], 'Wcontact_reference': wc[m],
                    'D_total': dtotal[m], 'D_ac': dac[m],
                    'Wfluctuating_minus_Dac_reference': wf[m]-dac[m],
                    'fluctuating_power': wf[m]/((b-a)*dt), **stats, **relations})
    static, dynamic = pd.DataFrame(static), pd.DataFrame(dynamic)
    additive = ['Wtotal', 'Wpos', 'Wneg', 'Wcontact_reference', 'D_total']
    for side in ('L', 'R'):
        pieces = [dynamic[(dynamic.side == side) & (dynamic.window == w)][additive].to_numpy() for w in ('early', 'late', 'evaluation')]
        close_gate(pieces[0]+pieces[1], pieces[2], atol=2e-18, rtol=1e-10, label='50ms work additivity')
    selected = {}
    for side in ('L', 'R'):
        sta = static[(static.side == side) & (static.window == 'evaluation')]
        dyn = dynamic[(dynamic.side == side) & (dynamic.window == 'evaluation')]
        selected[side] = sorted(set([1, 4, 5, 7, 35] + sta.nlargest(5, 'mean_abs_Q')['mode'].to_list()
            + sta.nlargest(5, 'mean_abs_q')['mode'].to_list()
            + dyn.assign(absW=dyn.Wfluctuating_force.abs()).nlargest(5, 'absW')['mode'].to_list()))
    validation['selected_modes'] = selected
    validation['50ms_split_additivity'] = True
    return static, dynamic, selected


def trend_and_cycles(raw, omega, zeta, dt, selected):
    rows, cycles, equal = [], [], []
    for s, side in enumerate(('L', 'R')):
        for m1 in range(1, raw.shape[2]+1):
            m = m1-1
            force, q, v, dq = [raw[:, s, m, k] for k in (2, 4, 5, 6)]
            # q_(19999) is the endpoint of the last INCLUDED interval.
            # It is not the excluded solver interval n=19999.
            qend = np.r_[q, q[-1]+dq[-1]]
            for cutoff in (20., 30., 40.):
                sos = signal.butter(4, cutoff, btype='highpass', fs=1/dt, output='sos')
                qhigh = signal.sosfiltfilt(sos, qend)
                fhigh = signal.sosfiltfilt(sos, force)
                vhigh = signal.sosfiltfilt(sos, v)
                for label, a, b in [('evaluation_full_edges', 10000, 19999),
                                     ('common_trim_25ms', 10000, 17499)]:
                    sl = slice(a, b)
                    w = float(np.sum(fhigh[sl]*np.diff(qhigh)[sl]))
                    dac = float(np.sum(2*zeta[s]*omega[s, m]*vhigh[sl]**2*dt))
                    rows.append({'side': side, 'mode': m1, 'cutoff_Hz': cutoff, 'window': label,
                        'start_s': a*dt, 'end_s': b*dt, 'W_highpass': w, 'D_highpass_reference': dac,
                        'q_highpass_RMS': np.std(qhigh[sl]), 'power_highpass': w/((b-a)*dt)})
            # Equal 4999-interval records, identical mean removal and Hann FFT.
            for label, a, b in [('early_equal', 10000, 14999), ('late_equal', 15000, 19999)]:
                _, wm, wf = work_decomposition(force[a:b], dq[a:b])
                equal.append({'side': side, 'mode': m1, 'window': label, 'count': b-a,
                    'duration_s': (b-a)*dt, 'q_AC_RMS': np.std(q[a:b]),
                    'q_linear_detrended_RMS': np.std(signal.detrend(q[a:b])),
                    'qdot_AC_RMS': np.std(v[a:b]), 'Wfluctuating_force': float(wf),
                    'power_fluctuating': float(wf/((b-a)*dt)), 'Wmean_force': float(wm)})
            if m1 not in selected[side]:
                continue
            # Reject unstable peaks; this is not a forced cycle assignment.
            stats = phase_diagnostic(q[10000:], v[10000:], force[10000:], dt)
            frequency = stats['q_peak_Hz']
            if stats['q_spectral_share'] < .65 or frequency < 40:
                cycles.append({'side': side, 'mode': m1, 'accepted': False, 'reason': 'q not sufficiently narrowband'})
                continue
            detrended = signal.detrend(q[10000:])
            peak, _ = signal.find_peaks(detrended, distance=max(1, int(.65/frequency/dt)),
                                        prominence=.5*np.std(detrended))
            periods = np.diff(peak)*dt
            if len(periods) < 3 or np.std(periods)/np.mean(periods) > .15:
                cycles.append({'side': side, 'mode': m1, 'accepted': False, 'reason': 'fewer than 3 stable peak intervals or period CV >15%'})
                continue
            for a0, b0 in zip(peak[:-1], peak[1:]):
                a, b = a0+10000, b0+10000
                cycles.append({'side': side, 'mode': m1, 'accepted': True, 'reason': 'stable q peaks after linear detrend',
                    'start_s': a*dt, 'end_s': b*dt, 'period_s': (b-a)*dt,
                    'amplitude_half_range': .5*np.ptp(detrended[a0:b0+1]),
                    'Wfluid': float(np.sum(force[a:b]*dq[a:b])),
                    'Wcontact_reference': float(np.sum(raw[a:b, s, m, 3]*dq[a:b])),
                    'D_total': float(np.sum(2*zeta[s]*omega[s, m]*v[a:b]**2*dt))})
    return pd.DataFrame(rows), pd.DataFrame(cycles), pd.DataFrame(equal)


def spatial_analysis(raw, flow, geom, static, selected, validation):
    points, modes, surfaces, nxs, sign = geom
    ni, nj = min(nxs), min(s.shape[1] for s in surfaces)
    grid = surfaces[0][:ni, :nj]
    p = points[0][grid]
    monitor = np.unravel_index(np.argmin((p[..., 0]-9.2)**2+(p[..., 2]-8.5)**2), grid.shape)
    # Fixed material points in the reconstructed mean configuration.
    qmean = raw[10000:, 0, :, 4].mean(axis=0)
    mean_coords = displace(points[0], modes[0], qmean)[grid]
    minx = flow.xMin_mm.iloc[10000:19999].mean()
    targets = {'monitor_fixed': monitor,
               'minimum_near_AP_center_fixed': np.unravel_index(np.argmin((mean_coords[..., 0]-minx)**2+(mean_coords[..., 2]-8.5)**2), grid.shape),
               'AP_center_material_fixed': (monitor[0], int(np.argmin(np.abs(p[monitor[0], :, 2]-8.5))))}
    contributions, fields = [], []
    for s, side in enumerate(('L', 'R')):
        sta = static[(static.side == side) & (static.window == 'evaluation')].sort_values('mode')
        qbar = sta.mean_q.to_numpy()
        field = np.einsum('m,mnc->nc', qbar, modes[s], optimize=True)*1e3
        ffield = np.einsum('m,mnc->nc', sta.qstatic_fluid, modes[s], optimize=True)*1e3
        cfield = np.einsum('m,mnc->nc', sta.qstatic_contact, modes[s], optimize=True)*1e3
        for i in range(ni):
            for j in range(nj):
                node = surfaces[s][i, j]
                fields.append({'side': side, 'i': i, 'j': j, 'node': int(node),
                    'x_initial_mm': points[s][node, 0], 'z_initial_mm': points[s][node, 2],
                    'mean_ux_mm': field[node, 0], 'mean_uy_mm': field[node, 1],
                    'mean_opening_y_mm': (-sign if s == 0 else sign)*field[node, 1],
                    'static_fluid_uy_mm': ffield[node, 1], 'static_contact_uy_mm': cfield[node, 1]})
        for name, (i, j) in targets.items():
            node = surfaces[s][i, j]
            for m in range(raw.shape[2]):
                contributions.append({'side': side, 'point': name, 'node': int(node), 'i': i, 'j': j,
                    'mode': m+1, 'mean_ux_mm': modes[s][m, node, 0]*qbar[m]*1e3,
                    'mean_uy_mm': modes[s][m, node, 1]*qbar[m]*1e3,
                    'opening_uy_mm': (-sign if s == 0 else sign)*modes[s][m, node, 1]*qbar[m]*1e3,
                    'opening_shape_per_q': (-sign if s == 0 else sign)*modes[s][m, node, 1]*1e3,
                    'static_fluid_uy_mm': modes[s][m, node, 1]*sta.qstatic_fluid.iloc[m]*1e3,
                    'static_contact_uy_mm': modes[s][m, node, 1]*sta.qstatic_contact.iloc[m]*1e3})
        node = surfaces[s][monitor]
        observed = flow['uyL_mm' if s == 0 else 'uyR_mm'].to_numpy()[:len(raw)]
        reconstructed = raw[:, s, :, 4] @ modes[s][:, node, 1]*1e3
        validation[f'{side}_fixed_monitor_replay'] = close_gate(reconstructed, observed, atol=1e-10, rtol=1e-9, label='monitor replay')
    # Verify and identify moving nearest-separation probes at selected frames.
    moving = []
    for n in validation['AP_selected_frames']:
        coords = displace(points[0], modes[0], raw[n, 0, :, 4])
        pg = coords[grid]
        i, j = np.unravel_index(np.argmin((pg[..., 0]-flow.xSep_mm.iloc[n])**2+(pg[..., 2]-8.5)**2), grid.shape)
        node = grid[i, j]
        close_gate(np.array([coords[node, 0], coords[node, 1]-points[0][node, 1]]),
                   np.array([flow.sep_node_x_mm.iloc[n], flow.sep_node_uyL_mm.iloc[n]]),
                   atol=1e-10, rtol=1e-9, label='moving separation probe')
        moving.append({'step': n, 'node': int(node), 'i': int(i), 'j': int(j)})
    validation['moving_sep_probes'] = moving
    validation['fixed_material_targets'] = {k: [int(v[0]), int(v[1])] for k, v in targets.items()}
    return pd.DataFrame(contributions), pd.DataFrame(fields)


def make_figures(out, coverage, frames, static, dynamic, equal, flow, raw, spatial, selected):
    os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='kanaya-three-checks-mpl-'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def save(fig, name):
        fig.tight_layout(); fig.savefig(out/'figures'/name, dpi=160); plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for step, g in coverage.groupby('step'):
        if step not in (0, 10000, 12500, 15000, 17500, 19500): continue
        label = f'{step*1e-5:.3f}s'
        for ax, y in zip(axs.flat, ('A_reconstructed', 'common_span_mm', 'mean_gap_common_mm', 'open_span_mm')):
            ax.plot(g.x, g[y], '.-', ms=2, label=label); ax.set(xlabel='x [mm]', ylabel=y)
            ax.grid(alpha=.25)
    axs[0, 0].set_yscale('symlog', linthresh=1.)
    axs[0, 0].set_ylim(0, 400); axs[0, 0].legend(fontsize=8)
    save(fig, '01_ap_coverage.png')
    fig, axs = plt.subplots(2, 2, figsize=(12, 7))
    for col, step in enumerate((10000, 17500)):
        x, sect, _ = frames[step]
        k = int(np.nanargmin([r['A'] for r in sect[1:-1]]))+1
        for p in sect[k]['profile']:
            lo, hi, yl0, yl1, yr0, yr1, h0, h1 = p
            axs[0, col].plot([lo, hi], [yl0, yl1], color='C0')
            axs[0, col].plot([lo, hi], [yr0, yr1], color='C1')
            axs[1, col].plot([lo, hi], [h0, h1], color='C2')
        axs[0, col].set(title=f't={step*1e-5:.3f}s, x={x[k]:.3f}mm, plane {k}', ylabel='yL / yR [mm]')
        axs[1, col].axhline(0, color='k', lw=.5)
        axs[1, col].set(xlabel='physical z [mm]', ylabel='signed local gap [mm]')
    save(fig, '02_ap_gap_profile.png')
    modes_to_plot = selected['L']
    st = static[(static.side == 'L') & (static.window == 'evaluation')].set_index('mode').loc[modes_to_plot]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    x = np.arange(len(st))
    axs[0, 0].bar(x-.2, st.mean_q, .4, label='measured mean q')
    axs[0, 0].bar(x+.2, st.qstatic, .4, label='mean Q / omega^2')
    axs[0, 0].set_xticks(x, st.index); axs[0, 0].legend(fontsize=8)
    axs[0, 0].set(xlabel='mode', ylabel='normalized q')
    axs[0, 1].bar(x, st.static_residual_normalized, label='static')
    axs[0, 1].plot(x, st.motion_corrected_normalized, 'ko', ms=3, label='+ mean inertia/damping')
    axs[0, 1].set_xticks(x, st.index); axs[0, 1].legend(fontsize=8)
    axs[0, 1].set(xlabel='mode', ylabel='normalized force residual')
    sp = spatial[spatial.side == 'L']
    for ax, key in zip(axs[1], ('mean_opening_y_mm', 'mean_ux_mm')):
        im = ax.scatter(sp.x_initial_mm, sp.z_initial_mm, c=sp[key], s=5, cmap='coolwarm')
        fig.colorbar(im, ax=ax, label=key); ax.set(xlabel='initial x [mm]', ylabel='initial z [mm]')
    save(fig, '03_static_modal_balance.png')
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for ax, window in zip(axs.flat, ('transient', 'early', 'late', 'evaluation')):
        d = dynamic[(dynamic.side == 'L') & (dynamic.window == window)].set_index('mode').loc[modes_to_plot]
        for offset, key, label in [(-.27, 'Wtotal', 'total'), (-.09, 'Wmean_force', 'mean force'),
                                   (.09, 'Wfluctuating_force', 'fluctuating'), (.27, 'D_ac', 'D AC reference')]:
            ax.bar(x+offset, d[key]*1e6, .18, label=label)
        ax.set_xticks(x, d.index); ax.set(title=window, xlabel='mode', ylabel='normalized work [micro-J equivalent]')
        ax.axhline(0, color='k', lw=.5)
    axs[0, 0].legend(fontsize=7)
    save(fig, '04_dynamic_work.png')
    fig, axs = plt.subplots(3, 1, figsize=(12, 10))
    ev = flow.iloc[10000:19999]
    axs[0].plot(ev.time_s, ev.Amin_mm2, label='Amin'); axs[0].axhline(.01, color='r', label='closure threshold')
    axs[0].set(ylabel='area [mm^2]'); axs[0].legend()
    for m in (1, 4, 5, 7, 35):
        a = equal[(equal.side == 'L') & (equal['mode'] == m)]
        axs[1].plot([.125, .175], a.q_AC_RMS, 'o-', label=f'mode {m}')
        axs[2].plot([.125, .175], a.power_fluctuating*1e6, 'o-', label=f'mode {m}')
    axs[1].set(ylabel='q AC RMS [normalized]'); axs[1].legend(ncol=5)
    axs[2].set(xlabel='time [s]', ylabel='fluctuating work rate [micro-W equivalent]')
    axs[2].axhline(0, color='k', lw=.5)
    save(fig, '05_amplitude_and_gap.png')


def markdown_table(frame, columns):
    # Avoid an optional tabulate dependency.
    result = ['|'+'|'.join(columns)+'|', '|'+'|'.join(['---']*len(columns))+'|']
    for _, row in frame[columns].iterrows():
        values = [f'{v:.6g}' if isinstance(v, (float, np.floating)) else str(v) for v in row]
        result.append('|'+'|'.join(values)+'|')
    return '\n'.join(result)


def build_report(out, coverage, static, dynamic, trends, cycles, equal, contributions,
                 spatial, flow, meta, validation):
    ev = flow.iloc[10000:19999]
    ap = pd.DataFrame(validation['AP_minimum_comparisons'])
    ap_ev = ap[ap.step >= 10000]
    ap_ev = ap_ev.assign(time_s=ap_ev.step*meta['dt_s'],
        mean_gap_common_mm=ap_ev.minimum_area_mm2/ap_ev.common_span_mm,
        width_loss_percent=100*(1-ap_ev.common_span_fraction_of_frame_max),
        width_log_fraction_percent=100*ap_ev.reference_to_min_log_width_fraction)
    ap_ev.loc[ap_ev.width_loss_percent.abs() < 1e-10, 'width_loss_percent'] = 0.
    selection = pd.DataFrame(validation['AP_width_selection_sensitivity'])
    selection_eval = selection[selection.step >= 10000]
    selection_changed = selection_eval[selection_eval.different_plane]
    sta = static[(static.side == 'L') & (static.window == 'evaluation')]
    dyn = dynamic[(dynamic.side == 'L') & (dynamic.window == 'evaluation')]
    monitor = contributions[(contributions.side == 'L') & (contributions.point == 'monitor_fixed')]
    measured = float(monitor.mean_uy_mm.sum())
    fluid_prediction = float(monitor.static_fluid_uy_mm.sum())
    contact_prediction = float(monitor.static_contact_uy_mm.sum())
    predicted = fluid_prediction+contact_prediction
    qs_error = abs(predicted-measured)/abs(measured)
    work_keys = ['Wtotal', 'Wmean_force', 'Wfluctuating_force', 'Wpos', 'Wneg', 'D_total', 'D_ac', 'Wcontact_reference']
    totals = dynamic[(dynamic.side == 'L') & dynamic.window.isin(['transient', 'early', 'late', 'evaluation'])].groupby('window')[work_keys].sum()
    eq = equal[equal.side == 'L']
    eq_a = eq[eq.window == 'early_equal'].set_index('mode')
    eq_b = eq[eq.window == 'late_equal'].set_index('mode')
    modes = validation['selected_modes']['L']
    amplitude = pd.DataFrame({'mode': modes,
        'q_RMS_early': eq_a.loc[modes].q_AC_RMS,
        'q_RMS_late': eq_b.loc[modes].q_AC_RMS,
        'RMS_late_over_early': eq_b.loc[modes].q_AC_RMS/eq_a.loc[modes].q_AC_RMS,
        'detrended_RMS_ratio': eq_b.loc[modes].q_linear_detrended_RMS/eq_a.loc[modes].q_linear_detrended_RMS,
        'power_early_micro': eq_a.loc[modes].power_fluctuating*1e6,
        'power_late_micro': eq_b.loc[modes].power_fluctuating*1e6})
    total_rate_a, total_rate_b = float(eq_a.power_fluctuating.sum()), float(eq_b.power_fluctuating.sum())
    aggregate_rms_a = float(np.sqrt(np.sum(eq_a.q_AC_RMS**2)))
    aggregate_rms_b = float(np.sqrt(np.sum(eq_b.q_AC_RMS**2)))
    area_early = flow.iloc[10000:14999].Amin_mm2.to_numpy()
    area_late = flow.iloc[15000:19999].Amin_mm2.to_numpy()
    amin = float(ev.Amin_mm2.min()); close_threshold = .01
    st_table = sta[sta['mode'].isin(modes)].copy()
    st_table['static_error_percent'] = st_table.static_residual_normalized*100
    st_table['motion_corrected_percent'] = st_table.motion_corrected_normalized*100
    summary_spatial = contributions.groupby(['side', 'point'])[['mean_ux_mm', 'mean_uy_mm', 'static_fluid_uy_mm', 'static_contact_uy_mm']].sum().reset_index()
    spatial_table = monitor[monitor['mode'].isin([1, 4, 5, 7, 35])]
    top_open = monitor.assign(abs_y=monitor.mean_uy_mm.abs()).nlargest(8, 'abs_y')
    trend_totals = trends[trends.side == 'L'].groupby(['window', 'cutoff_Hz'])[['W_highpass', 'D_highpass_reference', 'power_highpass']].sum().reset_index()
    trend_top = trends[(trends.side == 'L') & (trends.window == 'common_trim_25ms')].copy()
    ranks = {int(fc): g.nlargest(5, 'W_highpass')['mode'].to_list() for fc, g in trend_top.groupby('cutoff_Hz')}
    unstable = []
    for mode, group in trend_top.groupby('mode'):
        # Suppress effectively zero (<1e-12 normalized J) numerical signs.
        if group.W_highpass.max() > 1e-12 and group.W_highpass.min() < -1e-12:
            unstable.append(int(mode))
    full_trim = trend_totals[trend_totals.window == 'evaluation_full_edges'].power_highpass.to_numpy()
    cut_trim = trend_totals[trend_totals.window == 'common_trim_25ms'].power_highpass.to_numpy()
    cycle_good = cycles[(cycles.side == 'L') & cycles.accepted]
    cycle_summary = cycle_good.groupby('mode').agg(cycles=('period_s', 'size'),
        period_min_s=('period_s', 'min'), period_max_s=('period_s', 'max'),
        amplitude_min=('amplitude_half_range', 'min'), amplitude_max=('amplitude_half_range', 'max'),
        Wfluid_sum=('Wfluid', 'sum'), D_total_sum=('D_total', 'sum')).reset_index()
    phase_table = dyn[dyn['mode'].isin([1, 4, 5, 7, 35])][['mode', 'q_peak_Hz', 'qdot_peak_Hz', 'Q_peak_Hz', 'phase_accepted', 'phase_deg']]
    relation = dyn[dyn['mode'].isin([1, 4, 5, 7, 35])].merge(
        monitor[['mode', 'opening_shape_per_q']], on='mode')
    force_rows = pd.DataFrame(validation['decimated_force_regions_mean'])
    fpos = float(force_rows[force_rows.pressure_sign == 1].opening_force_N.sum())
    fneg = float(force_rows[force_rows.pressure_sign == -1].opening_force_N.sum())
    max_area_error = max(r['area']['max_abs'] for r in validation['AP_reconstruction'])
    max_amin_error = max(r['Amin_flow']['max_abs'] for r in validation['AP_reconstruction'])
    max_x_error = max(r['x']['max_abs'] for r in validation['AP_reconstruction'])
    max_g3_error = max(r['G3']['max_abs'] for r in validation['AP_reconstruction'])
    source_count = len(validation['stored_summary_checks'])
    before, after = validation['git_before'], validation['git_after']
    same_left_right = dynamic.pivot(index=['window', 'mode'], columns='side', values='Wtotal')
    lr_error = float(np.max(np.abs(same_left_right.L-same_left_right.R)))
    text = f'''# Kanaya型周期的閉鎖：既存1ランの3項目後処理検証

1. **大幅なAP幅欠落は不支持だが、最狭位置選択への局所的な幅の影響は支持。** 評価窓{len(ap_ev)}フレームで共通AP幅は{ap_ev.common_span_mm.min():.6f}〜{ap_ev.common_span_mm.max():.6f} mm（最大幅の{100*ap_ev.common_span_fraction_of_frame_max.min():.3f}〜100%）。ただし面積最小とA/S_common最小の断面が{len(selection_changed)}/{len(selection_eval)}フレームで異なり、位置差は最大{selection_changed.x_shift_mm.abs().max():.6f} mm。全{len(ap)}フレーム×50断面を丸め内で再現し、高精度Aminとの差は最大{max_amin_error:.3g} mm²。
2. **平均開大の準静的説明は、現行モデル内では支持。** 固定モニタ点の片側平均uyは{measured:.6f} mm、平均荷重/ω²の予測は{predicted:.6f} mm（差{100*qs_error:.3f}%）。モード1の開大と4・5の閉鎖方向寄与が大きく相殺する。物理的な材料剛性・荷重面の妥当性まで証明する結果ではない。
3. **変動力からの仕事供給は正。ただし全体の単調な自励成長は未判定。** 等長49.99 ms窓の全100モード合計変動仕事率は{total_rate_a*1e6:.3f}→{total_rate_b*1e6:.3f} µW相当。AminのAC RMSは{np.std(area_early):.4f}→{np.std(area_late):.4f} mm²だが、固定点uyのAC RMSは{np.std(flow.uyL_mm.iloc[10000:14999]):.5f}→{np.std(flow.uyL_mm.iloc[15000:19999]):.5f} mm。評価窓の最小面積{amin:.6f} mm²は閉鎖閾値{close_threshold} mm²をなお{amin-close_threshold:.6f} mm²上回る。

## 意思決定

| 診断結果 | 次の1手（まだ実行しない） |
|---|---|
| 面積最小が共通AP幅縮小に支配される | 全般的な大幅欠落は不支持だが、拮抗する最狭候補では幅の影響が支配的になる時刻あり。下記の剥離開始候補の選択指標だけを変える1条件比較を選ぶ。 |
| 幅の問題は小さく、平均荷重で平均開大を説明できる | 準静的説明は支持。ただし位置選択への幅の影響が確認されたため、今回は材料・初期gap・荷重面を変える比較を重ねない。物理的に「過大」かは実験の変位・材料値との照合待ち。 |
| 平均開大は妥当だが、変動流体仕事が不足/負 | 全モード一括の供給不足は不支持。個別モードの不足・流れ方向位相・圧力回復の因果は未判定であり、今回は別条件を追加しない。 |
| 仕事は正、振幅が成長中 | 7・35等は増幅、1・4は減幅し、モニタ点はほぼ横ばい。最狭面積の平均も上昇する。時間延長の必要性は残るが、単調な閉鎖への接近は確認できないため、今回は延長を次の1条件には選ばない。 |
| データ/再構成不一致 | 必須ゲートは全て合格。以前の単一スペクトル位相値は今回の帯域一致判定を満たさず、原因推定には採用しない。 |

**次の比較計算は1条件だけ：剥離開始候補の選択指標をAからA/S_commonへ変更する幾何選択感度比較。** 将来の別変更では`channelSections.cpp`で既に計算されている共通AP幅spanを読取り可能に保持し、`separationSelection.cpp`の下流検索の開始候補だけを、有効な内部断面のargmin(A/S_common)へ替える。全体の物理的Aminとその閉鎖判定は現行のargmin(A)のまま、各断面の実A/G3、流体係数の式、50断面、既存の面積比条件はそのままとする。span=0は未定義として候補から除く。今回のareaRatio=1では検索開始候補がそのまま剥離点となるため、変更因子を位置選択に限定できる。

凍結した既存形状では{len(selection_changed)}/{len(selection_eval)}フレームで剥離候補が約{selection_changed.x_shift_mm.abs().min():.3f}〜{selection_changed.x_shift_mm.abs().max():.3f} mm上流へ移る。これは現行ブレンド幅0.5 mmの約半分であり、小さい面積差でも負圧割当の位置に識別可能な差を作り得る。同じ形状のAmin/Rkは維持され、剥離点まで積分するLg/Rvは候補が上流へ移るフレームでは小さくなる方向。ブレンド帯も上流へ動くため、正負の力・平均開大・変動仕事の変化で位置選択の影響を調べる。再連成後の力の増減や閉鎖成立は保証しない。

これは**位置選択に対する人工的な識別比較**であり、A/Sを実面積に置換したり、欠けたAP域を補間したりする物理修正ではない。異なるAP分布を平均gap一つへ縮約する副作用があり、これを正しい剥離則と断言しない。Ps=1000 Pa、初期gap=0、soft材料・100モード・ζ=0.06、接触、荷重面、回復幅0.5 mm、dt、20,000ステップ等は全て固定する。モデル本体への変更が必要な別タスクであり、**このタスクでは変更・実行していない。**

## 入力・数値整合性

- 作業場所：`{ROOT}`。実作業HEAD：`{before['head'].strip()}`、ブランチ：`{before['branch'].strip()}`。本番時HEADとの違いは`manifest.json`のハッシュで区別した。
- 連成計算の起動0回、ソルバの再コンパイル0回、物理/数値モデル変更0件、既存出力の上書き0件。既存`analyze.py`も実行・編集していない。
- `complete=true`、20,000ステップに対して積分19,999区間、除外1区間、左右100モード、dt={meta['dt_s']:.8g} s、実ランプ0.05 sを検証。仕事区間は[0,0.19999] s、最後の荷重時刻は0.19998 s。除外されたn=19999を仕事・形状検証へ戻していない。
- `modal_native.bin`は{validation['native_bytes']:,} bytes。読取専用little-endian memmap。NaN/infなし、`Qpos+Qneg=Qfluid`最大絶対差{validation['Q_split']['max_abs']:.3g}、q連続性最大差{validation['q_continuity']['max_abs']:.3g}。既存選択モードCSVとも照合した。
- `consistency.csv`の39,998有効行全てvalid=1、最大規格化残差{validation['consistency_scaled_max']:.6g}<1。除外行は左右のn=19999だけ。
- 左右×全100モード×保存4窓について、Wpos/Wneg/Wfluid/Wcontact/Dの{source_count}配列照合が合格（最大絶対差{max(r['max_abs'] for r in validation['stored_summary_checks']):.3g}）。50ms分割の生仕事・減衰和もevaluationと一致。窓ごとに平均を再定義するWmean/Wfluctuatingは分割加法性を要求しない。
- DATはarea/section_x/gap_cubed/pressureがstep+50列、separationが9列、4,000行、同じn=5k。全DATのAmin/xMin/xSepを高精度flowと照合した。面積最小選択は現行通り内部断面1〜48で行い、入口0・出口49を最小選択へ混入させない。
- 元VTUから(raw Z,X,Y)へ軸変換し、全節点等分質量で各モードを正規化、NASトポロジーから64×69表面を独立生成。`surfp_output.csv`の全IDとも一致。nxsup=42、右y鏡映とuy符号反転後、初期gap剛体補正を同順で実施（今回の補正は片側{validation['geometry']['initial_gap_translation_per_side_mm']:.6g} mm）。x,z変形も含める。
- VTU3時刻はn=15000+20kの**更新後q_(n+1)**から再現、座標誤差最大{max(r['max_abs'] for r in validation['VTU_time_offset']):.3g} mm。誤ってq_nを使う対照では最大{max(r['wrong_q_n_max_abs_mm'] for r in validation['VTU_time_offset']):.6g} mmずれ、時刻ゲートが識別できた。
- 合成テスト{validation['synthetic_test_count']}件合格。一定幅・欠落AP行・真のくびれ・部分接触・交点なし/NaN・面積/平均gapの順位逆転・定数力・正逆位相・平均開大のみ・時刻ずれ・静的/非定常補正を検証。実ソルバを起動するテストはない。

## A：共通AP幅と真のくびれ

再構成した{len(ap)}フレームは`validation.json/AP_selected_frames`に列挙した。初期、0.05、0.10、0.125、0.15、0.175、0.195 s、Amin/xSep極値近傍を含む。厳密なDAT照合のため極値stepを最寄りの5刻みへ合わせた。全時刻のAmin極値自体は高精度flowから集計した。AP幅はこのサンプルに限定した結論で、全19,999時刻の幅を計算したわけではない。

再構成誤差はA最大{max_area_error:.6g} mm²、x最大{max_x_error:.6g} mm、G3最大{max_g3_error:.6g} mm⁴。通常DATがC++ defaultfloatの有効6桁なので、各値について`0.50001×10^(floor(log10(|value|))-5)+2e-10`を丸め許容値にした。全セルが丸め内、NaN位置も一致。高精度flowのAmin差{max_amin_error:.3g} mm²は丸めより十分小さい。

{markdown_table(ap_ev, ['time_s', 'minimum_plane', 'minimum_area_mm2', 'common_span_mm', 'open_span_mm', 'mean_gap_common_mm', 'width_loss_percent'])}

参照幅は同一フレーム内の最大共通AP幅（約17 mm）。最狭断面でも共通幅の損失は最大{ap_ev.width_loss_percent.max():.3f}%、A/S_commonは{ap_ev.mean_gap_common_mm.min():.6f}〜{ap_ev.mean_gap_common_mm.max():.6f} mmで有限のgapを持つ。最狭点のS_openは{ap_ev.open_span_mm.min():.6f}〜{ap_ev.open_span_mm.max():.6f} mmで、部分接触があってもS_commonが消えるわけではない。

さらに各最狭点の2断面以上上流にある、最大共通幅の99.5%以上の最近接正常断面を参照した。A比を`(S/Sref)×[(A/S)/(Aref/Sref)]`へ分解し、幅比{ap_ev.span_ratio.min():.6f}〜{ap_ev.span_ratio.max():.6f}、平均gap比{ap_ev.gap_ratio.min():.6f}〜{ap_ev.gap_ratio.max():.6f}、面積比{ap_ev.area_ratio.min():.6f}〜{ap_ev.area_ratio.max():.6f}。同参照に対する対数面積減少の幅寄与は{ap_ev.width_log_fraction_percent.min():.3f}〜{ap_ev.width_log_fraction_percent.max():.3f}%（参照依存の幾何分解であり、因果寄与率ではない）。全比較値はvalidation.jsonへ保存した。

**判定：共通AP幅の大幅欠落による面積過小評価は不支持。ただし最狭/剥離位置の選択への局所的な影響は支持、非閉鎖の原因としては未判定。** 例えば上流参照断面との差を分解すると、幅の減少が面積差の最大{ap_ev.width_log_fraction_percent.max():.3f}%を占める時刻がある。絶対幅の損失が小さくても、平坦な面積極小の順位には大きく効く。面積最小と共通幅で割った平均gap最小を別々に検索した結果、以下の時刻では選択位置が変わった（他の{len(selection_eval)-len(selection_changed)}フレームは同じ）。

{markdown_table(selection_changed.assign(time_s=selection_changed.step*meta['dt_s']), ['time_s', 'area_min_plane', 'gap_min_plane', 'x_shift_mm', 'area_min_A_mm2', 'gap_min_A_mm2', 'area_min_mean_gap_mm', 'minimum_mean_gap_mm'])}

この検索は全有効内部断面に対して行い、95%などの幅閾値は追加していない。幅の小さい端部を都合よく削って順位を変えた結果ではない。全フレームの検索結果はvalidation.json/AP_width_selection_sensitivityに保存。A/S_commonは接触部分も含む正gapの幅平均、A/S_openは開口部分だけの平均であり、同一量ではない。A/Sの順位比較は実Aや現行ソルバの位置選択を書き換えない。

入口/出口や最狭点の隣接±1〜3断面を含む全50断面は`ap_coverage.csv`と下図に保存。交線がない場合はarea=NaN、span=0とし、閉鎖に数えない。区間は左右のz端点で分割し、穴を補間せず重複を二重計数しない。

50断面のA/G3と剥離選択はLg/Rv/Rkへ波及し得る。最狭点の幅損失が小さくても、他断面を通る積分係数の正しさまで証明しない。今回それらの係数や流体式は変更・再計算していない。

![AP coverage](figures/01_ap_coverage.png)
![Local gap profiles at two times](figures/02_ap_gap_profile.png)

## B：平均荷重と準静的変形

`rungeStep.cpp`のNewmark式と`simulation.cpp`の最終荷重を確認した。モデル式はqddot+2ζω qdot+ω²q=Qfluid+Qcontact、ω=2πf。周波数ファイルの減衰0ではなく、呼出しで指定するζ=0.06を使用する。荷重Q_nは**更新後のq_(n+1),v_(n+1),a_(n+1)**の式に入る。保存速度からγ=0.55の関係で加速度を復元した離散式残差は最大{validation['recorded_newmark_equation']['max_abs']:.3g}（n=0〜19997。次の保存速度がない最後の有効区間はこの追加ゲートのみ対象外）。時間発展の再実行はしていない。

`static_balance.csv`には左右100モード、全6窓（allと指定5窓）を保存。mean(qddot)は窓内速度に2次精度の中央差分、両端は2次片側差分を適用して平均した。実際のNewmark右端荷重と左端平均のO(dt)ずれや微分ノイズが残るため、運動補正後の残差をゼロに強制していない。規模が1e-8未満の一般化力について相対残差をNaNにし、絶対残差を提示する。

{markdown_table(st_table, ['mode', 'mean_q', 'qstatic', 'static_error_percent', 'motion_corrected_percent'])}

qはmmではなく質量正規化一般化座標。モード1の準静的差は約{100*abs(sta[sta['mode'] == 1].static_residual_normalized.iloc[0]):.3f}%、4・5も小さい。各窓の残留運動はCSVを参照。100モードから平均|Q|・平均|q|・変動仕事上位を選び、重点表示モードは{modes}とした。

固定材料点での平均変位と平均流体/接触荷重の予測（mm）：

{markdown_table(summary_spatial[summary_spatial.side == 'L'], ['point', 'mean_ux_mm', 'mean_uy_mm', 'static_fluid_uy_mm', 'static_contact_uy_mm'])}

モニタは初期座標のx=9.2,z=8.5最近傍から決めた固定点（node={int(monitor.node.iloc[0])}, i={int(monitor.i.iloc[0])}, j={int(monitor.j.iloc[0])}）。今回AP中央指定点は同じ点になる。最狭近傍点も平均形状から一度選んだ固定材料点で比較した。一方、診断のnear-separationは毎時刻再選択される点であり、選択フレームでnodeが{len(set(r['node'] for r in validation['moving_sep_probes']))}種類に切り替わることを再現した。同じ材料点の波形とは解釈しない。

モニタ点での主要な加法寄与（opening_uy>0が開大、左右片側）：

{markdown_table(spatial_table, ['mode', 'mean_ux_mm', 'opening_uy_mm', 'static_fluid_uy_mm', 'static_contact_uy_mm'])}

全100モードの|y寄与|上位：

{markdown_table(top_open, ['mode', 'opening_uy_mm', 'mean_ux_mm'])}

モード1の開大を4・5が大きく打ち消すため、個別寄与を正味変位に対する単純な割合として足してはいけない。正圧の広い上流荷重はforce_regions.csvから再集計し、評価窓の片側開大力は正圧{fpos*1e3:.6f} mN、負圧{fneg*1e3:.6f} mN、正味{(fpos+fneg)*1e3:.6f} mN。20ステップ間引き集計と元の毎ステップ統計は0.1%内で一致。接触は予測反復で採用された荷重で、流体左端荷重と同時刻と断言しない。

{markdown_table(force_rows, ['pressure_sign', 'zone', 'opening_force_N', 'applied_area_mm2'])}

**判定：モデル内の平均変形は平均流体荷重と主な低次モードで説明でき、平均接触荷重の変位寄与は小さい。** ただしモニタの平均x変位は{monitor.mean_ux_mm.sum():.6f} mmあり、固定材料点が流れ方向へ大きく移る。線形モードの大変形適用、材料値、荷重対象面の物理的妥当性は未検証。全100モードの等分質量Gramの非対角最大は{validation['geometry']['full_gram_max_off_diagonal']:.6f}（以前の選択5モードの値より大きい）。独立した物理モードエネルギー保存は仮定しない。

![Mean displacement and modal balance](figures/03_static_modal_balance.png)

## C：変動仕事・振幅・時間発展

各窓・各モードでWtotal=ΣQ_nΔq、Wmean=mean(Q)ΣΔq、Wfluct=Σ(Q−mean(Q))Δqを**別々に計算して恒等性を検証**した。接触仕事は採用予測荷重の参考離散仕事であり散逸ではない。D_total=Σ2ζω v_n²dt、D_ac=Σ2ζω(v_n−mean(v))²dtを併記した。Wfluct−D_acをエネルギー保存残差とは呼ばない。

片側・全100モードの集計（µJ相当、正規化モデル量）：

{markdown_table((totals*1e6).reset_index(), ['window', *work_keys])}

全評価窓ではWtotal={totals.loc['evaluation', 'Wtotal']*1e6:.6f} µJ相当、Wmean={totals.loc['evaluation', 'Wmean_force']*1e6:.6f}、Wfluct={totals.loc['evaluation', 'Wfluctuating_force']*1e6:.6f}。以前の集計値を入力にせずbinから再集計して保存summaryと一致した。負の全仕事だけから「流体供給なし」とは言えない。

前後半比較は同じ4999区間・49.99 ms、同じ平均除去・同じHannスペクトル解析。前半[0.10,0.14999]、後半[0.15,0.19999] sで、前半の最後1区間だけ等長比較から外す。通常のearly窓[0.10,0.15]は5000区間として別途維持し、early+late=evaluationを検証した。

{markdown_table(amplitude, ['mode', 'RMS_late_over_early', 'detrended_RMS_ratio', 'power_early_micro', 'power_late_micro'])}

全モードqのAC RMS二乗和平方根は{aggregate_rms_a:.6g}→{aggregate_rms_b:.6g}（比{aggregate_rms_b/aggregate_rms_a:.3f}）。非直交なので物理的総振幅とは呼ばない。1・4が減少し、5・7・35等が増加する。線形トレンド除去の感度も上表に併記。50ms窓は45 Hzモードで約2周期しかなく、短い窓と拍・過渡の影響が大きい。

補助解析は全時系列へ4次Butterworth高域通過を前後適用（zero phase、20/30/40 Hz）。低域トレンドを差し引いたq′の端点差を用いてΣQ′Δq′を計算し、速度も同じフィルタで処理した。20,30,40 Hzは片方向フィルタのカットオフで、前後適用の振幅応答は二乗される。端点を含む0.10〜0.19999 sと、終端25 msを捨てた全カットオフ共通の0.10〜0.17499 sを比較。始端は少なくとも0.10 s離れている。25 ms除去は端点影響の完全消失を保証しないため、以下は感度診断である。

{markdown_table(trend_totals.assign(W_highpass_micro=trend_totals.W_highpass*1e6,
    D_highpass_micro=trend_totals.D_highpass_reference*1e6,
    power_highpass_micro=trend_totals.power_highpass*1e6), ['window', 'cutoff_Hz', 'W_highpass_micro', 'D_highpass_micro', 'power_highpass_micro'])}

共通トリム窓の全モード仕事は全カットオフで正、上位5順位は{ranks}。有意仕事の符号がカットオフで反転するモードは{unstable if unstable else 'なし（|W|>1e-12を有意とした）'}。端点込み/トリム後の仕事率の相対差は最大{100*np.max(np.abs(full_trim-cut_trim)/np.maximum(abs(cut_trim), 1e-30)):.3f}%だが、窓自体も異なるため純粋なフィルタ誤差とはみなさない。正の供給と上位順位の安定性までを支持し、厳密な自励仕事・成長率とは断言しない。

主要スペクトルの帯域一致・集中度・参考Welchコヒーレンスを要求した結果：

{markdown_table(phase_table, ['mode', 'q_peak_Hz', 'qdot_peak_Hz', 'Q_peak_Hz', 'phase_accepted', 'phase_deg'])}

Hann窓は線形トレンド除去後、ピーク±1ビンのパワー比≥0.65、Qと速度の主要ビン差≤1、参考coherence≥0.7の場合だけ位相を出す。短い記録のcoherence自体も因果証拠ではない。主要5モードは条件を満たさず、単一位相を採用しない。これは以前のreportの一つのqピークで計算した位相値と定義が違うためであり、元データや符号の矛盾ではない。

周期別評価はqの線形トレンドを除き、ピーク間隔≥主ピーク周期の65%、prominence≥0.5 RMS、3周期以上かつ周期CV≤15%を満たす場合だけ実施した。採用周期も正負仕事が混在し、振幅は単調でない。

{markdown_table(cycle_summary, ['mode', 'cycles', 'period_min_s', 'period_max_s', 'amplitude_min', 'amplitude_max', 'Wfluid_sum', 'D_total_sum'])}

条件を満たさなかったモードはcycle_work.csvに理由を記録した。ピークで切った周期は、未処理qの状態が厳密に閉じる周期とは限らず、Wmeanの残留も含む。周期仕事の正負だけで成長を断言しない。

同時刻のAmin/Ug/xSepとモードq・qdotの相関、および負圧仕事とモニタでのy形状を照合：

{markdown_table(relation, ['mode', 'Wneg', 'opening_shape_per_q', 'Amin_mm2_q_correlation', 'Amin_mm2_qdot_correlation', 'Ug_m3_s_q_correlation', 'xSep_mm_q_correlation'])}

例えばopening_shape_per_qの符号は正のqがその固定点を開くか閉じるかを示すだけで、モード座標の符号自体は任意。Wneg>0をそのまま「閉鎖運動への供給」と解釈しない。x,z変位による流路位置の変化や他モードとの相殺も含む。相関は同時刻・平均除去の観測量であり因果ではない。共通AP被覆が異なる局所断面の見かけの位相を混ぜないため、16フレームから流れ方向の位相遅れは推定しなかった。

![Work decomposition](figures/04_dynamic_work.png)
![Amplitude and area](figures/05_amplitude_and_gap.png)

## 原因の考察と限界

このランで直接確認できるのは、**平均流体荷重に支えられた有限の開口のまわりで、複数モードに変動仕事が入り、閉鎖に足りる一体的な変形には達していない**こと。評価窓Aminの平均は{ev.Amin_mm2.mean():.6f} mm²、前半/後半平均は{np.mean(area_early):.6f}/{np.mean(area_late):.6f} mm²、最小値は{np.min(area_early):.6f}→{np.min(area_late):.6f} mm²。後半のAC RMS増大は観測されるが、平均は低下せず、最小値も閾値の{amin/close_threshold:.1f}倍。平均から閾値までの距離を後半AC RMSで割ると{(np.mean(area_late)-close_threshold)/np.std(area_late):.2f}倍で、ガウス仮定や線形外挿による閉鎖予測には使わない。

したがって「大幅なAP幅不足だけで小面積になった」「全体の流体仕事が負だから供給がない」という説明は不支持。一方で、幅の数%の減少でも面積最小の順位と現行areaRatio=1の剥離位置を変える局面があり、無視できるとは言えない。平均開大を作る現行の剛性・荷重の組合せとモード相殺は有力なモデル内説明だが、どちらが実物と違うかは1ランから確定できない。まず観測できた位置選択の感度だけを比較する。成長中のモードがあるため、より長い時間で閉鎖に至る可能性も排除できない。

左右のWtotal差は最大{lr_error:.3g}で鏡映対称。報告は主に左側を表示したがCSVは左右全100モードを保持する。q・ω²q・Qの単位と等分節点質量近似、全Gramの非直交性、採用接触荷重の予測段階、有限時間/フィルタ端点、最終区間除外、線形モードの大変形を限界とする。新しい外部文献値や選択5モードのnpzを入力に用いていない。

## 保護監査・成果物

既存runの全出力・元ソース・使用したVTU/NAS/周波数等を`manifest.json/existing_files_before`にパス、解決先、サイズ、mtime、SHA256として記録し、解析後に全ハッシュとmtimeを再照合した。元`/home/kajiyama/code/simulation-main`のGit HEAD/status/diff/indexとinput/output/result/build/tools/diagnosticsの一覧・サイズ・mtime集約も作業前後で一致。元worktreeの既存未コミット変更はそのまま。

実装worktreeの作業前status（追加スクリプトとテストを含む）：

```text
{before['status'].rstrip()}
```

作業後status：

```text
{after['status'].rstrip()}
```

追加コードは`tools/kanaya_diagnostics/three_checks.py`と`tests/test_kanaya_three_checks.py`の2件。commit/pushは行っていない。Releaseビルド要件は計算本体を変更していないため既存の検証済みビルドを保持し、今回は指定書に従い再ビルド・連成計算を実行しなかった。

実行コマンド（既存出力があれば上書きを拒否）：

```bash
OPENBLAS_NUM_THREADS=1 python3 -B -m unittest discover -s tests -p test_kanaya_three_checks.py -v
OPENBLAS_NUM_THREADS=1 python3 -B tools/kanaya_diagnostics/three_checks.py
# 再解析時は --output-name postprocess_3checks_NEW を指定する
```

レポート絶対パス：`{out/'report.md'}`

```text
postprocess_3checks/
├── report.md
├── manifest.json
├── validation.json
├── ap_coverage.csv
├── static_balance.csv
├── dynamic_work.csv
├── trend_sensitivity.csv
├── cycle_work.csv
├── equal_windows.csv
├── spatial_contributions.csv
├── mean_displacement.csv
└── figures/
    ├── 01_ap_coverage.png
    ├── 02_ap_gap_profile.png
    ├── 03_static_modal_balance.png
    ├── 04_dynamic_work.png
    └── 05_amplitude_and_gap.png
```
'''
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output-name', default='postprocess_3checks')
    args = parser.parse_args()
    root = args.root.resolve()
    if Path(args.output_name).name != args.output_name:
        parser.error('output-name must be a single new subdirectory')
    base = root/'output/kanaya_diagnostic'; out = base/args.output_name
    if out.exists():
        raise FileExistsError(f'refusing to overwrite {out}; use a new --output-name')
    started = time.time()
    before_git, before_original = git_state(root), original_snapshot()
    print(f'cwd={Path.cwd()}\nroot={root}\nbranch={before_git["branch"].strip()}\nHEAD={before_git["head"].strip()}\nstatus={before_git["status"]}', flush=True)
    test = subprocess.run([sys.executable, '-B', '-m', 'unittest', 'discover', '-s', 'tests',
                           '-p', 'test_kanaya_three_checks.py', '-v'], cwd=root, capture_output=True, text=True)
    print(test.stdout+test.stderr, flush=True)
    if test.returncode: raise RuntimeError('synthetic tests failed: no physical conclusions')
    existing = [p for folder in ('output', 'src', 'include', 'input') for p in (root/folder).rglob('*') if p.is_file()]
    existing += [root/'CMakeLists.txt', root/'AGENTS.md', root/'tools/kanaya_diagnostics/analyze.py',
                 root/'docs/plans/KANAYA_THREE_CHECKS_POSTPROCESS.md']
    input_manifest = dict(line.split(' = ', 1) for line in (root/'output/txt/manifest.txt').read_text().splitlines() if ' = ' in line)
    existing += [Path(input_manifest[side+'_'+kind]) for side in ('left', 'right')
                 for kind in ('mode_vtu', 'frequency', 'surface_nas')]
    print(f'Hashing {len(existing)} existing inputs/outputs/source files...', flush=True)
    before = inventory(existing)
    raw, meta, flow, summary, manifest, omega, zeta, dat, validation = read_run(root)
    print('Native / summary / DAT / discrete-equation gates passed.', flush=True)
    geom = geometry_setup(root, manifest, raw, meta, validation)
    print('Geometry and VTU time-offset gates passed.', flush=True)
    coverage, frames = ap_analysis(raw, flow, dat, geom, validation)
    print(f'AP reconstruction passed for {len(frames)} frames.', flush=True)
    static, dynamic, selected = modal_analysis(raw, meta, flow, omega, zeta, validation)
    trends, cycles, equal = trend_and_cycles(raw, omega, zeta, meta['dt_s'], selected)
    contributions, spatial = spatial_analysis(raw, flow, geom, static, selected, validation)
    print('Static / dynamic / filter / fixed and moving probe checks passed.', flush=True)
    # No output directory is created until every mandatory numeric gate passes.
    after = inventory(existing)
    assert before == after, 'existing source or output changed during postprocessing'
    assert before_original == original_snapshot(), 'original worktree changed'
    assert before_git == git_state(root), 'existing Git state changed unexpectedly'
    out.mkdir(); (out/'figures').mkdir()
    for name, frame in [('ap_coverage', coverage), ('static_balance', static), ('dynamic_work', dynamic),
                        ('trend_sensitivity', trends), ('cycle_work', cycles), ('equal_windows', equal),
                        ('spatial_contributions', contributions), ('mean_displacement', spatial)]:
        frame.to_csv(out/(name+'.csv'), index=False, mode='x')
    make_figures(out, coverage, frames, static, dynamic, equal, flow, raw, spatial, selected)
    validation.update({'synthetic_tests_passed': True, 'synthetic_test_log': test.stdout+test.stderr,
                       'synthetic_test_count': int(re.search(r'Ran (\d+) tests', test.stdout+test.stderr)[1]),
                       'all_mandatory_gates_passed': True, 'solver_runs': 0, 'core_model_edits': 0,
                       'existing_file_overwrites': 0, 'existing_file_hashes_preserved': True,
                       'original_worktree_preserved': True, 'git_before': before_git, 'git_after': git_state(root),
                       'elapsed_s': time.time()-started})
    json_write(out/'validation.json', validation)
    json_write(out/'manifest.json', {'root': str(root), 'output': str(out), 'started_unix': started,
        'finished_unix': time.time(), 'git_before': before_git, 'git_after': git_state(root),
        'original_worktree_before': before_original, 'original_worktree_after': original_snapshot(),
        'existing_files_before': before, 'existing_files_after_identical': True,
        'script_sha256': sha(Path(__file__)), 'test_sha256': sha(root/'tests/test_kanaya_three_checks.py'),
        'run': {key: meta[key] for key in ('dt_s', 'nmode', 'integrated_intervals', 'excluded_intervals', 'actual_ramp_time_s')},
        'units': {'coordinates': 'mm', 'area': 'mm^2', 'G3': 'mm^4', 'span': 'mm',
                  'modal_q': 'sampled-mass-normalized, not mm', 'work': 'algebraic SI work / normalized modal J equivalent'},
        'limits': ['no solver execution or rebuild', 'finite 0.2s record', 'contact predicted-stage force',
                   'nonorthogonal sampled-mass Gram', '25ms filter boundary sensitivity',
                   'AP diagnostic is conditional on selected replay frames'],
        'reproduce': f'OPENBLAS_NUM_THREADS=1 python3 -B tools/kanaya_diagnostics/three_checks.py --output-name NEW_NAME'})
    # A separately defined report builder uses only this run's computed results.
    report = build_report(out, coverage, static, dynamic, trends, cycles, equal, contributions,
                          spatial, flow, meta, validation)
    with (out/'report.md').open('x') as stream: stream.write(report)
    print(f'Completed: {out}/report.md', flush=True)


if __name__ == '__main__':
    main()
