#!/usr/bin/env python3
"""Batch main-peak-associated airflow threshold Jitter (default 3%).

CSV and CLI relative input paths resolve from the repository root, not the CWD.
See tools/README_jitter_batch.md for the manifest, statuses and area conventions.
"""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from zoneinfo import ZoneInfo

sys.dont_write_bytecode = True
import analyze_airflow_metrics as core

np = core.np
REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_COLUMNS = ('case_id', 'flow_path', 'run_path', 'area_path', 'stiffness_ratio',
                    'pressure_pa', 'label', 'dt_s', 'start_s', 'end_s', 'notes')
SUMMARY_COLUMNS = (
    'case_id', 'label', 'stiffness_ratio', 'pressure_pa', 'start_s', 'end_s', 'status',
    'event_method', 'flow_threshold_pct', 'peak_selector', 'n_acf_regions', 'n_main_peaks', 'n_with_3pct_event',
    'n_without_3pct_event', 'n_events', 'n_intervals', 'n_jitter_terms', 'n_rap_terms',
    'n_ppq5_terms', 'jitter_local_pct', 'rap_pct', 'ppq5_pct', 'mean_interval_ms',
    'sd_interval_ms', 'min_interval_ms', 'max_interval_ms', 'median_interval_ms',
    'period_cv_pct', 'first_event_s', 'last_event_s', 'event_span_fraction',
    'n_multiple_crossing_branches', 'n_reached_zero', 'n_reached_ten_pct',
    'jitter_status', 'rap_status', 'ppq5_status', 'area_status', 'area_reason',
    'n_area_aligned_events', 'area_event_valid_fraction', 'area_lag_median_ms',
    'area_lag_iqr_ms', 'area_min_median_mm2', 'settings_id', 'warnings', 'notes')
FAILURE_COLUMNS = ('case_id', 'status', 'path', 'phase', 'reason')
EVENT_COLUMNS = ('case_id', 'parent_peak_id', 'event_id', 'event_time_s', 'peak_time_s',
                 'peak_raw_time_s', 'Upeak_m3_s', 'threshold_m3_s', 'sample_left',
                 'sample_right', 'value_left', 'value_right', 'selection_rule')
PERIOD_COLUMNS = ('case_id', 'interval_id', 'event_start', 'event_end', 'start_s',
                  'end_s', 'interval_s', 'next_T_s', 'jitter_contribution_pct')
AREA_COLUMNS = ('case_id', 'parent_peak_id', 'has_3pct_event', 't_Ug_3pct_s',
                'peak_raw_time_s', 'branch_end_s', 'status', 'reason', 't_Amin_min_s',
                'Amin_min_mm2', 'min_section_index', 'min_section_column_zero_based',
                'lag_ms', 'Amin_at_3pct_mm2', 'amin_at_event_rule',
                'left_window_edge_region', 'right_censored')


def threshold_pct(case):
    return core.flow_threshold_percent(case.get('flow_threshold_pct', 3.))


def threshold_columns(columns, pct):
    suffix = core.flow_threshold_suffix(pct)
    return tuple(c.replace('3pct', suffix) for c in columns)


def threshold_label(pct):
    return core.flow_threshold_suffix(pct)[:-3].replace('p', '.')


def threshold_rule(pct):
    return (f'last real {threshold_label(pct)}% of Upeak downcross before first exact zero; '
            'without zero, before first observed global main-branch minimum')


def repo_path(value):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else REPO_ROOT / p).resolve()


def numeric(value, default=None, name='value'):
    if value is None or str(value).strip() in ('', 'NA'):
        return default
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be numeric') from exc
    if not math.isfinite(result):
        raise ValueError(f'{name} must be finite')
    return result


def load_cases(args):
    """Structural CSV errors fail before writing; row errors remain visible cases."""
    if args.manifest:
        with repo_path(args.manifest).open(newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames or []
            if len(set(header)) != len(header) or 'case_id' not in header:
                raise ValueError('manifest needs a unique header including case_id')
            if set(header) - set(MANIFEST_COLUMNS):
                raise ValueError('unknown manifest columns: ' + str(sorted(set(header) - set(MANIFEST_COLUMNS))))
            if not {'flow_path', 'run_path'} & set(header):
                raise ValueError('manifest needs flow_path or run_path')
            source_rows = list(reader)
        if not source_rows:
            raise ValueError('manifest has no cases')
    else:
        source_rows = []
        for index, value in enumerate(args.paths, 1):
            p = repo_path(value)
            slug = re.sub('[^A-Za-z0-9_-]+', '_', p.stem).strip('_')[:32] or 'run'
            suffix = hashlib.sha256(str(p).encode()).hexdigest()[:8]
            source_rows.append({'case_id': f'case_{index:03d}_{slug}_{suffix}',
                                'flow_path' if p.is_file() else 'run_path': value})
    cases = []
    seen = set()
    for row in source_rows:
        identity = row.get('case_id') or ''
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*(?:\.[A-Za-z0-9_-]+)*', identity):
            raise ValueError('unsafe/blank case_id (ASCII letters/digits, _, -, single dots only): ' + repr(identity))
        if identity in seen:
            raise ValueError('duplicate case_id: ' + identity)
        seen.add(identity)
        case = dict(row, case_id=identity, flow_threshold_pct=core.flow_threshold_percent(args.flow_threshold_pct), peak_selector=args.peak_selector, label=row.get('label') or identity,
                    notes=row.get('notes') or '', error=None)
        try:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('malformed CSV row: missing or extra fields')
            if bool(row.get('flow_path')) == bool(row.get('run_path')):
                raise ValueError('specify exactly one of flow_path and run_path')
            for name, default in (('dt_s', args.dt), ('start_s', args.start),
                                  ('end_s', args.end), ('stiffness_ratio', None), ('pressure_pa', None)):
                case[name] = numeric(row.get(name), default, name)
            if case['dt_s'] <= 0 or not 0 <= case['start_s'] < case['end_s']:
                raise ValueError('require dt_s > 0 and 0 <= start_s < end_s')
        except ValueError as exc:
            case['error'] = str(exc)
        cases.append(case)
    return cases


def resolve_flow(case):
    if case.get('flow_path'):
        path = repo_path(case['flow_path'])
    else:
        directory = repo_path(case['run_path'])
        if not directory.is_dir():
            raise ValueError('run_path is not a directory: ' + str(directory))
        matches = sorted({p.resolve() for p in directory.rglob('airflow_vt.dat') if p.is_file()})
        if len(matches) != 1:
            raise ValueError(f'run_path must contain exactly one airflow_vt.dat; found {len(matches)}: {directory}')
        path = matches[0]
    if not path.is_file() or path.name != 'airflow_vt.dat':
        raise ValueError('expected an existing airflow_vt.dat file: ' + str(path))
    return path


def read_flow(case):
    flow = core.read_flow(case['resolved_flow'], case['dt_s'])
    if np.any(flow['x'] < 0):
        raise ValueError('peak-associated-closure requires nonnegative Ug')
    if not flow['t'][0] <= case['start_s'] < case['end_s'] <= flow['t'][-1] + 1 / flow['fs'] + 1e-12:
        raise ValueError('flow does not cover the requested [start_s,end_s) window')
    core.closure_archive_provenance(flow, None)
    return flow


def read_area(case, flow):
    area = core.read_minimum_area(case['resolved_area'], case['dt_s'])
    meta = area['area_metadata']
    if meta['input_columns'] != 51 or meta['selected_area_columns_zero_based'] != list(range(2, 50)):
        raise ValueError('area requires step + 50 sections; interior DAT columns 2..49')
    if not area['t'][0] <= case['start_s'] < case['end_s'] <= area['t'][-1] + 1 / area['fs'] + 1e-12:
        raise ValueError('area does not cover the requested window')
    core.closure_archive_provenance(flow, area)
    return area


def preflight(cases):
    for case in cases:
        case['area_status'] = 'AREA_NOT_AVAILABLE'
        case['area_reason'] = ''
        if case['error']:
            continue
        try:
            case['resolved_flow'] = resolve_flow(case)
            flow = read_flow(case)
            provenance = core.closure_archive_provenance(flow, None)
            case.update(flow_sha256=flow['sha256'], flow_stride=flow['stride'],
                        flow_fs_hz=flow['fs'], flow_first_s=float(flow['t'][0]),
                        flow_last_s=float(flow['t'][-1]), archive_id=provenance.get('archive_id'),
                        archive_metrics_sha256=provenance.get('archive_metrics_sha256'))
        except (ValueError, OSError, RuntimeError, OverflowError) as exc:
            case['error'] = str(exc)
            continue
        try:
            area_path = repo_path(case['area_path']) if case.get('area_path') else case['resolved_flow'].with_name('area.dat')
            if case.get('area_path') or area_path.is_file():
                if area_path.name != 'area.dat' or area_path.parent != case['resolved_flow'].parent:
                    raise ValueError('area must be the area.dat in the SAME dat directory as the flow')
                case['resolved_area'] = area_path
                area = read_area(case, flow)
                case.update(area_status='AVAILABLE', area_sha256=area['sha256'],
                            area_stride=area['stride'], area_fs_hz=area['fs'])
        except (ValueError, OSError, RuntimeError, OverflowError) as exc:
            case.update(area_status='AREA_FAILED', area_reason=str(exc))
        print(f"{case['case_id']}: READY [{case['start_s']:g},{case['end_s']:g}) s "
              f"dt={case['dt_s']:g}, stride={case['flow_stride']}, {case['area_status']} -> {case['resolved_flow']}", flush=True)
    for case in cases:
        if case['error']:
            print(f"{case['case_id']}: FAILED preflight: {case['error']}", flush=True)


def area_alignment(case, flow, result):
    """Independent native-sample area observations; no feedback to any flow event."""
    suffix = core.flow_threshold_suffix(threshold_pct(case))
    has_event, flow_time, at_event = 'has_'+suffix+'_event', 't_Ug_'+suffix+'_s', 'Amin_at_'+suffix+'_mm2'
    area = None
    if case['area_status'] == 'AVAILABLE':
        try:
            area = read_area(case, flow)
            if area['sha256'] != case['area_sha256']:
                raise ValueError('area changed since preflight')
            matrix = np.loadtxt(case['resolved_area'], comments='#', ndmin=2)
            if core.sha256(case['resolved_area']) != area['sha256']:
                raise ValueError('area changed while reading section indices')
            values = matrix[:, 2:50]
            if area['area_metadata']['input_unit'] == 'm2':
                values = values * 1e6
            # read_minimum_area already validates finite interior coverage and signs.
            winners = np.argmin(np.where(np.isfinite(values), values, np.inf), axis=1)
        except (ValueError, OSError, RuntimeError, OverflowError) as exc:
            case.update(area_status='AREA_FAILED', area_reason=str(exc))
            area = None
    events = {e['peak_id']: e for e in result['series']['main_'+suffix]}
    cycles = {c['main_peak_id']: c for c in result['cycles'] if c['main_peak_id']}
    rows = []
    for peak in result['peaks']:
        if peak['peak_role'] != 'MAIN':
            continue
        identity = peak['event_id']
        event = events.get(identity)
        left = peak['raw_time_s']
        right = result['branches'][identity]['branch_end_s']
        row = dict(case_id=case['case_id'], parent_peak_id=identity,
                   peak_raw_time_s=left, branch_end_s=right, status=case['area_status'],
                   reason=case['area_reason'], left_window_edge_region=peak['search_region_id'] == 0,
                   right_censored=cycles[identity]['right_censored'])
        row.update({has_event: event is not None, flow_time: event['event_time_s'] if event else None,
                    'flow_threshold_pct': threshold_pct(case)})
        rows.append(row)
        if area is None:
            continue
        at, ax = area['t'], area['x']
        ids = np.flatnonzero((at >= left) & (at < right))
        if not len(ids):
            row.update(status='AREA_NO_NATIVE_SAMPLE_IN_BRANCH', reason='no area sample in main branch')
            continue
        minimum = float(ax[ids].min())
        ties = ids[ax[ids] == minimum]
        k = int(ties[0])  # Exact ties: earliest native time, then lowest section column.
        plateau_end = k
        while plateau_end + 1 <= ids[-1] and ax[plateau_end + 1] == minimum:
            plateau_end += 1
        columns = np.flatnonzero(values[k] == minimum) + 2
        row.update(status='ALIGNED', t_Amin_min_s=float(at[k]), Amin_min_mm2=minimum,
                   min_section_index=int(winners[k]) + 1, min_section_column_zero_based=int(winners[k]) + 2,
                   min_section_ties=';'.join(str(int(c) - 1) for c in columns),
                   min_equal_sample_count=len(ties), min_plateau_start_s=float(at[k]),
                   min_plateau_end_s=float(at[plateau_end]), min_plateau_samples=plateau_end - k + 1,
                   min_section_switch_count=int(np.count_nonzero(np.diff(winners[ids]))),
                   rows_with_undefined_section=int((~np.isfinite(values[ids])).any(axis=1).sum()),
                   area_partial_left=bool(left < at[0]), area_partial_right=bool(right > at[-1] + 1 / area['fs'] + 1e-12))
        if event is not None:
            et = event['event_time_s']
            row['lag_ms'] = float(1000 * (at[k] - et))
            j = int(np.searchsorted(at, et))
            if j < len(at) and at[j] == et:
                row.update({at_event: float(ax[j]), 'amin_at_event_rule': 'EXACT_NATIVE_AMIN_SAMPLE'})
            elif 0 < j < len(at):
                row.update({at_event: float(ax[j-1] + (ax[j]-ax[j-1]) * (et-at[j-1]) / (at[j]-at[j-1])),
                            'amin_at_event_rule': 'LINEAR_NATIVE_AMIN_BRACKET'})
            else:
                row.update({at_event: None, 'amin_at_event_rule': 'OUTSIDE_SAMPLED_RANGE_NO_EXTRAPOLATION'})
    if area is not None and core.sha256(case['resolved_area']) != area['sha256']:
        raise RuntimeError('area changed during alignment')
    return rows


def initial_summary(case, settings_id):
    pct = threshold_pct(case)
    row = {**{k: None for k in threshold_columns(SUMMARY_COLUMNS, pct)},
            **{k: case.get(k) for k in ('case_id', 'label', 'stiffness_ratio', 'pressure_pa', 'start_s', 'end_s', 'notes')},
            'status': 'PENDING', 'event_method': 'peak-associated-closure/main_'+core.flow_threshold_suffix(pct),
            'flow_threshold_pct': pct, 'peak_selector': case.get('peak_selector','prominence'), 'settings_id': settings_id}
    for key in ('stiffness_ratio', 'pressure_pa', 'start_s', 'end_s'):
        try:
            row[key] = numeric(row.get(key))
        except ValueError:
            row[key] = None
    return row


def analyze_case(case, settings_id):
    pct = threshold_pct(case)
    suffix = core.flow_threshold_suffix(pct)
    method = 'main_'+suffix
    flow = read_flow(case)
    if flow['sha256'] != case['flow_sha256']:
        raise ValueError('flow changed since preflight')
    # These are the same three calls as the single-run implementation. Area is
    # deliberately omitted here; the Ug detector never consumes area information.
    base = core.acf_guided_extract(flow, None, case['start_s'], case['end_s'])
    result = core.peak_associated_extract(flow, case['start_s'], case['end_s'], base, flow_threshold_pct=pct, peak_selector=case.get('peak_selector','prominence'))
    events = result['series'][method]
    metric, periods, terms = core.all_event_metrics(events, method, case['start_s'], case['end_s'])
    summary = initial_summary(case, settings_id)
    main_peaks = [p for p in result['peaks'] if p['peak_role'] == 'MAIN']
    summary.update(status='COMPLETED' if all(metric[k] is not None for k in core.METRIC_NAMES) else 'INSUFFICIENT_EVENTS',
                   n_acf_regions=len(result['cycles']), n_main_peaks=len(main_peaks),
                   n_events=len(events), n_intervals=metric['n_intervals'], n_jitter_terms=metric['n_jitter_terms'],
                   n_rap_terms=metric['n_RAP_terms'], n_ppq5_terms=metric['n_PPQ5_terms'],
                   jitter_local_pct=metric['jitter_local_percent'], rap_pct=metric['rap_percent'], ppq5_pct=metric['ppq5_percent'],
                   period_cv_pct=metric['period_cv_percent'], first_event_s=metric['first_event_s'],
                   last_event_s=metric['last_event_s'], event_span_fraction=metric['event_span_fraction'],
                   n_multiple_crossing_branches=sum(r['peak_role'] == 'MAIN' and r['threshold_name'] == suffix and r['reach_count'] > 1 for r in result['reach_history']),
                   n_reached_zero=len(result['series']['main_zero']), n_reached_ten_pct=len(result['series']['main_10pct']))
    summary.update({'n_with_'+suffix+'_event': len(events), 'n_without_'+suffix+'_event': len(main_peaks)-len(events)})
    for target, source in (('mean_interval_ms', 'mean_period_s'), ('sd_interval_ms', 'std_interval_s'),
                           ('min_interval_ms', 'min_interval_s'), ('max_interval_ms', 'max_interval_s'),
                           ('median_interval_ms', 'median_period_s')):
        summary[target] = metric[source] * 1000 if metric[source] is not None else None
    for key, field in zip(('jitter_status', 'rap_status', 'ppq5_status'), core.METRIC_NAMES):
        summary[key] = 'DEFINED' if metric[field] is not None else 'INSUFFICIENT_EVENTS'
    try:
        aligned = area_alignment(case, flow, result)
    except (ValueError, OSError, RuntimeError, OverflowError, ArithmeticError) as exc:
        case.update(area_status='AREA_FAILED', area_reason=str(exc))
        aligned = area_alignment(case, flow, result)
    valid = [r for r in aligned if r['has_'+suffix+'_event'] and r['status'] == 'ALIGNED']
    lags = [r['lag_ms'] for r in valid]
    summary.update(area_status=case['area_status'], area_reason=case['area_reason'], n_area_aligned_events=len(valid),
                   area_event_valid_fraction=len(valid) / len(events) if events and case['area_status'] == 'AVAILABLE' else None,
                   area_lag_median_ms=float(np.median(lags)) if lags else None,
                   area_lag_iqr_ms=float(np.percentile(lags, 75) - np.percentile(lags, 25)) if lags else None,
                   area_min_median_mm2=float(np.median([r['Amin_min_mm2'] for r in valid])) if valid else None)
    warnings = [case['area_status']] if case['area_status'] != 'AVAILABLE' else []
    if summary['n_without_'+suffix+'_event']:
        warnings.append('MAIN_BRANCHES_WITHOUT_'+suffix.upper()+'_EVENT')
    if summary['n_multiple_crossing_branches']:
        warnings.append('MULTIPLE_REAL_DOWNCROSSINGS_FIXED_RULE')
    if summary['status'] == 'INSUFFICIENT_EVENTS':
        warnings.append('INSUFFICIENT_EVENTS')
    if main_peaks:
        warnings.append('FINAL_BRANCH_RIGHT_CENSORED')
    summary['warnings'] = ';'.join(warnings)
    event_rows = [dict(e, case_id=case['case_id'], parent_peak_id=e['peak_id']) for e in events]
    jitter_terms = {r['first_interval_id']: r for r in terms if r['metric'] == 'jitter_local_percent'}
    period_rows = []
    for i, p in enumerate(periods):
        period_rows.append(dict(p, case_id=case['case_id'], next_T_s=periods[i+1]['interval_s'] if i+1 < len(periods) else None,
                                jitter_contribution_pct=jitter_terms[i]['normalized_term_percent'] if i in jitter_terms else None))
    if core.sha256(case['resolved_flow']) != flow['sha256']:
        raise RuntimeError('flow changed during analysis')
    return {'summary': [summary], 'events_'+suffix: event_rows, 'periods_'+suffix: period_rows,
            'metric_terms_'+suffix: [dict(r, case_id=case['case_id']) for r in terms],
            'area_alignment': aligned, 'peaks': [dict(p, case_id=case['case_id']) for p in result['peaks']]}


def csv_columns(rows, required=()):
    return list(dict.fromkeys([*required, *(k for row in rows for k in row)]))


def csv_records(rows):
    return [dict(r, **{k: 'NA' for k in ('stiffness_ratio', 'pressure_pa') if k in r and r[k] is None}) for r in rows]


def atomic_csv(path, rows, required=()):
    with tempfile.NamedTemporaryFile(prefix='.' + path.name, dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
    try:
        core.write_csv(temporary, csv_records(rows), csv_columns(rows, required))
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def resolved_row(case, identifiers):
    pct = threshold_pct(case)
    return dict(case_id=case['case_id'], requested_flow_path=case.get('flow_path'), requested_run_path=case.get('run_path'),
                flow_path=str(case.get('resolved_flow', '')), area_path=str(case.get('resolved_area', '')),
                dt_s=case.get('dt_s'), start_s=case.get('start_s'), end_s=case.get('end_s'),
                flow_sha256=case.get('flow_sha256'), flow_stride=case.get('flow_stride'), fs_hz=case.get('flow_fs_hz'),
                flow_first_s=case.get('flow_first_s'), flow_last_s=case.get('flow_last_s'),
                area_sha256=case.get('area_sha256'), area_stride=case.get('area_stride'), area_fs_hz=case.get('area_fs_hz'),
                area_status=case['area_status'], area_reason=case['area_reason'], preflight_error=case['error'],
                source_run_directory=str(case['resolved_flow'].parent) if case.get('resolved_flow') else '',
                archive_id=case.get('archive_id'), archive_metrics_sha256=case.get('archive_metrics_sha256'),
                label=case['label'], stiffness_ratio=case.get('stiffness_ratio'), pressure_pa=case.get('pressure_pa'),
                flow_threshold_pct=pct, event_method='peak-associated-closure/main_'+core.flow_threshold_suffix(pct),
                primary_selection_rule=threshold_rule(pct), peak_selector=case.get('peak_selector','prominence'), **identifiers)


def comparison_groups(summaries):
    groups = {}
    for row in summaries:
        if row.get('stiffness_ratio') is None:
            continue
        pressure = row.get('pressure_pa')
        key = (pressure, row.get('start_s'), row.get('end_s'), (row['settings_id'], row.get('flow_threshold_pct', 3.), row.get('peak_selector','prominence')), row['case_id'] if pressure is None else None)
        groups.setdefault(key, []).append(row)
    return groups


def comparison_plot(out, summaries):
    groups = comparison_groups(summaries)
    if not any(r.get('jitter_local_pct') is not None for rows in groups.values() for r in rows):
        return
    (out / 'figures').mkdir(exist_ok=True)
    fig, ax = core.plt.subplots(figsize=(9, 5))
    for i, (key, rows) in enumerate(groups.items(), 1):
        pressure, start, end, _, unknown_case = key
        rows = sorted(rows, key=lambda r: r['stiffness_ratio'])
        xs = [r['stiffness_ratio'] for r in rows]
        ys = [r['jitter_local_pct'] if r.get('jitter_local_pct') is not None else np.nan for r in rows]
        label = f'P={pressure if pressure is not None else "NA"} Pa; [{start},{end}) s; group {i}'
        if unknown_case:
            label += '; ' + unknown_case
        style = 'o-' if pressure is not None and len(set(xs)) == len(xs) else 'o'
        ax.plot(xs, ys, style, ms=5, lw=.8, label=label)
    pct = summaries[0].get('flow_threshold_pct', 3.)
    ax.set(xlabel='Stiffness ratio (manifest)', ylabel=f'Main-associated {threshold_label(pct)}% Jitter local [%]')
    ax.grid(alpha=.15)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / 'figures/jitter_by_stiffness_ratio.png', dpi=150)
    core.plt.close(fig)


def report_text(summaries, cases, state, pct=3.):
    suffix = core.flow_threshold_suffix(pct)
    label = threshold_label(pct)
    def esc(value):
        return str(value if value is not None else 'NA').replace('|', '\\|').replace('\n', ' ')
    def number(value):
        return f'{value:.6f}' if isinstance(value, (int, float)) and math.isfinite(value) else 'NA'
    lines = [f'# Main-peak-associated {label}% airflow Jitter batch', '', 'Batch state: ' + state + '.', '',
             'The peak detector, crossing-selection rule and all-interval metric kernel are reused. Primary: peak-associated-closure/main_' + suffix + '.',
             threshold_rule(pct) + '. No missing events, long gaps, or short gaps are filled, split or dropped.',
             'Jitter/RAP/PPQ5 use all N−1/N−2/N−4 terms. Empty values mean undefined, never 0%. Manual validation is pending.', '',
             f'| Case | Label | Ratio | Ps [Pa] | Window [s] | Status | Jitter [%] | RAP [%] | PPQ5 [%] | Events | Main peaks | No {label}% | Mean [ms] | Max [ms] |',
             '|---|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in summaries:
        fields = [r['case_id'], r['label'], r.get('stiffness_ratio'), r.get('pressure_pa'),
                  f"[{r.get('start_s')},{r.get('end_s')})", r['status'],
                  *[number(r.get(k)) for k in ('jitter_local_pct', 'rap_pct', 'ppq5_pct')],
                  r.get('n_events'), r.get('n_main_peaks'), r.get('n_without_'+suffix+'_event'),
                  number(r.get('mean_interval_ms')), number(r.get('max_interval_ms'))]
        lines.append('| ' + ' | '.join(esc(x) for x in fields) + ' |')
    windows = {(r.get('start_s'), r.get('end_s')) for r in summaries}
    lines += ['', f'Comparison groups: {len(comparison_groups(summaries))}; distinct windows: {len(windows)}.',
              'Different pressure/window/settings never share a trend line; unknown pressure is a separate case. Unknown ratio is listed in input order only.',
              '`event_span_fraction` is only (last−first event)/(end−start), not search coverage, closure fraction or Jitter time coverage.', '',
              '## Area observations', '',
              'Area is optional and never participates in flow event selection. Only same-directory area.dat is accepted.',
              'Use finite interior DAT columns 2..49 (48 of 50 sections), convert declared units to mm². NaN is not zero.',
              'Native branch minimum: earliest exact minimum sample, then lowest tied section. min_section_index is 0-based among the 50 sections (DAT column minus 1).',
              f'lag_ms = 1000*(area minimum time − {label}% flow time); positive means area minimum later. Amin_at_{suffix} is a declared linear interpolation of native Amin brackets, without extrapolation/resampling.',
              f'Main peaks without {label}% events also have area rows with blank flow time/lag. Equal minima, plateaus, section switches, undefined sections and truncated branches are recorded. Proximity alone does not prove closure.', '',
              '## Inputs and warnings', '']
    for case, row in zip(cases, summaries):
        lines.append('- ' + esc(case['case_id']) + ': ' + esc(case.get('resolved_flow', case.get('flow_path') or case.get('run_path'))) + '; ' + esc(row.get('warnings')))
        if case.get('area_reason'):
            lines.append('  Area: ' + esc(case['area_reason']))
    lines += ['', '## Files', '',
              '- summary.csv: one row per requested case, in input order; insufficient/failed cases remain visible.',
              '- inputs_resolved.csv: resolved inputs, SHA256, actual sampling, mode/rule/settings and code identifiers.',
              '- failures.csv: explicit failures/skips; always has a header.',
              f'- events_{suffix}.csv / periods_{suffix}.csv: all events/intervals, next_T_s and every J_i (last interval has no next term).',
              f'- metric_terms_{suffix}.csv: every Jitter/RAP/PPQ5 numerator and normalized term.',
              '- area_alignment.csv / peaks.csv: independent area observations and all main/secondary clear peaks.',
              '- cases/<case_id>/: atomically completed numerical case results. No raw waveform copies or default per-case PNG/JSON.',
              '- figures/jitter_by_stiffness_ratio.png: only when numeric manifest ratios and defined Jitter exist.', '']
    return '\n'.join(lines)


def save_state(out, cases, summaries, failures, aggregate, identifiers, state):
    pct = threshold_pct(cases[0])
    suffix = core.flow_threshold_suffix(pct)
    atomic_csv(out / 'summary.csv', summaries, threshold_columns(SUMMARY_COLUMNS, pct))
    atomic_csv(out / 'inputs_resolved.csv', [resolved_row(c, identifiers) for c in cases])
    atomic_csv(out / 'failures.csv', failures, FAILURE_COLUMNS)
    required = {'events_'+suffix: EVENT_COLUMNS, 'periods_'+suffix: PERIOD_COLUMNS, 'area_alignment': threshold_columns(AREA_COLUMNS, pct),
                'metric_terms_'+suffix: ('case_id', 'metric', 'term_id'), 'peaks': ('case_id', 'event_id', 'peak_role')}
    for name, rows in aggregate.items():
        atomic_csv(out / (name + '.csv'), rows, required[name])
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', prefix='.report_', dir=out, delete=False) as f:
        f.write(report_text(summaries, cases, state, pct))
        temporary = Path(f.name)
    temporary.replace(out / 'report.md')


def selection_event_differences(reference, variant):
    """Time-first comparison only; this never changes either metric event vector."""
    pairs = core.optimal_event_assignment([e['event_time_s'] for e in reference],
                                          [e['event_time_s'] for e in variant], 1e-12)
    used_a = {a for a,b,dt in pairs}; used_b = {b for a,b,dt in pairs}
    rows = []
    def record(a, b, kind):
        rows.append(dict(change=kind, A_event_id=a['event_id'] if a else None,
            variant_event_id=b['event_id'] if b else None,
            A_region_id=a['search_region_id'] if a else None,
            variant_region_id=b['search_region_id'] if b else None,
            A_peak_id=a['peak_id'] if a else None, variant_peak_id=b['peak_id'] if b else None,
            A_event_time_s=a['event_time_s'] if a else None,
            variant_event_time_s=b['event_time_s'] if b else None,
            time_change_ms=1000*(b['event_time_s']-a['event_time_s']) if a and b else None,
            time_equivalence_tolerance_s=2e-12))
    for a,b,dt in pairs:
        ea,eb=reference[a],variant[b]
        record(ea,eb,'SAME_TIME_SAME_PARENT' if ea['peak_id']==eb['peak_id'] else 'SAME_TIME_REASSIGNED')
    remaining_b={e['search_region_id']:j for j,e in enumerate(variant) if j not in used_b}
    for i,a in enumerate(reference):
        if i in used_a:continue
        j=remaining_b.pop(a['search_region_id'],None)
        if j is None:record(a,None,'MISSING')
        else:
            record(a,variant[j],'TIME_CHANGED_IN_REGION');used_b.add(j)
    for j,b in enumerate(variant):
        if j not in used_b:record(None,b,'ADDED')
    return rows


def selection_regression(case, pct, result, metric, baseline_root):
    suffix=core.flow_threshold_suffix(pct)
    folder=baseline_root/('allRatio_001' if pct==3 else 'allRatio_'+suffix+'_001')
    def rows(name):
        with (folder/name).open(newline='') as f:return [r for r in csv.DictReader(f) if r['case_id']==case['case_id']]
    prior=rows('events_'+suffix+'.csv');current=result['series']['main_'+suffix]
    if [(r['event_id'],r['parent_peak_id'],float(r['event_time_s'])) for r in prior]!=[(r['event_id'],r['peak_id'],r['event_time_s']) for r in current]:
        raise RuntimeError('A baseline event regression failed: '+case['case_id']+' '+suffix)
    saved=rows('summary.csv');inputs=rows('inputs_resolved.csv')
    if len(saved)!=1 or len(inputs)!=1:raise ValueError('baseline needs one matching summary/input row')
    if inputs[0]['flow_sha256']!=case['flow_sha256'] or float(inputs[0]['dt_s'])!=case['dt_s'] or float(inputs[0]['start_s'])!=case['start_s'] or float(inputs[0]['end_s'])!=case['end_s']:
        raise ValueError('baseline input/dt/window mismatch')
    for key in core.METRIC_NAMES:
        column={'jitter_local_percent':'jitter_local_pct','rap_percent':'rap_pct','ppq5_percent':'ppq5_pct'}[key]
        value=float(saved[0][column]) if saved[0][column] else None
        if metric[key]!=value:raise RuntimeError('A baseline metric regression failed: '+key)
    prior_periods=rows('periods_'+suffix+'.csv')
    if [float(p['interval_s']) for p in prior_periods]!=list(np.diff([e['event_time_s'] for e in current])):
        raise RuntimeError('A baseline interval regression failed')
    return dict(case_id=case['case_id'],flow_threshold_pct=pct,n_events=len(current),n_intervals=len(prior_periods),
                ids_parent_ids_times_intervals_metrics='EXACT',max_time_difference_s=0.)


def run_peak_selection_comparison(args,cases,out):
    """Frozen A/B/C study using the same branch and metric kernels, including single-region interventions."""
    if args.peak_selector!='prominence':raise ValueError('--compare-peak-selectors always evaluates all three selectors; leave --peak-selector at its default')
    if any(c['error'] for c in cases):raise ValueError('comparison requires all manifest cases to pass preflight')
    percentages=list(dict.fromkeys(args.comparison_thresholds))
    baseline_root=repo_path(args.baseline_root)
    for pct in percentages:
        suffix=core.flow_threshold_suffix(pct)
        folder=baseline_root/('allRatio_001' if pct==3 else 'allRatio_'+suffix+'_001')
        for name in ('summary.csv','inputs_resolved.csv','events_'+suffix+'.csv','periods_'+suffix+'.csv'):
            if not (folder/name).is_file():raise ValueError('missing preserved A baseline: '+str(folder/name))
    out.mkdir(parents=True,exist_ok=False);(out/'figures').mkdir()
    (out/'report.md').write_text('# Peak-selection sensitivity\n\nIN_PROGRESS; incomplete study is not a completed result.\n')
    selectors={'A':'prominence','B':'clear-flow','C':'raw-flow'}
    tables={n:[] for n in ('comparison','main_peak_changes','events_diff','region_effects','single_region_impacts',
        'events','periods','metric_terms','peaks','search_regions','raw_peak_candidates','boundary_suspicions','regression_checks')}
    sources={};figure_data={}
    def remember(path):
        path=Path(path);sources[str(path)]=core.sha256(path)
    remember(__file__);remember(core.__file__)
    if args.manifest:remember(repo_path(args.manifest))
    for pct in percentages:
        suffix=core.flow_threshold_suffix(pct);folder=baseline_root/('allRatio_001' if pct==3 else 'allRatio_'+suffix+'_001')
        for name in ('summary.csv','inputs_resolved.csv','events_'+suffix+'.csv','periods_'+suffix+'.csv'):remember(folder/name)
    for index,case in enumerate(cases,1):
        print(f"Sensitivity {index}/{len(cases)}: {case['case_id']}",flush=True)
        flow=read_flow(case);remember(flow['path'])
        if flow['sha256']!=case['flow_sha256']:raise RuntimeError('input changed since preflight')
        start,end=case['start_s'],case['end_s'];base=core.acf_guided_extract(flow,None,start,end)
        regions=[r for r in base['regions'] if r['signal']=='Ug'];raw=[p for p in base['peak_candidates'] if p['method']=='Ug_max' and p['inside_window']]
        context=dict(case_id=case['case_id'],stiffness_ratio=case['stiffness_ratio'],pressure_pa=case['pressure_pa'])
        tables['raw_peak_candidates'].extend(dict(p,**context) for p in raw)
        tables['search_regions'].extend(dict(r,**context) for r in regions)
        # A single extraction with e.g. 5% already contains the fixed 3%/10% references.
        custom=[p for p in percentages if p not in (3,10)];extract_pcts=custom or [3.]
        runs={};metrics={};mains={}
        for label,selector in selectors.items():
            for p in extract_pcts:
                result=core.peak_associated_extract(flow,start,end,base,flow_threshold_pct=p,peak_selector=selector)
                for pct in percentages:
                    suffix=core.flow_threshold_suffix(pct)
                    if 'main_'+suffix in result['series'] and (label,pct) not in runs:runs[label,pct]=result
            result=runs[label,percentages[0]]
            mains[label]={p['search_region_id']:p for p in result['peaks'] if p['peak_role']=='MAIN'}
            tables['peaks'].extend(dict(p,selector=label,**context) for p in result['peaks'])
            for pct in percentages:
                result=runs[label,pct];suffix=core.flow_threshold_suffix(pct);es=result['series']['main_'+suffix]
                row,periods,terms=core.all_event_metrics(es,'main_'+suffix,start,end)
                metrics[label,pct]=row
                tables['events'].extend(dict(e,selector=label,flow_threshold_pct=pct,**context) for e in es)
                tables['periods'].extend(dict(e,selector=label,flow_threshold_pct=pct,**context) for e in periods)
                tables['metric_terms'].extend(dict(e,selector=label,flow_threshold_pct=pct,**context) for e in terms)
                if label=='A':tables['regression_checks'].append(selection_regression(case,pct,result,row,baseline_root))
        differences={};effects={};changed={}
        for label in ('B','C'):
            changed[label]=[r['region_id'] for r in regions if (mains['A'].get(r['region_id']) or {}).get('event_id')!=(mains[label].get(r['region_id']) or {}).get('event_id')]
        for pct in percentages:
            old_result=runs['A',pct];old_events=old_result['series']['main_'+core.flow_threshold_suffix(pct)]
            for label in selectors:
                result=runs[label,pct];es=result['series']['main_'+core.flow_threshold_suffix(pct)]
                diff=selection_event_differences(old_events,es);differences[label,pct]=diff
                tables['events_diff'].extend(dict(d,selector=label,flow_threshold_pct=pct,**context) for d in diff)
                oldmap={e['search_region_id']:e for e in old_events};newmap={e['search_region_id']:e for e in es}
                for reg in regions:
                    rid=reg['region_id'];pa,pb=mains['A'].get(rid),mains[label].get(rid);ea,eb=oldmap.get(rid),newmap.get(rid)
                    old_end=old_result['branches'][pa['event_id']]['branch_end_s'] if pa else None
                    new_end=result['branches'][pb['event_id']]['branch_end_s'] if pb else None
                    status='NOT_REACHED_BOTH' if not ea and not eb else 'ADDED' if not ea else 'MISSING' if not eb else 'SAME_TIME' if abs(ea['event_time_s']-eb['event_time_s'])<=2e-12 else 'TIME_CHANGED'
                    row=dict(**context,selector=label,flow_threshold_pct=pct,region_id=rid,region_start_s=reg['start_s'],region_end_s=reg['end_s'],
                        A_peak_id=pa['event_id'] if pa else None,variant_peak_id=pb['event_id'] if pb else None,
                        peak_changed=(pa or {}).get('event_id')!=(pb or {}).get('event_id'),
                        A_branch_end_s=old_end,variant_branch_end_s=new_end,branch_end_changed=old_end!=new_end,
                        A_event_time_s=ea['event_time_s'] if ea else None,variant_event_time_s=eb['event_time_s'] if eb else None,
                        event_change=status,event_time_change_ms=1000*(eb['event_time_s']-ea['event_time_s']) if ea and eb else None)
                    effects[label,pct,rid]=row
                    if row['peak_changed'] or row['branch_end_changed'] or status not in ('NOT_REACHED_BOTH','SAME_TIME'):tables['region_effects'].append(row)
                metric=metrics[label,pct];old=metrics['A',pct]
                deltas={key+'_delta_pp':metric[key]-old[key] if metric[key] is not None and old[key] is not None else None for key in core.METRIC_NAMES}
                tables['comparison'].append(dict(metric,**context,selector=label,peak_selector=selectors[label],flow_threshold_pct=pct,
                    n_acf_regions=len(regions),n_main_peaks=len(mains[label]),n_changed_main_regions=0 if label=='A' else len(changed[label]),
                    n_changed_existing_main_peaks=sum(rid in mains['A'] and rid in mains[label] and mains['A'][rid]['event_id']!=mains[label][rid]['event_id'] for rid in mains['A']),
                    n_added_main_peaks=len(set(mains[label])-set(mains['A'])),n_missing_main_peaks=len(set(mains['A'])-set(mains[label])),
                    n_event_added=sum(d['change']=='ADDED' for d in diff),n_event_missing=sum(d['change']=='MISSING' for d in diff),
                    n_event_time_changed=sum(d['change']=='TIME_CHANGED_IN_REGION' for d in diff),n_event_time_same_reassigned=sum(d['change']=='SAME_TIME_REASSIGNED' for d in diff),
                    **deltas,diagnostic_only=label!='A'))
        # Change one observed region at a time, including recomputing the previous branch.
        # These effects are not additive and are not physical causal attribution.
        for label in ('B','C'):
            for rid in changed[label]:
                pa,pb=mains['A'].get(rid),mains[label].get(rid)
                if pb is None:raise RuntimeError('unexpected disappearance of all candidates')
                single_runs={}
                for p in extract_pcts:
                    single=core.peak_associated_extract(flow,start,end,base,flow_threshold_pct=p,main_peak_overrides={rid:pb['event_id']})
                    for pct in percentages:
                        if 'main_'+core.flow_threshold_suffix(pct) in single['series'] and pct not in single_runs:single_runs[pct]=single
                change=dict(**context,selector=label,region_id=rid,region_start_s=regions[rid]['start_s'],region_end_s=regions[rid]['end_s'],
                    A_peak_id=pa['event_id'] if pa else None,variant_peak_id=pb['event_id'],
                    A_peak_time_s=pa['event_time_s'] if pa else None,variant_peak_time_s=pb['event_time_s'],
                    A_Upeak_m3_s=pa['sample_value'] if pa else None,variant_Upeak_m3_s=pb['sample_value'],
                    A_prominence=pa['prominence'] if pa else None,variant_prominence=pb['prominence'],
                    A_clarity_accepted=pa['clarity_accepted'] if pa else None,variant_clarity_accepted=pb['clarity_accepted'],
                    variant_clarity_rejection_reasons=pb['detector_rejection_reasons'])
                for pct in percentages:
                    suffix=core.flow_threshold_suffix(pct);single=single_runs[pct]
                    row,_,_=core.all_event_metrics(single['series']['main_'+suffix],'main_'+suffix,start,end)
                    original=metrics['A',pct];delta=row['jitter_local_percent']-original['jitter_local_percent'] if row['jitter_local_percent'] is not None and original['jitter_local_percent'] is not None else None
                    diff=selection_event_differences(runs['A',pct]['series']['main_'+suffix],single['series']['main_'+suffix])
                    joint=metrics[label,pct]['jitter_local_percent']-original['jitter_local_percent'] if metrics[label,pct]['jitter_local_percent'] is not None and original['jitter_local_percent'] is not None else None
                    tables['single_region_impacts'].append(dict(**context,selector=label,region_id=rid,flow_threshold_pct=pct,
                        replacement_peak_id=pb['event_id'],jitter_local_percent=row['jitter_local_percent'],single_region_delta_pp=delta,
                        n_events=row['n_events'],n_added=sum(d['change']=='ADDED' for d in diff),n_missing=sum(d['change']=='MISSING' for d in diff),
                        n_time_changed=sum(d['change']=='TIME_CHANGED_IN_REGION' for d in diff),joint_selector_delta_pp=joint))
                    effect=effects[label,pct,rid]
                    change.update({'A_'+suffix+'_event_s':effect['A_event_time_s'],'variant_'+suffix+'_event_s':effect['variant_event_time_s'],
                        suffix+'_event_change':effect['event_change'],suffix+'_single_region_delta_pp':delta,
                        suffix+'_joint_selector_delta_pp':joint})
                tables['main_peak_changes'].append(change)
        # Morphological diagnostics use real extrema and existing T0 only as a scale;
        # no warning changes event adoption or the ACF partition.
        raw_by_region={r['region_id']:[p for p in raw if p['search_region_id']==r['region_id']] for r in regions}
        for label in selectors:
            ordered=sorted(mains[label].values(),key=lambda p:p['event_time_s'])
            for a,b in zip(ordered[:-1],ordered[1:]):
                dt=b['event_time_s']-a['event_time_s'];T=regions[a['search_region_id']]['search_scale_s'];lo,hi=a['sample_id'],b['sample_id']
                trough=float(flow['x'][lo:hi+1].min());trough_ratio=trough/min(a['sample_value'],b['sample_value']) if min(a['sample_value'],b['sample_value'])>0 else None
                flags=[]
                if dt<.55*T:flags.append('CLOSE_MAIN_EVENTS_RELATIVE_TO_T0')
                if trough_ratio is not None and trough_ratio>.5:flags.append('ADJACENT_MAINS_WITHOUT_DEEP_INTERVENING_TROUGH')
                if flags:tables['boundary_suspicions'].append(dict(**context,selector=label,kind=';'.join(flags),peak_id=a['event_id'],other_peak_id=b['event_id'],time_s=a['event_time_s'],other_time_s=b['event_time_s'],region_id=a['search_region_id'],other_region_id=b['search_region_id'],interval_over_T0=dt/T,intervening_trough_ratio=trough_ratio,diagnostic_only=True))
            for rid,p in mains[label].items():
                near=[q for j in (rid-1,rid+1) for q in raw_by_region.get(j,[]) if q['sample_value']>p['sample_value'] and abs(q['event_time_s']-p['event_time_s'])<.55*regions[rid]['search_scale_s']]
                if near:
                    q=max(near,key=lambda q:q['sample_value'])
                    tables['boundary_suspicions'].append(dict(**context,selector=label,kind='SMALL_MAIN_NEAR_LARGER_PEAK_ACROSS_REGION_BOUNDARY',peak_id=p['event_id'],other_peak_id=q['event_id'],time_s=p['event_time_s'],other_time_s=q['event_time_s'],region_id=rid,other_region_id=q['search_region_id'],interval_over_T0=abs(q['event_time_s']-p['event_time_s'])/regions[rid]['search_scale_s'],diagnostic_only=True))
        if case['stiffness_ratio'] in (.45,.55,.60):figure_data[case['case_id']]=(case,flow,base,runs,mains)
        for name,rows in tables.items():atomic_csv(out/(name+'.csv'),rows,('case_id',))
        print(case['case_id'],'A/B/C main peaks',*[len(mains[k]) for k in selectors],flush=True)
    for path,digest in sources.items():
        if core.sha256(path)!=digest:raise RuntimeError('protected input/baseline/source changed during comparison: '+path)
    atomic_csv(out/'sources.csv',[dict(path=p,sha256=d) for p,d in sorted(sources.items())])
    selection_comparison_figures(out,tables,figure_data,percentages)
    (out/'report.md').write_text(selection_comparison_report(tables,percentages,len(cases)))
    print('Saved sensitivity study:',out,'rows=',len(tables['comparison']))
    return 0


def selection_comparison_figures(out,tables,figure_data,percentages):
    styles={'A':('#1764ab','^'),'B':('#d87813','s'),'C':('#15816b','D')}
    rows=tables['comparison']
    fig,axes=core.plt.subplots(1,len(percentages),figsize=(7*len(percentages),4.8),squeeze=False)
    for ax,pct in zip(axes[0],percentages):
        for label,(color,marker) in styles.items():
            rr=sorted((r for r in rows if r['selector']==label and r['flow_threshold_pct']==pct),key=lambda r:r['stiffness_ratio'])
            ax.plot([r['stiffness_ratio'] for r in rr],[r['jitter_local_percent'] if r['jitter_local_percent'] is not None else np.nan for r in rr],marker=marker,c=color,lw=1,ms=4,label=label)
        ax.set(xlabel='Stiffness ratio',ylabel='Jitter local [%]',title=f'{pct:g}% downcross');ax.grid(alpha=.15);ax.legend(title='Peak selector')
    fig.tight_layout();fig.savefig(out/'figures/jitter_comparison.png',dpi=170);core.plt.close(fig)
    windows=[]
    def overlay(ax,case,flow,base,runs,mains,pct,left,right):
        t,x=flow['t'],flow['x'];mask=(t>=left)&(t<right);ax.plot(t[mask],x[mask]*1e6,'k-',lw=.8)
        for reg in base['regions']:
            if reg['signal']=='Ug' and left<=reg['start_s']<right:ax.axvline(reg['start_s'],c='#bbbbbb',ls='--',lw=.6)
        for label,(color,marker) in styles.items():
            pp=[p for p in mains[label].values() if left<=p['event_time_s']<right]
            ax.scatter([p['event_time_s'] for p in pp],[p['sample_value']*1e6 for p in pp],s={'A':28,'B':58,'C':90}[label],
                marker=marker,facecolors=color if label=='A' else 'none',edgecolors=color,label=label+' main',zorder=4)
            es=[e for e in runs[label,pct]['series']['main_'+core.flow_threshold_suffix(pct)] if left<=e['event_time_s']<right]
            ax.scatter([e['event_time_s'] for e in es],[e['threshold_m3_s']*1e6 for e in es],s=30,marker={'A':'x','B':'+','C':'o'}[label],
                c=color if label!='C' else None,edgecolors=color if label=='C' else None,facecolors='none' if label=='C' else None,label=label+f' {pct:g}%',zorder=5)
        ax.set(xlim=(left,right),xlabel='Time [s]',ylabel='Ug [mL/s]');ax.ticklabel_format(axis='x',useOffset=False);ax.legend(fontsize=7,ncol=3)
    if figure_data:
        fig,axes=core.plt.subplots(len(figure_data),len(percentages),figsize=(8*len(percentages),3.7*len(figure_data)),squeeze=False)
        figb,axb=core.plt.subplots(len(figure_data),1,figsize=(15,3.5*len(figure_data)),squeeze=False)
        for i,(case_id,(case,flow,base,runs,mains)) in enumerate(figure_data.items()):
            impacts=[r for r in tables['single_region_impacts'] if r['case_id']==case_id and r['single_region_delta_pp'] is not None]
            strongest=max(impacts,key=lambda r:abs(r['single_region_delta_pp'])) if impacts else None
            reg=next(r for r in base['regions'] if r['signal']=='Ug' and r['region_id']==strongest['region_id']) if strongest else next(r for r in base['regions'] if r['signal']=='Ug' and r['start_s']>=.74)
            scale=reg['search_scale_s'];left=max(case['start_s'],reg['start_s']-2*scale);right=min(case['end_s'],reg['end_s']+2*scale)
            for ax,pct in zip(axes[i],percentages):
                overlay(ax,case,flow,base,runs,mains,pct,left,right);ax.set_title(f"r={case['stiffness_ratio']:g}, {pct:g}%: A/B/C peaks and crossings",fontsize=10)
                windows.append(dict(case_id=case_id,figure='representative_peak_changes.png',flow_threshold_pct=pct,start_s=left,end_s=right,region_id=reg['region_id'],selection='largest absolute single-region effect'))
            flags=[r for r in tables['boundary_suspicions'] if r['case_id']==case_id and r['selector']=='A']
            boundary=next((r for r in flags if r['kind']=='SMALL_MAIN_NEAR_LARGER_PEAK_ACROSS_REGION_BOUNDARY'),flags[0] if flags else None)
            middle=(boundary['time_s']+boundary['other_time_s'])/2 if boundary else (left+right)/2
            l=max(case['start_s'],middle-2*scale);r=min(case['end_s'],middle+2*scale)
            overlay(axb[i,0],case,flow,base,runs,mains,percentages[0],l,r)
            axb[i,0].set_title(f"r={case['stiffness_ratio']:g}: ACF-boundary diagnostic; regions unchanged",fontsize=10)
            windows.append(dict(case_id=case_id,figure='boundary_examples.png',flow_threshold_pct=percentages[0],start_s=l,end_s=r,selection=boundary['kind'] if boundary else 'no warning; representative context'))
        fig.tight_layout();fig.savefig(out/'figures/representative_peak_changes.png',dpi=170);core.plt.close(fig)
        figb.tight_layout();figb.savefig(out/'figures/boundary_examples.png',dpi=170);core.plt.close(figb)
    atomic_csv(out/'figure_windows.csv',windows,('case_id','figure','start_s','end_s'))


def selection_comparison_report(tables,percentages,n_cases):
    rows=tables['comparison'];changes=tables['main_peak_changes'];solo=tables['single_region_impacts']
    def f(v):return 'UNDEFINED' if v is None else f'{v:.6f}'
    def lookup(ratio,pct,label):return next((r for r in rows if r['stiffness_ratio']==ratio and r['flow_threshold_pct']==pct and r['selector']==label),None)
    text=['# 主ピーク選択基準によるJitter感度検証','',
        f'{n_cases}条件、900 Pa、[0.50,1.00)秒。A=現行prominence最大、B=10%判定通過候補内で流量最大、C=全raw極大の流量最大（補助）。B/Cの流量同値は早い補間時刻、次に候補IDで固定。Aの同値規則は従来どおり。','',
        'ACF領域・10%明瞭さ・MAD床は不変。Cだけは明瞭さを通らない実raw極大も候補に含む。選択後は各方式の主ピークから次主ピークまでを再探索し、同じ線形補間・最初の0まで/0なしは最初の実測全体最小点までの最後の下降交差を使用。全隣接間隔と全Jitter/RAP/PPQ5項を使用し、欠落補完・分割・長短間隔の除外なし。','',
        f'Aの保存済み結果との回帰：{len(tables["regression_checks"])}組の全イベントID・親ID・時刻・全間隔・J/RAP/PPQ5が完全一致（時刻差0）。','',
        '## 全条件比較','',
        '| 剛性比 | 閾値 | A Jitter [%] | B Jitter [%] | B−A [pp] | C Jitter [%] | C−A [pp] | A/B/Cイベント数 |',
        '|---:|---:|---:|---:|---:|---:|---:|---|']
    for ratio in sorted({r['stiffness_ratio'] for r in rows}):
        for pct in percentages:
            a,b,c=[lookup(ratio,pct,k) for k in ('A','B','C')]
            text.append(f"| {ratio:.2f} | {pct:g}% | {f(a['jitter_local_percent'])} | {f(b['jitter_local_percent'])} | {f(b['jitter_local_percent_delta_pp'])} | {f(c['jitter_local_percent'])} | {f(c['jitter_local_percent_delta_pp'])} | {a['n_events']}/{b['n_events']}/{c['n_events']} |")
    text+=['','## 主ピーク変更とイベント変化','']
    for label in ('B','C'):
        unique=[r for r in rows if r['selector']==label and r['flow_threshold_pct']==percentages[0]]
        text.append(f"{label}：既存主ピークの置換{sum(r['n_changed_existing_main_peaks'] for r in unique)}領域、主なし領域への実raw主ピーク追加{sum(r['n_added_main_peaks'] for r in unique)}領域。")
        for pct in percentages:
            rr=[r for r in rows if r['selector']==label and r['flow_threshold_pct']==pct and r['jitter_local_percent_delta_pp'] is not None]
            if rr:
                worst=max(rr,key=lambda r:abs(r['jitter_local_percent_delta_pp']))
                text.append(f"- {pct:g}%：最大|ΔJitter|はr={worst['stiffness_ratio']:g}の{f(worst['jitter_local_percent_delta_pp'])} pp。追加{sum(r['n_event_added'] for r in rr)}・欠落{sum(r['n_event_missing'] for r in rr)}・同領域内時刻変更{sum(r['n_event_time_changed'] for r in rr)}件。")
    text+=['','[comparison.csv](comparison.csv)にRAP/PPQ5・平均/標準偏差/最小/最大間隔・全差項数・主ピーク変更数を保存。[events_diff.csv](events_diff.csv)はまず時刻差2 ps以内を全列で1対1対応し、親/領域の移動を区別。残りは同領域を時刻変更として対応し、未対応を追加/欠落とした。時刻を動かしたり評価対象を間引く処理ではない。','',
        '[main_peak_changes.csv](main_peak_changes.csv)は主ピークID・時刻・流量・prominence・明瞭さ・拒否理由・各閾値イベントを記録。[region_effects.csv](region_effects.csv)には主ピーク自身が同じでも次主ピーク/下降枝の境界が変わった領域を含む。','',
        '## 0.55→0.60の差','',
        '| 閾値 | 方式 | J(0.55) [%] | J(0.60) [%] | J(0.55)−J(0.60) [pp] | 大小関係維持 |','|---:|---|---:|---:|---:|---|']
    for pct in percentages:
        for label in ('A','B','C'):
            a,b=lookup(.55,pct,label),lookup(.60,pct,label)
            if a and b and a['jitter_local_percent'] is not None and b['jitter_local_percent'] is not None:
                text.append(f"| {pct:g}% | {label} | {f(a['jitter_local_percent'])} | {f(b['jitter_local_percent'])} | {f(a['jitter_local_percent']-b['jitter_local_percent'])} | {'0.55 > 0.60' if a['jitter_local_percent']>b['jitter_local_percent'] else '反転/同値'} |")
    text+=['','## 個別領域の影響','',
        '下表は1領域だけをB/Cの観測ピークへ置換し、前後の下降枝も含めて全イベント列を再計算した診断。全領域を同時変更したΔとは別で、効果は非線形・非加算。物理的因果の寄与量とは扱わない。','',
        '| 条件 | 閾値 | 置換先 | 領域 | 単独置換ΔJitter [pp] |','|---|---:|---|---:|---:|']
    for r in sorted((r for r in solo if r['single_region_delta_pp'] is not None),key=lambda r:abs(r['single_region_delta_pp']),reverse=True)[:12]:
        text.append(f"| {r['case_id']} | {r['flow_threshold_pct']:g}% | {r['selector']} | {r['region_id']} | {f(r['single_region_delta_pp'])} |")
    text+=['','全件は[single_region_impacts.csv](single_region_impacts.csv)。raw候補の拒否理由・局所振幅・prominence等は[raw_peak_candidates.csv](raw_peak_candidates.csv)、採用されたC主ピークの明瞭さは[peaks.csv](peaks.csv)。','',
        '## ACF領域の妥当性と解釈','',
        'ACF分割は変更していない。隣領域の大きなraw峰が0.55T0以内、主ピーク間隔が0.55T0未満、または隣接主ピーク間の谷が両ピークの低い方の50%より高い箇所を、目視用の疑いとして記録した。診断フラグは選択やJitterの採否には一切使わない。同じ物理パルスという断定ではない。','',
        '[boundary_suspicions.csv](boundary_suspicions.csv)に時刻・候補ID・領域を保存。領域内で最大流量を選べても、目視上の1パルス1主ピーク対応は保証されない。Cがraw候補で空領域に主を追加しても、存在しないイベントの補間ではない。','',
        '![Jitter比較](figures/jitter_comparison.png)','',
        '[0.45/0.55/0.60の主ピーク・交差比較](figures/representative_peak_changes.png) / [領域境界の例](figures/boundary_examples.png) / [図の時間窓](figure_windows.csv)。','',
        'Aは主解析として維持し、Bを主要感度、Cを補助感度として併記する。値が小さくなることだけで選択法を変更しない。結果が選択依存する条件、領域境界の疑いを人手で検証してから恒久変更を判断する。振動安定化の物理的原因は今回評価していない。CPP・ソルバ・元波形・過去結果は非改変。','']
    return '\n'.join(text)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--manifest', type=Path, help='CSV list, relative to repository root; NOT the single-run dt manifest')
    source.add_argument('--paths', nargs='+', help='explicit airflow_vt.dat files or run directories; metadata remains NA')
    p.add_argument('--dt', type=float, default=1e-5, help='solver seconds per step; CSV dt_s overrides')
    p.add_argument('--start', type=float, default=.5, help='inclusive default start; CSV start_s overrides')
    p.add_argument('--end', type=float, default=1., help='exclusive default end; CSV end_s overrides')
    p.add_argument('--flow-threshold-pct', type=core.flow_threshold_percent, default=3., help='flow downcross level as percent of Upeak (0 < pct < 100); peak clarity remains 10%%')
    p.add_argument('--peak-selector',choices=core.PEAK_SELECTORS,default='prominence',help='A default: prominence; B diagnostic: clear-flow; C auxiliary diagnostic: raw-flow')
    p.add_argument('--compare-peak-selectors',action='store_true',help='separate A/B/C sensitivity study; original batch outputs are protected')
    p.add_argument('--comparison-thresholds',nargs='+',type=core.flow_threshold_percent,default=[3.,5.],help='percentages for --compare-peak-selectors; 10 is also supported')
    p.add_argument('--baseline-root',type=Path,default=Path('analysis_runs/jitter_batch_20261108'),help='parent of saved allRatio_001 / allRatio_<pct>pct_001 baseline batches')
    p.add_argument('--out-dir', type=Path, help='NEW directory; relative to repository root')
    p.add_argument('--dry-run', action='store_true', help='preflight all inputs, print list, write nothing and run no detector')
    policy = p.add_mutually_exclusive_group()
    policy.add_argument('--continue-on-error', dest='continue_on_error', action='store_true', default=True, help='continue after case failures (default)')
    policy.add_argument('--fail-fast', dest='continue_on_error', action='store_false', default=argparse.SUPPRESS, help='stop after first case failure, mark later cases NOT_RUN')
    return p


def run(args):
    cases = load_cases(args)
    preflight(cases)
    if args.dry_run:
        bad = sum(c['error'] is not None for c in cases)
        print(f'DRY RUN: {len(cases)-bad} ready, {bad} failed; no output written')
        return 1 if bad == len(cases) else 2 if bad else 0
    stamp = datetime.now(ZoneInfo('Asia/Tokyo')).strftime('%Y%m%dT%H%M%S%f')
    out = repo_path(args.out_dir) if args.out_dir else REPO_ROOT / 'analysis_runs/jitter_batch_20261108' / ('batch_' + stamp)
    protected = REPO_ROOT / 'analysis_runs/manual_run_archive'
    if out.is_relative_to(protected.resolve()) or any(p.name == 'manual_run_archive' for p in out.parents):
        raise ValueError('manual_run_archive is read-only')
    for case in cases:
        if case.get('resolved_flow') and (out.is_relative_to(case['resolved_flow'].parent) or case['resolved_flow'].is_relative_to(out)):
            raise ValueError('output must not contain or reside in source dat directory')
    if out.exists():
        raise ValueError('requires a NEW out-dir; existing results are never overwritten: ' + str(out))
    if args.compare_peak_selectors:
        return run_peak_selection_comparison(args,cases,out)
    config = json.dumps(core.GLOBAL_CYCLE_CONFIG, sort_keys=True)
    pct = threshold_pct(cases[0])
    suffix = core.flow_threshold_suffix(pct)
    settings_id = hashlib.sha256((core.sha256(core.__file__) + config + threshold_rule(pct)+args.peak_selector).encode()).hexdigest()
    try:
        commit = subprocess.check_output(['git', '-C', str(REPO_ROOT), 'rev-parse', 'HEAD'], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = 'UNAVAILABLE'
    identifiers = dict(source_commit=commit, core_sha256=core.sha256(core.__file__), batch_sha256=core.sha256(__file__),
                       settings_id=settings_id, fixed_settings=config)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'cases').mkdir()
    summaries = [initial_summary(c, settings_id) for c in cases]
    failures = []
    aggregate = {name: [] for name in ('events_'+suffix, 'periods_'+suffix, 'metric_terms_'+suffix, 'area_alignment', 'peaks')}
    save_state(out, cases, summaries, failures, aggregate, identifiers, 'IN_PROGRESS')
    halted = False
    interrupted = False
    for index, case in enumerate(cases):
        row = summaries[index]
        try:
            if halted:
                row.update(status='NOT_RUN_INTERRUPTED' if interrupted else 'NOT_RUN_FAIL_FAST', warnings='earlier case halted batch')
                failures.append(dict(case_id=case['case_id'], status=row['status'], path=str(case.get('resolved_flow', '')), phase='NOT_RUN', reason=row['warnings']))
            elif case['error']:
                raise ValueError(case['error'])
            else:
                print(f"Analyzing {index+1}/{len(cases)}: {case['case_id']}", flush=True)
                result = analyze_case(case, settings_id)
                with tempfile.TemporaryDirectory(prefix='.' + case['case_id'] + '_', dir=out / 'cases') as directory:
                    temporary = Path(directory)
                    for name, rows in result.items():
                        required = threshold_columns(SUMMARY_COLUMNS, pct) if name == 'summary' else threshold_columns(AREA_COLUMNS, pct) if name == 'area_alignment' else EVENT_COLUMNS if name == 'events_'+suffix else PERIOD_COLUMNS if name == 'periods_'+suffix else ('case_id',)
                        core.write_csv(temporary / (name + '.csv'), csv_records(rows), csv_columns(rows, required))
                    core.write_csv(temporary / 'inputs_resolved.csv', csv_records([resolved_row(case, identifiers)]))
                    temporary.rename(out / 'cases' / case['case_id'])
                summaries[index] = result['summary'][0]
                for name in aggregate:
                    aggregate[name].extend(result[name])
        except KeyboardInterrupt:
            row.update(status='INTERRUPTED', warnings='KeyboardInterrupt; incomplete case is not successful')
            failures.append(dict(case_id=case['case_id'], status='INTERRUPTED', path=str(case.get('resolved_flow', '')), phase='ANALYSIS', reason='KeyboardInterrupt'))
            interrupted = halted = True
        except Exception as exc:
            row.update(status='FAILED', warnings=str(exc), area_status=case['area_status'], area_reason=case['area_reason'])
            failures.append(dict(case_id=case['case_id'], status='FAILED', path=str(case.get('resolved_flow', case.get('flow_path') or case.get('run_path'))), phase='PREFLIGHT' if case['error'] else 'ANALYSIS', reason=str(exc)))
            halted = not args.continue_on_error
        save_state(out, cases, summaries, failures, aggregate, identifiers, 'IN_PROGRESS')
    state = 'INTERRUPTED' if interrupted else 'COMPLETE_WITH_FAILURES' if failures else 'COMPLETE'
    save_state(out, cases, summaries, failures, aggregate, identifiers, state)
    comparison_plot(out, summaries)
    print('\ncase_id\tstatus\tevents\tJitter[%]\tRAP[%]\tPPQ5[%]')
    for r in summaries:
        values = [f'{r[k]:.6f}' if r.get(k) is not None else 'NA' for k in ('jitter_local_pct', 'rap_pct', 'ppq5_pct')]
        print(f"{r['case_id']}\t{r['status']}\t{r.get('n_events')}\t" + '\t'.join(values))
    good = sum(r['status'] in ('COMPLETED', 'INSUFFICIENT_EVENTS') for r in summaries)
    print(f'Saved: {out}\nProcessed: {good}; failed/not run: {len(cases)-good}')
    return 130 if interrupted else 1 if good == 0 else 2 if good < len(cases) else 0


def main(argv=None):
    p = parser()
    try:
        return run(p.parse_args(argv))
    except (ValueError, OSError, RuntimeError, OverflowError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
