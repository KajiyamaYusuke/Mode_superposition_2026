#!/usr/bin/env python3
"""Prepare, but never execute, the requested independent 1 s cases."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'analysis_runs/12conditions_20261026_ZivIsl/cases/L12_R6_P900/input/param.txt'
SOURCE = ROOT / 'input/M5_test/261103'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def prepare(root, solver):
    rows = list(csv.DictReader((root / 'mode_inventory.csv').open()))
    if len(rows) != 18:
        raise ValueError('expected exactly 18 inventory rows')
    base = BASE.read_text()
    if base.count('# number of step nstep\n50000') != 1:
        raise ValueError('unexpected reference nstep')
    base = base.replace('# number of step nstep\n50000', '# number of step nstep\n100000')
    base += '\nwriteDeformationVtu = 0\n'
    paths = ('leftFrequencyFile', 'leftModeFile', 'leftSurfaceNasFile', 'rightSurfaceNasFile')
    for key in paths:
        if len(re.findall(r'^' + key + r'\s*=\s*/.+$', base, re.M)) != 1:
            raise ValueError('missing absolute input path: '+key)
    if not re.search(r'^pressureRampTimeSec\s*=\s*0\.10\s*$', base, re.M):
        raise ValueError('unexpected pressure ramp')
    if not re.search(r'# time step for time integration dt\s*\n1\.0E-5', base):
        raise ValueError('unexpected dt')
    plan = []
    for row in rows:
        if row['geometry_check'] != 'POINTS_IDENTICAL' or row['mode_order_check'] not in ('ROUNDED_FREQUENCIES_MATCH; identity unverified', 'USER_CONFIRMED_EXPORT_ORDER_19_20'):
            raise ValueError('unvalidated mode: '+str(row))
        token = row['body_token']
        freq, mode = SOURCE / (token+'.txt'), SOURCE / (token+'.vtu')
        if f'txt:{digest(freq)};vtu:{digest(mode)}' != row['raw_hashes']:
            raise ValueError('input changed: '+token)
        for pressure in (900, 1400):
            caseid = 'r'+f"{float(row['ratio']):.2f}".replace('.', 'p')+f'_P{pressure}'
            params = base
            params = re.sub(r'(?<=# subglottal pressure Ps \(Pa\)\n)900(?=\s*\n)', str(pressure), params)
            for key, val in (('rightFrequencyFile', freq), ('rightModeFile', mode)):
                params, count = re.subn(r'^'+key+r'\s*=.*$', key+' = '+str(val), params, flags=re.M)
                if count != 1:
                    raise ValueError('unexpected param key '+key)
            if params.count('# subglottal pressure Ps (Pa)\n'+str(pressure)) != 1:
                raise ValueError('pressure replacement failed')
            plan.append((row, pressure, caseid, params))
    cases = root / 'cases'
    if cases.exists():
        raise FileExistsError('cases directory already exists; refusing overwrite')
    cases.mkdir()
    for row, pressure, caseid, params in plan:
        folder = cases / caseid / 'input'
        folder.mkdir(parents=True)
        (folder / 'param.txt').write_text(params)
    config = {'execution_authorized': False, 'reason': 'preflight gates remain incomplete', 'solver': str(solver.resolve()), 'solver_sha256': digest(solver), 'reference_param_sha256': digest(BASE), 'writeDeformationVtu': '0 (user requested)', 'surface': 'M5_surface_T3_d2.nas', 'case_count': len(plan), 'dt_s': 1e-5, 'nstep': 100000}
    (root / 'prepared_config.json').write_text(json.dumps(config, indent=2)+'\n')
    with (root/'prepared_cases.csv').open('x', newline='') as f:
        w=csv.DictWriter(f,fieldnames=['case_id','ratio','pressure_pa','parameter_file','parameter_sha256','frequency_sha256','mode_sha256','status']); w.writeheader()
        for row, pressure, caseid, _ in plan:
            path=cases/caseid/'input/param.txt'; freq=SOURCE/(row['body_token']+'.txt'); mode=SOURCE/(row['body_token']+'.vtu')
            w.writerow(dict(case_id=caseid,ratio=row['ratio'],pressure_pa=pressure,parameter_file=str(path),parameter_sha256=digest(path),frequency_sha256=digest(freq),mode_sha256=digest(mode),status='PREFLIGHT_INCOMPLETE'))
    print(f'Prepared {len(plan)} cases; no solver was started.')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sweep-dir',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True)
    a=p.parse_args(); prepare(a.sweep_dir.resolve(),a.solver.resolve())
