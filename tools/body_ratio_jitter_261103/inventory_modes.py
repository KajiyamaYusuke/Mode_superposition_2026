#!/usr/bin/env python3
"""Conservative inventory only; never certify an input or start a solver run."""
import argparse
import csv
from decimal import Decimal
import hashlib
import math
from pathlib import Path
import re

import numpy as np


TOKEN = re.compile(r"b(\d+)(?:p(\d+))?c(\d+)\Z")
ARRAY = re.compile(rb'<DataArray\b[^>]*Name="Displacement_field,_([XYZ])-component_@_([^\"]+)"[^>]*>(.*?)</DataArray>', re.S)


def parse_token(token):
    match = TOKEN.fullmatch(token)
    if not match:
        raise ValueError(f"invalid body token: {token}")
    integer, fraction, cover = match.groups()
    return Decimal(integer + ('.' + fraction if fraction else '')), int(cover)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--left-mode', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--confirmed-export-order', action='store_true',
                        help='Record the user-confirmed COMSOL export order for b1p8c3')
    parser.add_argument('--confirmed-filename-materials', action='store_true',
                        help='Record the user-confirmed material values in filenames')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    left = args.left_mode.read_bytes()
    left_points = re.search(rb'<Points>.*?</Points>', left, re.S)
    if left_points is None:
        raise ValueError('left VTU has no Points')
    geometry_hash = sha(left_points.group())
    pairs = {}
    for path in args.source_dir.iterdir():
        if path.suffix.lower() not in ('.txt', '.vtu', '.vtk'):
            continue
        try:
            body, cover = parse_token(path.stem)
        except ValueError:
            continue
        pairs.setdefault((body, cover), {}).setdefault(path.suffix.lower(), []).append(path)
    fields = ['ratio','body_left_kpa','body_right_kpa','body_token','cover_kpa','frequency_txt','mode_vtk','raw_format','npoints','nmodes','frequency_1_hz','frequency_20_hz','raw_hashes','prepared_hashes','geometry_check','mode_order_check','status','notes']
    rows = []
    for i in range(2, 20):
        ratio = Decimal(i) / 20
        body = ratio * 12
        group = pairs.get((body, 3), {})
        txt, vtu = group.get('.txt', []), group.get('.vtu', [])
        row = dict.fromkeys(fields, '')
        row.update(ratio=str(ratio), body_left_kpa='12', body_right_kpa=str(body), cover_kpa='3', status='INPUT_BLOCKED')
        notes = []
        if len(txt) != 1 or len(vtu) != 1 or group.get('.vtk'):
            row['notes'] = 'missing or duplicate frequency/mode files'
            rows.append(row)
            continue
        row.update(body_token=txt[0].stem, frequency_txt=str(txt[0].resolve()), mode_vtk=str(vtu[0].resolve()))
        tb, vb = txt[0].read_bytes(), vtu[0].read_bytes()
        row['raw_hashes'] = f'txt:{sha(tb)};vtu:{sha(vb)}'
        row['raw_format'] = 'XML VTU' if b'<VTKFile' in vb[:500] else 'UNKNOWN'
        point_count = re.search(rb'NumberOfPoints="(\d+)"', vb)
        points = re.search(rb'<Points>.*?</Points>', vb, re.S)
        npoints = int(point_count.group(1)) if point_count else 0
        row['npoints'] = npoints
        row['geometry_check'] = 'POINTS_IDENTICAL' if points and sha(points.group()) == geometry_hash else 'MISMATCH_OR_MISSING'
        lines = tb.decode('utf-8-sig').splitlines()
        header = next((line for line in lines if line.startswith('% Model:')), '')
        notes.append('COMSOL ' + header)
        if args.confirmed_filename_materials:
            notes.append('material parameters confirmed by user from filename')
        elif row['body_token'] not in header:
            notes.append('model header does not match filename; material identity unverified')
        frequencies = []
        for line in lines:
            if line.strip() and not line.lstrip().startswith('%'):
                values = [float(x) for x in line.split()]
                if len(values) < 3 or not all(math.isfinite(x) for x in values[:3]):
                    raise ValueError(f'invalid frequency row in {txt[0]}')
                frequencies.append(values)
        if frequencies:
            row['frequency_1_hz'] = frequencies[0][0]
        if len(frequencies) >= 20:
            row['frequency_20_hz'] = frequencies[19][0]
        arrays = ARRAY.findall(vb)
        axes = {'X': [], 'Y': [], 'Z': []}
        valid = True
        for axis, name, data in arrays:
            values = np.fromstring(data.decode('ascii'), sep=' ')
            if len(values) != npoints or not np.isfinite(values).all():
                valid = False
            axes[axis.decode()].append(name.decode())
        row['nmodes'] = len(axes['X'])
        names = axes['X']
        if len(names) < 20 or axes['Y'] != names or axes['Z'] != names or not valid:
            notes.append('missing, nonfinite or mismatched mode arrays')
        hz = [re.match(r'([0-9.]+)_Hz', name) for name in names]
        if any(m is None for m in hz):
            row['mode_order_check'] = 'INVALID_NAMES'
        else:
            rounded = [Decimal(m.group(1)) for m in hz]
            mismatch = (len(frequencies) < 20 or len(rounded) < 20 or
                        any(abs(Decimal(str(frequencies[j][0])) - rounded[j]) > Decimal(5).scaleb(rounded[j].as_tuple().exponent - 1) for j in range(min(20, len(rounded), len(frequencies)))))
            if mismatch:
                row['mode_order_check'] = 'MISMATCH'
            elif len(set(rounded)) != len(rounded):
                if args.confirmed_export_order and row['body_token'] == 'b1p8c3' and names[18:20] == ['109.17_Hz_(1)', '109.17_Hz_(2)'] and len(set(rounded[:18])) == 18 and rounded[18] not in rounded[:18]:
                    row['mode_order_check'] = 'USER_CONFIRMED_EXPORT_ORDER_19_20'
                    notes.append('COMSOL export follows eigenmode order; user confirmed (1)=19, (2)=20')
                else:
                    row['mode_order_check'] = 'MODE_ORDER_AMBIGUOUS'
                    notes.append('repeated rounded frequency; exact mode correspondence not established')
            else:
                row['mode_order_check'] = 'ROUNDED_FREQUENCIES_MATCH; identity unverified'
        notes.append('full provenance, geometry/NAS and solver-read gates not certified')
        row['notes'] = '; '.join(notes)
        rows.append(row)
        print(row['body_token'], row['geometry_check'], row['mode_order_check'], flush=True)
    with (args.output_dir / 'mode_inventory.csv').open('x', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (args.output_dir / 'cases_manifest.csv').open('x', newline='') as out:
        names = ['ratio','body_right_kpa','pressure_pa','case_id','param_path','run_dir','left_inputs_hash','right_inputs_hash','mesh_hash','solver_hash','physics_config_hash','nstep','dt_s','analysis_window','status','reason']
        writer = csv.DictWriter(out, fieldnames=names)
        writer.writeheader()
        for row in rows:
            for pressure in (900, 1400):
                writer.writerow(dict(ratio=row['ratio'], body_right_kpa=row['body_right_kpa'], pressure_pa=pressure, case_id='r'+f"{Decimal(row['ratio']):.2f}".replace('.', 'p')+f'_P{pressure}', right_inputs_hash=row['raw_hashes'], nstep=100000, dt_s='1e-5', analysis_window='[0.50,1.00)', status='INPUT_BLOCKED', reason='preflight not passed; '+row['notes']))


if __name__ == '__main__':
    main()
