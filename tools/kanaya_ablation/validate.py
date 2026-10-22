#!/usr/bin/env python3
"""Fixed-shape user-main equivalence, parser tests and modern regression gates."""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from run_matrix import ROOT, OUT, COMPILER, dump, sha, environment, parameters, run_case


def extract(source, signature):
    # Ignore braces inside comments/strings when locating the body; preserve
    # the original bytes in the exported method.
    clean = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                   lambda m: ' '*len(m[0]), source)
    start = clean.index(signature)
    begin = clean.index('{', start)
    depth = 1; end = begin+1
    while depth:
        if clean[end] == '{': depth += 1
        if clean[end] == '}': depth -= 1
        end += 1
    return source[start:end]


def fixed():
    directory = OUT/'validation/fixed_final'; directory.mkdir(parents=True)
    refs = OUT/'references'
    for name in ['Geometry.h', 'ForceCalculator.h', 'ModeData.h', 'State.h',
                 'Displacement.h', 'SimulationParams.h', 'TimeIntegrator.h']:
        shutil.copyfile(refs/('user_main_'+name), directory/name)
    fc = (refs/'user_main_forceCalculator.cpp').read_text()
    state = (refs/'user_main_state.cpp').read_text()
    signatures = ['ForceCalculator::ForceCalculator(', 'void ForceCalculator::calcForce(',
                  'void ForceCalculator::f2mode(', 'double ForceCalculator::findMinHarea(',
                  'int ForceCalculator::findNsep(']
    bodies = [extract(fc, name) for name in signatures]+[extract(state, 'void State::calcArea(')]
    source = '#include <bits/stdc++.h>\n#include "ForceCalculator.h"\n'+'\n\n'.join(bodies)+'\n'
    (directory/'user_methods.cpp').write_text(source)
    binary = directory/'user_main_frozen'
    command = [str(COMPILER/'bin/icpx'), '-O3', '-std=c++17', '-I'+str(directory),
               '-I'+str(ROOT/'tests'), str(directory/'user_methods.cpp'),
               str(ROOT/'tests/legacy_user_reference_driver.cpp'), str(refs/'user_main_rungeStep.cpp'),
               '-o', str(binary)]
    subprocess.run(command, check=True, env=environment())
    subprocess.run([str(binary), str(directory/'user_snapshots.csv')], check=True, env=environment())
    subprocess.run([str(ROOT/'build/switchable/test_legacy_glottal_methods'), '--snapshots',
                    str(directory/'adapter_snapshots.csv')], check=True, env=environment(), stdout=subprocess.DEVNULL)
    ref = pd.read_csv(directory/'user_snapshots.csv')
    adapter = pd.read_csv(directory/'adapter_snapshots.csv')
    assert list(ref.columns) == list(adapter.columns) and ref.shape == adapter.shape
    comparisons = []
    for name in ref.columns:
        error = np.abs(ref[name]-adapter[name])
        scaled = error/(1e-12+1e-11*np.maximum(ref[name].abs(), adapter[name].abs()))
        assert scaled.max() <= 1, f'user-main frozen equivalence failed {name}: {scaled.max()}'
        comparisons.append({'field': name, 'max_abs': float(error.max()), 'max_scaled': float(scaled.max())})
    # Frozen-geometry, SAME minimum A, SAME U, SAME circuit pressure. The
    # original coefficient statements are exported verbatim to a tiny oracle.
    statements = []
    for name in ['lis', 'lg_m', 'Lg1', 'Rk1', 'Rv1']:
        match = re.search(r'double '+name+r'\s*=([^;]+);', fc)
        assert match is not None
        statements.append(match[0])
    oracle = '#include <bits/stdc++.h>\nstruct G{double xsup;};\nint main(){ G geom{12.}; double rho=1.13,mu=1.83e-5,beta=.3175,lg=17.; std::cout<<std::setprecision(17); for(double min_area_m2:{1e-6,3e-6,1e-5}){\n'+ '\n'.join(statements)+'\nstd::cout<<min_area_m2<<" "<<Lg1<<" "<<Rv1<<" "<<Rk1<<"\\n";}}\n'
    (directory/'coefficient_oracle.cpp').write_text(oracle)
    subprocess.run([str(COMPILER/'bin/icpx'), '-O3', '-std=c++17', str(directory/'coefficient_oracle.cpp'), '-o', str(directory/'coefficient_oracle')], env=environment(), check=True)
    old_coeff = np.fromstring(subprocess.check_output([str(directory/'coefficient_oracle')], env=environment(), text=True), sep=' ').reshape(-1, 4)
    records = []
    for area, lg, rv, rk in old_coeff:
        lg_new = 1.13*.5*.008/area
        rv_new = 12*1.83e-5*.017**2*.008/area**3
        ug, previous, dt = .00025, .0002, 1e-5
        for method, inertance, viscous in [('legacy', lg, rv), ('modern_constant_gap_to_8mm', lg_new, rv_new)]:
            # NR and analytic positive root for the same frozen airway state.
            b=viscous+(inertance+.1+10.)/dt
            c=-(inertance+.1+10.)*previous/dt-1000.
            guess=ug
            for iterations in range(100):
                residual=rk*abs(guess)*guess+b*guess+c
                if abs(residual)<1e-9: break
                guess-=residual/(2*rk*abs(guess)+b)
            analytic=-2*c/(b+np.sqrt(b*b-4*rk*c))
            assert abs(residual)<1e-9 and abs(guess-analytic)<1e-14
            records.append({'Amin_m2': area, 'frozen_Ug_m3_s': ug, 'method': method,
                'Lg': inertance, 'Rv': viscous, 'Rk': rk, 'Lg_old_over_modern': lg/lg_new,
                'Rv_old_over_modern': rv/rv_new, 'NR_residual_at_frozen_U': rk*ug*ug+b*ug+c,
                'NR_solution': guess, 'NR_residual_at_solution': residual, 'analytic_root': analytic})
    pd.DataFrame(records).to_csv(directory/'coefficient_comparison.csv', index=False)
    dump(directory/'validation.json', {'passed': True, 'rows': len(ref), 'columns': comparisons,
        'source_sha256': {str(p.relative_to(refs)): sha(p) for p in refs.glob('user_main_*')},
        'verbatim_method_signatures': signatures+['void State::calcArea('],
        'verbatim_body_sha256': sha(directory/'user_methods.cpp'),
        'source_URL': 'https://github.com/KajiyamaYusuke/Mode_superposition/tree/main',
        'limits': 'Frozen Ug/Pg, geometry/pressure/load/projection/work and Newmark; not the old coupled acoustic boundary/contact algorithm.'})
    print('PASS user-confirmed main frozen area, pressure, full nodal force, modal Q/work, Newmark and SI coefficients')


def compare_case(ref, case, skip_vtu=False):
    required = ['regression_state.bin']
    required += ['dat/'+p.name for p in (ref/'dat').glob('*.dat')]
    required += ['kanaya_diagnostic/'+name for name in ['modal_native.bin', 'modal_energy_summary.csv',
        'flow_geometry.csv', 'force_regions.csv', 'contact_regions.csv', 'consistency.csv']]
    required += ['wav/test_sound.wav']
    if not skip_vtu:
        required += ['result/'+p.name for p in (ref/'result').glob('*.vtu')]
    differences = []
    for name in required:
        if sha(ref/name) != sha(case/name): differences.append(name)
    return {'passed': not differences, 'bitwise_identical_files': len(required)-len(differences),
            'different': differences}


def modern():
    ref = OUT/'validation/untouched_80/output'
    results = {}
    for case in ['modern_implicit_80', 'modern_explicit_80']:
        results[case] = compare_case(ref, OUT/'validation'/case/'output')
        assert results[case]['passed'], results
    dump(OUT/'validation/modern_regression.json', results)
    print('PASS untouched / implicit modern / explicit modern are bitwise identical', results)


def open_modern():
    ref=OUT/'validation/untouched_open_80/output'
    results={}
    for case in ['modern_open_final_80','modern_open_no_vtu_80']:
        results[case]=compare_case(ref,OUT/'validation'/case/'output',skip_vtu='no_vtu' in case)
        assert results[case]['passed'], results
        flow=pd.read_csv(OUT/'validation'/case/'output/kanaya_diagnostic/flow_geometry.csv')
        assert flow.Ug_m3_s.abs().max()>0, 'open-flow regression did not exercise flow'
    dump(OUT/'validation/modern_open_regression.json',results)
    print('PASS open-flow current final binary and VTU suppression are bitwise identical',results)


def parsers():
    directory = OUT/'validation/parser_final'; directory.mkdir(parents=True)
    # Parser-only binary does not initialize geometry or execute a solver.
    code = '#include "SimulationParams.h"\nint main(int c,char**v){SimulationParams p;std::string e;if(c!=2)return 2;if(!p.loadFromFile(v[1],e)||!p.validate(e)){std::cerr<<e;return 1;}return 0;}\n'
    (directory/'parser.cpp').write_text(code)
    binary = directory/'parser'
    subprocess.run([str(COMPILER/'bin/icpx'), '-O3', '-std=c++17', '-I'+str(ROOT/'include/simulation'),
        '-I'+str(ROOT/'include/fluid'), str(directory/'parser.cpp'), str(ROOT/'src/simulation/readParams.cpp'),
        '-o', str(binary)], env=environment(), check=True)
    examples = [('default', '', True), ('modern', 'glottalMethodProfile=modern', True),
        ('all', 'glottalMethodProfile=legacy_fluid_all', True),
        ('coefficient', 'glottalMethodProfile=custom\nlegacyFlowCoefficients=1', True),
        ('bad_flag', 'glottalMethodProfile=custom\nlegacyArea=2', False),
        ('conflict', 'glottalMethodProfile=modern\nlegacyArea=1', False),
        ('bad_profile', 'glottalMethodProfile=unknown', False),
        ('incompatible', 'glottalMethodProfile=custom\nlegacySurfacePressureLoad=1', False),
        ('contact', 'glottalMethodProfile=custom\nlegacyContact=1', False),
        ('integrator', 'glottalMethodProfile=custom\nlegacyIntegrator=1', False),
        ('tract', 'vocalTractLengthM=0.175', True),
        ('bad_tract', 'vocalTractLengthM=-1', False)]
    results = []
    for name, extra, accepted in examples:
        path = directory/(name+'.txt'); path.write_text(parameters(80, extra=extra))
        r = subprocess.run([str(binary), str(path)], env=environment(), capture_output=True, text=True)
        assert (r.returncode == 0) == accepted, (name, r.stderr)
        results.append({'name': name, 'expected_accept': accepted, 'returncode': r.returncode, 'message': r.stderr})
    dump(directory/'validation.json', {'passed': True, 'cases': results})
    print('PASS parser acceptance/rejection, old positional format and profiles')


def mapping():
    directory=OUT/'validation/mapping'; directory.mkdir(parents=True)
    binary=ROOT/'build/switchable/test_legacy_glottal_methods'
    for option,name in [('--mapping','mapping_comparison.csv'),('--snapshots','adapter_snapshots.csv')]:
        subprocess.run([str(binary),option,str(directory/name)],env=environment(),check=True)
    reference=pd.read_csv(OUT/'validation/fixed_final/user_snapshots.csv')
    current=pd.read_csv(directory/'adapter_snapshots.csv')
    delta=np.abs(reference.to_numpy()-current.to_numpy())
    scaled=delta/(1e-12+1e-11*np.maximum(np.abs(reference.to_numpy()),np.abs(current.to_numpy())))
    assert scaled.max()<=1
    frame=pd.read_csv(directory/'mapping_comparison.csv')
    assert len(frame)==36 and np.isfinite(frame.select_dtypes('number')).all().all()
    flat=frame[(frame['case']==0)&(frame.experiment==0)].select_dtypes('number').to_numpy()
    assert np.allclose(flat[0],flat[1],atol=1e-14,rtol=1e-12)
    dump(directory/'validation.json',{'passed':True,'rows':len(frame),
        'original_method_equivalence_max_scaled':float(scaled.max()),
        'experiments':{'0':'same uniform pressure / common nodal support; area-vector formula only',
                       '1':'same uniform pressure / native load support',
                       '2':'native pressure rules / SAME material rows and frozen Ug/Pg; not 50-plane modern'},
        'flat_same_pressure_load_bitwise_or_tolerance_equal':True,
        'linear_physical_x_pressure_interpolation_asserted':True})
    print('PASS separated same-pressure mapping / native-pressure comparison and original-main equivalence')


def boundary():
    directory=OUT/'validation/boundary_fix';directory.mkdir(parents=True)
    subprocess.run([str(ROOT/'build/switchable/test_legacy_glottal_methods'),'--snapshots',
                    str(directory/'adapter_snapshots.csv')],env=environment(),check=True)
    ref=pd.read_csv(OUT/'validation/fixed_final/user_snapshots.csv').to_numpy()
    current=pd.read_csv(directory/'adapter_snapshots.csv').to_numpy()
    scaled=np.abs(ref-current)/(1e-12+1e-11*np.maximum(np.abs(ref),np.abs(current)))
    assert scaled.max()<=1
    binary=ROOT/'build/switchable/test_kanaya_diagnostics'
    with ThreadPoolExecutor(2) as pool:
        a=pool.submit(run_case,'validation/modern_boundary_open_80',binary,steps=80,threads=1,
                      cpus='0',extra='initialGapMm = 0.1\nglottalMethodProfile = modern\n')
        b=pool.submit(run_case,'validation/legacy_boundary_2000',binary,steps=2000,threads=2,
                      cpus='2,4',extra='writeDeformationVtu = 0\nglottalMethodProfile = legacy_fluid_all\n')
        a.result();b.result()
    results={'modern':compare_case(OUT/'validation/untouched_open_80/output',
                                   OUT/'validation/modern_boundary_open_80/output'),
             'legacy':compare_case(OUT/'benchmark/legacy_t4/output',
                                   OUT/'validation/legacy_boundary_2000/output',skip_vtu=True)}
    assert all(v['passed'] for v in results.values())
    assert sha(OUT/'benchmark/legacy_t4/output/csv/glottal_methods.csv')==sha(
        OUT/'validation/legacy_boundary_2000/output/csv/glottal_methods.csv')
    dump(directory/'validation.json',dict(passed=True,regressions=results,
        original_frozen_equivalence_max_scaled=float(scaled.max()),
        fix='machine-equal x at abrupt legacy zero boundary selects upstream endpoint; 32 epsilon relative coordinate tie only',
        no_finite_length_blend=True,modern_and_all_legacy_unchanged=True))
    print('PASS boundary tie regression, original frozen equivalence, modern open and legacy coupled bitwise invariance')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('action', choices=['fixed', 'modern', 'open-modern', 'parser','mapping','boundary'])
    args = parser.parse_args()
    {'fixed': fixed, 'modern': modern, 'open-modern':open_modern, 'parser': parsers,
     'mapping':mapping,'boundary':boundary}[args.action]()
