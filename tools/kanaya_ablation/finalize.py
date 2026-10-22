#!/usr/bin/env python3
"""Final regression, manifests and preservation audit (no physical solver changes)."""
import argparse
import json
import subprocess
from pathlib import Path
from run_matrix import ROOT, OUT, dump, sha, git, protection
from validate import compare_case


def regression():
    reference=OUT/'validation/pristine_modern_full/output'
    modern=OUT/'R0_modern/output'
    result=compare_case(reference,modern,skip_vtu=True)
    assert result['passed'],result
    result.update(pristine_revision='bc62cf8f3a2c01fb1ed0a33f6d4a6864f6c1d2d5',
                  nmode=20,nstep=20000,VTU_note='VTU coordinates checked bitwise in 80-step closed/open regressions; primary VTU intentionally suppressed')
    dump(OUT/'validation/full_modern_regression.json',result)
    print('PASS full-length pristine modern / R0 bitwise invariance',result)


def finish():
    cases=['R0_modern','R1_legacy_fluid_all','R3_modern_coefficients',
           'R6_modern_surface_mapping_tie_safe','S1_legacy_zeta010','S2_legacy_vt0175']
    validation=json.loads((OUT/'comparison_validation.json').read_text())
    assert all(validation[name]['passed'] for name in cases)
    assert json.loads((OUT/'validation/full_modern_regression.json').read_text())['passed']
    before=json.loads((OUT/'pre_task_protection.json').read_text())
    after=protection();assert before==after,'protected original changed'
    dump(OUT/'final_preservation_audit.json',{'passed':True,'snapshot':after})
    records=[]
    keys=['legacy_area','legacy_flow_coefficients','legacy_wall_viscosity','legacy_separation',
          'legacy_post_separation','legacy_surface_pressure_load','legacy_contact','legacy_integrator']
    for name in cases:
        path=OUT/name
        manifest=json.loads((path/'manifest.json').read_text())
        assert manifest['exit_code']==0
        effective=manifest['effective_parameters']
        record={k:manifest[k] for k in ['name','binary_sha256','parameters_sha256','elapsed_s','threads','affinity']}
        record['flags']={k:effective[k] for k in keys}
        record['zeta']=effective['zetaL'];record['vocal_tract_length_m']=effective['vocal_tract_length_m']
        records.append(record)
        dump(path/'analysis_manifest.json',dict(passed=True,source_case_manifest_sha256=sha(path/'manifest.json'),
            compatibility_contract={
                'legacy_surface_pressure_load':'requires legacy area; material-row direct pressure and historical load crop',
                'legacy_area':'matching topology and reflected x,z,y checked each step; 42 actual material rows, not 50 resamples',
                'hybrid_physical_mapping':'requires monotonic row x; abrupt pressure endpoint resolves only 32-epsilon coordinate equality upstream',
                'legacy_contact_integrator':'unsupported true rejected; current contact and Newmark retained'},
            validation=validation[name],work_interval='Q_n * (q_(n+1)-q_n)',excluded_final_step=19999,
            equal_windows={'early_steps':[10000,14998],'late_steps':[15000,19998]},
            output_VTU_disabled='waveform/force/modal full-data outputs retained; short invariance tested'))
    failed=json.loads((OUT/'R6_modern_surface_mapping/manifest.json').read_text())
    dump(OUT/'run_matrix.json',dict(completed_physical_cases=records,
        failed_excluded_attempt={'case':failed['name'],'exit_code':failed['exit_code'],
            'first_Qfluid_asymmetry_step':1028,'stopped_step':1030,
            'reason':'one-ulp coordinate ambiguity at discontinuous legacy wall-pressure boundary',
            'replacement':'R6_modern_surface_mapping_tie_safe'},
        decisions=[
            'Both primary cases sustain open oscillation but no evaluation-window closure: do not claim no motion or successful periodic closure.',
            'Select only R3 coefficients and R6 mapping to discriminate mechanism sensitivity; do not run full R2–R6 or all combinations.',
            'Run two explicitly separate historical-condition controls: zeta .10 and vocal tract .175 m; never combine these changes.',
            'No contact/integrator modification or large pressure/mode sweep.'],
        omitted_cases={'R2':'area/mapping dependency; full matrix not warranted by primary lack of closure',
                       'R4':'wall/separation combination not independently isolated in this finite experiment',
                       'R5':'historical load crop strongly limits downstream-pressure leverage; not executed',
                       'combinations':'no factorial/interaction expansion executed'},
        measured_schedule=json.loads((OUT/'benchmark_final_selection.json').read_text())['winner']))
    sources=[p for directory in ['src','include','tests','tools/kanaya_ablation']
             for p in (ROOT/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    sources.append(ROOT/'CMakeLists.txt')
    dump(OUT/'provenance.json',dict(worktree=str(ROOT),git=git(ROOT),
        base_revision='bc62cf8f3a2c01fb1ed0a33f6d4a6864f6c1d2d5',
        user_reference=json.loads((OUT/'references/user_reference.json').read_text()),
        plan_sha256=sha(ROOT/'docs/plans/KANAYA_LEGACY_SWITCH_ABLATION_20MODES_1600PA.md'),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},
        build_flags=(ROOT/'build/switchable/CMakeFiles/simulation.dir/flags.make').read_text(),
        compiler='/opt/intel/oneapi/compiler/2025.1/bin/icpx',
        limitations=['Historical successful output and binary unavailable; equation equivalence is not full coupled replay.',
                     'Original historical runtime parameters/input hashes cannot be recovered from current checked-in example parameters.',
                     'Open oscillations over .2 s do not establish an asymptotic attractor or perturbation stability.',
                     'Per-mode damping work is a normalized reference; sampled mass Gram matrix is non-diagonal.']))
    print('PASS preservation and final provenance/analysis manifests')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['regression','finish']);a=p.parse_args()
    {'regression':regression,'finish':finish}[a.action]()
