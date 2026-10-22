#!/usr/bin/env python3
"""Read-only preservation audit and append-only persistent archive."""
import argparse
import json
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from next_common import ROOT,OUT,PREVIOUS,ARCHIVE,BUILD,sha,dump,git,environment
from run_matrix import protection
from validate import extract
from compare_runs import numeric_gate


def audit():
    first=json.loads((OUT/'before.json').read_text())
    assert first['originals']==protection(),'protected original worktree changed'
    allowed={'CMakeLists.txt','include/fluid/FlowModel.h','include/forces/ForceCalculator.h','include/simulation/SimulationParams.h',
             'src/fluid/flowModel.cpp','src/simulation/readParams.cpp','src/simulation/simulation.cpp',
             'tools/kanaya_ablation/compare_runs.py'}
    inspected=0
    for entry in first['files']:
        p=ROOT/entry['path'];assert p.exists(),f'previous file deleted: {p}'
        if entry['path'] in allowed:continue
        stat=p.stat()
        assert (stat.st_size,stat.st_mtime_ns)==(entry['bytes'],entry['mtime_ns']),f'previous file changed: {p}'
        inspected+=1
    before_source=(ARCHIVE/'previous/src/fluid/flowModel.cpp').read_text()
    now_source=(ROOT/'src/fluid/flowModel.cpp').read_text()
    assert extract(before_source,'void FlowModel::calcFlowStep(')==extract(now_source,'void FlowModel::calcFlowStep(')
    production_manifest=json.loads((OUT/'S3_legacy_upstream_circuit/manifest.json').read_text())
    assert sha(BUILD/'test_kanaya_diagnostics')==production_manifest['binary_sha256']
    postrun_changes=[p for p,digest in production_manifest['source_sha256'].items() if sha(ROOT/p)!=digest]
    assert set(postrun_changes)<=set(['include/forces/ForceCalculator.h'])
    paths=['src/fluid/legacyChannelArea.cpp','src/fluid/legacyFlowCoefficients.cpp',
           'src/fluid/legacyWallPressure.cpp','src/forces/legacySurfacePressureLoad.cpp',
           'src/structure/state.cpp','src/structure/rungeStep.cpp','src/structure/modedata.cpp',
           'src/forces/forceCalculatorContact.cpp','src/forces/modalProjector.cpp','src/fluid/channelSections.cpp']
    assert all(sha(ROOT/p)==sha(ARCHIVE/'previous'/p) for p in paths)
    current=git(ROOT)
    assert current['head']==first['git']['head'] and current['branch']==first['git']['branch']
    assert current['staged']==first['git']['staged']
    production=[];short=[]
    for path in OUT.rglob('started.json'):
        data=json.loads(path.read_text())
        (production if data['steps']==20000 else short).append(data['name'])
    assert production==['S3_legacy_upstream_circuit'],production
    assert json.loads((OUT/'r0_aperture/validation.json').read_text())['passed']
    assert json.loads((OUT/'s2_vs_s3/validation.json').read_text())['passed']
    # Validate hidden-state reconstruction also against directly observed
    # short-run Uu/Pu, not only the full-length boundary pressure records.
    replay=pd.read_csv(OUT/'s2_vs_s3/S2_acoustic_states_reconstructed.csv',float_precision='round_trip').iloc[:300]
    measured=pd.read_csv(OUT/'validation/s2_omitted_v2_300/output/csv/acoustic_trace.csv',float_precision='round_trip')
    fields=[c for c in replay if c.startswith(('Uu_','Pu_accum_','Pu_physical_'))]
    gate=numeric_gate(replay[fields],measured[fields],atol=1e-10,rtol=1e-9)
    dump(OUT/'s2_vs_s3/extended_validation.json',dict(passed=True,S2_Uu_Pu_short_direct_comparison=gate,
        reconstructed_fields=len(fields),short_rows=300))
    tests={}
    for label,command in [('cpp',['ctest','--test-dir',str(BUILD),'--output-on-failure']),
                           ('old_analysis',['python3','-B','tests/test_kanaya_ablation_analysis.py']),
                           ('aperture_analysis',['python3','-B','tests/test_r0_aperture_analysis.py']),
                           ('whitespace',['git','diff','--check'])]:
        result=subprocess.run(command,cwd=ROOT,env=environment(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (OUT/'validation'/(label+'_final.log')).write_text(result.stdout)
        assert result.returncode==0,result.stdout
        tests[label]={'passed':True,'command':command}
    dump(OUT/'preservation_audit.json',dict(passed=True,protected_worktrees_unchanged=True,
        existing_output_and_build_files_unchanged=True,checked_unchanged_files=inspected,
        allowed_task_source_edits=sorted(allowed),legacy_structure_contact_geometry_hashes_unchanged=paths,
        production_flow_NR_kernel_body_byte_identical=True,git_before=first['git'],git_after=current,
        final_binary_byte_identical_to_production=True,postrun_source_visibility_fix_only=postrun_changes,
        production_coupled_runs=len(production),production_cases=production,short_smoke_cases=short,
        R0_S2_new_coupled_runs=0,tests=tests,persistent_archive=str(ARCHIVE)))
    sources=[p for folder in ['src','include','tests','tools/kanaya_ablation','docs/plans']
             for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    sources.append(ROOT/'CMakeLists.txt')
    dump(OUT/'provenance.json',dict(worktree=str(ROOT),git=current,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},
        plan=json.loads((OUT/'plan_provenance.json').read_text()),
        legacy_source_commit='58cd92281aa6a36c732baa680d7817acbc98ddb2',
        compiler_flags=(BUILD/'CMakeFiles/simulation.dir/flags.make').read_text(),
        validation_gate_sha256=sha(OUT/'validation/preproduction_gate.json'),
        numerical_diagnostics='original NR loop unchanged; external post-step NR terms; sentinel iteration count -1'))
    print('PASS original/worktree/output preservation, single production count, source scope and all tests')


def archive():
    assert json.loads((OUT/'preservation_audit.json').read_text())['passed']
    destination=ARCHIVE/'next';destination.mkdir(exist_ok=False)
    rows=[];omitted=[]
    def copy(source,target):
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        digest=sha(source);assert sha(target)==digest
        rows.append(dict(source=str(source),archive=str(target.relative_to(ARCHIVE)),bytes=source.stat().st_size,sha256=digest))
    for folder in ['src','include','tests','tools','docs']:
        for source in (ROOT/folder).rglob('*'):
            if source.is_file() and '__pycache__' not in source.parts:copy(source,destination/'source'/source.relative_to(ROOT))
    for name in ['CMakeLists.txt','AGENTS.md','.gitignore','input/param.txt']:
        copy(ROOT/name,destination/'source'/name)
    for source in OUT.rglob('*'):
        if not source.is_file():continue
        if source.suffix=='.vtu':
            omitted.append(dict(path=str(source),bytes=source.stat().st_size,reason='generated short-smoke deformation VTU; hashes in test manifests, all-mode source retained'))
            continue
        copy(source,destination/'output/kanaya_ablation_next'/source.relative_to(OUT))
    for name in ['simulation','test_kanaya_diagnostics','kanaya_section_replay','test_upstream_acoustics']:
        copy(BUILD/name,destination/'binaries'/name)
    (destination/'task_diff.patch').write_text(git(ROOT)['diff'])
    with (ARCHIVE/'next_archive_manifest.json').open('x') as f:
        json.dump(dict(copied=rows,omitted=omitted,copied_bytes=sum(r['bytes'] for r in rows),
            previous_archive_untouched=True,source='uncommitted current worktree with prior implementation retained'),f,ensure_ascii=False,indent=2)
    print('Persistent final archive complete:',destination,'bytes',sum(r['bytes'] for r in rows))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['audit','archive']);a=p.parse_args()
    {'audit':audit,'archive':archive}[a.action]()
