#!/usr/bin/env python3
"""Isolated, resumable Kanaya ablation workflow. Never writes to source worktrees."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(os.environ.get('KANAYA_CASE_ROOT', ROOT/'output/kanaya_ablation')).resolve()
COMPILER = Path(os.environ.get('KANAYA_COMPILER_ROOT',
    Path(shutil.which('icpx') or 'icpx').resolve().parent.parent))
PROTECTED = [Path(p).resolve() for p in os.environ.get('KANAYA_PROTECTED_ROOTS', str(ROOT)).split(os.pathsep) if p]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f: json.dump(value, f, ensure_ascii=False, indent=2)


def git(root):
    return {k: subprocess.check_output(['git', *v], cwd=root, text=True) for k, v in
            {'head': ['rev-parse', 'HEAD'], 'branch': ['branch', '--show-current'],
             'status': ['status', '--short'], 'diff': ['diff', '--binary'],
             'staged': ['diff', '--cached', '--binary']}.items()}


def protection():
    results = {}
    for root in PROTECTED:
        trees = {}
        for folder in ['src', 'include', 'input', 'output', 'result', 'build', 'tools']:
            rows = []
            for base, dirs, names in os.walk(root/folder):
                dirs.sort()
                for name in sorted(names):
                    p = Path(base)/name
                    st = p.stat()
                    rows.append([str(p.relative_to(root)), st.st_size, st.st_mtime_ns])
            trees[folder] = {'count': len(rows), 'metadata_sha256': hashlib.sha256(json.dumps(rows).encode()).hexdigest()}
        files = [p for folder in ['src', 'include'] for p in (root/folder).rglob('*') if p.is_file()]
        files += [p for p in [root/'CMakeLists.txt', root/'input/param.txt'] if p.exists()]
        results[str(root)] = {'git': git(root), 'trees': trees,
                             'source_input_sha256': {str(p): sha(p) for p in files}}
    return results


def environment(threads=1, cpus=None, diagnostics=True):
    e = os.environ.copy()
    e['PATH'] = str(COMPILER/'bin')+os.pathsep+e['PATH']
    e['LD_LIBRARY_PATH'] = str(COMPILER/'lib')+os.pathsep+e.get('LD_LIBRARY_PATH', '')
    e.update(OMP_NUM_THREADS=str(threads), OMP_PROC_BIND='close', OMP_PLACES='cores',
             OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1', SIMULATION_KANAYA_DIAGNOSTICS='1' if diagnostics else '0')
    return e


def parameters(nstep=20000, zeta=.06, extra=''):
    template = Path(os.environ.get('KANAYA_PARAMETER_TEMPLATE', PROTECTED[0]/'input/param.txt')).resolve()
    lines = template.read_text().splitlines()
    indices = [i for i, s in enumerate(lines) if s.strip() and not s.lstrip().startswith('#') and '=' not in s]
    for index, value in [(0, '20'), (2, str(nstep)), (5, str(zeta)), (6, str(zeta)), (15, '1600')]:
        lines[indices[index]] = value
    # Replace only path keys; preserve the actual protected physical input.
    for i, line in enumerate(lines):
        if '=' not in line or line.lstrip().startswith('#'): continue
        key, value = line.split('=', 1)
        if key.strip().lower() in ['leftfrequencyfile', 'rightfrequencyfile', 'leftmodefile',
                                   'rightmodefile', 'leftsurfacenasfile', 'rightsurfacenasfile']:
            lines[i] = key.strip()+' = '+str((template.parent/value.strip()).resolve())
    return '\n'.join(lines)+'\n'+extra+'\n'


def run_case(name, binary, *, steps=20000, threads=1, cpus=None, extra='', zeta=.06, diagnostics=True):
    case = OUT/name
    if case.exists(): raise FileExistsError(f'refusing existing case {case}')
    (case/'work').mkdir(parents=True)
    param = case/'params.txt'
    param.write_text(parameters(steps, zeta, extra))
    env = environment(threads, cpus, diagnostics)
    env['SIMULATION_RUN_DIR'] = str(case/'output')
    env['OMP_DISPLAY_AFFINITY'] = 'TRUE'
    command = [str(binary), '--smoke', str(param)]
    if cpus is not None: command = ['taskset', '-c', cpus]+command
    command = ['/usr/bin/time', '-v', '-o', str(case/'resource_usage.txt')]+command
    record = {'name': name, 'binary': str(binary), 'binary_sha256': sha(binary), 'parameters_sha256': sha(param),
              'nstep': steps, 'threads': threads, 'affinity': cpus,
              'environment': {k: v for k, v in env.items() if k.startswith(('OMP_', 'SIMULATION_', 'MKL_', 'OPENBLAS_'))},
              'command': command, 'started_unix': time.time(), 'git': git(ROOT)}
    dump(case/'started.json', record)
    print(f'Start {name}: threads={threads}, cpus={cpus}, steps={steps}', flush=True)
    with (case/'run.log').open('x') as log:
        result = subprocess.run(command, cwd=case/'work', env=env, stdout=log, stderr=subprocess.STDOUT)
    record.update(exit_code=result.returncode, elapsed_s=time.time()-record['started_unix'])
    manifest_txt=case/'output/txt/manifest.txt'
    if manifest_txt.exists():
        record['effective_parameters']=dict(line.split(' = ',1) for line in manifest_txt.read_text().splitlines() if ' = ' in line)
        paths=[Path(v) for k,v in record['effective_parameters'].items()
               if k in ['left_mode_vtu','right_mode_vtu','left_frequency','right_frequency','left_surface_nas','right_surface_nas']]
        record['model_inputs']={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in set(paths)}
    record['source_sha256']={str(p.relative_to(ROOT)):sha(p) for folder in ['src','include'] for p in (ROOT/folder).rglob('*') if p.is_file()}
    if binary.parent.name=='modern_baseline':
        # This binary was built before the ablation edits. Do not mislabel
        # the current workspace sources as the sources that produced it.
        revision='bc62cf8f3a2c01fb1ed0a33f6d4a6864f6c1d2d5'
        record['workspace_source_sha256_at_capture']=record['source_sha256']
        files=subprocess.check_output(['git','ls-tree','-r','--name-only',revision,'src','include'],
                                      cwd=ROOT,text=True).splitlines()
        record['source_sha256']={name:hashlib.sha256(subprocess.check_output(
            ['git','show',revision+':'+name],cwd=ROOT)).hexdigest() for name in files}
        record['binary_source_revision']=revision
        original=json.loads((OUT/'validation/untouched_80/manifest.json').read_text())
        assert record['binary_sha256']==original['binary_sha256'], 'pristine reference binary changed'
    record['legacy_reference_sha256']={p.name:sha(p) for p in (OUT/'references').glob('user_main_*')}
    record['legacy_source_commit']=json.loads((OUT/'references/user_reference.json').read_text())
    record['outputs'] = {str(p.relative_to(case)): {'bytes': p.stat().st_size, 'sha256': sha(p)}
                         for p in (case/'output').rglob('*') if p.is_file() and 'result' not in p.parts}
    dump(case/'manifest.json', record)
    print(f'Finish {name}: code={result.returncode}, wall={record["elapsed_s"]:.3f}s', flush=True)
    if result.returncode: raise RuntimeError(f'case failed: {case}; inspect run.log')
    return record


def prepare():
    dump(OUT/'pre_task_protection.json', protection())
    dump(OUT/'start.json', {'cwd': str(Path.cwd()), 'git': git(ROOT),
                           'CPU': subprocess.check_output(['lscpu'], text=True),
                           'allowed': sorted(os.sched_getaffinity(0)),
                           'cpu_max': Path('/sys/fs/cgroup/cpu.max').read_text(),
                           'actual_run_parameter_template': parameters()})
    (ROOT/'input').mkdir(exist_ok=True)
    (ROOT/'input/ablation_reference.txt').write_text(parameters())
    references = OUT/'references'; references.mkdir()
    for revision, names in [('6baf527', ['src/forceCalculator.cpp', 'src/state.cpp', 'src/simulation.cpp',
                                        'include/SimulationParams.h', 'input/param.txt']),
                             ('bc62cf8', ['src/structure/rungeStep.cpp'])]:
        for name in names:
            content = subprocess.check_output(['git', 'show', revision+':'+name], cwd=ROOT)
            (references/(revision+'_'+name.replace('/', '_'))).write_bytes(content)
    print('Protected original worktrees and prepared independent inputs.', flush=True)


def build(directory='modern_baseline', diagnostics=True):
    commands = [['cmake', '-S', str(ROOT), '-B', str(ROOT/'build'/directory), '-DCMAKE_BUILD_TYPE=Release',
                 '-DENABLE_KANAYA_DIAGNOSTICS='+('ON' if diagnostics else 'OFF')],
                ['cmake', '--build', str(ROOT/'build'/directory), '-j', str(len(os.sched_getaffinity(0)))]]
    for cmd in commands: subprocess.run(cmd, cwd=ROOT, env=environment(), check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'baseline-build', 'build', 'audit'])
    args = parser.parse_args()
    if args.action == 'prepare': prepare()
    elif args.action == 'baseline-build': build()
    elif args.action == 'build': build('switchable')
    else:
        before = json.loads((OUT/'pre_task_protection.json').read_text())
        after = protection()
        assert before == after, 'protected original changed'
        dump(OUT/'final_preservation_audit.json', {'passed': True, 'snapshot': after})
        print('All protected source/input/output/build trees and Git states unchanged.')


if __name__ == '__main__': main()
