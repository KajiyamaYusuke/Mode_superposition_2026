"""Isolated R0 replay/S3 workflow; old outputs and binaries are read-only."""
import json
import os
from pathlib import Path
import subprocess
import time
from run_matrix import ROOT, OUT as PREVIOUS, sha, dump, environment, git

OUT=Path(os.environ.get('KANAYA_NEXT_ROOT',ROOT/'output/kanaya_ablation_next')).resolve()
ARCHIVE=Path(os.environ.get('KANAYA_ARCHIVE_ROOT',ROOT/'output/kanaya_archive')).resolve()
BUILD=Path(os.environ.get('KANAYA_BUILD_DIR',ROOT/'build/diagnostic')).resolve()
BIN=BUILD/'test_kanaya_diagnostics'
OLD_BIN=ARCHIVE/'previous/build/switchable/test_kanaya_diagnostics'
S2=PREVIOUS/'S2_legacy_vt0175'


def physical_parameters(path):
    tokens=[s.strip() for s in Path(path).read_text().splitlines()
            if s.strip() and not s.lstrip().startswith('#') and '=' not in s]
    return dict(nmode=int(tokens[0]),nstep=int(tokens[2]),dt_s=float(tokens[4]),
                zetaL=float(tokens[5]),zetaR=float(tokens[6]),ps_Pa=float(tokens[15]),
                rho_kg_m3=float(tokens[16]),mu_Pa_s=float(tokens[17]),mass_kg=float(tokens[18]),
                c_sound_m_s=float(tokens[19]),iforce=int(tokens[20]))


def copied_parameters(source,steps=None,extra=''):
    lines=Path(source).read_text().splitlines()
    indices=[i for i,s in enumerate(lines) if s.strip() and not s.lstrip().startswith('#') and '=' not in s]
    if steps is not None:lines[indices[2]]=str(steps)
    return '\n'.join(lines)+'\n'+extra


def run(name,*,source=S2/'params.txt',steps=300,extra='',binary=BIN,trace=True,threads=4,cpus='0,2,4,6'):
    case=OUT/name
    if case.exists():raise FileExistsError(case)
    (case/'work').mkdir(parents=True)
    param=case/'params.txt';param.write_text(copied_parameters(source,steps,extra))
    env=environment(threads)
    env.update(SIMULATION_RUN_DIR=str(case/'output'),OMP_DISPLAY_AFFINITY='TRUE',BLIS_NUM_THREADS='1')
    if trace:env['SIMULATION_ACOUSTIC_TRACE']='1'
    else:env.pop('SIMULATION_ACOUSTIC_TRACE',None)
    command=['/usr/bin/time','-v','-o',str(case/'resource_usage.txt'),'taskset','-c',cpus,str(binary),'--smoke',str(param)]
    record=dict(name=name,binary=str(binary),binary_sha256=sha(binary),parameters_sha256=sha(param),
        physical_parameters=physical_parameters(param),
        parent_parameters=str(source),parent_parameters_sha256=sha(source),
        steps=steps,threads=threads,affinity=cpus,command=command,started_unix=time.time(),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for folder in ['src','include'] for p in (ROOT/folder).rglob('*') if p.is_file()})
    dump(case/'started.json',record)
    print('Start',name,'steps',steps,'threads',threads,flush=True)
    with (case/'run.log').open('x') as log:
        result=subprocess.run(command,cwd=case/'work',env=env,stdout=log,stderr=subprocess.STDOUT)
    record.update(exit_code=result.returncode,elapsed_s=time.time()-record['started_unix'])
    path=case/'output/txt/manifest.txt'
    if path.exists():record['effective_parameters']=dict(s.split(' = ',1) for s in path.read_text().splitlines() if ' = ' in s)
    record['outputs']={str(p.relative_to(case)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in (case/'output').rglob('*') if p.is_file()}
    record['model_inputs']=json.loads((S2/'manifest.json').read_text())['model_inputs']
    record['legacy_source_commit']='58cd92281aa6a36c732baa680d7817acbc98ddb2'
    dump(case/'manifest.json',record)
    print('Finish',name,'exit',result.returncode,'wall',record['elapsed_s'],flush=True)
    if result.returncode:raise RuntimeError(f'{name} failed; inspect preserved run.log')
    return record


def prepare():
    OUT.mkdir(parents=True,exist_ok=False)
    before=json.loads((ARCHIVE/'before.json').read_text())
    dump(OUT/'before.json',before)
    affinity=sorted(os.sched_getaffinity(0))
    topology=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True)
    distinct={}
    for line in topology.splitlines():
        if line.startswith('#'):continue
        cpu,core,socket=map(int,line.split(','))
        if cpu in affinity:distinct.setdefault((socket,core),cpu)
    cpus=list(distinct.values())[:8]
    assert len(cpus)==8,'recheck scheduling; physical-core exposure changed'
    dump(OUT/'benchmark_decision.json',dict(threads=8,affinity=','.join(map(str,cpus)),
        basis='previous measured fastest single-job 8 threads; current topology and cpuset rechecked; no benchmark sweep repeated',
        allowed=affinity,topology=topology,cpu_max=Path('/sys/fs/cgroup/cpu.max').read_text(),
        previous_benchmark_sha256=sha(PREVIOUS/'benchmark_final_selection.json'),
        prior_single_modern_s={'4':49.199830,'8':47.167546,'16':51.317304}))
    dump(OUT/'plan_provenance.json',dict(repository='KajiyamaYusuke/Mode_superposition_2026',
        branch='codex/kanaya-force-energy-implementation-20261021',remote_blob='22e0ee741950d54bbe705c03116e33664c4dd5e9',
        local_git_blob=subprocess.check_output(['git','hash-object','docs/plans/KANAYA_R0_RESIDUAL_GAP_AND_S3_LEGACY_UPSTREAM_20261021.md'],cwd=ROOT,text=True).strip(),
        persistent_archive=str(ARCHIVE),archive_manifest_sha256=sha(ARCHIVE/'archive_manifest.json')))
    # Test fixture inputs are copied from the ACTUAL S2, not production param.
    fixtures=OUT/'validation';fixtures.mkdir()
    (fixtures/'s2_params.txt').write_text(copied_parameters(S2/'params.txt',300))
    (fixtures/'s3_params.txt').write_text(copied_parameters(S2/'params.txt',300,'upstreamAcousticProfile = legacy_main\n'))
    (fixtures/'bad_profile.txt').write_text(copied_parameters(S2/'params.txt',3,'upstreamAcousticProfile = nonsense\n'))


if __name__=='__main__':prepare()
