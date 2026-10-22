#!/usr/bin/env python3
"""Pre-production stop gates for the single S3 production case."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import json
import subprocess
import numpy as np
import pandas as pd
from next_common import ROOT,OUT,PREVIOUS,BUILD,BIN,OLD_BIN,S2,run,environment,dump,sha
from validate import compare_case,extract


def resume_smoke(name,**kwargs):
    path=OUT/name/'manifest.json'
    if path.exists():
        value=json.loads(path.read_text())
        assert value['exit_code']==0 and value['binary_sha256']==sha(kwargs.get('binary',BIN))
        return value
    return run(name,**kwargs)


def smoke():
    cases=[('validation/s2_omitted_v2_300',{}),
           ('validation/s2_explicit_v2_300',dict(extra='upstreamAcousticProfile = modern\n')),
           ('validation/s3_legacy_v2_300',dict(extra='upstreamAcousticProfile = legacy_main\n')),
           ('validation/s3_open_v2_300',dict(extra='upstreamAcousticProfile = legacy_main\ninitialGapMm = 0.1\n'))]
    with ThreadPoolExecutor(4) as pool:
        futures=[pool.submit(resume_smoke,name,threads=2,cpus=','.join(map(str,[4*i,4*i+2])),**kw)
                 for i,(name,kw) in enumerate(cases)]
        [f.result() for f in futures]
    comparisons={}
    for name in ['s2_omitted_v2_300','s2_explicit_v2_300']:
        comparisons[name]=compare_case(OUT/'validation/s2_old_300/output',OUT/'validation'/name/'output',skip_vtu=True)
        assert comparisons[name]['passed'],comparisons
    modern_source=PREVIOUS/'validation/modern_boundary_open_80/params.txt'
    resume_smoke('validation/r0_open_80',source=modern_source,steps=80,trace=False)
    comparisons['modern_open']=compare_case(PREVIOUS/'validation/modern_boundary_open_80/output',OUT/'validation/r0_open_80/output')
    assert comparisons['modern_open']['passed'],comparisons
    for name,ng in [('s2_omitted_v2_300',10),('s2_explicit_v2_300',10),('s3_legacy_v2_300',3),('s3_open_v2_300',3)]:
        path=OUT/'validation'/name
        frame=pd.read_csv(path/'output/csv/acoustic_trace.csv',float_precision='round_trip')
        assert len(frame)==300 and np.isfinite(frame.select_dtypes('number')).all().all()
        assert len([c for c in frame if c.startswith('Uu_')])==ng+1
        assert len([c for c in frame if c.startswith('Pu_accum_')])==ng+2
        active=frame.NR_active==1
        if name=='s3_legacy_v2_300' and not active.any():
            assert frame.Ug_m3_s.eq(0).all() # zero-gap startup remains closed during 3 ms
        else:
            assert active.any() and frame[active].NR_converged.eq(1).all()
            assert frame[active].NR_residual_Pa.abs().max()<1e-8
        meta=json.loads((path/'manifest.json').read_text())['effective_parameters']
        assert int(meta['subglottal_sections'])==ng and float(meta['vocal_tract_length_m'])==.175
    invalid=OUT/'validation/bad_profile_output';invalid.mkdir()
    env=environment();env['SIMULATION_RUN_DIR']=str(invalid)
    bad=subprocess.run([str(BIN),'--smoke',str(OUT/'validation/bad_profile.txt')],env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/'validation/bad_profile.log').write_text(bad.stdout)
    assert bad.returncode!=0 and 'unknown upstreamAcousticProfile' in bad.stdout
    oracle=coefficient_oracle()
    dump(OUT/'validation/preproduction_gate.json',dict(passed=True,regressions=comparisons,
        partition_bounds_checked=[3,10],bad_profile_rejected=True,coefficient_oracle=oracle,
        production_coupled_runs=0,short_smokes={'accepted_300_steps':4,'accepted_80_steps':1},
        extra_test_condition='S3 initialGapMm=.1 ONLY in open-flow validation; production inherits S2 zero gap',
        NR_iteration_count='not instrumented (-1); preserves original NR loop rounding'))
    print('PASS S3 pre-production gate: bitwise defaults/modern, 3/10 partitions, old coefficients, NR, invalid profile')


def coefficient_oracle():
    # Compile the ORIGINAL initialization method with original headers, not
    # a second hand transcription of its nine coefficient formulas.
    directory=OUT/'validation/original_acoustic_oracle';directory.mkdir()
    refs=PREVIOUS/'references'
    import shutil
    for name in ['Geometry.h','ForceCalculator.h','ModeData.h','State.h','Displacement.h','SimulationParams.h']:
        shutil.copyfile(refs/('user_main_'+name),directory/name)
    source=(refs/'user_main_forceCalculator.cpp').read_text()
    methods='\n\n'.join(extract(source,sig) for sig in ['ForceCalculator::ForceCalculator(','void ForceCalculator::initialize('])
    code='#include <bits/stdc++.h>\n#define private public\n#include "ForceCalculator.h"\n#undef private\n'+methods+'''\n
int main(){ Geometry g;g.nPoints=2;g.nsurfl=7;g.nsurfz=5;g.nxsup=7;g.zmax=17.;
ModeData md;md.nModes=20;State st;st.nSteps=300;SimulationParams p;
p.rho=1.13;p.mu=1.83e-5;p.ps=1600.;p.c_sound=343.;
ForceCalculator f(g,md,st,p);f.initialize();
std::cout<<"ORACLE "<<std::setprecision(17)<<f.Lui<<" "<<f.Cui<<" "<<f.Lu<<" "<<f.Cu<<" "<<f.R2<<" "<<f.La<<" "<<f.Ca<<" "<<f.Lr<<" "<<f.Rr<<"\\n";}
'''
    (directory/'oracle.cpp').write_text(code)
    command=['/opt/intel/oneapi/compiler/2025.1/bin/icpx','-O3','-std=c++17','-I'+str(directory),str(directory/'oracle.cpp'),'-o',str(directory/'oracle')]
    subprocess.run(command,env=environment(),check=True)
    output=subprocess.check_output([str(directory/'oracle')],env=environment(),text=True)
    (directory/'output.log').write_text(output)
    values=np.fromstring(next(s for s in output.splitlines() if s.startswith('ORACLE '))[7:],sep=' ')
    actual=pd.read_csv(OUT/'validation/s3_legacy_v2_300/output/csv/acoustic_coefficients.csv',float_precision='round_trip')
    assert len(values)==9 and np.allclose(actual.value,values,rtol=1e-14,atol=1e-25)
    actual['original_value']=values;actual['absolute_error']=np.abs(actual.value-values)
    actual.to_csv(directory/'coefficient_comparison.csv',index=False)
    return dict(passed=True,maximum_relative_error=float(np.max(np.abs(actual.value-values)/np.abs(values))),
                source_sha256=sha(refs/'user_main_forceCalculator.cpp'),constructor_and_initialize_verbatim=True)


def production():
    gate=json.loads((OUT/'validation/preproduction_gate.json').read_text());assert gate['passed']
    decision=json.loads((OUT/'benchmark_decision.json').read_text())
    run('S3_legacy_upstream_circuit',steps=20000,extra='upstreamAcousticProfile = legacy_main\n',
        threads=decision['threads'],cpus=decision['affinity'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['smoke','production']);a=p.parse_args()
    {'smoke':smoke,'production':production}[a.action]()
