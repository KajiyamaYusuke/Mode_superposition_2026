#!/usr/bin/env python3
"""Run isolated sweep cases after immutable-input and completion checks."""
import argparse
import csv
import datetime
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def save_status(path, data):
    tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n')
    tmp.replace(path)


def inputs(text):
    result={}
    for name in ('leftFrequencyFile','rightFrequencyFile','leftModeFile','rightModeFile','leftSurfaceNasFile','rightSurfaceNasFile'):
        m=re.findall(r'^'+name+r'\s*=\s*(.+?)\s*$',text,re.M)
        if len(m)!=1 or not Path(m[0]).is_absolute() or not Path(m[0]).is_file(): raise ValueError('bad input '+name)
        result[name]=Path(m[0])
    return result


def audit(rows, config):
    if len(rows)!=36: raise ValueError('expected 36 cases')
    solver=Path(config['solver']);
    if sha(solver)!=config['solver_sha256']: raise ValueError('solver binary changed')
    canonical=set(); combos=set(); left=set()
    for row in rows:
        path=Path(row['parameter_file']); text=path.read_text(); p=inputs(text)
        if sha(path)!=row['parameter_sha256'] or sha(p['rightFrequencyFile'])!=row['frequency_sha256'] or sha(p['rightModeFile'])!=row['mode_sha256']: raise ValueError('input changed: '+row['case_id'])
        if '# number of step nstep\n100000\n' not in text or '# time step for time integration dt\n1.0E-5\n' not in text or text.count('writeDeformationVtu = 0')!=1 or '# subglottal pressure Ps (Pa)\n'+row['pressure_pa']+'\n' not in text: raise ValueError('invalid settings: '+row['case_id'])
        token=re.fullmatch(r'b(\d+)(?:p(\d+))?c(\d+)',p['rightFrequencyFile'].stem)
        if not token or p['rightModeFile'].stem!=p['rightFrequencyFile'].stem: raise ValueError('right mode/frequency token mismatch')
        value=Decimal(token.group(1)+(('.'+token.group(2)) if token.group(2) else ''))
        if value!=12*Decimal(row['ratio']) or token.group(3)!='3' or row['pressure_pa'] not in ('900','1400'): raise ValueError('ratio or cover mismatch: '+row['case_id'])
        surface=ROOT/'input/M5_test/M5_surface_T3_d2.nas'
        if p['leftSurfaceNasFile']!=surface or p['rightSurfaceNasFile']!=surface: raise ValueError('unexpected surface')
        normalized=re.sub(r'(?<=# subglottal pressure Ps \(Pa\)\n)\d+','PS',text)
        normalized=re.sub(r'^(rightFrequencyFile|rightModeFile)\s*=.*$',r'\1 = RIGHT',normalized,flags=re.M)
        canonical.add(hashlib.sha256(normalized.encode()).hexdigest()); combos.add((row['ratio'],row['pressure_pa']))
        left.add(tuple(sha(p[k]) for k in ('leftFrequencyFile','leftModeFile','leftSurfaceNasFile','rightSurfaceNasFile')))
    expected={(str(Decimal(i)/20),str(pressure)) for i in range(2,20) for pressure in (900,1400)}
    if len(canonical)!=1 or combos!=expected or len(left)!=1: raise ValueError('case configurations, ratio set or left inputs differ')
    if shutil.disk_usage(STUDY).free<5*1024**3: raise ValueError('less than 5 GiB free')
    return solver


def validate(run, pressure):
    file=run/'dat/airflow_vt.dat'
    with file.open() as f: header=f.readline()
    if not header.startswith('#') or 'airflow' not in header: raise ValueError('invalid flow header')
    a=np.loadtxt(file,comments='#',ndmin=2); expected=np.arange(0,100000,5)
    if a.shape!=(20000,2) or not np.array_equal(a[:,0],expected) or not np.isfinite(a).all(): raise ValueError('flow samples missing or nonfinite')
    y=a[:,1]; dif=np.abs(np.diff(y)); maxflow=float(np.max(np.abs(y))); jump=float(np.max(dif)); p99=float(np.quantile(dif,.99))
    # A conservative gross-failure screen. Smaller jumps remain visible in QC plots.
    jump_limit=max(.001,20*p99)
    if maxflow>.01 or jump>jump_limit: raise ValueError('flow magnitude or jump requires inspection')
    for name in ('dat/area.dat','dat/displace.dat','dat/pressure_vt.dat','wav/test_sound.wav','txt/params_used.txt','txt/manifest.txt'):
        if not (run/name).is_file() or (run/name).stat().st_size==0: raise ValueError('missing output '+name)
    for name in ('dat/displace.dat','dat/displace_xy.dat','dat/pressure_vt.dat'):
        scalar=np.loadtxt(run/name,comments='#',ndmin=2)
        if scalar.shape[0] != 20000 or not np.isfinite(scalar).all():
            raise ValueError('missing or nonfinite scalar output '+name)
        if name.endswith('pressure_vt.dat'):
            if not np.array_equal(scalar[:,0],expected): raise ValueError('invalid pressure sample steps')
        elif not np.allclose(scalar[:,0],expected*1e-5,rtol=0,atol=1e-9):
            raise ValueError('invalid displacement sample times')
    manifest=(run/'txt/manifest.txt').read_text()
    for line in ('nstep = 100000','dt_s = 1e-05','nmode = 20','write_deformation_vtu = 0'):
        if line not in manifest: raise ValueError('manifest mismatch: '+line)
    # Legacy cosine pressure is recorded in the solver's initialization log,
    # not the manifest of the current solver revision.
    stdout=(run.parent/'stdout.log').read_text()
    match=re.search(r'^\s*ps\s*=\s*([\d.eE+-]+)\s*\[Pa\]',stdout,re.M)
    if not match or float(match.group(1))!=float(pressure): raise ValueError('solver pressure log mismatch')
    return dict(flow_samples=len(y),first_step=int(a[0,0]),last_step=int(a[-1,0]),max_abs_flow_m3_s=maxflow,max_sample_jump_m3_s=jump,p99_sample_jump_m3_s=p99,gross_jump_limit_m3_s=jump_limit,flow_sha256=sha(file))


def run_case(row,solver,threads):
    case=Path(row['parameter_file']).parent.parent; status=case/'status.json'; run=case/'output'
    if shutil.disk_usage(STUDY).free<5*1024**3: raise RuntimeError('less than 5 GiB free')
    if status.exists() or run.exists() or (case/'stdout.log').exists() or (case/'stderr.log').exists(): raise FileExistsError('existing case preserved: '+str(case))
    text=Path(row['parameter_file']).read_text(); inp=inputs(text); hashes={name:sha(path) for name,path in inp.items()}
    data=dict(case_id=row['case_id'],simulation_status='RUNNING',jitter_status='NOT_ANALYZED',parameter_sha256=row['parameter_sha256'],solver_sha256=sha(solver),input_sha256=hashes,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),threads=threads)
    save_status(status,data); env=os.environ.copy(); env.update(OMP_NUM_THREADS=str(threads),OMP_PROC_BIND='close',OMP_PLACES='cores',SIMULATION_RUN_DIR=str(run)); started=time.monotonic()
    try:
        with (case/'stdout.log').open('x') as out,(case/'stderr.log').open('x') as err:
            process=subprocess.Popen([str(solver),row['parameter_file']],cwd=ROOT/'build',env=env,stdout=out,stderr=err)
            position=0; pending=b''; previous_step=-5
            while process.poll() is None:
                time.sleep(10)
                path=run/'dat/airflow_vt.dat'
                if not path.exists(): continue
                with path.open('rb') as stream:
                    stream.seek(position);chunk=stream.read();position=stream.tell()
                pieces=(pending+chunk).split(b'\n');pending=pieces.pop()
                try:
                    for line in pieces:
                        if not line.strip() or line.startswith(b'#'): continue
                        parts=line.split()
                        if len(parts)!=2: raise ValueError('invalid streamed flow row')
                        step,value=int(parts[0]),float(parts[1])
                        if step!=previous_step+5 or not np.isfinite(value) or abs(value)>.01:
                            raise ValueError('invalid or divergent streamed flow')
                        previous_step=step
                except Exception:
                    process.terminate()
                    try: process.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        process.kill();process.wait()
                    raise
            data['return_code']=process.returncode
        if process.returncode!=0: raise RuntimeError('solver returned '+str(process.returncode))
        if sha(Path(row['parameter_file'])) != row['parameter_sha256'] or any(sha(path)!=hashes[name] for name,path in inp.items()): raise ValueError('input changed during run')
        if sha(run/'txt/params_used.txt') != row['parameter_sha256']: raise ValueError('recorded parameter file differs')
        data['validation']=validate(run,row['pressure_pa']); data['simulation_status']='SIMULATION_COMPLETED'
    except Exception as exc:
        data['simulation_status']='RUN_FAILED';data['reason']=str(exc)
    data['elapsed_s']=time.monotonic()-started;data['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save_status(status,data);print(row['case_id'],data['simulation_status'],data.get('reason',''),flush=True)
    if data['simulation_status']!='SIMULATION_COMPLETED': raise RuntimeError('case stopped: '+row['case_id'])


def verified_completed(row, solver):
    case=Path(row['parameter_file']).parent.parent
    status=json.loads((case/'status.json').read_text())
    if status.get('simulation_status')!='SIMULATION_COMPLETED' or status.get('parameter_sha256')!=row['parameter_sha256'] or status.get('solver_sha256')!=sha(solver):
        raise ValueError('completed case identity/status mismatch: '+row['case_id'])
    current={name:sha(path) for name,path in inputs(Path(row['parameter_file']).read_text()).items()}
    if current != status.get('input_sha256'):
        raise ValueError('completed input hashes changed: '+row['case_id'])
    now=validate(case/'output',row['pressure_pa'])
    if now['flow_sha256']!=status['validation']['flow_sha256']:
        raise ValueError('completed flow hash changed: '+row['case_id'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');g=p.add_mutually_exclusive_group(required=True);g.add_argument('--representatives',action='store_true');g.add_argument('--remaining',action='store_true');g.add_argument('--pressure-900',action='store_true',help='Run only the 900 Pa series after validating its representative');p.add_argument('--resume',action='store_true');p.add_argument('--threads',type=int,default=4);a=p.parse_args()
    if not a.execute or a.threads<=0: p.error('specify --execute and positive --threads')
    rows=list(csv.DictReader((STUDY/'prepared_cases.csv').open()));config=json.loads((STUDY/'prepared_config.json').read_text());solver=audit(rows,config)
    representative_ids=('r0p50_P900','r0p50_P1400')
    representatives=[]
    for ident in representative_ids:
        matches=[r for r in rows if r['case_id']==ident]
        if len(matches)!=1: raise ValueError('representative missing')
        representatives.append(matches[0])
    if a.pressure_900:
        verified_completed(representatives[0],solver)
        selected=[r for r in rows if r['pressure_pa']=='900']
    elif a.remaining:
        for row in representatives: verified_completed(row,solver)
        selected=[r for r in rows if r['case_id'] not in representative_ids]
    else:
        selected=representatives
    for row in selected:
        status=Path(row['parameter_file']).parent.parent/'status.json'
        if a.resume and status.exists():
            verified_completed(row,solver)
            print(row['case_id'],'verified completed; skipping',flush=True)
            continue
        run_case(row,solver,a.threads)

if __name__=='__main__': main()
