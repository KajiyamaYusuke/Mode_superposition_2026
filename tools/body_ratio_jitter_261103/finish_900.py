#!/usr/bin/env python3
"""Monitor the explicitly requested 900 Pa sweep, then analyze saved cases."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import time
import run_sweep

ROOT=Path(__file__).resolve().parents[2]
STUDY=run_sweep.STUDY


def main():
    rows=[r for r in csv.DictReader((STUDY/'prepared_cases.csv').open()) if r['pressure_pa']=='900']
    if len(rows)!=18: raise ValueError('expected 18 cases at 900 Pa')
    while True:
        completed=0
        for row in rows:
            path=Path(row['parameter_file']).parent.parent/'status.json'
            if not path.exists():continue
            try: status=json.loads(path.read_text())
            except json.JSONDecodeError: continue
            state=status.get('simulation_status')
            if state=='RUN_FAILED':raise RuntimeError(row['case_id']+' failed: '+status.get('reason',''))
            if state=='SIMULATION_COMPLETED':completed+=1
            elif state!='RUNNING':raise RuntimeError('unexpected case state: '+str(state))
        if completed==18:break
        time.sleep(30)
    allrows=list(csv.DictReader((STUDY/'prepared_cases.csv').open()))
    config=json.loads((STUDY/'prepared_config.json').read_text());solver=run_sweep.audit(allrows,config)
    for row in rows:run_sweep.verified_completed(row,solver)
    print('All 18 900 Pa simulations passed numerical validation.',flush=True)
    for row in rows:
        case=Path(row['parameter_file']).parent.parent
        if not (case/'jitter_analysis').exists():
            subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/analyze_saved.py'),str(case)],cwd=ROOT,check=True)
    versions=set()
    for row in rows:
        case=Path(row['parameter_file']).parent.parent
        data=json.loads((case/'jitter_analysis/metrics.json').read_text())
        versions.add(json.dumps({'code':data.get('analysis_code_sha256'), 'config':data.get('analysis_config_sha256')},sort_keys=True))
    if len(versions)!=1:raise RuntimeError('analysis versions differ across cases; preserving outputs')
    subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/summarize_900.py')],cwd=ROOT,check=True)
    print('900 Pa analysis complete.',flush=True)

if __name__=='__main__':main()
