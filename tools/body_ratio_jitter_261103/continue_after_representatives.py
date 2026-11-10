#!/usr/bin/env python3
"""Wait for representatives; validate before sequential remaining simulations."""
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
import run_sweep

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'

def main():
    rows=list(csv.DictReader((STUDY/'prepared_cases.csv').open()))
    config=json.loads((STUDY/'prepared_config.json').read_text())
    for case in ('r0p50_P900','r0p50_P1400'):
        status=STUDY/'cases'/case/'status.json'
        while True:
            try:
                data=json.loads(status.read_text())
            except (FileNotFoundError, json.JSONDecodeError):
                time.sleep(20)
                continue
            if data.get('simulation_status')!='RUNNING': break
            time.sleep(20)
        if data.get('simulation_status')!='SIMULATION_COMPLETED':raise RuntimeError(case+' representative failed: '+str(data))
    solver=run_sweep.audit(rows,config)
    for row in rows:
        if row['case_id'] in ('r0p50_P900','r0p50_P1400'): run_sweep.verified_completed(row,solver)
    subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/check_representatives.py')],cwd=ROOT,check=True)
    print('Representative numerical and flow checks passed; starting remaining cases.',flush=True)
    subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/run_sweep.py'),'--execute','--remaining','--resume','--threads','4'],cwd=ROOT,check=True)
    print('All cases completed; analyzing saved flow.',flush=True)
    for row in rows:
        case=Path(row['parameter_file']).parent.parent
        if not (case/'jitter_analysis').exists():
            subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/analyze_saved.py'),str(case)],cwd=ROOT,check=True)
    subprocess.run([sys.executable,str(ROOT/'tools/body_ratio_jitter_261103/plot_jitter_sweep.py')],cwd=ROOT,check=True)
    print('All cases and preliminary analyses completed.',flush=True)

if __name__=='__main__':main()
