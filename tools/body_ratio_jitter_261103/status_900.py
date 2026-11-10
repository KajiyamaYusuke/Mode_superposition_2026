#!/usr/bin/env python3
"""Show read-only progress of the 900 Pa body-ratio cases."""
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]/'analysis_runs/jitter_body_ratio_261103_900_1400_1s'


def main():
    rows=[r for r in csv.DictReader((ROOT/'prepared_cases.csv').open()) if r['pressure_pa']=='900']
    print('case_id,simulation_status,last_saved_time_s')
    for row in rows:
        case=Path(row['parameter_file']).parent.parent;path=case/'status.json';state='NOT_STARTED';t=''
        if path.exists():
            try:state=json.loads(path.read_text()).get('simulation_status','UNKNOWN')
            except (OSError,json.JSONDecodeError):state='STATUS_UPDATING'
        flow=case/'output/dat/airflow_vt.dat'
        if flow.exists():
            with flow.open('rb') as f:
                f.seek(max(0,flow.stat().st_size-8192));lines=f.read().splitlines()
            for line in reversed(lines):
                fields=line.split()
                if len(fields)==2:
                    try:t=f'{int(fields[0])*1e-5:.5f}';break
                    except ValueError:pass
        print(f"{row['case_id']},{state},{t}")

if __name__=='__main__':main()
