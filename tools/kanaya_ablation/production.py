#!/usr/bin/env python3
"""Run only selected physical cases with the measured, nonoversubscribed schedule."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from benchmark import BIN, NO_VTU, OLD, check_gate
from run_matrix import OUT, ROOT, run_case, dump


def run_batch(cases):
    check_gate()
    selection=json.loads((OUT/'benchmark_final_selection.json').read_text())
    assert selection['passed']
    winner=selection['winner'];threads=winner['threads_per_job']
    jobs=winner['jobs']
    if jobs==2:
        affinities=[r['cpus'] for r in winner['job_records']]
    else:
        affinities=[','.join(str(v) for v in range(0,2*threads,2))]
    results=[]
    for start in range(0,len(cases),jobs):
        chunk=cases[start:start+jobs]
        with ThreadPoolExecutor(len(chunk)) as pool:
            pending=[]
            for i,(name,extra,zeta) in enumerate(chunk):
                path=OUT/name/'manifest.json'
                if path.exists():
                    record=json.loads(path.read_text());assert record['exit_code']==0
                    results.append(record);continue
                pending.append(pool.submit(run_case,name,BIN,steps=20000,
                    threads=threads,cpus=affinities[i],extra=NO_VTU+extra,zeta=zeta))
            results.extend(p.result() for p in pending)
    return results


def custom(**changes):
    flags=['legacyArea','legacyFlowCoefficients','legacyWallViscosity',
           'legacySeparation','legacyPostSeparation','legacySurfacePressureLoad']
    assert set(changes)<=set(flags)
    return 'glottalMethodProfile = custom\n'+''.join(f'{k} = {int(changes.get(k,True))}\n' for k in flags)


def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['primary','targeted','weak','fluid','baseline'])
    args=p.parse_args()
    if args.stage=='primary':
        cases=[('R0_modern','',.06),('R1_legacy_fluid_all',OLD,.06)]
    elif args.stage=='targeted':
        cases=[('R3_modern_coefficients',custom(legacyFlowCoefficients=False),.06),
               ('R6_modern_surface_mapping',custom(legacySurfacePressureLoad=False),.06)]
    elif args.stage=='weak':
        # Explicit single-factor comparisons against R1, not simultaneous
        # changes to fluid mechanics, damping and acoustics.
        cases=[('S1_legacy_zeta010',OLD,.10),
               ('S2_legacy_vt0175',OLD+'vocalTractLengthM = 0.175\n',.06)]
    elif args.stage=='fluid':
        cases=[('R2_modern_area_mapping',custom(legacyArea=False,legacySurfacePressureLoad=False),.06),
               ('R3_modern_coefficients',custom(legacyFlowCoefficients=False),.06),
               ('R4_modern_wall_separation',custom(legacyWallViscosity=False,legacySeparation=False),.06),
               ('R5_modern_post_separation',custom(legacyPostSeparation=False),.06),
               ('R6_modern_surface_mapping',custom(legacySurfacePressureLoad=False),.06)]
    else:
        check_gate()
        selection=json.loads((OUT/'benchmark_final_selection.json').read_text())
        single=selection['fastest_single_modern']
        record=run_case('validation/pristine_modern_full',ROOT/'build/modern_baseline/test_kanaya_diagnostics',
                        threads=single['threads'],cpus=single['cpus'])
        dump(OUT/'pristine_full_run.json',{'case':record['name'],'passed':record['exit_code']==0})
        return
    records=run_batch(cases)
    dump(OUT/(args.stage+'_runs.json'),{'stage':args.stage,'cases':[r['name'] for r in records],
         'measured_schedule':json.loads((OUT/'benchmark_final_selection.json').read_text())['winner']})


if __name__=='__main__':main()
