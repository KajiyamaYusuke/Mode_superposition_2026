#!/usr/bin/env python3
"""Measured CPU/affinity selection, including actual two-job throughput."""
from concurrent.futures import ThreadPoolExecutor
import json
import re
import time
import numpy as np
import pandas as pd
from run_matrix import ROOT, OUT, run_case, dump, sha

BIN=ROOT/'build/switchable/test_kanaya_diagnostics'
NO_VTU='writeDeformationVtu = 0\n'
OLD='glottalMethodProfile = legacy_fluid_all\n'


def check_gate():
    for path in ['validation/fixed_final/validation.json','validation/parser_final/validation.json',
                 'validation/modern_regression.json','validation/modern_open_regression.json']:
        values=json.loads((OUT/path).read_text())
        if 'passed' in values: assert values['passed']
        else: assert all(v['passed'] for v in values.values())


def metrics(record):
    path=OUT/record['name']
    text=(path/'resource_usage.txt').read_text()
    def value(pattern):
        return float(re.search(pattern,text)[1])
    timings=re.findall(r'^(calcArea|calcForce|f2mode|mode2uf|calcDis|output)\s*:\s*([0-9.eE+\-]+)\s*ms',
                       (path/'run.log').read_text(),re.M)
    return {'name':record['name'],'threads':record['threads'],'cpus':record['affinity'],
            'elapsed_s':record['elapsed_s'],'user_s':value(r'User time \(seconds\): ([0-9.]+)'),
            'system_s':value(r'System time \(seconds\): ([0-9.]+)'),
            'peak_rss_kb':value(r'Maximum resident set size \(kbytes\): (\d+)'),
            'CPU_percent':value(r'Percent of CPU this job got: (\d+)%'),
            'core_loop_reported_s':sum(float(v) for k,v in timings)/1000,
            'modal_sha256':sha(path/'output/kanaya_diagnostic/modal_native.bin'),
            'state_sha256':sha(path/'output/regression_state.bin')}


def pair(name,threads,legacy_second=False):
    cpu1=','.join(str(v) for v in range(0,8,2))
    cpu2=','.join(str(v) for v in range(8,16,2))
    # Each job stays on one half of the exposed physical cores. With one
    # thread the Intel runtime binds the process worker to the first core.
    if threads==1: cpu1='0';cpu2='8'
    start=time.perf_counter()
    with ThreadPoolExecutor(2) as pool:
        a=pool.submit(run_case,f'benchmark/{name}_a',BIN,steps=2000,threads=threads,cpus=cpu1,extra=NO_VTU)
        b=pool.submit(run_case,f'benchmark/{name}_b',BIN,steps=2000,threads=threads,cpus=cpu2,
                      extra=NO_VTU+(OLD if legacy_second else ''))
        records=[a.result(),b.result()]
    return {'name':name,'threads_per_job':threads,'jobs':2,'batch_wall_s':time.perf_counter()-start,
            'job_records':[metrics(r) for r in records]}


def main():
    check_gate()
    if (OUT/'benchmark_selection.json').exists(): raise FileExistsError('benchmark already completed')
    rows=[]
    physical=','.join(str(v) for v in range(0,16,2))
    for threads in [1,4,8,16,2]:
        name=f'benchmark/modern_t{threads}'
        if (OUT/name/'manifest.json').exists():
            record=json.loads((OUT/name/'manifest.json').read_text()); assert record['exit_code']==0
        else:
            cpus='0-15' if threads==16 else ','.join(str(v) for v in range(0,threads*2,2))
            record=run_case(name,BIN,steps=2000,threads=threads,cpus=cpus,extra=NO_VTU)
        rows.append(metrics(record))
    # schedule(static) reductions must preserve the same trajectory at every
    # tested thread count, not merely a matching rounded waveform.
    assert len({r['modal_sha256'] for r in rows})==1, 'thread configuration changed modal results'
    assert len({r['state_sha256'] for r in rows})==1, 'thread configuration changed full final state'
    for threads in [1,4]:
        name=f'benchmark/legacy_t{threads}'
        record=run_case(name,BIN,steps=2000,threads=threads,
                        cpus=','.join(str(v) for v in range(0,threads*2,2)),extra=NO_VTU+OLD)
        rows.append(metrics(record))
    pairs=[pair('pair_t1',1,True),pair('pair_t4',4,True)]
    flat=rows+[r for p in pairs for r in p['job_records']]
    pd.DataFrame(flat).to_csv(OUT/'benchmark_results.csv',index=False)
    serials=[]
    for threads in [1,4]:
        total=sum(r['elapsed_s'] for r in rows if r['name'] in
                  [f'benchmark/modern_t{threads}',f'benchmark/legacy_t{threads}'])
        serials.append({'name':f'serial_t{threads}','threads_per_job':threads,'jobs':1,'batch_wall_s':total})
    winner=min(pairs+serials,key=lambda p:p['batch_wall_s'])
    single_modern=min(rows[:5],key=lambda r:r['elapsed_s'])
    dump(OUT/'benchmark_selection.json',{'passed':True,'logical_CPUs':16,'exposed_physical_cores':8,
        'conditions':'same 20 modes / 1600 Pa / 2000 steps / full diagnostics; VTU disabled and numerically validated',
        'rows':rows,'pairs':pairs,'serial_pairs':serials,'winner':winner,
        'fastest_single_modern':single_modern,'thread_trajectories_bitwise_equal':True,
        'limits':'One sample per single configuration; includes model loading. This is a measured finite candidate set, not a universal optimality proof.'})
    print('Selected measured pair schedule:',winner,flush=True)


if __name__=='__main__':main()
