#!/usr/bin/env python3
"""Complete the 8-thread single-job / 16-logical-CPU two-job comparison."""
from concurrent.futures import ThreadPoolExecutor
import json
import time
import pandas as pd
from benchmark import BIN, NO_VTU, OLD, metrics, check_gate
from run_matrix import OUT, run_case, dump


def main():
    check_gate()
    initial=json.loads((OUT/'benchmark_selection.json').read_text())
    if (OUT/'benchmark_final_selection.json').exists():
        raise FileExistsError('final benchmark already completed')
    record=run_case('benchmark/legacy_t8',BIN,steps=2000,threads=8,
                    cpus='0,2,4,6,8,10,12,14',extra=NO_VTU+OLD)
    legacy8=metrics(record)
    start=time.perf_counter()
    # Both hardware threads of every exposed core are exercised, with one
    # sibling assigned to each job. This measures, rather than assumes, the
    # value of SMT for throughput. No job exceeds the 16 available CPUs.
    with ThreadPoolExecutor(2) as pool:
        a=pool.submit(run_case,'benchmark/pair_t8_a',BIN,steps=2000,threads=8,
                      cpus='0,2,4,6,8,10,12,14',extra=NO_VTU)
        b=pool.submit(run_case,'benchmark/pair_t8_b',BIN,steps=2000,threads=8,
                      cpus='1,3,5,7,9,11,13,15',extra=NO_VTU+OLD)
        records=[a.result(),b.result()]
    pair8={'name':'pair_t8','threads_per_job':8,'jobs':2,
           'batch_wall_s':time.perf_counter()-start,
           'job_records':[metrics(r) for r in records]}
    serial8={'name':'serial_t8','threads_per_job':8,'jobs':1,
             'batch_wall_s':legacy8['elapsed_s']+initial['fastest_single_modern']['elapsed_s']}
    # Use the actual t8 modern result, even if it was not the fastest sample.
    serial8['batch_wall_s']=legacy8['elapsed_s']+next(
        r['elapsed_s'] for r in initial['rows'] if r['name']=='benchmark/modern_t8')
    candidates=initial['pairs']+initial['serial_pairs']+[pair8,serial8]
    winner=min(candidates,key=lambda r:r['batch_wall_s'])
    rows=initial['rows']+[legacy8]+[r for p in initial['pairs']+[pair8] for r in p['job_records']]
    for prefix in ['modern','legacy']:
        selected=[r for r in rows if '/'+prefix+'_t' in r['name']]
        assert len({r['modal_sha256'] for r in selected})==1
        assert len({r['state_sha256'] for r in selected})==1
    pd.DataFrame(rows).to_csv(OUT/'benchmark_final_results.csv',index=False)
    dump(OUT/'benchmark_final_selection.json',dict(passed=True,candidates=candidates,winner=winner,
        fastest_single_modern=initial['fastest_single_modern'],rows=rows,
        thread_trajectories_bitwise_equal=True,
        limits=initial['limits'],logical_CPUs=16,exposed_physical_cores=8))
    print('Final measured schedule:',winner,flush=True)


if __name__=='__main__':main()
