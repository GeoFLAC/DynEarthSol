#!/usr/bin/env python3
"""Compare equally sized OpenMP/MPI CPU budgets with warmups and interleaved runs.

Open MPI and Linux CPU affinity are required for this pinned local comparison.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import sys
import re
import tempfile

from run_mpi import configure,execute,metrics,compare_physical_hdf5


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--serial',type=Path,required=True)
    p.add_argument('--mpi',type=Path,required=True)
    p.add_argument('--workers',type=int,nargs='+',default=[1,2,3,4,6,8])
    p.add_argument('--resolutions',type=int,nargs='+',default=[500,125])
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--workdir',type=Path)
    a=p.parse_args()
    available=sorted(os.sched_getaffinity(0))
    if (not a.workers or min(a.workers)<1 or max(a.workers)>len(available)
            or len(a.workers)!=len(set(a.workers)) or a.repeats<3):
        p.error('Need available CPU IDs, positive worker counts, and at least three repeats')
    if a.workdir:a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='scaling-',dir=a.workdir));print('Artifacts:',work,flush=True)
    serial,mpi=a.serial.resolve(),a.mpi.resolve()
    executables={key:dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                 for key,path in [('serial',serial),('mpi',mpi)]}
    output=dict(cpus_available=available,repeats=a.repeats,executables=executables,
                comparison='All physical HDF5 datasets exact; walltime excluded',measurements=[])
    for resolution in a.resolutions:
        cfg=configure(2,False,resolution,steps=1000,output_interval=1000)
        variants=[(mode,n) for n in a.workers for mode in ('openmp','mpi')]
        samples={v:[] for v in variants};reference=None
        def run(mode,n,label):
            directory=work/f'r{resolution}-{mode}{n}-{label}';directory.mkdir()
            with (directory/'model.cfg').open('w') as f:cfg.write(f)
            cpus=available[:n];cpu_list=','.join(map(str,cpus))
            env=dict(os.environ,OMP_DYNAMIC='FALSE',OMP_NUM_THREADS=str(n if mode=='openmp' else 1))
            binding=[sys.executable,str(Path(__file__).with_name('pin_cpu.py').resolve()),
                     '--cpus',cpu_list,'--record',str(directory),'--exe',str(serial if mode=='openmp' else mpi)]
            if mode=='openmp':
                env.update(OMP_PROC_BIND='close',OMP_PLACES=','.join('{'+str(c)+'}' for c in cpus))
                command=binding
            else:
                env.update(OMP_PROC_BIND='false',OMP_PLACES='cores')
                command=['mpiexec','--bind-to','none','-n',str(n)]+binding+['--mpi']
            log,elapsed=execute(command,directory,env)
            actual_threads=int(re.search(r'omp_threads=(\d+)',log)[1])
            if actual_threads!=(n if mode=='openmp' else 1):raise RuntimeError('Unexpected OpenMP thread count')
            masks=[json.loads(path.read_text()) for path in directory.glob('binding-*.json')]
            if len(masks)!=(n if mode=='mpi' else 1):raise RuntimeError('Missing CPU binding records')
            if mode=='mpi' and any(row['cpus']!=[cpus[row['rank']]] for row in masks):raise RuntimeError('Incorrect MPI CPU binding')
            timing=metrics(log,elapsed);timing.update(command=command,cpus=cpus,bindings=masks)
            return directory,timing
        # Warm each worker mode once. The one-thread parent provides the reference.
        for mode,n in variants:
            current,_=run(mode,n,'warmup')
            if reference is None:reference=current
            compare_physical_hdf5(reference,current,1)
        randomizer=random.Random(20261105+resolution)
        for repeat in range(a.repeats):
            order=list(variants);randomizer.shuffle(order)
            for mode,n in order:
                current,timing=run(mode,n,f'sample{repeat}')
                if repeat==0:compare_physical_hdf5(reference,current,1)
                samples[(mode,n)].append(timing)
            print(f'resolution={resolution}: completed measured round {repeat+1}/{a.repeats}',flush=True)
        for (mode,n),data in samples.items():
            record=dict(resolution=resolution,mode=mode,workers=n,
                        median_seconds=statistics.median(r['simulation_s'] for r in data),samples=data)
            output['measurements'].append(record)
            print(f"resolution={resolution} {mode} workers={n}: {record['median_seconds']:.6f} s",flush=True)
        (work/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Results:',work/'results.json',flush=True)


if __name__=='__main__':main()
