#!/usr/bin/env python3
"""Compare the preceding MPI backend with blocking/overlapped updated kernels.

Pinned CPU IDs, warmups and randomized interleaved repeats. Run without competing
builds/simulations. Reports include per-rank profiles and physical-output checks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import statistics
import sys
import tempfile
from run_mpi import configure, execute, metrics, compare_physical_hdf5


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before',type=Path,required=True)
    p.add_argument('--after',type=Path,required=True)
    p.add_argument('--serial',type=Path,required=True)
    p.add_argument('--workers',nargs='+',type=int,default=[3,4])
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args();cpus=sorted(os.sched_getaffinity(0))
    if a.repeats<3 or min(a.workers)<1 or max(a.workers)>len(cpus):p.error('Need >=3 repeats and available positive CPU counts')
    a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='overhead-',dir=a.workdir));print('Artifacts:',work,flush=True)
    report=dict(completed=False,cpus=cpus,repeats=a.repeats,measurements=[],executables={})
    for key in ['before','after','serial']:
        exe=getattr(a,key).resolve();setattr(a,key,exe)
        report['executables'][key]=dict(path=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    for label,resolution,moving,steps,interval in [('small-fixed',500,False,1000,1000),('large-fixed',125,False,1000,1000),('adaptive',1000,True,196,49)]:
        cfg=configure(2,moving,resolution,steps=steps,output_interval=interval)
        reference=work/(label+'-parent');reference.mkdir()
        with (reference/'model.cfg').open('w') as f:cfg.write(f)
        execute([str(a.serial),'model.cfg'],reference,dict(os.environ,OMP_NUM_THREADS='1'))
        variants=[(mode,n) for n in a.workers for mode in ['before','blocking','overlap']]
        samples={v:[] for v in variants}
        def run(mode,n,tag):
            directory=work/f'{label}-{mode}{n}-{tag}';directory.mkdir()
            with (directory/'model.cfg').open('w') as f:cfg.write(f)
            executable=a.before if mode=='before' else a.after
            env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='false',
                     DES_MPI_PARTITIONER='rcb',DES_MPI_OVERLAP='1' if mode=='overlap' else '0')
            command=['mpiexec','--bind-to','none','-n',str(n),sys.executable,str(Path(__file__).with_name('pin_cpu.py').resolve()),
                     '--mpi','--cpus',','.join(map(str,cpus[:n])),'--record',str(directory),'--exe',str(executable)]
            log,elapsed=execute(command,directory,env)
            masks=[json.loads(p.read_text()) for p in directory.glob('binding-*.json')]
            if len(masks)!=n or any(row['cpus']!=[cpus[row['rank']]] for row in masks):raise RuntimeError('Incorrect rank affinity')
            frames=sorted(int(p.name.split('.')[-2]) for p in reference.glob('result.save.*.vtkhdf'))
            for frame in frames:compare_physical_hdf5(reference,directory,frame)
            row=metrics(log,elapsed);row.update(command=command,bindings=masks,bit_exact=True,directory=str(directory))
            return row
        for mode,n in variants:run(mode,n,'warmup')
        rng=random.Random(20261106+resolution)
        for repeat in range(a.repeats):
            order=list(variants);rng.shuffle(order)
            for mode,n in order:samples[(mode,n)].append(run(mode,n,f'sample{repeat}'))
            print(label,'round',repeat+1,'complete',flush=True)
        for (mode,n),rows in samples.items():
            item=dict(case=label,mode=mode,workers=n,median_seconds=statistics.median(r['simulation_s'] for r in rows),samples=rows)
            report['measurements'].append(item);print(label,mode,n,item['median_seconds'],flush=True)
        (work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    report['completed']=True;(work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Results:',work/'results.json',flush=True)


if __name__=='__main__':main()
