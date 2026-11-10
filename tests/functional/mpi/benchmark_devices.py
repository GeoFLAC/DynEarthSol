#!/usr/bin/env python3
"""Measure checked OpenMP, CPU MPI and GPU MPI runs on one allocated host.

Variants are a JSON list of {name, kind, workers, executable, launcher, env}.
kind is openmp, mpi or gpu; MPI launchers accept '-n N executable args'.
Each variant supplies its own runtime environment, so distinct MPI ABIs never
share a library path. Warmups and every measured run are physically compared.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shlex
import statistics
import sys
import tempfile

from launch_sge import allocation
from run_mpi import configure, execute, metrics
from run_gpu_mpi import compare_fields, device_records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variants',type=Path,required=True)
    p.add_argument('--resolutions',type=int,nargs='+',default=[500,125])
    p.add_argument('--steps',type=int,default=1000)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--timeout',type=int,default=900)
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args()
    variants=json.loads(a.variants.read_text())
    cpus=sorted(os.sched_getaffinity(0))
    context=allocation() if os.environ.get('JOB_ID') else None
    if context and len(context['host_slots'])!=1:p.error('Use a single allocated host')
    budget=min(len(cpus),context['nslots'] if context else len(cpus))
    if a.repeats<3 or a.steps<1 or a.timeout<1:p.error('Use at least three repeats and positive steps/timeout')
    if not variants or variants[0]['kind']!='openmp' or variants[0]['workers']!=1:
        p.error('The first variant must be the one-thread CPU reference')
    names=[v['name'] for v in variants]
    if len(names)!=len(set(names)):p.error('Variant names must be distinct')
    for v in variants:
        if not re.fullmatch('[A-Za-z0-9_-]+',v['name']):p.error('Use simple variant names')
        if v['kind'] not in ('openmp','mpi','gpu') or not 1<=v['workers']<=budget:
            p.error('Invalid variant kind or CPU budget')
        v['executable']=str(Path(v['executable']).resolve())
        v['sha256']=hashlib.sha256(Path(v['executable']).read_bytes()).hexdigest()
        if v['kind']!='openmp' and not v.get('launcher'):p.error('MPI variants require explicit launchers')
    a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='devices-',dir=str(a.workdir))).resolve()
    print('Artifacts:',work,flush=True)
    report=dict(completed=False,all_passed=False,allocation=context,cpus_available=cpus,
                cpu_budget=budget,variants=variants,steps=a.steps,repeats=a.repeats,
                timing_note='Simulation includes initiation, output and management. MPI operation and communication timers overlap.',
                measurements=[],runs=[])
    def save():(work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    save()
    try:
        for resolution in a.resolutions:
            cfg=configure(2,False,resolution,steps=a.steps,output_interval=a.steps)
            # Triangle compiled with different host compilers can choose distinct
            # initial meshes. Every measured variant restarts the same CPU frame
            # zero, so geometry, fields and workload are identical at entry.
            seed=work/('r{}-initial-checkpoint'.format(resolution));seed.mkdir()
            initial=configure(2,False,resolution,steps=1,output_interval=a.steps)
            with (seed/'model.cfg').open('w') as stream:initial.write(stream)
            first=variants[0]
            initial_env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1')
            initial_env.update(first.get('env',{}))
            execute([first['executable'],'model.cfg'],seed,initial_env,timeout=a.timeout)
            cfg['sim'].update(is_restarting='yes',restarting_from_modelname=str(seed/'result'),restarting_from_frame='0')
            report.setdefault('initial_checkpoints',[]).append(dict(resolution=resolution,source=str(seed),
                                                                    sha256=hashlib.sha256((seed/'result.chkpt.000000.vtkhdf').read_bytes()).hexdigest()))
            reference=None
            samples={v['name']:[] for v in variants}
            def run(v,label):
                directory=work/('r{}-{}-{}'.format(resolution,v['name'],label));directory.mkdir()
                with (directory/'model.cfg').open('w') as stream:cfg.write(stream)
                n=v['workers'];kind=v['kind'];selected=cpus[:n]
                env=dict(os.environ,OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',
                         OMP_NUM_THREADS=str(n if kind=='openmp' else 1))
                env.update(v.get('env',{}))
                if kind=='openmp':env.update(OMP_PROC_BIND='close',OMP_PLACES=','.join('{'+str(c)+'}' for c in selected))
                else:env.update(OMP_PROC_BIND='false')
                binding=[sys.executable,str(Path(__file__).with_name('pin_cpu.py').resolve()),
                         '--cpus',','.join(map(str,selected)),'--record',str(directory),'--exe',v['executable']]
                if kind!='openmp':binding+=['--mpi']
                command=binding if kind=='openmp' else shlex.split(v['launcher'])+['-n',str(n)]+binding
                (directory/'command.json').write_text(json.dumps(command,indent=2)+'\n')
                log,wall=execute(command,directory,env,timeout=a.timeout)
                if 'Ending simulation.' not in log:raise RuntimeError('Incomplete simulation')
                masks=[json.loads(path.read_text()) for path in directory.glob('binding-*.json')]
                if len(masks)!=(1 if kind=='openmp' else n):raise RuntimeError('Missing affinity records')
                if any(row['cpus']!=(selected if kind=='openmp' else [selected[row['rank']]]) for row in masks):
                    raise RuntimeError('Unexpected CPU affinity')
                timing=metrics(log,wall)
                if 'simulation_s' not in timing:raise RuntimeError('Missing simulation timer')
                timing['phases_s']={key.lower():3600*int(h)+60*int(m)+float(s) for key,h,m,s in
                                    re.findall(r'(Execute|Initiate|Compute|Remesh|Output)\s*:\s*(\d+):(\d+):([\d.]+)',log)}
                checks=[]
                if reference:
                    # Restart writes the initial save frame, but does not emit
                    # another copy of its input checkpoint at frame zero.
                    checks=[compare_fields(reference,directory,frame,kind=output)
                            for output,frames in (('save',(0,1)),('chkpt',(1,))) for frame in frames]
                devices=device_records(log,n) if kind=='gpu' else []
                row=dict(resolution=resolution,name=v['name'],sample=label,timing=timing,
                         bindings=masks,devices=devices,comparisons=checks,directory=str(directory))
                report['runs'].append(row);save()
                return directory,timing
            for v in variants:
                current,_=run(v,'warmup')
                if reference is None:reference=current
            randomizer=random.Random(20261109+resolution)
            for repeat in range(a.repeats):
                order=list(variants);randomizer.shuffle(order)
                for v in order:
                    _,timing=run(v,'sample'+str(repeat));samples[v['name']].append(timing)
                print('resolution',resolution,'measured round',repeat+1,'PASS',flush=True)
            baseline=statistics.median(t['simulation_s'] for t in samples[variants[0]['name']])
            for v in variants:
                data=samples[v['name']];median=statistics.median(t['simulation_s'] for t in data)
                row=dict(resolution=resolution,name=v['name'],kind=v['kind'],workers=v['workers'],
                         median_seconds=median,speedup_from_cpu1=baseline/median,samples=data)
                report['measurements'].append(row)
                print(v['name'],'median',median,'s speedup',row['speedup_from_cpu1'],flush=True)
            save()
    except Exception as error:
        report['failure']=str(error);save();raise
    report.update(completed=True,all_passed=True);save()
    print('Results:',work/'results.json',flush=True)


if __name__=='__main__':main()
