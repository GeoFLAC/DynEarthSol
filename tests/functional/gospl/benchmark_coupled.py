#!/usr/bin/env python3
"""One controlled 3D coupled workload: CPU/OpenMP, CPU MPI and GPU MPI.

All variants restart the same initial checkpoint. Launcher/device startup,
coupling calls and output are timed separately; every measured output is checked.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tests/functional/mpi'))
from run_mpi import metrics


def timing_parts(text):
    result = {}
    for name, key in (('Execute', 'execute_s'), ('Initiate', 'initialize_s'),
                      ('Compute', 'compute_s'), ('Output', 'output_s')):
        match = re.search(r'\b'+name+r'\s*:\s*(\d+):(\d+):([\d.]+)', text)
        if match is None:
            raise RuntimeError('Missing timing category: '+name)
        result[key] = 3600*int(match[1])+60*int(match[2])+float(match[3])
    return result


def summarize(report):
    result = []
    for name in ('omp1', 'omp4', 'mpi3', 'gpu1', 'gpu3'):
        rows = [r for r in report['runs'] if r['label'].startswith('round') and r['label'].endswith('-'+name)]
        row = dict(name=name, kind=rows[0]['kind'], ranks=rows[0]['ranks'], threads=rows[0]['threads'],
                   median_wall_seconds=statistics.median(r['wall_seconds'] for r in rows),
                   min_wall_seconds=min(r['wall_seconds'] for r in rows),
                   max_wall_seconds=max(r['wall_seconds'] for r in rows), sample_labels=[r['label'] for r in rows])
        for key in ('execute_s', 'initialize_s', 'compute_s', 'output_s'):
            row['median_'+key] = statistics.median(r['timing'][key] for r in rows)
        row['median_worker_erosion_s'] = statistics.median(r['worker_seconds']['run_and_get_erosion'] for r in rows)
        result.append(row)
    report['measurements'] = result


def physical_cpus():
    result=[];seen=set()
    for cpu in sorted(os.sched_getaffinity(0)):
        path=Path('/sys/devices/system/cpu/cpu%d/topology'%cpu)
        key=((path/'physical_package_id').read_text(),(path/'core_id').read_text())
        if key not in seen:seen.add(key);result.append(cpu)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('cpu','gpu','seed-config','extension','gospl-python','workdir'):
        p.add_argument('--'+flag,type=Path,required=True)
    for flag in ('cpu-launcher','gpu-launcher','gospl-launcher'):
        p.add_argument('--'+flag,required=True)
    p.add_argument('--worker-env',action='append',default=[])
    p.add_argument('--gpu-env',action='append',default=[])
    p.add_argument('--resolution',type=float,default=1250)
    p.add_argument('--gospl-resolution',type=float,default=5000)
    p.add_argument('--steps',type=int,default=100)
    p.add_argument('--coupling-every',type=int,default=20)
    p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args()
    if a.steps%a.coupling_every or a.repeats<1:p.error('Use complete coupling intervals and positive repeats')
    cpus=physical_cpus()
    if len(cpus)<5:p.error('At least five physical CPUs are required')
    a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='coupled-cost-',dir=str(a.workdir))).resolve()
    print('Artifacts:',work,flush=True)
    report=dict(completed=False,all_passed=False,resolution=a.resolution,gospl_resolution=a.gospl_resolution,
                steps=a.steps,repeats=a.repeats,des_cpus=cpus[:4],worker_cpus=[cpus[4]],runs=[],measurements=[],
                executable_sha256={k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest() for k in ('cpu','gpu')})
    target=work/'results.json'
    def save():target.write_text(json.dumps(report,indent=2)+'\n')
    save()
    variants=[('omp1','cpu',1,1),('omp4','cpu',1,4),('mpi3','cpu',3,1),('gpu1','gpu',1,1),('gpu3','gpu',3,1)]
    seed=work/'seed'
    def run(label,kind,ranks,threads,seed_only=False):
        folder=work/label;folder.mkdir()
        c=configparser.ConfigParser(inline_comment_prefixes=('#',';'));c.optionxform=str;c.read(a.seed_config)
        c['sim'].update(modelname='result',max_steps=str(1 if seed_only else a.steps),
            output_step_interval=str(1 if seed_only else a.steps),output_time_interval_in_yr='1e20',
            max_time_in_yr='1e20',checkpoint_frame_interval='1',is_restarting='no' if seed_only else 'yes',
            restarting_from_modelname=str(seed/'result'),restarting_from_frame='0')
        c['mesh'].update(resolution=str(a.resolution),quality_check_step_interval='1000')
        c['control'].update(gospl_coupling_frequency=str(a.coupling_every),gospl_coupling_mode='steps',
                            gospl_mesh_resolution=str(a.gospl_resolution))
        with (folder/'model.cfg').open('w') as f:c.write(f)
        shutil.copy2(a.seed_config.parent/'gospl.yml',folder/'gospl.yml')
        if not seed_only:shutil.copy2(seed/'gospl_mesh.npz',folder/'gospl_mesh.npz')
        env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OMP_WAIT_POLICY='PASSIVE',OPENBLAS_NUM_THREADS='1')
        if kind=='gpu':
            for setting in a.gpu_env:
                key,value=setting.split('=',1);env[key]=value
        prefix=a.cpu_launcher if kind=='cpu' else a.gpu_launcher
        command=[sys.executable,str(ROOT/'gospl_driver/run_coupled.py'),
            '--gospl-python',str(a.gospl_python.resolve()),
            '--interface',str(a.extension.resolve()/'cpp_interface/gospl_python_interface.py'),
            '--gospl-launcher',a.gospl_launcher+' -n 1 taskset -c '+str(cpus[4]),
            '--workdir',str(folder),'--timeout','180','--run-timeout','600']
        for setting in a.worker_env:command+=['--worker-env',setting]
        command+=['--']+shlex.split(prefix)+['-n',str(ranks),str(a.gospl_python.resolve()),str(Path(__file__).resolve()),
            '--rank',','.join(map(str,cpus[:4])),str(threads),str(getattr(a,kind).resolve()),'model.cfg']
        start=time.perf_counter()
        with (folder/'driver.log').open('w') as log:
            result=subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=660,
                                  preexec_fn=lambda:os.sched_setaffinity(0,cpus[:5]))
        wall=time.perf_counter()-start
        if result.returncode:raise RuntimeError('Coupled benchmark failed: '+label)
        text=(folder/'des.log').read_text()
        if 'Ending simulation.' not in text:raise RuntimeError('Incomplete DES run: '+label)
        timing=metrics(text,wall)
        timing.update(timing_parts(text))
        calls=[json.loads(line) for line in (folder/'worker-calls.jsonl').read_text().splitlines()]
        check=dict(label=label,kind=kind,ranks=ranks,threads=threads,wall_seconds=wall,timing=timing,
                   worker_seconds={op:sum(r['seconds'] for r in calls if r['operation']==op)
                                   for op in {r['operation'] for r in calls}})
        if not seed_only:
            reference=work/'warmup-omp1'
            with (folder/'comparison.log').open('w') as log:
                verify=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--compare',str(reference),str(folder)],
                                      stdout=log,stderr=subprocess.STDOUT,timeout=180)
            if verify.returncode:raise RuntimeError('Physical comparison failed: '+label)
            check['comparison']=json.loads((folder/'comparison.json').read_text())
        report['runs'].append(check);save();print(label,'PASS',flush=True)
    try:
        run('seed','cpu',1,1,True)
        report['initial_checkpoint_sha256']=hashlib.sha256((seed/'result.chkpt.000000.vtkhdf').read_bytes()).hexdigest()
        for name,kind,ranks,threads in variants:run('warmup-'+name,kind,ranks,threads)
        for repeat in range(a.repeats):
            order=list(variants);random.Random(117+repeat).shuffle(order)
            for name,kind,ranks,threads in order:run('round%d-%s'%(repeat,name),kind,ranks,threads)
        summarize(report)
        report.update(completed=True,all_passed=True);save()
    except BaseException as error:
        report['failure']=str(error);save();raise


if __name__=='__main__':
    if sys.argv[1:2]==['--summarize-existing']:
        work=Path(sys.argv[2]);target=work/'results.json'
        report=json.loads(target.read_text())
        if not report['completed'] or not report['all_passed']:raise RuntimeError('Incomplete benchmark')
        for row in report['runs']:
            row['timing'].update(timing_parts((work/row['label']/'des.log').read_text()))
        summarize(report);target.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report['measurements'],indent=2))
    elif sys.argv[1:2]==['--rank']:
        cpus=list(map(int,sys.argv[2].split(',')));threads=int(sys.argv[3])
        rank=int(os.environ.get('OMPI_COMM_WORLD_LOCAL_RANK','0'))
        os.sched_setaffinity(0,cpus[:threads] if threads>1 else [cpus[rank]])
        os.execv(sys.argv[4],sys.argv[4:])
    elif sys.argv[1:2]==['--compare']:
        from run_gpu_mpi import compare_fields
        reference,current=map(Path,sys.argv[2:4])
        checks=[compare_fields(reference,current,frame,kind) for kind,frames in [('save',[0,1]),('chkpt',[1])] for frame in frames]
        (current/'comparison.json').write_text(json.dumps(checks,indent=2)+'\n')
    else:main()
