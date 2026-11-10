#!/usr/bin/env python3
"""Bounded real coupling checks with independent DES and GoSPL MPI worlds."""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import h5py
import numpy as np

from run_coupling_smoke import audit

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tests/functional/mpi'))
from run_gpu_mpi import compare_fields, device_records
from run_mpi import metrics


def check_protocol_failure(a,work):
    folder=work/'protocol-error';folder.mkdir()
    command=[sys.executable,str(ROOT/'gospl_driver/run_coupled.py'),
             '--gospl-python',str(a.gospl_python.resolve()),
             '--interface',str(a.extension.resolve()/'cpp_interface/gospl_python_interface.py'),
             '--gospl-launcher',a.gospl_launcher+' -n 2',
             '--workdir',str(folder),'--timeout','30','--run-timeout','10']
    for value in a.worker_env:command+=['--worker-env',value]
    command+=['--',str(a.gospl_python.resolve()),'-c',
              "from des_gospl_client import _rpc; _rpc('invalid_operation')"]
    with (folder/'driver.log').open('w') as log:
        result=subprocess.run(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=60)
    record=json.loads((folder/'coupling-launch.json').read_text())
    if result.returncode==0 or record['completed'] or record['des_exit'] in (None,0) or record['worker_exit']!=53:
        raise RuntimeError('Worker protocol failure did not stop and reap both process groups')
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu', type=Path, required=True)
    p.add_argument('--gpu', type=Path, required=True)
    p.add_argument('--seed', type=Path, required=True, help='Previous bounded CPU fixture directory')
    p.add_argument('--extension', type=Path, required=True)
    p.add_argument('--gospl-python', type=Path, required=True)
    p.add_argument('--cpu-launcher', required=True)
    p.add_argument('--gpu-launcher', required=True)
    p.add_argument('--gospl-launcher', required=True)
    p.add_argument('--worker-env', action='append', default=[])
    p.add_argument('--gpu-env', action='append', default=[])
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--resume', type=Path, help='Reuse completed fixtures and append only selected cases')
    p.add_argument('--cases', nargs='+')
    a = p.parse_args()
    a.workdir.mkdir(parents=True, exist_ok=True)
    work = a.resume.resolve() if a.resume else Path(tempfile.mkdtemp(prefix='mpi-coupling-', dir=str(a.workdir))).resolve()
    print('Artifacts:', work, flush=True)
    seed = a.seed.resolve()
    report = dict(completed=False, all_passed=False, scope='DES coupling boundary; GoSPL algorithms unchanged',
                  checks=[], comparisons=[], executables={kind:dict(path=str(exe.resolve()),
                  sha256=hashlib.sha256(exe.read_bytes()).hexdigest()) for kind,exe in [('cpu',a.cpu),('gpu',a.gpu)]})
    target = work / 'results.json'
    if a.resume:
        previous=json.loads(target.read_text())
        if previous['executables']!=report['executables']:p.error('Resume requires the same executable hashes')
        report=previous
        report.update(completed=False,all_passed=False,comparisons=[])
        report.pop('failure',None)
    def save(): target.write_text(json.dumps(report, indent=2) + '\n')
    save()
    variants = [('cpu-d1-g1','cpu',1,1,False), ('cpu-d3-g1','cpu',3,1,False),
                ('cpu-d3-g2','cpu',3,2,False), ('gpu-d2-g1','gpu',2,1,False),
                ('gpu-d2-g2','gpu',2,2,False), ('remesh-d3-g2','cpu',3,2,True),
                ('gpu-remesh-d2-g2','gpu',2,2,True), ('time-d3-g2','cpu',3,2,False)]
    if a.cases:
        if set(a.cases)-{v[0] for v in variants}:p.error('Unknown coupling case')
        variants=[v for v in variants if v[0] in a.cases]
    try:
        for name,kind,des_ranks,gospl_ranks,remesh in variants:
            folder = work / name
            folder.mkdir()
            cfg = configparser.ConfigParser(inline_comment_prefixes=('#',';'))
            cfg.optionxform = str
            cfg.read(seed / 'model.cfg')
            cfg['sim'].update(is_restarting='yes', restarting_from_modelname=str(seed/'result'),
                              restarting_from_frame='0')
            if remesh:
                cfg['mesh'].update(quality_check_step_interval='2', max_boundary_distortion='1e-9')
            if name=='time-d3-g2':
                cfg['control'].update(gospl_coupling_mode='time',gospl_coupling_interval_in_yr='30')
            with (folder/'model.cfg').open('w') as stream: cfg.write(stream)
            for file in ('gospl_mesh.npz','gospl.yml'):
                shutil.copy2(seed/file,folder/file)
            mpi = a.cpu_launcher if kind == 'cpu' else a.gpu_launcher
            des = shlex.split(mpi) + ['-n',str(des_ranks),str(getattr(a,kind).resolve()),'model.cfg']
            command = [sys.executable,str(ROOT/'gospl_driver/run_coupled.py'),
                '--gospl-python',str(a.gospl_python.resolve()),
                '--interface',str(ROOT/'tests/functional/gospl/audit_interface/gospl_python_interface.py'),
                '--gospl-launcher',a.gospl_launcher+' -n '+str(gospl_ranks),
                '--workdir',str(folder),'--audit-dir',str(folder/'audit'),'--timeout','120','--run-timeout','180']
            for value in a.worker_env: command += ['--worker-env',value]
            command += ['--'] + des
            env = dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
                       GOSPL_REAL_INTERFACE=str(a.extension.resolve()/'cpp_interface/gospl_python_interface.py'))
            if kind == 'gpu':
                for setting in a.gpu_env:
                    key,value=setting.split('=',1);env[key]=value
            start = time.perf_counter()
            with (folder/'driver.log').open('w') as log:
                result = subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=240)
            if result.returncode:raise RuntimeError('Coupled run failed: '+name)
            text = (folder/'des.log').read_text()
            if 'Ending simulation.' not in text:raise RuntimeError('DES did not finish: '+name)
            check = audit(folder,gospl_ranks,verify_seed=True)
            timing = metrics(text,time.perf_counter()-start)
            check.update(name=name,des_ranks=des_ranks,gospl_ranks=gospl_ranks,timing=timing)
            if kind == 'gpu':check['devices']=device_records(text,des_ranks)
            if remesh:
                if not timing['remeshes']:raise RuntimeError('Remesh coupling fixture did not remesh')
            elif des_ranks>1 or kind=='gpu':
                if timing['kernels']['Surface']['calls']!=4:raise RuntimeError('Idle coupling steps performed surface communication')
                if timing['kernels']['State']['calls']!=1:raise RuntimeError('Coupling unnecessarily redistributed full state')
            report['checks'].append(check);save();print(name,'coupling PASS',flush=True)
        # Compare DES rank counts and CPU/GPU while keeping the GoSPL rank count
        # fixed. GoSPL partition-dependent numerical equivalence is out of scope.
        pairs=[('cpu-d1-g1','cpu-d3-g1'),('cpu-d1-g1','gpu-d2-g1'),('cpu-d3-g2','gpu-d2-g2')]
        if (work/'time-d3-g2').exists():pairs.append(('cpu-d3-g2','time-d3-g2'))
        for left,right in pairs:
            comparisons=[compare_fields(work/left,work/right,frame,kind)
                         for kind,frames in [('save',range(5)),('chkpt',range(1,5))] for frame in frames]  # restart omits checkpoint zero
            report['comparisons'].append(dict(reference=left,current=right,fields=comparisons))
        if 'protocol_failure_propagation' not in report:
            report['protocol_failure_propagation']=check_protocol_failure(a,work)
        report.update(completed=True,all_passed=True);save()
    except BaseException as error:
        report['failure']=str(error);save();raise
    print('Results:',target,flush=True)


if __name__ == '__main__':main()
