#!/usr/bin/env python3
"""Compare coupled continuation checkpoints to an uninterrupted 3D run.

Checkpoints before the first exchange and between exchanges retain the pending
clock and velocity reference. GoSPL rank count stays fixed; DES ranks may change.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import pickle
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tests/functional/mpi'))


def verify(work,frames,total):
    import h5py
    import numpy as np
    from run_gpu_mpi import compare_fields
    def snapshot(folder,frame):
        path=folder/('result.chkpt.%06d.vtkhdf'%frame)
        with h5py.File(path) as file:
            token=file['provenance'].attrs['gospl_generation']
            if isinstance(token,bytes):token=token.decode()
            clock=float(file['VTKHDF/grid/time'][()].item())/(365.2422*86400)
            pending=float(file['VTKHDF/grid/gospl coupling state'][1])
        directory=Path(str(path)+'.gospl')/token
        meta=json.loads((directory/'manifest.json').read_text());states=[]
        for rank,digest in enumerate(meta['files']):
            payload=(directory/('rank-%04d.pkl'%rank)).read_bytes()
            if hashlib.sha256(payload).hexdigest()!=digest:raise RuntimeError('Checkpoint checksum changed')
            state=pickle.loads(payload)
            if not np.isclose(float(state['data']['tNow'])+pending,clock,rtol=1e-12,atol=1e-12):
                raise RuntimeError('DES/GoSPL clocks differ in checkpoint')
            states.append(state)
        return states
    reference=snapshot(work/'continuous',total)
    checks=[]
    physical={'hGlobal','hLocal','cumED','cumEDLocal','vSed','vSedLocal','vSedF','Eb','FAG','headG','duriHG','Gsoil'}
    for start in frames:
        folder=work/('resume%d'%start)
        comparisons=[compare_fields(work/'continuous',folder,frame,kind)
                     for frame in range(start+1,total+1) for kind in ('save','chkpt')]
        current=snapshot(folder,total);vector_errors=[]
        if len(current)!=len(reference):raise RuntimeError('GoSPL rank count changed')
        for rank,(old,new) in enumerate(zip(reference,current)):
            for name in physical & old['vectors'].keys():
                a,b=old['vectors'][name],new['vectors'][name]
                if a is None or b is None:
                    if a is not b:raise RuntimeError('Missing GoSPL vector '+name)
                    continue
                x,y=a['values'],b['values']
                if x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():raise RuntimeError('Invalid GoSPL state '+name)
                error=float(np.max(abs(x-y))) if x.size else 0.
                relative=error/(float(np.max(abs(x))) or 1.) if x.size else 0.
                if relative>1e-8:raise RuntimeError('GoSPL continuation differs: '+name+' '+str(relative))
                vector_errors.append(dict(rank=rank,name=name,relative_max=relative,absolute_max=error))
        checks.append(dict(restart_frame=start,comparisons=comparisons,gospl_vectors=vector_errors,
                           all_bit_exact=all(x['bit_exact'] for x in comparisons)))
    report=dict(completed=True,all_passed=True,checks=checks,total_physical_comparisons=sum(len(c['comparisons']) for c in checks))
    (work/'continuation-results.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify-existing',type=Path)
    p.add_argument('--frames',type=int,nargs='+',default=[1,5,7])
    p.add_argument('--steps',type=int,default=12)
    for flag in ('exe','seed','extension','gospl-python','workdir'):p.add_argument('--'+flag,type=Path)
    p.add_argument('--launcher');p.add_argument('--gospl-launcher')
    p.add_argument('--ranks',type=int,default=3);p.add_argument('--restart-ranks',type=int)
    p.add_argument('--gospl-ranks',type=int,default=2)
    p.add_argument('--worker-env',action='append',default=[])
    p.add_argument('--des-env',action='append',default=[])
    p.add_argument('--remesh',action='store_true')
    p.add_argument('--time-mode',action='store_true')
    a=p.parse_args()
    if a.verify_existing:
        result=verify(a.verify_existing,a.frames,a.steps)
        print('Continuation PASS:',result['total_physical_comparisons'],'physical comparisons');return
    if any(getattr(a,name) is None for name in ('exe','seed','extension','gospl_python','workdir','launcher','gospl_launcher')):
        p.error('Provide executable, seed, launchers and environment paths')
    if any(frame<=0 or frame>=a.steps for frame in a.frames):p.error('Restart frames must lie inside the run')
    a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='continuation-',dir=str(a.workdir))).resolve()
    print('Artifacts:',work,flush=True)
    (work/'case.json').write_text(json.dumps(dict(arguments=vars(a),
        executable_sha256=hashlib.sha256(a.exe.read_bytes()).hexdigest(),
        coupling_sources={path.name:hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (ROOT/'gospl_driver').glob('*.py')}),default=str,indent=2)+'\n')
    for start in [None]+a.frames:
        folder=work/('continuous' if start is None else 'resume%d'%start);folder.mkdir()
        c=configparser.ConfigParser(inline_comment_prefixes=('#',';'));c.optionxform=str;c.read(a.seed/'model.cfg')
        c['sim'].update(max_steps=str(a.steps),output_step_interval='1',is_restarting='yes',
            restarting_from_modelname=str((a.seed if start is None else work/'continuous')/'result'),
            restarting_from_frame=str(start or 0))
        c['control'].update(gospl_coupling_frequency='3',gospl_coupling_mode='time' if a.time_mode else 'steps',
                            gospl_coupling_interval_in_yr='50')
        if a.remesh:c['mesh'].update(quality_check_step_interval='2',max_boundary_distortion='1e-9')
        with (folder/'model.cfg').open('w') as stream:c.write(stream)
        shutil.copy2(a.seed/'gospl.yml',folder/'gospl.yml')
        if start is None:shutil.copy2(a.seed/'gospl_mesh.npz',folder/'gospl_mesh.npz')
        ranks=a.ranks if start is None else (a.restart_ranks or a.ranks)
        command=[sys.executable,str(ROOT/'gospl_driver/run_coupled.py'),'--gospl-python',str(a.gospl_python.resolve()),
            '--interface',str(a.extension.resolve()/'cpp_interface/gospl_python_interface.py'),
            '--gospl-launcher',a.gospl_launcher+' -n '+str(a.gospl_ranks),
            '--workdir',str(folder),'--timeout','90','--run-timeout','180']
        for setting in a.worker_env:command+=['--worker-env',setting]
        command+=['--']+shlex.split(a.launcher)+['-n',str(ranks),str(a.exe.resolve()),'model.cfg']
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        for setting in a.des_env:
            key,value=setting.split('=',1);env[key]=value
        with (folder/'driver.log').open('w') as log:
            result=subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=240)
        if result.returncode:raise RuntimeError('Coupled continuation failed: '+str(folder))
        print(folder.name,'execution PASS',flush=True)
    result=verify(work,a.frames,a.steps)
    print('Continuation PASS:',result['total_physical_comparisons'],'physical comparisons',flush=True)


if __name__=='__main__':main()
