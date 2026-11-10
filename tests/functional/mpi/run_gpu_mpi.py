#!/usr/bin/env python3
"""Check complete GPU MPI timesteps against CPU outputs and one-GPU outputs.

Run through qsub on godo-s. Other GPU hosts may supply their own launcher. Full
physical HDF5 arrays are checked, including vectors/tensors and checkpoints;
integer data and topology must be exact, floating data finite within the normal
DES 1e-8 relative-field criterion. Timings here are not scaling measurements.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import sys
import tempfile

from run_mpi import configure, execute, metrics

ROOT = Path(__file__).resolve().parents[3]


def compare_fields(reference, current, frame, kind='save', rtol=1e-8, prefix='result'):
    import h5py
    import numpy as np
    name = '{}.{}.{:06d}.vtkhdf'.format(prefix,kind,frame)
    def arrays(path):
        found={}
        with h5py.File(str(path),'r') as file:
            def visit(name,obj):
                if isinstance(obj,h5py.Dataset) and name.split('/')[-1]!='walltime_sec':
                    found[name]=obj[()]
            file['VTKHDF'].visititems(visit)
        return found
    old,new=arrays(reference/name),arrays(current/name)
    if not old or old.keys()!=new.keys():raise RuntimeError('Physical dataset sets differ: '+name)
    differences=[]
    for field in old:
        a,b=np.asarray(old[field]),np.asarray(new[field])
        if a.shape!=b.shape or a.dtype!=b.dtype:raise RuntimeError('Shape/dtype differs: '+name+' '+field)
        exact=a.tobytes()==b.tobytes()
        if a.dtype.kind in 'fc':
            if not np.isfinite(a).all() or not np.isfinite(b).all():raise RuntimeError('Nonfinite field: '+field)
            maximum=float(np.max(np.abs(a-b))) if a.size else 0.0
            scale=float(np.max(np.abs(a))) if a.size else 0.0
            relative=maximum/(scale if scale else 1.0)
            if relative>rtol:raise RuntimeError('{} {} relative error {} > {}'.format(name,field,relative,rtol))
            if not exact:differences.append(dict(field=field,relative_max=relative,absolute_max=maximum))
        elif not exact:raise RuntimeError('Integer/topology field differs: '+name+' '+field)
    return dict(kind=kind,frame=frame,datasets=len(old),bit_exact=not differences,
                floating_differences=differences,rtol=rtol)


def device_records(log,ranks):
    rows=[dict(re.findall(r'(\w+)=(\S+)',line)) for line in
          re.findall(r'\[mpi\]\[device\] (.*)',log)]
    if len(rows)!=ranks or {int(row['rank']) for row in rows}!=set(range(ranks)):
        raise RuntimeError('Missing/duplicate GPU rank identities')
    for host in {row['host'] for row in rows}:
        local=[row for row in rows if row['host']==host]
        if len({row['uuid'] for row in local})!=len(local):raise RuntimeError('GPU UUID is shared by local ranks')
        if any(not re.fullmatch('[0-9a-f]{32}',row['uuid']) or
               row.get('halo_transport')!='host_staged' for row in local):
            raise RuntimeError('Missing physical GPU or unexpected transport')
    if 'kernel=GPU' not in log or '[mpi][kernel] name=Stress' not in log:
        raise RuntimeError('The complete compact GPU MPI timestep backend was not exercised')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,required=True)
    p.add_argument('--gpu',type=Path,required=True)
    p.add_argument('--cpu3d',type=Path)
    p.add_argument('--gpu3d',type=Path)
    p.add_argument('--launcher',default='mpiexec --bind-to none')
    p.add_argument('--ranks',type=int,nargs='+',default=[1,2])
    p.add_argument('--cases',nargs='+',choices=['fixed','remesh','hydraulic','averaged'],default=['fixed','remesh','hydraulic','averaged'])
    p.add_argument('--steps',type=int,default=76)
    p.add_argument('--output-interval',type=int,default=19)
    p.add_argument('--rtol',type=float,default=1e-8)
    p.add_argument('--restart',action='store_true',help='Resume remeshed GPU checkpoints and compare equivalent CPU restarts')
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args()
    if bool(a.cpu3d)!=bool(a.gpu3d):p.error('Provide both 3D executables')
    if min(a.ranks)<1 or len(set(a.ranks))!=len(a.ranks) or 1 not in a.ranks:p.error('Use distinct positive ranks including one GPU')
    if a.steps<1 or a.output_interval<1 or a.steps%a.output_interval:p.error('Steps must be a positive multiple of output interval')
    if a.restart and a.steps//a.output_interval<4:p.error('Restart checks require at least four output intervals')
    if not 0<a.rtol<=1e-8:p.error('rtol must be positive and no larger than 1e-8')
    a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='gpu-mpi-',dir=str(a.workdir))).resolve()
    print('Artifacts:',work,flush=True)
    executables={key:str(getattr(a,key).resolve()) for key in ('cpu','gpu','cpu3d','gpu3d') if getattr(a,key)}
    report=dict(completed=False,all_passed=False,scope='Full GPU MPI timesteps',checks=[],executables=executables,
                executable_sha256={key:hashlib.sha256(Path(path).read_bytes()).hexdigest() for key,path in executables.items()},
                environment={name:os.environ[name] for name in ('JOB_ID','NSLOTS','PE_HOSTFILE','LD_LIBRARY_PATH','PATH','CUDA_VISIBLE_DEVICES') if name in os.environ})
    target=work/'results.json'
    def save():target.write_text(json.dumps(report,indent=2)+'\n')
    save()
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1')
    env['PYTHONPATH']=str(ROOT)+os.pathsep+env.get('PYTHONPATH','')
    def run(name,cfg,executable,ranks=None):
        folder=work/name;folder.mkdir()
        with (folder/'model.cfg').open('w') as stream:cfg.write(stream)
        command=[executable,'model.cfg'];settings=dict(env)
        if ranks:command=shlex.split(a.launcher)+['-n',str(ranks)]+command
        if ranks==1 and not report['checks']:
            settings.update(NV_ACC_NOTIFY='1',NVCOMPILER_ACC_NOTIFY='1')
        (folder/'command.json').write_text(json.dumps(command,indent=2)+'\n')
        log,wall=execute(command,folder,settings)
        if 'Ending simulation.' not in log:raise RuntimeError('Run did not finish: '+name)
        return folder,log,metrics(log,wall)
    try:
        for dim in ((2,3) if a.cpu3d else (2,)):
            suffix='' if dim==2 else '3d'
            for case in a.cases:
                moving=case in ('remesh','averaged')
                cfg=configure(dim,moving,1000,steps=a.steps,output_interval=a.output_interval)
                if case=='hydraulic':
                    cfg['control']['has_hydraulic_diffusion']='yes'
                    cfg['mat'].update(porosity='[0.01]',hydraulic_perm='[1e-20]')
                if case=='averaged':
                    cfg['sim']['is_outputting_averaged_fields']='yes'
                    cfg['mesh']['quality_check_step_interval']=str(a.output_interval)
                label='{}-{}d'.format(case,dim)
                reference,_,baseline=run(label+'-cpu',cfg,executables['cpu'+suffix])
                if moving and not baseline['remeshes']:raise RuntimeError('Remesh case did not remesh: '+label)
                restart_reference=None
                restart_frame=a.steps//a.output_interval//2
                continuation=range(restart_frame+1,a.steps//a.output_interval+1)
                # Averaged output restarts reset the parent's averaging window;
                # use non-averaged MMG checkpoints for trajectory restart checks.
                if a.restart and case=='remesh':
                    cfg['sim'].update(is_restarting='yes',restarting_from_modelname=str(reference/'result'),
                                      restarting_from_frame=str(restart_frame))
                    restart_reference,_,_=run(label+'-cpu-restart',cfg,executables['cpu'+suffix])
                    for frame in continuation:
                        compare_fields(reference,restart_reference,frame,rtol=a.rtol)
                    cfg['sim']['is_restarting']='no'
                single=None
                for ranks in sorted(a.ranks):
                    current,log,timing=run(label+'-gpu'+str(ranks),cfg,executables['gpu'+suffix],ranks)
                    devices=device_records(log,ranks);checks=[]
                    for kind in ('save','chkpt'):
                        for frame in range(a.steps//a.output_interval+1):
                            checks.append(compare_fields(reference,current,frame,kind,a.rtol))
                            if single:compare_fields(single,current,frame,kind,a.rtol)
                    if timing['remeshes']!=baseline['remeshes']:raise RuntimeError('Remesh count differs: '+label)
                    if ranks==1:single=current
                    report['checks'].append(dict(case=label,ranks=ranks,devices=devices,timing=timing,comparisons=checks,
                                                 kernel_trace_lines=log.count('launch CUDA kernel')))
                    save();print(label,'ranks',ranks,'full physical outputs/checkpoints PASS',flush=True)
                    if restart_reference:
                        cfg['sim'].update(is_restarting='yes',restarting_from_modelname=str(current/'result'),
                                          restarting_from_frame=str(restart_frame))
                        resumed,resume_log,resume_timing=run(label+'-gpu'+str(ranks)+'-restart',cfg,
                                                             executables['gpu'+suffix],ranks)
                        resume_checks=[compare_fields(restart_reference,resumed,frame,kind,a.rtol)
                                       for kind in ('save','chkpt') for frame in continuation]
                        for frame in continuation:compare_fields(current,resumed,frame,rtol=a.rtol)
                        report['checks'].append(dict(case=label+'-restart',ranks=ranks,
                                                     devices=device_records(resume_log,ranks),timing=resume_timing,
                                                     comparisons=resume_checks,checkpoint_source=str(current)))
                        cfg['sim']['is_restarting']='no'
                        save();print(label,'ranks',ranks,'own checkpoint restart PASS',flush=True)
    except Exception as error:
        report['failure']=str(error);save();raise
    report.update(completed=True,all_passed=True);save()
    print('Passed. Results:',target,flush=True)
    return 0


if __name__=='__main__':raise SystemExit(main())
