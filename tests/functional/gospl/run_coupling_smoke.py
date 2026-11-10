#!/usr/bin/env python3
"""Bounded tutorial coupling check: real CPU/OpenMP and single-GPU DES runs.

Checks the DES/GoSPL boundary, not GoSPL's surface algorithms or scaling.
The GPU control resumes the CPU's initial checkpoint and reuses its surface mesh.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import h5py
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tests/functional/mpi'))
from run_gpu_mpi import compare_fields
YEAR_SECONDS=365.2422*86400

def audit(folder, gospl_ranks=1, verify_seed=False):
    audit_dir=folder/'audit'
    if (audit_dir/'rank-0000').exists():audit_dir=audit_dir/'rank-0000'
    rows=[json.loads(line) for line in (audit_dir/'calls.jsonl').read_text().splitlines()]
    steps=[row for row in rows if row['operation']=='step']
    if len(steps)!=4 or any(not row['returned'] for row in steps):raise RuntimeError('Expected four successful coupling events')
    if sum(row['operation']=='seed' for row in rows)!=1:raise RuntimeError('DES must seed GoSPL elevation once')
    if [row['mpi_size'] for row in rows if row['operation']=='create']!=[gospl_ranks]:raise RuntimeError('Unexpected GoSPL MPI world size')
    if any(row.get('status',0)!=0 for row in rows):raise RuntimeError('Bridge operation failed')
    if not rows or rows[-1]['operation']!='destroy' or rows[-1]['remaining_models']!=0:raise RuntimeError('GoSPL cleanup did not finish')
    seed_errors=[]
    if verify_seed:
        from scipy.spatial import cKDTree
        directories=list((folder/'audit').glob('rank-*')) or [audit_dir]
        if len(directories)!=gospl_ranks:raise RuntimeError('Missing GoSPL rank seed audit')
        for directory in directories:
            seed=np.load(directory/'seed.npz')
            for kind in ('owned','local'):
                distance,ids=cKDTree(seed['coords']).query(seed[kind+'_coords'],k=3)
                weight=1/np.maximum(distance,1e-20)
                expected=(weight*seed['elevation'][ids]).sum(axis=1)/weight.sum(axis=1)
                on=distance[:,0]<1e-20;expected[on]=seed['elevation'][ids[on,0]]
                error=float(np.max(abs(expected-seed[kind+'_height'])))
                if error>1e-9:raise RuntimeError('Incorrect GoSPL '+kind+' seed placement')
                seed_errors.append(error)
    previous=None
    previous_remesh=None
    events=[]
    for row in steps:
        number=row['event'];data=np.load(audit_dir/('event-%04d.npz'%number))
        coords,velocity,dh=data['coords'],data['velocity'],data['dh']
        if velocity.shape!=coords.shape or dh.shape!=(len(coords),):raise RuntimeError('Coupling array shapes differ')
        if not all(np.isfinite(x).all() for x in (coords,velocity,dh)):raise RuntimeError('Nonfinite boundary data')
        with h5py.File(folder/('result.save.%06d.vtkhdf'%number)) as file:
            grid=file['VTKHDF/grid'];points=grid['Points'][()];vel=grid['PointData/velocity'][()]
            clock=float(grid['FieldData/time_yr'][()].item())
            # Match unchanged horizontal positions, then verify the returned
            # erosion increment was added to each DES surface height once.
            from scipy.spatial import cKDTree
            flags=grid['PointData/bcflag'][()]
            top=np.flatnonzero((flags&32)!=0)
            distance,local=cKDTree(points[top,:2]).query(coords[:,:2])
            indices=top[local]
            if len(np.unique(indices))!=len(coords) or np.max(distance)>1e-7:raise RuntimeError('Surface point identities changed')
            error=float(np.max(np.abs(points[indices,2]-(coords[:,2]+dh))))
            if error>1e-9:raise RuntimeError('DES did not apply returned erosion exactly once')
            if not np.isclose(row['time_after'],clock,rtol=1e-12,atol=1e-12):raise RuntimeError('DES/GoSPL clocks use different year units')
            if not np.isclose(row['time_after']-row['time_before'],row['dt'],rtol=1e-12,atol=1e-12):raise RuntimeError('GoSPL advanced a different interval')
            with h5py.File(folder/('result.chkpt.%06d.vtkhdf'%number)) as checkpoint:
                remesh=int(checkpoint['VTKHDF/grid/nremesh'][()].item())
            instantaneous=previous is None or remesh!=previous_remesh
            expected=vel[indices]*YEAR_SECONDS if instantaneous else (coords-previous)/row['dt']
            scale=max(float(np.max(np.abs(expected))),1e-30)
            if float(np.max(np.abs(velocity-expected)))/scale>1e-8:raise RuntimeError('Surface velocity transfer differs')
            previous=points[indices].copy()
            previous_remesh=remesh
            events.append(dict(event=number,points=len(coords),time_yr=clock,height_application_error=error,
                               dh_max=float(np.max(np.abs(dh)))))
    return dict(passed=True,events=events,seed_placement_errors=seed_errors,call_seconds={kind:sum(r['call_s'] for r in rows if r['operation']==kind)
                for kind in ('create','seed','velocity','step','destroy')},
                audit_write_seconds=sum(r.get('audit_write_s',0) for r in rows))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cpu',type=Path,required=True)
    p.add_argument('--gpu',type=Path,required=True)
    p.add_argument('--extension',type=Path,required=True)
    p.add_argument('--workdir',type=Path,required=True)
    a=p.parse_args();a.workdir.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='coupling-',dir=str(a.workdir))).resolve()
    print('Artifacts:',work,flush=True)
    source=ROOT/'gospl_driver/examples/gaussian-weakzone-3d-with-gospl.cfg'
    cfg=configparser.ConfigParser(inline_comment_prefixes=('#',';'));cfg.optionxform=str;cfg.read(source)
    cfg['sim'].update(modelname='result',max_steps='8',max_time_in_yr='1e20',output_step_interval='2',
                      output_time_interval_in_yr='1e20',info_display_step_interval='2',checkpoint_frame_interval='1',is_restarting='no')
    # A uniform coarse box keeps this check independent of the tutorial's
    # refined-zone PLC construction and concentrates on the coupling boundary.
    cfg['mesh'].update(meshing_option='1',meshing_elem_shape='1',resolution='10000',quality_check_step_interval='1000')
    cfg['control'].update(surface_process_gospl_config_file='gospl.yml',gospl_coupling_frequency='2',
                          gospl_mesh_resolution='10000',gospl_mesh_perturbation='0')
    surface=yaml.safe_load((source.parent/'gospl_config_gaussian_weakzone_3D.yml').read_text())
    surface['time'].update(dt=1.,tout=1.,end=1e20);surface['domain']['seadepo']=False
    surface['output']['dir']='gospl-output'
    report=dict(completed=False,scope='DES coupling boundary only',checks=[],executables={},steps=8,
                coupling_events=4,des_resolution=10000,gospl_resolution=10000,
                note='Uniform coarse initialization; marine deposition is inactive below the surface and disabled for this bounded coupling check.')
    target=work/'results.json'
    def save():target.write_text(json.dumps(report,indent=2)+'\n')
    save()
    try:
        for kind,exe,threads in (('cpu',a.cpu,2),('gpu',a.gpu,1)):
            exe=exe.resolve();report['executables'][kind]=dict(path=str(exe),sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
            folder=work/kind;folder.mkdir()
            if kind=='gpu':
                cfg['sim'].update(is_restarting='yes',restarting_from_modelname=str(work/'cpu/result'),restarting_from_frame='0')
                shutil.copy2(work/'cpu/gospl_mesh.npz',folder/'gospl_mesh.npz')
            with (folder/'model.cfg').open('w') as stream:cfg.write(stream)
            (folder/'gospl.yml').write_text(yaml.safe_dump(surface))
            env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',
                     GOSPL_REAL_INTERFACE=str(a.extension.resolve()/'cpp_interface/gospl_python_interface.py'),
                     GOSPL_AUDIT_DIR=str(folder/'audit'),CUDA_VISIBLE_DEVICES='0')
            env['PYTHONPATH']=str(Path(__file__).with_name('audit_interface'))+os.pathsep+env.get('PYTHONPATH','')
            start=time.perf_counter()
            with (folder/'run.log').open('w') as log:
                result=subprocess.run([str(exe),'model.cfg'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
            elapsed=time.perf_counter()-start
            text=(folder/'run.log').read_text()
            if result.returncode or 'Ending simulation.' not in text or 'Traceback' in text:raise RuntimeError('Coupled run failed: '+kind)
            check=audit(folder);check.update(kind=kind,des_openmp_threads=threads,wall_seconds=elapsed)
            report['checks'].append(check);save();print(kind,'coupling boundary PASS',flush=True)
        report['comparisons']=[compare_fields(work/'cpu',work/'gpu',frame,kind)
                               for kind,frames in (('save',range(5)),('chkpt',range(1,5))) for frame in frames]
        report.update(completed=True,all_passed=True);save()
    except Exception as error:
        report['failure']=str(error);save();raise
    print('Results:',target,flush=True)

if __name__=='__main__':main()
