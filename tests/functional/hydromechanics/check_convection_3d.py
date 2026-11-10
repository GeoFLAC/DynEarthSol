#!/usr/bin/env python3
"""Three-dimensional porous-cell convection against a separable eigenmode.

The cos(pi*x) cos(pi*y) sin(pi*(z+1)) mode varies in all three directions.
Its horizontal wavenumber squared is 2*pi^2, vertical pi^2, and critical
Rayleigh number 4.5*pi^2. A small Gaussian seeds this and other modes; the
volume projection measures its growth independently of the other linear
modes. This tests linear 3D cells, not nonlinear steady planform selection.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import check_convection as cc
from convection_reference import critical_rayleigh,growth_rate,storage_growth_rate


def series(folder):
    model=cc.Dynearthsol(str(folder/'box'))
    rows=cc.scan_frames(str(folder/'box'))
    x=model.read_field(0,'coordinate');conn=model.read_field(0,'connectivity')
    assert x.shape[1]==3
    volume=np.abs(np.linalg.det(x[conn][:,1:]-x[conn][:,:1]))/math.factorial(3)
    weights=np.zeros(len(x))
    np.add.at(weights,conn,np.broadcast_to(volume[:,None]/conn.shape[1],conn.shape))
    assert abs(weights.sum()-1)<1e-10
    basis=8*weights*np.cos(math.pi*x[:,0])*np.cos(math.pi*x[:,1])*np.sin(math.pi*(x[:,2]+1))
    result=[]
    for row in rows:
        theta=model.read_field(row['frame'],'temperature')-(cc.T_TOP-cc.DT*x[:,2])
        assert np.isfinite(theta).all()
        assert np.isfinite(model.read_field(row['frame'],'pore pressure')).all()
        result.append([row['time'],float(basis@theta),float(abs(theta).max())])
    return result,len(x),len(conn)


def study(exe,out,spacing,ratio,steps,dt,seed,threads,fluid_bulk=1000.):
    out.mkdir(parents=True,exist_ok=False)
    critical=critical_rayleigh(1,1,transverse=1,width_y=1)
    ra=ratio*critical
    cfg=cc.config(spacing,ra/critical_rayleigh(1,1),seed,steps,dt,40)
    cfg['mesh']['meshing_elem_shape']='1'
    cfg['mat']['fluid_bulk_modulus']=repr([fluid_bulk])
    cfg['ic']['temp_anomaly_ycenter']='0.25'
    cfg['sim'].update(has_marker_output='false',has_initial_checkpoint='false',checkpoint_frame_interval='100000')
    with (out/'box.cfg').open('w') as stream:cfg.write(stream)
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (out/'run.log').open('w') as stream:
        subprocess.run([str(exe),'box.cfg'],cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,
                       check=True,timeout=14400)
    measured,nodes,elements=series(out)
    result=dict(dimension=3,spacing_m=spacing,nodes=nodes,elements=elements,
        rayleigh=ra,critical_rayleigh=critical,mode=[1,1,1],seed_K=seed,steps=steps,dt_s=dt,threads=threads,fluid_bulk_Pa=fluid_bulk,
        series=measured,executable_sha256=cc.digest(exe),runner_sha256=cc.digest(Path(__file__)),
        shared_runner_sha256=cc.digest(Path(cc.__file__)),reference_sha256=cc.digest(Path(cc.ROOT/'benchmarks/poromechanics/convection_reference.py')),
        input_sha256=cc.digest(out/'box.cfg'),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=cc.ROOT,text=True).strip(),
        scope='Rigid saturated Darcy-Boussinesq matrix; linear 3D mode with exact finite-storage reference; no nonlinear steady-cell claim')
    if seed:
        values=np.array(measured);late=values[len(values)//2:]
        fitted=float(np.polyfit(late[:,0],np.log(abs(late[:,1])),1)[0])
        value=lambda name:json.loads(cfg['mat'][name])[0]
        storage=value('porosity')/fluid_bulk+(value('biot_coeff')-value('porosity'))/value('bulk_modulus_s')+value('biot_coeff')**2/value('bulk_modulus')
        hydraulic_ratio=value('hydraulic_perm')/value('fluid_visc')/storage/(cc.LAMBDA/cc.CAPACITY)
        exact=storage_growth_rate(ra,hydraulic_ratio,cc.LAMBDA,cc.CAPACITY,1,1,transverse=1,width_y=1)
        result.update(fitted_rate_s_inv=fitted,reference_rate_s_inv=exact,
                      hydraulic_to_thermal_diffusivity=hydraulic_ratio,reference_model='finite-storage',
                      incompressible_rate_s_inv=growth_rate(ra,cc.LAMBDA,cc.CAPACITY,1,1,transverse=1,width_y=1),
                      relative_growth_error=abs(fitted/exact-1),correct_sign=bool(fitted*exact>0),
                      linearity_max_theta_over_dT=float(values[:,2].max()/cc.DT))
    else:
        result['rest_max_theta_over_dT']=max(row[2] for row in measured)/cc.DT
    result['outputs_sha256']={p.name:cc.digest(p) for p in sorted(out.glob('box.*'))
                             if '.save.' in p.name or p.suffix=='.info'}
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('series','outputs_sha256')},indent=2),flush=True)
    if seed:
        assert result['correct_sign']
        assert result['linearity_max_theta_over_dT']<.01
    else:assert result['rest_max_theta_over_dT']<1e-6
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--spacing',type=float,default=.0625)
    parser.add_argument('--ratio',type=float,default=1.4)
    parser.add_argument('--steps',type=int,default=400)
    parser.add_argument('--dt',type=float,default=500.)
    parser.add_argument('--seed',type=float,default=1e-4)
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--fluid-bulk',type=float,default=1000.)
    args=parser.parse_args()
    if (not math.isfinite(args.spacing) or args.spacing<=0 or args.ratio<=0 or args.ratio==1 or
            args.steps<40 or args.steps%40 or args.dt<=0 or args.seed<0 or args.threads<1 or args.fluid_bulk<=0):
        parser.error('Use positive spacing/ratio/dt/threads, ratio !=1, nonnegative seed, and steps divisible by40')
    if args.threads>int(os.environ.get('NSLOTS',args.threads)):
        parser.error('Threads exceed the scheduler allocation')
    study(args.executable.resolve(),args.output.resolve(),args.spacing,args.ratio,args.steps,args.dt,args.seed,args.threads,args.fluid_bulk)
