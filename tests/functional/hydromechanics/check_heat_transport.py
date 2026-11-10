#!/usr/bin/env python3
"""Check the integrated Ogata-Banks problem, mixed capacity and restart.

The steady pressure field fixes q. A deliberately enlarged hydraulic storage
keeps the explicit hydraulic solve cheap without changing q or the thermal
reference. These scaled coefficients isolate transport, not an application.
"""
import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.special import erfc, erfcx

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from Dynearthsol import Dynearthsol, scan_frames, YEAR2SEC


def exact(x,t,v,d):
    if t == 0:
        return np.where(x == 0,1.,0.)
    root = 2*math.sqrt(d*t)
    a,b = (x-v*t)/root,(x+v*t)/root
    return .5*(erfc(a)+np.exp(-a*a)*erfcx(b))


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def study(exe,out,spacings,steps,porosity,restart):
    out.mkdir(parents=True,exist_ok=False)
    (out/'runner.py').write_bytes(Path(__file__).read_bytes())
    result=dict(executable=str(exe),executable_sha256=digest(exe),runner_sha256=digest(Path(__file__)),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git','diff','HEAD','--binary'],cwd=ROOT))
    rho_s,cs,rho_f,cf,q,conductivity=2500.,2000.,1000.,4000.,1.,23500.
    capacity=(1-porosity)*rho_s*cs+porosity*rho_f*cf
    velocity=rho_f*cf*q/capacity
    diffusivity=conductivity/capacity
    duration=.4/velocity
    result['reference']=dict(bulk_capacity=capacity,velocity=velocity,diffusivity=diffusivity,
                            duration_s=duration,porosity=porosity,scaled_hydraulic_storage=True)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(folder,cfg,name):
        with (folder/f'{name}.cfg').open('w') as f:cfg.write(f)
        with (folder/f'{name}.log').open('w') as f:
            subprocess.run([str(exe),f'{name}.cfg'],cwd=folder,env=env,stdout=f,
                           stderr=subprocess.STDOUT,check=True,timeout=180)
    for h in spacings:
        for count in steps:
            folder=out/f'h{h:g}-n{count}';folder.mkdir()
            cfg=configparser.ConfigParser();cfg.optionxform=str
            cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg['sim'].update(modelname='column',max_steps=str(count),output_step_interval=str(count//4),
                output_time_interval_in_yr='1',max_time_in_yr=str(10*duration/YEAR2SEC))
            cfg['mesh'].update(xlength='1',ylength=str(2*h),zlength=str(2*h),resolution=str(h),meshing_elem_shape='1')
            cfg['control'].update(has_PT='false',is_quasi_static='true',inertial_scaling='1e6',
                use_global_velocity_scaling='false',has_moving_mesh='false',gravity='0',damping_option='0',
                fixed_dt=repr(duration/count),has_thermal_diffusion='true',has_fluid_thermal_advection='true',
                has_poroelastic_pressure_feedback='false')
            cfg['ic'].update(has_body_force_adjustment='false',temperature_option='4',temp_background='273.15',
                             temp_anomaly_amplitude='0',pore_pressure_option='1')
            for side in ('x0','x1','y0','y1','z0','z1'):
                cfg['bc'][f'vbc_{side}']='1';cfg['bc'][f'hbc_{side}']='0';cfg['bc'][f'tbc_{side}']='0'
            cfg['bc'].update(hbc_x0='2',hbc_x1='2',hbc_pressure='[2,1,0,0,0,0]',
                tbc_x0='1',tbc_val_x0='274.15',stress_bc_z1='0',stress_val_z1='0')
            cfg['mat'].update(rho0=f'[{rho_s}]',alpha='[0]',heat_capacity=f'[{cs}]',therm_cond=f'[{conductivity}]',
                fluid_rho0=f'[{rho_f}]',fluid_heat_capacity=f'[{cf}]',fluid_visc='[1]',hydraulic_perm='[1]',
                porosity=f'[{porosity}]',biot_coeff='[1]',fluid_bulk_modulus='[1e-6]',
                bulk_modulus='[1e9]',shear_modulus='[1e9]',retention_model='[0]')
            run(folder,cfg,'fresh')
            model=Dynearthsol(str(folder/'column'));rows=scan_frames(str(folder/'column'))
            x=model.read_field(0,'coordinate')[:,0]
            errors=[];prange=[]
            for row in rows:
                t=model.read_field(row['frame'],'temperature')-273.15
                p=model.read_field(row['frame'],'pore pressure')
                assert np.isfinite(t).all()
                assert t.min()>=-1e-10 and t.max()<=1+1e-10
                errors.append(float(np.max(np.abs(t-exact(x,row['time'],velocity,diffusivity)))))
                prange.append(float(np.max(np.abs(p-(2-x)))))
            item=dict(spacing_m=h,steps=count,dimension=model.ndims,nodes=len(x),
                elements=len(model.read_field(0,'connectivity')),max_error=errors,
                max_steady_pressure_error=max(prange),input_sha256=digest(folder/'fresh.cfg'))
            assert max(prange)<1e-10,prange
            if restart:
                quarter=rows[1]
                cfg['sim'].update(modelname='continued',is_restarting='true',
                    restarting_from_modelname=str(folder/'column'),restarting_from_frame=str(quarter['frame']))
                run(folder,cfg,'restart')
                resumed=Dynearthsol(str(folder/'continued'))
                fresh={r['steps']:r for r in rows}
                for row in scan_frames(str(folder/'continued')):
                    f=fresh[row['steps']];assert row['time']==f['time']
                    for key in ('temperature','pore pressure','velocity','stress','strain','coordinate'):
                        np.testing.assert_array_equal(model.read_field(f['frame'],key),resumed.read_field(row['frame'],key))
                item['physical_restart_exact']=True
            item['outputs_sha256']={p.name:digest(p) for p in folder.glob('column.*')
                if '.save.' in p.name or '.chkpt.' in p.name or p.suffix=='.info'}
            result['cases'].append(item)
            (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False))
            print(json.dumps(item|{'outputs_sha256':'recorded'}),flush=True)
    assert digest(exe)==result['executable_sha256']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--spacings',type=float,nargs='+',default=[.05,.025,.0125])
    p.add_argument('--steps',type=int,nargs='+',default=[2048])
    p.add_argument('--porosity',type=float,default=.3);p.add_argument('--restart',action='store_true')
    a=p.parse_args()
    if (not 0<a.porosity<=1 or any(not math.isfinite(h) or h<=0 for h in a.spacings) or
        any(n<4 or n%4 for n in a.steps)):p.error('Use finite positive spacings, positive step counts divisible by four and 0<porosity<=1')
    study(a.executable.resolve(),a.output.resolve(),a.spacings,a.steps,a.porosity,a.restart)
