#!/usr/bin/env python3
"""Qualify signed fluid/heat schedules by independent closed-box balances.

The analytical references are integrated fluid volume and sensible heat,
not a fitted spatial point-source profile. Fixed-mesh restart must preserve
every saved/checkpoint physical payload. A lifted column remeshed at every
step checks that point sources and a well screen move with the material, as
monitor points do, and continue exactly on restart.
"""

import argparse
import ast
import csv
import json
import math
import os
from pathlib import Path
import shutil
import struct
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'tests/functional'))
from Dynearthsol import Dynearthsol,scan_frames,YEAR2SEC
from pt_initial_checkpoint import fields


from des_case import digest, read_config, reset_boundaries, run_solver


def config():
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='fresh',max_steps='4',output_step_interval='1',
        max_time_in_yr=str(100/YEAR2SEC),output_time_interval_in_yr='1')
    cfg['mesh'].update(resolution='.25',largest_size='.25',smallest_size='.25')
    cfg['control'].update(gravity='0',has_PT='false',is_quasi_static='true',inertial_scaling='1e11',
        has_moving_mesh='false',has_thermal_diffusion='false',fixed_dt='1',
        has_poroelastic_pressure_feedback='false',has_fluid_thermal_advection='false')
    cfg['ic'].update(temperature_option='4',temp_background='300',temp_anomaly_amplitude='0',
                     has_body_force_adjustment='false',pore_pressure_option='2',pore_pressure_value='0')
    reset_boundaries(cfg['bc'], vbc='1', hbc='0', tbc='0')
    cfg['bc'].update(stress_bc_z1='0',stress_val_z1='0')
    cfg['mat'].update(rho0='[2500]',alpha='[0]',heat_capacity='[1000]',therm_cond='[0]',
                     hydraulic_perm='[0]',fluid_heat_capacity='[4000]',fluid_alpha='[0]')
    for group,rates in (('injection',[8e-8,-3e-8]),('heat_injection',[4000,-1500])):
        cfg[group]=dict(enabled='true',num_points='2',points_x='[0.37,0.61]',points_y='[0.43,0.57]',
            points_z='[-0.53,-0.41]',points_unit='m',rate_model='constant_rate',rate=json.dumps(rates),
            start_time_in_yr=json.dumps([.25/YEAR2SEC,.8/YEAR2SEC]),
            end_time_in_yr=json.dumps([2.75/YEAR2SEC,2.1/YEAR2SEC]))
    return cfg


def amount(t,rate):
    return float(np.dot(rate,np.maximum(0,np.minimum(t,[2.75,2.1])-[.25,.8])))


def check_balance(folder,cfg):
    model=Dynearthsol(str(folder/'fresh'));dim=model.ndims
    x=model.read_field(0,'coordinate');conn=model.read_field(0,'connectivity')
    cells=x[conn];vol=np.abs(np.linalg.det(cells[:,1:]-cells[:,0,None]))/math.factorial(dim)
    nodal=np.zeros(len(x));np.add.at(nodal,conn.ravel(),np.repeat(vol/(dim+1),dim+1))
    mat=lambda k:json.loads(cfg['mat'][k])[0]
    # DES records in-plane mean stress in 2D; its fixed-stress compliance is
    # 1/(K+G/3), rather than the three-dimensional 1/K.
    bulk=mat('bulk_modulus')+(mat('shear_modulus')/3 if dim==2 else 0)
    sb=mat('porosity')/mat('fluid_bulk_modulus')+(mat('biot_coeff')-mat('porosity'))/mat('bulk_modulus_s')+mat('biot_coeff')**2/bulk
    capacity=mat('rho0')*mat('heat_capacity')
    if cfg.getboolean('control','has_fluid_thermal_advection'):
        capacity=(1-mat('porosity'))*capacity+mat('porosity')*mat('fluid_rho0')*mat('fluid_heat_capacity')
    p0=model.read_field(0,'pore pressure');t0=model.read_field(0,'temperature')
    report=[]
    for row in scan_frames(str(folder/'fresh')):
        p=model.read_field(row['frame'],'pore pressure');temp=model.read_field(row['frame'],'temperature')
        assert np.isfinite(p).all() and np.isfinite(temp).all()
        fluid=float(np.dot(nodal*sb,p-p0));energy=float(np.dot(nodal*capacity,temp-t0))
        expected_fluid=amount(row['time'],[8e-8,-3e-8]) if cfg.getboolean('injection','enabled') else 0.
        if cfg.getint('bc','num_wells',fallback=0):expected_fluid-=5e-9*row['time']
        expected_energy=amount(row['time'],[4000,-1500]) if cfg.getboolean('heat_injection','enabled') else 0.
        np.testing.assert_allclose(fluid,expected_fluid,rtol=3e-12,atol=1e-18)
        np.testing.assert_allclose(energy,expected_energy,rtol=3e-10,atol=1e-6)
        report.append(dict(time_s=row['time'],fluid_m3=fluid,expected_fluid_m3=expected_fluid,
                           heat_J=energy,expected_heat_J=expected_energy))
    return report


def run(exe,folder,cfg,env,expect=0):
    folder.mkdir()
    with (folder/'fresh.cfg').open('w') as f:cfg.write(f)
    p = run_solver([str(exe), 'fresh.cfg'], log=folder / 'fresh.log', cwd=folder, env=env, timeout=180)
    if p.returncode!=expect:
        raise AssertionError((folder,p.returncode,expect,(folder/'fresh.log').read_text()[-2500:]))


def restart(exe,folder,cfg,env):
    cfg['sim'].update(modelname='resumed',is_restarting='true',
        restarting_from_modelname=str(folder/'fresh'),restarting_from_frame='1')
    with (folder/'restart.cfg').open('w') as f:cfg.write(f)
    run_solver([str(exe), 'restart.cfg'], log=folder / 'restart.log', cwd=folder, env=env, check=True, timeout=180)
    count=0;clock_rounding=[]
    for suffix in ('save','chkpt'):
        for p in sorted(folder.glob(f'resumed.{suffix}.*')):
            # The shared payload reader resolves the container from a bare
            # frame name; passing an HDF5 filename selects its native parser.
            frame=p.with_suffix('') if p.suffix=='.vtkhdf' else p
            a,b=fields(folder/frame.name.replace('resumed.','fresh.',1)),fields(frame)
            assert a.keys()==b.keys()
            for key in a:
                base=key.rsplit('/',1)[-1].lstrip('@').removesuffix('.layout')
                if base in ('walltime_sec','provenance') or '/provenance/' in key:continue
                if base=='reference_frame_time' and a[key]!=b[key]:
                    # The existing non-PT output scheduler rebases its clock at
                    # restart. Record its ulp-level regrouping difference; do
                    # not confuse that unresolved lifecycle issue with exact
                    # equality of physical state and actual saved times.
                    before,after=struct.unpack('d',a[key])[0],struct.unpack('d',b[key])[0]
                    np.testing.assert_allclose(before,after,rtol=4*np.finfo(float).eps,atol=0)
                    clock_rounding.append(dict(file=p.name,fresh=before,resumed=after))
                    continue
                assert a[key]==b[key],(p,key)
                count+=1
    assert count>0
    return dict(exact_physical_payloads=count,legacy_output_clock_rounding=clock_rounding)


def checkpoint_arrays(path,names):
    """Named arrays of a checkpoint, native or HDF5."""
    data=fields(path);found={}
    for name in names:
        key=next(k for k in data if not k.startswith('@') and k.split('/')[-1]==name)
        layout=data.get('@'+key+'.layout')
        kind=ast.literal_eval(layout.decode())[1] if layout else ('<i4' if name.endswith(('elements','wells')) else '<f8')
        found[name]=np.frombuffer(data[key],dtype=np.dtype(kind))
    return found


def material(exe,out,dim,env):
    """Sources are fixed in the material, as monitor points are. A column lifted
    rigidly and remeshed at every step carries a point source and a well screen;
    the positions their checkpoints give must follow the monitor points that
    start with them, the screen must keep its whole rate, and a restart must
    continue them exactly."""
    lift,dt,steps,point_rate,well_rate=1.,.1,4,1e-10,3e-10
    point,bottom,top=[.4,.4,-5],[.6,.6,-6.5],[.6,.6,-3.5]
    cfg=read_config(ROOT/'tests/functional/pt_remesh_restart.cfg')
    cfg['sim']['has_output_during_remeshing']='no'
    cfg['control'].update(fixed_dt=str(dt),has_thermal_diffusion='true',gravity='0',PT_absolute_tolerance='1e-8')
    cfg['bc'].update(vbc_val_z0=str(lift),num_wells='1',well_x=f'[{bottom[0]}]',well_y=f'[{bottom[1]}]',
                     well_z_bottom=f'[{bottom[2]}]',well_z_top=f'[{top[2]}]',well_rate=f'[{well_rate}]')
    for group,rate in (('injection',f'[{point_rate}]'),('heat_injection','[100]')):
        cfg[group]=dict(enabled='true',num_points='1',points_x=f'[{point[0]}]',points_y=f'[{point[1]}]',
                        points_z=f'[{point[2]}]',rate=rate,rate_model='constant_rate')
    starts=[point,bottom,top]
    cfg['monitor']=dict(enabled='true',step_interval='1',num_points='3',point_type='particle',
                        points_x=json.dumps([q[0] for q in starts]),points_y=json.dumps([q[1] for q in starts]),
                        points_z=json.dumps([q[2] for q in starts]))
    if dim==2:
        for section in ('injection','heat_injection','monitor'):del cfg[section]['points_y']
        del cfg['bc']['well_y']
    with (out/'remesh.cfg').open('w') as f:cfg.write(f)
    run_solver([sys.executable,str(ROOT/'tests/functional/pt_initial_checkpoint.py'),str(exe),str(out/'remesh.cfg'),
                str(out/'remesh'),'--steps',str(steps),'--restart-frame','1','--require-remesh'],
               log=out/'remesh.log',env=env,check=True,timeout=240)
    folder=out/'remesh';model=Dynearthsol(str(folder/'fresh'));nodes=dim+1
    axes=('x','y','z') if dim==3 else ('x','z')
    tracked={}
    for i in range(3):
        with (folder/f'monitor_point_{i}.csv').open() as stream:
            for row in csv.DictReader(stream):
                if int(row['row'])==0:
                    tracked[i,int(row['step'])]=np.array([float(row['query_'+a]) for a in axes])
    names=['fluid point elements','fluid point weights','heat point elements','heat point weights',
           'well piece elements','well piece starts','well piece ends']
    gap=0.
    for frame in model.frames:
        x=model.read_field(frame,'coordinate');conn=model.read_field(frame,'connectivity')
        a=checkpoint_arrays(folder/f'fresh.chkpt.{frame:06d}',names)
        where=lambda element,weights:weights[:nodes]@x[conn[element]]
        found=[where(a['fluid point elements'][0],a['fluid point weights']),
               where(a['heat point elements'][0],a['heat point weights']),
               where(a['well piece elements'][0],a['well piece starts'][:nodes]),
               where(a['well piece elements'][-1],a['well piece ends'][-nodes:])]
        for source,monitor in zip(found,(0,0,1,2)):
            gap=max(gap,float(np.max(np.abs(source-tracked[monitor,frame]))))
        if frame:
            b=model.read_field(frame,'fluid balance')
            np.testing.assert_allclose(b[2],(well_rate-point_rate)*dt,rtol=1e-9,atol=0)
    assert gap<=1e-9,('Sources left the material their monitor points follow',gap)
    start=np.array([point[0],point[2]] if dim==2 else point)
    moved=float(np.max(np.abs(found[0]-start)))
    assert moved>=.5*lift*dt*steps,('The source did not move with the lifted column',moved)
    return dict(actual_mesh_change=True,all_continuation_payloads_exact=True,
                max_source_monitor_gap_m=gap,point_displacement_m=moved)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--dimension',type=int,choices=(2,3),required=True)
    args=parser.parse_args();dim=args.dimension
    exe,out=args.executable.resolve(),args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(Path(__file__),out/'runner.py')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    result=dict(executable_sha256=digest(exe),checks={},runner_sha256=digest(Path(__file__)))
    modes=['combined','total-amount','half-step','subcycled','heat-only','fluid-only','mixed-capacity',
           'millimetres','centimetres','kilometres','point-and-well']
    if dim==2:modes.append('legacy-coordinate')
    for mode in modes:
        cfg=config()
        # The two-coordinate legacy alias is tested in 2D after probing a base run.
        if mode=='total-amount':
            for group,rates in (('injection',[8e-8,-3e-8]),('heat_injection',[4000,-1500])):
                cfg[group].update(rate_model='total_amount',total_amount=json.dumps((np.array(rates)*[2.5,1.3]).tolist()),
                                  duration_in_yr=json.dumps([2.5/YEAR2SEC,1.3/YEAR2SEC]))
                del cfg[group]['end_time_in_yr']
        if mode=='half-step':cfg['control']['fixed_dt']='.5';cfg['sim']['max_steps']='8'
        if mode=='subcycled':
            cfg['control']['has_thermal_diffusion']='true'
            cfg['mat'].update(hydraulic_perm='[1e-10]',therm_cond='[1e6]')
        if mode=='heat-only':cfg['injection']['enabled']='false';cfg['control']['has_hydraulic_diffusion']='false'
        if mode=='fluid-only':cfg['heat_injection']['enabled']='false'
        if mode=='mixed-capacity':cfg['control']['has_fluid_thermal_advection']='true'
        if mode=='point-and-well':
            cfg['bc'].update(num_wells='1',well_x='[0.5]',well_y='[0.5]',
                well_z_top='[-0.2]',well_z_bottom='[-0.8]',well_rate='[5e-9]')
        if mode in ('millimetres','centimetres','kilometres'):
            unit,scale={'millimetres':('mm',1000),'centimetres':('cm',100),'kilometres':('km',.001)}[mode]
            for group in ('injection','heat_injection'):
                cfg[group]['points_unit']=unit
                for name in ('points_x','points_y','points_z'):
                    cfg[group][name]=json.dumps((np.array(json.loads(cfg[group][name]))*scale).tolist())
        # The solver rejects ambiguous legacy-y and modern-z source coordinates in 2D.
        if dim==2:
            for group in ('injection','heat_injection'):
                del cfg[group]['points_y']
                if mode=='legacy-coordinate':cfg[group]['points_y']=cfg[group].pop('points_z')
        folder=out/mode;run(exe,folder,cfg,env)
        assert Dynearthsol(str(folder/'fresh')).ndims==dim
        records=check_balance(folder,cfg);payloads=restart(exe,folder,cfg,env)
        result['checks'][mode]=dict(balance=records,restart=payloads)
        print(f'PASS {dim}D {mode}: balance and exact physical restart; '+
              f'{len(payloads["legacy_output_clock_rounding"])} legacy output-clock roundings',flush=True)
    invalid={
        'zero-count':('injection','num_points','0'),
        'outside':('injection','points_x','[2,.5]'),
        'bad-unit':('injection','points_unit','feet'),
        'nan-rate':('injection','rate','[nan]'),
        'bad-duration':('injection','duration_in_yr','[0]'),
        'nan-time':('injection','start_time_in_yr','[nan]'),
        'bad-rate-model':('injection','rate_model','unknown'),
        'negative-heat-capacity':('mat','heat_capacity','[-1]')}
    for name,(section,key,value) in invalid.items():
        cfg=config()
        if dim==2:
            for group in ('injection','heat_injection'):del cfg[group]['points_y']
        if name=='bad-duration':del cfg['injection']['end_time_in_yr']
        cfg[section][key]=value
        run(exe,out/name,cfg,env,11);result['checks'][name]=True
    for source in ('injection','heat_injection'):
        cfg=config()
        if dim==2:
            for group in ('injection','heat_injection'):del cfg[group]['points_y']
        cfg[source]['points_x']='[0,0]'
        if source=='injection':cfg['bc'].update(hbc_x0='2',hbc_pressure='[0,0,0,0,0,0]')
        else:cfg['bc'].update(tbc_x0='1',tbc_val_x0='300')
        run(exe,out/(source+'-fixed-node'),cfg,env,51)
        result['checks'][source+'-fixed-node']=True
    result['checks']['remesh']=material(exe,out,dim,env)
    assert digest(exe)==result['executable_sha256']
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS: {len(result["checks"])} point-source cases in {dim}D')


if __name__=='__main__':main()
