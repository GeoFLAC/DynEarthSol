#!/usr/bin/env python3
"""Finite-density water, sensible heat and remap qualification.

The physical inventory M and numerical fixed-stress buffer R are distinct:
the split closes M+R, while sensible heat is (Cs+cf*M)*T. References use a
constant-compressibility EOS and Eulerian Biot porosity, with explicit
reservoir mass/enthalpy reactions and geometric loss/import at remeshing.

    check_mass.py EXECUTABLE OUTPUT [--cases mass heat remap lag transport]
"""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT),str(ROOT/'benchmarks/poromechanics'),str(ROOT/'tests/functional')]
import check_plastic_content as cell_tests
import check_fluid_remap as remap_tests
from mass_reference import Cell, Transport
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames


def write_summary(output, result):
    (output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def configuration():
    cfg = cell_tests.configuration('pt',0,elastic=True)
    cfg['sim'].update(max_time_in_yr='1e9',has_marker_output='false')
    cfg['mesh'].update(xlength='1',ylength='1',zlength='1',resolution='1',smallest_size='1',largest_size='1')
    cfg['control'].update(hydraulic_content_model='1',has_thermal_diffusion='false',
                          has_fluid_thermal_advection='false',is_using_mixed_stress='false',inertial_scaling='1e5')
    cfg['bc'].update(vbc_val_z1='0',surface_temperature='273',mantle_temperature='273')
    for side in ('x0','x1','y0','y1','z0','z1'):
        cfg['bc'][f'tbc_{side}'] = '0'
    cfg['ic'].update(temperature_option='4',temp_background='273',temp_anomaly_amplitude='0',
                     weakzone_option='0',has_body_force_adjustment='false')
    cfg['mat'].update(porosity_model='1',porosity='[.25]',biot_coeff='[.98]',rho0='[2650]',
                      bulk_modulus='[2e7]',shear_modulus='[1e7]',bulk_modulus_s='[1e9]',fluid_bulk_modulus='[1e8]',
                      alpha='[0]',fluid_alpha='[.0002]',heat_capacity='[1000]',fluid_heat_capacity='[4180]',
                      therm_cond='[0]',hydraulic_perm='[0]')
    return cfg


def lag_checks(exe, out, timeout=600):
    out.mkdir(parents=True,exist_ok=False)
    reference = Cell()
    reports = []
    for steps in (1000,2000,4000):
        folder = out/str(steps);folder.mkdir()
        cfg = configuration()
        cfg['sim'].update(max_steps=str(steps),output_step_interval=str(steps//10),checkpoint_frame_interval='100000')
        cfg['control']['fixed_dt'] = repr(.1/steps)
        cfg['bc']['vbc_val_z1'] = '-.01'
        model = cell_tests.run(exe,folder,cfg,timeout=timeout)
        history = []
        for frame,time in zip(model.frames,model.time):
            mass = float(model.read_field(frame,'fluid mass').sum())
            buffer = float(model.read_field(frame,'hydraulic split buffer').sum())
            pressure = float(model.read_field(frame,'pore pressure').mean())
            exact = reference.sealed_pressure(273.,-.01*time)
            assert abs(mass+buffer-reference.initial_mass)<2e-9
            history.append([time,pressure,exact,mass,buffer,reference.initial_mass])
        time,pressure,exact,mass,buffer,total = history[-1]
        error = abs(pressure/exact-1)
        reports.append(dict(steps=steps,dt=.1/steps,pressure_error=error,
                            physical_mass_defect_kg=mass-total,numerical_buffer_kg=buffer))
        np.savetxt(folder/'history.csv',history,delimiter=',',comments='',
                   header='time_s,p_Pa,exact_p_Pa,M_kg,R_kg,initial_M_kg')
    errors = [row['pressure_error'] for row in reports]
    assert all(1.8<errors[i]/errors[i+1]<2.2 for i in range(len(errors)-1)),errors
    assert errors[-1]<.005,errors
    return write_summary(out,dict(passed=True,cases=reports,
        scope='First-order fixed-stress splitting error; exact budget is M+R, not physical M alone'))


def shear_checks(exe, out, timeout=600):
    import plastic_heat as plastic
    out.mkdir(parents=True,exist_ok=False)
    seed = out/'seed';seed.mkdir()
    initial = configuration();initial['sim']['max_steps']='0'
    dimensions = cell_tests.run(exe,seed,initial,timeout=timeout).ndims
    cfg = plastic.configuration(dimensions,'pt1',.0025,40)
    cfg['control'].update(has_hydraulic_diffusion='true',hydraulic_content_model='1',
                          has_thermal_diffusion='true',has_fluid_thermal_advection='true',
                          has_thermal_pressurization='true',has_poroelastic_pressure_feedback='true')
    cfg['mat'].update(porosity_model='1',hydraulic_perm='[0]',porosity='[.25]',
                      fluid_rho0='[1000]',fluid_heat_capacity='[4180]',bulk_modulus_s='[1e20]')
    model = plastic.run(exe,out,cfg,timeout=timeout)
    def inspect(model, initial_energy=None):
        records = []
        for frame,time in zip(model.frames,model.time):
            mass = model.read_field(frame,'fluid mass')
            temperature = model.read_field(frame,'temperature')
            energy = float(np.sum((model.read_field(frame,'solid heat capacity')+4180*mass)*temperature))
            if initial_energy is None:initial_energy=energy
            expected = float(plastic.work(time))
            assert abs(energy-initial_energy-expected)<3e-5
            assert np.max(abs(model.read_field(frame,'pore pressure')))<1e-4
            assert np.max(abs(model.read_field(frame,'pending shear heat')))==0
            records.append([time,energy-initial_energy,expected])
        return records,initial_energy
    history,initial_energy = inspect(model)
    restart = copy.deepcopy(cfg)
    restart['sim'].update(is_restarting='true',restarting_from_modelname=model.modelname,
                          restarting_from_frame='20',max_steps='40')
    resumed = plastic.run(exe,out,restart,'resumed',timeout=timeout)
    exact = cell_tests.compare_continuation(out,model.modelname,resumed)
    restart['control']['fixed_dt'] = '.005'
    restart['sim']['max_steps'] = '30'
    changed = plastic.run(exe,out,restart,'changed',timeout=timeout)
    inspect(changed,initial_energy)
    np.savetxt(out/'history.csv',history,delimiter=',',comments='',header='time_s,stored_heat_J,plastic_work_J')
    return write_summary(out,dict(passed=True,exact_restart_payloads=exact,changed_dt_pass=True))


def transport_checks(exe, out, timeout=600):
    out.mkdir(parents=True,exist_ok=False)
    seed = out/'seed';seed.mkdir()
    cfg = configuration();cfg['sim']['max_steps']='0'
    dim = int(cell_tests.run(exe,seed,cfg,timeout=timeout).ndims)
    def setup(dt, hot, balanced=0.):
        cfg = configuration()
        cfg['sim'].update(max_steps='1',output_step_interval='1')
        cfg['mesh'].update(resolution='.5',smallest_size='.5',largest_size='.5')
        cfg['control'].update(fixed_dt=repr(dt),has_thermal_diffusion='true',
                              has_fluid_thermal_advection='true',has_thermal_pressurization='true')
        cfg['bc'].update(hbc_x0='2',hbc_x1='2',hbc_pressure='[20000,80000,0,0,0,0]')
        cfg['ic'].update(pore_pressure_option='2',pore_pressure_value='50000',temp_background='298',
                         temp_anomaly_amplitude='10',temp_anomaly_xcenter='.5',temp_anomaly_ycenter='.5',
                         temp_anomaly_zcenter='-.5',temp_anomaly_xsemi_axis='1',temp_anomaly_ysemi_axis='1',
                         temp_anomaly_zsemi_axis='1')
        cfg['mat'].update(hydraulic_perm='[1e-10]',permeability_tensor='[2,1.2,1.5,.1,.3,.2]',therm_cond='[2]')
        points = [[.25]*(dim-1)+[-.25],[.75]*(dim-1)+[-.75]]
        rates, starts, stops = [.001,-.0005],[.13*dt,0.],[.47*dt,dt]
        if balanced:
            points += [[.5]*(dim-1)+[-.5]]*2
            rates += [balanced,-balanced];starts += [0.,0.];stops += [dt,dt]
        source = dict(enabled='true',num_points=str(len(points)),points_unit='m',rate_model='constant_rate',
            points_x=json.dumps([p[0] for p in points]),points_z=json.dumps([p[-1] for p in points]),
            rate=json.dumps(rates),start_time_in_yr=json.dumps((np.array(starts)/YEAR2SEC).tolist()),
            end_time_in_yr=json.dumps((np.array(stops)/YEAR2SEC).tolist()))
        if dim==3:source['points_y']=json.dumps([p[1] for p in points])
        if hot:source['temperature']='[333]'
        cfg['injection']=source
        return cfg
    folder = out/'flux';folder.mkdir()
    model = cell_tests.run(exe,folder,setup(1e-5,True),timeout=timeout)
    reference = Transport(model)
    frame = model.frames[-1]
    p,t = (model.read_field(frame,k) for k in ('pore pressure','temperature'))
    mass,mp,mt,moving = reference.state(p,t,1.)
    rm,re = reference.rates(reference.p0,reference.t0,0.,reference.duration,True)
    m0 = model.read_field(0,'fluid mass')
    e0 = (reference.cs+4180*m0)*reference.t0
    predicted_mass = m0+reference.duration*rm
    predicted_heat = e0+reference.duration*re
    reservoir = np.where(reference.fixed,mass-predicted_mass,0.)
    mass_error = float(np.max(np.abs((mass+reference.lag*(p-reference.p0)-predicted_mass)[~reference.fixed])))
    energy_error = float(np.max(abs((reference.cs+4180*mass)*t-predicted_heat-4180*t*reservoir)))
    np.testing.assert_allclose(mass,model.read_field(frame,'fluid mass'),rtol=2e-13,atol=2e-11)
    assert mass_error<2e-10,mass_error
    assert energy_error<3e-5,energy_error
    errors = []
    # Co-located equal inflow/outflow at the local temperature leaves both
    # physical right-hand sides unchanged, while tightening the CFL bound.
    for label,balanced in (('subcycle',0.),('refined',200.)):
        folder = out/label;folder.mkdir()
        model = cell_tests.run(exe,folder,setup(1.,False,balanced),timeout=timeout)
        reference = Transport(model)
        expected_p,expected_t = reference.solve()
        actual_p,actual_t = (model.read_field(model.frames[-1],k) for k in ('pore pressure','temperature'))
        error = max(float(np.max(abs(actual_p-expected_p)))/60000.,float(np.max(abs(actual_t-expected_t)))/50.)
        errors.append(error)
        np.savetxt(folder/'profile.csv',np.column_stack((reference.x,actual_p,expected_p,actual_t,expected_t)),
                   delimiter=',',comments='',header=','.join((['x','z'] if dim==2 else ['x','y','z'])+
                                                           ['p_Pa','reference_p_Pa','T_K','reference_T_K']))
        assert np.all(np.isfinite(actual_p)) and np.all(actual_t>0)
    assert errors[0]<.01 and errors[1]<.5*errors[0],errors
    return write_summary(out,dict(passed=True,flux_mass_error_kg=mass_error,flux_energy_error_J=energy_error,
        time_errors=errors,scope='Accepted-path spatial transport and explicit CFL refinement; independent Radau reference'))

def mass_checks(exe, out, timeout=600):
    out.mkdir(parents=True,exist_ok=False)
    c = cell_tests
    ref=Cell(reference_temperature=273.)
    def config(method,drained):
        cfg=c.configuration(method,0,elastic=True)
        cfg['sim'].update(max_steps='100',has_marker_output='true')
        cfg['mesh'].update(resolution='1',smallest_size='1',largest_size='1')
        cfg['control'].update(hydraulic_content_model='1',fixed_dt='.01',is_using_mixed_stress='false',use_global_velocity_scaling=str(method=='gvs'))
        cfg['bc'].update(vbc_val_z1='-.01',surface_temperature='273',mantle_temperature='273')
        if drained:
            for s in ('x0','x1','y0','y1','z0','z1'):cfg['bc'][f'hbc_{s}']='2'
        cfg['mat'].update(porosity_model='1',biot_coeff='[.98]',bulk_modulus_s='[1e9]',fluid_bulk_modulus='[1e8]')
        cfg['ic'].update(temperature_option='4',temp_background='273',temp_anomaly_amplitude='0')
        return cfg
    reports=[]
    for method,drained,scaling in [('pt',False,1e5),('pt',True,1e5),('explicit',False,1e4),('gvs',False,1e6)]:
        folder=out/(method+('-drained' if drained else '-sealed'));folder.mkdir();cfg=config(method,drained);cfg['control']['inertial_scaling']=repr(scaling)
        m=c.run(exe,folder,cfg,timeout=timeout);dim=int(m.ndims);bulk=2e7+(1e7/3 if dim==2 else 0)
        pprev=0.;epsprev=0.;boundary=0.;worst=0.;budget=0.;history=[]
        for f,t in zip(m.frames,m.time):
            pressure=m.read_field(f,'pore pressure');mass=m.read_field(f,'fluid mass').sum();buf=m.read_field(f,'hydraulic split buffer').sum();balance=m.read_field(f,'fluid mass balance');boundary+=balance[3];eps=-.01*t
            budget=max(budget,abs(mass+buf-ref.initial_mass-boundary))
            if drained:
                expected=ref.mass(0,273,eps);worst=max(worst,abs(mass-expected));assert np.max(abs(pressure))==0
            elif f:
                lag=ref.density(pprev,273)*math.exp(epsprev)*.98**2/bulk if method=='pt' else 0.
                expected=brentq(lambda x:ref.mass(x,273,eps)+lag*(x-pprev)-ref.initial_mass,-1e7,1e7,xtol=1e-7)
                worst=max(worst,float(abs(pressure-expected).max()));pprev=expected;epsprev=eps
            np.testing.assert_allclose(mass,ref.mass(float(pressure.mean()),273,eps),rtol=2e-13,atol=2e-11)
            history.append([t,float(pressure.mean()),mass,buf,boundary])
        assert budget<2e-9,budget
        assert worst<(1e-8 if drained else 1e-4),worst
        resume=copy.deepcopy(cfg);resume['sim'].update(is_restarting='true',restarting_from_modelname=m.modelname,restarting_from_frame='50',max_steps='100')
        resumed=c.run(exe,folder,resume,'resumed',timeout=timeout);count=c.compare_continuation(folder,m.modelname,resumed)
        resume['control']['fixed_dt']='.025';resume['sim']['max_steps']='70';changed=c.run(exe,folder,resume,'changed',timeout=timeout)
        # A changed dt must retain the old extensive buffer without rescaling it.
        changed_boundary=sum(m.read_field(f,'fluid mass balance')[3] for f in m.frames if f<=50)
        for f in changed.frames[1:]:
            changed_boundary+=changed.read_field(f,'fluid mass balance')[3]
            actual=changed.read_field(f,'fluid mass').sum()+changed.read_field(f,'hydraulic split buffer').sum()
            assert abs(actual-ref.initial_mass-changed_boundary)<2e-9
        np.savetxt(folder/'history.csv',history,delimiter=',',comments='',header='time_s,mean_pressure_Pa,fluid_mass_kg,split_buffer_kg,boundary_in_kg')
        reports.append(dict(method=method,drained=drained,ndims=dim,scaling=scaling,budget_error_kg=float(budget),reference_error=float(worst),exact_restart_payloads=count,changed_dt_mass_pass=True))
        print('PASS',folder.name,flush=True)
    return write_summary(out,dict(passed=True,cases=reports))


def heat_checks(exe, out, timeout=600):
    out.mkdir(parents=True,exist_ok=False)
    c = cell_tests
    ref=Cell();duration=5.;power=3e7;inlet=323.;volume_rate=.01
    cfg=c.configuration('pt',0,elastic=True)
    cfg['sim'].update(max_steps='0',has_marker_output='true',max_time_in_yr='1e9')
    cfg['mesh'].update(resolution='1',smallest_size='1',largest_size='1')
    cfg['control'].update(hydraulic_content_model='1',fixed_dt='.05',is_using_mixed_stress='false',
        has_thermal_diffusion='true',has_fluid_thermal_advection='true',has_thermal_pressurization='true')
    cfg['bc'].update(vbc_val_z1='0',surface_temperature='273',mantle_temperature='273')
    for s in ('x0','x1','y0','y1','z0','z1'):cfg['bc'][f'tbc_{s}']='0'
    cfg['mat'].update(porosity_model='1',biot_coeff='[.98]',bulk_modulus_s='[1e9]',fluid_bulk_modulus='[1e8]',
        rho0='[2650]',alpha='[0]',fluid_alpha='[.0002]',heat_capacity='[1000]',fluid_heat_capacity='[4180]',therm_cond='[0]')
    cfg['ic'].update(temperature_option='4',temp_background='273',temp_anomaly_amplitude='0')
    seed=out/'mesh';seed.mkdir();m=c.run(exe,seed,cfg,timeout=timeout)
    x,conn=(m.read_field(0,key) for key in ('coordinate','connectivity'));dim=int(m.ndims)
    vol=np.abs(np.linalg.det(x[conn][:,1:]-x[conn][:,:1]))/math.factorial(dim)
    weights=np.zeros(len(x));np.add.at(weights,conn,np.broadcast_to(vol[:,None]/conn.shape[1],conn.shape))
    assert abs(vol.sum()-1)<1e-12
    def point_sources(strength):
        d=dict(enabled='true',num_points=str(len(x)),points_x=json.dumps(x[:,0].tolist()),points_z=json.dumps(x[:,-1].tolist()),
            points_unit='m',rate_model='constant_rate',rate=json.dumps((strength*weights).tolist()),start_time_in_yr='[0]')
        if dim==3:d['points_y']=json.dumps(x[:,1].tolist())
        return d
    reports=[]
    for kind in ('drained','sealed','mixing','expanding-mixing'):
        errors=[]
        for steps in (100,200):
            path=out/(kind+'-'+str(steps));path.mkdir();case=copy.deepcopy(cfg)
            case['sim']['max_steps']=str(steps);case['control']['fixed_dt']=repr(duration/steps)
            if kind!='sealed':
                for s in ('x0','x1','y0','y1','z0','z1'):case['bc'][f'hbc_{s}']='2'
            if 'mixing' in kind:
                case['injection']=point_sources(volume_rate)
                case['injection']['temperature']=f'[{inlet}]'
                expansion=0 if kind=='mixing' else .0002
                case['mat']['fluid_alpha']=repr([expansion])
                cell=Cell(fluid_expansion=expansion)
            else:
                case['heat_injection']=point_sources(power);cell=ref
            model=c.run(exe,path,case,timeout=timeout)
            supplied=well=carried=reaction=boundary=0.;worst=0.;pold=0.;told=273.
            history=[]
            for f,t in zip(model.frames,model.time):
                pn=model.read_field(f,'pore pressure');tn=model.read_field(f,'temperature')
                mass=float(model.read_field(f,'fluid mass').sum());buffer=float(model.read_field(f,'hydraulic split buffer').sum())
                cs=float(model.read_field(f,'solid heat capacity').sum())
                ledger=model.read_field(f,'sensible energy balance');water=model.read_field(f,'fluid mass balance')
                supplied+=ledger[1];well+=ledger[2];carried+=ledger[3];reaction+=ledger[4];boundary+=water[3]
                energy=float(np.sum((model.read_field(f,'solid heat capacity')+cell.fluid_heat*model.read_field(f,'fluid mass'))*tn))
                assert abs(energy-cell.energy(0.,273.)-supplied-well-carried-reaction)<2e-5
                np.testing.assert_allclose(mass,cell.mass(float(pn.mean()),float(tn.mean())),rtol=2e-12)
                np.testing.assert_allclose(tn,tn.mean(),rtol=0,atol=2e-8)
                if kind=='sealed' and f:
                    dt=duration/steps;bulk=2e7+(1e7/3 if dim==2 else 0)
                    lag=cell.density(pold,told)*.98**2/bulk
                    target=cell.energy(0.,273.)+power*t
                    solved=root(lambda v:[(cell.mass(v[0],v[1])+lag*(v[0]-pold)-cell.initial_mass)/cell.initial_mass,
                        (cell.energy(v[0],v[1])-target)/target],[pold,told])
                    assert solved.success,solved.message
                    pold,told=solved.x
                    worst=max(worst,abs(float(pn.mean())-pold)/1e5,abs(float(tn.mean())-told))
                if kind=='drained':
                    expected=brentq(lambda v:cell.drained_heat(v)-power*t,273.,400.,xtol=1e-10)
                    worst=max(worst,abs(float(tn.mean())-expected))
                if 'mixing' in kind:
                    expected=inlet-(inlet-273)*math.exp(-cell.fluid_heat*cell.rho0*volume_rate*t/(cell.solid_capacity+cell.fluid_heat*cell.initial_mass)) if kind=='mixing' else (273. if t==0 else brentq(lambda v:cell.injection_time(v,inlet,volume_rate)-t,273.,inlet-1e-8))
                    worst=max(worst,abs(float(tn.mean())-expected))
                history.append([t,float(pn.mean()),float(tn.mean()),mass,buffer,supplied,well,carried,reaction,energy])
            if kind=='sealed':assert worst<1e-5,worst
            else:assert worst<.1,worst
            resume=copy.deepcopy(case);resume['sim'].update(is_restarting='true',restarting_from_modelname=model.modelname,restarting_from_frame=str(steps//2))
            resumed=c.run(exe,path,resume,'resumed',timeout=timeout)
            count=c.compare_continuation(path,model.modelname,resumed)
            # Changing the prescribed reference boundary must not change the EOS
            # reference saved with the physical inventory. All actual T BCs are off.
            boundary_case=copy.deepcopy(resume);boundary_case['bc']['surface_temperature']='303'
            boundary_model=c.run(exe,path,boundary_case,'boundary',timeout=timeout)
            for frame in boundary_model.frames[1:]:
                for field in ('fluid mass','hydraulic split buffer','solid heat capacity','pore pressure','temperature','pending temperature increment'):
                    np.testing.assert_array_equal(model.read_field(frame,field),boundary_model.read_field(frame,field))
            incompatible=copy.deepcopy(resume);incompatible['sim']['modelname']='bad-law'
            incompatible['mat']['fluid_bulk_modulus']='[2e8]'
            with (path/'bad-law.cfg').open('w') as stream:incompatible.write(stream)
            with (path/'bad-law.log').open('w') as stream:
                failed=subprocess.run([str(exe),'bad-law.cfg'],cwd=path,env=dict(os.environ,OMP_NUM_THREADS='1'),stdout=stream,stderr=subprocess.STDOUT)
            assert failed.returncode==22,(path/'bad-law.log').read_text()
            np.savetxt(path/'history.csv',history,delimiter=',',comments='',header='time_s,p_Pa,T_K,M_kg,R_kg,Q_J,well_heat_J,boundary_heat_J,thermal_reaction_J,E_J')
            reports.append(dict(kind=kind,steps=steps,reference_error=worst,exact_restart_payloads=count))
            errors.append(worst);print('PASS',kind,steps,worst,flush=True)
        if kind!='sealed':assert errors[1]<.6*errors[0],errors
    return write_summary(out,dict(passed=True,cases=reports))


def remap_checks(exe, out, timeout=600):
    out.mkdir(parents=True,exist_ok=False)
    c = remap_tests
    reports=[]
    for kind in ('sealed','drained','thermal','explicit','gvs','domain-add','domain-remove'):
        folder=out/kind;folder.mkdir();cfg=c.configuration(kind)
        cfg['control'].update(hydraulic_content_model='1',inertial_scaling='1e5',surface_process_option='0',has_thermal_diffusion='false',has_fluid_thermal_advection='false')
        cfg['mat']['porosity_model']='1'
        for side in ('x0','x1','y0','y1','z0','z1'):cfg['bc'][f'tbc_{side}']='0'
        if kind=='thermal':cfg['control'].update(has_thermal_diffusion='true',has_fluid_thermal_advection='true')
        c.run(exe,folder,cfg,'input',timeout=timeout)
        model=Dynearthsol(str(folder/'flow'));rows=scan_frames(model.modelname)
        log=(folder/'input.log').read_text()
        by_frame={row['frame']:row for row in rows};pairs=[];last=None;before=None
        for event in re.finditer(r'Output #\s*(\d+)|Remeshing starts',log):
            if event.group(1) is None:
                before=last
            else:
                last=int(event.group(1))
                if before is not None:
                    pairs.append((by_frame[before],by_frame[last]));before=None
        assert pairs
        masses=re.findall(r'Fluid mass remap: removed/added/boundary \[kg\] = (\S+) (\S+) (\S+); numerical reset \[kg\] = (\S+)',log)
        heats=re.findall(r'Sensible heat remap: removed/added/carried/reaction \[J\] = (\S+) (\S+) (\S+) (\S+);',log)
        solids=re.findall(r'Solid remap: removed/added mass \[kg\] = (\S+) (\S+); capacity \[J/K\] = (\S+) (\S+)',log)
        buffers=re.findall(r'Numerical buffer remap: removed/added/reset \[kg\] = (\S+) (\S+) (\S+)',log)
        assert len(solids)==len(buffers)==len(masses)==len(heats)==len(pairs),(kind,len(masses),len(pairs))
        worst=0.;topologies=0
        def state(frame):
            mass=model.read_field(frame,'fluid mass');capacity=model.read_field(frame,'solid heat capacity');t=model.read_field(frame,'temperature')
            return dict(mass=float(mass.sum()),buffer=float(model.read_field(frame,'hydraulic split buffer').sum()),
                energy=float(np.sum((capacity+4180*mass)*t)),solid=float(model.read_field(frame,'solid mass').sum()),capacity=float(capacity.sum()),
                cells=model.read_field(frame,'connectivity'),nodes=model.read_field(frame,'coordinate'))
        for (before,after),mi,ei,si,bi in zip(pairs,masses,heats,solids,buffers):
            x,y=state(before['frame']),state(after['frame']);removed,added,boundary,reset=map(float,mi);eremoved,eadded,carried,reaction=map(float,ei)
            mass_error=abs(y['mass']-x['mass']+removed-added-boundary)/x['mass']
            heat_error=abs(y['energy']-x['energy']+eremoved-eadded-carried-reaction)/x['energy']
            assert mass_error<3e-12,(kind,mass_error)
            assert heat_error<3e-12,(kind,heat_error)
            sremoved,sadded,cremoved,cadded=map(float,si);bremoved,badded,breset=map(float,bi)
            assert abs(y['solid']-x['solid']+sremoved-sadded)<3e-12*x['solid']
            assert abs(y['capacity']-x['capacity']+cremoved-cadded)<3e-12*x['capacity']
            assert abs(y['buffer']-x['buffer']+bremoved-badded-breset)<1e-9
            worst=max(worst,mass_error,heat_error)
            topologies+=not np.array_equal(x['cells'],y['cells']) or not np.array_equal(x['nodes'],y['nodes'])
        assert topologies
        exact=0 if kind.startswith('domain') else c.continuation(exe,folder,cfg,rows)
        reports.append(dict(case=kind,remaps=len(pairs),changed_meshes=topologies,relative_inventory_defect=worst,exact_restart_payloads=exact))
        print('PASS',kind,worst,exact,flush=True)
    return write_summary(out,dict(passed=True,cases=reports))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--cases',nargs='+',choices=('mass','heat','remap','lag','shear','transport'),
                        default=['mass','heat','remap','lag','shear','transport'])
    parser.add_argument('--timeout',type=float,default=600,help='Wall-clock limit in seconds per solver case; does not change physical timesteps')
    args = parser.parse_args()
    if not (args.timeout>0 and math.isfinite(args.timeout)):parser.error('--timeout must be positive and finite')
    output = args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    checks = {'mass':mass_checks,'heat':heat_checks,'remap':remap_checks,
              'lag':lag_checks,'shear':shear_checks,'transport':transport_checks}
    results = {}
    for name in args.cases:
        results[name] = checks[name](args.executable.resolve(),output/name,timeout=args.timeout)
        print('PASS',name,flush=True)
    write_summary(output,dict(passed=True,checks=results))


if __name__ == '__main__':
    main()
