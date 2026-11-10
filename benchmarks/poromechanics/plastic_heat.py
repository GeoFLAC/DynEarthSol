#!/usr/bin/env python3
"""Adiabatic plastic heating in a prescribed, spin-free pure-shear cell.

The independent scalar reference is tau=min(G*gamma,c) and
Q=c*max(gamma-c/G,0), where gamma=2*rate*t. This is a material-work
qualification for earthquake heating, not a simulated earthquake cycle.
Stored heat plus the undeposited energy must equal that plastic work.
"""

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

import fluid_content as water
from Dynearthsol import Dynearthsol, scan_frames

RATE, BULK, SHEAR = .5, 2e10, 1e9
YIELD, CONFINING, CAPACITY, T0 = 1e7, 3e7, 2700*1000., 300.


def work(time, elastic=False):
    gamma = 2*RATE*np.asarray(time)
    return np.zeros_like(gamma) if elastic else YIELD*np.maximum(gamma-YIELD/SHEAR, 0)


def rsf_continuum(times, theta0, width):
    """Independent aging-law solution and work integral for this shear cell.

    The strength decreases after first yield for the stated initial state.
    Quadrature integrates boundary work; subtracting elastic storage gives
    irreversible work. This is a material path, not a fault rupture solution.
    """
    from scipy.integrate import quad
    from scipy.optimize import brentq

    speed = 2*RATE*width
    equilibrium = .01/speed

    def strength(time):
        theta = equilibrium+(theta0-equilibrium)*np.exp(-speed*time/.01)
        mu = np.tan(np.pi/6)+.01*np.log(speed/.001)+.015*np.log(.001*theta/.01)
        return (YIELD+CONFINING*mu)/np.sqrt(1+mu*mu)

    yield_time = brentq(lambda t: 2*SHEAR*RATE*t-strength(t), 0., 1., xtol=1e-15)
    times = np.asarray(times)
    tau = np.minimum(2*SHEAR*RATE*times, strength(times))
    energy = np.array([0. if t <= yield_time else
        2*RATE*quad(strength, yield_time, t, epsabs=1e-6, epsrel=1e-12)[0]
        -(strength(t)**2-strength(yield_time)**2)/(2*SHEAR) for t in times])
    return tau, energy


def configuration(dim, method, dt, steps, elastic=False):
    cfg = water.configuration(.5, dt, 1., 'pt' if method.startswith('pt') else method, 1e5, 0.)
    cfg['sim'].update(max_steps=str(steps), max_time_in_yr='1e9')
    cfg['mesh'].update(xlength='1', ylength='1', zlength='1', resolution='1',
                       largest_size='1', smallest_size='1')
    cfg['control'].update(fixed_dt=repr(dt), has_hydraulic_diffusion='false',
                          has_thermal_diffusion='false', has_shear_heating='true',
                          has_fluid_thermal_advection='false', is_using_mixed_stress='false',
                          PT_option=method[-1] if method.startswith('pt') else '1')
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        cfg['bc'].update({f'vbc_{side}': '1', f'vbc_val_{side}': '0',
                          f'hbc_{side}': '0', f'tbc_{side}': '0',
                          f'stress_bc_{side}': '0'})
    cfg['bc'].update(vbc_val_x1=repr(RATE), vbc_val_z0=repr(RATE),
                     hbc_flux='[0,0,0,0,0,0]', surface_temperature=repr(T0),
                     mantle_temperature=repr(T0))
    cfg['ic'].update(temperature_option='4', temp_background=repr(T0), temp_anomaly_amplitude='0',
                     initial_stress_option='1',
                     initial_stress=repr([-CONFINING]*dim+[0.]*(3 if dim == 3 else 1)))
    cfg['mat'].update(rheology_type='elastic' if elastic else 'elasto-plastic',
                      is_plane_strain=str(dim == 2).lower(), rho0='[2700]', heat_capacity='[1000]',
                      bulk_modulus=repr([BULK]), shear_modulus=repr([SHEAR]),
                      alpha='[0]', fluid_alpha='[0]', therm_cond='[0]',
                      cohesion0=repr([YIELD]), cohesion1=repr([YIELD]),
                      friction_angle0='[0]', friction_angle1='[0]',
                      dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e12')
    return cfg


def run(exe, folder, cfg, name='cell', timeout=180):
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    with (folder/f'{name}.log').open('w') as stream:
        subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, stdout=stream,
                       stderr=subprocess.STDOUT, timeout=timeout, check=True,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    return Dynearthsol(str(folder/name))


def analyze(model, elastic=False, initial_generated=0., reference=None, capacity=CAPACITY):
    records = scan_frames(model.modelname)
    dim = int(model.ndims)
    x, conn = (model.read_field(records[0]['frame'], key) for key in ('coordinate', 'connectivity'))
    volume = abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    total_volume = volume.sum()
    generated = initial_generated
    times = np.array([row['time'] for row in records])
    targets = work(times, elastic) if reference is None else reference(times)
    worst_balance = worst_work = 0.
    rows = []
    for i, row in enumerate(records):
        frame, time = row['frame'], row['time']
        pending = model.read_field(frame, 'pending shear heat').sum()
        power = model.read_field(frame, 'shear heating') @ volume
        temperature = model.read_field(frame, 'temperature')
        stored = capacity*(temperature[conn].mean(axis=1)-T0) @ volume
        if i:
            assert row['steps'] == records[i-1]['steps']+1
            generated += power*(time-records[i-1]['time'])
        target = float(targets[i])*total_volume
        scale = YIELD*total_volume
        worst_balance = max(worst_balance, abs(stored+pending-generated)/scale)
        worst_work = max(worst_work, abs(generated-target)/scale)
        rows.append(dict(time=time, step=row['steps'], stored=stored, pending=float(pending),
                         generated=float(generated), reference=target, power=float(power)))
    assert worst_balance < 2e-11, ('heat+pending != accepted work', worst_balance)
    assert worst_work < 2e-11, ('plastic work reference', worst_work)
    if elastic:
        assert all(r['power'] == 0 and r['pending'] == 0 and r['stored'] == 0 for r in rows)
    return dict(dimensions=dim, nodes=len(x), elements=len(conn), saved_frames=len(rows),
                balance_error_over_yield_volume=worst_balance,
                work_error_over_yield_volume=worst_work, rows=rows)


def continuation(base, resumed, extra_fields=()):
    original = {r['steps']: r for r in scan_frames(base.modelname)}
    fields = ('coordinate', 'connectivity', 'temperature', 'velocity', 'strain', 'stress',
              'plastic strain', 'plastic strain-rate', 'viscosity', 'shear heating', 'pending shear heat') + tuple(extra_fields)
    count = 0
    for row in scan_frames(resumed.modelname):
        old = original[row['steps']]
        assert old['time'] == row['time']
        for key in fields:
            np.testing.assert_array_equal(base.read_field(old['frame'], key), resumed.read_field(row['frame'], key))
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--dimensions', type=int, choices=(2, 3), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--methods', nargs='+', default=['pt0', 'pt1', 'explicit', 'gvs'])
    parser.add_argument('--dt', type=float, default=.005)
    parser.add_argument('--steps', type=int, default=40)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    results = []
    for method in args.methods:
        for elastic in (False, True):
            name = method+('-elastic' if elastic else '-plastic')
            folder = args.output/name
            folder.mkdir()
            cfg = configuration(args.dimensions, method, args.dt, args.steps, elastic)
            base = run(args.exe.resolve(), folder, cfg)
            result = analyze(base, elastic)
            result.update(case=name)
            if not elastic:
                restart = copy.deepcopy(cfg)
                middle = args.steps//2
                restart['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'cell'),
                                      restarting_from_frame=str(middle))
                resumed = run(args.exe.resolve(), folder, restart, 'resumed')
                result['exact_continuation_payloads'] = continuation(base, resumed)
                restart['control']['fixed_dt'] = repr(args.dt*.4)
                changed = run(args.exe.resolve(), folder, restart, 'changed-step')
                result['changed_step'] = analyze(changed, initial_generated=result['rows'][middle]['generated'])
            results.append(result)
            print(name, result['balance_error_over_yield_volume'], result['work_error_over_yield_volume'], flush=True)
            (args.output/'summary.json').write_text(json.dumps(results, indent=2)+'\n')
    shutil.copy2(__file__, args.output/'runner.py')
    (args.output/'runner.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()+'\n')


if __name__ == '__main__':
    main()
