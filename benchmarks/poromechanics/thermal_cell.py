#!/usr/bin/env python3
"""Uniform heating of sealed fixed-strain, sealed free-stress and drained cells.

Nodal heat sources are proportional to independently integrated nodal volume,
so temperature stays uniform without modifying a saved field. Separate heat,
water and thermoelastic references distinguish the fixed-capacity model from
the retained-fluid enthalpy balance. The drained case records the latter's
known defect explicitly; it is not an energy qualification of open flow.
"""

import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

import fluid_content as water
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

ROOT = water.ROOT
T0, RATE, STOP, END = 273., 40., .25, 1.
K, G, ALPHA, PHI, KF, KS = 2e7, 1e7, .98, .25, 1e8, 1e9
EXPANSION = 2e-4
STORAGE = PHI/KF+(ALPHA-PHI)/KS
BETA = PHI*EXPANSION
FLUID = 1000.*4180.
CAPACITY = (1-PHI)*2650.*1000.+PHI*FLUID


sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import run_solver

def configuration(method, boundary, dt):
    cfg = water.configuration(.5, dt, END, method, 1e5, water.INFLOW)
    cfg['sim'].update(max_steps=str(round(END/dt)), max_time_in_yr='1e9')
    cfg['control'].update(fixed_dt=repr(dt), has_thermal_pressurization='true',
                          has_thermal_diffusion='true', has_fluid_thermal_advection='true')
    if boundary == 'drained':
        # A material cell with every vertex drained: no unresolved pressure
        # boundary layer is used as an approximation to the uniform reference.
        cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        cfg['bc'].update({f'hbc_{side}': '2' if boundary == 'drained' else '0',
                          f'tbc_{side}': '0', f'vbc_{side}': '3'})
        if boundary == 'free':
            cfg['bc'][f'vbc_{side}'] = '1' if side.endswith('0') else '0'
    cfg['bc'].update(hbc_flux='[0,0,0,0,0,0]', stress_bc_z1='0',
                     surface_temperature=repr(T0), mantle_temperature=repr(T0))
    cfg['ic'].update(temperature_option='4', temp_background=repr(T0), temp_anomaly_amplitude='0')
    cfg['mat'].update(bulk_modulus=f'[{K}]', shear_modulus=f'[{G}]', biot_coeff=f'[{ALPHA}]',
                      porosity=f'[{PHI}]', fluid_bulk_modulus=f'[{KF}]', bulk_modulus_s=f'[{KS}]',
                      hydraulic_perm='[0]', alpha='[0]', fluid_alpha=f'[{EXPANSION}]',
                      heat_capacity='[1000]', fluid_heat_capacity='[4180]', therm_cond='[0]')
    return cfg


def run(exe, folder, cfg, name, env):
    cfg['sim']['modelname'] = name
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    run_solver([str(exe), f'{name}.cfg'], log=folder / f'{name}.log', cwd=folder, env=env, check=True, timeout=600)
    return Dynearthsol(str(folder/name))


def volumes(model):
    x, conn = (model.read_field(model.frames[0], key) for key in ('coordinate', 'connectivity'))
    v = np.abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(int(model.ndims))
    lumped = np.zeros(len(x))
    np.add.at(lumped, conn, np.broadcast_to((v/conn.shape[1])[:, None], conn.shape))
    return x, conn, v, lumped


def heating(cfg, x, lumped):
    source = dict(enabled='true', num_points=str(len(x)), points_unit='m',
                  points_x=json.dumps(x[:, 0].tolist()), points_z=json.dumps(x[:, -1].tolist()),
                  rate_model='constant_rate', rate=json.dumps((RATE*CAPACITY*lumped).tolist()),
                  start_time_in_yr=json.dumps([0.]*len(x)),
                  duration_in_yr=json.dumps([STOP/YEAR2SEC]*len(x)))
    if x.shape[1] == 3:
        source['points_y'] = json.dumps(x[:, 1].tolist())
    cfg['heat_injection'] = source


def analyze(folder, method, boundary, modelname='heated'):
    model = Dynearthsol(str(folder/modelname))
    rows = scan_frames(str(folder/modelname))
    x, conn, v, lumped = volumes(model)
    dim, volume = int(model.ndims), float(v.sum())
    bulk_mean = K+(G/3 if dim == 2 else 0)
    lag = ALPHA**2/bulk_mean if method == 'pt' else 0.
    pressure_scale = BETA*RATE*STOP/STORAGE
    history, water_boundary, carried_heat = [], 0., 0.
    previous_mean, reference_pressure = T0, 0.
    maxima = dict(temperature=0., heat=0., water=0., pressure_discrete_reference=0., stress=0.)
    for row in rows:
        frame, time = row['frame'], row['time']
        theta = RATE*min(time, STOP)
        temperature = model.read_field(frame, 'temperature')
        p = model.read_field(frame, 'pore pressure')
        strain = model.read_field(frame, 'strain')[:, :dim].sum(axis=1)
        stress = model.read_field(frame, 'stress')[:, :dim]
        mean_temperature = float(lumped @ temperature/volume)
        mean_pressure = float(v @ p[conn].mean(axis=1)/volume)
        strain_integral = float(v @ strain)
        pressure_integral = mean_pressure*volume
        physical = STORAGE*pressure_integral+ALPHA*strain_integral-BETA*float(lumped @ (temperature-T0))
        dp = (water.checkpoint_field(model, frame, 'pending pore pressure increment', len(x))
              if method == 'pt' else np.zeros(len(x)))
        pending = -lag*float(v @ dp[conn].mean(axis=1))
        b = model.read_field(frame, 'fluid balance')
        water_boundary += b[3]
        carried_heat -= FLUID*.5*(previous_mean+mean_temperature)*b[3]
        energy = CAPACITY*float(lumped @ (temperature-T0))
        supplied = CAPACITY*volume*theta
        if boundary == 'fixed':
            exact_pressure = BETA*theta/STORAGE
            if frame:
                reference_pressure = (BETA*theta+lag*reference_pressure)/(STORAGE+lag)
            expected_stress = -ALPHA*(mean_pressure + float(v @ dp[conn].mean(axis=1)/volume)
                                      if method == 'pt' else mean_pressure)
        elif boundary == 'free':
            exact_pressure = BETA*theta/(STORAGE+ALPHA**2/bulk_mean)
            reference_pressure = exact_pressure
            expected_stress = 0.
        else:
            exact_pressure = reference_pressure = expected_stress = 0.
        maxima['temperature'] = max(maxima['temperature'], float(np.max(np.abs(temperature-(T0+theta))))/(RATE*STOP))
        maxima['heat'] = max(maxima['heat'], abs(energy-supplied)/(CAPACITY*volume*RATE*STOP))
        maxima['water'] = max(maxima['water'], abs(physical+pending-water_boundary)/(BETA*volume*RATE*STOP))
        maxima['pressure_discrete_reference'] = max(maxima['pressure_discrete_reference'],
                                                   float(np.max(np.abs(p-reference_pressure)))/pressure_scale)
        maxima['stress'] = max(maxima['stress'], float(np.max(np.abs(stress-expected_stress)))/pressure_scale)
        # Reference masses plus the measured linear Biot fluid-content change.
        # No pressure work or skeleton work enters this rigid, p=0 drain cell.
        retained_energy_change = energy+FLUID*float(np.sum(v*(STORAGE*p[conn].mean(axis=1)+ALPHA*strain-
                                                           BETA*(temperature[conn].mean(axis=1)-T0))*
                                                        temperature[conn].mean(axis=1)))
        enthalpy_residual = retained_energy_change+carried_heat-supplied
        history.append([time, mean_temperature, mean_pressure, exact_pressure, physical,
                        pending, water_boundary, energy, supplied, carried_heat, enthalpy_residual])
        previous_mean = mean_temperature
    for key, error in maxima.items():
        assert error < 2e-8, (boundary, method, key, error)
    history = np.asarray(history)
    sample = history[int(np.argmin(abs(history[:, 0]-.2)))]
    np.savetxt(folder/'history.csv', history, delimiter=',', comments='',
               header='time_s,temperature_K,pressure_Pa,physical_reference_Pa,physical_content_m3,pending_split_m3,boundary_in_m3,fixed_capacity_heat_J,supplied_J,outflow_enthalpy_J,retained_enthalpy_residual_J')
    report = dict(method=method, boundary=boundary, dimension=dim, frames=len(rows), nodes=len(x),
                  maximum_errors=maxima, final_pressure_Pa=float(history[-1, 2]),
                  final_pressure_error_over_reference=(abs(history[-1, 2]/history[-1, 3]-1) if boundary != 'drained' else None),
                  heating_pressure_error_at_point2_s=(abs(sample[2]/sample[3]-1) if boundary != 'drained' else None),
                  fixed_capacity_heat_error=maxima['heat'],
                  retained_enthalpy_residual_J=float(history[-1, 10]) if boundary == 'drained' else None,
                  retained_enthalpy_relative_error=float(abs(history[-1, 10])/history[-1, 8]) if boundary == 'drained' else None,
                  drained_defect_reference_J=-.5*FLUID*BETA*(RATE*STOP)**2*volume if boundary == 'drained' else None,
                  full_open_energy_qualified=False)
    if boundary == 'drained':
        np.testing.assert_allclose(report['retained_enthalpy_residual_J'], report['drained_defect_reference_J'], rtol=1e-7)
    water.write_json(folder/'analysis.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--methods', nargs='+', choices=('pt', 'explicit'), default=['pt', 'explicit'])
    parser.add_argument('--boundaries', nargs='+', choices=('fixed', 'free', 'drained'), default=['fixed', 'free', 'drained'])
    parser.add_argument('--dt', type=float, default=.01)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if not 0 < args.dt <= END:
        parser.error('--dt must be positive and no larger than one second')
    exe, output = args.executable.resolve(), args.output.resolve()
    if args.replay:
        manifest = json.loads((output/'manifest.json').read_text())
        for name, sha in manifest['artifacts_sha256'].items():
            assert water.digest(output/name) == sha, name
    else:
        output.mkdir(parents=True, exist_ok=False)
        manifest = dict(executable=str(exe), executable_sha256=water.digest(exe), command=sys.argv,
                        runner_sha256=water.digest(__file__), cases=[])
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        for method in args.methods:
            for boundary in args.boundaries:
                if method == 'explicit' and boundary == 'free':
                    # A free specimen has an inertial transient; it cannot use
                    # the quasistatic pressure reference without equilibration.
                    continue
                name = f'{method}-{boundary}'; folder = output/name; folder.mkdir()
                cfg = configuration(method, boundary, args.dt)
                steps = cfg['sim']['max_steps']; cfg['sim']['max_steps'] = '1'
                seed = run(exe, folder, cfg, 'seed', env)
                assert int(seed.ndims) == args.ndims
                x, _, _, lumped = volumes(seed)
                heating(cfg, x, lumped); cfg['sim']['max_steps'] = steps
                run(exe, folder, cfg, 'heated', env)
                assert water.digest(exe) == manifest['executable_sha256']
                manifest['cases'].append(dict(folder=name, method=method, boundary=boundary))
                # Retain the raw failure evidence if an analytical check fails.
                manifest['artifacts_sha256'] = {str(p.relative_to(output)): water.digest(p)
                                                for p in sorted(output.glob('*/*')) if p.is_file()
                                                and p.name not in ('analysis.json', 'history.csv')}
                water.write_json(output/'manifest.json', manifest)
                print(json.dumps(analyze(folder, method, boundary)), flush=True)
    summary = dict(cases=[analyze(output/c['folder'], c['method'], c['boundary']) for c in manifest['cases']],
                   full_open_energy_qualified=False)
    water.write_json(output/'summary.json', summary)


if __name__ == '__main__':
    main()
