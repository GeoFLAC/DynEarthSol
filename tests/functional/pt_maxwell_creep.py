#!/usr/bin/env python3
"""Check confined Maxwell creep under fixed traction against a discrete reference.

The independent homogeneous reference solves prescribed axial stress with zero
lateral strain. Deviatoric stress uses Crank-Nicolson; mean stress is K*strain.
This tests continued physical creep, not just the first accepted load step.
"""
import argparse
import configparser
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Dynearthsol import Dynearthsol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('config', type=Path, help='Base dry elastic load example')
    parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), required=True)
    parser.add_argument('--option', type=int, choices=(0, 1), default=0)
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--timeout', type=float, default=240)
    args = parser.parse_args()
    if args.steps < 3 or args.timeout <= 0:
        parser.error('Require at least three steps and a positive timeout')
    executable, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config)
    cfg['sim'].update(modelname='creep', max_steps=str(args.steps),
                      output_step_interval='1', is_restarting='no')
    cfg['control'].update(PT_option=str(args.option), PT_use_running_scale='yes',
                          fixed_dt='0.01', has_moving_mesh='no', gravity='0')
    for name in ('PT_relative_tolerance', 'PT_absolute_tolerance'):
        cfg['control'].pop(name, None)
    cfg['mat'].update(rheology_type='maxwell', bulk_modulus='[1e7]', shear_modulus='[1e7]',
                      min_viscosity='1e6', max_viscosity='1e6',
                      is_plane_strain='yes' if args.ndims == 2 else 'no')
    with (output / 'input.cfg').open('w') as stream:
        cfg.write(stream)
    with (output / 'run.log').open('w') as stream:
        subprocess.run([str(executable), 'input.cfg'], cwd=output,
                       stdout=stream, stderr=subprocess.STDOUT, check=True,
                       timeout=args.timeout)
    model = Dynearthsol(str(output / 'creep'))
    assert model.ndims == args.ndims
    assert list(model.frames) == list(range(args.steps + 1))
    bulk = shear = 1e7
    viscosity, dt = 1e6, 0.01
    load = cfg.getfloat('bc', 'stress_val_z1')
    half_step = dt * shear / (2 * viscosity)
    decay = (1 - half_step) / (1 + half_step)
    deviatoric_modulus = 4 * shear / (3 * (1 + half_step))
    velocity_scale = abs(load) / (bulk + deviatoric_modulus) / dt
    mean = deviatoric = strain = 0.0
    maximum = dict(stress_zz=0.0, stress_xx=0.0, strain_zz=0.0, velocity=0.0)
    history = []
    for frame in model.frames[1:]:
        increment = (load - mean - decay * deviatoric) / (bulk + deviatoric_modulus)
        mean += bulk * increment
        deviatoric = decay * deviatoric + deviatoric_modulus * increment
        strain += increment
        stress = model.read_field(frame, 'stress')
        total_strain = model.read_field(frame, 'strain')
        velocity = model.read_field(frame, 'velocity')
        coordinate = model.read_field(frame, 'coordinate')
        assert all(np.isfinite(field).all() for field in (stress, total_strain, velocity, coordinate))
        reference_velocity = (coordinate[:, args.ndims-1] + 1) * increment / dt
        errors = dict(
            stress_zz=float(np.max(np.abs(stress[:, args.ndims-1] - load)) / abs(load)),
            stress_xx=float(np.max(np.abs(stress[:, 0] - (mean-deviatoric/2))) / abs(load)),
            strain_zz=float(np.max(np.abs(total_strain[:, args.ndims-1] - strain)) / (abs(load)/bulk)),
            velocity=float(np.max(np.abs(velocity[:, args.ndims-1] - reference_velocity)) / velocity_scale))
        assert np.isfinite(list(errors.values())).all(), (frame, errors)
        for name, error in errors.items():
            maximum[name] = max(maximum[name], error)
        history.append(dict(step=int(frame), strain_zz=strain,
                            stress_xx=mean-deviatoric/2, errors=errors))
    # Maxwell motion must persist beyond loading; a static-equilibrium shortcut
    # that ignores physical-time relaxation would fail both of these checks.
    assert np.max(np.abs(model.read_field(2, 'velocity'))) > 1e-4
    assert np.max(np.abs(model.read_field(args.steps, 'strain') - model.read_field(1, 'strain'))) > 1e-4
    assert np.isfinite(list(maximum.values())).all()
    (output / 'analytic.json').write_text(json.dumps(dict(max_errors=maximum, history=history), indent=2) + '\n')
    assert max(maximum.values()) <= 5e-5, maximum
    print(f'PASS: {args.ndims}D PT_option={args.option}, {args.steps} Maxwell creep steps; max error={max(maximum.values()):.6g}')


if __name__ == '__main__':
    main()
