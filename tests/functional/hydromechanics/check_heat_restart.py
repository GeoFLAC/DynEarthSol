#!/usr/bin/env python3
"""Check pending plastic heat, aging EP-RSF work and checkpoint guards."""

import argparse
import copy
import json
import itertools
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks/poromechanics'), str(ROOT/'tests/functional')]
import plastic_heat as heat
from pt_initial_checkpoint import fields

HEAT_STATE = ('shear heat version', 'pending shear heat', 'shear heating')


def checkpoint_copy(prefix, frame, target, drop=(), replace=None):
    """Modify a copied fixture, preserving every unrelated serialized field."""
    replace = replace or {}
    for suffix in ('info', 'manifest'):
        source = Path(f'{prefix}.{suffix}')
        if source.exists():
            shutil.copy2(source, Path(f'{target}.{suffix}'))
    for suffix in ('save', 'chkpt'):
        source = Path(f'{prefix}.{suffix}.{frame:06d}')
        dest = Path(f'{target}.{suffix}.{frame:06d}')
        if source.exists():
            shutil.copy2(source, dest)
        else:
            shutil.copy2(Path(str(source)+'.vtkhdf'), Path(str(dest)+'.vtkhdf'))
        if suffix != 'chkpt':
            continue
        if dest.exists():
            data = dest.read_bytes()
            lines = data[:4096].split(b'\0', 1)[0].splitlines()
            header_size = int(lines[1].split(b'\t')[1])
            payloads = fields(dest)
            for name in drop:
                payloads.pop(name, None)
            payloads.update(replace)
            header = lines[0]+b'\n'
            body = b''
            for name, value in payloads.items():
                header += name.encode()+b'\t'+str(header_size+len(body)).encode()+b'\n'
                body += value
            assert len(header) < header_size
            dest.write_bytes(header.ljust(header_size, b'\0')+body)
        else:
            import h5py
            import numpy as np
            with h5py.File(str(dest)+'.vtkhdf', 'r+') as stream:
                for name in drop:
                    if name in stream:
                        del stream[name]
                for name, value in replace.items():
                    dataset = stream[name]
                    dataset[...] = np.frombuffer(value, dtype=dataset.dtype).reshape(dataset.shape)


def reject(exe, folder, cfg, name, reason):
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    result = subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120,
                            env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    (folder/f'{name}.log').write_text(result.stdout)
    assert result.returncode == 22, (name, result.returncode, result.stdout[-2000:])
    assert reason in result.stdout, (name, result.stdout[-2000:])


def rsf_reference(times, theta0, widths):
    """Aging Euler solution and scalar boundary work minus elastic storage."""
    a, b, distance, speed0 = .01, .015, .01, .001
    speed = 2*heat.RATE*widths
    intervals = np.diff(times)
    equilibrium = distance/speed
    factors = np.vstack((np.ones(len(widths)), np.cumprod(1-intervals[:, None]*speed/distance, axis=0)))
    theta = equilibrium+(theta0-equilibrium)*factors
    mu = np.tan(np.pi/6)+a*np.log(speed/speed0)+b*np.log(speed0*theta/distance)
    strength = (heat.YIELD+heat.CONFINING*mu)/np.sqrt(1+mu*mu)
    tau, energy = np.zeros_like(theta), np.zeros_like(theta)
    for i, dt in enumerate(intervals, 1):
        tau[i] = np.minimum(tau[i-1]+2*heat.SHEAR*heat.RATE*dt, strength[i])
        external = (tau[i-1]+tau[i])*heat.RATE*dt
        stored = (tau[i]*tau[i]-tau[i-1]*tau[i-1])/(2*heat.SHEAR)
        energy[i] = energy[i-1]+np.maximum(0., external-stored)
    return theta, mu, tau, energy


def rsf_case(exe, folder, dim, dt=.005, steps=12, coupled=False, method='pt1'):
    folder.mkdir()
    cfg = heat.configuration(dim, method, dt, steps)
    cfg['mat'].update(rheology_type='elasto-plastic-rate-state-friction',
                      direct_a='[.01]', evolution_b='[.015]', state_var_model='1',
                      characteristic_distance='[.01]', characteristic_velocity='[.001]',
                      friction_angle0='[30]', friction_angle1='[30]')
    cfg['control']['rsf_slip_rate_projection_option'] = '1'
    capacity = heat.CAPACITY
    if coupled:
        cfg['control'].update(has_hydraulic_diffusion='true', has_thermal_diffusion='true',
                              has_fluid_thermal_advection='true', has_thermal_pressurization='true')
        cfg['mat'].update(porosity='[.2]', biot_coeff='[.8]', bulk_modulus_s='[1e11]',
                          fluid_bulk_modulus='[2.2e9]', hydraulic_perm='[0]',
                          fluid_rho0='[1000]', fluid_heat_capacity='[4180]', fluid_alpha='[2e-4]')
        capacity = .8*2700*1000+.2*1000*4180
    model = heat.run(exe, folder, cfg)
    times = np.array([row['time'] for row in heat.scan_frames(model.modelname)])
    theta0 = model.read_field(model.frames[0], 'friction state variable')
    x, conn = (model.read_field(model.frames[0], key) for key in ('coordinate', 'connectivity'))
    volume = abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    facets = []
    for face in itertools.combinations(range(dim+1), dim):
        points = x[conn][:, face, :]
        area = np.linalg.norm(points[:, 1]-points[:, 0], axis=1) if dim == 2 else \
            np.linalg.norm(np.cross(points[:, 1]-points[:, 0], points[:, 2]-points[:, 0]), axis=1)/2
        facets.append(area)
    widths = dim*volume/np.max(facets, axis=0)
    theta, mu, tau, energy = rsf_reference(times, theta0, widths)
    result = heat.analyze(model, capacity=capacity,
        reference=lambda ts: rsf_reference(ts, theta0, widths)[3] @ volume / volume.sum())
    for i, frame in enumerate(model.frames[1:], 1):
        np.testing.assert_allclose(model.read_field(frame, 'friction state variable'), theta[i], rtol=2e-13)
        np.testing.assert_allclose(model.read_field(frame, 'dynamic friction coefficient'), mu[i], rtol=2e-13)
        stress = model.read_field(frame, 'stress')
        np.testing.assert_allclose((stress[:, 0]-stress[:, dim-1])/2, tau[i], rtol=2e-12)
    restart = copy.deepcopy(cfg)
    middle = steps//2
    restart['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'cell'),
                          restarting_from_frame=str(middle))
    resumed = heat.run(exe, folder, restart, 'resumed')
    result['exact_continuation_payloads'] = heat.continuation(
        model, resumed, ('friction state variable', 'dynamic friction coefficient'))
    restart['control']['fixed_dt'] = repr(.4*dt)
    changed = heat.run(exe, folder, restart, 'changed')
    result['changed_timestep'] = heat.analyze(changed, initial_generated=energy[middle] @ volume, capacity=capacity,
        reference=lambda ts: rsf_reference(np.r_[times[:middle], ts], theta0, widths)[3][middle:] @ volume / volume.sum())
    result['element_widths_m'] = widths.tolist()
    result['coupled'] = coupled
    result['heat_capacity_J_m3_K'] = capacity
    if coupled:
        storage = .2/2.2e9+(.8-.2)/1e11
        beta = .2*2e-4
        lag = .8**2/(heat.BULK+(heat.SHEAR/3 if dim == 2 else 0)) if method.startswith('pt') else 0.
        epsilon0 = model.read_field(model.frames[0], 'strain')[:, :dim].sum(axis=1)
        balances = []
        for frame in model.frames:
            pressure = model.read_field(frame, 'pore pressure')[conn].mean(axis=1)
            temperature = model.read_field(frame, 'temperature')[conn].mean(axis=1)
            epsilon = model.read_field(frame, 'strain')[:, :dim].sum(axis=1)-epsilon0
            pending = heat.water.checkpoint_field(model, frame, 'pending pore pressure increment', len(x))
            balance = (storage*pressure+.8*epsilon-beta*(temperature-heat.T0)-lag*pending[conn].mean(axis=1)) @ volume
            balances.append(float(balance))
        assert max(abs(value) for value in balances) < 1e-15, balances
        result['split_fluid_budget_m3'] = balances
    result['scope'] = 'Nonzero a/b, evolving aging state and friction; scalar boundary work minus elastic storage.'
    (folder/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--study', action='store_true', help='Run only the prescribed RSF shear-cell study.')
    parser.add_argument('--coupled', action='store_true', help='Include sealed porosity and fluid thermal expansion.')
    parser.add_argument('--dt', type=float, default=.005)
    parser.add_argument('--steps', type=int, default=12)
    args = parser.parse_args()
    exe, folder = args.executable.resolve(), args.output.resolve()
    if args.study:
        result = rsf_case(exe, folder, args.ndims, args.dt, args.steps, args.coupled)
        shutil.copy2(__file__, folder/'checker.py')
        print(json.dumps({key: result[key] for key in ('dimensions', 'saved_frames', 'coupled',
                          'balance_error_over_yield_volume', 'work_error_over_yield_volume')}, indent=2))
        return
    folder.mkdir(parents=True, exist_ok=False)
    cfg = heat.configuration(args.ndims, 'pt1', .005, 6)
    base = heat.run(exe, folder, cfg)
    reference = heat.analyze(base)
    restart = copy.deepcopy(cfg)
    restart['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'cell'),
                          restarting_from_frame='3')
    resumed = heat.run(exe, folder, restart, 'resumed')
    exact = heat.continuation(base, resumed)
    changed = copy.deepcopy(restart)
    changed['control']['fixed_dt'] = '.002'
    changed_model = heat.run(exe, folder, changed, 'changed')
    changed_result = heat.analyze(changed_model, initial_generated=reference['rows'][3]['generated'])

    checkpoint_copy(folder/'cell', 0, folder/'legacy-zero', drop=HEAT_STATE)
    zero = copy.deepcopy(restart)
    zero['sim'].update(restarting_from_modelname=str(folder/'legacy-zero'), restarting_from_frame='0')
    legacy_zero = heat.run(exe, folder, zero, 'zero-resumed')
    legacy_exact = heat.continuation(base, legacy_zero)

    for key in HEAT_STATE:
        name = 'missing-'+key.replace(' ', '-')
        checkpoint_copy(folder/'cell', 3, folder/name, drop=(key,))
        broken = copy.deepcopy(restart)
        broken['sim']['restarting_from_modelname'] = str(folder/name)
        reject(exe, folder, broken, name+'-run', 'Incomplete shear-heat checkpoint')
    checkpoint_copy(folder/'cell', 3, folder/'legacy-power', drop=HEAT_STATE)
    legacy = copy.deepcopy(restart)
    legacy['sim']['restarting_from_modelname'] = str(folder/'legacy-power')
    reject(exe, folder, legacy, 'legacy-power-run', 'no recoverable pending energy')
    disabled = copy.deepcopy(restart)
    disabled['control']['has_shear_heating'] = 'false'
    reject(exe, folder, disabled, 'disable-pending', 'Cannot disable shear heating')
    checkpoint_copy(folder/'cell', 3, folder/'unknown-version',
                    replace={'shear heat version': struct.pack('i', 2)})
    unknown = copy.deepcopy(restart)
    unknown['sim']['restarting_from_modelname'] = str(folder/'unknown-version')
    reject(exe, folder, unknown, 'unknown-version-run', 'Unsupported shear-heat checkpoint version')
    result = dict(dimensions=args.ndims, exact_continuation_payloads=exact,
                  legacy_zero_exact_payloads=legacy_exact, rejected_schemas=6,
                  reference=reference, changed_timestep=changed_result,
                  rsf=rsf_case(exe, folder/'rsf', args.ndims),
                  scope='Plastic and aging EP-RSF work and checkpoint ownership; this is not a complete earthquake cycle.')
    (folder/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    shutil.copy2(__file__, folder/'checker.py')
    print(json.dumps({k: result[k] for k in ('dimensions', 'exact_continuation_payloads',
                                           'legacy_zero_exact_payloads', 'rejected_schemas')}, indent=2))


if __name__ == '__main__':
    main()
