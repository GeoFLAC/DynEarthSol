#!/usr/bin/env python3
"""Terzaghi loading/unloading in fixed, prescribed material regimes.

An initial compressive load is reduced by half at Tv=0.5; the
column then drains through Tv=2.5. Signed step superposition supplies both
pressure and settlement. Mesh and timestep sequences are separate.
"""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

import fluid_content as water
from Dynearthsol import Dynearthsol, scan_frames
from terzaghi_reference import TerzaghiReference

REGIMES = {
    'soft-fluid': dict(bulk=2e7, shear=1e7, porosity=.25, fluid_bulk=1e6,
                       grain_bulk=1e9, biot=.98, load=1e4, permeability=1e-11),
    'water': dict(bulk=2e7, shear=1e7, porosity=.25, fluid_bulk=2.2e9,
                  grain_bulk=1e9, biot=.98, load=1e4, permeability=1e-11),
    'soft-solid': dict(bulk=2e5, shear=1e5, porosity=.4, fluid_bulk=2.2e9,
                       grain_bulk=1e9, biot=.9998, load=100., permeability=1e-14),
}


def run(exe, folder, cfg, name):
    cfg = copy.deepcopy(cfg)
    cfg['sim']['modelname'] = name
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    with (folder/f'{name}.log').open('w') as stream:
        subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, check=True,
                       stdout=stream, stderr=subprocess.STDOUT, timeout=600,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    return Dynearthsol(str(folder/name))


def configuration(regime, spacing, count):
    props = dict(height=1., viscosity=.001, **REGIMES[regime])
    ref = TerzaghiReference(**props)
    clock = 1/ref.diffusivity
    cfg = water.configuration(spacing, clock/count, 2.5*clock, 'pt', 1e5, water.INFLOW)
    cfg['sim'].update(max_steps=str(count//2), max_time_in_yr='1e6', output_step_interval=str(count//20))
    cfg['control'].update(PT_relative_tolerance='1e-11', PT_absolute_tolerance='1e-9')
    cfg['bc'].update(hbc_z0='0', hbc_flux='[0,0,0,0,0,0]', stress_val_z1=repr(-props['load']))
    names = dict(bulk='bulk_modulus', shear='shear_modulus', porosity='porosity', fluid_bulk='fluid_bulk_modulus',
                 grain_bulk='bulk_modulus_s', biot='biot_coeff', permeability='hydraulic_perm')
    for name, parameter in names.items():
        cfg['mat'][parameter] = json.dumps([props[name]])
    return cfg, ref, clock


def analyze(folder, regime, spacing, count):
    _, ref, clock = configuration(regime, spacing, count)
    first, second = (Dynearthsol(str(folder/name)) for name in ('loaded', 'released'))
    x, conn = (first.read_field(first.frames[0], key) for key in ('coordinate', 'connectivity'))
    dim = int(first.ndims)
    volume = np.abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    area = volume.sum()/ref.height
    records = [(first, row, False) for row in scan_frames(str(folder/'loaded'))]
    records += [(second, row, True) for row in scan_frames(str(folder/'released'))[1:]]
    history, truncation = [], 0.
    for model, row, released in records:
        frame, time = row['frame'], row['time']
        np.testing.assert_array_equal(model.read_field(frame, 'coordinate'), x)
        p = model.read_field(frame, 'pore pressure')
        strain = model.read_field(frame, 'strain')
        events = [(0., ref.load)]+([(.5*clock, -.5*ref.load)] if released else [])
        expected, settlement = ref.load_history(-x[:, -1], time, events)
        shorter, _ = ref.load_history(-x[:, -1], time, events, terms=256)
        truncation = max(truncation, float(np.max(abs(shorter-expected)))/ref.load)
        measured = -float(volume @ strain[:, dim-1])/area
        error = float(np.max(abs(p-expected[0]))/ref.load)
        serror = abs(measured-settlement[0])/ref.drained_settlement
        history.append([time, time/clock, float(p.max()), float(p.min()), measured, settlement[0], error, serror])
    a = np.asarray(history)
    assert np.isfinite(a).all() and truncation < 1e-10
    np.testing.assert_allclose(a[-1, 1], 2.5, rtol=1e-12)
    np.savetxt(folder/'history.csv', a, delimiter=',', comments='',
               header='time_s,Tv,pmax_Pa,pmin_Pa,settlement_m,reference_settlement_m,pressure_error_over_initial_load,settlement_error_over_initial_drained')
    result = dict(regime=regime, properties=REGIMES[regime], dimension=dim, spacing_m=spacing,
                  steps_per_Tv=count, clock_s=clock, frames=len(records), nodes=len(x), elements=len(conn),
                  max_pressure_error=float(a[1:, 6].max()), max_settlement_error=float(a[1:, 7].max()),
                  final_pressure_error=float(a[-1, 6]), final_settlement_error=float(a[-1, 7]),
                  series_truncation_over_load=truncation)
    water.write_json(folder/'analysis.json', result)
    return result


def case(exe, folder, regime, spacing, count):
    folder.mkdir(parents=True, exist_ok=False)
    cfg, _, _ = configuration(regime, spacing, count)
    first = run(exe, folder, cfg, 'loaded')
    shifted = copy.deepcopy(cfg)
    shifted['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'loaded'),
                          restarting_from_frame=str(first.frames[-1]), max_steps=str(5*count//2))
    shifted['bc']['stress_val_z1'] = repr(-REGIMES[regime]['load']/2)
    run(exe, folder, shifted, 'released')
    return analyze(folder, regime, spacing, count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--regimes', nargs='+', choices=tuple(REGIMES), default=list(REGIMES))
    parser.add_argument('--spacings', nargs='+', type=float, default=[.2, .1, .05])
    parser.add_argument('--counts', nargs='+', type=int, default=[1000, 2000, 4000])
    parser.add_argument('--time-spacing', type=float, default=.1)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if min(args.spacings+[args.time_spacing]) <= 0 or any(n <= 0 or n % 20 for n in args.counts):
        parser.error('spacings must be positive and counts must be positive multiples of 20')
    exe, output = args.executable.resolve(), args.output.resolve()
    if args.replay:
        report = json.loads((output/'summary.json').read_text())
        for name, expected in report['artifacts_sha256'].items():
            assert water.digest(output/name) == expected, name
        for previous in report['cases']:
            name = f"{previous['regime']}/h{previous['spacing_m']:g}-n{previous['steps_per_Tv']}"
            assert analyze(output/name, previous['regime'], previous['spacing_m'], previous['steps_per_Tv']) == previous
        print('PASS: material/load references and artifact hashes replayed')
        return
    output.mkdir(parents=True, exist_ok=False)
    report = dict(executable_sha256=water.digest(exe), runner_sha256=water.digest(__file__),
                  reference_sha256=water.digest(Path(__file__).with_name('terzaghi_reference.py')), cases=[], passed=False)
    plan = sorted({(h, max(args.counts)) for h in args.spacings} |
                  {(args.time_spacing, n) for n in args.counts})
    for regime in args.regimes:
        for spacing, count in plan:
            result = case(exe, output/regime/f'h{spacing:g}-n{count}', regime, spacing, count)
            report['cases'].append(result)
            print(json.dumps(result), flush=True)
            water.write_json(output/'summary.json', report)
    assert water.digest(exe) == report['executable_sha256']
    report['artifacts_sha256'] = {str(p.relative_to(output)): water.digest(p) for p in sorted(output.rglob('*'))
        if p.is_file() and p.name not in ('analysis.json', 'history.csv', 'summary.json')}
    report['passed'] = True
    water.write_json(output/'summary.json', report)


if __name__ == '__main__':
    main()
