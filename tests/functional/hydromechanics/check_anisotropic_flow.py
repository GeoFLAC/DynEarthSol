#!/usr/bin/env python3
"""Check two diffusion modes across a prescribed anisotropic permeability jump.

All four x/z faces drain at zero pressure; 3D is a sealed extrusion.
Strain-to-pressure feedback is disabled on a fixed mesh; pressure can still
load the solid, so the timestep must also resolve elastic waves. The coefficient
changes through a restart, not a continuous state-dependent material law.
The analytical modes have different horizontal wave numbers, so exchanging
the directional diffusivities cannot masquerade as the correct decay.
"""
import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames
from pt_initial_checkpoint import fields

P0, MU = 1000., 1.002e-3
STORAGE = .3/2.17e9 + .7/37e9
PERMEABILITY = (1e-15, 2e-15)
TENSORS = ([4, 1, 1, 0, 0, 0], [1, 1, 9, 0, 0, 0])
INITIAL = np.array([1., .35])


from des_case import digest, read_config, run_solver, source_revision


def basis(x):
    return np.column_stack([np.sin(n*np.pi*x[:, 0])*np.sin(np.pi*(x[:, -1]+1))
                            for n in (1, 2)])


def amplitudes(time, switch):
    decay = np.zeros(2)
    for duration, permeability, tensor in zip((min(time, switch), max(0., time-switch)),
                                             PERMEABILITY, TENSORS):
        decay += duration*permeability/(MU*STORAGE)*np.pi**2*(tensor[0]*np.array([1., 4.])+tensor[2])
    return INITIAL*np.exp(-decay)


def configuration(h, dt, dim):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='seed', max_steps='1', output_step_interval='1',
                      max_time_in_yr='1', output_time_interval_in_yr='1',
                      has_initial_checkpoint='yes', checkpoint_frame_interval='1')
    cfg['mesh'].update(meshing_option='1', meshing_elem_shape='1', xlength='1',
                       ylength=repr(h if dim == 3 else 1.), zlength='1', resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h), remeshing_option='0')
    cfg['control'].update(has_PT='false', is_quasi_static='false', fixed_dt=repr(dt),
                          has_poroelastic_pressure_feedback='false', has_thermal_diffusion='false',
                          gravity='0', has_moving_mesh='false', damping_option='0',
                          inertial_scaling='1e5', characteristic_speed='1e-5')
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        cfg['bc'].update({f'vbc_{side}': '3', f'vbc_val_{side}': '0',
                          f'hbc_{side}': '0' if side[0] == 'y' else '2'})
    cfg['bc'].update(hbc_pressure='[0,0,0,0,0,0]', stress_bc_z1='0')
    cfg['ic'].update(pore_pressure_option='2', pore_pressure_value='0',
                     has_body_force_adjustment='false', weakzone_option='0')
    cfg['mat'].update(porosity='[0.3]', biot_coeff='[1]', fluid_bulk_modulus='[2.17e9]',
                      bulk_modulus_s='[37e9]', fluid_visc=repr([MU]),
                      hydraulic_perm=repr([PERMEABILITY[0]]), permeability_tensor=repr(TENSORS[0]))
    return cfg


def run(executable, folder, cfg, name):
    cfg['sim']['modelname'] = name
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    run_solver([str(executable), f'{name}.cfg'], log=folder / f'{name}.log', cwd=folder, check=True, timeout=180, env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    return Dynearthsol(str(folder/name))


def case(executable, folder, h, steps, dim):
    folder.mkdir()
    cfg = configuration(h, 1/steps, dim)
    seed = run(executable, folder, cfg, 'seed')
    assert seed.format == 'binary', 'Use a binary executable for the analytical initial-state seed'
    assert seed.ndims == dim
    x, conn = seed.read_field(0, 'coordinate'), seed.read_field(0, 'connectivity')
    modes = basis(x)
    initial = P0*(modes@INITIAL)
    for d in (0, dim-1):
        initial[np.isclose(x[:, d], x[:, d].min()) | np.isclose(x[:, d], x[:, d].max())] = 0.
    seed.read_header(0)
    with Path(seed.get_fn(0)).open('r+b') as stream:
        stream.seek(seed.field_pos['pore pressure']); stream.write(initial.astype(np.float64).tobytes())
    cfg['sim'].update(is_restarting='true', restarting_from_modelname='seed', restarting_from_frame='0',
                      output_step_interval=str(steps//4), max_steps=str(steps//2))
    first = run(executable, folder, cfg, 'first')
    switch = scan_frames(str(folder/'first'))[-1]['time']
    cfg['mat'].update(hydraulic_perm=repr([PERMEABILITY[1]]), permeability_tensor=repr(TENSORS[1]))
    cfg['sim'].update(restarting_from_modelname='first', restarting_from_frame='2', max_steps=str(steps))
    second = run(executable, folder, cfg, 'second')
    cfg['sim'].update(restarting_from_modelname='second', restarting_from_frame='3')
    run(executable, folder, cfg, 'resumed')
    return analyze(folder, h, steps, dim)


def analyze(folder, h, steps, dim):
    first, second = Dynearthsol(str(folder/'first')), Dynearthsol(str(folder/'second'))
    x, conn = first.read_field(0, 'coordinate'), first.read_field(0, 'connectivity')
    modes = basis(x)
    switch = scan_frames(str(folder/'first'))[-1]['time']
    for kind, frame in (('save', 3), ('save', 4), ('chkpt', 4)):
        a, b = fields(folder/f'second.{kind}.{frame:06d}'), fields(folder/f'resumed.{kind}.{frame:06d}')
        assert a.keys() == b.keys()
        assert all(a[k] == b[k] for k in a if k not in ('walltime_sec', 'provenance'))
    volumes = np.abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    weights = np.zeros(len(x))
    np.add.at(weights, conn, np.broadcast_to((volumes/(dim+1))[:, None], conn.shape))
    rows = []
    for name, model in (('first', first), ('second', second)):
        for row in scan_frames(str(folder/name)):
            if name == 'second' and row['frame'] == 2: continue
            p = model.read_field(row['frame'], 'pore pressure')
            assert np.isfinite(p).all()
            assert np.isfinite(model.read_field(row['frame'], 'velocity')).all()
            assert np.isfinite(model.read_field(row['frame'], 'strain')).all()
            exact = amplitudes(row['time'], switch)
            fitted = np.linalg.lstsq(modes, p/P0, rcond=None)[0]
            content = STORAGE*float(weights@p)
            reference_content = STORAGE*P0*4/np.pi**2*exact[0]*(h if dim == 3 else 1.)
            rows.append(dict(time=row['time'], phase=name, frame=row['frame'],
                             pressure_error_over_initial_scale=float(np.max(np.abs(p/P0-modes@exact))),
                             modal_amplitudes=fitted.tolist(), reference_amplitudes=exact.tolist(),
                             modal_error=float(np.max(np.abs(fitted-exact))),
                             content_error_over_initial=float(abs(content-reference_content)/
                                 (STORAGE*P0*4/np.pi**2*(h if dim == 3 else 1.)))))
    outputs = {p.name: digest(p) for p in folder.iterdir() if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.cfg'}
    return dict(spacing_m=h, steps=steps, nodes=len(x), elements=len(conn), switch_time=switch,
                frames=rows, maximum_pressure_error=max(r['pressure_error_over_initial_scale'] for r in rows),
                maximum_modal_error=max(r['modal_error'] for r in rows), exact_restart=True, files_sha256=outputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=2)
    parser.add_argument('--spacings', type=float, nargs='+', default=[.1, .05, .025])
    parser.add_argument('--steps', type=int, default=20000)
    parser.add_argument('--binary-manifest', type=Path)
    args = parser.parse_args()
    if args.steps < 4 or args.steps % 4 or not all(np.isfinite(h) and 0 < h <= .25 for h in args.spacings):
        parser.error('Use positive spacings <=0.25 and steps divisible by four')
    executable, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = dict(executable_sha256=digest(executable), runner_sha256=digest(__file__), dimension=args.dimension,
                  source_commit=source_revision(ROOT), cases=[])
    if args.binary_manifest:
        q = json.loads(args.binary_manifest.read_text()); assert report['executable_sha256'] in q['binaries'].values()
        report.update(solver_commit=q['implementation_commit'], binary_manifest_sha256=digest(args.binary_manifest))
    (output/'runner.py').write_bytes(Path(__file__).read_bytes())
    for h in args.spacings:
        record = case(executable, output/f'h{h:g}', h, args.steps, args.dimension)
        report['cases'].append(record)
        (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
        print(h, record['maximum_pressure_error'], record['maximum_modal_error'], flush=True)
    assert digest(executable) == report['executable_sha256']


if __name__ == '__main__':
    main()
