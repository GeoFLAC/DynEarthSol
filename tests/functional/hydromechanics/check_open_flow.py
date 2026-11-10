#!/usr/bin/env python3
"""Drain a smooth pressure pulse through two open faces and account for fluid.

The independently specified reference is
  p(x,t) = P sin(pi x/L) exp(-D pi^2 t/L^2).
The side faces are sealed; both ends hold zero pressure. Exact stored fluid
is 2 S P A L/pi times the same exponential. Check pressure and storage on
refined 3D meshes, and integrate every saved physical-step boundary volume
against fluid content computed independently from the mesh and pressure.
Constant properties and disabled mechanical feedback isolate transport.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

import check_theis as ct

L, WIDTH, P, END, DT = 10., 1., 1000., 10., .05
RATE = ct.DIFFUSIVITY*math.pi**2/L**2


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(executable, folder, cfg, name):
    with (folder/f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    with (folder/f'{name}.log').open('w') as stream:
        subprocess.run([str(executable), f'{name}.cfg'], cwd=folder, check=True,
                       stdout=stream, stderr=subprocess.STDOUT, timeout=600,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))


def run(executable, folder, spacing):
    folder.mkdir()
    cfg = ct.config(spacing)
    cfg['sim'].update(modelname='seed', max_steps='1', output_step_interval='1',
                      checkpoint_frame_interval='1', has_initial_checkpoint='true')
    cfg['mesh'].update(xlength=str(L), ylength=str(WIDTH), zlength=str(WIDTH))
    cfg['control']['fixed_dt'] = str(DT)
    for key in list(cfg['bc']):
        if key.startswith(('well_', 'num_well')):
            del cfg['bc'][key]
    cfg['bc'].update(num_wells='0', hbc_x0='2', hbc_x1='2',
                     hbc_pressure='[0,0,0,0,0,0]')
    execute(executable, folder, cfg, 'seed')
    seed = ct.Dynearthsol(str(folder/'seed'))
    assert seed.ndims == 3 and seed.format == 'binary'
    x = seed.read_field(0, 'coordinate')
    conn = seed.read_field(0, 'connectivity')
    initial = P*np.sin(math.pi*x[:, 0]/L)
    initial[np.isclose(x[:, 0], 0) | np.isclose(x[:, 0], L)] = 0.
    seed.read_header(0)
    with Path(seed.get_fn(0)).open('r+b') as stream:
        stream.seek(seed.field_pos['pore pressure'])
        stream.write(initial.astype(np.float64).tobytes())
    np.testing.assert_array_equal(seed.read_field(0, 'pore pressure'), initial)
    cfg['sim'].update(modelname='pulse', is_restarting='true',
                      restarting_from_modelname='seed', restarting_from_frame='0',
                      max_steps=str(round(END/DT)))
    execute(executable, folder, cfg, 'pulse')
    model = ct.Dynearthsol(str(folder/'pulse'))
    rows = ct.scan_frames(str(folder/'pulse'))
    assert len(rows) == round(END/DT)+1
    cells = x[conn]
    volume = np.abs(np.linalg.det(cells[:, 1:]-cells[:, :1]))/6
    initial_stored = ct.STORAGE*float(volume @ initial[conn].mean(axis=1))
    exact_initial_stored = 2*ct.STORAGE*P*WIDTH**2*L/math.pi
    cumulative_boundary = 0.
    history = []
    for row in rows[1:]:
        frame, time = row['frame'], row['time']
        p = model.read_field(frame, 'pore pressure')
        assert np.isfinite(p).all()
        np.testing.assert_array_equal(model.read_field(frame, 'coordinate'), x)
        stored = ct.STORAGE*float(volume @ p[conn].mean(axis=1))
        balance = model.read_field(frame, 'fluid balance')
        np.testing.assert_array_equal(balance[1:3], 0.)
        cumulative_boundary += balance[3]
        exact_pressure = initial*math.exp(-RATE*time)
        exact_stored = exact_initial_stored*math.exp(-RATE*time)
        history.append(dict(time_s=time, stored_fluid_m3=stored,
                            exact_stored_fluid_m3=exact_stored,
                            cumulative_boundary_m3=float(cumulative_boundary),
                            global_budget_error_m3=abs(stored-initial_stored-cumulative_boundary),
                            pressure_error_over_amplitude=float(np.max(np.abs(p-exact_pressure)))/P,
                            storage_error_over_initial=abs(stored-exact_stored)/exact_initial_stored))
    worst_budget = max(r['global_budget_error_m3'] for r in history)
    assert worst_budget/exact_initial_stored < 1e-10
    assert cumulative_boundary < 0 and history[-1]['stored_fluid_m3'] < initial_stored
    result = dict(spacing_m=spacing, nodes=len(x), elements=len(conn),
                  max_pressure_error_over_amplitude=max(r['pressure_error_over_amplitude'] for r in history),
                  max_storage_error_over_initial=max(r['storage_error_over_initial'] for r in history),
                  max_global_budget_error_m3=worst_budget,
                  global_budget_error_over_initial=worst_budget/exact_initial_stored,
                  seed_save_sha256=digest(folder/'seed.save.000000'),
                  seed_checkpoint_sha256=digest(folder/'seed.chkpt.000000'),
                  input_sha256=digest(folder/'pulse.cfg'), history=history,
                  outputs_sha256={p.name: digest(p) for p in sorted(folder.glob('pulse.*'))
                                  if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.info'})
    print(folder.name, 'pressure error', result['max_pressure_error_over_amplitude'],
          'relative budget error', result['global_budget_error_over_initial'], flush=True)
    return result


def study(executable, output, spacings, manifest):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = digest(executable)
    solver = json.loads(manifest.read_text()) if manifest else None
    if solver:
        assert binary_hash in solver['binaries'].values()
        shutil.copy2(manifest, output/'solver-qualification.json')
    for source in (Path(__file__), Path(ct.__file__)):
        shutil.copy2(source, output/source.name)
    summary = dict(executable=str(executable), executable_sha256=binary_hash,
                   solver_commit=solver['implementation_commit'] if solver else None,
                   runner_sha256=digest(Path(__file__)),
                   reference=dict(length_m=L, area_m2=WIDTH**2, amplitude_Pa=P,
                                  diffusivity_m2_s=ct.DIFFUSIVITY, decay_rate_s_inv=RATE,
                                  storage_per_Pa=ct.STORAGE, duration_s=END, dt_s=DT), cases=[])
    for h in spacings:
        summary['cases'].append(run(executable, output/f'h{h:g}', h))
        (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    ordered = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    summary['passed'] = (ordered[-1]['max_pressure_error_over_amplitude'] < .01 and
                         ordered[-1]['max_pressure_error_over_amplitude'] < ordered[0]['max_pressure_error_over_amplitude'])
    assert digest(executable) == binary_hash
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    assert summary['passed'], 'Pressure accuracy/convergence target not met; keep the raw evidence'
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--spacings', type=float, nargs='+', default=[1., .5, .25])
    parser.add_argument('--binary-manifest', type=Path)
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.spacings,
          args.binary_manifest.resolve() if args.binary_manifest else None)
