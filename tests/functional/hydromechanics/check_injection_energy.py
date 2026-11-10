#!/usr/bin/env python3
"""DES heat ledger of a hot-water doublet in a sealed, insulated box.

Fluid at T_in is injected at (L/4, -L/2) and the same rate extracted at
(3L/4, -L/2) in a 2D box whose sides are sealed and insulated. The hot
front reaches the producer, so the extracted water warms. No heat or fluid
crosses the boundary, so by energy conservation alone (independent of the
form DES integrates)

    E(t) - E(0) = rho_f c_f q INT_0^t (T_in - T_prod(s)) ds,

with E = SUM C A_n T_n the stored sensible heat and T_prod the temperature
at the producer, interpolated with the element's shape functions. The
integral uses every saved step (trapezoid rule). The ledger error is
reported against the injected heat rho_f c_f q (T_in - T0) t. The fluid is
slightly compressible, so storage changes while the pressure settles.

Thermal advection uses the Darcy flux of the pressure at the start of each
step, so in the first step of a source switched on from rest there is none:
the injection node saturates and part of that step's heat is lost. That
start deficit is reported in steps of injected heat; the ledger error is
measured after it, and the total error includes it.

    check_injection_energy.py EXECUTABLE OUTPUT [--spacings 0.2 0.1] [--steps 400]
"""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(Path(__file__).resolve().parent)]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
import check_injection_temperature as hot

L, Q, DURATION = 8.0, 1e-4, 6e5
WELLS = ((L/4, -L/2), (3*L/4, -L/2))


from des_case import reset_boundaries, run_solver, source_diff, source_revision

def config(h, steps):
    cfg = hot.config(h, steps)
    cfg['mesh'].update(xlength=repr(L), zlength=repr(L))
    cfg['sim'].update(output_step_interval='1')
    cfg['control'].update(fixed_dt=repr(DURATION/steps))
    cfg['bc'].update(hbc_x1='0', hbc_z0='0')
    (xi, zi), (xp, zp) = WELLS
    cfg['injection'].update(num_points='2', points_x=f'[{xi!r},{xp!r}]', points_z=f'[{zi!r},{zp!r}]',
                            rate=f'[{Q!r},{-Q!r}]', start_time_in_yr='[0,0]',
                            temperature=f'[{hot.TIN!r},{hot.TIN!r}]')
    return cfg


def shape_weights(coord, conn, point):
    """Element containing the point and its linear shape-function weights."""
    for e, nodes in enumerate(conn):
        a, b, c = coord[nodes]
        m = np.array([[b[0]-a[0], c[0]-a[0]], [b[1]-a[1], c[1]-a[1]]])
        s, t = np.linalg.solve(m, np.asarray(point) - a)
        w = np.array([1 - s - t, s, t])
        if w.min() > -1e-9:
            return nodes, w
    raise ValueError(f'{point} is outside the mesh')


def ledger(folder):
    model = Dynearthsol(str(folder/'hot')); rows = scan_frames(str(folder/'hot'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    area = hot.lumped_area(coord, conn)
    nodes, w = shape_weights(coord, conn, WELLS[1])
    times, stored, produced = [], [], []
    for row in rows:
        T = model.read_field(row['frame'], 'temperature')
        times.append(row['time']); stored.append(hot.CAPACITY*float(area @ T)); produced.append(float(w @ T[nodes]))
    times, stored, produced = map(np.array, (times, stored, produced))
    net = np.concatenate([[0.0], np.cumsum(np.diff(times)*.5*(produced[1:] + produced[:-1]))])
    expected = hot.FLUID*Q*(hot.TIN*times - net)
    gained = stored - stored[0]
    deviation = gained - expected
    # Thermal advection uses the Darcy flux of the pressure at the start of the
    # step, which is zero in the first step of a source switched on from rest;
    # the injection node saturates and part of that step's heat is lost. Report
    # it in steps of injected heat, and the ledger after it separately.
    one_step = hot.FLUID*Q*(hot.TIN - hot.T0)*times[1]
    drift = np.abs(deviation[2:] - deviation[1])/(hot.FLUID*Q*(hot.TIN - hot.T0)*times[2:])
    return dict(nodes=len(coord), frames=len(rows), breakthrough=float(produced[-1] - hot.T0),
                producer_rise=float((produced[-1] - hot.T0)/(hot.TIN - hot.T0)),
                start_deficit_steps=float(-deviation[1]/one_step),
                ledger_error=float(drift.max()), final_ledger_error=float(drift[-1]),
                total_error=float(abs(deviation[-1])/(hot.FLUID*Q*(hot.TIN - hot.T0)*times[-1])),
                stored_gain=float(gained[-1]), injected=float(hot.FLUID*Q*hot.TIN*times[-1]),
                extracted=float(hot.FLUID*Q*net[-1]))


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=hot.digest(exe), runner_sha256=hot.digest(__file__),
                  source_commit=source_revision(ROOT),
                  steps=steps, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'hot.cfg').open('w') as f: config(h, steps).write(f)
        run_solver([str(exe), 'hot.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=40000)
        item = dict(spacing_m=h, input_sha256=hot.digest(folder/'hot.cfg'), **ledger(folder))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps(item), flush=True)
    return result


def pulse_study(executable, output, ndims):
    """Bounded mixing with zero Darcy flux; not a retained-fluid energy test.

    The source sits on a vertex, so only that nodal capacity is heated.
    Its independent ODE gives T = Tin - (Tin-T0)*exp(-rho*cf*amount/Cnode).
    A large coefficient exposes a time-averaged stability bound; its stored
    volume is deliberately outside the small-storage physical qualification.
    """
    output.mkdir(parents=True, exist_ok=False)
    cfg = hot.config(1, 4, ndims)
    cfg['mesh'].update(xlength='1', ylength='1', zlength='1', resolution='1')
    cfg['sim'].update(has_initial_checkpoint='true', checkpoint_frame_interval='1')
    cfg['control']['has_thermal_diffusion'] = 'false'
    cfg['mat'].update(hydraulic_perm='[0]', fluid_bulk_modulus='[1e3]')
    reset_boundaries(cfg['bc'], vbc='3', hbc='0')
    if ndims == 3:
        cfg['injection']['points_y'] = '[0]'
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')

    def run(name, steps):
        cfg['sim'].update(modelname=name, max_steps=str(steps), output_step_interval=str(max(1, steps//10)))
        cfg['control']['fixed_dt'] = repr(1/steps)
        with (output/f'{name}.cfg').open('w') as stream:
            cfg.write(stream)
        run_solver([str(executable), f'{name}.cfg'], log=output / f'{name}.log', cwd=output, env=env, check=True, timeout=120)
        return Dynearthsol(str(output/name))

    cfg['injection']['enabled'] = 'false'
    seed = run('seed', 1)
    assert seed.ndims == ndims
    x, conn = seed.read_field(0, 'coordinate'), seed.read_field(0, 'connectivity')
    nodal_volume = hot.lumped_area(x, conn)
    node = int(np.argmin(np.linalg.norm(x, axis=1)))
    assert np.linalg.norm(x[node]) < 1e-12
    amount = float(2.25*hot.CAPACITY*nodal_volume[node]/hot.FLUID)
    exact = hot.TIN - (hot.TIN-hot.T0)*np.exp(-2.25)
    cfg['injection']['enabled'] = 'true'
    cases = {}
    for name, steps, start, duration, rate_model in (
            ('constant', 1, 0., 1., 'constant_rate'),
            ('early', 1, 0., .01, 'constant_rate'),
            ('late', 1, .99, .01, 'constant_rate'),
            ('amount', 1, 0., .01, 'total_amount'),
            ('refined1000', 1000, 0., .01, 'constant_rate'),
            ('refined10000', 10000, 0., .01, 'constant_rate')):
        cfg['injection'].update(rate_model=rate_model, rate=json.dumps([amount/duration]),
                                total_amount=json.dumps([amount]), start_time_in_yr=json.dumps([start/YEAR2SEC]),
                                duration_in_yr=json.dumps([duration/YEAR2SEC]))
        model = run(name, steps)
        frames = scan_frames(str(output/name))
        temperatures = [model.read_field(row['frame'], 'temperature') for row in frames]
        low, high = min(float(t.min()) for t in temperatures), max(float(t.max()) for t in temperatures)
        final = float(temperatures[-1][node])
        cases[name] = dict(minimum=low, maximum=high, final=final,
                           bounded=low >= hot.T0-1e-10 and high <= hot.TIN+1e-10,
                           final_error_over_jump=abs(final-exact)/(hot.TIN-hot.T0))
    report = dict(executable_sha256=hot.digest(executable), ndims=ndims, cases=cases,
                  ode_final_temperature=float(exact), injected_amount=amount,
                  scope='Boundedness and time convergence of the implemented fixed-capacity mixing ODE; no changing-storage energy claim.')
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    assert all(case['bounded'] for case in cases.values()), report
    assert cases['refined10000']['final_error_over_jump'] < .005
    assert cases['refined10000']['final_error_over_jump'] < .2*cases['refined1000']['final_error_over_jump']
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1])
    p.add_argument('--steps', type=int, default=400)
    p.add_argument('--pulse', action='store_true', help='Test bounded inlet pulses and mixing-ODE time convergence')
    p.add_argument('--ndims', type=int, choices=(2, 3), default=2, help='Dimension for --pulse')
    a = p.parse_args()
    if a.pulse:
        pulse_study(a.executable.resolve(), a.output.resolve(), a.ndims)
    else:
        if a.ndims != 2:
            p.error('The doublet ledger is 2D; --ndims 3 requires --pulse')
        study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
