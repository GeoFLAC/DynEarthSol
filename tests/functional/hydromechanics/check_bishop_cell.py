#!/usr/bin/env python3
"""DES sealed Bishop cell: the partly saturated undrained response under compression.

A homogeneous, sealed, almost impermeable block at uniform suction is
compressed uniaxially at a prescribed rate. Without flow the implemented
Bishop rate equation, (S0 + phi (1 - Sr) dSe/dp) dp/dt = -alpha chi deps/dt
with chi = Se, gives the exact pressure along the strain path,

    dp/d eps = -alpha Se(p) / (S0 + phi (1 - Sr) dSe/dp),

with no skeleton-stiffness term: the fixed-stress stabilisation must cancel
as the step shrinks. This tests the rate equation, not conservation of a
product state alpha*chi*eps, whose derivative would add alpha*eps*dchi/dp.
Two states: moderate suction (chi*Pi' > 0) and strong suction, where
Pi' = Se + p dSe/dp < 0 for the Gardner law (alpha_p*suction > 1). The pressure of every PT run, at several step counts,
and of the explicit strain-form run is compared with this ODE integrated
along each run's own vertical strain.

    check_bishop_cell.py EXECUTABLE OUTPUT [--ndims 2] [--steps 10 20 40 80] [--explicit]
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics')]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from unsaturated_reference import Retention

MAT = dict(bulk=1e4, shear=1e4, biot=1.0, porosity=.3, fluid_bulk=2.2e9, grain_bulk=1e20, permeability=1e-30,
           viscosity=1e-3, sr=.1)
STATES = dict(moderate=(-5000.0, 1e-4), strong=(-2000.0, 1e-3))   # initial pressure (Pa), Gardner alpha (1/Pa)
STRAIN, DURATION = -1e-3, 1000.0                   # final vertical strain, end time (s)
S0 = MAT['porosity']/MAT['fluid_bulk'] + (MAT['biot'] - MAT['porosity'])/MAT['grain_bulk']


from des_case import digest, read_config, reset_boundaries, run_solver, source_diff, source_revision


def exact(state, strain):
    """Pressure along a strain path from the initial state (strain increasing in magnitude)."""
    p0, alpha = STATES[state]
    law = Retention('gardner', alpha, sr=MAT['sr'])
    def rate(eps, p):
        se = float(law.se(p[0])); dse = alpha*se if p[0] < 0 else 0.0
        return [-MAT['biot']*se/(S0 + MAT['porosity']*(1 - MAT['sr'])*dse)]
    sol = solve_ivp(rate, (0.0, float(np.min(strain))), [p0], method='DOP853', rtol=1e-12, atol=1e-9,
                    dense_output=True)
    return sol.sol(np.asarray(strain, float))[0]


def config(ndims, steps, state, pt_option=0):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    explicit = steps == 0
    cfg['sim'].update(modelname='cell', max_time_in_yr=repr(DURATION/YEAR2SEC),
                      max_steps=str(10**8 if explicit else steps), output_step_interval=str(10**8 if explicit else max(1, steps//10)),
                      output_time_interval_in_yr=repr(DURATION/10/YEAR2SEC) if explicit else '1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'].update(meshing_elem_shape='1', xlength='1', ylength='1', zlength='1', resolution='0.5',
                       largest_size='0.5', smallest_size='0.5')
    control = dict(gravity='0', has_thermal_diffusion='false', has_poroelastic_pressure_feedback='true',
                   PT_max_iter='100000')
    if explicit:
        control.update(has_PT='false', inertial_scaling='1e5', fixed_dt='0', use_global_velocity_scaling='false')
    else:
        control.update(fixed_dt=repr(DURATION/steps), PT_option=str(pt_option))
    cfg['control'].update(control)
    bc = cfg['bc']
    reset_boundaries(bc, vbc='1', vbc_val='0', hbc='0')
    bc['vbc_val_z1'] = repr(STRAIN/DURATION)              # the cell is 1 m high
    bc.update(stress_bc_z1='0', hbc_pressure='[0,0,0,0,0,0]')
    p0, alpha = STATES[state]
    cfg['ic'].update(pore_pressure_option='2', pore_pressure_value=repr(p0), has_body_force_adjustment='false')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']!r}]", shear_modulus=f"[{MAT['shear']!r}]",
                      biot_coeff=f"[{MAT['biot']!r}]", porosity=f"[{MAT['porosity']!r}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']!r}]", bulk_modulus_s=f"[{MAT['grain_bulk']!r}]",
                      hydraulic_perm=f"[{MAT['permeability']!r}]", fluid_visc=f"[{MAT['viscosity']!r}]",
                      retention_model='[3]', retention_alpha=f'[{alpha!r}]',
                      residual_saturation=f"[{MAT['sr']!r}]", retention_n='[2]')
    return cfg


def compare(folder, state):
    model = Dynearthsol(str(folder/'cell')); rows = scan_frames(str(folder/'cell'))
    zz = model.read_field(0, 'coordinate').shape[1] - 1
    out = []
    for row in rows[1:]:
        eps = float(np.mean(model.read_field(row['frame'], 'strain')[:, zz]))
        p = model.read_field(row['frame'], 'pore pressure')
        ref = float(exact(state, [eps])[0])
        out.append(dict(time=row['time'], strain=eps, pressure_mean=float(p.mean()), pressure_spread=float(np.ptp(p)),
                        pressure_ref=ref))
    change = abs(out[-1]['pressure_ref'] - STATES[state][0])
    worst = max(abs(f['pressure_mean'] - f['pressure_ref']) for f in out)/change
    return dict(frames=out, max_error_over_change=worst, final_strain=out[-1]['strain'],
                final_pressure=out[-1]['pressure_mean'], final_ref=out[-1]['pressure_ref'])


def study(exe, out, ndims, steps, explicit, states, pt_option=0):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=source_revision(ROOT),
                  ndims=ndims, material=MAT, states=STATES, pt_option=pt_option, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for state, n in [(s, n) for s in states for n in list(steps) + ([0] if explicit else [])]:
        folder = out/f"{state}-{f'pt-{n}' if n else 'explicit'}"; folder.mkdir()
        with (folder/'cell.cfg').open('w') as f: config(ndims, n, state, pt_option).write(f)
        run_solver([str(exe), 'cell.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=20000)
        item = dict(state=state, mode='explicit' if n == 0 else 'PT', steps=n, **compare(folder, state))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    return result


def restart_checks(exe, out, ndims, state, steps=40):
    """Restart halfway: with the same step it must reproduce the fresh run exactly
    (the frozen increments come from the checkpointed pending increment); with the
    step divided by four it must still follow the exact rate equation. DES keeps
    the checkpointed step until its next slow update (10 steps), so the second
    run is a genuine variable-step sequence; it is compared at its last saved
    frame, on its own strain path."""
    sys.path.insert(0, str(ROOT/'tests'/'functional'))
    from pt_initial_checkpoint import fields
    env = dict(os.environ, OMP_NUM_THREADS='1')
    half = steps//2
    out.mkdir(parents=True, exist_ok=False)
    fresh = config(ndims, steps, state)
    fresh['sim'].update(output_step_interval=str(half), checkpoint_frame_interval='1')
    run = lambda folder, cfg, name: (folder.mkdir(exist_ok=True), cfg.write((folder/f'{name}.cfg').open('w')),
        subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, env=env, stdout=(folder/f'{name}.log').open('w'),
                       stderr=subprocess.STDOUT, check=True, timeout=20000))
    fresh['sim']['modelname'] = 'fresh'
    run(out, fresh, 'fresh')
    report = {}
    for name, factor in (('same', 1), ('finer', 4)):
        cfg = config(ndims, steps, state)
        cfg['sim'].update(modelname=name, is_restarting='yes', restarting_from_modelname=str(out/'fresh'),
                          restarting_from_frame='1', max_steps=str(half + half*factor),
                          output_step_interval=str(half),   # a different interval stops restart output
                          checkpoint_frame_interval='1')
        cfg['control']['fixed_dt'] = repr(DURATION/steps/factor)
        run(out, cfg, name)
    last = lambda name: out/f"{name}.save.{scan_frames(str(out/name))[-1]['frame']:06d}"
    a, b = fields(last('fresh')), fields(last('same'))
    assert scan_frames(str(out/'fresh'))[-1]['time'] == scan_frames(str(out/'same'))[-1]['time']
    report['same_step_identical_payloads'] = sum(a[k] == b[k] for k in a if k not in ('walltime_sec', 'provenance') and k in b)
    report['same_step_different_payloads'] = sum(a[k] != b[k] for k in a if k not in ('walltime_sec', 'provenance') and k in b)
    zz = ndims - 1                                      # strain is (xx, zz, xz) in 2D, (xx, yy, zz, ...) in 3D
    for name in ('fresh', 'finer'):
        # A restarted run has no frame 0, so read its last frame's payloads directly.
        row = scan_frames(str(out/name))[-1]
        f = fields(out/f"{name}.save.{row['frame']:06d}")
        strain = np.frombuffer(f['strain'], dtype=np.float64).reshape(-1, 3 if ndims == 2 else 6)
        eps = float(np.mean(strain[:, zz]))
        p = float(np.mean(np.frombuffer(f['pore pressure'], dtype=np.float64)))
        ref = float(exact(state, [eps])[0])
        report[f'{name}_time'] = row['time']
        report[f'{name}_error_over_change'] = abs(p - ref)/abs(ref - STATES[state][0])
    (out/'restart.json').write_text(json.dumps(report, indent=1))
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, default=2)
    p.add_argument('--steps', type=int, nargs='*', default=[10, 20, 40, 80])
    p.add_argument('--explicit', action='store_true')
    p.add_argument('--states', nargs='+', default=list(STATES))
    p.add_argument('--restart', action='store_true', help='run only the restart and changed-step checks')
    p.add_argument('--pt-option', type=int, default=0, help='control.PT_option of the PT runs')
    a = p.parse_args()
    if a.restart:
        for state in a.states:
            print(state, json.dumps(restart_checks(a.executable.resolve(), a.output.resolve()/state, a.ndims, state)))
        sys.exit(0)
    study(a.executable.resolve(), a.output.resolve(), a.ndims, a.steps, a.explicit, a.states, a.pt_option)
