#!/usr/bin/env python3
"""Rate-and-state friction under a velocity step, at a fixed pore pressure.

The material-strength cell (check_pore_strength.py) imposes normal
extensions +e and -e on every node, so stress follows the constitutive law
alone: principal shear t = 2 G e until yield, then the Mohr-Coulomb load
t_y = c cos(phi) + (C - q) sin(phi) at confinement C and Bishop pressure q.
With EP-RSF, tan(phi) = mu follows

    mu = mu_0 + a ln(V/V_0) + b ln(V_0 theta/D_c),   d(theta)/dt = 1 - V theta/D_c   (aging law),

with slip rate V = 2 w eps_II (rsf_slip_rate_projection_option = 1; w is the
element's minimum altitude). The cell loads at rate e1 with V_0 = V_1, so it
starts and stays at steady state, then restarts from its last checkpoint at
10 e1. Exact responses:

- discrete: DES's forward-Euler state update, theta_{n+1} = theta_n +
  (1 - V theta_n/D_c) dt, then mu and t_n = min(t_{n-1} + 2 G de, t_y(mu_n));
  DES must match this to round-off (one state update per physical step);
- continuum: theta = D_c/V_2 + (theta_1 - D_c/V_2) exp(-V_2 (t - t_1)/D_c),
  approached at first order in the step.

    check_rsf_step.py EXECUTABLE OUTPUT [--refinements 1 2 4]
"""
import argparse
import configparser
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_pore_strength as ps

cell = ps.cell
G, COHESION, C, P, BIOT = ps.G, ps.COHESION, 1e6, 8e5, 1.0
A, B, DC = .01, .015, 1e-3
MU0 = math.tan(ps.PHI)
RATE1, RATIO = ps.STRAIN/ps.STEPS, 10.0             # strain rate per second (1 m cell), step ratio
WIDTH = 1/math.sqrt(2)                              # minimum altitude of the cell's two right triangles
V1 = 2*WIDTH*RATE1                                  # eps_II = e for normal rates +e, -e
T1, T2 = 64.0, 64.0                                 # seconds before and after the step


from des_case import run_solver, source_revision

def configs(refine):
    dt = 1.0/refine
    n1, n2 = int(T1*refine), int(T2*refine)
    base = ps.config(2, 'elasto-plastic-rate-state-friction', C, P, BIOT, 0)
    base['sim'].update(max_steps=str(n1), output_step_interval='1', checkpoint_frame_interval='1',
                       has_initial_checkpoint='true', max_time_in_yr='1')
    base['control'].update(fixed_dt=repr(dt), rsf_slip_rate_projection_option='1')
    base['monitor'] = dict(output_state_variable='true')
    base['bc'].update(vbc_val_x1=repr(RATE1), vbc_val_z1=repr(-RATE1))
    base['mat'].update(direct_a=f'[{A!r}]', evolution_b=f'[{B!r}]', characteristic_velocity=f'[{V1!r}]',
                       characteristic_distance=f'[{DC!r}]', state_var_model='1')
    step = configparser.ConfigParser(); step.optionxform = str
    step.read_dict(base)
    step['sim'].update(modelname='step', max_steps=str(n1 + n2), is_restarting='yes',
                       restarting_from_modelname='cell', restarting_from_frame=str(n1))
    step['bc'].update(vbc_val_x1=repr(RATIO*RATE1), vbc_val_z1=repr(-RATIO*RATE1))
    return base, step, dt, n1, n2


def discrete(dt, n1, n2):
    """DES's update, step by step: time, theta, mu and principal shear t after each step."""
    theta, t_load, out = DC/V1, 0.0, []
    for n in range(1, n1 + n2 + 1):
        rate = RATE1 if n <= n1 else RATIO*RATE1
        V = 2*WIDTH*rate
        theta = theta + (1 - V*theta/DC)*dt
        mu = MU0 + A*math.log(V/V1) + B*math.log(V1*theta/DC)
        phi = math.atan(mu)
        t_y = COHESION*math.cos(phi) + (C - BIOT*P)*math.sin(phi)
        t_load = min(t_load + 2*G*rate*dt, t_y)
        out.append((n*dt, theta, mu, t_load))
    return np.array(out)


def continuum(t):
    V2 = 2*WIDTH*RATIO*RATE1
    theta = np.where(t <= T1, DC/V1, DC/V2 + (DC/V1 - DC/V2)*np.exp(-V2*(t - T1)/DC))
    V = np.where(t <= T1, V1, V2)
    return MU0 + A*np.log(V/V1) + B*np.log(V1*theta/DC)


def read(folder, model_name, after=0.0):
    model = cell.Dynearthsol(str(folder/model_name)); rows = cell.scan_frames(str(folder/model_name))
    out = []
    for row in rows:
        if row['time'] <= after + 1e-9: continue      # a restart also saves its starting frame
        s = model.read_field(row['frame'], 'stress')
        t = float(np.mean((s[:, :2].max(axis=1) - s[:, :2].min(axis=1))/2))
        mu = model.read_field(row['frame'], 'dynamic friction coefficient')
        theta = model.read_field(row['frame'], 'friction state variable')
        out.append((row['time'], float(np.mean(theta)), float(np.mean(mu)), t,
                    float(np.ptp(mu)), float(np.ptp(theta))))
    return np.array(out)


def study(exe, out, refinements):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=cell.digest(exe), runner_sha256=cell.digest(Path(__file__)),
                  source_commit=source_revision(cell.ROOT),
                  a=A, b=B, dc=DC, v1=V1, ratio=RATIO, confinement=C, pressure=P, cases=[])
    for refine in refinements:
        folder = out/f'r{refine:g}'; folder.mkdir()
        base, step, dt, n1, n2 = configs(refine)
        for name, cfg in (('cell', base), ('step', step)):
            with (folder/f'{name}.cfg').open('w') as f: cfg.write(f)
            run_solver([str(exe), f'{name}.cfg'], log=folder / f'{name}.log', cwd=folder, env=env, check=True, timeout=3600)
        des = np.concatenate([read(folder, 'cell'), read(folder, 'step', T1)])
        ref = discrete(dt, n1, n2)
        assert len(des) == len(ref) and np.allclose(des[:, 0], ref[:, 0], rtol=0, atol=1e-9), (len(des), len(ref))
        exact = continuum(des[:, 0])
        after = des[:, 0] > T1
        item = dict(refine=refine, dt=dt, frames=len(des),
                    theta_vs_discrete=float(np.max(np.abs(des[:, 1]/ref[:, 1] - 1))),
                    mu_vs_discrete=float(np.max(np.abs(des[:, 2] - ref[:, 2]))),
                    stress_vs_discrete=float(np.max(np.abs(des[:, 3] - ref[:, 3]))/np.max(ref[:, 3])),
                    element_spread=float(max(des[:, 4].max(), des[:, 5].max()/des[:, 1].min())),
                    mu_vs_continuum=float(np.max(np.abs(des[after, 2] - exact[after]))),
                    mu_jump=float(A*math.log(RATIO)), mu_drop=float((A - B)*math.log(RATIO)),
                    series=dict(time=des[:, 0].tolist(), mu=des[:, 2].tolist(), stress=des[:, 3].tolist()))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'series'}), flush=True)
    shutil.copy2(__file__, out/'runner.py')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--refinements', type=float, nargs='+', default=[1, 2, 4])
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.refinements)
