#!/usr/bin/env python3
"""Rain on a partly saturated slope builds a perched water table until it slides.

The slope of slope_application.py starts at a uniform suction (Gardner
retention, Bishop chi = S_e). Rain enters the surface at a constant rate; the
bedrock and crest end are sealed and the toe end drains. Failure is read at a
fixed basal plastic strain and compared with the bedrock-pressure loading of
the same slope and with a fill-the-pore-space time estimate.

    rainfall_slope.py EXECUTABLE OUTPUT [--spacings 0.5 0.25] [--steps 960]
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np

import slope_application as sa
from unsaturated_reference import Retention

RAIN, SUCTION, ALPHA_G, SR = 1e-5, 5e3, 2.4e-4, .05      # m/s, Pa, 1/Pa (Gardner), residual saturation
COHESION, FRICTION = 1e3, math.radians(30.0)             # weaker than the pressure-driven slope: fails before saturating
U_POINT, U_DAVIS, U_STAR = sa.estimates(COHESION, FRICTION)
DURATION = 1.2e5                                            # s


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import source_revision

se = Retention('gardner', ALPHA_G, sr=SR).se


def config(h, steps, loaded):
    cfg = sa.config(h, steps, loaded, cohesion=COHESION, friction=FRICTION)
    cfg['sim'].update(max_time_in_yr=repr(2*DURATION/sa.YEAR2SEC), output_step_interval=str(max(1, steps//240)))
    cfg['control']['fixed_dt'] = repr(DURATION/steps)
    bc = cfg['bc']
    bc.update(hbc_z0='0', hbc_z1='3', hbc_x1='2', hbc_x0='0', hbc_pressure='[0,0,0,0,0,0]',
              hbc_flux=f'[0,0,0,0,0,{-RAIN!r}]', num_hbc_schedule='1', hbc_schedule_time_in_yr='[0]', hbc_schedule_ratio='[1]')
    cfg['ic'].update(pore_pressure_option='2', pore_pressure_value=repr(-SUCTION))
    for key in ('retention_model', 'retention_alpha', 'residual_saturation'):
        cfg['mat'].pop(key, None)
    cfg['mat'].update(retention_model='[3, 3]', retention_alpha=f'[{ALPHA_G!r}, {ALPHA_G!r}]',
                      residual_saturation=f'[{SR!r}, {SR!r}]')
    return cfg


def surface_pressure(folder):
    """Largest surface pore pressure at each saved step. A positive value means rain is being forced into a
    saturated surface: DES has no runoff or ponding, so results after it are not physical."""
    model = sa.Dynearthsol(str(folder/'slope')); rows = sa.scan_frames(str(folder/'slope'))
    x = model.read_field(0, 'coordinate')
    surface = np.abs(x[:, 1] + x[:, 0]*math.tan(sa.BETA)) < 1e-9
    return [float(model.read_field(r['frame'], 'pore pressure')[surface].max()) for r in rows if r['time'] > 0]


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    se0 = float(se(-SUCTION)); s0 = SR + (1 - SR)*se0
    # Lower bound on the time to failure: rain fills the pore space above the initial saturation up to the
    # perched height at which the base pressure (slope-parallel seepage, u = rho_f g h_w cos^2 beta) reaches
    # the pressure-driven collapse value; no water leaves.
    result = dict(executable=str(exe), executable_sha256=sa.digest(exe), runner_sha256=sa.digest(__file__),
                  slope_runner_sha256=sa.digest(sa.__file__),
                  source_commit=source_revision(sa.ROOT),
                  rain=RAIN, suction=SUCTION, alpha=ALPHA_G, sr=SR, se0=se0, s0=s0, bishop0=-se0*SUCTION,
                  cohesion=COHESION, friction_deg=math.degrees(FRICTION),
                  u_point=U_POINT, u_davis=U_DAVIS, u_star=U_STAR, steps=steps, cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        sa.polygon(folder/'slope.poly', h)
        for name, loaded in (('seed', False), ('slope', True)):
            with (folder/f'{name}.cfg').open('w') as f: config(h, steps, loaded).write(f)
            with (folder/f'{name}.log').open('w') as f:
                code = subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, env=env, stdout=f,
                                      stderr=subprocess.STDOUT, timeout=80000).returncode
            if code not in (0, 54) or (code == 54 and not loaded):
                raise RuntimeError(f'{folder}/{name}: exit {code}')
        frames, nodes, elements = sa.analyse(folder)
        for f, ps in zip(frames, surface_pressure(folder)): f['surface_pressure_max'] = ps
        item = dict(spacing_m=h, nodes=nodes, elements=elements, lost_equilibrium=code == 54, frames=frames)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        last = frames[-1]
        print(json.dumps(dict(spacing_m=h, lost_equilibrium=code == 54, last_time=last['time'],
                              last_base_pressure=last['base_pressure'],
                              last_basal_plastic_fraction=last['basal_plastic_fraction'])), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.5, 0.25])
    p.add_argument('--steps', type=int, default=960)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
