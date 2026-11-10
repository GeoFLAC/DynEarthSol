#!/usr/bin/env python3
"""Check rate-controlled wells with the production executable.

Usage: check_wells.py EXECUTABLE OUTPUT_DIRECTORY
In a sealed domain the explicit pressure update conserves fluid exactly:
    sum_n hmass_n (p_n - p0_n) + sum_e alpha V_e (s_e - s0_e) / K_b = -pumped volume,
where s is the recorded mean total stress. A centred well and a well on a
domain edge, which several elements share, follow a schedule that pumps, shuts
in and injects, with events between time steps. The output directory is new.
"""

import configparser
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

DURATION = 2.0e4                  # s
EVENTS = (0.33 * DURATION, 0.61 * DURATION)
MAT = dict(bulk=1e8, shear=1e8, biot=0.8, porosity=0.25, fluid_bulk=2.2e9,
           perm=1e-12, visc=1e-3)
MAT['grain_bulk'] = MAT['bulk'] / (1 - MAT['biot'])


def config(dim, edge):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(ROOT / f'benchmarks/poromechanics/1d-consolidation-des{dim}d_traction.cfg')
    L, H, Q = 100.0, 10.0, (1e-3 if dim == 3 else 1e-4)
    t = [0, EVENTS[0], EVENTS[0], EVENTS[1], EVENTS[1]]
    cfg['sim'].update(modelname='wells', max_time_in_yr=repr(DURATION / YEAR2SEC),
                      output_time_interval_in_yr=repr(DURATION / 8 / YEAR2SEC))
    cfg['mesh'] = dict(meshing_option='1', xlength=repr(L), ylength=repr(L), zlength=repr(H),
                       resolution='5', remeshing_option='0')
    cfg['control'].update(gravity='0', PT_option='1', PT_absolute_tolerance='1e-6')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    bc.update(hbc_z1='0', stress_bc_z1='0', stress_val_z1='0', num_wells='1',
              well_x=f'[{0.0 if edge else L / 2}]', well_y=f'[{0.0 if edge else L / 2}]',
              well_z_top='[0]', well_z_bottom=f'[{-H}]', well_rate=f'[{Q}]',
              num_well_schedule='5',
              well_schedule_time_in_yr='[' + ', '.join(repr(x / YEAR2SEC) for x in t) + ']',
              well_schedule_ratio='[1, 1, 0, 0, -0.5]')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']}]", shear_modulus=f"[{MAT['shear']}]",
                      biot_coeff=f"[{MAT['biot']}]", porosity=f"[{MAT['porosity']}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']}]",
                      bulk_modulus_s=f"[{MAT['grain_bulk']}]",
                      hydraulic_perm=f"[{MAT['perm']}]", fluid_visc=f"[{MAT['visc']}]")
    return cfg, Q


def pumped(Q, t):
    """Exact volume of the schedule [1 until t1, 0 until t2, -0.5 after]."""
    return Q * (min(t, EVENTS[0]) - 0.5 * max(0.0, t - EVENTS[1]))


def balance(case, dim, Q):
    model = Dynearthsol(str(case / 'wells'))
    times = {r['frame']: r['time'] for r in scan_frames(str(case / 'wells'))}
    alpha, K, G = MAT['biot'], MAT['bulk'], MAT['shear']
    Kb = K if dim == 3 else K + G / 3
    S = (MAT['porosity'] / MAT['fluid_bulk'] + (alpha - MAT['porosity']) / MAT['grain_bulk']
         + alpha**2 / Kb)
    conn = model.read_field(model.frames[0], 'connectivity')

    def state(frame):
        xyz = model.read_field(frame, 'coordinate')
        cells = xyz[conn]
        vol = np.abs(np.linalg.det(cells[:, 1:] - cells[:, :1])) / math.factorial(dim)
        hmass = np.zeros(len(xyz))
        np.add.at(hmass, conn, (S * vol / (dim + 1))[:, None])
        stress = model.read_field(frame, 'stress')
        return hmass, vol, model.read_field(frame, 'pore pressure'), stress[:, :dim].sum(1) / dim

    h0, v0, p0, s0 = state(model.frames[0])
    worst = 0.0
    for frame in model.frames[1:]:
        h, v, p, s = state(frame)
        stored = np.sum(h0 * (p - p0)) + np.sum(alpha * v0 * (s - s0) / Kb)
        scale = abs(Q) * EVENTS[0]
        worst = max(worst, abs(stored + pumped(Q, times[frame])) / scale)
    return worst, len(model.frames)


def main():
    exe, output = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    dim = 3 if exe.name.endswith('3d') else 2
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    summary, ok = {}, True
    for edge in (False, True):
        name = 'edge' if edge else 'centre'
        case = output / name
        case.mkdir()
        cfg, Q = config(dim, edge)
        with (case / 'input.cfg').open('w') as stream:
            cfg.write(stream)
        with (case / 'run.log').open('w') as stream:
            subprocess.run([str(exe), 'input.cfg'], cwd=case, env=env, check=True,
                           stdout=stream, stderr=subprocess.STDOUT, timeout=900)
        worst, frames = balance(case, dim, Q)
        summary[name] = dict(frames=frames, max_balance_error_over_pumped=worst)
        ok &= worst < 1e-4
    print(json.dumps(summary, indent=2))
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
