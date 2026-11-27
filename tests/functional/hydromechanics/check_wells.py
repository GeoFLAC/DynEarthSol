#!/usr/bin/env python3
"""Check rate-controlled wells with the production executable.

Usage: check_wells.py EXECUTABLE OUTPUT_DIRECTORY [--ndims 3]
In a sealed domain the fluid-volume budget must close within 1e-9 of the
scheduled pumped volume:
    sum_n hmass_n (p_n - p0_n) + sum_e alpha V_e (s_e - s0_e) / K_b = -pumped volume,
where s is the recorded mean total stress. A centred well and a well on a
domain edge, which several elements share, follow a schedule that pumps, shuts
in and injects, with events between time steps. OUTPUT_DIRECTORY must not exist.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

DURATION = 2.0e4                  # s
EVENTS = (0.33 * DURATION, 0.61 * DURATION)
MAT = dict(bulk=1e8, shear=1e8, biot=0.8, porosity=0.25, fluid_bulk=2.2e9,
           permeability=1e-12, viscosity=1e-3)
MAT['grain_bulk'] = MAT['bulk'] / (1 - MAT['biot'])
BALANCE_RTOL = 1e-9


from check_fluid_remap import configuration as remesh_configuration
from des_case import material_options, mean_stress_modulus, read_config, restart_check, run_des, write_json
from des_test_utils import checkpoint_arrays, lumped_volumes, simplex_volumes

def config(dim, edge):
    cfg = read_config(ROOT / f'benchmarks/poromechanics/1d-consolidation-des{dim}d_traction.cfg')
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
    cfg['mat'].update(material_options(MAT))
    return cfg, Q


def pumped(Q, t):
    """Exact volume of the schedule [1 until t1, 0 until t2, -0.5 after]."""
    return Q * (min(t, EVENTS[0]) - 0.5 * max(0.0, t - EVENTS[1]))


def balance(case, dim, Q):
    model = Dynearthsol(str(case / 'wells'))
    assert model.ndims == dim, 'the executable and --ndims disagree'
    times = {r['frame']: r['time'] for r in scan_frames(str(case / 'wells'))}
    alpha, K, G = MAT['biot'], MAT['bulk'], MAT['shear']
    Kb = mean_stress_modulus(dim, K, G)
    S = (MAT['porosity'] / MAT['fluid_bulk'] + (alpha - MAT['porosity']) / MAT['grain_bulk']
         + alpha**2 / Kb)
    first = model.frames[0]
    xyz, conn = (model.read_field(first, key) for key in ('coordinate', 'connectivity'))
    vol = simplex_volumes(xyz, conn)
    mass = lumped_volumes(xyz, conn, S * vol)

    def state(frame):
        stress = model.read_field(frame, 'stress')
        return model.read_field(frame, 'pore pressure'), stress[:, :dim].sum(1) / dim

    p0, s0 = state(first)
    worst = 0.0
    for frame in model.frames[1:]:
        p, s = state(frame)
        stored = np.sum(mass * (p - p0)) + np.sum(alpha * vol * (s - s0) / Kb)
        scale = abs(Q) * EVENTS[0]
        worst = max(worst, abs(stored + pumped(Q, times[frame])) / scale)
    return worst, len(model.frames)


def trimmed_screen(exe, folder, dim):
    cfg = remesh_configuration('domain-remove')
    rate, dt = 1e-8, cfg.getfloat('control', 'fixed_dt')
    cfg['sim']['has_output_during_remeshing'] = 'no'
    cfg['bc'].update(num_wells='1', well_x='[0.5]', well_y='[0.5]',
                     well_z_bottom='[-1]', well_z_top='[-0.2]', well_rate=f'[{rate}]')
    restart_check(exe, cfg, folder, steps=3, restart_frame=2, require_remesh=True)
    model = Dynearthsol(str(folder/'fresh'))
    assert model.ndims == dim
    for frame in model.frames[1:]:
        balance = model.read_field(frame, 'fluid balance')
        np.testing.assert_allclose(balance[2], rate*dt, rtol=BALANCE_RTOL, atol=0)
        keys = ('well piece elements', 'well piece starts')
        checkpoint = checkpoint_arrays(folder/f'fresh.chkpt.{frame:06d}', keys)
        conn = model.read_field(frame, 'connectivity')
        xyz = model.read_field(frame, 'coordinate')
        z = checkpoint[keys[1]][:dim+1] @ xyz[conn[checkpoint[keys[0]][0]], -1]
        assert z < -1., ('the screen must reach the receding base', frame, z)
        np.testing.assert_allclose(z, xyz[:, -1].min(), rtol=0, atol=1e-10)
    return dict(remeshed=True, exact_restart=True, well_rate=rate)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    a = p.parse_args()
    exe, output, dim = a.executable.resolve(), a.output.resolve(), a.ndims
    output.mkdir(parents=True, exist_ok=False)
    summary, ok = {}, True
    for edge in (False, True):
        name = 'edge' if edge else 'centre'
        case = output / name
        case.mkdir()
        cfg, Q = config(dim, edge)
        run_des(exe, case, 'input', cfg, timeout=900)
        worst, frames = balance(case, dim, Q)
        summary[name] = dict(frames=frames, max_balance_error_over_pumped=worst)
        ok &= worst < BALANCE_RTOL
    summary['trimmed_screen'] = trimmed_screen(exe, output/'trimmed-screen', dim)
    print(json.dumps(summary, indent=2))
    write_json(output / 'summary.json', summary)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
