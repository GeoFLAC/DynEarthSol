#!/usr/bin/env python3
"""DES Theis pumping and recovery in a confined aquifer, on refined meshes.

A quarter box (x, y >= 0) with sealed top, base and sides carries a fully
penetrating well on its corner edge, pumping a quarter of the total rate Q.
The skeleton is stiff and mechanical pressure feedback is off, so pressure
obeys S dp/dt = (kappa/mu) lap p with the fixed-stress storage
S = phi/K_f + (alpha - phi)/K_s + alpha^2/K. With gravity off, Theis gives
the pressure drawdown s_p = Q mu/(4 pi kappa B) E1(r^2 S mu/(4 kappa t)).
The well pumps, then shuts in; recovery is the superposed solution. The box
is wide enough that the drawdown at its far sides stays negligible.

    check_theis.py EXECUTABLE OUTPUT [--spacings 5 2.5 1.25]
"""
import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'benchmarks'/'poromechanics'))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from theis_reference import scheduled_drawdown

L, B, Q = 60.0, 5.0, 1e-4               # box side, thickness (m), total rate (m3/s)
PUMP, END, DT = 80.0, 160.0, 2.0        # shut-in time, end time, physical step (s)
MAT = dict(bulk=1e12, shear=1e12, biot=1.0, porosity=.25, fluid_bulk=2.2e9, grain_bulk=1e20, visc=1e-3)
STORAGE = (MAT['porosity']/MAT['fluid_bulk'] + (MAT['biot']-MAT['porosity'])/MAT['grain_bulk']
           + MAT['biot']**2/MAT['bulk'])
DIFFUSIVITY = 1.0                       # m2/s
PERM = DIFFUSIVITY*STORAGE*MAT['visc']
T = PERM*B/MAT['visc']                  # transmissivity in pressure units, m3/(Pa s)
S = STORAGE*B                           # storativity in pressure units, m/Pa
SCALE = Q/(4*math.pi*T)                 # pressure scale of the Theis solution, Pa
SCHEDULE = [(0.0, Q), (PUMP, 0.0)]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def config(h):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    steps = int(round(END/DT))
    cfg['sim'].update(modelname='theis', max_steps=str(steps), output_step_interval=str(steps//8),
                      max_time_in_yr=repr(2*END/YEAR2SEC), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='4')
    cfg['mesh'] = dict(meshing_option='1', meshing_elem_shape='1', xlength=repr(L), ylength=repr(L),
                       zlength=repr(B), resolution=repr(h), remeshing_option='0')
    cfg['control'].update(gravity='0', has_PT='false', is_quasi_static='true', inertial_scaling='1e6',
                          use_global_velocity_scaling='false', has_moving_mesh='false', damping_option='0',
                          fixed_dt=repr(DT), has_poroelastic_pressure_feedback='false',
                          has_thermal_diffusion='false')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    t = [0.0, PUMP, PUMP]
    bc.update(stress_bc_z1='0', stress_val_z1='0', num_wells='1', well_x='[0.0]', well_y='[0.0]',
              well_z_top='[0]', well_z_bottom=f'[{-B}]', well_rate=f'[{Q/4!r}]',
              num_well_schedule='3', well_schedule_time_in_yr='['+', '.join(repr(x/YEAR2SEC) for x in t)+']',
              well_schedule_ratio='[1, 1, 0]')
    cfg['ic'].update(has_body_force_adjustment='false')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']}]", shear_modulus=f"[{MAT['shear']}]",
                      biot_coeff=f"[{MAT['biot']}]", porosity=f"[{MAT['porosity']}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']}]", bulk_modulus_s=f"[{MAT['grain_bulk']}]",
                      hydraulic_perm=f'[{PERM!r}]', fluid_visc=f"[{MAT['visc']}]", retention_model='[0]')
    return cfg


def study(exe, out, spacings):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  reference=dict(L=L, B=B, Q=Q, pump_until_s=PUMP, end_s=END, dt_s=DT, storage_per_Pa=STORAGE,
                                 diffusivity=DIFFUSIVITY, permeability=PERM, transmissivity=T, storativity=S,
                                 pressure_scale_Pa=SCALE), cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'theis.cfg').open('w') as f: config(h).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'theis.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=3600)
        model = Dynearthsol(str(folder/'theis')); rows = scan_frames(str(folder/'theis'))
        coord = model.read_field(0, 'coordinate')
        r = np.hypot(coord[:, 0], coord[:, 1])
        window = (r >= 10) & (r <= L/2)        # same radii on every mesh: no axis, no far field
        frames = []
        for row in rows[1:]:
            t = row['time']
            p = model.read_field(row['frame'], 'pore pressure')
            p0 = model.read_field(0, 'pore pressure')
            des = p0 - p                          # drawdown: pressure decrease
            ref = scheduled_drawdown(r[window], np.full(window.sum(), t), SCHEDULE, T, S)
            far = scheduled_drawdown(np.array([L]), np.array([t]), SCHEDULE, T, S)[0]
            frames.append(dict(time=t, max_error_over_scale=float(np.max(np.abs(des[window]-ref)))/SCALE,
                               rms_error_over_scale=float(np.sqrt(np.mean((des[window]-ref)**2)))/SCALE,
                               far_side_reference_over_scale=far/SCALE))
        strain = max(float(np.abs(model.read_field(f, 'strain')).max()) for f in model.frames)
        item = dict(spacing_m=h, nodes=len(r), window=[10, L/2], frames=frames,
                    max_error_over_scale=max(f['max_error_over_scale'] for f in frames),
                    rms_error_over_scale=max(f['rms_error_over_scale'] for f in frames),
                    max_abs_strain=strain, input_sha256=digest(folder/'theis.cfg'))
        assert strain < 1e-6, ('stiff skeleton must stay rigid', strain)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[5, 2.5, 1.25])
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings)
