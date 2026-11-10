#!/usr/bin/env python3
"""DES coupled gravity drainage of a partly saturating column (U9).

The column of benchmarks/poromechanics/drainage_reference.py: 1 m high,
laterally confined, saturated and hydrostatic with its water table at the
top; at t = 0+ the base pressure drops to zero (a 1 ms ramp of the scheduled
boundary pressure), the top stays sealed and traction-free. Water drains,
the top desaturates and the column settles under the Bishop stress and its
changing weight. Pressure at the nodes and the settlement (integrated layer
strain change) are compared with the independent reference at 300, 1200,
3600 and 7200 s, for PT (fixed-stress with the Bishop correction) and the
explicit strain form.

    check_unsaturated_drainage.py EXECUTABLE OUTPUT [--ndims 2] [--spacings 0.05 0.025 0.0125]
                                  [--steps 480] [--explicit]
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics')]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from drainage_reference import Drainage, liakopoulos_like

PROPS = liakopoulos_like()
REF = Drainage(**PROPS)
TIMES = (300.0, 1200.0, 3600.0, 7200.0)
DURATION, G, RHO = 7200.0, 9.81, 1000.0
RAMP = 1e-3                                          # s, base pressure ramp to zero


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config(ndims, h, steps):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    explicit = steps == 0
    cfg['sim'].update(modelname='drain', max_time_in_yr=repr(DURATION/YEAR2SEC), max_steps=str(10**9 if explicit else steps),
                      output_step_interval=str(10**9 if explicit else steps//24),
                      output_time_interval_in_yr=repr(300.0/YEAR2SEC) if explicit else '1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    H = PROPS['height']
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(h), ylength=repr(h), zlength=repr(H), resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h))
    control = dict(gravity=repr(G), has_thermal_diffusion='false', has_poroelastic_pressure_feedback='true',
                   PT_max_iter='100000')
    if explicit:
        # No boundary moves, so the explicit pseudo-mass needs an explicit characteristic speed.
        control.update(has_PT='false', inertial_scaling='1e5', fixed_dt='0', use_global_velocity_scaling='false',
                       characteristic_speed='1e-5', dt_fraction='0.5')
    else:
        control.update(fixed_dt=repr(DURATION/steps))
    cfg['control'].update(control)
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    bc.update(hbc_z0='2', hbc_z1='0', hbc_pressure=f'[0,0,0,0,{RHO*G*H!r},0]', num_hbc_schedule='2',
              hbc_schedule_time_in_yr=f'[0, {RAMP/YEAR2SEC!r}]', hbc_schedule_ratio='[1, 0]')
    cfg['ic'].update(water_table_z='0', pore_pressure_option='0')
    p = PROPS
    cfg['mat'].update(rho0='[2650]', bulk_modulus=f"[{p['bulk']!r}]", shear_modulus=f"[{p['shear']!r}]",
                      biot_coeff=f"[{p['biot']!r}]", porosity=f"[{p['porosity']!r}]",
                      fluid_bulk_modulus=f"[{p['fluid_bulk']!r}]", bulk_modulus_s=f"[{p['grain_bulk']!r}]",
                      hydraulic_perm=f"[{p['permeability']!r}]", fluid_visc=f"[{p['viscosity']!r}]",
                      fluid_rho0=f'[{RHO!r}]', fluid_alpha='[0]', retention_model='[1]',
                      retention_alpha=f"[{p['vg_alpha']!r}]", retention_n=f"[{p['vg_n']!r}]",
                      residual_saturation=f"[{p['sr']!r}]")
    return cfg


def layer_strain(coord, conn, ezz, h):
    depth = -coord[conn][:, :, -1].mean(axis=1)
    n = int(round(PROPS['height']/h)); layer = np.minimum((depth/h).astype(int), n - 1)
    count = np.bincount(layer, minlength=n).astype(float)
    assert count.min() > 0, 'the column mesh must be structured in layers of height h'
    return np.bincount(layer, weights=ezz, minlength=n)/count


def compare(folder, h, reference):
    model = Dynearthsol(str(folder/'drain')); rows = scan_frames(str(folder/'drain'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    zz = coord.shape[1] - 1; zn = coord[:, -1] + PROPS['height']
    e0 = layer_strain(coord, conn, model.read_field(0, 'strain')[:, zz], h)
    scale = RHO*G*PROPS['height']
    out = []
    for ref in reference:
        row = min(rows[1:], key=lambda r: abs(r['time'] - ref['time']))
        assert abs(row["time"] - ref["time"]) < 1.0, (row["time"], ref["time"])   # explicit runs pick their own step
        p = model.read_field(row['frame'], 'pore pressure')
        e = layer_strain(coord, conn, model.read_field(row['frame'], 'strain')[:, zz], h) - e0
        settle = float(np.sum(e)*h)
        pref = np.interp(zn, ref['z'], ref['p'])
        out.append(dict(time=ref['time'], pressure_error=float(np.max(np.abs(p - pref)))/scale,
                        settlement=settle, settlement_ref=ref['settlement'],
                        settlement_error=abs(settle/ref['settlement'] - 1), p_min=float(p.min())))
    return out


def study(exe, out, ndims, spacings, steps, explicit):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    reference = REF.solve(200, TIMES)
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  ndims=ndims, steps=steps, properties=PROPS, reference_grid=200,
                  reference=[dict(time=r['time'], settlement=r['settlement'], p_top=float(r['p'][-1]),
                                  saturation_top=float(r['saturation'][-1])) for r in reference], cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    runs = [(h, steps) for h in spacings] + ([(spacings[0], 0)] if explicit else [])
    for h, n in runs:
        folder = out/(f'h{h:g}' + ('-explicit' if n == 0 else '')); folder.mkdir()
        with (folder/'drain.cfg').open('w') as f: config(ndims, h, n).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'drain.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=40000)
        frames = compare(folder, h, reference)
        item = dict(spacing_m=h, mode='explicit' if n == 0 else 'PT', steps=n, frames=frames,
                    pressure_error=max(f['pressure_error'] for f in frames),
                    settlement_error=max(f['settlement_error'] for f in frames), input_sha256=digest(folder/'drain.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, default=2)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.05, 0.025, 0.0125])
    p.add_argument('--steps', type=int, default=480)
    p.add_argument('--explicit', action='store_true', help='also run the explicit strain form on the first mesh')
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.ndims, a.spacings, a.steps, a.explicit)
