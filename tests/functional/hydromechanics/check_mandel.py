#!/usr/bin/env python3
"""DES Mandel problem with a rigid frictionless platen, on refined meshes.

Quarter of the plane-strain strip: 0 <= x <= a with symmetry at x = 0, and
-b <= z <= 0 with symmetry at the mid-plane z = -b. The top carries the rigid
frictionless platen (bc.vbc_z1 = 8) loaded by the uniform traction -P, so the
force per unit length on the half width is F = P a. The side x = a is drained
(p = 0) and traction-free; the platen and the symmetry planes are sealed.
Gravity is off and pressure feedback is on (two-way coupling); PT is used for
the instantaneous load. The pressure and the uniform plate settlement are
compared with benchmarks/poromechanics/mandel_reference.py, which uses the
same material values.

    check_mandel.py EXECUTABLE OUTPUT [--spacings 0.1 0.05 0.025]
"""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'benchmarks'/'poromechanics'))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from mandel_reference import Mandel

A, HEIGHT, P = 1.0, 1.0, 1e5
MAT = dict(bulk=1e7, shear=1e7, biot=1.0, porosity=.3, fluid_bulk=2.17e9, grain_bulk=37e9,
           permeability=1.03e-12, viscosity=1.002e-3)
REF = Mandel(MAT['bulk'], MAT['shear'], MAT['biot'], MAT['porosity'], MAT['fluid_bulk'],
             MAT['grain_bulk'], MAT['permeability'], MAT['viscosity'], A, P*A)
TC = A**2/REF.c                      # consolidation time a^2/c


from des_case import digest, read_config, run_solver, source_diff, source_revision


def config(h, steps, duration):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='mandel', max_steps=str(steps), output_step_interval=str(max(1, steps//40)),
                      output_time_interval_in_yr='1e9', max_time_in_yr=repr(2*duration/YEAR2SEC),
                      checkpoint_frame_interval='10')
    cfg['mesh'].update(meshing_option='1', xlength=repr(A), ylength='1', zlength=repr(HEIGHT),
                       resolution=repr(h), largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(gravity='0', has_PT='true', PT_max_iter='100000', PT_relative_tolerance='1e-9',
                          PT_absolute_tolerance='1e-8', is_quasi_static='true', inertial_scaling='1e11',
                          has_moving_mesh='false', fixed_dt=repr(duration/steps), has_thermal_diffusion='false')
    bc = cfg['bc']
    bc.update(vbc_x0='1', vbc_val_x0='0', vbc_x1='0', vbc_z0='1', vbc_val_z0='0', vbc_z1='8',
              hbc_x0='0', hbc_x1='2', hbc_z0='0', hbc_z1='0', hbc_pressure='[0,0,0,0,0,0]',
              stress_bc_z1='3', stress_val_z1=repr(-P))
    cfg['ic'].update(has_body_force_adjustment='false', pore_pressure_option='0')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']}]", shear_modulus=f"[{MAT['shear']}]",
                      biot_coeff=f"[{MAT['biot']}]", porosity=f"[{MAT['porosity']}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']}]", bulk_modulus_s=f"[{MAT['grain_bulk']}]",
                      hydraulic_perm=f"[{MAT['permeability']}]", fluid_visc=f"[{MAT['viscosity']}]",
                      retention_model='[0]')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'mandel')); rows = scan_frames(str(folder/'mandel'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    x = coord[:, 0]
    tri = coord[conn]
    area = .5*np.abs((tri[:, 1, 0]-tri[:, 0, 0])*(tri[:, 2, 1]-tri[:, 0, 1])
                     - (tri[:, 2, 0]-tri[:, 0, 0])*(tri[:, 1, 1]-tri[:, 0, 1]))
    top = np.abs(coord[:, 1]) < 1e-9
    frames = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        ref = REF.pressure(x, np.full_like(x, t))
        ezz = model.read_field(row['frame'], 'strain')[:, 1]
        settle = -float(area @ ezz/area.sum())*HEIGHT   # each column integrates to the plate motion
        vz = model.read_field(row['frame'], 'velocity')[top, 1]
        frames.append(dict(time=t, t_over_tc=t/TC, max_error_over_p0=float(np.max(np.abs(p-ref)))/REF.p0,
                           centre_des_over_p0=float(np.mean(p[x < 1e-9]))/REF.p0,
                           centre_ref_over_p0=float(REF.pressure(0.0, t))/REF.p0,
                           settlement_m=settle, settlement_ref_m=float(REF.settlement(HEIGHT, t)),
                           platen_vz_spread=float(np.ptp(vz)/max(np.abs(vz).max(), 1e-300))))
    return model, rows, frames


def study(exe, out, spacings, steps, duration):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  source_commit=source_revision(ROOT),
                  reference=dict(c=REF.c, p0=REF.p0, B=REF.B, nu=REF.nu, nu_u=REF.nu_u, tc=TC, F=P*A,
                                 half_width=A, half_height=HEIGHT), steps=steps, duration_s=duration, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'mandel.cfg').open('w') as f: config(h, steps, duration).write(f)
        run_solver([str(exe), 'mandel.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=7200)
        _, rows, frames = compare(folder)
        settle = max(abs(f['settlement_m']-f['settlement_ref_m'])/f['settlement_ref_m'] for f in frames)
        item = dict(spacing_m=h, frames=frames, max_error_over_p0=max(f['max_error_over_p0'] for f in frames),
                    max_settlement_error=settle, max_platen_vz_spread=max(f['platen_vz_spread'] for f in frames),
                    centre_peak_des_over_p0=max(f['centre_des_over_p0'] for f in frames),
                    centre_peak_ref_over_p0=max(f['centre_ref_over_p0'] for f in frames),
                    input_sha256=digest(folder/'mandel.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.1, 0.05, 0.025])
    p.add_argument('--steps', type=int, default=200)
    p.add_argument('--duration', type=float, default=None, help='seconds; default one consolidation time a^2/c')
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps, a.duration or TC)
