#!/usr/bin/env python3
"""DES Rudnicki problem: a continuous fluid point source in a poroelastic solid.

One octant (x, y >= 0, z <= 0) of a large box: symmetry planes x = 0, y = 0 and
z = 0 (normal displacement zero, sealed); far faces with zero normal
displacement and p = 0. A scheduled point source at the corner node injects
Q/8 at a constant rate from t = 0. The mesh is refined near the source
(mesh.meshing_option = 2) and coarse outside. Gravity is off; pressure
feedback is on and PT keeps mechanical equilibrium. Pressure and the radial
and hoop strains are compared with benchmarks/poromechanics/rudnicki_reference.py
for 2.5 <= r <= 8 m while the diffusion front stays well inside the box.

    check_rudnicki.py EXECUTABLE OUTPUT [--spacings 1 0.7 0.5] [--steps 120]
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
from rudnicki_reference import Rudnicki

L, ZONE, Q, DURATION = 60.0, 10.0, 1e-4, 240.0      # box, refined zone (m), total rate, end time
BAND = (2.5, 8.0)
MAT = dict(bulk=2e8, shear=1e8, biot=.8, porosity=.25, fluid_bulk=2.2e9, grain_bulk=1e9,
           permeability=1e-13, viscosity=1e-3)
REF = Rudnicki(MAT['bulk'], MAT['shear'], MAT['biot'], MAT['porosity'], MAT['fluid_bulk'],
               MAT['grain_bulk'], MAT['permeability'], MAT['viscosity'], Q)


from des_case import digest, read_config, run_solver, source_diff, source_revision


def config(h, steps):
    cfg = read_config(ROOT / 'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    cfg['sim'].update(modelname='rudnicki', max_steps=str(steps), output_step_interval=str(max(1, steps//12)),
                      output_time_interval_in_yr='1e9', max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      checkpoint_frame_interval='6')
    zone = f'[0, {ZONE/L!r}]'
    cfg['mesh'] = dict(meshing_option='2', xlength=repr(L), ylength=repr(L), zlength=repr(L), resolution=repr(h),
                       refined_zonex=zone, refined_zoney=zone, refined_zonez=zone,
                       largest_size='1000', smallest_size='0.01', remeshing_option='0')
    cfg['control'].update(gravity='0', has_PT='true', PT_max_iter='100000', PT_relative_tolerance='1e-9',
                          PT_absolute_tolerance='1e-8', is_quasi_static='true', inertial_scaling='1e11',
                          has_moving_mesh='false', fixed_dt=repr(DURATION/steps), has_thermal_diffusion='false')
    bc = cfg['bc']
    for side in ('x0', 'y0', 'z1'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    for side in ('x1', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '2'
    bc.update(hbc_pressure='[0,0,0,0,0,0]', stress_bc_z1='0', stress_val_z1='0')
    cfg['ic'].update(has_body_force_adjustment='false', pore_pressure_option='0')
    cfg['injection'] = dict(enabled='true', num_points='1', points_x='[0]', points_y='[0]', points_z='[0]',
                            points_unit='m', rate_model='constant_rate', rate=f'[{Q/8!r}]',
                            start_time_in_yr='[0]', end_time_in_yr=f'[{2*DURATION/YEAR2SEC!r}]')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']}]", shear_modulus=f"[{MAT['shear']}]",
                      biot_coeff=f"[{MAT['biot']}]", porosity=f"[{MAT['porosity']}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']}]", bulk_modulus_s=f"[{MAT['grain_bulk']}]",
                      hydraulic_perm=f"[{MAT['permeability']}]", fluid_visc=f"[{MAT['viscosity']}]",
                      retention_model='[0]')
    return cfg


def polar_strains(coord, conn, strain):
    """Radial and hoop strain at element centroids from the Cartesian strain (xx, yy, zz, xy, xz, yz)."""
    c = coord[conn].mean(axis=1)
    r = np.linalg.norm(c, axis=1); n = c/r[:, None]
    e = np.zeros((len(c), 3, 3))
    e[:, 0, 0], e[:, 1, 1], e[:, 2, 2] = strain[:, 0], strain[:, 1], strain[:, 2]
    e[:, 0, 1] = e[:, 1, 0] = strain[:, 3]; e[:, 0, 2] = e[:, 2, 0] = strain[:, 4]; e[:, 1, 2] = e[:, 2, 1] = strain[:, 5]
    err = np.einsum('ei,eij,ej->e', n, e, n)
    return r, err, (np.trace(e, axis1=1, axis2=2) - err)/2


def compare(folder):
    model = Dynearthsol(str(folder/'rudnicki')); rows = scan_frames(str(folder/'rudnicki'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.linalg.norm(coord, axis=1); node = (rn >= BAND[0]) & (rn <= BAND[1])
    frames = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        re, err, ett = polar_strains(coord, conn, model.read_field(row['frame'], 'strain'))
        elem = (re >= BAND[0]) & (re <= BAND[1])
        _, rr, tt = REF.strains(re[elem], t)
        pref = REF.pressure(rn[node], t)
        pscale, escale = np.abs(pref).max(), max(np.abs(rr).max(), np.abs(tt).max())
        frames.append(dict(time=t, pressure_error=float(np.max(np.abs(p[node]-pref))/pscale),
                           pressure_rms=float(np.sqrt(np.mean((p[node]-pref)**2))/pscale),
                           radial_rms=float(np.sqrt(np.mean((err[elem]-rr)**2))/escale),
                           hoop_rms=float(np.sqrt(np.mean((ett[elem]-tt)**2))/escale)))
    return frames, len(coord), len(conn)


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  source_commit=source_revision(ROOT),
                  reference=dict(c=REF.c, S=REF.S, A=REF.A, rate=Q, box=L, refined_zone=ZONE, band=BAND),
                  steps=steps, duration_s=DURATION, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'rudnicki.cfg').open('w') as f: config(h, steps).write(f)
        run_solver([str(exe), 'rudnicki.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=20000)
        frames, nodes, elements = compare(folder)
        item = dict(spacing_m=h, nodes=nodes, elements=elements, frames=frames,
                    pressure_error=max(f['pressure_error'] for f in frames),
                    pressure_rms=max(f['pressure_rms'] for f in frames),
                    radial_rms=max(f['radial_rms'] for f in frames), hoop_rms=max(f['hoop_rms'] for f in frames),
                    input_sha256=digest(folder/'rudnicki.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[1, 0.7, 0.5])
    p.add_argument('--steps', type=int, default=120)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
