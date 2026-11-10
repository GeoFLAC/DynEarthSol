#!/usr/bin/env python3
"""DES Cryer problem: a poroelastic sphere under a sudden external pressure.

One octant (x, y >= 0, z <= 0) of a sphere of radius R is meshed from a .poly
file. The three flat faces are symmetry planes (normal velocity zero, sealed);
the curved surface carries the BOUNDZ1 flag, so it is drained (hbc_z1 = 2,
p = 0) and loaded by the normal stress -P (bc.stress_bc_z1 = 4). Gravity is
off and PT handles the instantaneous load. The pressure is compared with
benchmarks/poromechanics/cryer_reference.py, which uses the same material.

    check_cryer.py EXECUTABLE OUTPUT [--spacings 0.2 0.1 0.07]
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
from cryer_reference import Cryer

R, P = 1.0, 1e5
MAT = dict(bulk=1e7, shear=1e7, biot=1.0, porosity=.3, fluid_bulk=2.17e9, grain_bulk=37e9,
           permeability=1.03e-12, viscosity=1.002e-3)
REF = Cryer(MAT['bulk'], MAT['shear'], MAT['biot'], MAT['porosity'], MAT['fluid_bulk'], MAT['grain_bulk'],
            MAT['permeability'], MAT['viscosity'], R, P)
TC = R**2/REF.c
BOUNDX0, BOUNDY0, BOUNDZ0, BOUNDZ1 = 1, 4, 16, 32


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_octant_poly(path, h):
    """Octant of a sphere: a subdivided spherical triangle plus three quarter disks."""
    n = max(2, math.ceil(.5*math.pi*R/h))
    A, B, C = np.array([R, 0, 0.]), np.array([0, R, 0.]), np.array([0, 0, -R])
    index, points = {}, [np.zeros(3)]                       # point 0 is the centre
    for i in range(n+1):
        for j in range(n+1-i):
            k = n-i-j
            v = i*A + j*B + k*C
            index[i, j] = len(points); points.append(R*v/np.linalg.norm(v))
    facets = []
    for i in range(n):
        for j in range(n-i):
            facets.append((BOUNDZ1, [index[i, j], index[i+1, j], index[i, j+1]]))
            if i+j < n-1:
                facets.append((BOUNDZ1, [index[i+1, j], index[i+1, j+1], index[i, j+1]]))
    arc_z = [index[i, n-i] for i in range(n+1)]             # k = 0: plane z = 0, from B to A
    arc_y = [index[i, 0] for i in range(n+1)]               # j = 0: plane y = 0, from C to A
    arc_x = [index[0, j] for j in range(n+1)]               # i = 0: plane x = 0, from C to B
    facets += [(BOUNDZ0, [0] + arc_z), (BOUNDY0, [0] + arc_y), (BOUNDX0, [0] + arc_x)]
    with path.open('w') as f:
        f.write(f'# Octant of a sphere, R = {R}, surface subdivision {n}\n{len(points)} 3 0 0\n')
        for k, p in enumerate(points):
            f.write(f'{k} {p[0]:.16e} {p[1]:.16e} {p[2]:.16e}\n')
        f.write(f'{len(facets)} 1\n')
        for marker, v in facets:
            f.write(f'1 0 {marker}\n{len(v)} ' + ' '.join(map(str, v)) + '\n')
        f.write('0\n1\n')
        f.write(f'0 {R/4:.16e} {R/4:.16e} {-R/4:.16e} 0 {h**3/(6*math.sqrt(2)):.16e}\n')
    return n


def config(h, steps, duration):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    cfg['sim'].update(modelname='cryer', max_steps=str(steps), output_step_interval=str(max(1, steps//40)),
                      output_time_interval_in_yr='1e9', max_time_in_yr=repr(2*duration/YEAR2SEC),
                      checkpoint_frame_interval='10')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='octant.poly', xlength=repr(R), ylength=repr(R),
                       zlength=repr(R), resolution=repr(h), remeshing_option='0')
    cfg['control'].update(gravity='0', has_PT='true', PT_max_iter='100000', PT_relative_tolerance='1e-9',
                          PT_absolute_tolerance='1e-8', is_quasi_static='true', inertial_scaling='1e11',
                          has_moving_mesh='false', fixed_dt=repr(duration/steps), has_thermal_diffusion='false')
    bc = cfg['bc']
    for side in ('x0', 'y0', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    for side in ('x1', 'y1', 'z1'):
        bc[f'vbc_{side}'], bc[f'hbc_{side}'] = '0', '0'
    bc.update(hbc_z1='2', hbc_pressure='[0,0,0,0,0,0]', stress_bc_z1='4', stress_val_z1=repr(-P))
    cfg['ic'].update(has_body_force_adjustment='false', pore_pressure_option='0')
    cfg['mat'].update(bulk_modulus=f"[{MAT['bulk']}]", shear_modulus=f"[{MAT['shear']}]",
                      biot_coeff=f"[{MAT['biot']}]", porosity=f"[{MAT['porosity']}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']}]", bulk_modulus_s=f"[{MAT['grain_bulk']}]",
                      hydraulic_perm=f"[{MAT['permeability']}]", fluid_visc=f"[{MAT['viscosity']}]",
                      retention_model='[0]')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'cryer')); rows = scan_frames(str(folder/'cryer'))
    coord = model.read_field(0, 'coordinate')
    r = np.linalg.norm(coord, axis=1)
    centre = int(np.argmin(r))
    frames = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        ref = REF.pressure(r, np.full_like(r, t))
        frames.append(dict(time=t, t_over_tc=t/TC, max_error_over_p0=float(np.max(np.abs(p-ref)))/REF.p0,
                           rms_error_over_p0=float(np.sqrt(np.mean((p-ref)**2)))/REF.p0,
                           centre_des_over_p0=float(p[centre])/REF.p0,
                           centre_ref_over_p0=float(REF.pressure(0.0, t))/REF.p0))
    return frames, len(r), float(r[centre])


def study(exe, out, spacings, steps, duration):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  reference=dict(c=REF.c, p0=REF.p0, beta=REF.beta, tc=TC, radius=R, pressure=P),
                  steps=steps, duration_s=duration, cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        n = write_octant_poly(folder/'octant.poly', h)
        with (folder/'cryer.cfg').open('w') as f: config(h, steps, duration).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'cryer.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=14400)
        frames, nodes, r0 = compare(folder)
        item = dict(spacing_m=h, surface_subdivision=n, nodes=nodes, centre_node_radius=r0, frames=frames,
                    max_error_over_p0=max(f['max_error_over_p0'] for f in frames),
                    rms_error_over_p0=max(f['rms_error_over_p0'] for f in frames),
                    centre_peak_des_over_p0=max(f['centre_des_over_p0'] for f in frames),
                    centre_peak_ref_over_p0=max(f['centre_ref_over_p0'] for f in frames),
                    input_sha256=digest(folder/'cryer.cfg'), poly_sha256=digest(folder/'octant.poly'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1, 0.07])
    p.add_argument('--steps', type=int, default=200)
    p.add_argument('--duration', type=float, default=None, help='seconds; default 0.5 R^2/c')
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps, a.duration or .5*TC)
