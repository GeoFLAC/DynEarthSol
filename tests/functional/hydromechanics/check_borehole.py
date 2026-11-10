#!/usr/bin/env python3
"""DES fluid injection into a poroelastic borehole (plane strain).

A quarter annulus a = 1 m, R = 20 m (the Kirsch geometry, refined out to
5 m) is unstressed at p = 0. At t = 0+ the wall takes fluid at p_w: pore
pressure p_w (a 1 ms ramp of the scheduled boundary pressure) and normal
traction -p_w. The outer arc is sealed and traction-free (DES's input
check allows prescribed pressures on only one of two intersecting box
faces, which the annulus maps to x1 and z1); the symmetry cuts are
rollers and sealed. Gravity and heat are off; pressure feedback
is on. Compared with benchmarks/poromechanics/borehole_reference.py,
which includes the finite-annulus coupling: nodal pressure and element
radial/hoop total stress for 1.05 a <= r <= 4 a at four times, relative to
p_w.

    check_borehole.py EXECUTABLE OUTPUT [--spacings 0.2 0.1 0.05] [--steps 144]
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
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics'), str(HERE)]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from borehole_reference import Borehole
import check_kirsch as ck

A, R, RING = ck.A, ck.R, 5.0
PW, DURATION, RAMP = 1e6, 36.0, 1e-3
MAT = dict(bulk=2e8, shear=1e8, biot=1.0, porosity=.2, fluid_bulk=2e8, grain_bulk=1e20,
           permeability=1e-12, viscosity=1e-3)
REF = Borehole(radius=A, outer=R, bulk=MAT['bulk'], shear=MAT['shear'], biot=MAT['biot'],
               storage_modulus=1/(MAT['porosity']/MAT['fluid_bulk']),
               mobility=MAT['permeability']/MAT['viscosity'], wall_pressure=PW, sealed=True)
BAND = (1.05*A, 4*A)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config(h, steps):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='hole', max_steps=str(steps), max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval=str(steps//4), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='annulus.poly', xlength=repr(R), ylength='1', zlength=repr(R),
                       resolution=repr(h), remeshing_option='0', min_angle='28', max_steiner_factor='0')
    cfg['control'].update(gravity='0', has_thermal_diffusion='false', has_moving_mesh='false',
                          fixed_dt=repr(DURATION/steps), has_poroelastic_pressure_feedback='true',
                          is_using_mixed_stress='false', PT_option='1', PT_max_iter='200000',
                          PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc.update({f'vbc_{side}': '0', f'vbc_val_{side}': '0', f'stress_bc_{side}': '0', f'stress_val_{side}': '0',
                   f'hbc_{side}': '0'})
    bc.update(vbc_x0='1', vbc_z0='1', stress_bc_x1='4', stress_val_x1=repr(-PW), hbc_x1='2',
              hbc_pressure=f'[0,{PW!r},0,0,0,0]', num_hbc_schedule='2',
              hbc_schedule_time_in_yr=f'[0, {RAMP/YEAR2SEC!r}]', hbc_schedule_ratio='[0, 1]')
    cfg['ic'].update(pore_pressure_option='0', has_body_force_adjustment='false', weakzone_option='0')
    m = MAT
    cfg['mat'].update(rheology_type='elastic', is_plane_strain='true', alpha='[0]',
                      bulk_modulus=f"[{m['bulk']!r}]", shear_modulus=f"[{m['shear']!r}]", biot_coeff=f"[{m['biot']!r}]",
                      porosity=f"[{m['porosity']!r}]", fluid_bulk_modulus=f"[{m['fluid_bulk']!r}]",
                      bulk_modulus_s=f"[{m['grain_bulk']!r}]", hydraulic_perm=f"[{m['permeability']!r}]",
                      fluid_visc=f"[{m['viscosity']!r}]", fluid_alpha='[0]', retention_model='[0]')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'hole')); rows = scan_frames(str(folder/'hole'))
    x = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.linalg.norm(x, axis=1)
    centre = x[conn].mean(axis=1); re = np.linalg.norm(centre, axis=1)
    nx, nz = (centre/re[:, None]).T
    nodes, elems = (rn >= BAND[0]) & (rn <= BAND[1]), (re >= BAND[0]) & (re <= BAND[1])
    out = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        s = model.read_field(row['frame'], 'stress')
        srr = nx**2*s[:, 0] + nz**2*s[:, 1] + 2*nx*nz*s[:, 2]
        stt = nz**2*s[:, 0] + nx**2*s[:, 1] - 2*nx*nz*s[:, 2]
        fn, fe = REF.fields(rn[nodes], t), REF.fields(re[elems], t)
        out.append(dict(time=t, pressure_error=float(np.max(np.abs(p[nodes] - fn['p'])))/PW,
                        pressure_rms=float(np.sqrt(np.mean((p[nodes] - fn['p'])**2)))/PW,
                        radial_error=float(np.max(np.abs(srr[elems] - fe['srr'])))/PW,
                        hoop_error=float(np.max(np.abs(stt[elems] - fe['stt'])))/PW,
                        hoop_rms=float(np.sqrt(np.mean((stt[elems] - fe['stt'])**2)))/PW))
    return out, len(x), len(conn)


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  steps=steps, material=MAT, wall_pressure=PW, diffusivity=REF.c, cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        ck.polygon(folder/'annulus.poly', h, refined_radius=RING)
        with (folder/'hole.cfg').open('w') as f: config(h, steps).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'hole.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=40000)
        frames, nodes, elements = compare(folder)
        item = dict(spacing_m=h, nodes=nodes, elements=elements, frames=frames, input_sha256=digest(folder/'hole.cfg'),
                    mesh_sha256=digest(folder/'annulus.poly'),
                    pressure_error=max(f['pressure_error'] for f in frames),
                    stress_error=max(max(f['radial_error'], f['hoop_error']) for f in frames))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1, 0.05])
    p.add_argument('--steps', type=int, default=144)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
