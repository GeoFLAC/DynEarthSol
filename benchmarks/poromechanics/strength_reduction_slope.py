#!/usr/bin/env python3
"""Factor of safety of a homogeneous embankment by strength reduction.

H = 10 m, 2:1 face, foundation H deep; c = 10 kPa, phi = 20 deg, psi = 0, dry.
Gravity is applied elastically by two PT steps from an unstressed seed. Each
trial restarts with c/F and atan(tan(phi)/F) and runs explicitly for DURATION
seconds on the fixed mesh; it is at rest if its peak speed at the end is below
that three quarters through. Reference: Bishop's method (bishop_reference.py).

    strength_reduction_slope.py EXECUTABLE OUTPUT [--spacings 1 0.5] [--factors 1.0 1.2 ...]
"""
import argparse
import configparser
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
import bishop_reference as br

G_ACC = 9.81
RHO = br.GAMMA/G_ACC
E, NU = 1e8, 0.3
BULK, SHEAR = E/(3*(1 - 2*NU)), E/(2*(1 + NU))
WIDTH = br.TOE + br.RUN + br.CREST
GRAVITY_DT = 0.02                         # s per metre of mesh spacing: PT gravity step, below the explicit limit
DURATION, FRAMES = 2000.0, 200            # s of mass-scaled explicit time per trial, saved frames
SCALING, DAMPING = 1e5, 0.8               # inertial scaling; DES's default local damping factor
FACTORS = [1.0, 1.2, 1.3, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.45]


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def polygon(path, h):
    """Counter-clockwise outline: base, right side, crest, face, toe platform, left side."""
    def seg(a, b):
        n = max(1, math.ceil(math.dist(a, b)/h))
        return [(a[0] + (b[0] - a[0])*k/n, a[1] + (b[1] - a[1])*k/n) for k in range(n)]
    corners = [(0.0, br.BASE), (WIDTH, br.BASE), (WIDTH, 0.0), (br.TOE + br.RUN, 0.0), (br.TOE, -br.H), (0.0, -br.H)]
    marks = [16, 2, 32, 32, 32, 1]                      # z0, x1, surface (crest, face, toe), x0
    points, markers = [], []
    for i, (a, mk) in enumerate(zip(corners, marks)):
        piece = seg(a, corners[(i + 1) % len(corners)])
        points += piece; markers += [mk]*len(piece)
    with path.open('w') as f:
        f.write(f'{len(points)} 2 0 0\n')
        for i, (x, z) in enumerate(points): f.write(f'{i} {x:.17g} {z:.17g}\n')
        f.write(f'{len(points)} 1\n')
        for i, mk in enumerate(markers): f.write(f'{i} {i} {(i + 1) % len(points)} {mk}\n')
        f.write(f'0\n1\n0 {WIDTH/2} {br.BASE/2} 0 {h*h/2}\n')


def config(h, name, factor=None):
    """seed: no gravity, zero stress, zero steps; gravity: two elastic PT steps; otherwise the explicit reduced-strength run."""
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname=name, output_step_interval='1', output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='1' if name == 'gravity' else '100000', has_initial_checkpoint='yes',
                      max_steps={'seed': '0', 'gravity': '2'}.get(name, '1000000000'),
                      max_time_in_yr=repr((DURATION if factor else 10*GRAVITY_DT*h)/YEAR2SEC))
    if name != 'seed':
        cfg['sim'].update(is_restarting='yes', restarting_from_modelname='seed' if name == 'gravity' else 'gravity',
                          restarting_from_frame='0' if name == 'gravity' else '2')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='slope.poly', xlength=repr(WIDTH), ylength='1',
                       zlength=repr(-br.BASE), resolution=repr(h), remeshing_option='0', min_angle='28', max_steiner_factor='0')
    cfg['control'].update(gravity='0' if name == 'seed' else repr(G_ACC), fixed_dt=repr(GRAVITY_DT*h), has_thermal_diffusion='false',
                          has_hydraulic_diffusion='false', has_poroelastic_pressure_feedback='false',
                          has_moving_mesh='false', is_using_mixed_stress='false',
                          PT_option='1', PT_max_iter='200000', PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8')
    if factor:
        cfg['sim'].update(output_step_interval='1000000000', output_time_interval_in_yr=repr(DURATION/FRAMES/YEAR2SEC))
        cfg['control'].update(has_PT='false', fixed_dt='0', inertial_scaling=repr(SCALING), characteristic_speed='1e-5',
                              dt_fraction='0.5', damping_factor=repr(DAMPING))
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc.update({f'vbc_{side}': '0', f'vbc_val_{side}': '0', f'stress_bc_{side}': '0', f'hbc_{side}': '0'})
    bc.update(vbc_x0='1', vbc_x1='1', vbc_z0='3')
    cfg['ic'].update(initial_stress_option='1' if name == 'seed' else '0', initial_stress='[0,0,0]', pore_pressure_option='0',
                     has_body_force_adjustment='false', weakzone_option='0', mattype_option='0')
    c, phi = br.COHESION, br.FRICTION
    if factor: c, phi = c/factor, math.atan(math.tan(phi)/factor)
    cfg['mat'].update(rheology_type='elastic' if factor is None else 'elasto-plastic', is_plane_strain='true', num_materials='1',
                      alpha='[0]', rho0=f'[{RHO!r}]', bulk_modulus=f'[{BULK!r}]', shear_modulus=f'[{SHEAR!r}]',
                      cohesion0=f'[{c!r}]', cohesion1=f'[{c!r}]', friction_angle0=f'[{math.degrees(phi)!r}]',
                      friction_angle1=f'[{math.degrees(phi)!r}]', dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e9')
    return cfg


def run(exe, folder, name, cfg, threads=1):
    with (folder/f'{name}.cfg').open('w') as f: cfg.write(f)
    with (folder/f'{name}.log').open('w') as f:
        return subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, env=dict(os.environ, OMP_NUM_THREADS=str(threads)),
                              stdout=f, stderr=subprocess.STDOUT, timeout=80000).returncode


def analyse(folder, name):
    """Largest displacement (saved velocities integrated by the trapezoid rule), crest-edge displacement, peak speed and
    plastic extent per frame."""
    model = Dynearthsol(str(folder/name)); rows = scan_frames(str(folder/name))
    first = rows[0]['frame']
    x0 = model.read_field(first, 'coordinate'); conn = model.read_field(first, 'connectivity')
    edge = int(np.argmin(np.hypot(x0[:, 0] - br.TOE - br.RUN, x0[:, 1])))
    out = []
    u = np.zeros_like(x0); v0 = model.read_field(first, 'velocity'); t0 = rows[0]['time']
    for row in rows[1:]:
        v1 = model.read_field(row['frame'], 'velocity')
        u += (v0 + v1)/2*(row['time'] - t0); v0, t0 = v1, row['time']
        mag = np.linalg.norm(u, axis=1)
        pls = model.read_field(row['frame'], 'plastic strain')
        out.append(dict(time=row['time'], max_displacement=float(mag.max()), crest_displacement=float(mag[edge]),
                        max_speed=float(np.linalg.norm(model.read_field(row['frame'], 'velocity'), axis=1).max()),
                        plastic_fraction=float(np.mean(pls > 0)), max_plastic_strain=float(pls.max())))
    return out, len(x0), len(conn)


def deceleration(frames):
    """Mean peak speed of the last five frames over that of five frames three quarters through; below 1 the slope is coming to rest."""
    v = [f['max_speed'] for f in frames]
    q = 3*len(v)//4
    return float(np.mean(v[-5:])/np.mean(v[q - 2:q + 3]))


def at_rest(frames):
    return deceleration(frames) < 1


def study(exe, out, spacings, factors, workers):
    out.mkdir(parents=True, exist_ok=False)
    fs_bishop, circle = br.embankment().critical(br.starts())
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  reference_sha256=digest(br.__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  height=br.H, run=br.RUN, toe=br.TOE, crest=br.CREST, base=br.BASE, gamma=br.GAMMA, cohesion=br.COHESION,
                  friction_deg=math.degrees(br.FRICTION), young=E, poisson=NU, duration=DURATION, frames=FRAMES,
                  scaling=SCALING, damping=DAMPING, fs_bishop=fs_bishop, bishop_circle=list(br.embankment().circle(*circle)), cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        polygon(folder/'slope.poly', h)
        for name in ('seed', 'gravity'):
            code = run(exe, folder, name, config(h, name))
            if code: raise RuntimeError(f'{folder}/{name}: exit {code}')
        names = [f'srf{f:.3f}' for f in factors]
        with ThreadPoolExecutor(workers) as pool:
            codes = list(pool.map(lambda nf: run(exe, folder, nf[0], config(h, nf[0], nf[1])), zip(names, factors)))
        case = dict(spacing_m=h, trials=[])
        for name, f, code in zip(names, factors, codes):
            if code: raise RuntimeError(f'{folder}/{name}: exit {code}')
            frames, case['nodes'], case['elements'] = analyse(folder, name)
            case['trials'].append(dict(factor=f, at_rest=at_rest(frames), deceleration=deceleration(frames), frames=frames))
            print(json.dumps(dict(spacing_m=h, factor=f, at_rest=at_rest(frames), max_displacement=frames[-1]['max_displacement'],
                                  final_speed=frames[-1]['max_speed'])), flush=True)
        result['cases'].append(case)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[1.0, 0.5])
    p.add_argument('--factors', type=float, nargs='+', default=FACTORS)
    p.add_argument('--workers', type=int, default=4)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.factors, a.workers)
