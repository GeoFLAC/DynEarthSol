#!/usr/bin/env python3
"""A soil layer on inclined bedrock slides as the bedrock pore pressure rises.

A layer of vertical thickness D on bedrock inclined at beta, with strong ends.
The base is fixed and supplies water, the surface drains, and gravity is
applied as a load by restarting an unstressed seed. The base pressure then
rises slowly. Failure is read at a fixed basal plastic strain and compared
with the infinite-slope factor of safety and Davis's non-associated estimate.

    slope_application.py EXECUTABLE OUTPUT [--spacings 0.5 0.25] [--steps 404] [--cohesion PA] [--friction DEG]
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

BETA, D, L = math.radians(20.0), 2.0, 40.0
G_ACC, RHO_S, RHO_F, PHI_N = 9.81, 2650.0, 1000.0, .3
GAMMA = (RHO_S*(1 - PHI_N) + RHO_F*PHI_N)*G_ACC
COHESION, FRICTION = 2e3, math.radians(35.0)
MAT = dict(bulk=5e7, shear=5.17e6, biot=1.0, fluid_bulk=2.2e9, permeability=1e-11, viscosity=1e-3)
SETTLE, RAMP = 20.0, 2000.0              # s: gravity settles, then the base pressure ramps
SIGMA_N, TAU = GAMMA*D*math.cos(BETA)**2, GAMMA*D*math.sin(BETA)*math.cos(BETA)
U_STAR = SIGMA_N - (TAU - COHESION)/math.tan(FRICTION)     # base pore pressure at FS = 1
# Pointwise Mohr-Coulomb yield on the poroelastic path of the infinite slope. Gravity turn-on gives
# slope-parallel effective stress K0 sigma_n (zero strain along the slope); a pore-pressure rise u
# lowers the slope-normal effective stress by u and the slope-parallel one by K0 u, so
# s = sigma_n - u scales both while the shear stays tau.
NU = (3*MAT['bulk'] - 2*MAT['shear'])/(2*(3*MAT['bulk'] + MAT['shear']))
K0 = NU/(1 - NU)


def _pointwise_yield():
    from scipy.optimize import brentq
    f = lambda s: (math.hypot((1 - K0)/2*s, TAU) - (1 + K0)/2*s*math.sin(FRICTION) - COHESION*math.cos(FRICTION))
    return SIGMA_N - brentq(f, 1.0, SIGMA_N)


U_POINT = _pointwise_yield()


def estimates(cohesion, friction):
    """Pointwise onset, Davis non-associated and classical base pore pressures for another strength."""
    from scipy.optimize import brentq
    f = lambda s: (math.hypot((1 - K0)/2*s, TAU) - (1 + K0)/2*s*math.sin(friction) - cohesion*math.cos(friction))
    point = SIGMA_N - brentq(f, 1.0, SIGMA_N)
    davis = SIGMA_N - (TAU - cohesion*math.cos(friction))/math.sin(friction)
    classical = SIGMA_N - (TAU - cohesion)/math.tan(friction)
    return point, davis, classical
# Davis (1968): a non-associated shear band fails as if c* = c cos(psi) cos(phi)/(1 - sin(psi) sin(phi)),
# tan(phi*) = cos(psi) sin(phi)/(1 - sin(psi) sin(phi)); here psi = 0.
PHI_DAVIS, C_DAVIS = math.atan(math.sin(FRICTION)), COHESION*math.cos(FRICTION)
U_DAVIS = SIGMA_N - (TAU - C_DAVIS)/math.tan(PHI_DAVIS)
P_MAX = 1.5*U_STAR


def factor_of_safety(u):
    return (COHESION + (SIGMA_N - u)*math.tan(FRICTION))/TAU


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


END = 4.0                                 # m of strong material at each end (no corner yielding)


def polygon(path, h):
    """Parallelogram: surface z = -x tan(beta), base D below it, vertical ends; counter-clockwise.
    Vertical constraint edges at x = END and L - END split it into strong end regions (attribute 1)
    and the Mohr-Coulomb middle (attribute 0)."""
    t = math.tan(BETA); m = math.ceil(D/h)
    pieces = [(0.0, END), (END, L - END), (L - END, L)]
    xs = np.unique(np.concatenate([np.linspace(a, b, math.ceil((b - a)/math.cos(BETA)/h) + 1) for a, b in pieces]))
    n = len(xs) - 1
    base = [(x, -x*t - D) for x in xs]
    right = [(L, -L*t - D + D*k/m) for k in range(1, m)]
    surface = [(x, -x*t) for x in xs[::-1]]
    left = [(0.0, -D*k/m) for k in range(1, m)]
    points = base + right + surface + left
    markers = [16]*n + [2]*m + [32]*n + [1]*m          # z0 base, x1, z1 surface, x0
    edges = [(i, (i + 1) % len(points), mk) for i, mk in enumerate(markers)]
    for xe in (END, L - END):                           # interior vertical constraints, unmarked
        ib = int(np.argmin(np.abs(xs - xe))); it = len(base) + len(right) + (n - ib)
        ids = [ib]
        for k in range(1, m):
            ids.append(len(points)); points.append((xe, -xe*t - D + D*k/m))
        ids.append(it)
        edges += [(a, b, 0) for a, b in zip(ids[:-1], ids[1:])]
    with path.open('w') as f:
        f.write(f'{len(points)} 2 0 0\n')
        for i, (x, z) in enumerate(points): f.write(f'{i} {x:.17g} {z:.17g}\n')
        f.write(f'{len(edges)} 1\n')
        for i, (a, b, mk) in enumerate(edges): f.write(f'{i} {a} {b} {mk}\n')
        regions = [(END/2, 1), (L/2, 0), (L - END/2, 1)]
        f.write(f'0\n{len(regions)}\n')
        for i, (x, attr) in enumerate(regions):
            f.write(f'{i} {x} {-x*t - D/2} {attr} {h*h/2}\n')


def config(h, steps, loaded, dilation=0.0, cohesion=COHESION, friction=FRICTION):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    duration = SETTLE + RAMP
    cfg['sim'].update(modelname='slope' if loaded else 'seed', max_steps=str(steps if loaded else 0),
                      max_time_in_yr=repr(2*duration/YEAR2SEC), output_step_interval='1', output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='yes')
    if loaded:
        cfg['sim'].update(is_restarting='yes', restarting_from_modelname='seed', restarting_from_frame='0')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='slope.poly', xlength=repr(L), ylength='1',
                       zlength=repr(L*math.tan(BETA) + D), resolution=repr(h), remeshing_option='0',
                       min_angle='28', max_steiner_factor='0')
    cfg['control'].update(gravity=repr(G_ACC) if loaded else '0', fixed_dt=repr(duration/steps), has_thermal_diffusion='false',
                          has_moving_mesh='false', is_using_mixed_stress='false', has_poroelastic_pressure_feedback='true',
                          PT_option='1', PT_max_iter='200000', PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc.update({f'vbc_{side}': '0', f'vbc_val_{side}': '0', f'stress_bc_{side}': '0', f'hbc_{side}': '0'})
    bc.update(vbc_x0='1', vbc_x1='1', vbc_z0='3', hbc_z0='2', hbc_z1='2',
              hbc_pressure=f'[0,0,0,0,{P_MAX!r},0]', num_hbc_schedule='3',
              hbc_schedule_time_in_yr=f'[0, {SETTLE/YEAR2SEC!r}, {(SETTLE + RAMP)/YEAR2SEC!r}]',
              hbc_schedule_ratio='[0, 0, 1]')
    cfg['ic'].update(initial_stress_option='0' if loaded else '1', initial_stress='[0,0,0]', pore_pressure_option='0',
                     has_body_force_adjustment='false', weakzone_option='0')
    m = MAT
    cfg['mat'].update(rheology_type='elasto-plastic', is_plane_strain='true', alpha='[0]', rho0=f'[{RHO_S!r}]',
                      fluid_rho0=f'[{RHO_F!r}]', porosity=f'[{PHI_N!r}]', bulk_modulus=f"[{m['bulk']!r}]",
                      shear_modulus=f"[{m['shear']!r}]", biot_coeff=f"[{m['biot']!r}]", fluid_bulk_modulus=f"[{m['fluid_bulk']!r}]",
                      bulk_modulus_s='[1e20]', hydraulic_perm=f"[{m['permeability']!r}]", fluid_visc=f"[{m['viscosity']!r}]",
                      fluid_alpha='[0]', retention_model='[0]',
                      cohesion0=f'[{cohesion!r}]', cohesion1=f'[{cohesion!r}]',
                      friction_angle0=f'[{math.degrees(friction)!r}]', friction_angle1=f'[{math.degrees(friction)!r}]',
                      dilation_angle0=f'[{dilation!r}]', dilation_angle1=f'[{dilation!r}]', max_tension='1e9')
    # Material 1 (both ends): the same rock, strong enough to stay elastic, so the corners cannot yield.
    cfg['mat']['num_materials'] = '2'
    cfg['ic']['mattype_option'] = '0'
    for key, value in list(cfg['mat'].items()):
        if value.startswith('[') and value.endswith(']') and ',' not in value:
            cfg['mat'][key] = f'[{value[1:-1]}, {value[1:-1]}]'
    cfg['mat'].update(cohesion0=f'[{cohesion!r}, 1e7]', cohesion1=f'[{cohesion!r}, 1e7]')
    return cfg


def analyse(folder):
    """Base pressure, mid-slope plastic strain at the base and surface slip, per saved step."""
    model = Dynearthsol(str(folder/'slope')); rows = scan_frames(str(folder/'slope'))
    first = rows[0]['frame']                       # a restart from the gravity checkpoint has no frame 0
    x = model.read_field(first, 'coordinate'); conn = model.read_field(first, 'connectivity')
    centre = x[conn].mean(axis=1)
    depth = -centre[:, 1] - centre[:, 0]*math.tan(BETA)          # vertical depth below the surface
    mid = (centre[:, 0] > L/3) & (centre[:, 0] < 2*L/3)
    basal = mid & (depth > D - 0.6)
    surface_nodes = (np.abs(x[:, 1] + x[:, 0]*math.tan(BETA)) < 1e-9) & (x[:, 0] > L/3) & (x[:, 0] < 2*L/3)
    base_nodes = (np.abs(x[:, 1] + x[:, 0]*math.tan(BETA) + D) < 1e-9) & (x[:, 0] > L/3) & (x[:, 0] < 2*L/3)
    down = np.array([math.cos(BETA), -math.sin(BETA)])
    out = []
    for row in rows:
        if row['time'] <= 0: continue
        p = model.read_field(row['frame'], 'pore pressure')
        pls = model.read_field(row['frame'], 'plastic strain')
        v = model.read_field(row['frame'], 'velocity')
        out.append(dict(time=row['time'], base_pressure=float(np.mean(p[base_nodes])),
                        basal_plastic_fraction=float(np.mean(pls[basal] > 0)),
                        basal_plastic_strain=float(np.mean(pls[basal])),
                        surface_velocity=float(np.mean(v[surface_nodes] @ down))))
    return out, len(x), len(conn)


def study(exe, out, spacings, steps, dilation=0.0, cohesion=COHESION, friction=FRICTION, elastic_gravity=False):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  beta_deg=math.degrees(BETA), thickness=D, length=L, cohesion=COHESION, friction_deg=math.degrees(FRICTION),
                  gamma=GAMMA, sigma_n=SIGMA_N, tau=TAU, u_star=U_STAR, u_point=estimates(cohesion, friction)[0], u_davis=estimates(cohesion, friction)[1], k0=K0, p_max=P_MAX, steps=steps,
                  dilation_deg=dilation, cohesion_Pa=cohesion, friction_deg_used=math.degrees(friction),
                  elastic_gravity=elastic_gravity, cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        polygon(folder/'slope.poly', h)
        stages = [('seed', config(h, steps, False, dilation, cohesion, friction))]
        if elastic_gravity:
            # Gravity first with elastic material (one step, checkpointed), then the plastic run restarts from it.
            grav = config(h, steps, True, dilation, cohesion, friction)
            grav['sim'].update(modelname='gravity', max_steps='1', checkpoint_frame_interval='1')
            grav['mat']['rheology_type'] = 'elastic'
            stages.append(('gravity', grav))
        load = config(h, steps, True, dilation, cohesion, friction)
        if elastic_gravity:
            load['sim'].update(restarting_from_modelname='gravity', restarting_from_frame='1')
        stages.append(('slope', load))
        for name, cfg in stages:
            with (folder/f'{name}.cfg').open('w') as f: cfg.write(f)
            with (folder/f'{name}.log').open('w') as f:
                code = subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, env=env, stdout=f,
                                      stderr=subprocess.STDOUT, timeout=40000).returncode
            if code not in (0, 54) or (code == 54 and name != 'slope'):
                raise RuntimeError(f'{folder}/{name}: exit {code}')
        frames, nodes, elements = analyse(folder)
        item = dict(spacing_m=h, nodes=nodes, elements=elements, lost_equilibrium=code == 54, frames=frames)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        last = frames[-1]
        print(json.dumps(dict(spacing_m=h, lost_equilibrium=code == 54, last_base_pressure=last['base_pressure'],
                              u_star=U_STAR, last_basal_plastic_fraction=last['basal_plastic_fraction'])), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.5, 0.25])
    p.add_argument('--steps', type=int, default=404)
    p.add_argument('--dilation', type=float, default=0.0, help='dilation angle (deg); equal to friction for associated flow')
    p.add_argument('--cohesion', type=float, default=COHESION, help='Pa')
    p.add_argument('--friction', type=float, default=math.degrees(FRICTION), help='deg')
    p.add_argument('--elastic-gravity', action='store_true', help='apply gravity with elastic material first')
    a = p.parse_args()
    print(f'infinite slope: sigma_n {SIGMA_N:.0f} Pa, tau {TAU:.0f} Pa, FS(0) {factor_of_safety(0):.3f}, u* {U_STAR:.0f} Pa, '
          f'pointwise yield u {U_POINT:.0f} Pa (K0 {K0:.3f}), Davis non-associated u {U_DAVIS:.0f} Pa')
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps, a.dilation, a.cohesion, math.radians(a.friction),
          a.elastic_gravity)
