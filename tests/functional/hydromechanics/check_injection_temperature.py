#!/usr/bin/env python3
"""DES hot-water injection from a line source: injection temperature (Avdonin).

A 2D quarter plane 0 <= x, -z <= L receives fluid at T_in from a point
source at its corner (a line in the plane-strain sense), Q/4 of the total
rate q per unit thickness; the symmetry sides are sealed and insulated and
the far sides drained (p = 0) with natural advective outflow. The layer
starts at T0, gravity is off and the skeleton is stiff without pressure
feedback. Storage is chosen so the flow is quasi-steady over the thermal
time while the explicit pressure substeps stay affordable. The injected
fluid carries injection.temperature: an injecting node gains
rho_f c_f Q (T_in - T).

Checks: nodal temperature against benchmarks/poromechanics/avdonin_reference.py
at four times, relative to T_in - T0; and the stored heat
sum C A_n (T_n - T0) against the injected rho_f c_f (Q/4)(T_in - T0) t,
which holds while the warm region stays inside the domain.

In 3D (--ndims 3) the quarter plane is a slab one element thick with
insulated, sealed faces and the source at mid-thickness on its corner edge,
which acts as a line source of q/4 per unit thickness: a sanity check.

    check_injection_temperature.py EXECUTABLE OUTPUT [--spacings 0.2 0.1 0.05] [--steps 500] [--ndims 2]
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
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics')]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from avdonin_reference import Avdonin

L, Q, DURATION = 8.0, 2e-5, 5e5                     # quarter side (m), total rate per thickness (m2/s), end (s)
T0, TIN = 283.0, 333.0
MAT = dict(rho_solid=2650.0, cp_solid=900.0, conductivity=2.0, porosity=.25, rho_fluid=1000.0, cp_fluid=4180.0,
           permeability=1e-13, viscosity=1e-3, fluid_bulk=2.5e7, bulk=1e12, shear=1e12)
CAPACITY = (1 - MAT['porosity'])*MAT['rho_solid']*MAT['cp_solid'] + MAT['porosity']*MAT['rho_fluid']*MAT['cp_fluid']
FLUID = MAT['rho_fluid']*MAT['cp_fluid']
REF = Avdonin(rate=Q, conductivity=MAT['conductivity'], capacity=CAPACITY, fluid_capacity=FLUID, t0=T0, t_in=TIN)
TIMES = (1.25e5, 2.5e5, 3.75e5, 5e5)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config(h, steps, ndims=2):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='hot', max_steps=str(steps), max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval=str(steps//4), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(L), ylength=repr(h if ndims == 3 else L), zlength=repr(L), resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(gravity='0', has_PT='false', is_quasi_static='true', inertial_scaling='1e6',
                          use_global_velocity_scaling='false', has_moving_mesh='false', damping_option='0',
                          fixed_dt=repr(DURATION/steps), has_poroelastic_pressure_feedback='false',
                          has_thermal_diffusion='true', has_fluid_thermal_advection='true')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'], bc[f'tbc_{side}'] = '1', '0', '0', '0'
    bc.update(hbc_x1='2', hbc_z0='2', hbc_pressure='[0,0,0,0,0,0]', stress_bc_z1='0',
              surface_temperature=repr(T0), mantle_temperature=repr(T0))
    cfg['ic'].update(temperature_option='4', temp_background=repr(T0), temp_anomaly_amplitude='0',
                     pore_pressure_option='0', has_body_force_adjustment='false', weakzone_option='0')
    cfg['injection'] = dict(enabled='true', num_points='1', points_x='[0]', points_z='[0]', points_unit='m',
                            rate_model='constant_rate', rate=f'[{Q/4*(h if ndims == 3 else 1)!r}]', start_time_in_yr='[0]',
                            temperature=f'[{TIN!r}]')
    if ndims == 3: cfg['injection']['points_y'] = f'[{h/2!r}]'
    m = MAT
    cfg['mat'].update(rho0=f"[{m['rho_solid']!r}]", heat_capacity=f"[{m['cp_solid']!r}]",
                      therm_cond=f"[{m['conductivity']!r}]", alpha='[0]', porosity=f"[{m['porosity']!r}]",
                      fluid_rho0=f"[{m['rho_fluid']!r}]", fluid_heat_capacity=f"[{m['cp_fluid']!r}]", fluid_alpha='[0]',
                      hydraulic_perm=f"[{m['permeability']!r}]", fluid_visc=f"[{m['viscosity']!r}]",
                      fluid_bulk_modulus=f"[{m['fluid_bulk']!r}]", bulk_modulus=f"[{m['bulk']!r}]",
                      shear_modulus=f"[{m['shear']!r}]", biot_coeff='[1]', bulk_modulus_s='[1e20]', retention_model='[0]')
    return cfg


def lumped_area(coord, conn):
    """Lumped nodal area (2D) or volume (3D)."""
    c = coord[conn]
    edges = c[:, 1:] - c[:, :1]
    size = np.abs(np.linalg.det(edges))/(2 if coord.shape[1] == 2 else 6)
    w = np.zeros(len(coord))
    np.add.at(w, conn, size[:, None]/conn.shape[1])
    return w


def compare(folder):
    model = Dynearthsol(str(folder/'hot')); rows = scan_frames(str(folder/'hot'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    r = np.hypot(coord[:, 0], coord[:, -1]); area = lumped_area(coord, conn)
    thickness = np.ptp(coord[:, 1]) if coord.shape[1] == 3 else 1.0
    out = []
    for target in TIMES:
        row = min(rows[1:], key=lambda x: abs(x['time'] - target))
        assert abs(row['time'] - target) < 1e-6*DURATION
        T = model.read_field(row['frame'], 'temperature')
        band = r > 0.25                                     # the singular source neighbourhood is excluded
        err = np.abs(T[band] - REF.temperature(r[band], target))/(TIN - T0)
        stored = CAPACITY*float(area @ (T - T0))
        injected = FLUID*Q/4*thickness*(TIN - T0)*target
        out.append(dict(time=target, temperature_error=float(err.max()), temperature_rms=float(np.sqrt(np.mean(err**2))),
                        stored_heat=stored, injected_heat=injected, heat_error=abs(stored/injected - 1),
                        far_side=float(np.max(T[(coord[:, 0] > L - 1e-9) | (coord[:, -1] < -L + 1e-9)]) - T0)))
    return out, len(coord)


def study(exe, out, spacings, steps, ndims=2):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  steps=steps, ndims=ndims, material=MAT, capacity=CAPACITY, nu=REF.nu, front_m=REF.front(DURATION), cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'hot.cfg').open('w') as f: config(h, steps, ndims).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'hot.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=40000)
        frames, nodes = compare(folder)
        item = dict(spacing_m=h, nodes=nodes, frames=frames, input_sha256=digest(folder/'hot.cfg'),
                    temperature_error=max(f['temperature_error'] for f in frames),
                    heat_error=max(f['heat_error'] for f in frames))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1, 0.05])
    p.add_argument('--steps', type=int, default=500)
    p.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps, a.ndims)
