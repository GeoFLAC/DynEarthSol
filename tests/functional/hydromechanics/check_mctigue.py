#!/usr/bin/env python3
"""DES McTigue column: thermal pressurization under a sudden surface heating.

A laterally confined poroelastic column -L <= z <= 0 at uniform temperature
T0 has its top raised to T0 + dT at t = 0; the top is drained (p = 0) and
traction-free, the base and sides are sealed, insulated and fixed normally.
Heat conducts without advection; gravity is off; PT keeps equilibrium and
the thermal pressurization source drives the pressure. While both diffusion
lengths stay well inside the column it behaves as McTigue's half-space:
benchmarks/poromechanics/mctigue_reference.py.

DES applies a fixed boundary temperature at initialization, so a surface
that starts hot carries no thermal stress or pressurization at its nodes.
To apply the step at t = 0+, as in the reference, a zero-step run at T0
writes the initial checkpoint and the column restarts from it with the
surface at T0 + dT.

Temperature and pressure are compared at the nodes, vertical strain at the
element layers and the surface heave from the integrated layer strain.

    check_mctigue.py EXECUTABLE OUTPUT [--ndims 2] [--spacings 0.1 0.05 0.025] [--steps 200]
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
from mctigue_reference import McTigue

L, T0, DT, DURATION = 1.0, 273.0, 10.0, 4e4     # column (m), start and step temperature (K), end time (s)
MAT = dict(bulk=1e9, shear=1e9, biot=.8, porosity=.2, fluid_bulk=2.2e9, grain_bulk=5e9, permeability=1e-19,
           viscosity=1e-3, alpha_solid=3e-5, alpha_fluid=2.07e-4, rho=2650.0, cp=1000.0, conductivity=2.65)
KAPPA = MAT['conductivity']/(MAT['rho']*MAT['cp'])
REF = McTigue(MAT['bulk'], MAT['shear'], MAT['biot'], MAT['porosity'], MAT['fluid_bulk'], MAT['grain_bulk'],
              MAT['permeability'], MAT['viscosity'], MAT['alpha_solid'], MAT['alpha_fluid'], KAPPA, DT)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config(ndims, h, steps, start=None):
    """The heated run; with start=None, the zero-step run at T0 that it restarts from."""
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='mctigue' if start else 'start', max_steps=str(steps if start else 0),
                      max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval=str(max(1, steps//8)), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='yes')
    if start:
        cfg['sim'].update(is_restarting='yes', restarting_from_modelname=str(start), restarting_from_frame='0')
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(h), ylength=repr(h), zlength=repr(L), resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(gravity='0', fixed_dt=repr(DURATION/steps), PT_max_iter='100000',
                          has_thermal_diffusion='true', has_thermal_pressurization='true',
                          has_fluid_thermal_advection='false', has_poroelastic_pressure_feedback='true')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'], bc[f'tbc_{side}'] = '1', '0', '0', '0'
    bc.update(hbc_z1='2', hbc_pressure='[0,0,0,0,0,0]', tbc_z1='1', tbc_val_z1=repr(T0 + DT if start else T0),
              surface_temperature=repr(T0), mantle_temperature=repr(T0))
    cfg['ic'].update(temperature_option='4', temp_background=repr(T0), temp_anomaly_amplitude='0',
                     pore_pressure_option='0', has_body_force_adjustment='false')
    cfg['mat'].update(rho0=f"[{MAT['rho']!r}]", heat_capacity=f"[{MAT['cp']!r}]", therm_cond=f"[{MAT['conductivity']!r}]",
                      alpha=f"[{MAT['alpha_solid']!r}]", fluid_alpha=f"[{MAT['alpha_fluid']!r}]",
                      bulk_modulus=f"[{MAT['bulk']!r}]", shear_modulus=f"[{MAT['shear']!r}]",
                      biot_coeff=f"[{MAT['biot']!r}]", porosity=f"[{MAT['porosity']!r}]",
                      fluid_bulk_modulus=f"[{MAT['fluid_bulk']!r}]", bulk_modulus_s=f"[{MAT['grain_bulk']!r}]",
                      hydraulic_perm=f"[{MAT['permeability']!r}]", fluid_visc=f"[{MAT['viscosity']!r}]",
                      retention_model='[0]')
    return cfg


def layer_strain(coord, conn, ezz, h):
    """Mean vertical strain of each element layer, from the top down, and layer mid-depths."""
    depth = -coord[conn][:, :, -1].mean(axis=1)
    layer = np.floor(depth/h).astype(int)
    n = int(round(L/h))
    count = np.bincount(layer, minlength=n).astype(float)
    assert len(count) == n and count.min() > 0, 'the column mesh must be structured in layers of height h'
    return np.bincount(layer, weights=ezz, minlength=n)/count, (np.arange(n) + .5)*h


def compare(folder, h):
    base = Dynearthsol(str(folder/'start'))
    model = Dynearthsol(str(folder/'mctigue')); rows = scan_frames(str(folder/'mctigue'))
    coord = base.read_field(0, 'coordinate'); conn = base.read_field(0, 'connectivity')
    assert np.all(base.read_field(0, 'temperature') == T0), 'the column must start uniformly at T0'
    depth = -coord[:, -1]; zz = coord.shape[1] - 1
    e0, _ = layer_strain(coord, conn, base.read_field(0, 'strain')[:, zz], h)
    frames = []
    for row in [r for r in rows if r['time'] > 0]:
        t = row['time']
        T = model.read_field(row['frame'], 'temperature') - T0
        p = model.read_field(row['frame'], 'pore pressure')
        e, mid = layer_strain(coord, conn, model.read_field(row['frame'], 'strain')[:, zz], h)
        e = e - e0
        pref, Tref, eref = REF.pressure(depth, t), REF.temperature(depth, t), REF.strain(mid, t)
        heave = float(np.sum(e)*h)
        frames.append(dict(time=t, temperature_error=float(np.max(np.abs(T - Tref)))/DT,
                           pressure_error=float(np.max(np.abs(p - pref))/np.abs(pref).max()),
                           strain_error=float(np.max(np.abs(e - eref))/np.abs(eref).max()),
                           heave_m=heave, heave_ref_m=REF.heave(t), heave_error=abs(heave/REF.heave(t) - 1),
                           peak_pressure_Pa=float(p.max()), peak_pressure_ref_Pa=float(pref.max())))
    return frames, len(coord)


def study(exe, out, ndims, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  ndims=ndims, steps=steps, duration_s=DURATION, material=MAT,
                  reference=dict(kappa=KAPPA, c=REF.c, gamma=REF.gamma, S=REF.S, dT=DT), cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'start.cfg').open('w') as f: config(ndims, h, steps).write(f)
        with (folder/'mctigue.cfg').open('w') as f: config(ndims, h, steps, folder/'start').write(f)
        for name in ('start', 'mctigue'):
            with (folder/f'{name}.log').open('w') as f:
                subprocess.run([str(exe), f'{name}.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                               check=True, timeout=20000)
        frames, nodes = compare(folder, h)
        item = dict(spacing_m=h, nodes=nodes, frames=frames, input_sha256=digest(folder/'mctigue.cfg'),
                    **{k: max(f[k] for f in frames) for k in
                       ('temperature_error', 'pressure_error', 'strain_error', 'heave_error')})
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, default=2)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.1, 0.05, 0.025])
    p.add_argument('--steps', type=int, default=200)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.ndims, a.spacings, a.steps)
