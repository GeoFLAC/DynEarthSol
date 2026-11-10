#!/usr/bin/env python3
"""Check Boussinesq buoyancy against Horton-Rogers-Lapwood linear stability.

A unit square porous box is heated from below: closed hydraulic boundaries,
fixed top/bottom temperatures, insulated sides. The conductive profile plus
a small Gaussian seed starts each run (ic.temperature_option 5). The growth
rate of the cos(pi*x)*sin(pi*z) mode is fitted and compared with
benchmarks/poromechanics/convection_reference.py. A seed-free run measures
the discrete rest state. Storage is chosen so that pressure diffusion is
~1e3 times faster than heat diffusion, the incompressible limit of the
reference. Mechanics is stiff and does not feed back on the pressure.
"""
import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'benchmarks'/'poromechanics'))
from Dynearthsol import Dynearthsol, scan_frames
from convection_reference import rayleigh, growth_rate, critical_rayleigh

RHO, CF, RHO_S, CS, LAMBDA, G, ALPHA, DT, MU = 1000., 1000., 1000., 1000., 1., 10., 1e-3, 10., 1e-3
T_TOP = 273.
CAPACITY = RHO_S*CS  # equals RHO*CF, so independent of porosity


from des_case import digest, read_config, reset_boundaries, run_solver, source_diff, source_revision


def config(h, ratio, seed, steps, dt, outputs):
    k = ratio*critical_rayleigh(1, 1)/rayleigh(RHO, G, ALPHA, DT, 1, 1, CF, MU, LAMBDA)
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='box', max_steps=str(steps), output_step_interval=str(steps//outputs),
                      output_time_interval_in_yr='1e9', max_time_in_yr='1e9', checkpoint_frame_interval='10')
    cfg['mesh'].update(xlength='1', ylength='1', zlength='1', resolution=str(h),
                       largest_size=str(h), smallest_size=str(h))
    cfg['control'].update(has_PT='false', is_quasi_static='true', inertial_scaling='1e6',
        use_global_velocity_scaling='false', has_moving_mesh='false', gravity=str(G), damping_option='0',
        fixed_dt=repr(dt), has_thermal_diffusion='true', has_fluid_thermal_advection='true',
        has_poroelastic_pressure_feedback='false')
    cfg['ic'].update(has_body_force_adjustment='false', temperature_option='5',
        temp_anomaly_amplitude=repr(seed), temp_anomaly_xcenter='0.25', temp_anomaly_zcenter='-0.5',
        temp_anomaly_xsemi_axis='0.25', temp_anomaly_ysemi_axis='0.25', temp_anomaly_zsemi_axis='0.25')
    reset_boundaries(cfg['bc'], vbc='1', hbc='0', tbc='0')
    cfg['bc'].update(surface_temperature=repr(T_TOP), mantle_temperature=repr(T_TOP+DT),
                     tbc_z0='1', tbc_val_z0=repr(T_TOP+DT), tbc_z1='1', tbc_val_z1=repr(T_TOP),
                     stress_bc_z1='0', stress_val_z1='0')
    cfg['mat'].update(rho0=f'[{RHO_S}]', alpha='[0]', heat_capacity=f'[{CS}]', therm_cond=f'[{LAMBDA}]',
        fluid_rho0=f'[{RHO}]', fluid_heat_capacity=f'[{CF}]', fluid_alpha=f'[{ALPHA}]', fluid_visc=f'[{MU}]',
        hydraulic_perm=f'[{k!r}]', porosity='[0.5]', biot_coeff='[1]', fluid_bulk_modulus='[1e3]',
        bulk_modulus='[1e9]', shear_modulus='[1e9]', bulk_modulus_s='[1e12]', retention_model='[0]')
    return cfg


def lumped_area(coord, conn):
    w = np.zeros(len(coord))
    for tri in conn:
        (x0, z0), (x1, z1), (x2, z2) = coord[tri]
        np.add.at(w, tri, abs((x1-x0)*(z2-z0) - (x2-x0)*(z1-z0))/6)
    return w


def amplitudes(folder):
    model = Dynearthsol(str(folder/'box')); rows = scan_frames(str(folder/'box'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    assert coord.shape[1] == 2, 'the projection is two-dimensional'
    x, z = coord[:, 0], coord[:, 1]
    w = lumped_area(coord, conn)*4*np.cos(math.pi*x)*np.sin(math.pi*(z+1))
    out = []
    for row in rows:
        theta = model.read_field(row['frame'], 'temperature') - (T_TOP - DT*z)
        assert np.isfinite(theta).all()
        out.append((row['time'], float(w @ theta), float(np.abs(theta).max())))
    return out


def study(exe, out, h, ratios, steps, dt, seed):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(Path(__file__)),
                  source_commit=source_revision(ROOT),
                  spacing=h, steps=steps, dt=dt, seed_K=seed, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for ratio, amp in [(r, seed) for r in ratios] + [(ratios[-1], 0.0)]:
        folder = out/f'ra{ratio:g}-seed{amp:g}'; folder.mkdir()
        cfg = config(h, ratio, amp, steps, dt, 40)
        with (folder/'box.cfg').open('w') as f: cfg.write(f)
        run_solver([str(exe), 'box.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=3600)
        series = amplitudes(folder)
        item = dict(ra_over_critical=ratio, seed_K=amp, series=series)
        if amp:
            t, a = np.array([s[0] for s in series]), np.array([s[1] for s in series])
            late = slice(len(t)//2, None)
            fit = float(np.polyfit(t[late], np.log(np.abs(a[late])), 1)[0])
            exact = growth_rate(ratio*critical_rayleigh(1, 1), LAMBDA, CAPACITY, 1, 1)
            item.update(fitted_rate=fit, exact_rate=exact, relative_error=abs(fit-exact)/abs(exact),
                        sign_ok=bool(np.sign(fit) == np.sign(exact)))
        else:
            item.update(rest_max_theta_over_dT=series[-1][2]/DT)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2))
        print(json.dumps({k: v for k, v in item.items() if k != 'series'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacing', type=float, default=0.05)
    p.add_argument('--ratios', type=float, nargs='+', default=[0.5, 1.5, 2.0])
    p.add_argument('--steps', type=int, default=400)
    p.add_argument('--dt', type=float, default=1000.)
    p.add_argument('--seed', type=float, default=1e-3)
    p.add_argument('--tolerance', type=float, default=None,
                   help='fail if any fitted growth rate differs from the reference by more than this fraction')
    a = p.parse_args()
    r = study(a.executable.resolve(), a.output.resolve(), a.spacing, a.ratios, a.steps, a.dt, a.seed)
    seeded = [c for c in r['cases'] if c['seed_K']]
    assert all(c['sign_ok'] for c in seeded), 'growth/decay disagrees with the onset criterion'
    if a.tolerance is not None:
        assert all(c['relative_error'] <= a.tolerance for c in seeded), 'growth rate outside tolerance'
