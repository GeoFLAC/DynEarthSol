#!/usr/bin/env python3
"""Verify explicit unsaturated transport in DES (verification levels U6, U7).

Usage: check_unsaturated_column.py EXECUTABLE OUTPUT_DIRECTORY

- transient: Gardner column against Srivastava and Yeh (1991). A water table at
  the base (their initial state for zero flux) receives a top infiltration
  flux. Checks: the error falls with resolution, the per-step fluid balance
  closes, and pressure goes negative (the unsaturated path runs).
- steady: van Genuchten and Brooks-Corey columns reach the exact steady
  profile for a top flux. Storage does not enter the steady state, so a
  compressible pore fluid shortens the explicit run.
- modes: the 5 cm transient without PT (dynamic relaxation, inertial_scaling
  1e5, with and without global velocity scaling) matches the PT result. The
  transport does not depend on how mechanics reaches equilibrium.
- restart: exact checkpoint continuation of a partly unsaturated column, on a
  fixed and on a remeshed mesh, with the pt_initial_checkpoint.py comparator.
"""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT / 'benchmarks' / 'poromechanics')]
from unsaturated_reference import Retention, srivastava_yeh, steady_pressure

G, RHO, MU = 9.81, 1000.0, 1.002e-3
GAMMA = RHO * G
YEAR = 365.2422 * 86400
KS, ALPHA_HEAD = 1 / 3.6e5, 10.0  # Srivastava and Yeh example: 1 cm/h, 0.1/cm
PERM = KS * MU / GAMMA


from des_case import read_config, restart_check, retention_options, run_model, write_json

def column(ndims, res, tend, dt, *, retention, table=-1.0, flux=0.9 * KS, fluid_bulk=2.17e9, frames=4):
    """A laterally confined column, z in [-1, 0], water table at `table`,
    fixed p = 0 at the base and an inward flux at the top."""
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='column', max_time_in_yr=repr(tend / YEAR),
                      output_time_interval_in_yr=repr(tend / frames / YEAR),
                      max_steps='100000000', output_step_interval='100000000',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    width = repr(2 * res)
    cfg['mesh'].update(xlength=width, ylength=width, zlength='1', resolution=repr(res),
                       largest_size=repr(res), smallest_size=repr(res))
    # Above the water table the soil weighs less than the saturated lithostatic
    # initial stress assumes; the initial PT equilibrium removes that imbalance.
    cfg['control'].update(gravity=repr(G), fixed_dt=repr(dt), has_thermal_diffusion='false',
                          PT_max_iter='100000')
    cfg['bc'].update(vbc_x1='1', hbc_x1='0', hbc_z0='2', hbc_z1='3',
                     hbc_flux=f'[0,0,0,0,0,{-flux!r}]', hbc_pressure='[0,0,0,0,0,0]')
    cfg['ic']['water_table_z'] = repr(table)
    cfg['mat'].update(porosity='[0.4]', hydraulic_perm=f'[{PERM!r}]', fluid_visc=f'[{MU!r}]',
                      fluid_bulk_modulus=f'[{fluid_bulk!r}]', **retention)
    return cfg


def run(exe, cfg, case):
    case.mkdir(parents=True)
    return run_model(exe, case, cfg, 'column', timeout=1800)


# A node's content is rounded at ~eps times its pore volume on each update, with
# random sign. Relative to one step's inflow, that floor grows with nodes and
# substeps, so it is largest for the 2.5 cm columns in 3D.
BALANCE_TOL = 1e-8


def balance_error(model, inflow):
    """Largest per-step balance residual, relative to the larger of the step's
    terms and its imposed inflow (all terms vanish at a steady state)."""
    worst = 0.0
    for frame in model.frames[1:]:
        b = model.read_field(frame, 'fluid balance')
        worst = max(worst, abs(b[4]) / max(np.abs(b[:4]).max(), inflow))
    return worst


def transient(exe, out, ndims):
    """U6: Gardner column against Srivastava-Yeh at T = alpha*Ks*t/dtheta = 0.5."""
    ret = retention_options(Retention('gardner', ALPHA_HEAD / GAMMA))
    sy = dict(length=1.0, alpha=ALPHA_HEAD, ks=KS, dtheta=0.4, qa=0.0, qb=0.9 * KS)
    # The step stays fixed while the spatial mesh is refined.
    # This is an isothermal flow check, so unused heat diffusion is disabled.
    dt = 60.0 if ndims == 2 else 30.0
    rows = []
    for res in (0.05, 0.025):
        model = run(exe, column(ndims, res, 7200.0, dt, retention=ret), out / f'transient-{res}')
        assert model.ndims == ndims, 'the executable and --ndims disagree'
        err = []
        for i, frame in enumerate(model.frames):
            z = model.read_field(frame, 'coordinate')[:, -1] + 1.0
            p = model.read_field(frame, 'pore pressure')
            k = np.exp(ALPHA_HEAD * p / GAMMA)
            err.append(float(np.max(np.abs(k - srivastava_yeh(z, model.time[i], **sy)))))
        inflow = 0.9 * KS * (2 * res)**(ndims - 1) * dt
        rows.append(dict(resolution=res, nodes=len(z), max_error=err, final_time=model.time[-1],
                         p_min=float(p.min()), balance=balance_error(model, inflow)))
    coarse, fine = rows
    assert coarse['p_min'] < -5e3 and fine['p_min'] < -5e3, 'column did not stay unsaturated'
    assert fine['max_error'][-1] <= 2.5e-2, fine
    for a, b in zip(coarse['max_error'][1:], fine['max_error'][1:]):
        assert a / b >= 2.2, (coarse['max_error'], fine['max_error'])
    assert max(r['balance'] for r in rows) <= BALANCE_TOL, rows
    return rows


def modes(exe, out, ndims, reference):
    """The 5 cm transient without PT, as tectonic and cycle runs are set up."""
    ret = retention_options(Retention('gardner', ALPHA_HEAD / GAMMA))
    sy = dict(length=1.0, alpha=ALPHA_HEAD, ks=KS, dtheta=0.4, qa=0.0, qb=0.9 * KS)
    rows = {}
    for gvs in ('false', 'true'):
        cfg = column(ndims, 0.05, 7200.0, 0.0, retention=ret)
        cfg['control'].update(has_PT='false', inertial_scaling='1e5', use_global_velocity_scaling=gvs)
        cfg['ic']['has_body_force_adjustment'] = 'false'
        model = run(exe, cfg, out / f'modes-gvs-{gvs}')
        assert model.ndims == ndims, 'the executable and --ndims disagree'
        frame = model.frames[-1]
        z = model.read_field(frame, 'coordinate')[:, -1] + 1.0
        k = np.exp(ALPHA_HEAD * model.read_field(frame, 'pore pressure') / GAMMA)
        rows[gvs] = float(np.max(np.abs(k - srivastava_yeh(z, model.time[-1], **sy))))
        assert abs(rows[gvs] - reference) <= 0.1 * reference, (rows, reference)
    return rows


def steady(exe, out, ndims):
    """U7: steady van Genuchten and Brooks-Corey profiles for a flux of Ks/2.

    In 2D the error must fall with resolution: at second order for van
    Genuchten and at about first order for Brooks-Corey, whose air-entry kink
    limits it. In 3D the short steps at 2.5 cm would make 2e6 s of explicit
    run too long, so one 5 cm run must stay within 6% of the profile range.
    """
    laws = {'van_genuchten': (Retention('van_genuchten', 1 / 2e3, 2.0, 0.1), 3.0),
            'brooks_corey': (Retention('brooks_corey', 1 / 1e3, 1.5, 0.05), 1.5)}
    rows = {}
    for name, (ret, rate) in laws.items():
        rows[name] = []
        # A compressible fluid shortens the explicit run; it changes no steady
        # profile.
        dt = 300.0 if ndims == 2 else 120.0
        for res in ((0.05, 0.025) if ndims == 2 else (0.05,)):
            cfg = column(ndims, res, 2e6, dt, retention=retention_options(ret), flux=0.5 * KS,
                         fluid_bulk=4e4, frames=2)
            model = run(exe, cfg, out / f'steady-{name}-{res}')
            assert model.ndims == ndims, 'the executable and --ndims disagree'
            frame = model.frames[-1]
            z = model.read_field(frame, 'coordinate')[:, -1] + 1.0
            p = model.read_field(frame, 'pore pressure')
            exact = steady_pressure(ret, z, flux=0.5 * KS, perm=PERM, mu=MU, rho=RHO, g=G)
            inflow = 0.5 * KS * (2 * res)**(ndims - 1) * dt
            rows[name].append(dict(resolution=res, max_error_pa=float(np.max(np.abs(p - exact))),
                                   p_range=float(-exact.min()), balance=balance_error(model, inflow)))
        assert max(r['balance'] for r in rows[name]) <= BALANCE_TOL, rows[name]
        if ndims == 2:
            coarse, fine = rows[name]
            assert coarse['max_error_pa'] / fine['max_error_pa'] >= rate, rows[name]
            assert fine['max_error_pa'] <= 0.02 * fine['p_range'], rows[name]
        else:
            assert rows[name][0]['max_error_pa'] <= 0.06 * rows[name][0]['p_range'], rows[name]
    return rows


def restart(exe, out):
    """Exact continuation of a partly unsaturated column, fixed and remeshed."""
    result = {}
    for mode in ('fixed', 'remesh'):
        cfg = read_config(ROOT / 'tests/functional/pt_remesh_restart.cfg')
        cfg['control']['PT_absolute_tolerance'] = '1e-8'
        if mode == 'fixed':
            cfg['control']['has_moving_mesh'] = 'false'
        cfg['ic']['water_table_z'] = '-0.5'
        cfg['mat'].update(retention_options(Retention('van_genuchten', 2e-4, 1.6, 0.1)))
        restart_check(exe, cfg, out / f'restart-{mode}', 4, 1, require_remesh=mode == 'remesh', timeout=600)
        result[mode] = 'exact'
    return result


def main():
    exe, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    ndims = 3 if exe.name.endswith('3d') else 2
    summary = dict(ndims=ndims, transient=transient(exe, out, ndims))
    summary['modes'] = modes(exe, out, ndims, summary['transient'][0]['max_error'][-1])
    summary.update(steady=steady(exe, out, ndims), restart=restart(exe, out))
    write_json(out / 'summary.json', summary)
    print(json.dumps(summary, indent=2))
    print('PASS: unsaturated transient (PT and dynamic relaxation), steady profiles, exact restart and remeshing')


if __name__ == '__main__':
    main()
