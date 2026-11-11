#!/usr/bin/env python3
"""DES steady convection beyond onset: the Nusselt number of the square cell.

The Horton-Rogers-Lapwood box of check_convection.py is run at Ra above
onset until its heat transport stops changing. Storage does not enter the
steady state, so a more compressible pore fluid (pressure diffusing ~20
times faster than heat instead of ~2000) and a larger seed shorten the run.

Nu = 1 + rho_f c_f H <q_z T>/(lambda dT) over the box, with the element
Darcy flux q = -(k/mu)(grad p + rho_f (1 - alpha (T - T0)) g e_z) rebuilt
from the saved pressure and temperature; a wall estimate from the top
element layer is reported alongside. The reference is the independent
steady solution in benchmarks/poromechanics/convection_steady.py.

    check_convection_nusselt.py EXECUTABLE OUTPUT [--ra 100 200] [--spacings 0.05 0.025]
"""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_convection as cc
from convection_reference import critical_rayleigh
from convection_steady import nusselt

FLUID_BULK = 1e1                     # Pa; storage only sets the transient
DURATION, OUTPUTS = 1.6e6, 32        # s: steady by about 0.5e6 s at h = 0.05 m


from des_case import digest, run_solver, source_diff, source_revision


def config(h, ra, dt):
    steps = int(round(DURATION/dt))
    cfg = cc.config(h, ra/critical_rayleigh(1, 1), 0.5, steps, dt, OUTPUTS)
    cfg['mat']['fluid_bulk_modulus'] = f'[{FLUID_BULK!r}]'
    cfg['sim']['checkpoint_frame_interval'] = '100000'
    return cfg


def element_gradients(coord, conn):
    """Constant gradients of linear triangles: (nelem, 3 nodes, 2) and areas."""
    x, z = coord[conn][:, :, 0], coord[conn][:, :, 1]
    det = (x[:, 1]-x[:, 0])*(z[:, 2]-z[:, 0]) - (x[:, 2]-x[:, 0])*(z[:, 1]-z[:, 0])
    g = np.empty(conn.shape + (2,))
    g[:, 0, 0], g[:, 1, 0], g[:, 2, 0] = z[:, 1]-z[:, 2], z[:, 2]-z[:, 0], z[:, 0]-z[:, 1]
    g[:, 0, 1], g[:, 1, 1], g[:, 2, 1] = x[:, 2]-x[:, 1], x[:, 0]-x[:, 2], x[:, 1]-x[:, 0]
    return g/det[:, None, None], np.abs(det)/2


def nusselt_des(folder, ra):
    model = cc.Dynearthsol(str(folder/'box')); rows = cc.scan_frames(str(folder/'box'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    grad, area = element_gradients(coord, conn)
    k = ra/critical_rayleigh(1, 1)*critical_rayleigh(1, 1)/cc.rayleigh(cc.RHO, cc.G, cc.ALPHA, cc.DT, 1, 1, cc.CF, cc.MU, cc.LAMBDA)
    top = np.abs(coord[conn][:, :, 1]).min(axis=1) < 1e-9       # elements touching the top
    out = []
    for row in rows[1:]:
        p = model.read_field(row['frame'], 'pore pressure'); T = model.read_field(row['frame'], 'temperature')
        Te = T[conn].mean(axis=1)
        dp = np.einsum('ejd,ej->ed', grad, p[conn]); dT = np.einsum('ejd,ej->ed', grad, T[conn])
        rho = cc.RHO*(1 - cc.ALPHA*(Te - cc.T_TOP))
        qz = -(k/cc.MU)*(dp[:, 1] + rho*cc.G)
        nu = 1 + cc.RHO*cc.CF*float(area @ (qz*Te))/(cc.LAMBDA*cc.DT)
        wall = float(-(area[top] @ dT[top, 1])/area[top].sum())/cc.DT
        out.append(dict(time=row['time'], nu=nu, nu_wall=wall))
    return out


def study(exe, out, ras, spacings, dt):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=source_revision(cc.ROOT),
                  dt=dt, duration_s=DURATION, fluid_bulk_Pa=FLUID_BULK,
                  reference={f'{ra:g}': nusselt(ra, cache=out/f'reference-ra{ra:g}.npz') for ra in ras}, cases=[])
    result['reference_sha256'] = {f'{ra:g}': digest(out/f'reference-ra{ra:g}.npz') for ra in ras}
    (out/'source.diff').write_bytes(source_diff(cc.ROOT))
    for ra in ras:
        ref = result['reference'][f'{ra:g}']['nu_extrapolated']
        for h in spacings:
            folder = out/f'ra{ra:g}-h{h:g}'; folder.mkdir()
            with (folder/'box.cfg').open('w') as f: config(h, ra, dt).write(f)
            run_solver([str(exe), 'box.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=40000)
            series = nusselt_des(folder, ra)
            last = series[-1]
            item = dict(ra=ra, spacing_m=h, series=series, nu=last['nu'], nu_wall=last['nu_wall'], nu_ref=ref,
                        error=abs(last['nu']/ref - 1), wall_error=abs(last['nu_wall']/ref - 1),
                        steady_change=abs(series[-1]['nu'] - series[-5]['nu'])/last['nu'])
            result['cases'].append(item)
            (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
            print(json.dumps({k: v for k, v in item.items() if k != 'series'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ra', type=float, nargs='+', default=[100.0, 200.0])
    p.add_argument('--spacings', type=float, nargs='+', default=[0.05, 0.025])
    p.add_argument('--dt', type=float, default=1000.0)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.ra, a.spacings, a.dt)
