#!/usr/bin/env python3
"""DES on Elder's heater problem against the independent transient reference.

The symmetric half of Elder's box, 2 m wide and 1 m high: the base is held
10 K hot for 0 <= x <= 1 m (bc.tbc_z0 on segment flag 16) and insulated
beyond it (slant boundary n0, flag 64, with zero normal velocity); the top is
cold, the symmetry axis x = 0 and the far wall are insulated, and every
boundary is impermeable. The fluid and matrix constants are those of
check_convection.py, so kappa = 1e-6 m^2/s, sigma = 1 and tau = 1e-6 t/s.
Ra = 400 sets the permeability. The rigid matrix does not feed back on the
pressure; the fluid bulk modulus sets how much faster pressure diffuses than
heat (the reference is incompressible).

Each case writes a frame at every tau of elder_reference.TIMES, and the
temperature is compared node by node with the finest reference grid,
interpolated bilinearly. The heat content of the layer is compared too.

    check_elder.py EXECUTABLE REFERENCE.npz OUTPUT [--spacings 0.0625 0.03125]
                   [--dt 50] [--fluid-bulk 100] [--threads 8]
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_convection as cc
from convection_reference import rayleigh

RA, WIDTH, HEIGHT = 400.0, 2.0, 1.0
KAPPA = cc.LAMBDA/(cc.RHO*cc.CF)


def poly(h):
    """Half box with the heater and the insulated base on separate boundary flags."""
    area = math.sqrt(3)/4*h*h
    return '\n'.join([
        '#### node coordinates ####',
        '5 2 0 0',
        '0 0 0', f'1 0 {-HEIGHT}', f'2 1 {-HEIGHT}', f'3 {WIDTH} {-HEIGHT}', f'4 {WIDTH} 0',
        '#### segments ####',
        '5 1',
        '0 0 1 1      # x0: symmetry axis',
        '1 1 2 16     # z0: heater',
        '2 2 3 64     # n0: insulated base',
        '3 3 4 2      # x1: far wall',
        '4 4 0 32     # z1: cold top',
        '#### holes, must be 0 ####',
        '0',
        '#### regions ####',
        '1',
        f'0 1 {-HEIGHT/2} 0 {area!r}',
        ''])


def config(h, dt, fluid_bulk, outputs):
    k = RA/rayleigh(cc.RHO, cc.G, cc.ALPHA, cc.DT, 1, HEIGHT, cc.CF, cc.MU, cc.LAMBDA)
    steps = [round(t*HEIGHT**2/KAPPA/dt) for t in outputs]
    assert all(abs(s*dt - t*HEIGHT**2/KAPPA) < 1e-9*s*dt for s, t in zip(steps, outputs)), 'dt must divide each time'
    interval = math.gcd(*steps)
    cfg = cc.config(h, 1.0, 0.0, steps[-1], dt, 1)
    cfg['sim'].update(modelname='elder', output_step_interval=str(interval), checkpoint_frame_interval='100000')
    cfg['mesh'].update(meshing_option='90', poly_filename='elder.poly', xlength=str(WIDTH))
    cfg['ic'].update(temperature_option='4', temp_background=repr(cc.T_TOP), temp_anomaly_amplitude='0')
    cfg['bc'].update(tbc_z0='1', tbc_val_z0=repr(cc.T_TOP + cc.DT), tbc_z1='1', tbc_val_z1=repr(cc.T_TOP),
                     tbc_x0='0', tbc_x1='0', vbc_n0='1', vbc_val_n0='0')
    cfg['mat'].update(hydraulic_perm=f'[{k!r}]', fluid_bulk_modulus=f'[{fluid_bulk!r}]')
    return cfg, interval


def compare(folder, reference, interval, dt):
    model = cc.Dynearthsol(str(folder/'elder')); rows = cc.scan_frames(str(folder/'elder'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    w = cc.lumped_area(coord, conn)
    n = max(int(k.split('_')[1]) for k in reference.files if k.startswith('theta_'))
    xr, zr, snaps = reference[f'x_{n}'], reference[f'z_{n}'], reference[f'theta_{n}']
    points = np.column_stack((coord[:, 1] + HEIGHT, coord[:, 0]))
    by_step = {r['steps']: r for r in rows}
    out = []
    for tau, ref in zip(reference['times'], snaps):
        row = by_step[round(tau*HEIGHT**2/KAPPA/dt)]
        assert abs(row['time']*KAPPA/HEIGHT**2 - tau) < 1e-9
        theta = (model.read_field(row['frame'], 'temperature') - cc.T_TOP)/cc.DT
        exact = RegularGridInterpolator((zr, xr), ref)(points)
        diff = theta - exact
        hr = np.full(len(xr), xr[1]-xr[0]); hr[[0, -1]] /= 2
        hz = np.full(len(zr), zr[1]-zr[0]); hz[[0, -1]] /= 2
        out.append(dict(tau=float(tau), max_abs_difference=float(np.abs(diff).max()),
                        rms_difference=float(math.sqrt(w @ diff**2/w.sum())),
                        heat_des=float(w @ theta), heat_reference=float(hz @ ref @ hr),
                        theta_min=float(theta.min()), theta_max=float(theta.max())))
    return out, len(coord), len(conn)


def study(exe, reference_file, out, spacings, dt, fluid_bulk, threads):
    out.mkdir(parents=True, exist_ok=False)
    reference = np.load(reference_file)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=cc.digest(exe), runner_sha256=cc.digest(Path(__file__)),
                  reference_sha256=cc.digest(reference_file), ra=RA, dt_s=dt, fluid_bulk_Pa=fluid_bulk,
                  kappa_m2_s=KAPPA, threads=threads, cases=[])
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        (folder/'elder.poly').write_text(poly(h))
        cfg, interval = config(h, dt, fluid_bulk, reference['times'])
        with (folder/'elder.cfg').open('w') as f: cfg.write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'elder.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=200000)
        series, nodes, elements = compare(folder, reference, interval, dt)
        item = dict(spacing_m=h, nodes=nodes, elements=elements, series=series)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps(dict(spacing_m=h, nodes=nodes, max_abs=[round(s['max_abs_difference'], 4) for s in series],
                              heat=[(round(s['heat_des'], 5), round(s['heat_reference'], 5)) for s in series])),
              flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('reference', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.0625, 0.03125])
    p.add_argument('--dt', type=float, default=50.0)
    p.add_argument('--fluid-bulk', type=float, default=100.0)
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    study(a.executable.resolve(), a.reference.resolve(), a.output.resolve(), a.spacings, a.dt, a.fluid_bulk, a.threads)
