#!/usr/bin/env python3
"""DES Bishop column (unsaturated verification level U8).

A laterally confined elastic column 0 <= z <= H (DES z = z - H) starts drained
and at rest with a hydrostatic water table at TABLE0. The base pressure is
then lowered to that of a table at TABLE1; the sealed top and sides let the
column drain to the new hydrostatic state. Between the two drained states the
Bishop effective stress sigma' = sigma + alpha*Se*p and the bulk weight
rho_s*(1-phi) + rho_f*phi*S change; the exact vertical displacement is
benchmarks/poromechanics/unsaturated_reference.py:bishop_column.

The 'layered' case puts van Genuchten below z = INTERFACE and Gardner above,
inside the zone that drains, so the interface nodes are shared by two laws; its reference integrates the same
formula numerically, law by layer.

Storage does not enter either drained state, so a compressible pore fluid
shortens the explicit drainage. Strains are differences from the drained
starting frame; element-layer strains are integrated from the fixed base.

    check_bishop_column.py EXECUTABLE OUTPUT [--ndims 2] [--spacings 1 0.5 0.25]
"""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks'/'poromechanics')]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from scipy.integrate import cumulative_trapezoid
from unsaturated_reference import Retention, bishop_column, bishop_pressure

G, RHO, MU = 9.81, 1000.0, 1e-3
GAMMA = RHO*G
H, TABLE0, TABLE1 = 10.0, 9.0, 5.0
COLUMN = dict(height=H, phi=0.4, biot=0.9, bulk=2e7, shear=1e7)
GRAIN = COLUMN['bulk']/(1 - COLUMN['biot'])      # consistent with the Biot coefficient
FLUID_BULK, PERM = 1e5, 1e-10                    # compressible pore fluid; permeable column
LAWS = dict(van_genuchten=(Retention('van_genuchten', 1/2e4, 1.6, 0.1), 1),
            brooks_corey=(Retention('brooks_corey', 1/1e4, 0.8, 0.05), 2),
            gardner=(Retention('gardner', 0.5/GAMMA), 3))
LAYERED = ('van_genuchten', 'gardner')          # below and above z = INTERFACE
INTERFACE = 7.0                                  # between TABLE1 and TABLE0
DURATION = 4e5                                   # end time (s)


from des_case import digest, read_config, run_solver, source_diff, source_revision


def retention_lists(laws):
    rets = [LAWS[k][0] for k in laws]
    out = dict(retention_model='[' + ', '.join(str(LAWS[k][1]) for k in laws) + ']',
               retention_alpha='[' + ', '.join(repr(r.alpha) for r in rets) + ']',
               residual_saturation='[' + ', '.join(repr(r.sr) for r in rets) + ']',
               retention_n='[' + ', '.join(repr(r.n if r.n is not None else 2.0) for r in rets) + ']')
    return out


def layered_column(z):
    """bishop_column for the two-law column, by integration on a fine grid."""
    zf = np.unique(np.concatenate((np.linspace(0, H, 400001), z)))
    upper = zf >= INTERFACE
    dpi, dse = np.empty_like(zf), np.empty_like(zf)
    for side, key in ((~upper, LAYERED[0]), (upper, LAYERED[1])):
        ret = LAWS[key][0]
        p0, p1 = -GAMMA*(zf[side] - TABLE0), -GAMMA*(zf[side] - TABLE1)
        dpi[side] = bishop_pressure(ret, p1) - bishop_pressure(ret, p0)
        dse[side] = (1 - ret.sr)*(ret.se(p1) - ret.se(p0))
    above = cumulative_trapezoid(dse[::-1], zf[::-1], initial=0)[::-1]    # -int_z^H dSe (reversed grid)
    dsig = GAMMA*COLUMN['phi']*above
    eps = (dsig + COLUMN['biot']*dpi)/(COLUMN['bulk'] + 4*COLUMN['shear']/3)
    return np.interp(z, zf, cumulative_trapezoid(eps, zf, initial=0))


def reference(law, z):
    if law == 'layered': return layered_column(z)
    return bishop_column(LAWS[law][0], z, table0=TABLE0, table1=TABLE1, rho=RHO, g=G, **COLUMN)


def config(ndims, h, law, steps):
    laws = LAYERED if law == 'layered' else (law,)
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    dt, hold = DURATION/steps, max(1, steps//10)    # first output frame is the drained start
    cfg['sim'].update(modelname='bishop', max_steps=str(steps), max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval=str(hold), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(h), ylength=repr(h), zlength=repr(H), resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(gravity=repr(G), fixed_dt=repr(dt), has_thermal_diffusion='false', PT_max_iter='100000',
                          has_poroelastic_pressure_feedback='true')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'] = '1', '0'
    bc.update(hbc_x1='0', hbc_z0='2', hbc_z1='0', hbc_pressure=f'[0,0,0,0,{GAMMA*TABLE0!r},0]',
              num_hbc_schedule='3',
              hbc_schedule_time_in_yr='[0, ' + ', '.join(repr(k*dt/YEAR2SEC) for k in (hold, hold+1)) + ']',
              hbc_schedule_ratio=f'[1, 1, {TABLE1/TABLE0!r}]')
    cfg['ic'].update(water_table_z=repr(TABLE0 - H), pore_pressure_option='0')
    n = len(laws)
    if n > 1:
        cfg['ic'].update(mattype_option='1', num_mattype_layers='2', layer_mattypes='[1, 0]',
                         mattype_layer_depths=f'[{(H - INTERFACE)/H!r}]')

    def each(value):
        return '[' + ', '.join([repr(value)]*n) + ']'
    cfg['mat'].update(num_materials=str(n), rho0=each(2650.0), bulk_modulus=each(COLUMN['bulk']),
                      shear_modulus=each(COLUMN['shear']), biot_coeff=each(COLUMN['biot']),
                      porosity=each(COLUMN['phi']), bulk_modulus_s=each(GRAIN), fluid_bulk_modulus=each(FLUID_BULK),
                      hydraulic_perm=each(PERM), fluid_visc=each(MU), fluid_rho0=each(RHO), fluid_alpha=each(0.0),
                      **retention_lists(laws))
    return cfg


def layer_strain(coord, conn, ezz, h):
    """Volume-weighted mean vertical strain of each element layer, base upward."""
    z = coord[conn][:, :, -1].mean(axis=1) + H
    layer = np.floor(z/h).astype(int)
    n = int(round(H/h))
    vol = np.bincount(layer, minlength=n).astype(float)   # equal-volume elements in a structured column
    assert len(vol) == n and vol.min() > 0, 'the column mesh must be structured in layers of height h'
    return np.bincount(layer, weights=ezz, minlength=n)/vol


def compare(folder, h, law):
    model = Dynearthsol(str(folder/'bishop')); rows = scan_frames(str(folder/'bishop'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    base = rows[1]                                  # step steps//10, before the base pressure drops
    last = rows[-1]
    z = coord[:, -1] + H
    p = model.read_field(last['frame'], 'pore pressure')
    drained = float(np.max(np.abs(p - (-GAMMA*(z - TABLE1)))))/(GAMMA*(TABLE0 - TABLE1))
    zz = coord.shape[1] - 1                         # strain is (xx, zz, xz) in 2D, (xx, yy, zz, ...) in 3D
    e0 = layer_strain(coord, conn, model.read_field(base['frame'], 'strain')[:, zz], h)
    e1 = layer_strain(coord, conn, model.read_field(last['frame'], 'strain')[:, zz], h)
    zn = np.arange(len(e0) + 1)*h
    u = np.concatenate(([0], np.cumsum((e1 - e0)*h)))
    ref = reference(law, zn)
    scale = float(np.max(np.abs(ref)))
    return dict(spacing_m=h, law=law, baseline_step=base['steps'], final_time=last['time'],
                drained_residual=drained, top_des_m=float(u[-1]), top_ref_m=float(ref[-1]),
                max_error=float(np.max(np.abs(u - ref)))/scale, z=zn.tolist(), u_des=u.tolist(), u_ref=ref.tolist())


def study(exe, out, ndims, spacings, steps, laws):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=source_revision(ROOT),
                  ndims=ndims, steps=steps, duration_s=DURATION, column=COLUMN, tables=[TABLE0, TABLE1], cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for law in laws:
        for h in spacings:
            folder = out/f'{law}-h{h:g}'; folder.mkdir()
            with (folder/'bishop.cfg').open('w') as f: config(ndims, h, law, steps).write(f)
            run_solver([str(exe), 'bishop.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=20000)
            item = compare(folder, h, law)
            item['input_sha256'] = digest(folder/'bishop.cfg')
            result['cases'].append(item)
            (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
            print(json.dumps({k: v for k, v in item.items() if k not in ('z', 'u_des', 'u_ref')}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, default=2)
    p.add_argument('--spacings', type=float, nargs='+', default=[1, 0.5, 0.25])
    p.add_argument('--steps', type=int, default=40)
    p.add_argument('--laws', nargs='+', default=list(LAWS) + ['layered'])
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.ndims, a.spacings, a.steps, a.laws)
