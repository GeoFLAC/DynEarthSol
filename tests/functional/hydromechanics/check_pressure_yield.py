#!/usr/bin/env python3
"""DES pressure-driven Mohr-Coulomb yielding of a confined saturated layer.

The layer of benchmarks/poromechanics/pressure_yield_reference.py: 2 m
high, laterally confined, fixed base, top loaded by a surcharge q, saturated
with its water table at the top and the DES legacy isotropic initial stress.
The pore pressure at the top rises linearly to 95% of the pressure at which
the yielded top would reach the Mohr-Coulomb apex (beyond it no equilibrium
exists: the layer would lift off); the base is sealed and the permeability
high, so the pressure rises almost uniformly.
Each element's yield pressure is its own pressure rise (mean of its nodes)
at the first saved step with plastic strain, which removes any diffusion
lag; it is compared with the closed-form onset at the element's depth. The
step bounds the measurement to one pressure increment.

    check_pressure_yield.py EXECUTABLE OUTPUT [--ndims 2] [--spacings 1 0.5 0.25] [--steps 1000]
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
from pressure_yield_reference import PressureYield

H, Q, DURATION = 2.0, 2e4, 2000.0                    # height (m), surcharge (Pa), ramp (s)
MAT = dict(bulk=1e8, shear=6e7, biot=1.0, porosity=.3, fluid_bulk=2.2e9, grain_bulk=1e20, permeability=1e-11,
           viscosity=1e-3, rho_solid=2650.0, rho_fluid=1000.0, cohesion=1e4, friction=20.0)
RHO_BULK = MAT['rho_solid']*(1 - MAT['porosity']) + MAT['rho_fluid']*MAT['porosity']
REF = PressureYield(bulk=MAT['bulk'], shear=MAT['shear'], biot=MAT['biot'], cohesion=MAT['cohesion'],
                    friction_deg=MAT['friction'], rho_bulk=RHO_BULK, rho_fluid=MAT['rho_fluid'], surcharge=Q)
DP = 0.95*float(REF.apex(0.0))                     # final pressure rise (Pa)


from des_case import digest, read_config, run_solver, source_diff, source_revision


def config(ndims, h, steps, mode='pt', scaling=1e5, pt_option=None):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='layer', max_steps=str(steps), max_time_in_yr=repr(2*DURATION/YEAR2SEC),
                      output_step_interval='1', output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='100000', has_initial_checkpoint='no')
    cfg['mesh'].update(meshing_elem_shape='1', xlength=repr(h), ylength=repr(h), zlength=repr(H), resolution=repr(h),
                       largest_size=repr(h), smallest_size=repr(h))
    cfg['control'].update(gravity='9.81', fixed_dt=repr(DURATION/steps), has_thermal_diffusion='false',
                          has_poroelastic_pressure_feedback='true', PT_max_iter='100000')
    if pt_option is not None:
        cfg['control']['PT_option'] = str(pt_option)
    if mode != 'pt':
        # Explicit dynamic relaxation or GVS at the stable step; output keeps the PT cadence.
        cfg['control'].update(has_PT='false', fixed_dt='0', inertial_scaling=repr(scaling), characteristic_speed='1e-5',
                              dt_fraction='0.5', use_global_velocity_scaling=str(mode == 'gvs').lower())
        cfg['sim'].update(max_steps=str(10**9), output_step_interval=str(10**9),
                          output_time_interval_in_yr=repr(DURATION/steps/YEAR2SEC))
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'], bc[f'hbc_{side}'] = '1', '0', '0'
    bc.update(stress_bc_z1='3', stress_val_z1=repr(-Q), hbc_z1='2', hbc_pressure=f'[0,0,0,0,0,{DP!r}]', num_hbc_schedule='2',
              hbc_schedule_time_in_yr=f'[0, {DURATION/YEAR2SEC!r}]', hbc_schedule_ratio='[0, 1]')
    cfg['ic'].update(water_table_z='0', pore_pressure_option='0', weakzone_option='0')
    m = MAT
    cfg['mat'].update(rheology_type='elasto-plastic', rho0=f"[{m['rho_solid']!r}]", bulk_modulus=f"[{m['bulk']!r}]",
                      shear_modulus=f"[{m['shear']!r}]", biot_coeff=f"[{m['biot']!r}]", porosity=f"[{m['porosity']!r}]",
                      fluid_bulk_modulus=f"[{m['fluid_bulk']!r}]", bulk_modulus_s=f"[{m['grain_bulk']!r}]",
                      hydraulic_perm=f"[{m['permeability']!r}]", fluid_visc=f"[{m['viscosity']!r}]",
                      fluid_rho0=f"[{m['rho_fluid']!r}]", fluid_alpha='[0]', retention_model='[0]',
                      cohesion0=f"[{m['cohesion']!r}]", cohesion1=f"[{m['cohesion']!r}]",
                      friction_angle0=f"[{m['friction']!r}]", friction_angle1=f"[{m['friction']!r}]",
                      dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e9')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'layer')); rows = scan_frames(str(folder/'layer'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    depth = -coord[conn][:, :, -1].mean(axis=1)
    p0 = model.read_field(0, 'pore pressure')[conn].mean(axis=1)
    yielded = np.full(len(conn), np.nan); before = np.full(len(conn), np.nan)
    previous = np.zeros(len(conn))
    for row in rows[1:]:
        dp = model.read_field(row['frame'], 'pore pressure')[conn].mean(axis=1) - p0
        new = np.isnan(yielded) & (model.read_field(row['frame'], 'plastic strain') > 0)
        yielded[new], before[new] = dp[new], previous[new]
        previous = dp
    done = ~np.isnan(yielded)
    ref = REF.onset(depth)
    reachable = ref < .95*DP
    assert np.all(done[reachable]), 'elements below the final pressure did not yield'
    assert not np.any(done & (ref > 1.05*DP)), 'elements yielded beyond the analytical reach'
    err = np.abs(yielded[done] - ref[done])/ref[done]
    step = (yielded[done] - before[done])/ref[done]          # one-increment bound
    return dict(elements=int(len(conn)), yielded=int(done.sum()), max_error=float(err.max()),
                mean_error=float(err.mean()), max_increment=float(step.max()),
                onset_top_des=float(np.min(yielded[done])), onset_top_ref=float(REF.onset(depth[done].min())),
                profile=dict(depth=depth[done].tolist(), des=yielded[done].tolist(), ref=ref[done].tolist()))


def study(exe, out, ndims, spacings, steps, mode='pt', scaling=1e5, pt_option=None):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=source_revision(ROOT),
                  ndims=ndims, steps=steps, mode=mode, scaling=scaling, pt_option=pt_option, material=MAT, rho_bulk=RHO_BULK, kappa=REF.kappa, surcharge=Q, dp_max=DP, cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'layer.cfg').open('w') as f: config(ndims, h, steps, mode, scaling, pt_option).write(f)
        run_solver([str(exe), 'layer.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=40000)
        item = dict(spacing_m=h, **compare(folder), input_sha256=digest(folder/'layer.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'profile'}), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--ndims', type=int, default=2)
    p.add_argument('--spacings', type=float, nargs='+', default=[0.2, 0.1, 0.05])
    p.add_argument('--steps', type=int, default=1000, help='PT steps, or output frames otherwise')
    p.add_argument('--mode', choices=('pt', 'explicit', 'gvs'), default='pt')
    p.add_argument('--scaling', type=float, default=1e5, help='inertial scaling for explicit and GVS')
    p.add_argument('--pt-option', type=int, choices=(0, 1), default=None, help='PT_option (default: the base input)')
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.ndims, a.spacings, a.steps, a.mode, a.scaling, a.pt_option)
