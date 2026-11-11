#!/usr/bin/env python3
"""DES Hantush-Jacob pumping from a leaky confined aquifer, on refined meshes.

A quarter box (x, y >= 0) holds a confined aquifer of thickness B under an
aquitard one element (h) thick whose top is held at p = 0. The aquitard
leakance k'_z/(mu b') is the same on every mesh, so the leakage factor
lambda^2 = k B b'/k'_z is fixed; its storage is a thousandth of the aquifer's.
A fully penetrating well on the corner edge pumps a quarter of Q through the
aquifer, then shuts in. The skeleton is stiff and pressure feedback is off,
as in check_theis.py.

Hantush and Jacob neglect vertical resistance in the aquifer, aquitard
storage and horizontal aquitard flow. Vertical permeability is ANISOTROPY
times horizontal in both layers, so the aquifer's vertical pressure drop is
B^2/(2 ANISOTROPY lambda^2) = 0.14% of the drawdown. The aquitard pressure
is linear in z, which one element represents exactly; its storage S'b'/(S B)
and horizontal transmissivity (b'/lambda)^2/ANISOTROPY shrink with h.
Drawdown is compared with benchmarks/poromechanics/hantush_reference.py.

    check_hantush.py EXECUTABLE OUTPUT [--spacings 2.5 1.25 0.625]
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
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from hantush_reference import scheduled_drawdown
import theis_reference

L, B, Q, LAMBDA = 75.0, 2.5, 1e-4, 15.0   # box side, aquifer thickness, total rate, leakage factor (m)
ANISOTROPY = 10.0                         # vertical / horizontal permeability
PUMP, END, DT = 200.0, 320.0, 2.0         # shut-in time, end time, physical step (s)
WINDOW = (10.0, 30.0)
AQUIFER = dict(porosity=.25, biot=1.0)
AQUITARD = dict(porosity=1e-4, biot=.3)
SOLID = dict(bulk=1e12, shear=1e12, fluid_bulk=2.2e9, grain_bulk=1e20, visc=1e-3)


from des_case import digest, read_config, reset_boundaries, run_solver, source_diff, source_revision

def storage(layer):
    return (layer['porosity']/SOLID['fluid_bulk'] + (layer['biot']-layer['porosity'])/SOLID['grain_bulk']
            + layer['biot']**2/SOLID['bulk'])


STORAGE, STORAGE_AQUITARD = storage(AQUIFER), storage(AQUITARD)
DIFFUSIVITY = 1.0                         # aquifer, m2/s
PERM = DIFFUSIVITY*STORAGE*SOLID['visc']
T = PERM*B/SOLID['visc']                  # transmissivity in pressure units, m3/(Pa s)
S = STORAGE*B                             # storativity in pressure units, m/Pa
SCALE = Q/(4*math.pi*T)
SCHEDULE = [(0.0, Q), (PUMP, 0.0)]


def aquitard_perm(h):
    """Scalar k' of an aquitard h thick whose vertical leakance is ANISOTROPY*k'/h = k B/lambda^2."""
    return PERM*B*h/LAMBDA**2/ANISOTROPY


def config(h):
    cfg = read_config(ROOT / 'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    steps = int(round(END/DT))
    cfg['sim'].update(modelname='hantush', max_steps=str(steps), output_step_interval=str(steps//8),
                      max_time_in_yr=repr(2*END/YEAR2SEC), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='4')
    cfg['mesh'] = dict(meshing_option='1', meshing_elem_shape='1', xlength=repr(L), ylength=repr(L),
                       zlength=repr(B + h), resolution=repr(h), remeshing_option='0')
    cfg['control'].update(gravity='0', has_PT='false', is_quasi_static='true', inertial_scaling='1e6',
                          use_global_velocity_scaling='false', has_moving_mesh='false', damping_option='0',
                          fixed_dt=repr(DT), has_poroelastic_pressure_feedback='false',
                          has_thermal_diffusion='false')
    bc = cfg['bc']
    reset_boundaries(bc, vbc='1', vbc_val='0', hbc='0')
    bc['hbc_z1'] = '2'
    t = [0.0, PUMP, PUMP]
    bc.update(hbc_pressure='[0,0,0,0,0,0]', stress_bc_z1='0', stress_val_z1='0', num_wells='1',
              well_x='[0.0]', well_y='[0.0]', well_z_top=f'[{-h!r}]', well_z_bottom=f'[{-(B + h)!r}]',
              well_rate=f'[{Q/4!r}]', num_well_schedule='3',
              well_schedule_time_in_yr='['+', '.join(repr(x/YEAR2SEC) for x in t)+']',
              well_schedule_ratio='[1, 1, 0]')
    cfg['ic'].update(has_body_force_adjustment='false', pore_pressure_option='0', mattype_option='1',
                     num_mattype_layers='2', layer_mattypes='[1, 0]', mattype_layer_depths=f'[{h/(B + h)!r}]')

    def pair(a, b):
        return f'[{a!r}, {b!r}]'
    cfg['mat'].update(num_materials='2', rho0='[2650, 2650]', fluid_rho0='[1000, 1000]', fluid_alpha='[0, 0]',
                      bulk_modulus=pair(SOLID['bulk'], SOLID['bulk']), shear_modulus=pair(SOLID['shear'], SOLID['shear']),
                      biot_coeff=pair(AQUIFER['biot'], AQUITARD['biot']),
                      porosity=pair(AQUIFER['porosity'], AQUITARD['porosity']),
                      fluid_bulk_modulus=pair(SOLID['fluid_bulk'], SOLID['fluid_bulk']),
                      bulk_modulus_s=pair(SOLID['grain_bulk'], SOLID['grain_bulk']),
                      hydraulic_perm=pair(PERM, aquitard_perm(h)), fluid_visc=pair(SOLID['visc'], SOLID['visc']),
                      permeability_tensor=f'[1, 1, {ANISOTROPY!r}, 0, 0, 0]', retention_model='[0, 0]')
    return cfg


def compare(folder, h):
    model = Dynearthsol(str(folder/'hantush')); rows = scan_frames(str(folder/'hantush'))
    coord = model.read_field(0, 'coordinate')
    r = np.hypot(coord[:, 0], coord[:, 1])
    # Aquifer nodes below the interface, same radii on every mesh: no axis, no far field.
    window = (r >= WINDOW[0]) & (r <= WINDOW[1]) & (coord[:, 2] < -h + 1e-9)
    p0 = model.read_field(0, 'pore pressure')
    frames = []
    for row in rows[1:]:
        t = row['time']
        des = p0 - model.read_field(row['frame'], 'pore pressure')
        ref = scheduled_drawdown(r[window], np.full(window.sum(), t), SCHEDULE, T, S, LAMBDA)
        theis = theis_reference.scheduled_drawdown(r[window], np.full(window.sum(), t), SCHEDULE, T, S)
        frames.append(dict(time=t, max_error_over_scale=float(np.max(np.abs(des[window]-ref)))/SCALE,
                           rms_error_over_scale=float(np.sqrt(np.mean((des[window]-ref)**2)))/SCALE,
                           leakage_effect_over_scale=float(np.max(np.abs(theis-ref)))/SCALE))
    return model, frames, len(r)


def study(exe, out, spacings):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=source_revision(ROOT),
                  reference=dict(L=L, B=B, Q=Q, leakage_factor=LAMBDA, anisotropy=ANISOTROPY, pump_until_s=PUMP, end_s=END, dt_s=DT,
                                 storage_per_Pa=STORAGE, aquitard_storage_per_Pa=STORAGE_AQUITARD,
                                 diffusivity=DIFFUSIVITY, permeability=PERM, transmissivity=T, storativity=S,
                                 pressure_scale_Pa=SCALE, window=WINDOW), cases=[])
    (out/'source.diff').write_bytes(source_diff(ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'hantush.cfg').open('w') as f: config(h).write(f)
        run_solver([str(exe), 'hantush.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=20000)
        model, frames, nodes = compare(folder, h)
        strain = max(float(np.abs(model.read_field(f, 'strain')).max()) for f in model.frames)
        item = dict(spacing_m=h, nodes=nodes, frames=frames,
                    aquitard_storage_ratio=STORAGE_AQUITARD*h/(STORAGE*B), aquitard_transmissivity_ratio=(h/LAMBDA)**2/ANISOTROPY,
                    max_error_over_scale=max(f['max_error_over_scale'] for f in frames),
                    rms_error_over_scale=max(f['rms_error_over_scale'] for f in frames),
                    leakage_effect_over_scale=max(f['leakage_effect_over_scale'] for f in frames),
                    max_abs_strain=strain, input_sha256=digest(folder/'hantush.cfg'))
        assert strain < 1e-6, ('stiff skeleton must stay rigid', strain)
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[2.5, 1.25, 0.625])
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings)
