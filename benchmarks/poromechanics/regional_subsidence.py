#!/usr/bin/env python3
"""Application: seasonal pumping and land subsidence in a soft confined aquifer.

A 50 m confined aquifer of soft basin sediment is pumped by a well field at
0.02 m^3/s for six months, rests for six months, and repeats once. DES models
a quarter of a 6 km square slab (symmetry rollers at x = 0 and y = 0, roller
base, free top, far sides held at p = 0 with fixed normal displacement) with a
full-height well on its corner edge, and PT mechanics. Gravity is off: the
fields are changes from the initial state, which is what the averaged-aquifer
model (averaged_aquifer_reference.py, generalised Bear-Corapcioglu) gives.
Its finite-domain solution (fixed outer pressure and displacement, which
brings in the uniform areal strain h(t) of the closed basin) is superposed
over the pumping schedule, on the circle of the same area as the square.
The infinite-domain line well misses h(t) and departs by up to ~20% in
the rest seasons.

Subsidence is B times the vertical strain of the slab; the observation wells
and the compared band stay inside the diffusion front of the far sides.

    regional_subsidence.py EXECUTABLE OUTPUT [--spacings 100 50 25] [--steps 104]
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
sys.path[:0] = [str(ROOT), str(Path(__file__).resolve().parent)]
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames
from averaged_aquifer_reference import AveragedAquifer

L, B, Q = 3000.0, 50.0, 0.02                  # quarter-domain side, thickness (m), well-field rate (m^3/s)
SEASON = 0.5*YEAR2SEC                         # pumping and resting periods (s)
SCHEDULE = [(0.0, Q), (SEASON, 0.0), (2*SEASON, Q), (3*SEASON, 0.0)]
DURATION = 4*SEASON
PROPS = dict(thickness=B, bulk=5e6, shear=3e6, biot=0.9, porosity=0.3, fluid_bulk=2.2e9, grain_bulk=3.7e10,
             permeability=1e-12, viscosity=1e-3)
BAND = (250.0, 1500.0)                        # compared radii (m)
WELLS = (250.0, 500.0, 1000.0)                # observation radii (m)
REF = AveragedAquifer(**PROPS)


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


OUTER = 2*L/math.sqrt(math.pi)                # radius of the circle with the quarter square's area (m)
WELL_RADIUS = 1.0                             # m; negligible at the compared radii


def reference(radius, t):
    """Superposed finite-domain fields for the schedule at one time t."""
    radius = np.atleast_1d(np.asarray(radius, float))
    out = {k: np.zeros(radius.shape) for k in ('pressure', 'areal_strain', 'vertical_strain', 'subsidence')}
    previous = 0.0
    for start, rate in SCHEDULE:
        # A change that has not yet acted (a frame on the change, up to rounding) adds nothing.
        if t - start > 1e-9*SEASON and rate != previous:
            part = REF.finite_well(radius, [t - start], rate - previous, well_radius=WELL_RADIUS, outer_radius=OUTER)
            for k in out:
                out[k] += np.asarray(part[k]).reshape(-1)
        previous = rate
    return out


def config(h, steps):
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(ROOT/'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    cfg['sim'].update(modelname='basin', max_steps=str(steps), output_step_interval=str(steps//16),
                      max_time_in_yr=repr(2*DURATION/YEAR2SEC), output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='4')
    # Structured: the unstructured mesher fails on a slab 60 times wider than thick.
    cfg['mesh'] = dict(meshing_option='1', meshing_elem_shape='1', xlength=repr(L), ylength=repr(L),
                       zlength=repr(B), resolution=repr(h), remeshing_option='0')
    # Adaptive dynamic relaxation, as in the qualified slab study; the accelerated scheme
    # stalls near 1e-2 of the initial residual on a 5 km domain.
    cfg['control'].update(gravity='0', has_PT='true', PT_option='1', PT_max_iter='100000',
                          PT_relative_tolerance='1e-7', PT_absolute_tolerance='1e-6',
                          is_quasi_static='true', inertial_scaling='1e11',
                          has_moving_mesh='false', fixed_dt=repr(DURATION/steps), has_thermal_diffusion='false')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'] = '1', '0'
    # Linear interpolation between entries: a repeated time makes each change a step.
    times, ratios = [0.0], [1.0]
    for start, rate in SCHEDULE[1:]:
        times += [start, start]; ratios += [ratios[-1], rate/Q]
    bc.update(vbc_z1='0', stress_bc_z1='0', stress_val_z1='0', hbc_x0='0', hbc_y0='0', hbc_z0='0', hbc_z1='0',
              hbc_x1='2', hbc_y1='2', hbc_pressure='[0,0,0,0,0,0]', num_wells='1', well_x='[0]', well_y='[0]',
              well_z_top='[0]', well_z_bottom=f'[{-B}]', well_rate=f'[{Q/4!r}]',
              num_well_schedule=str(len(times)),
              well_schedule_time_in_yr='[' + ', '.join(repr(x/YEAR2SEC) for x in times) + ']',
              well_schedule_ratio='[' + ', '.join(repr(r) for r in ratios) + ']')
    cfg['ic'].update(has_body_force_adjustment='false', pore_pressure_option='0')
    cfg['mat'].update(bulk_modulus=f"[{PROPS['bulk']!r}]", shear_modulus=f"[{PROPS['shear']!r}]",
                      biot_coeff=f"[{PROPS['biot']!r}]", porosity=f"[{PROPS['porosity']!r}]",
                      fluid_bulk_modulus=f"[{PROPS['fluid_bulk']!r}]", bulk_modulus_s=f"[{PROPS['grain_bulk']!r}]",
                      hydraulic_perm=f"[{PROPS['permeability']!r}]", fluid_visc=f"[{PROPS['viscosity']!r}]",
                      retention_model='[0]')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'basin')); rows = scan_frames(str(folder/'basin'))
    coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
    rn = np.hypot(coord[:, 0], coord[:, 1])
    centre = coord[conn].mean(axis=1); re = np.hypot(centre[:, 0], centre[:, 1])
    node = (rn >= BAND[0]) & (rn <= BAND[1]); elem = (re >= BAND[0]) & (re <= BAND[1])
    frames = []
    for row in rows[1:]:
        t = row['time']
        p = model.read_field(row['frame'], 'pore pressure')
        sub = -B*model.read_field(row['frame'], 'strain')[:, 2]
        rp, rs = reference(rn[node], t)['pressure'], reference(re[elem], t)['subsidence']
        wells = {}
        for r in WELLS:
            near = np.argsort(np.abs(re - r))[:8]               # the elements closest to the observation radius
            wells[f'{r:g}'] = dict(des=float(sub[near].mean()), ref=float(reference(re[near], t)['subsidence'].mean()))
        frames.append(dict(time=t, pressure_error=float(np.max(np.abs(p[node] - rp))/np.abs(rp).max()) if np.abs(rp).max() > 0 else 0.0,
                           subsidence_error=float(np.max(np.abs(sub[elem] - rs))/np.abs(rs).max()) if np.abs(rs).max() > 0 else 0.0,
                           peak_pressure_ref=float(np.abs(rp).max()), peak_subsidence_ref=float(np.abs(rs).max()),
                           wells=wells))
    return frames, len(coord)


def study(exe, out, spacings, steps):
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    result = dict(executable=str(exe), executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  steps=steps, duration_s=DURATION, schedule=SCHEDULE, properties=PROPS,
                  reference=dict(diffusivity=REF.diffusivity, specific_storage=REF.specific_storage,
                                 front_m=math.sqrt(4*REF.diffusivity*DURATION)), cases=[])
    (out/'source.diff').write_bytes(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=ROOT))
    for h in spacings:
        folder = out/f'h{h:g}'; folder.mkdir()
        with (folder/'basin.cfg').open('w') as f: config(h, steps).write(f)
        with (folder/'run.log').open('w') as f:
            subprocess.run([str(exe), 'basin.cfg'], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                           check=True, timeout=40000)
        frames, nodes = compare(folder)
        item = dict(spacing_m=h, nodes=nodes, frames=frames,
                    pressure_error=max(f['pressure_error'] for f in frames),
                    subsidence_error=max(f['subsidence_error'] for f in frames), input_sha256=digest(folder/'basin.cfg'))
        result['cases'].append(item)
        (out/'summary.json').write_text(json.dumps(result, indent=1, allow_nan=False))
        print(json.dumps({k: v for k, v in item.items() if k != 'frames'}), flush=True)
    assert digest(exe) == result['executable_sha256']
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--spacings', type=float, nargs='+', default=[100.0, 50.0, 25.0])
    p.add_argument('--steps', type=int, default=104)
    a = p.parse_args()
    study(a.executable.resolve(), a.output.resolve(), a.spacings, a.steps)
