#!/usr/bin/env python3
"""Pump a thin 3D DES slab and compare it with the averaged-aquifer reference.

Usage: pumping_slab.py OUTPUT_DIRECTORY [--analyze-only]
A quarter slab (symmetry rollers at x=0 and y=0, roller base, free top) has a
full-height well on its corner edge. While the pressure front stays well inside
the domain, the slab approaches the plane-stress model of
averaged_aquifer_reference.py, up to O((B/r)^2). Nodal pressure and element
areal/vertical strains are compared with AveragedAquifer.line_well.
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

from averaged_aquifer_reference import AveragedAquifer

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

L, B, Q = 400.0, 10.0, 0.01               # m, m, m^3/s for the full well
PROPS = dict(thickness=B, bulk=2e8, shear=1e8, biot=0.8, porosity=0.25,
             fluid_bulk=2.2e9, grain_bulk=3.7e10, permeability=1e-12, viscosity=1e-3)
BAND = (15.0, 120.0)                      # compared radii, m
DURATION = 1.0e4                          # s: sqrt(4 D t) stays below 0.3 L


def config(resolution=20.0):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(ROOT / 'benchmarks/poromechanics/1d-consolidation-des3d_traction.cfg')
    cfg['sim'].update(modelname='slab', max_time_in_yr=repr(DURATION / YEAR2SEC),
                      output_time_interval_in_yr=repr(DURATION / 4 / YEAR2SEC))
    cfg['mesh'] = dict(meshing_option='1', xlength=repr(L), ylength=repr(L), zlength=repr(B),
                       resolution=repr(resolution), remeshing_option='0')
    cfg['control'].update(gravity='0', PT_option='1', PT_absolute_tolerance='1e-6')
    bc = cfg['bc']
    for side in ('x0', 'x1', 'y0', 'y1', 'z0'):
        bc[f'vbc_{side}'], bc[f'vbc_val_{side}'] = '1', '0'
    bc.update(vbc_z1='0', stress_bc_z1='0', stress_val_z1='0', hbc_x0='0', hbc_y0='0',
              hbc_z0='0', hbc_z1='0', hbc_x1='1', hbc_y1='1', num_wells='1',
              well_x='[0]', well_y='[0]', well_z_top='[0]', well_z_bottom=f'[{-B}]',
              well_rate=f'[{Q / 4}]')
    cfg['mat'].update(bulk_modulus=f"[{PROPS['bulk']}]", shear_modulus=f"[{PROPS['shear']}]",
                      biot_coeff=f"[{PROPS['biot']}]", porosity=f"[{PROPS['porosity']}]",
                      fluid_bulk_modulus=f"[{PROPS['fluid_bulk']}]",
                      bulk_modulus_s=f"[{PROPS['grain_bulk']}]",
                      hydraulic_perm=f"[{PROPS['permeability']}]",
                      fluid_visc=f"[{PROPS['viscosity']}]")
    return cfg


def compare(case):
    model = Dynearthsol(str(case / 'slab'))
    times = {r['frame']: r['time'] for r in scan_frames(str(case / 'slab'))}
    aq = AveragedAquifer(**PROPS)
    xyz = model.read_field(model.frames[0], 'coordinate')
    conn = model.read_field(model.frames[0], 'connectivity')
    r_node = np.hypot(xyz[:, 0], xyz[:, 1])
    centre = xyz[conn].mean(axis=1)
    r_elem = np.hypot(centre[:, 0], centre[:, 1])
    node_in = (r_node >= BAND[0]) & (r_node <= BAND[1])
    elem_in = (r_elem >= BAND[0]) & (r_elem <= BAND[1])
    rows = []
    for frame in model.frames[1:]:
        t = times[frame]
        p = model.read_field(frame, 'pore pressure')
        eps = model.read_field(frame, 'strain')
        ref_p = aq.line_well(r_node[node_in], t, Q)['pressure']
        ref_e = aq.line_well(r_elem[elem_in], t, Q)
        areal, vertical = eps[elem_in, 0] + eps[elem_in, 1], eps[elem_in, 2]
        rows.append(dict(
            time_s=t, front_over_L=math.sqrt(4 * aq.diffusivity * t) / L,
            pressure_error=float(np.max(np.abs(p[node_in] - ref_p)) / np.max(np.abs(ref_p))),
            pressure_rms=float(np.sqrt(np.mean((p[node_in] - ref_p)**2)) / np.max(np.abs(ref_p))),
            areal_rms=float(np.sqrt(np.mean((areal - ref_e['areal_strain'])**2)) /
                            np.max(np.abs(ref_e['areal_strain']))),
            vertical_rms=float(np.sqrt(np.mean((vertical - ref_e['vertical_strain'])**2)) /
                               np.max(np.abs(ref_e['vertical_strain']))),
            plot=dict(r_node=r_node[node_in], p=p[node_in], r_elem=r_elem[elem_in],
                      areal=areal, vertical=vertical)))
    return aq, rows




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--analyze-only', '--plot-only', dest='plot_only', action='store_true', help='Check saved runs')
    parser.add_argument('--resolution', type=float, default=20.0, help='Uniform mesh spacing in m')
    args = parser.parse_args()
    if not math.isfinite(args.resolution) or args.resolution <= 0:
        parser.error('--resolution must be finite and positive')
    output = args.output.resolve()
    case = output / 'slab-3d'
    if not args.plot_only:
        output.mkdir(parents=True, exist_ok=False)
        case.mkdir()
        with (case / 'input.cfg').open('w') as stream:
            config(args.resolution).write(stream)
        provenance = dict(executable=str(ROOT/'dynearthsol3d'),
                          executable_sha256=hashlib.sha256((ROOT/'dynearthsol3d').read_bytes()).hexdigest(),
                          config_sha256=hashlib.sha256((case/'input.cfg').read_bytes()).hexdigest(),
                          resolution_m=args.resolution, omp_threads=1)
        (output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
        with (case / 'run.log').open('w') as stream:
            subprocess.run([str(ROOT / 'dynearthsol3d'), 'input.cfg'], cwd=case, check=True,
                           env=dict(os.environ, OMP_NUM_THREADS='1'),
                           stdout=stream, stderr=subprocess.STDOUT, timeout=3600)
    aq, rows = compare(case)
    summary = [{k: v for k, v in row.items() if k != 'plot'} for row in rows]
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    # This qualifies a coarse thin-slab comparison, not spatial/time convergence.
    limits = dict(pressure_rms=0.03, areal_rms=0.08, vertical_rms=0.08, front_over_L=0.3)
    failed = {key: max(row[key] for row in summary) for key, limit in limits.items()
              if max(row[key] for row in summary) > limit}
    if failed:
        raise SystemExit(f'Thin-slab acceptance limits exceeded: {failed}')


if __name__ == '__main__':
    main()
