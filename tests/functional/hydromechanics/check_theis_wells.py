#!/usr/bin/env python3
"""Qualify several fully penetrating wells against Theis superposition.

Uses a full 3D confined slab, two extracting wells and an injecting well,
followed by recovery. Pressure feedback is disabled and material properties
are constant: this tests transport, not deformation or wellbore physics.
Independent volume integration and separate single-well runs check the
closed budget and discrete superposition as well as the analytical error.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

import check_theis as ct
from theis_reference import scheduled_drawdown

L, PUMP, END = 160., 40., 80.
WELLS = ((60., 80., ct.Q), (100., 80., .6*ct.Q), (80., 100., -.25*ct.Q))
SCALE = sum(abs(q) for _, _, q in WELLS)/(4*math.pi*ct.T)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference(xy, time, wells):
    result = np.zeros(len(xy))
    for x, y, q in wells:
        radius = np.linalg.norm(xy-np.array([x, y]), axis=1)
        result += scheduled_drawdown(radius, np.full(len(xy), time),
                                     [(0., q), (PUMP, 0.)], ct.T, ct.S)
    return result


def config(spacing, dt, wells):
    cfg = ct.config(spacing)
    steps = int(round(END/dt))
    cfg['sim'].update(modelname='multiwell', max_steps=str(steps),
                      output_step_interval=str(steps//8), checkpoint_frame_interval='1')
    cfg['mesh'].update(xlength=str(L), ylength=str(L))
    cfg['control']['fixed_dt'] = str(dt)
    cfg['bc'].update(num_wells=str(len(wells)),
                     well_x=str([x for x, _, _ in wells]),
                     well_y=str([y for _, y, _ in wells]),
                     well_rate=str([q for _, _, q in wells]),
                     well_z_top=str([0.]*len(wells)),
                     well_z_bottom=str([-ct.B]*len(wells)),
                     well_schedule_time_in_yr=str([0., PUMP/ct.YEAR2SEC, PUMP/ct.YEAR2SEC]))
    return cfg


def run(executable, folder, spacing, dt, wells):
    folder.mkdir()
    with (folder/'input.cfg').open('w') as stream:
        config(spacing, dt, wells).write(stream)
    with (folder/'run.log').open('w') as stream:
        subprocess.run([str(executable), 'input.cfg'], cwd=folder, check=True,
                       stdout=stream, stderr=subprocess.STDOUT, timeout=600,
                       env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
    model = ct.Dynearthsol(str(folder/'multiwell'))
    rows = ct.scan_frames(str(folder/'multiwell'))
    assert len(rows) == 9
    x = model.read_field(0, 'coordinate')
    assert model.ndims == 3
    conn = model.read_field(0, 'connectivity')
    cells = x[conn]
    volume = np.abs(np.linalg.det(cells[:, 1:]-cells[:, :1]))/6
    p0 = model.read_field(0, 'pore pressure')
    radii = np.stack([np.hypot(x[:, 0]-wx, x[:, 1]-wy) for wx, wy, _ in wells])
    band = (radii.min(axis=0) >= 10) & (radii.min(axis=0) <= 30)
    assert band.any()
    boundary = ((x[:, 0] == 0) | (x[:, 0] == L) | (x[:, 1] == 0) | (x[:, 1] == L))
    frames = []
    for row in rows[1:]:
        frame, time = row['frame'], row['time']
        p = model.read_field(frame, 'pore pressure')
        np.testing.assert_array_equal(model.read_field(frame, 'coordinate'), x)
        assert np.isfinite(p).all()
        exact = reference(x[band, :2], time, wells)
        error = (p0[band]-p[band]-exact)/SCALE
        stored = ct.STORAGE*float(volume @ np.mean((p-p0)[conn], axis=1))
        expected = -sum(q for _, _, q in wells)*min(time, PUMP)
        boundary_tail = reference(x[boundary, :2], time, wells)
        frames.append(dict(time_s=time, max_error_over_scale=float(np.abs(error).max()),
                           rms_error_over_scale=float(np.sqrt(np.mean(error**2))),
                           stored_m3=stored, exact_stored_m3=expected,
                           absolute_volume_error_m3=abs(stored-expected),
                           far_reference_over_scale=float(np.max(np.abs(boundary_tail)))/SCALE))
    result = dict(spacing_m=spacing, dt_s=dt, wells=wells, nodes=len(x), elements=len(conn),
                  frames=frames, max_error_over_scale=max(f['max_error_over_scale'] for f in frames),
                  max_volume_error_m3=max(f['absolute_volume_error_m3'] for f in frames),
                  outputs_sha256={p.name: digest(p) for p in sorted(folder.glob('multiwell.*'))
                                  if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.info'},
                  input_sha256=digest(folder/'input.cfg'))
    assert result['max_volume_error_m3'] < sum(abs(q) for _, _, q in wells)*PUMP*1e-10
    assert max(f['far_reference_over_scale'] for f in frames) < 1e-5
    print(folder.name, 'pressure/scale', result['max_error_over_scale'],
          'volume error m3', result['max_volume_error_m3'], flush=True)
    return result


def study(executable, output, spacings, manifest):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = digest(executable)
    solver = None
    if manifest:
        solver = json.loads(manifest.read_text())
        assert binary_hash in solver['binaries'].values()
    summary = dict(executable=str(executable), executable_sha256=binary_hash,
                   solver_commit=solver['implementation_commit'] if solver else None,
                   binary_manifest_sha256=digest(manifest) if manifest else None,
                   analyzer_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ct.ROOT, text=True).strip(),
                   runner_sha256=digest(Path(__file__)),
                   reference=dict(box_side_m=L, thickness_m=ct.B, pump_until_s=PUMP,
                                  end_s=END, pressure_scale_Pa=SCALE, storage_per_Pa=ct.STORAGE,
                                  transmissivity=ct.T, storativity=ct.S), cases={})
    for source in (Path(__file__), Path(ct.__file__), ct.ROOT/'benchmarks/poromechanics/theis_reference.py'):
        shutil.copy2(source, output/source.name)
    if manifest:
        shutil.copy2(manifest, output/'solver-qualification.json')

    def save():
        (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')

    for spacing in spacings:
        key = f'joint-h{spacing:g}-dt1'
        summary['cases'][key] = run(executable, output/key, spacing, 1., WELLS)
        save()
    finest = min(spacings)
    key = f'joint-h{finest:g}-dt0.125'
    summary['cases'][key] = run(executable, output/key, finest, .125, WELLS)
    joint = ct.Dynearthsol(str(output/key/'multiwell'))
    individual = []
    for i, well in enumerate(WELLS):
        name = f'well{i+1}-h{finest:g}'
        summary['cases'][name] = run(executable, output/name, finest, .125, [well])
        individual.append(ct.Dynearthsol(str(output/name/'multiwell')))
        save()
    worst = 0.
    for frame in joint.frames:
        expected = sum(m.read_field(frame, 'pore pressure') for m in individual)
        actual = joint.read_field(frame, 'pore pressure')
        worst = max(worst, float(np.max(np.abs(expected-actual)))/SCALE)
    assert worst < 2e-12, worst
    summary.update(max_discrete_superposition_error_over_scale=worst,
                   analytical_error_below_one_percent=summary['cases'][key]['max_error_over_scale'] < .01,
                   mesh_error_decreases=summary['cases'][f'joint-h{finest:g}-dt1']['max_error_over_scale'] <
                       summary['cases'][f'joint-h{max(spacings):g}-dt1']['max_error_over_scale'])
    summary['passed'] = summary['analytical_error_below_one_percent'] and summary['mesh_error_decreases']
    assert digest(executable) == binary_hash
    save()
    assert summary['passed'], 'Analytical accuracy/convergence target not met; keep the raw evidence'
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--spacings', type=float, nargs='+', default=[5., 2.5])
    parser.add_argument('--binary-manifest', type=Path)
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.spacings,
          args.binary_manifest.resolve() if args.binary_manifest else None)
