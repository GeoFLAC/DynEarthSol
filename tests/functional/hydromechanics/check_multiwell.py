#!/usr/bin/env python3
"""Verify two-well superposition and independent sealed-volume closure.

Usage: check_multiwell.py EXECUTABLE NEW_OUTPUT_DIRECTORY
One well pumps and one injects, then both shut in and reverse their rates.
Three complete DES runs use identical fixed meshes, timesteps and schedules.
Rates keep strains below 1e-5: DES retains its objective stress rotation
even on a fixed mesh, so exact finite-amplitude superposition is not its
constitutive model. This checks the small-strain coupled/source limit, not
infinite-domain Theis accuracy or nonlinear/evolving-property superposition.
"""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from check_wells import ROOT, balance, config

sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

RATE_SCALE = 1e-3


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    exe, output = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    dim = 3 if exe.name.endswith('3d') else 2
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    original_hash = digest(exe)
    cases, records = {}, {}
    for label, wells in [('pump', [0]), ('inject', [1]), ('combined', [0, 1])]:
        folder = output / label
        folder.mkdir()
        cfg, rate = config(dim, False)
        rate *= RATE_SCALE
        cfg['mesh'] = dict(meshing_option='1', meshing_elem_shape='1', xlength='40',
                           ylength='40', zlength='10', resolution='10', remeshing_option='0')
        cfg['control'].update(fixed_dt='37', PT_absolute_tolerance='1e-8')
        cfg['sim'].pop('output_time_interval_in_yr', None)
        cfg['sim'].update(max_steps='600', max_time_in_yr=repr(600 * 37 / YEAR2SEC),
                          output_step_interval='100')
        rates = [rate, -0.35 * rate]
        cfg['bc'].update(num_wells=str(len(wells)),
                         well_x=json.dumps([10 if w == 0 else 30 for w in wells]),
                         well_y=json.dumps([20] * len(wells)),
                         well_z_top=json.dumps([0] * len(wells)),
                         well_z_bottom=json.dumps([-10] * len(wells)),
                         well_rate=json.dumps([rates[w] for w in wells]))
        with (folder / 'input.cfg').open('w') as stream:
            cfg.write(stream)
        records[label] = dict(executable_sha256=original_hash,
                              input_sha256=digest(folder / 'input.cfg'),
                              test_sha256=digest(Path(__file__)), rate_scale=RATE_SCALE,
                              net_rate=sum(rates[w] for w in wells), status='running')
        (output / 'provenance.json').write_text(json.dumps(records, indent=2) + '\n')
        print(f'Running {dim}D {label}', flush=True)
        with (folder / 'run.log').open('w') as stream:
            subprocess.run([str(exe), 'input.cfg'], cwd=folder, env=env,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=180)
        if digest(exe) != original_hash:
            raise RuntimeError('Executable changed during multi-well verification')
        model = Dynearthsol(str(folder / 'wells'))
        frames = scan_frames(str(folder / 'wells'))
        np.testing.assert_allclose([f['time'] for f in frames], np.arange(7) * 3700,
                                   rtol=0, atol=1e-8)
        residual, count = balance(folder, dim, records[label]['net_rate'])
        assert residual < 1e-9, (label, residual)
        strain = max(float(np.max(np.abs(model.read_field(frame, 'strain')))) for frame in model.frames)
        assert strain < 1e-5, ('Superposition fixture must stay in the small-strain limit', label, strain)
        for frame in model.frames:
            b = model.read_field(frame, 'fluid balance')
            assert abs(b[4]) <= 1e-9 * max(np.abs(b[:4]).max(), 1e-30) + 1e-20, (label, frame, b)
        cases[label] = model
        records[label].update(status='complete', max_content_balance_error_over_volume=residual,
                              frames=count, max_abs_strain=strain)
        (output / 'provenance.json').write_text(json.dumps(records, indent=2) + '\n')

    combined, pump, inject = (cases[name] for name in ('combined', 'pump', 'inject'))
    for frame in combined.frames:
        for name in ('coordinate', 'connectivity'):
            np.testing.assert_array_equal(combined.read_field(frame, name), pump.read_field(frame, name))
            np.testing.assert_array_equal(combined.read_field(frame, name), inject.read_field(frame, name))
    errors = {}
    for name in ('pore pressure', 'stress', 'strain', 'fluid balance'):
        a = np.asarray([pump.read_field(f, name) for f in pump.frames])
        b = np.asarray([inject.read_field(f, name) for f in inject.frames])
        c = np.asarray([combined.read_field(f, name) for f in combined.frames])
        scale = max(float(np.max(np.abs(a) + np.abs(b))), 1e-30)
        errors[name] = float(np.max(np.abs(c - (a + b)))) / scale
        assert errors[name] < 1e-8, (name, errors[name])
    summary = dict(dimension=dim, fields_max_superposition_error=errors, cases=records,
                   notes=['Small-strain elastic limit, fixed geometry/properties, sealed boundaries.',
                          'Schedule events fall between physical timesteps.',
                          'Injection-only closure uses a positive normalization volume.',
                          'This is not an analytical radial-flow accuracy test.'])
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
