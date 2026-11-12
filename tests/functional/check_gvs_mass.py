#!/usr/bin/env python3
"""Check fixed-mesh GVS mass from the independent velocity-update identity.

Usage: check_gvs_mass.py EXECUTABLE OUTPUT_DIRECTORY
Runs both wave-speed ceilings with a sudden dry traction, one saved frame
per step. Checks the mass actually used, rather than only saved dt or vmax.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, scan_frames
from pt_initial_checkpoint import fields
from des_test_utils import read_config, run_des


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    executable, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    executable_hash = hashlib.sha256(executable.read_bytes()).hexdigest()
    summary = {}
    for ceiling in ('shear', 'bulk'):
        folder = output / ceiling
        folder.mkdir()
        cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
        cfg['sim'].update(modelname='fresh', max_steps='64', max_time_in_yr='1e6',
                          output_step_interval='1', checkpoint_frame_interval='1',
                          has_initial_checkpoint='yes')
        cfg['sim'].pop('output_time_interval_in_yr', None)
        cfg['control'].update(has_PT='false', has_hydraulic_diffusion='false',
                              has_thermal_diffusion='false', gravity='0',
                              use_global_velocity_scaling='true',
                              mass_scaling_reference_speed=ceiling,
                              inertial_scaling='1e6', characteristic_speed='1e-10',
                              dt_fraction='.1', fixed_dt='0')
        cfg['bc'].update(vbc_x1='1', stress_val_z1='-1000')
        cfg['ic']['has_body_force_adjustment'] = 'false'
        cfg['mat'].update(bulk_modulus='[2e7]', shear_modulus='[1e7]', alpha='[0]')
        cfg['control']['has_moving_mesh'] = 'false'
        timeout = 60
        run_des(executable, folder, 'input', cfg, timeout)
        model = Dynearthsol(str(folder / 'fresh'))
        rows = scan_frames(str(folder / 'fresh'))
        times = {r['frame']: r['time'] for r in rows}
        assert len(model.frames) == 65
        x = model.read_field(0, 'coordinate')
        conn = model.read_field(0, 'connectivity')
        dim = model.ndims
        cells = x[conn]
        volume = np.abs(np.linalg.det(cells[:, 1:] - cells[:, :1])) / math.factorial(dim)
        lumped_volume = np.zeros(len(x))
        np.add.at(lumped_volume, conn, np.broadcast_to((volume / (dim + 1))[:, None], conn.shape))
        free = np.ones_like(x, dtype=bool)
        for d in range(dim - 1):
            free[:, d] = ~(np.isclose(x[:, d], 0) | np.isclose(x[:, d], 1))
        free[:, -1] = ~np.isclose(x[:, -1], -1)
        old_velocity = model.read_field(0, 'velocity')
        scale_min, scale_max, max_error = float('inf'), 0., 0.
        for previous, frame in zip(model.frames, model.frames[1:]):
            scale = max(1e-10, np.max(np.linalg.norm(old_velocity[conn].mean(axis=1), axis=1)))
            wave_speed = min(scale * 1e6, math.sqrt((1e7 if ceiling == 'shear' else 2e7) / 2650))
            mass = lumped_volume * 2e7 / wave_speed**2
            checkpoint = fields(folder / f'fresh.chkpt.{previous:06d}')
            # Saved dt is the next interval. Subtracting two large accumulated
            # clocks loses precision after GVS shortens the step sharply.
            dt = float(np.frombuffer(checkpoint['dt'], dtype=np.float64)[0])
            np.testing.assert_allclose(np.frombuffer(checkpoint['max_global_vel_mag'], dtype=np.float64)[0],
                                       scale, atol=0, rtol=2e-14)
            assert abs(times[frame] - times[previous] - dt) < 8*np.finfo(float).eps*max(1, abs(times[frame]))
            velocity = model.read_field(frame, 'velocity')
            force = model.read_field(frame, 'force')
            actual = velocity - old_velocity
            expected = dt * force / mass[:, None]
            scale_error = max(np.max(np.abs(actual[free])), np.max(np.abs(expected[free])), 1e-30)
            error = np.max(np.abs((actual - expected)[free])) / scale_error
            max_error = max(max_error, float(error))
            assert error < 2e-12, f'{ceiling} step {frame}: mass/velocity identity error {error}'
            scale_min, scale_max = min(scale_min, scale), max(scale_max, scale)
            old_velocity = velocity
        assert scale_max / scale_min > 100, 'GVS did not change its scale enough to exercise mass updates'
        cfg['sim'].update(modelname='resumed', is_restarting='true',
                          restarting_from_modelname='fresh', restarting_from_frame='16')
        run_des(executable, folder, 'restart', cfg, timeout)
        resumed = Dynearthsol(str(folder / 'resumed'))
        fresh_steps = {r['steps']: r for r in rows}
        for row in scan_frames(str(folder / 'resumed')):
            fresh = fresh_steps[row['steps']]
            assert row['time'] == fresh['time']
            names = ('coordinate','velocity','stress','strain','temperature')
            for field in names:
                np.testing.assert_array_equal(model.read_field(fresh['frame'], field),
                                              resumed.read_field(row['frame'], field))
        summary[ceiling] = dict(dimension=int(dim), frames=len(model.frames),
                                max_relative_velocity_identity_error=max_error,
                                velocity_scale_min=scale_min, velocity_scale_max=scale_max,
                                exact_restart=True)
    assert hashlib.sha256(executable.read_bytes()).hexdigest() == executable_hash
    (output / 'summary.json').write_text(json.dumps(dict(executable_sha256=executable_hash,
                                                        moving_mesh=False,
                                                        cases=summary), indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
