#!/usr/bin/env python3
"""Check explicit/GVS fluid content, wet-law equivalence and restart.

Usage: check_explicit.py EXECUTABLE OUTPUT_DIRECTORY
The sealed fixed mesh has no hydraulic flux or thermal source. Its nodal
fluid-content change must equal minus the projected volumetric strain.
"""

import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tests/functional'))
sys.path.insert(0, str(ROOT / 'benchmarks/poromechanics'))
from Dynearthsol import Dynearthsol, scan_frames
from pt_initial_checkpoint import fields
from unsaturated_reference import Retention


def run(executable, folder, cfg, name):
    with (folder / f'{name}.cfg').open('w') as stream:
        cfg.write(stream)
    with (folder / f'{name}.log').open('w') as stream:
        subprocess.run([str(executable), f'{name}.cfg'], cwd=folder, check=True,
                       stdout=stream, stderr=subprocess.STDOUT, timeout=60,
                       env=dict(os.environ, OMP_NUM_THREADS='1',
                                OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))


def main():
    executable, output = (Path(p).resolve() for p in sys.argv[1:])
    output.mkdir(parents=True, exist_ok=False)
    executable_hash = hashlib.sha256(executable.read_bytes()).hexdigest()
    summary = {}
    for method in ('explicit', 'gvs-shear', 'gvs-bulk'):
        baseline = None
        for law, initial in [(0, 1000), (1, 1000), (2, 1000), (3, 1000), (1, -1e4)]:
            name = f'{method}-law{law}-p{initial:g}'
            folder = output / name
            folder.mkdir()
            cfg = configparser.ConfigParser()
            cfg.optionxform = str
            cfg.read(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg['sim'].update(modelname='fresh', max_steps='64', max_time_in_yr='1',
                              output_step_interval='1', checkpoint_frame_interval='1')
            cfg['control'].update(has_PT='false', is_quasi_static='true',
                                  has_thermal_diffusion='false', gravity='0',
                                  has_moving_mesh='false', fixed_dt='0', dt_fraction='.5',
                                  inertial_scaling='1e5', characteristic_speed='1e-5',
                                  use_global_velocity_scaling=str(method != 'explicit').lower(),
                                  mass_scaling_reference_speed=method.split('-')[-1]
                                      if method != 'explicit' else 'shear')
            cfg['bc'].update(vbc_x1='1', stress_val_z1='-1000')
            for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
                cfg['bc'][f'hbc_{side}'] = '0'
            cfg['ic'].update(has_body_force_adjustment='false',
                             pore_pressure_value=str(initial), pore_pressure_option='2')
            cfg['mat'].update(bulk_modulus='[2e7]', shear_modulus='[1e7]',
                              alpha='[0]', hydraulic_perm='[0]', retention_model=f'[{law}]',
                              retention_alpha='[2e-4]', retention_n='[1.6]',
                              residual_saturation='[0.1]')
            run(executable, folder, cfg, 'fresh')
            model = Dynearthsol(str(folder / 'fresh'))
            rows = scan_frames(str(folder / 'fresh'))
            assert len(rows) == 65
            x = model.read_field(0, 'coordinate')
            conn = model.read_field(0, 'connectivity')
            dim = model.ndims
            volume = np.abs(np.linalg.det(x[conn][:, 1:] - x[conn][:, :1])) / math.factorial(dim)
            measures = []
            for omit in range(dim+1):
                face = np.delete(x[conn], omit, axis=1)
                measures.append(np.linalg.norm(face[:, 1]-face[:, 0], axis=1) if dim == 2 else
                                np.linalg.norm(np.cross(face[:, 1]-face[:, 0],
                                                       face[:, 2]-face[:, 0]), axis=1)/2)
            height = np.min(dim*volume/np.max(measures, axis=0))
            nodal_volume = np.zeros(len(x))
            np.add.at(nodal_volume, conn, np.broadcast_to((volume/(dim+1))[:, None], conn.shape))
            p0 = model.read_field(0, 'pore pressure')
            np.testing.assert_array_equal(p0, initial)
            eps0 = model.read_field(0, 'strain')[:, :dim].sum(axis=1)
            storage = .3/2.17e9 + .7/37e9
            modulus = 2e7 + 4e7/3 + 1/storage
            expected_dt = .5*height*math.sqrt(2e7/modulus)/5
            np.testing.assert_allclose(rows[1]['time'], expected_dt, atol=0, rtol=2e-14)
            retention = Retention('van_genuchten', 2e-4, 1.6, .1)
            max_content_error = 0.
            for row in rows:
                frame = row['frame']
                p = model.read_field(frame, 'pore pressure')
                eps = model.read_field(frame, 'strain')[:, :dim].sum(axis=1) - eps0
                projection = np.zeros(len(x))
                np.add.at(projection, conn, np.broadcast_to((volume*eps/(dim+1))[:, None], conn.shape))
                content = storage*(p-p0) + projection/nodal_volume
                if initial < 0:
                    content += .3*(retention.sat(p)-retention.sat(p0))
                else:
                    assert np.min(p) > 0, 'Wet-law equivalence left the saturated range'
                max_content_error = max(max_content_error, float(np.max(np.abs(content))))
                for field in ('velocity', 'stress', 'strain', 'pore pressure'):
                    assert np.isfinite(model.read_field(frame, field)).all()
                if baseline is not None and initial > 0:
                    for suffix in ('save', 'chkpt'):
                        a = fields(baseline / f'fresh.{suffix}.{frame:06d}')
                        b = fields(folder / f'fresh.{suffix}.{frame:06d}')
                        assert a.keys() == b.keys()
                        assert all(a[k] == b[k] for k in a if k not in ('provenance', 'walltime_sec'))
            assert max_content_error < 2e-13, (name, max_content_error)
            if law == 0:
                baseline = folder
            cfg['sim'].update(modelname='resumed', is_restarting='true',
                              restarting_from_modelname='fresh', restarting_from_frame='16')
            run(executable, folder, cfg, 'resumed')
            resumed = Dynearthsol(str(folder / 'resumed'))
            fresh_by_step = {r['steps']: r for r in rows}
            for row in scan_frames(str(folder / 'resumed')):
                original = fresh_by_step[row['steps']]
                assert row['time'] == original['time']
                for field in ('coordinate', 'velocity', 'stress', 'strain', 'pore pressure', 'temperature'):
                    np.testing.assert_array_equal(model.read_field(original['frame'], field),
                                                  resumed.read_field(row['frame'], field))
                if row['steps'] > 16:
                    a = fields(folder / f'fresh.chkpt.{original["frame"]:06d}')
                    b = fields(folder / f'resumed.chkpt.{row["frame"]:06d}')
                    for key in ('dt', 'pending pore pressure increment', 'physical hydraulic step',
                                'fluid balance', 'max_global_vel_mag'):
                        assert a[key] == b[key], (name, key)
            summary[name] = dict(dimension=dim, frames=len(rows), exact_physical_restart=True,
                                 saturated_payloads_exact=initial > 0,
                                 max_local_content_error=max_content_error)
            (output / 'summary.json').write_text(json.dumps(dict(
                executable_sha256=executable_hash, cases=summary), indent=2) + '\n')
            print(name, 'PASS', flush=True)
    assert hashlib.sha256(executable.read_bytes()).hexdigest() == executable_hash


if __name__ == '__main__':
    main()
