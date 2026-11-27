#!/usr/bin/env python3
"""Check explicit/GVS fluid content, wet-law equivalence and restart.

Usage: check_explicit.py EXECUTABLE OUTPUT_DIRECTORY
The sealed fixed mesh has no hydraulic flux or thermal source. Its nodal
fluid-content change must equal minus the projected volumetric strain,
weighted by the preceding element Bishop coefficient when unsaturated.
"""

import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tests/functional'))
sys.path.insert(0, str(ROOT / 'benchmarks/poromechanics'))
from Dynearthsol import Dynearthsol, scan_frames
from unsaturated_reference import Retention


from des_case import digest, read_config, reset_boundaries, retention_options, write_json
from des_test_utils import fields, lumped_volumes, mesh_measures, record_name, run_des, runtime_key, simplex_volumes

RETENTION = Retention('van_genuchten', 2e-4, 1.6, .1)


def main():
    executable, output = (Path(p).resolve() for p in sys.argv[1:])
    output.mkdir(parents=True, exist_ok=False)
    executable_hash = digest(executable)
    summary = {}
    for method in ('explicit', 'gvs-shear', 'gvs-bulk'):
        baseline = None
        for law, initial in [(0, 1000), (1, 1000), (2, 1000), (3, 1000), (1, -1e4)]:
            name = f'{method}-law{law}-p{initial:g}'
            folder = output / name
            folder.mkdir()
            cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
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
            reset_boundaries(cfg['bc'], hbc='0')
            cfg['ic'].update(has_body_force_adjustment='false',
                             pore_pressure_value=str(initial), pore_pressure_option='2')
            cfg['mat'].update(bulk_modulus='[2e7]', shear_modulus='[1e7]', alpha='[0]', hydraulic_perm='[0]')
            # Every law takes the parameters of RETENTION, the reference of the unsaturated case.
            cfg['mat'].update(retention_options(RETENTION), retention_model=f'[{law}]')
            run_des(executable, folder, 'fresh', cfg, timeout=60)
            model = Dynearthsol(str(folder / 'fresh'))
            rows = scan_frames(str(folder / 'fresh'))
            assert len(rows) == 65
            x, height, nodal_volume = mesh_measures(model)
            conn = model.read_field(0, 'connectivity')
            dim = int(model.ndims)
            volume = simplex_volumes(x, conn)
            p0 = model.read_field(0, 'pore pressure')
            np.testing.assert_array_equal(p0, initial)
            eps0 = model.read_field(0, 'strain')[:, :dim].sum(axis=1)
            storage = .3/2.17e9 + .7/37e9
            modulus = 2e7 + 4e7/3 + 1/storage
            expected_dt = .5*height*math.sqrt(2e7/modulus)/5
            np.testing.assert_allclose(rows[1]['time'], expected_dt, atol=0, rtol=2e-14)
            previous_p, previous_eps = p0, np.zeros(len(conn))
            mechanical_content = np.zeros(len(x))
            max_content_error = 0.
            for row in rows:
                frame = row['frame']
                p = model.read_field(frame, 'pore pressure')
                eps = model.read_field(frame, 'strain')[:, :dim].sum(axis=1) - eps0
                content = storage*(p-p0) + lumped_volumes(x, conn, volume*eps)/nodal_volume
                if initial < 0:
                    # The rate law contains alpha*chi*d(eps_v), not alpha*d(eps_v)
                    # or d(alpha*chi*eps_v). This zero-mobility fixture takes one
                    # hydraulic interval per saved mechanical step.
                    chi = RETENTION.se(previous_p[conn]).mean(axis=1)
                    change = volume*chi*(eps-previous_eps)/(dim+1)
                    np.add.at(mechanical_content, conn, np.broadcast_to(change[:, None], conn.shape))
                    content = storage*(p-p0) + mechanical_content/nodal_volume
                    content += .3*(RETENTION.sat(p)-RETENTION.sat(p0))
                    previous_p, previous_eps = p, eps
                else:
                    assert np.min(p) > 0, 'Wet-law equivalence left the saturated range'
                max_content_error = max(max_content_error, float(np.max(np.abs(content))))
                for field in ('velocity', 'stress', 'strain', 'pore pressure'):
                    assert np.isfinite(model.read_field(frame, field)).all()
                if baseline is not None and initial > 0:
                    for suffix in ('save', 'chkpt'):
                        a = fields(baseline / f'fresh.{suffix}.{frame:06d}')
                        b = fields(folder / f'fresh.{suffix}.{frame:06d}')
                        # A retentive law adds four diagnostics even on a wet branch;
                        # every other field must equal the saturated baseline byte for byte.
                        diagnostics = {'saturation', 'effective saturation',
                                       'relative permeability', 'bulk density'}
                        assert a.keys() <= b.keys()
                        assert all(record_name(k) in diagnostics for k in b.keys()-a.keys())
                        assert all(a[k] == b[k] for k in a if not runtime_key(k))
                    for field in ('saturation', 'effective saturation', 'relative permeability'):
                        np.testing.assert_array_equal(model.read_field(frame, field), 1.0)
            assert max_content_error < 2e-13, (name, max_content_error)
            if law == 0:
                baseline = folder
            cfg['sim'].update(modelname='resumed', is_restarting='true',
                              restarting_from_modelname='fresh', restarting_from_frame='16')
            run_des(executable, folder, 'resumed', cfg, timeout=60)
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
            write_json(output / 'summary.json', dict(executable_sha256=executable_hash, cases=summary))
            print(name, 'PASS', flush=True)
    assert digest(executable) == executable_hash


if __name__ == '__main__':
    main()
