#!/usr/bin/env python3
"""Check explicit and GVS fluid content and restart.

Usage: check_explicit.py EXECUTABLE OUTPUT_DIRECTORY

On a sealed fixed mesh with zero hydraulic mobility, the pressure storage
change must equal minus the projected volumetric strain.
"""

import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tests/functional'))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

from des_case import digest, read_config, reset_boundaries, write_json
from des_test_utils import checkpoint_arrays, fields, lumped_volumes, mesh_measures, run_des, simplex_volumes


def stress_coupling(executable, output, dim, method, rheology):
    """An isochoric pressure ramp and decaying shear select the viscous EVP branch."""
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    dt, viscosity, shear, biot, pressure = .01, 1e5, 1e7, .6, 4e4
    cfg['sim'].update(modelname='ramp', max_steps='2')
    cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
    cfg['control'].update(has_PT='false', gravity='0', has_thermal_diffusion='false',
                          fixed_dt=repr(dt), is_using_mixed_stress='false',
                          use_global_velocity_scaling=str(method == 'gvs'))
    initial = [5000., -5000.] + [0.] * (1 if dim == 2 else 4)
    cfg['ic'].update(has_body_force_adjustment='false', weakzone_option='0',
                     initial_stress_option='1', initial_stress=json.dumps(initial),
                     pore_pressure_option='2', pore_pressure_value='0')
    reset_boundaries(cfg['bc'], vbc='3', vbc_val='0', stress_bc='0', hbc='2')
    cfg['bc'].update(hbc_pressure=json.dumps([pressure]*6), num_hbc_schedule='2',
                     hbc_schedule_time_in_yr=json.dumps([0., 2*dt/YEAR2SEC]),
                     hbc_schedule_ratio='[0,1]')
    cfg['mat'].update(rheology_type=rheology, is_plane_strain=str(dim == 2),
                      min_viscosity=repr(viscosity), max_viscosity=repr(viscosity),
                      biot_coeff=f'[{biot}]', hydraulic_perm='[0]',
                      cohesion0='[1e7]', cohesion1='[1e7]')
    output.mkdir()
    run_des(executable, output, 'input', cfg, timeout=60)
    model = Dynearthsol(str(output/'ramp'))
    rows = scan_frames(model.modelname)
    assert model.ndims == dim and [row['steps'] for row in rows] == [0, 1, 2]
    relax = (1 - dt*shear/(2*viscosity))/(1 + dt*shear/(2*viscosity))
    errors = []
    for row in rows[1:]:
        frame, step = row['frame'], row['steps']
        p = pressure*step/2
        np.testing.assert_allclose(model.read_field(frame, 'pore pressure'), p, rtol=1e-12, atol=0)
        np.testing.assert_array_equal(model.read_field(frame, 'strain'), 0.)
        effective = np.zeros(len(initial)) if rheology == 'viscous' else np.array(initial)*relax**step
        total = effective.copy()
        total[:dim] -= biot*p
        s = model.read_field(frame, 'stress')
        np.testing.assert_allclose(s, np.broadcast_to(total, s.shape), rtol=0, atol=1e-10)
        errors.append(float(np.max(np.abs(s-total))))
        if dim == 2:
            yy = checkpoint_arrays(output/f'ramp.chkpt.{frame:06d}', ['stressyy'])['stressyy']
            np.testing.assert_allclose(yy, -biot*p, rtol=0, atol=1e-10)
    return max(errors)

def main():
    executable, output = (Path(p).resolve() for p in sys.argv[1:])
    output.mkdir(parents=True, exist_ok=False)
    executable_hash = digest(executable)
    summary = {}
    initial = 1000
    for method in ('explicit', 'gvs-shear', 'gvs-bulk'):
        folder = output / method
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
        max_content_error = 0.
        for row in rows:
            frame = row['frame']
            p = model.read_field(frame, 'pore pressure')
            eps = model.read_field(frame, 'strain')[:, :dim].sum(axis=1) - eps0
            content = storage*(p-p0) + lumped_volumes(x, conn, volume*eps)/nodal_volume
            assert np.min(p) > 0
            max_content_error = max(max_content_error, float(np.max(np.abs(content))))
            for field in ('velocity', 'stress', 'strain', 'pore pressure'):
                assert np.isfinite(model.read_field(frame, field)).all()
        assert max_content_error < 2e-13, (method, max_content_error)
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
                    assert a[key] == b[key], (method, key)
        summary[method] = dict(dimension=dim, frames=len(rows), exact_physical_restart=True,
                             max_local_content_error=max_content_error)
        print(method, 'PASS', flush=True)
    for method in ('explicit', 'gvs'):
        rheologies = ['viscous', 'maxwell', 'elasto-visco-plastic']
        if method == 'gvs':
            rheologies.append('elasto-visco-plastic-rsf')
        for rheology in rheologies:
            name = f'{method}-{rheology}'
            error = stress_coupling(executable, output/name, dim, method, rheology)
            summary[name] = dict(dimension=dim, max_stress_error_pa=error)
            print(name, 'PASS', flush=True)
    assert digest(executable) == executable_hash
    write_json(output / 'summary.json', dict(executable_sha256=executable_hash, cases=summary))


if __name__ == '__main__':
    main()
