#!/usr/bin/env python3
"""Check changed restart controls against independent timestep/mass identities.

The tiny elastic cell tests exact continuation with unchanged controls and the
first explicit interval after changing mechanical method, inertia, scaling or
fixed dt. --legacy-executable adds restarts from the PT checkpoints it writes.
"""

import argparse
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, scan_frames
from des_test_utils import compare_continuation, digest, fields, mesh_measures, read_config, roller_free, run_des

BULK, SHEAR, DENSITY = 2e7, 1e7, 2650.


def config():
    cfg = read_config(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(max_time_in_yr='1e6', output_step_interval='1',
                      output_time_interval_in_yr='1e6', checkpoint_frame_interval='1',
                      has_initial_checkpoint='true', is_outputting_averaged_fields='false')
    cfg['control'].update(gravity='0', has_thermal_diffusion='false',
                          has_hydraulic_diffusion='false', has_moving_mesh='false',
                          use_global_velocity_scaling='false', inertial_scaling='1e5',
                          characteristic_speed='1e-5', dt_fraction='.2', fixed_dt='1',
                          mass_scaling_reference_speed='shear', is_using_mixed_stress='false')
    cfg['bc'].update(vbc_x1='1', stress_val_z1='-1000')
    cfg['ic'].update(has_body_force_adjustment='false', weakzone_option='0')
    cfg['mat'].update(bulk_modulus=f'[{BULK}]', shear_modulus=f'[{SHEAR}]',
                      rho0=f'[{DENSITY}]', alpha='[0]')
    return cfg


def run(executable, folder, cfg, name, steps, source=None):
    cfg = copy.deepcopy(cfg)
    cfg['sim'].update(modelname=name, max_steps=str(steps),
                      is_restarting=str(source is not None).lower())
    if source is not None:
        cfg['sim'].update(restarting_from_modelname=str(source), restarting_from_frame='1')
    run_des(executable, folder, name, cfg)
    return Dynearthsol(str(folder/name))


def expected_dt(cfg, height):
    c = cfg['control']
    fixed = c.getfloat('fixed_dt')
    if fixed:
        return fixed
    speed = (c.getfloat('characteristic_speed')*c.getfloat('inertial_scaling')
             if c.getboolean('is_quasi_static') else math.sqrt(BULK/DENSITY))
    return .5*height/speed*c.getfloat('dt_fraction')


def check_interval(folder, cfg, model):
    rows = scan_frames(str(folder/'resumed'))
    x, height, lumped = mesh_measures(model)
    expected = expected_dt(cfg,height)
    checkpoint = fields(folder/f'resumed.chkpt.{model.frames[1]:06d}')
    saved = float(np.frombuffer(checkpoint['dt'],dtype=np.float64)[0])
    actual = rows[1]['time']-rows[0]['time']
    result = dict(expected_dt=expected,first_interval=actual,saved_next_dt=saved,
                  height=height,dimension=int(model.ndims))
    np.testing.assert_allclose([actual,saved],expected,rtol=2e-13,atol=2e-15)
    if not cfg['control'].getboolean('has_PT'):
        c = cfg['control']
        rho = BULK/(c.getfloat('characteristic_speed')*c.getfloat('inertial_scaling'))**2 \
              if c.getboolean('is_quasi_static') else DENSITY
        free = roller_free(x)
        before = model.read_field(model.frames[0],'velocity')
        after = model.read_field(model.frames[1],'velocity')
        force = model.read_field(model.frames[1],'force')
        expected_increment = expected*force/(rho*lumped[:,None])
        error = float(np.max(np.abs((after-before-expected_increment)[free])))
        scale = max(float(np.max(np.abs(before))),float(np.max(np.abs(after))),1e-30)
        assert error <= 64*np.finfo(float).eps*scale+2e-12*float(np.max(np.abs(expected_increment)))
        result['velocity_mass_identity_absolute_error'] = error
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--legacy-executable',type=Path)
    args = parser.parse_args()
    executable, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    report = dict(executable_sha256=digest(executable),cases=[],passed=False)
    def save():
        (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    sources = {}
    for method in ('pt','explicit'):
        folder = output/(method+'-unchanged');folder.mkdir()
        cfg = config();cfg['control']['has_PT'] = str(method=='pt').lower()
        if method == 'explicit':
            cfg['control']['fixed_dt'] = '0'
        run(executable,folder,cfg,'fresh',4)
        resumed = run(executable,folder,cfg,'resumed',4,folder/'fresh')
        payloads = compare_continuation(folder,'fresh',resumed)
        report['cases'].append(dict(name=method+'-unchanged',exact_payloads=payloads))
        sources[method] = folder/'fresh',cfg
        save()
    changes = [('pt-to-explicit','pt',dict(has_PT='false',fixed_dt='0')),
               ('pt-to-physical','pt',dict(has_PT='false',is_quasi_static='false',fixed_dt='0')),
               ('scaling','explicit',dict(inertial_scaling='1e6')),
               ('velocity-scale','explicit',dict(characteristic_speed='2e-5')),
               ('fraction','explicit',dict(dt_fraction='.1')),
               ('explicit-fixed','explicit',dict(fixed_dt='.0125')),
               ('pt-fixed','pt',dict(fixed_dt='.5')),
               ('pt-fixed-increase','pt',dict(fixed_dt='2'))]
    if args.legacy_executable:
        folder = output/'legacy-source';folder.mkdir()
        cfg = config()
        run(args.legacy_executable.resolve(),folder,cfg,'fresh',1)
        sources['legacy'] = folder/'fresh',cfg
        changes += [('legacy-to-explicit','legacy',dict(has_PT='false',fixed_dt='0')),
                    ('legacy-to-physical','legacy',dict(has_PT='false',is_quasi_static='false',fixed_dt='0'))]
    for name,kind,change in changes:
        folder = output/name;folder.mkdir()
        source,original = sources[kind]
        cfg = copy.deepcopy(original);cfg['control'].update(change)
        model = run(executable,folder,cfg,'resumed',3,source)
        try:
            result = check_interval(folder,cfg,model)
        except AssertionError:
            report['failure'] = name;save();raise
        report['cases'].append(dict(name=name,**result));save()
    assert digest(executable) == report['executable_sha256']
    report['passed'] = True;save()
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
