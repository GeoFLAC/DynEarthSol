#!/usr/bin/env python3
"""Check linear Biot fluid content across real remeshes and exact continuation.

Fixed saturated material properties only. Across each remesh the mobile and
matrix contents change only by what the domain change removes and adds, as
DES reports it; a remesh that is not meant to change the domain changes its
volume by at most 1e-6, and the moving-domain controls, a rigid translation,
add or remove the uniform excess content of the swept volume. This does not
check total retained mass or energy.
"""
import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames


from des_case import read_config, reset_boundaries, run_des, write_json
from des_test_utils import checkpoint_file, fields, physical_differences, runtime_key, simplex_volumes

def configuration(case):
    cfg = read_config(ROOT / 'tests/functional/pt_remesh_restart.cfg')
    cfg['sim'].update(modelname='flow', max_steps='3', has_marker_output='true',
                      has_output_during_remeshing='true')
    cfg['control']['PT_absolute_tolerance'] = '1e-8'
    if case == 'sealed':
        cfg['bc']['hbc_z1'] = '0'
    if case in ('explicit', 'gvs'):
        cfg['control'].update(has_PT='false', fixed_dt='.001', inertial_scaling='1e5',
                              use_global_velocity_scaling=str(case == 'gvs'))
    if case == 'matrix':
        cfg['mat'].update(matrix_storage='[1e-9]', matrix_exchange='[1e-12]')
    if case in ('uncoupled', 'uncoupled-thermal'):
        cfg['control']['has_poroelastic_pressure_feedback'] = 'false'
    if case in ('thermal', 'uncoupled-thermal'):
        cfg['control']['has_thermal_pressurization'] = 'true'
        cfg['mat'].update(alpha='[3e-5]', fluid_alpha='[2e-4]')
        cfg['ic'].update(temperature_option='4', temp_background='273',
                         temp_anomaly_amplitude='10', temp_anomaly_xcenter='.5',
                         temp_anomaly_ycenter='.5', temp_anomaly_zcenter='-5',
                         temp_anomaly_xsemi_axis='1', temp_anomaly_ysemi_axis='1',
                         temp_anomaly_zsemi_axis='3')
    if case.startswith('domain'):
        cfg['sim']['max_steps'] = '1'
        cfg['mesh'].update(xlength='1', ylength='1', zlength='1', resolution='1',
                           smallest_size='.25', largest_size='1', remeshing_option='1',
                           min_quality='.01', max_boundary_distortion='1e-9')
        cfg['control'].update(gravity='0', fixed_dt='.001', has_thermal_diffusion='false')
        cfg['ic'].update(weakzone_option='0', has_body_force_adjustment='false',
                         pore_pressure_option='2', pore_pressure_value='100000')
        reset_boundaries(cfg['bc'], vbc='1', vbc_val='0', hbc='0', stress_bc='0')
        velocity = '.1' if case == 'domain-add' else '-.1'
        cfg['bc'].update(vbc_val_z0=velocity, vbc_val_z1=velocity)
    return cfg


def inventories(model, frame, cfg):
    dim = int(model.ndims)
    x, cells = (model.read_field(frame, key) for key in ('coordinate', 'connectivity'))
    volume = simplex_volumes(x, cells)
    value = lambda key: json.loads(cfg['mat'][key])[0]
    phi, alpha = value('porosity'), value('biot_coeff')
    storage = phi/value('fluid_bulk_modulus')+(alpha-phi)/value('bulk_modulus_s')
    feedback = cfg.getboolean('control', 'has_poroelastic_pressure_feedback', fallback=True)
    mean_bulk = value('bulk_modulus')+(value('shear_modulus')/3 if dim == 2 else 0)
    if not feedback:
        storage += alpha**2/mean_bulk
    pressure = model.read_field(frame, 'pore pressure')[cells].mean(axis=1)
    strain = model.read_field(frame, 'strain')[:, :dim].sum(axis=1)
    mobile = float(volume@(storage*pressure+(alpha*strain if feedback else 0)))
    if cfg.getboolean('control', 'has_thermal_pressurization', fallback=False):
        beta = phi*value('fluid_alpha')+(alpha-phi)*value('alpha')
        if not feedback:
            beta -= alpha*value('bulk_modulus')*value('alpha')/mean_bulk
        temperature = model.read_field(frame, 'temperature')[cells].mean(axis=1)
        mobile -= float(volume@(beta*(temperature-cfg.getfloat('bc', 'surface_temperature'))))
    matrix = (float(volume@model.read_field(frame, 'matrix pressure')*value('matrix_storage'))
              if cfg.has_option('mat', 'matrix_storage') else 0.)
    return dict(volume=float(volume.sum()), mobile=mobile, matrix=matrix,
                max_pressure=float(np.max(np.abs(pressure))), nodes=len(x), cells=len(cells))


def continuation(exe, folder, cfg, rows):
    checkpoint = next(r['frame'] for r in rows if r['steps'] == 1)
    cfg['sim'].update(modelname='continued', is_restarting='true',
        restarting_from_modelname=str(folder/'flow'), restarting_from_frame=str(checkpoint))
    run_des(exe, folder, 'restart', cfg)
    count = 0
    for row in scan_frames(str(folder/'continued')):
        frame = row['frame']
        assert any(r['frame'] == frame and r['steps'] == row['steps'] and r['time'] == row['time'] for r in rows)
        for suffix in ('save', 'chkpt'):
            if frame == checkpoint and suffix == 'chkpt':
                continue
            a, b = (folder/f'{name}.{suffix}.{frame:06d}' for name in ('flow', 'continued'))
            assert (checkpoint_file(a) is None) == (checkpoint_file(b) is None)
            if checkpoint_file(a) is None:
                continue
            old, new = fields(a), fields(b)
            assert old.keys() == new.keys()
            different = physical_differences(old, new)
            assert not different, (frame, suffix, different)
            count += sum(not runtime_key(key) for key in old)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    records = []
    for case in ('drained', 'sealed', 'explicit', 'gvs', 'matrix', 'thermal',
                 'uncoupled', 'uncoupled-thermal', 'domain-add', 'domain-remove'):
        folder = output/case
        folder.mkdir()
        cfg = configuration(case)
        run_des(exe, folder, 'input', cfg)
        model = Dynearthsol(str(folder/'flow'))
        rows = scan_frames(model.modelname)
        pairs = [(a, b) for a, b in zip(rows, rows[1:]) if a['steps'] == b['steps'] and a['time'] == b['time']]
        remeshes = [(a, b) for a, b in pairs
                    if not np.array_equal(*(model.read_field(r['frame'], 'coordinate') for r in (a, b)))]
        log = (folder/'input.log').read_text()
        mobile_changes = [added-removed+reaction for removed, added, reaction in
                          (map(float, m) for m in re.findall(
                              r'Biot content remap: removed/added/Dirichlet \[m\^\d\] = (\S+) (\S+) (\S+);', log))]
        matrix_changes = [added-removed for removed, added in (map(float, m) for m in re.findall(
                              r'Matrix content remap: removed/added \[m\^\d\] = (\S+) (\S+)', log))]
        assert remeshes and len(remeshes) == len(mobile_changes), (case, len(remeshes), len(mobile_changes))
        defects = []
        for i, (a, b) in enumerate(remeshes):
            before, after = (inventories(model, r['frame'], cfg) for r in (a, b))
            scale = max(abs(before['mobile']), 1e-20)
            if case.startswith('domain'):
                swept = before['mobile']*(after['volume']/before['volume']-1)
                assert abs(mobile_changes[i]-swept) < 1e-10*scale, (case, a['frame'], mobile_changes[i], swept)
            else:
                assert abs(after['volume']/before['volume']-1) < 1e-6, (case, a['frame'], before['volume'], after['volume'])
            defect = after['mobile']-before['mobile']-mobile_changes[i]
            defects.append(abs(defect)/scale)
            assert defects[-1] < 1e-10, (case, a['frame'], b['frame'], defect)
            matrix_change = matrix_changes[i] if matrix_changes else 0.
            np.testing.assert_allclose(after['matrix']-before['matrix'], matrix_change,
                                       rtol=0, atol=1e-10*abs(before['matrix'])+1e-15)
        if case.startswith('domain'):
            before, after = (inventories(model, r['frame'], cfg) for r in remeshes[-1])
            np.testing.assert_allclose(after['volume']-before['volume'],
                                      1e-4 if case == 'domain-add' else -1e-4, rtol=1e-7)
            exact = 0
        else:
            exact = continuation(exe, folder, cfg, rows)
        records.append(dict(case=case, maximum_relative_remap_defect=max(defects),
                            exact_restart_payloads=exact, remeshes=len(remeshes)))
        print('PASS fluid remap', case, flush=True)
    write_json(output/'summary.json', dict(passed=True, cases=records,
        scope='Fixed-property saturated linear Biot inventories; not evolving-storage/unsaturated qualification'))


if __name__ == '__main__':
    main()
