#!/usr/bin/env python3
"""Replay and qualify the closed-column physical-inertia wave study.

Refinement data stay in the evidence record. The guided review selects one
finest-grid case only. No material, phase or amplitude fit is performed.
"""

import argparse
import configparser
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import prepare_figures as pf


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_case(study, count, consumed=None):
    """Reconstruct the reference from the actual deck and initial velocity."""
    summary = json.loads((study / 'summary.json').read_text())
    case = summary['cases'][str(count)]
    folder = study / f'steps-{count}'

    def record(path):
        if consumed is not None:
            consumed[str(path.resolve())] = digest(path)

    record(study / 'summary.json')
    assert digest(folder/'wave.cfg') == case['input_sha256']
    assert digest(folder/'seed.save.000000') == case['seed_save_sha256']
    assert digest(folder/'seed.chkpt.000000') == case['seed_checkpoint_sha256']
    for name, sha in case.get('outputs_sha256', {}).items():
        assert digest(folder/name) == sha, folder/name
    for path in folder.iterdir():
        if path.name.startswith(('wave.save.', 'wave.chkpt.', 'wave.info', 'wave.manifest',
                                 'wave.cfg', 'seed.save.000000', 'seed.chkpt.000000')):
            record(path)
    cfg = configparser.ConfigParser(); cfg.read(folder/'wave.cfg')
    for key in ('has_PT', 'is_quasi_static', 'has_thermal_diffusion',
                'has_moving_mesh', 'use_global_velocity_scaling'):
        assert not cfg.getboolean('control', key)
    assert cfg.getint('control', 'damping_option') == 0
    assert cfg.getfloat('control', 'gravity') == 0
    for side in ('x0', 'x1', 'y0', 'y1', 'z0', 'z1'):
        assert cfg.getint('bc', f'hbc_{side}') == 0
    value = lambda key: json.loads(cfg['mat'][key])[0]
    assert value('retention_model') == value('hydraulic_perm') == value('alpha') == 0
    bulk, shear, phi, alpha = (value(k) for k in ('bulk_modulus', 'shear_modulus', 'porosity', 'biot_coeff'))
    storage = phi/value('fluid_bulk_modulus')+(alpha-phi)/value('bulk_modulus_s')
    density = (1-phi)*value('rho0')+phi*value('fluid_rho0')
    speed = math.sqrt((bulk+4*shear/3+alpha**2/storage)/density)
    length = cfg.getfloat('mesh', 'zlength')
    period = 2*length/speed
    model = pf.Dynearthsol(str(folder/'wave'))
    rows = pf.scan_frames(str(folder/'wave'))
    x = model.read_field(0, 'coordinate')
    conn = model.read_field(0, 'connectivity')
    dim = model.ndims
    volume = np.abs(np.linalg.det(x[conn][:,1:]-x[conn][:,:1]))/math.factorial(dim)
    nodal_volume = np.zeros(len(x))
    np.add.at(nodal_volume, conn, np.broadcast_to((volume/(dim+1))[:,None], conn.shape))
    amplitude = float(model.read_field(0, 'velocity')[:,-1].max())
    pressure_amplitude = alpha*amplitude/(storage*speed)
    phase_z = np.pi*(x[:,-1]+length)/length
    np.testing.assert_allclose(model.read_field(0, 'velocity')[:,-1], amplitude*np.sin(phase_z),
                               rtol=1e-13, atol=1e-20)
    velocities, pressures, measured = [], [], []
    for row in rows:
        frame, time = row['frame'], row['time']
        v = model.read_field(frame, 'velocity')[:,-1]
        p = model.read_field(frame, 'pore pressure')
        strain = model.read_field(frame, 'strain')
        assert np.isfinite(v).all() and np.isfinite(p).all() and np.isfinite(strain).all()
        projection = np.zeros(len(x))
        np.add.at(projection, conn, np.broadcast_to(
            (volume*strain[:,:dim].sum(axis=1)/(dim+1))[:,None], conn.shape))
        content = alpha*projection/nodal_volume+storage*p
        velocities.append(v); pressures.append(p)
        measured.append(dict(time_s=time,
            velocity_error_over_amplitude=float(np.max(np.abs(v-amplitude*np.sin(phase_z)*np.cos(2*np.pi*time/period)))/amplitude),
            pressure_error_over_amplitude=float(np.max(np.abs(p+pressure_amplitude*np.cos(phase_z)*np.sin(2*np.pi*time/period)))/pressure_amplitude),
            max_local_fluid_content_error=float(np.max(np.abs(content))),
            local_content_error_over_strain_amplitude=float(np.max(np.abs(content))/(alpha*amplitude/speed))))
    for key in measured[0]:
        np.testing.assert_allclose([r[key] for r in measured], [r[key] for r in case['measured']],
                                   rtol=1e-10, atol=1e-14 if key != 'time_s' else period*1e-13)
    times = np.array([r['time'] for r in rows])
    np.testing.assert_allclose(times[-1], period, rtol=1e-12, atol=0)
    return dict(coordinate=x, times=times, velocity=np.array(velocities), pressure=np.array(pressures),
                velocity_amplitude=amplitude, pressure_amplitude=pressure_amplitude, period_s=period,
                dimension=dim, spacing_m=case['spacing_m'], steps=count, measured=measured,
                maximum={k:max(r[k] for r in measured) for k in measured[0] if k != 'time_s'})


def qualify(output, spatial):
    consumed, dimensions = {}, {}
    for dim in (2, 3):
        study = output/f'wave-{dim}d'
        summary = json.loads((study/'summary.json').read_text())
        assert digest(study/'undrained_wave.py') == summary['analyzer_sha256']
        assert digest(output/f'qualified-binaries/dynearthsol{dim}d') == summary['executable_sha256']
        cases = [load_case(study, n, consumed) for n in (1024, 2048, 4096, 8192)]
        differences = []
        for a, b in zip(cases, cases[1:]):
            np.testing.assert_array_equal(a['coordinate'], b['coordinate'])
            np.testing.assert_allclose(a['times'], b['times'], rtol=0, atol=b['period_s']*1e-12)
            differences.append([float(np.max(np.abs(a[key]-b[key]))/b[f'{key}_amplitude'])
                                for key in ('velocity', 'pressure')])
        orders = np.log2(np.array(differences[:-1])/np.array(differences[1:]))
        assert ((orders[:,0]>.9) & (orders[:,0]<1.1)).all(), orders
        assert ((orders[:,1]>1.8) & (orders[:,1]<2.2)).all(), orders
        finest = cases[-1]
        for key in ('velocity_error_over_amplitude', 'pressure_error_over_amplitude'):
            assert finest['maximum'][key] < .003, finest['maximum']
        assert finest['maximum']['local_content_error_over_strain_amplitude'] < 5e-12
        for c in summary['cases'].values():
            assert c['restart_physical_fields_exact'] and c['incompatible_ordering_rejected']
        assert summary['automatic_step']['matches_undrained_wave_limit']
        spatial_cases = [load_case(p,4096,consumed) for p in spatial[dim]]+[cases[-2]]
        h = [c['spacing_m'] for c in spatial_cases]
        errors = [c['maximum']['pressure_error_over_amplitude'] for c in spatial_cases]
        h_orders = [math.log(a/b)/math.log(x/y) for a,b,x,y in zip(errors,errors[1:],h,h[1:])]
        assert all(1.8 < o < 2.2 for o in h_orders), h_orders
        dimensions[str(dim)] = dict(finest_maximum=finest['maximum'], temporal_pair_differences=differences,
            temporal_orders_velocity_pressure=orders.tolist(), spatial_pressure_errors=errors,
            spacings_m=h, spatial_pressure_orders=h_orders, exact_physical_restart=True,
            incompatible_ordering_rejected=True, automatic_step=summary['automatic_step'])
    report_paths = [p for p in output.glob('compare-*.json') if not p.name.startswith('compare-stamped-')]
    reports = [json.loads(p.read_text()) for p in report_paths]
    verification = json.loads((output/'verification.json').read_text())
    assert verification['passed'] and not any(r['differences'] for r in reports)
    old_restart = json.loads((output/'old-checkpoint-rejection.json').read_text())
    assert all(old_restart[str(d)]['rejected_exit'] == 22 for d in (2,3))
    stamped_path = output/'stamped-restart-final-verification.json'
    assert json.loads(stamped_path.read_text())['passed']
    for dim in (2,3):
        converted_path = output/f'stamped-restart-final-wave-{dim}d/summary.json'
        converted = json.loads(converted_path.read_text())
        assert converted['automatic_step']['initial_ordering_conversion_reselects_step']
        consumed[str(converted_path.resolve())] = digest(converted_path)
    consumed[str(stamped_path.resolve())] = digest(stamped_path)
    for path in report_paths+[output/'verification.json',output/'old-checkpoint-rejection.json']:
        consumed[str(path.resolve())] = digest(path)
    for dim in (2,3):
        log = (output/f'assembly-final-{dim}d.log').read_text()
        assert '0 failed' in log
        consumed[str((output/f'assembly-final-{dim}d.log').resolve())] = digest(output/f'assembly-final-{dim}d.log')
    result = dict(passed=True, scope='Closed saturated elastic physical-inertia wave, one period; no GVS/PT or earthquake-cycle qualification.',
        dimensions=dimensions, selected=dict(dimension=3,study='wave-3d',steps=8192,spacing_m=.0125),
        regression_checks=len(verification['checks']), files_compared=sum(r['files_compared'] for r in reports),
        physical_payloads_compared=sum(r['physical_payloads_compared'] for r in reports),
        differences=0, analyzer_sha256=digest(Path(__file__)), consumed_files=consumed)
    (output/'wave-qualification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'consumed_files'},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path,help='Evidence group containing wave-2d and wave-3d')
    parser.add_argument('--spatial-2d',type=Path,nargs=2,required=True,help='Coarse then middle study, both at 4096 steps')
    parser.add_argument('--spatial-3d',type=Path,nargs=2,required=True)
    args = parser.parse_args()
    qualify(args.output.resolve(),{2:args.spatial_2d,3:args.spatial_3d})
