#!/usr/bin/env python3
"""Check a manual coupled explicit/GVS-to-physical-inertia restart.

Stress, strain, pressure and velocity must survive the switch. Independent
undrained wave-speed and velocity/force identities check the first physical
interval and its mixture mass. Artificial kinetic energy is reported apart
from physical kinetic energy; the prescribed mass change supplies no heat.
This is not an automatic earthquake-event or radiated-energy qualification.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks/poromechanics'), str(ROOT/'tests/functional')]
from Dynearthsol import Dynearthsol, scan_frames
import fluid_content as water
from des_test_utils import (compare_continuation, fields, mesh_measures, physical_differences, record_name,
                            roller_free, runtime_key, simplex_volumes)

K, G, PHI, ALPHA, KF, KS = 2e7, 1e7, .25, .8, 1e8, 1e9
DENSITY = (1-PHI)*2650+PHI*1000
STORAGE = PHI/KF+(ALPHA-PHI)/KS
UNDRAINED = K+4*G/3+ALPHA**2/STORAGE


sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import check_hashes, digest, file_hashes, provenance, reset_boundaries, run_model, write_json


def configuration(method, scaling):
    cfg = water.configuration(.5, .01, 1., method, scaling, water.INFLOW)
    cfg['sim'].update(max_steps='32', max_time_in_yr='1e6')
    cfg['control'].update(fixed_dt='0', damping_option='0', is_using_mixed_stress='false')
    reset_boundaries(cfg['bc'], hbc='0')
    cfg['bc'].update(hbc_flux='[0,0,0,0,0,0]', stress_val_z1='-1000')
    return cfg


def analyze(folder, method, scaling):
    warm, dynamic = (Dynearthsol(str(folder/name)) for name in ('warm', 'physical'))
    assert np.max(abs(warm.read_field(32, 'pore pressure'))) > 1
    assert np.max(abs(warm.read_field(32, 'strain'))) > 1e-9
    before = fields(folder/'warm.save.000032')
    after = fields(folder/'physical.save.000032')
    keys = {key for key in before.keys() | after.keys() if not runtime_key(key) and record_name(key) != 'dt_sec'}
    different = [key for key in physical_differences(before, after) if record_name(key) != 'dt_sec']
    assert not different, different
    x, height, lumped = mesh_measures(dynamic)
    expected = .2*min(.5*height/math.sqrt(K/2650), height/(5*math.sqrt(UNDRAINED/DENSITY)))
    rows = scan_frames(str(folder/'physical'))
    interval = rows[1]['time']-rows[0]['time']
    np.testing.assert_allclose(interval, expected, rtol=2e-12)
    # Restart emits an initial save, then checkpoints after accepted steps.
    checkpoint = fields(folder/'physical.chkpt.000033')
    np.testing.assert_allclose(np.frombuffer(checkpoint['dt'], dtype=np.float64)[0], expected, rtol=2e-13)
    free = roller_free(x)
    v0 = dynamic.read_field(dynamic.frames[0], 'velocity')
    v1 = dynamic.read_field(dynamic.frames[1], 'velocity')
    force = dynamic.read_field(dynamic.frames[1], 'force')
    dv = expected*force/(DENSITY*lumped[:, None])
    error = float(np.max(abs((v1-v0-dv)[free]))/
                  max(np.max(abs(dv[free])), np.max(abs(v0[free])), 1e-30))
    assert error < 1e-12, error

    # Preserved velocity means continuous physical kinetic energy. Artificial
    # mass energy changes by prescription and must not become supplied heat.
    conn = dynamic.read_field(dynamic.frames[0], 'connectivity')
    saved = fields(folder/'warm.chkpt.000032')
    vmax = float(np.frombuffer(saved['max_global_vel_mag'], dtype=np.float64)[0])
    old_dt = float(np.frombuffer(saved['dt'], dtype=np.float64)[0])
    scale = max(1e-5, vmax) if method.startswith('gvs') else 1e-5
    speed = scale*scaling
    if method.startswith('gvs'):
        speed = min(speed, math.sqrt((G if method == 'gvs-shear' else K)/2650))
    artificial = K/speed**2
    velocity_integral = float(lumped @ np.sum(v0*v0, axis=1))
    old_energy, physical_energy = .5*artificial*velocity_integral, .5*DENSITY*velocity_integral
    assert old_energy > 0 and physical_energy > 0 and old_dt != expected
    dim = int(dynamic.ndims)
    volume = simplex_volumes(x, conn)
    budget = []
    for row in rows:
        frame = row['frame']
        p = dynamic.read_field(frame, 'pore pressure')
        strain = dynamic.read_field(frame, 'strain')[:, :dim].sum(axis=1)
        stored, deformed = STORAGE*p[conn].mean(axis=1), ALPHA*strain
        content = float(volume @ (stored+deformed))
        scale = float(volume @ (abs(stored)+abs(deformed)))
        budget.append(abs(content)/max(scale, 1e-30))
    assert max(budget) < 1e-11, budget
    resumed = Dynearthsol(str(folder/'continuation/resumed'))
    exact = compare_continuation(folder/'continuation', str(folder/'physical'), resumed)
    return dict(method=method, scaling=scaling, dimension=dim, restored_save_payloads_exact=len(keys),
                old_artificial_dt_s=old_dt, expected_physical_dt_s=expected, first_physical_interval_s=interval,
                velocity_mass_identity_error=error, maximum_relative_closed_content_error=max(budget),
                artificial_density_kg_m3=artificial, physical_density_kg_m3=DENSITY,
                artificial_kinetic_energy_J=old_energy, physical_kinetic_energy_J=physical_energy,
                numerical_mass_energy_change_J=physical_energy-old_energy, exact_physical_restart_payloads=exact)


def case(exe, folder, method, scaling):
    folder.mkdir(exist_ok=False)
    cfg = configuration(method, scaling)
    run_model(exe, folder, cfg, 'warm', timeout=120)
    physical = copy.deepcopy(cfg)
    physical['sim'].update(is_restarting='true', restarting_from_modelname=str(folder/'warm'),
                           restarting_from_frame='32', max_steps='38')
    physical['control'].update(is_quasi_static='false', use_global_velocity_scaling='false',
                               mass_scaling_reference_speed='shear')
    run_model(exe, folder, physical, 'physical', timeout=120)
    resumed = copy.deepcopy(physical)
    resumed['sim'].update(restarting_from_modelname=str(folder/'physical'), restarting_from_frame='34')
    sub = folder/'continuation'
    sub.mkdir()
    run_model(exe, sub, resumed, 'resumed', timeout=120)
    return analyze(folder, method, scaling)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    if args.replay:
        report = json.loads((output/'summary.json').read_text())
        check_hashes(report['artifacts_sha256'], output)
        for previous in report['cases']:
            name = f"{previous['method']}-S{previous['scaling']:g}"
            assert analyze(output/name, previous['method'], previous['scaling']) == previous, name
        print('PASS: inertia switch references, restart payloads and artifact hashes replayed')
        return
    output.mkdir(parents=True, exist_ok=False)
    report = dict(provenance(exe, __file__), cases=[], passed=False)
    for method in ('explicit', 'gvs-shear', 'gvs-bulk'):
        for scaling in (1e4, 1e5, 1e6):
            result = case(exe, output/f'{method}-S{scaling:g}', method, scaling)
            report['cases'].append(result)
            print(json.dumps(result), flush=True)
            write_json(output/'summary.json', report)
    assert digest(exe) == report['executable_sha256']
    report['artifacts_sha256'] = file_hashes((p for p in output.rglob('*') if p.name != 'summary.json'), output)
    report['passed'] = True
    write_json(output/'summary.json', report)


if __name__ == '__main__':
    main()
