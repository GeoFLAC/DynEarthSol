#!/usr/bin/env python3
"""Recompute the heat-transport qualification from frozen outputs and hashes."""

import argparse
import configparser
import importlib.util
import json
from pathlib import Path

import numpy as np

import consolidation_refinement as cr
import prepare_figures as pf


def load_case(study, spacing, steps, consumed, executable=None):
    def record(path):
        consumed[str(path.resolve())] = cr.digest(path)
    record(study/'summary.json')
    summary = json.loads((study/'summary.json').read_text())
    item = next(r for r in summary['cases'] if r['spacing_m'] == spacing and r['steps'] == steps)
    folder = study/f'h{spacing:g}-n{steps}'
    record(folder/'fresh.cfg'); record(study/'runner.py'); record(study/'source.diff')
    assert cr.digest(folder/'fresh.cfg') == item['input_sha256']
    assert cr.digest(study/'runner.py') == summary['runner_sha256']
    binary = Path(summary['executable']) if executable is None else Path(executable)
    assert cr.digest(binary) == summary['executable_sha256']
    record(binary)
    for name, sha in item['outputs_sha256'].items():
        path = folder/name
        assert cr.digest(path) == sha, path
        record(path)
    spec = importlib.util.spec_from_file_location('frozen_heat_reference', study/'runner.py')
    reference = importlib.util.module_from_spec(spec); spec.loader.exec_module(reference)
    cfg = configparser.ConfigParser(); cfg.read(folder/'fresh.cfg')
    model = pf.Dynearthsol(str(folder/'column'))
    x = model.read_field(0, 'coordinate')[:, 0]
    rows = pf.scan_frames(str(folder/'column'))
    times = np.array([r['time'] for r in rows])
    theta = np.array([model.read_field(r['frame'], 'temperature')-273.15 for r in rows])
    ref = summary['reference']
    exact = np.array([reference.exact(x, t, ref['velocity'], ref['diffusivity']) for t in times])
    assert np.isfinite(theta).all() and theta.min() >= -1e-10 and theta.max() <= 1+1e-10
    errors = np.max(np.abs(theta-exact), axis=1)
    np.testing.assert_allclose(errors, item['max_error'], rtol=0, atol=1e-12)
    np.testing.assert_allclose(times[-1], ref['duration_s'], rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(times[-1], steps*cfg.getfloat('control', 'fixed_dt'), rtol=1e-12)
    for row in rows:
        np.testing.assert_allclose(model.read_field(row['frame'], 'pore pressure'), 2-x, rtol=0, atol=1e-10)
    assert item['physical_restart_exact']
    for path in folder.glob('continued.*'):
        if '.save.' in path.name or path.suffix == '.info':
            record(path)
    resumed = pf.Dynearthsol(str(folder/'continued'))
    by_step = {r['steps']: r for r in rows}
    for row in pf.scan_frames(str(folder/'continued')):
        fresh = by_step[row['steps']]
        assert row['time'] == fresh['time']
        for key in ('temperature', 'pore pressure', 'velocity', 'stress', 'strain', 'coordinate'):
            np.testing.assert_array_equal(model.read_field(fresh['frame'], key), resumed.read_field(row['frame'], key))
    return dict(coordinate=x, times=times, theta=theta, exact=exact, errors=errors,
                reference=ref, item=item, folder=folder, exact_function=reference.exact)


def qualify(group):
    consumed = {}
    verification = json.loads((group/'verification.json').read_text())
    checks = verification['checks']
    assert len(checks) == 94 and all(r['exit'] == 0 for r in checks.values())
    consumed[str((group/'verification.json').resolve())] = cr.digest(group/'verification.json')
    dimensions = {}
    for dim, spacings in ((2, [.025, .0125, .00625, .003125]), (3, [.0125, .00625, .003125])):
        cases = [load_case(group/f'heat-spatial-{dim}d', h, 8192, consumed) for h in spacings]
        errors = [float(r['errors'].max()) for r in cases]
        orders = np.log2(np.array(errors[:-1])/errors[1:])
        assert errors[-1] < .014 and np.all(np.diff(errors) < 0)
        assert orders[-1] > .95
        dimensions[str(dim)] = dict(spacings_m=spacings, max_errors_over_jump=errors,
                                   spatial_orders=orders.tolist(), restart_exact=True)
    cases = [load_case(group/'heat-temporal-2d', .0125, n, consumed) for n in (512, 1024, 2048, 4096)]
    differences = [float(np.max(np.abs(a['theta']-b['theta']))) for a, b in zip(cases[:-1], cases[1:])]
    temporal_orders = np.log2(np.array(differences[:-1])/differences[1:])
    assert np.all((temporal_orders > .95) & (temporal_orders < 1.1))
    limits = {}
    for dim in (2, 3):
        case = load_case(group/f'heat-fluid-limit-{dim}d', .00625, 4096, consumed)
        assert case['reference']['porosity'] == 1 and max(case['errors']) < .03
        limits[str(dim)] = float(max(case['errors']))
    large = load_case(group/'heat-large-fixed-step-2d', .003125, 128, consumed)
    assert max(large['errors']) < .014
    inputs = json.loads((group/'heat-input-checks.json').read_text())
    assert inputs['passed'] and len(inputs['checks']) == 26 and all(inputs['checks'].values())
    consumed[str((group/'heat-input-checks.json').resolve())] = cr.digest(group/'heat-input-checks.json')
    result = dict(passed=True, dimensions=dimensions,
        temporal=dict(steps=[512, 1024, 2048, 4096], max_pair_differences=differences,
                      orders=temporal_orders.tolist()),
        fluid_only_capacity_limit_max_error=limits,
        large_fixed_step_max_error=float(max(large['errors'])),
        legacy_regression_checks=len(checks),
        selected=dict(study='heat-spatial-3d', spacing_m=.003125, steps=8192),
        consumed_files=consumed, input_checks=inputs,
        scope='Constant-capacity divergence-free transport with graph diffusion; scaled hydraulic storage; no general evolving-storage enthalpy claim.')
    cr.write_json(group/'heat-qualification.json', result)
    return {k: v for k, v in result.items() if k not in ('consumed_files', 'input_checks')}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('group', type=Path)
    print(json.dumps(qualify(parser.parse_args().group.resolve()), indent=2))
