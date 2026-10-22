#!/usr/bin/env python3
"""Check fixed-load PT equilibrium, force output, restart and a changed load."""
import argparse
import configparser
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Dynearthsol import Dynearthsol


def stress_tensor(s):
    if len(s) == 3:
        return np.array([[s[0], s[2]], [s[2], s[1]]])
    return np.array([[s[0], s[3], s[4]], [s[3], s[1], s[5]], [s[4], s[5], s[2]]])


def stress_components(s):
    if len(s) == 2:
        return np.array([s[0, 0], s[1, 1], s[0, 1]])
    return np.array([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[0, 2], s[1, 2]])


def box_force(coordinate, connectivity, stress, traction, free):
    """Independent simplex integration, using saved physical stress and geometry."""
    ndims = coordinate.shape[1]
    force = np.zeros_like(coordinate)
    facets = {}
    for nodes, sigma in zip(connectivity, stress):
        vertices = np.column_stack((np.ones(ndims + 1), coordinate[nodes]))
        gradient = np.linalg.inv(vertices)[1:, :].T
        volume = abs(np.linalg.det(vertices)) / math.factorial(ndims)
        force[nodes] -= volume * (gradient @ stress_tensor(sigma).T)
        for face in itertools.combinations(nodes, ndims):
            key = tuple(sorted(int(node) for node in face))
            facets[key] = facets.get(key, 0) + 1
    for face, count in facets.items():
        vertices = coordinate[list(face)]
        if count != 1 or not np.all(np.isclose(vertices[:, -1], 0)):
            continue
        area = np.linalg.norm(vertices[1] - vertices[0])
        if ndims == 3:
            area = np.linalg.norm(np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])) / 2
        force[list(face), -1] += traction * area / ndims
    return np.where(free, force, 0)


def before_rotation(coordinate, connectivity, stress, velocity, dt):
    """Invert one linear Jaumann update independently of the production routine."""
    ndims = coordinate.shape[1]
    result = stress.copy()
    for element, nodes in enumerate(connectivity):
        vertices = np.column_stack((np.ones(ndims + 1), coordinate[nodes]))
        gradient = velocity[nodes].T @ np.linalg.inv(vertices)[1:, :].T
        if ndims == 2:
            spin = (gradient[1, 0] - gradient[0, 1]) / 2
            operator = np.array([[0, 0, -2*spin], [0, 0, 2*spin], [spin, -spin, 0]])
        else:
            spin = (gradient - gradient.T) / 2
            operator = np.column_stack([
                stress_components(stress_tensor(basis) @ spin - spin @ stress_tensor(basis))
                for basis in np.eye(6)])
        result[element] = np.linalg.solve(np.eye(len(stress[element])) + dt*operator, stress[element])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('config', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--option', type=int, choices=(0, 1), default=0)
    parser.add_argument('--timeout', type=float, default=60,
                        help='Maximum seconds per simulation process (default: 60)')
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error('--timeout must be positive')
    executable = args.executable.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(args.config)
    assert not cfg.getboolean('control', 'has_moving_mesh')
    assert not cfg.getboolean('control', 'has_hydraulic_diffusion')
    assert cfg.getint('control', 'damping_option') == 1
    assert cfg.getfloat('control', 'gravity') == 0
    cfg['control']['PT_option'] = str(args.option)
    for phase in ('fresh', 'resumed'):
        cfg['sim'].update(modelname=phase, is_restarting=str(phase == 'resumed').lower())
        if phase == 'resumed':
            cfg['sim'].update(restarting_from_modelname=str(output / 'fresh'), restarting_from_frame='1')
        with (output / f'{phase}.cfg').open('w') as stream:
            cfg.write(stream)
        with (output / f'{phase}.log').open('w') as stream:
            subprocess.run([str(executable), f'{phase}.cfg'], cwd=output,
                           stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=args.timeout)
    des = Dynearthsol(str(output / 'fresh'))
    try:
        des.read_field(0, 'PT residual force')
    except KeyError:
        pass
    else:
        raise AssertionError('Unsolved initial frame must not claim a PT residual')
    reports = {int(step): dict(iterations=int(iterations), residual=float(residual),
                              initial=float(initial), threshold=float(threshold),
                              warm_initial=float(warm), evaluations=int(evaluations))
               for step, iterations, residual, initial, threshold, warm, evaluations in re.findall(
                   r'PT step=(\d+) status=converged iterations=(\d+) residual=(\S+) initial=(\S+) threshold=(\S+) '
                   r'warm_initial=(\S+) evaluations=(\d+) stage=before_physical_updates',
                   (output / 'fresh.log').read_text())}
    assert set(reports) == {1, 2, 3}, 'Missing complete PT diagnostic records'
    assert reports[1]['initial'] > 1, 'Example must exercise a nonzero force imbalance'
    if cfg.getboolean('control', 'PT_use_running_scale', fallback=False):
        for step in (2, 3):
            assert reports[step]['iterations'] == 1, 'Held equilibrium must not require relaxation updates'
            np.testing.assert_array_equal(des.read_field(step, 'velocity'), 0)
            np.testing.assert_array_equal(des.read_field(step, 'strain'), des.read_field(1, 'strain'))
    checked = {}
    physical_end_rms = {}
    for frame in des.frames[1:]:
        force = des.read_field(frame, 'force')
        velocity = des.read_field(frame, 'velocity')
        flags = des.read_field(frame, 'bcflag')
        residual = des.read_field(frame, 'PT residual force')
        free = np.ones(force.shape, dtype=bool)
        free[(flags & 3) != 0, 0] = False       # roller x sides
        if des.ndims == 3:
            free[(flags & 12) != 0, 1] = False  # roller y sides
        free[(flags & 16) != 0, :] = False     # clamped base
        factor = np.where(np.abs(velocity) > 1e-13,
                          1 - cfg.getfloat('control', 'damping_factor') * np.sign(force) * np.sign(velocity), 1)
        expected = np.where(free, force / factor, 0)
        np.testing.assert_allclose(residual, expected, rtol=1e-12, atol=1e-14)
        assert np.all(residual[~free] == 0)
        rms = float(np.sqrt(np.sum(residual * residual) / free.sum()))
        report = reports[des.steps[des.frames.index(frame)]]
        assert report['iterations'] <= report['evaluations'] <= report['iterations'] + 2
        reported, threshold = report['residual'], report['threshold']
        np.testing.assert_allclose(rms, reported, rtol=1e-12, atol=1e-14)
        assert rms <= threshold and np.isfinite(residual).all()
        coordinate = des.read_field(frame, 'coordinate')
        connectivity = des.read_field(frame, 'connectivity')
        stress = des.read_field(frame, 'stress')
        traction = cfg.getfloat('bc', 'stress_val_z1')
        trial_stress = before_rotation(coordinate, connectivity, stress, velocity,
                                       cfg.getfloat('control', 'fixed_dt'))
        independent = box_force(coordinate, connectivity, trial_stress, traction, free)
        np.testing.assert_allclose(independent, residual, rtol=1e-10, atol=5e-11)
        physical = box_force(coordinate, connectivity, stress, traction, free)
        physical_end_rms[frame] = float(np.sqrt(np.sum(physical*physical) / free.sum()))
        checked[int(frame)] = dict(report, rms=rms, physical_end_rms=physical_end_rms[frame])
    if cfg.getboolean('control', 'PT_use_running_scale', fallback=False):
        for step in (2, 3):
            np.testing.assert_allclose(reports[step]['initial'], physical_end_rms[step-1], rtol=1e-10, atol=5e-11)
            assert reports[step]['initial'] == reports[step]['residual']
        assert reports[2]['warm_initial'] > reports[2]['threshold']
        assert reports[2]['evaluations'] == 2 and reports[3]['evaluations'] == 1
        assert reports[3]['warm_initial'] == reports[3]['initial']
        for report in reports.values():
            np.testing.assert_allclose(report['threshold'], 1e-6 * reports[1]['warm_initial'], rtol=1e-14)
    resumed = Dynearthsol(str(output / 'resumed'))
    for frame in resumed.frames:
        np.testing.assert_array_equal(des.read_field(frame, 'PT residual force'),
                                      resumed.read_field(frame, 'PT residual force'))
    # A cached equilibrium must not suppress a new load after restarting from
    # a held step with zero velocity. Re-evaluate the complete physical problem.
    cfg['sim'].update(modelname='reloaded', is_restarting='yes', restarting_from_frame='2')
    load = 2 * cfg.getfloat('bc', 'stress_val_z1')
    cfg['bc']['stress_val_z1'] = str(load)
    with (output / 'reloaded.cfg').open('w') as stream:
        cfg.write(stream)
    with (output / 'reloaded.log').open('w') as stream:
        subprocess.run([str(executable), 'reloaded.cfg'], cwd=output,
                       stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=args.timeout)
    reloaded = Dynearthsol(str(output / 'reloaded'))
    np.testing.assert_allclose(reloaded.read_field(3, 'stress')[:, des.ndims-1], load, rtol=1e-5)
    assert np.max(np.abs(reloaded.read_field(3, 'strain') - des.read_field(2, 'strain'))) > 1e-4
    (output / 'residual-check.json').write_text(json.dumps(checked, indent=2) + '\n')
    print(f'PASS: PT_option={args.option}, {des.ndims}D force/RMS, held equilibrium, restart and changed load')


if __name__ == '__main__':
    main()
