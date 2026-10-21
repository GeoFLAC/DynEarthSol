#!/usr/bin/env python3
"""Check fixed-load PT equilibrium, force output, restart and a changed load."""
import argparse
import configparser
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Dynearthsol import Dynearthsol


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
                              initial=float(initial), threshold=float(threshold))
               for step, iterations, residual, initial, threshold in re.findall(
                   r'PT step=(\d+) status=converged iterations=(\d+) residual=(\S+) initial=(\S+) threshold=(\S+)',
                   (output / 'fresh.log').read_text())}
    assert reports[1]['initial'] > 1, 'Example must exercise a nonzero force imbalance'
    if cfg.getboolean('control', 'PT_use_running_scale', fallback=False):
        for step in (2, 3):
            assert reports[step]['iterations'] == 1, 'Held equilibrium must not require relaxation updates'
            np.testing.assert_array_equal(des.read_field(step, 'velocity'), 0)
            np.testing.assert_array_equal(des.read_field(step, 'strain'), des.read_field(1, 'strain'))
    checked = {}
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
        reported, threshold = report['residual'], report['threshold']
        np.testing.assert_allclose(rms, reported, rtol=1e-12, atol=1e-14)
        assert rms <= threshold and np.isfinite(residual).all()
        checked[int(frame)] = {'rms': rms, 'threshold': threshold, 'iterations': report['iterations']}
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
