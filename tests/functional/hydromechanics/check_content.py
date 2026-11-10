#!/usr/bin/env python3
"""Integrate saturated constant-property content independently of its residual.

All accepted steps must be saved. Fixed geometry, fixed porosity, no
retention, no thermal source and no wells are required. The initial pending
PT increment comes from the run or its source checkpoint; it is not assumed
zero at a nonzero restart. --measure retains a failing parent diagnostic.
"""
import argparse
import configparser
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'benchmarks/poromechanics')]
from Dynearthsol import Dynearthsol, scan_frames
import fluid_content as water


def analyze(config, until=None):
    config = config.resolve()
    cfg = configparser.ConfigParser(); cfg.optionxform = str
    cfg.read(config)
    assert cfg.getboolean('control', 'has_hydraulic_diffusion', fallback=False)
    assert cfg.getboolean('control', 'has_poroelastic_pressure_feedback', fallback=False)
    assert not cfg.getboolean('control', 'has_moving_mesh', fallback=False)
    assert not cfg.getboolean('control', 'has_thermal_pressurization', fallback=False)
    assert cfg.getint('mat', 'porosity_model', fallback=0) == 0
    assert not any(json.loads(cfg.get('mat', 'retention_model', fallback='[0]')))

    def prop(name):
        values = json.loads(cfg['mat'][name])
        assert len(set(values)) == 1, ('nonuniform coefficient', name)
        return values[0]

    phi, alpha = prop('porosity'), prop('biot_coeff')
    storage = phi/prop('fluid_bulk_modulus')+(alpha-phi)/prop('bulk_modulus_s')
    model = Dynearthsol(str(config.parent/cfg['sim']['modelname']))
    rows = scan_frames(model.modelname)
    if until is not None:
        rows = [row for row in rows if row['time'] <= until]
    assert rows
    dim = int(model.ndims)
    lag = alpha**2/(prop('bulk_modulus')+(prop('shear_modulus')/3 if dim == 2 else 0)) \
        if cfg.getboolean('control', 'has_PT', fallback=False) else 0.
    x = model.read_field(rows[0]['frame'], 'coordinate')
    conn = model.read_field(rows[0]['frame'], 'connectivity')
    volume = abs(np.linalg.det(x[conn][:, 1:]-x[conn][:, :1]))/math.factorial(dim)
    previous = model.read_field(rows[0]['frame'], 'pore pressure')
    if lag:
        # A restart's initial save may have no corresponding new checkpoint.
        prefix, frame = Path(model.modelname), rows[0]['frame']
        checkpoint = Path(f'{prefix}.chkpt.{frame:06d}')
        if checkpoint.is_file() or Path(str(checkpoint)+'.vtkhdf').is_file():
            source = model
        else:
            assert cfg.getboolean('sim', 'is_restarting', fallback=False)
            prefix = config.parent/cfg['sim']['restarting_from_modelname']
            frame = cfg.getint('sim', 'restarting_from_frame')
            source = Dynearthsol(str(prefix))
        previous = previous+water.checkpoint_field(source, frame, 'pending pore pressure increment', len(previous))
    cumulative, initial, error, scale, residual = 0., None, 0., 0., 0.
    history = []
    for i, row in enumerate(rows):
        if i:
            assert row['steps']-rows[i-1]['steps'] == 1, 'save every accepted step for a boundary budget'
        frame = row['frame']
        np.testing.assert_array_equal(model.read_field(frame, 'coordinate'), x)
        pressure = model.read_field(frame, 'pore pressure')
        strain = model.read_field(frame, 'strain')[:, :dim].sum(axis=1)
        physical = float(volume @ (storage*pressure[conn].mean(axis=1)+alpha*strain))
        pending = float(lag*(volume @ (pressure-previous)[conn].mean(axis=1)))
        balance = model.read_field(frame, 'fluid balance')
        assert balance[2] == 0, 'this check excludes well sources'
        if i:
            cumulative += float(balance[3])
        else:
            initial = physical+pending
        defect = physical+pending-initial-cumulative
        error, scale = max(error, abs(defect)), max(scale, abs(cumulative))
        residual = max(residual, float(abs(balance[4])))
        history.append([row['time'], physical, pending, cumulative, defect])
        previous = pressure
    assert scale > 0, 'a nonzero boundary-volume scale is required'
    return dict(frames=len(rows), common_end_s=rows[-1]['time'], maximum_content_defect_m3=error,
                max_net_boundary_volume_m3=scale, normalized_content_defect=error/scale,
                largest_reported_step_residual_m3=residual), np.asarray(history)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', type=Path)
    parser.add_argument('--until', type=float)
    parser.add_argument('--output', type=Path, help='New JSON report and matching CSV history')
    parser.add_argument('--measure', action='store_true', help='Report an old failing budget without treating it as qualified')
    args = parser.parse_args()
    report, history = analyze(args.config, args.until)
    report.update(passed=report['normalized_content_defect'] < 1e-9,
                  config_sha256=water.digest(args.config), analyzer_sha256=water.digest(__file__))
    if args.output:
        assert not args.output.exists(), args.output
        water.write_json(args.output, report)
        np.savetxt(args.output.with_suffix('.csv'), history, delimiter=',', comments='',
                   header='time_s,physical_content_m3,pending_content_m3,boundary_net_m3,defect_m3')
    print(json.dumps(report, indent=2))
    if not args.measure:
        assert report['passed'], 'physical content + pending reservoir does not match boundary exchange'


if __name__ == '__main__':
    main()
