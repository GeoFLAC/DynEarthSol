#!/usr/bin/env python3
"""Run and check the DES consolidation columns against the Biot reference."""
import argparse
import configparser
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import numpy as np
from terzaghi_reference import TerzaghiReference

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import read_config, run_solver, source_revision

def load_study(study, record_source, active_runner=None):
    """Verify the archived and optional active runner before recording a study."""
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    if active_runner is not None:
        assert hashlib.sha256(active_runner.read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def make_config(dim):
    cfg = read_config(ROOT / f'benchmarks/poromechanics/1d-consolidation-des{dim}d_traction.cfg')
    cfg['control'].update(PT_option='1', PT_absolute_tolerance='1e-6', gravity='10')
    return cfg


def run_cases(output, timeout):
    environment = dict(os.environ, OMP_NUM_THREADS='1')
    runs = {}
    for dim in (2, 3):
        executable = ROOT / f'dynearthsol{dim}d'
        case = output / f'terzaghi-{dim}d'
        case.mkdir()
        cfg = make_config(dim)
        with (case / 'input.cfg').open('w') as stream:
            cfg.write(stream)
        command = [str(executable), 'input.cfg']
        executable_sha256 = hashlib.sha256(executable.read_bytes()).hexdigest()
        print(f'Running {dim}D Terzaghi in {case}', flush=True)
        run_solver(command, log=case / 'run.log', cwd=case, env=environment, check=True, timeout=timeout)
        if hashlib.sha256(executable.read_bytes()).hexdigest() != executable_sha256:
            raise RuntimeError('Executable changed during the run; rebuild before benchmarking')
        runs[str(dim)] = dict(command=command, omp_threads=1,
                             executable_sha256=executable_sha256,
                             input_sha256=hashlib.sha256((case / 'input.cfg').read_bytes()).hexdigest())
    return runs


def reference_from_config(cfg):
    def material(name):
        return float(json.loads(cfg['mat'][name])[0])
    return TerzaghiReference(height=cfg.getfloat('mesh', 'zlength'),
                            load=-cfg.getfloat('bc', 'stress_val_z1'),
                            bulk=material('bulk_modulus'), shear=material('shear_modulus'),
                            porosity=material('porosity'), fluid_bulk=material('fluid_bulk_modulus'),
                            grain_bulk=material('bulk_modulus_s'), biot=material('biot_coeff'),
                            permeability=material('hydraulic_perm'), viscosity=material('fluid_visc'))


def load_case(output, dim):
    case = output / f'terzaghi-{dim}d'
    cfg = configparser.ConfigParser()
    cfg.read(case / 'input.cfg')
    model = Dynearthsol(str(case / cfg['sim']['modelname']))
    if model.ndims != dim or len(model.frames) < 2:
        raise RuntimeError(f'{dim}D case has incomplete or unexpected output')
    frame_times = {row['frame']: row['time']
                   for row in scan_frames(str(case / cfg['sim']['modelname']))}
    times = np.asarray([frame_times[frame] for frame in model.frames])
    expected_end = cfg.getfloat('sim', 'max_time_in_yr') * YEAR2SEC
    if times[-1] < expected_end - 1e-6:
        raise RuntimeError(f'{dim}D case stopped before its requested end time')
    reference = reference_from_config(cfg)
    coordinate = model.read_field(model.frames[0], 'coordinate')
    connectivity = model.read_field(model.frames[0], 'connectivity')
    cells = coordinate[connectivity]
    volume = np.abs(np.linalg.det(cells[:, 1:] - cells[:, :1])) / math.factorial(dim)
    area = np.sum(volume) / reference.height
    depth = -coordinate[:, -1]
    initial_strain = model.read_field(model.frames[0], 'strain')
    density = json.loads(cfg['mat']['fluid_rho0'])[0]
    fluid_weight = density * cfg.getfloat('control', 'gravity')
    excess, settlement = [], []
    for frame in model.frames:
        current_coordinate = model.read_field(frame, 'coordinate')
        np.testing.assert_array_equal(current_coordinate, coordinate)
        pressure = model.read_field(frame, 'pore pressure')
        # Compare perturbations about the same material reference hydrostat
        # used by the solver; these columns have one uniform fluid density.
        excess.append(pressure + fluid_weight * coordinate[:, -1])
        delta_strain = model.read_field(frame, 'strain') - initial_strain
        settlement.append(-np.sum(volume * delta_strain[:, dim - 1]) / area)
    excess = np.asarray(excess)
    settlement = np.asarray(settlement)
    exact = reference.pressure(depth, times)
    exact_settlement = reference.settlement(times)
    active = times > 0  # DES frame zero is pre-load, while the reference is t=0+.
    pressure_error = np.max(np.abs(excess - exact), axis=1) / reference.load
    settlement_error = np.abs(settlement - exact_settlement) / reference.drained_settlement
    truncation = np.max(np.abs(reference.pressure(depth, times[active], terms=256) - exact[active]))
    if not np.isfinite(excess).all() or not np.isfinite(settlement).all():
        raise RuntimeError(f'{dim}D output contains nonfinite values')
    if truncation > 1e-10 * reference.load:
        raise RuntimeError('Analytical series truncation exceeds its error budget')
    bottom = np.isclose(depth, reference.height)
    bottom_pressure = excess[:, bottom].mean(axis=1)
    np.savetxt(output / f'terzaghi-{dim}d-history.csv',
               np.column_stack((times, bottom_pressure, reference.pressure([reference.height], times)[:, 0],
                                settlement, exact_settlement, pressure_error, settlement_error)),
               delimiter=',', header='time_s,bottom_excess_Pa,reference_bottom_Pa,settlement_m,reference_settlement_m,max_pressure_error_over_load,settlement_error_over_drained', comments='')
    return dict(reference=reference, coordinate=coordinate, connectivity=connectivity,
                depth=depth, times=times, pressure=excess, settlement=settlement,
                pressure_error=pressure_error, settlement_error=settlement_error,
                metrics=dict(frames=len(times), final_time_s=float(times[-1]),
                             nodes=len(coordinate), elements=len(connectivity),
                             max_pressure_error_over_load=float(pressure_error[active].max()),
                             max_settlement_error_over_drained=float(settlement_error[active].max()),
                             final_settlement_mm=float(settlement[-1] * 1000),
                             reference_final_settlement_mm=float(exact_settlement[-1] * 1000),
                             series_256_vs_512_max_difference_Pa=float(truncation)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--analyze-only', '--plot-only', dest='analyze_only', action='store_true',
                        help='Check saved runs without executing DES')
    parser.add_argument('--timeout', type=float, default=180)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.analyze_only:
        if not output.is_dir():
            parser.error('--analyze-only requires an existing run directory')
        runs = json.loads((output/'summary.json').read_text()).get('runs', {})
    else:
        output.mkdir(parents=True, exist_ok=False)
        runs = run_cases(output, args.timeout)
    summary = dict(source_commit=source_revision(ROOT),
                   runs=runs, terzaghi={str(d): load_case(output, d)['metrics'] for d in (2, 3)})
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary['terzaghi'], indent=2))

if __name__ == '__main__':
    main()
