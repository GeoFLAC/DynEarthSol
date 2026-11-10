#!/usr/bin/env python3
"""Check Bishop history when hydraulic substeps enter suction from saturation.

A sealed affine cell expands at a fixed rate. Nonzero pressure history is
consumed while wet, then pressure crosses zero inside a physical step.
Two permeabilities refine hydraulic subcycling without changing uniform
flow. Compare each physical interval with an independent scalar ODE for
the positive fixed-stress stabilizer and the actual preceding increment.
High residual saturation resolves the transition at small imposed strain.
"""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

import check_bishop_cell as cell
sys.path.insert(0, str(cell.ROOT/'tests/functional'))
from pt_initial_checkpoint import fields

DT, RATE, STEPS = .1, 1e-4, 24
PHI, RESIDUAL, ALPHA = .3, .999, 1.
S0 = PHI/2.2e9+(1-PHI)/37e9


def macro_reference(pressure, previous_increment, stabilizer):
    """Integrate (SB+phi*S') p_dot = -chi*eps_dot + L*D/dt."""
    storage = S0+stabilizer
    forcing = stabilizer*previous_increment/DT
    span = DT
    if pressure > 0:
        velocity = (-RATE+forcing)/storage
        if pressure+velocity*DT >= 0:
            return pressure+velocity*DT
        span = DT+pressure/velocity
        pressure = 0.

    def rhs(time, value):
        se = math.exp(ALPHA*min(float(value[0]), 0.))
        return [(-RATE*se+forcing)/(storage+PHI*(1-RESIDUAL)*ALPHA*se)]

    solution = solve_ivp(rhs, (0., span), [pressure], method='DOP853', rtol=1e-12, atol=1e-13)
    assert solution.success
    return float(solution.y[0, -1])


def substeps(coordinate, connectivity, permeability, stabilizer):
    """Independent simplex mass and absolute-row operator bound."""
    dim = coordinate.shape[1]
    vertices = np.concatenate((np.ones((*connectivity.shape, 1)), coordinate[connectivity]), axis=2)
    volume = np.abs(np.linalg.det(vertices))/math.factorial(dim)
    gradient = np.linalg.inv(vertices)[:, 1:, :].swapaxes(1, 2)
    mass, operator = np.zeros(len(coordinate)), np.zeros(len(coordinate))
    for e, nodes in enumerate(connectivity):
        mass[nodes] += (S0+stabilizer)*volume[e]/(dim+1)
        operator[nodes] += permeability/.001*volume[e]*np.abs(gradient[e]@gradient[e].T).sum(axis=1)
    return math.ceil(DT/float(np.min(mass/operator)))


def study(executable, output, ndims):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = cell.digest(executable)
    stabilizer = 1/(1e4+(1e4/3 if ndims == 2 else 0))
    report = dict(executable=str(executable), executable_sha256=binary_hash, ndims=ndims,
                  reference=dict(dt_s=DT, expansion_rate_per_s=RATE, steps=STEPS, initial_pressure_Pa=1.,
                                 porosity=PHI, residual_saturation=RESIDUAL, retention_alpha_per_Pa=ALPHA,
                                 inverse_biot_modulus_per_Pa=S0, stabilizer_per_Pa=stabilizer), cases=[])
    for source in (Path(__file__), Path(cell.__file__)):
        shutil.copy2(source, output/source.name)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for permeability in (1e-5, 1e-4):
        folder = output/f'k{permeability:g}'; folder.mkdir()
        cfg = cell.config(ndims, STEPS, 'moderate', 1)
        cfg['sim'].update(max_steps=str(STEPS), max_time_in_yr='1', output_step_interval='1', checkpoint_frame_interval='1')
        cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
        cfg['control'].update(fixed_dt=str(DT), has_PT='true', PT_option='1', has_moving_mesh='false')
        cfg['bc'].update(vbc_val_z1='0', vbc_val_x1=str(RATE))
        cfg['ic']['pore_pressure_value'] = '1'
        cfg['mat'].update(bulk_modulus='[1e4]', shear_modulus='[1e4]', bulk_modulus_s='[37e9]',
                          fluid_bulk_modulus='[2.2e9]', porosity=f'[{PHI}]', alpha='[0]', retention_alpha=f'[{ALPHA}]',
                          residual_saturation=f'[{RESIDUAL}]', hydraulic_perm=f'[{permeability}]',
                          is_plane_strain=str(ndims == 2).lower())
        with (folder/'input.cfg').open('w') as stream: cfg.write(stream)
        with (folder/'run.log').open('w') as stream:
            subprocess.run([str(executable), 'input.cfg'], cwd=folder, env=env, check=True,
                           stdout=stream, stderr=subprocess.STDOUT, timeout=120)
        model = cell.Dynearthsol(str(folder/'cell'))
        assert model.ndims == ndims and model.nnode_list[0] == 2**ndims
        pressures, history = [], []
        for row in cell.scan_frames(str(folder/'cell')):
            pressure = model.read_field(row['frame'], 'pore pressure')
            assert np.ptp(pressure) < 1e-8
            strain = model.read_field(row['frame'], 'strain')[:, :ndims].sum(axis=1)
            np.testing.assert_allclose(strain, RATE*row['time'], rtol=1e-10, atol=1e-14)
            pressures.append(float(pressure.mean()))
            if len(pressures) < 2: continue
            previous = pressures[-2]-pressures[-3] if len(pressures) > 2 else 0.
            exact = macro_reference(pressures[-2], previous, stabilizer)
            history.append(dict(step=row['steps'], time_s=row['time'], start_Pa=pressures[-2],
                                previous_increment_Pa=previous, pressure_Pa=pressures[-1],
                                reference_Pa=exact, error_Pa=pressures[-1]-exact))
        entry = next(r for r in history if r['start_Pa'] > 0 > r['pressure_Pa'])
        assert entry['previous_increment_Pa'] != 0
        before = next(r for r in cell.scan_frames(str(folder/'cell')) if r['steps'] == entry['step']-1)
        cfg['sim'].update(modelname='resumed', is_restarting='true', restarting_from_modelname=str(folder/'cell'),
                          restarting_from_frame=str(before['frame']))
        with (folder/'restart.cfg').open('w') as stream: cfg.write(stream)
        with (folder/'restart.log').open('w') as stream:
            subprocess.run([str(executable), 'restart.cfg'], cwd=folder, env=env, check=True,
                           stdout=stream, stderr=subprocess.STDOUT, timeout=120)
        fresh = {r['steps']:r for r in cell.scan_frames(str(folder/'cell'))}
        checked = 0
        for row in cell.scan_frames(str(folder/'resumed')):
            if row['steps'] == before['steps']: continue
            for suffix in ('save', 'chkpt'):
                a = fields(folder/f"cell.{suffix}.{fresh[row['steps']]['frame']:06d}")
                b = fields(folder/f"resumed.{suffix}.{row['frame']:06d}")
                assert a.keys() == b.keys()
                for key in a:
                    if key in ('provenance', 'walltime_sec'): continue
                    assert a[key] == b[key], (row['steps'], suffix, key)
                    checked += 1
        count = substeps(model.read_field(0, 'coordinate'), model.read_field(0, 'connectivity'), permeability, stabilizer)
        case = dict(permeability_m2=permeability, independent_substeps=count, entry=entry,
                    max_error_Pa=max(abs(r['error_Pa']) for r in history), history=history,
                    restart_before_entry_frame=before['frame'], exact_restart_payloads=checked,
                    input_sha256=cell.digest(folder/'input.cfg'),
                    restart_input_sha256=cell.digest(folder/'restart.cfg'),
                    outputs_sha256={p.name:cell.digest(p) for p in sorted(folder.iterdir())
                                    if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.info'})
        report['cases'].append(case)
        print(permeability, count, entry, flush=True)
    assert cell.digest(executable) == binary_hash
    coarse, fine = report['cases']
    report['passed'] = bool(fine['independent_substeps'] >= 9*coarse['independent_substeps'] and
                            fine['max_error_Pa'] < 1e-5 and fine['max_error_Pa'] < .25*coarse['max_error_Pa'])
    report['runner_sha256'] = cell.digest(Path(__file__))
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert report['passed'], 'Wet-transition history does not converge under hydraulic subcycling'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.ndims)
