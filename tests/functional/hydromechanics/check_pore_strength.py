#!/usr/bin/env python3
"""Check effective-stress Mohr-Coulomb strength on an affine material path.

Normal extensions +e and -e give zero volume change and principal shear
stress t=2*G*e. At zero dilation the mean stress remains -C. The yield
load follows directly from the Mohr circle:
t_y = c*cos(phi) + (C - alpha*Se*p)*sin(phi).
All nodes have prescribed velocities; no equilibrium or mesh error can
obscure the pressure level. EP-RSF is checked in its constant-friction
limit (a=b=0), independently of state/slip-rate accuracy.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

import check_bishop_cell as cell
from unsaturated_reference import Retention

G, COHESION, PHI, STRAIN, STEPS = 1e9, 1e5, np.pi/6, .0005, 64
CASES = (
    ('dry-effective', 2e5, 0., 1., 0),
    ('pressurized', 1e6, 8e5, 1., 0),
    ('dry-total', 1e6, 0., 1., 0),
    ('biot', 1e6, 8e5, .6, 0),
    ('suction', 2e5, -2e5, 1., 3),
)


def config(ndims, rheology, confinement, pressure, biot, retention):
    cfg = cell.config(ndims, STEPS, 'moderate', 1)
    cfg['sim'].update(max_steps=str(STEPS), max_time_in_yr='1', output_step_interval='4',
                      checkpoint_frame_interval='1', has_initial_checkpoint='true')
    cfg['mesh'].update(resolution='1', largest_size='1', smallest_size='1')
    cfg['control'].update(fixed_dt='1', has_PT='true', PT_option='1', has_moving_mesh='false',
                          has_poroelastic_pressure_feedback='false', is_using_mixed_stress='false')
    cfg['bc'].update(vbc_val_x1=repr(STRAIN/STEPS), vbc_val_z1=repr(-STRAIN/STEPS))
    nstr = ndims*(ndims+1)//2
    initial = [-confinement]*ndims + [0.]*(nstr-ndims)
    cfg['ic'].update(weakzone_option='0', initial_stress_option='1', initial_stress=json.dumps(initial),
                      pore_pressure_value=repr(pressure))
    cfg['mat'].update(rheology_type=rheology, bulk_modulus='[1e9]', shear_modulus='[1e9]',
                      alpha='[0]', hydraulic_perm='[0]', retention_model=f'[{retention}]',
                      retention_alpha='[1e-5]', residual_saturation='[.1]',
                      biot_coeff=f'[{biot}]', is_plane_strain=str(ndims == 2).lower(),
                      cohesion0='[1e5]', cohesion1='[1e5]', friction_angle0='[30]', friction_angle1='[30]',
                      dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e9',
                      direct_a='[0]', evolution_b='[0]', state_var_model='1')
    return cfg


def compare(folder, ndims, confinement, pressure, biot, retention):
    model = cell.Dynearthsol(str(folder/'cell'))
    assert model.nnode_list[0] == 2**ndims
    se = float(Retention('gardner', 1e-5, sr=.1).se(pressure)) if retention else 1.
    q = biot*se*pressure
    yield_load = float(COHESION*np.cos(PHI)+(confinement-q)*np.sin(PHI))
    worst, final_load = 0., None
    for row in cell.scan_frames(str(folder/'cell')):
        frame = row['frame']; strain = model.read_field(frame, 'strain')
        imposed = STRAIN*row['time']/STEPS
        target = np.zeros(strain.shape[1]); target[0] = imposed; target[ndims-1] = -imposed
        np.testing.assert_allclose(strain, np.tile(target, (len(strain), 1)), atol=1e-16)
        np.testing.assert_allclose(model.read_field(frame, 'pore pressure'), pressure, rtol=0, atol=1e-8)
        load = min(2*G*imposed, yield_load)
        expected = np.zeros_like(strain); expected[:, :ndims] = -confinement
        expected[:, 0] += load; expected[:, ndims-1] -= load
        actual = model.read_field(frame, 'stress')
        worst = max(worst, float(np.max(np.abs(actual-expected)))/yield_load)
        plastic = max(0., imposed-yield_load/(2*G))
        worst = max(worst, float(np.max(np.abs(model.read_field(frame, 'plastic strain')-plastic)))*2*G/yield_load)
        final_load = float(np.mean((actual[:, :ndims].max(axis=1)-actual[:, :ndims].min(axis=1))/2))
        if 'rate-state' in folder.name:
            np.testing.assert_allclose(model.read_field(frame, 'dynamic friction coefficient'), np.tan(PHI), rtol=1e-14)
    return dict(passed=worst < 1e-10, max_error_over_yield_load=worst,
                expected_yield_load_Pa=yield_load, final_load_Pa=final_load, bishop_pressure_Pa=q,
                frames=len(model.frames))


def study(executable, output, ndims, show_failure):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = cell.digest(executable)
    sources = (Path(__file__), Path(cell.__file__), cell.ROOT/'benchmarks/poromechanics/unsaturated_reference.py')
    for source in sources: shutil.copy2(source, output/source.name)
    report = dict(executable=str(executable), executable_sha256=binary_hash, ndims=ndims,
                  source_sha256={p.name:cell.digest(p) for p in sources}, cases={})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for rheology in ('elasto-plastic', 'elasto-plastic-rate-state-friction'):
        for name, confinement, pressure, biot, retention in CASES:
            folder = output/f'{rheology}-{name}'; folder.mkdir()
            cfg = config(ndims, rheology, confinement, pressure, biot, retention)
            with (folder/'input.cfg').open('w') as stream: cfg.write(stream)
            with (folder/'run.log').open('w') as stream:
                subprocess.run([str(executable), 'input.cfg'], cwd=folder, env=env,
                               stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=120)
            case = compare(folder, ndims, confinement, pressure, biot, retention)
            report['cases'][folder.name] = case
            print(folder.name, json.dumps(case), flush=True)
    assert cell.digest(executable) == binary_hash
    report['passed'] = all(c['passed'] for c in report['cases'].values())
    report['files_sha256'] = {str(p.relative_to(output)):cell.digest(p) for p in sorted(output.rglob('*'))
                              if p.is_file() and (p.suffix in ('.cfg', '.info') or '.save.' in p.name or '.chkpt.' in p.name)}
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert report['passed'] or show_failure, 'Pressure-dependent strength fails the reference'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--show-failure', action='store_true', help='Record a failed qualification without accepting it')
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.ndims, args.show_failure)
