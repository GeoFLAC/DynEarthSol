#!/usr/bin/env python3
"""Check effective-stress Mohr-Coulomb strength on an affine material path.

Normal extensions +e and -e give zero volume change and principal shear
stress t=2*G*e. At zero dilation the mean stress remains -C. The yield
load follows directly from the Mohr circle:
t_y = c*cos(phi) + (C - alpha*p)*sin(phi).
All nodes have prescribed velocities; no equilibrium or mesh error can
obscure the pressure level. EP-RSF is checked in its constant-friction
limit (a=b=0), independently of state/slip-rate accuracy.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, scan_frames
from des_case import digest, file_hashes, read_config, reset_boundaries, run_des, write_json

G, COHESION, PHI, STRAIN, STEPS = 1e9, 1e5, np.pi/6, .0005, 64
CASES = (
    ('dry-effective', 2e5, 0., 1.),
    ('pressurized', 1e6, 8e5, 1.),
    ('dry-total', 1e6, 0., 1.),
    ('biot', 1e6, 8e5, .6),
)


def config(ndims, rheology, confinement, pressure, biot):
    cfg = read_config(ROOT/'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='cell', max_steps=str(STEPS), max_time_in_yr='1',
                      output_step_interval='4', output_time_interval_in_yr='1e9',
                      checkpoint_frame_interval='1', has_initial_checkpoint='true')
    cfg['mesh'].update(meshing_elem_shape='1', resolution='1', largest_size='1', smallest_size='1')
    cfg['control'].update(fixed_dt='1', has_PT='true', PT_option='1', has_moving_mesh='false',
                          has_poroelastic_pressure_feedback='false', is_using_mixed_stress='false',
                          gravity='0', has_thermal_diffusion='false', PT_max_iter='100000')
    bc = cfg['bc']
    reset_boundaries(bc, vbc='1', vbc_val='0', hbc='0')
    bc.update(vbc_val_x1=repr(STRAIN/STEPS), vbc_val_z1=repr(-STRAIN/STEPS),
              stress_bc_z1='0', hbc_pressure='[0,0,0,0,0,0]')
    nstr = ndims*(ndims+1)//2
    initial = [-confinement]*ndims + [0.]*(nstr-ndims)
    cfg['ic'].update(weakzone_option='0', initial_stress_option='1', initial_stress=json.dumps(initial),
                      pore_pressure_option='2', pore_pressure_value=repr(pressure), has_body_force_adjustment='false')
    cfg['mat'].update(rheology_type=rheology, bulk_modulus='[1e9]', shear_modulus='[1e9]',
                      alpha='[0]', hydraulic_perm='[0]', biot_coeff=f'[{biot}]', is_plane_strain=str(ndims == 2).lower(),
                      porosity='[0.3]', fluid_bulk_modulus='[2.2e9]', bulk_modulus_s='[1e20]', fluid_visc='[1e-3]',
                      cohesion0='[1e5]', cohesion1='[1e5]', friction_angle0='[30]', friction_angle1='[30]',
                      dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e9',
                      direct_a='[0]', evolution_b='[0]', state_var_model='1')
    return cfg


def compare(folder, ndims, confinement, pressure, biot):
    model = Dynearthsol(str(folder/'cell'))
    assert model.nnode_list[0] == 2**ndims
    q = biot*pressure
    yield_load = float(COHESION*np.cos(PHI)+(confinement-q)*np.sin(PHI))
    worst, final_load = 0., None
    for row in scan_frames(str(folder/'cell')):
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
                expected_yield_load_Pa=yield_load, final_load_Pa=final_load, effective_pressure_Pa=q,
                frames=len(model.frames))


def study(executable, output, ndims, show_failure):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = digest(executable)
    sources = (Path(__file__),)
    for source in sources: shutil.copy2(source, output/source.name)
    report = dict(executable=str(executable), executable_sha256=binary_hash, ndims=ndims,
                  source_sha256=file_hashes(sources), cases={})
    for rheology in ('elasto-plastic', 'elasto-plastic-rate-state-friction'):
        for name, confinement, pressure, biot in CASES:
            folder = output/f'{rheology}-{name}'; folder.mkdir()
            run_des(executable, folder, 'input', config(ndims, rheology, confinement, pressure, biot))
            case = compare(folder, ndims, confinement, pressure, biot)
            report['cases'][folder.name] = case
            print(folder.name, json.dumps(case), flush=True)
    assert digest(executable) == binary_hash
    report['passed'] = all(c['passed'] for c in report['cases'].values())
    report['files_sha256'] = file_hashes((p for p in output.rglob('*')
                                          if p.suffix in ('.cfg', '.info') or '.save.' in p.name or '.chkpt.' in p.name), output)
    write_json(output/'qualification.json', report)
    assert report['passed'] or show_failure, 'Pressure-dependent strength fails the reference'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--ndims', type=int, choices=(2, 3), default=2)
    parser.add_argument('--show-failure', action='store_true', help='Record a failed qualification without accepting it')
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.ndims, args.show_failure)
