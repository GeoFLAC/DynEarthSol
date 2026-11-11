#!/usr/bin/env python3
"""Check scheduled flux/pressure, discrete balance and exact restart in DES.

Usage: check_boundaries.py EXECUTABLE OUTPUT_DIRECTORY
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames


from des_case import read_config, run_solver

def mesh_tags(exe, output):
    """Imported tag names do not determine which pressure faces intersect."""
    output.mkdir(parents=True, exist_ok=False)
    result = {}
    for name, top, right, active, values, expected in (
            ('disjoint', 1, 32, (1, 5), (0, 1000, 0, 0, 0, 0), 0),
            ('shared', 32, 32, (1, 5), (0, 1000, 0, 0, 0, 0), 11),
            ('shared-axis', 1, 1, (0, 1), (0, 1000, 0, 0, 0, 0), 11),
            ('equal-shared', 32, 32, (1, 5), (0, 1000, 0, 0, 0, 1000), 0)):
        folder = output/name; folder.mkdir()
        poly = folder/'boundary.poly'
        poly.write_text(f'''4 2 0 0
0 0 -1
1 1 -1
2 1 0
3 0 0
4 1
0 0 1 16
1 1 2 {right}
2 2 3 {top}
3 3 0 2
0
1
0 .5 -.5 0 1
''')
        cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
        cfg['sim'].update(modelname='tags', max_steps='2', output_step_interval='1',
                          has_initial_checkpoint='yes', checkpoint_frame_interval='1')
        cfg['mesh'].update(meshing_option='90', poly_filename=str(poly), xlength='1',
                           ylength='1', zlength='1', resolution='1', largest_size='1',
                           smallest_size='1', remeshing_option='0')
        cfg['control'].update(has_PT='false', is_quasi_static='true', gravity='0',
                              fixed_dt='.01', has_thermal_diffusion='false', has_moving_mesh='false',
                              has_poroelastic_pressure_feedback='false',
                              inertial_scaling='1e5', characteristic_speed='1e-5')
        for i, side in enumerate(('x0', 'x1', 'y0', 'y1', 'z0', 'z1')):
            cfg['bc'].update({f'hbc_{side}': '2' if i in active else '0',
                              f'vbc_{side}': '3', f'vbc_val_{side}': '0'})
        cfg['bc'].update(hbc_pressure=json.dumps(values), stress_bc_z1='0')
        cfg['ic'].update(pore_pressure_option='2', pore_pressure_value='0',
                         has_body_force_adjustment='false', weakzone_option='0')
        with (folder/'input.cfg').open('w') as stream: cfg.write(stream)
        done = run_solver([str(exe), 'input.cfg'], log=folder / 'run.log', cwd=folder, timeout=60, env=dict(os.environ, OMP_NUM_THREADS='1'))
        result[name] = dict(returncode=done.returncode, expected=expected)
        if done.returncode == 0:
            model = Dynearthsol(str(folder/'tags')); x = model.read_field(0, 'coordinate')
            assert model.ndims == 2
            for frame in model.frames:
                p = model.read_field(frame, 'pore pressure')
                if name in ('disjoint', 'equal-shared'):
                    np.testing.assert_array_equal(p[np.isclose(x[:, 0], 0)], 1000)
                    np.testing.assert_array_equal(p[np.isclose(x[:, 0], 1)], values[5])
        elif done.returncode == 11:
            result[name]['actual_node_diagnostic'] = 'Node ' in (folder/'run.log').read_text()
    # A changed condition on restart must be checked against the saved mesh.
    cfg.read(output/'disjoint/input.cfg')
    cfg['sim'].update(modelname='bad-restart', is_restarting='yes',
                      restarting_from_modelname=str(output/'disjoint/tags'), restarting_from_frame='1', max_steps='3')
    cfg['bc']['hbc_x0'] = '2'
    with (output/'restart.cfg').open('w') as stream: cfg.write(stream)
    done = run_solver([str(exe), 'restart.cfg'], log=output / 'restart.log', cwd=output, timeout=60, env=dict(os.environ, OMP_NUM_THREADS='1'))
    result['restart-shared'] = dict(returncode=done.returncode, expected=11,
                                    actual_node_diagnostic='Node ' in (output/'restart.log').read_text())
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    assert all(c['returncode'] == c['expected'] for c in result.values()), result
    assert all(c.get('actual_node_diagnostic', True) for c in result.values()), result
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--mesh-tags-only', action='store_true', help='2D imported pressure-tag conflicts')
    args = parser.parse_args()
    exe, output = args.executable.resolve(), args.output.resolve()
    if args.mesh_tags_only:
        mesh_tags(exe, output)
        return
    output.mkdir(parents=True, exist_ok=False)
    result = {}
    for mode in ('flux', 'flux-side', 'pressure'):
        cfg = read_config(ROOT / 'tests/functional/pt_remesh_restart.cfg')
        cfg['control'].update(gravity='0', has_moving_mesh='false', fixed_dt='0.01',
                              PT_absolute_tolerance='1e-8')
        cfg['bc'].update(stress_val_z1='0', hbc_z1='2' if mode == 'pressure' else '3',
                         hbc_flux='[0,0,0,0,0,-1e-7]', hbc_pressure='[0,0,0,0,0,1000]',
                         num_hbc_schedule='4',
                         hbc_schedule_time_in_yr='[0, 3e-10, 3e-10, 6e-10]',
                         hbc_schedule_ratio='[1,1,-0.5,0]')
        if mode == 'flux-side':
            cfg['bc'].update(hbc_z1='0', hbc_x0='3', hbc_flux='[-1e-7,0,0,0,0,0]')
        cfg['mat']['permeability_tensor'] = '[2,1,1,0,0.3,0]'
        source = output / f'{mode}.cfg'
        with source.open('w') as f:
            cfg.write(f)
        case = output / mode
        subprocess.run([sys.executable, str(ROOT/'tests/functional/pt_initial_checkpoint.py'),
                        str(exe), str(source), str(case), '--restart-frame', '1', '--steps', '4'],
                       env=dict(os.environ, OMP_NUM_THREADS='1'), check=True, timeout=240)
        model = Dynearthsol(str(case/'fresh'))
        worst = 0.0
        times = {r['frame']: r['time'] for r in scan_frames(str(case/'fresh'))}
        previous = 0.0
        for frame in model.frames[1:]:
            b = model.read_field(frame, 'fluid balance')
            scale = max(np.abs(b[:4]).max(), 1e-30)
            worst = max(worst, abs(b[4])/scale)
            assert abs(b[4]) <= 1e-9*scale + 1e-25, b
            a, stop = 3e-10*YEAR2SEC, 6e-10*YEAR2SEC
            if mode.startswith('flux'):
                def integral(t):
                    u = min(max(t-a, 0), stop-a)
                    return min(t,a) - 0.5*u + 0.25*u*u/(stop-a)
                area = 10 if mode == 'flux-side' else 1
                imposed = 1e-7*area*(integral(times[frame])-integral(previous))
                np.testing.assert_allclose(b[3], imposed, rtol=1e-12, atol=1e-20)
                np.testing.assert_allclose(b[0]+b[1], imposed, rtol=1e-9, atol=1e-20)
            else:
                xyz = model.read_field(frame, 'coordinate')
                top = np.isclose(xyz[:, -1], 0)
                p = model.read_field(frame, 'pore pressure')[top]
                t = times[frame]
                a, b = 3e-10*YEAR2SEC, 6e-10*YEAR2SEC
                ratio = 1 if t < a else (-0.5*(b-t)/(b-a) if t < b else 0)
                np.testing.assert_allclose(p, 1000*ratio, atol=1e-10, rtol=1e-12)
            previous = times[frame]
        result[mode] = dict(max_relative_step_balance_error=worst, exact_restart=True)
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    if model.ndims == 2:
        mesh_tags(exe, output/'mesh-tags')


if __name__ == '__main__':
    main()
