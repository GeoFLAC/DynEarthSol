#!/usr/bin/env python3
"""Affine plane-strain material paths sampled from the Salencon solution.

All four nodes have prescribed velocities. Thus no free mechanical residual
or mesh resolution can obscure the constitutive update. Large and refined
increments are reported separately; a small-increment pass does not qualify
the large-increment return or the spatial tunnel problem.
"""

import argparse
import json
import os
from pathlib import Path
import shutil

import numpy as np

import check_salencon as cs


from des_case import run_solver

def patch_config(poly):
    poly.write_text('4 2 0 0\n0 0 0\n1 1 0\n2 1 -1\n3 0 -1\n'
                    '4 1\n0 0 1 16\n1 1 2 2\n2 2 3 32\n3 3 0 1\n0\n1\n0 .5 -.5 0 2\n')
    cfg = cs.config(1.)
    cfg['mesh'].update(poly_filename=str(poly), xlength='1', zlength='1', resolution='1')
    cfg['control'].update(PT_option='1', is_using_mixed_stress='false')
    for side in ('x0', 'x1', 'z0', 'z1'):
        cfg['bc'].update({f'vbc_{side}':'1', f'vbc_val_{side}':'0',
                          f'stress_bc_{side}':'0', f'stress_val_{side}':'0'})
    return cfg


def run(executable, cfg, folder, name):
    with (folder/f'{name}.cfg').open('w') as stream: cfg.write(stream)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    run_solver([str(executable), f'{name}.cfg'], log=folder / f'{name}.log', cwd=folder, env=env, check=True, timeout=90)


def study(executable, output, manifest):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = cs.ck.digest(executable)
    solver = json.loads(manifest.read_text()) if manifest else None
    if solver:
        cs.ck.check_binary(solver, binary_hash)
        shutil.copy2(manifest, output/'solver-qualification.json')
    sources = [Path(__file__), Path(cs.__file__), Path(cs.ck.__file__),
               cs.ROOT/'benchmarks/poromechanics/salencon_reference.py']
    for path in sources: shutil.copy2(path, output/path.name)
    report = dict(executable=str(executable), executable_sha256=binary_hash,
                  solver_commit=solver['implementation_commit'] if solver else None,
                  source_sha256={p.name:cs.ck.digest(p) for p in sources}, cases=[])
    for dilation, radius in ((0., 1.05), (30., 1.1)):
        group = output/f'psi{dilation:g}'; group.mkdir()
        cfg = patch_config(group/'patch.poly')
        cfg['mat'].update(dilation_angle0=f'[{dilation}]', dilation_angle1=f'[{dilation}]')
        run(executable, cfg, group, 'seed')
        target = cs.state(radius, cs.A, cs.R, cs.P, cs.C, cs.F, cs.K, cs.G, dilation=dilation)
        expected = np.array([target[0], target[1], 0.])
        for steps in (1, 16, 64):
            folder = group/f'n{steps}'; folder.mkdir()
            cfg['sim'].update(modelname='patch', max_steps=str(steps), is_restarting='true',
                              restarting_from_modelname=str(group/'seed'), restarting_from_frame='0')
            cfg['bc'].update(vbc_val_x1=str(float(target[4])/steps), vbc_val_z1=str(-float(target[5])/steps))
            run(executable, cfg, folder, 'patch')
            model = cs.ck.Dynearthsol(str(folder/'patch')); frame = model.frames[-1]
            assert len(model.read_field(frame, 'coordinate')) == 4
            strain = model.read_field(frame, 'strain')
            np.testing.assert_allclose(strain, np.tile([target[4], target[5], 0.], (2, 1)), rtol=1e-12, atol=1e-18)
            stress = model.read_field(frame, 'stress')
            syy = np.frombuffer(cs.ck.fields(folder/f'patch.chkpt.{frame:06d}')['stressyy'], dtype=np.float64)
            error = float(max(np.max(np.abs(stress-expected)), np.max(np.abs(syy-target[2])))/cs.P)
            case = dict(dilation_degrees=dilation, reference_radius_m=radius, load_steps=steps,
                        error_over_P=error, expected_stress_Pa=expected.tolist(),
                        expected_syy_Pa=float(target[2]), actual_stress_Pa=stress.tolist(), actual_syy_Pa=syy.tolist())
            report['cases'].append(case)
            print(dilation, steps, error, flush=True)
    assert cs.ck.digest(executable) == binary_hash
    report['passed_for_small_increments'] = all(c['error_over_P'] < 1e-10 for c in report['cases'] if c['load_steps'] > 1)
    report['large_increment_qualified'] = all(c['error_over_P'] < 1e-10 for c in report['cases'] if c['load_steps'] == 1)
    report['files_sha256'] = {str(p.relative_to(output)):cs.ck.digest(p) for p in sorted(output.rglob('*'))
                              if p.is_file() and ('.save.' in p.name or '.chkpt.' in p.name or p.suffix in ('.cfg', '.poly', '.info'))}
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert report['passed_for_small_increments'], 'Refined constitutive path fails the analytical state'
    assert report['large_increment_qualified'], 'Uncut shear return fails the large-increment analytical state'


def edge_state(dilation, extension):
    """Two active faces at equal in-plane maximum principals; plane strain."""
    lam = cs.K - 2*cs.G/3
    stiffness = lam*np.ones((3, 3)) + 2*cs.G*np.eye(3)
    trial = -cs.P + stiffness@np.array([extension, extension, 0.])
    factor = lambda angle: (1 + np.sin(np.deg2rad(angle)))/(1 - np.sin(np.deg2rad(angle)))
    nphi, npsi = factor(cs.F), factor(dilation)
    normals = np.array([[nphi, 0., -1.], [0., nphi, -1.]])
    flow = np.array([[npsi, 0., -1.], [0., npsi, -1.]])
    multipliers = np.linalg.solve(normals@stiffness@flow.T,
                                   normals@trial - 2*cs.C*np.sqrt(nphi))
    assert np.all(multipliers >= 0)
    result = trial - stiffness@flow.T@multipliers
    assert result[2] < min(result[:2])
    np.testing.assert_allclose(normals@result, 2*cs.C*np.sqrt(nphi), rtol=1e-14)
    np.testing.assert_allclose(result[0], result[1], rtol=1e-14)
    return result


def edges(executable, output, manifest, show_failure):
    """Isolate the edge return without smoothing or free mechanical nodes."""
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = cs.ck.digest(executable)
    solver = json.loads(manifest.read_text()) if manifest else None
    if solver:
        cs.ck.check_binary(solver, binary_hash)
        shutil.copy2(manifest, output/'solver-qualification.json')
    sources = [Path(__file__), Path(cs.__file__), Path(cs.ck.__file__)]
    for source in sources: shutil.copy2(source, output/source.name)
    report = dict(executable=str(executable), executable_sha256=binary_hash,
                  solver_commit=solver['implementation_commit'] if solver else None,
                  source_sha256={p.name:cs.ck.digest(p) for p in sources}, cases=[])
    cfg = patch_config(output/'patch.poly')
    run(executable, cfg, output, 'seed')
    extension = 1.2e-5
    for dilation in (0., 30.):
        target = edge_state(dilation, extension)
        for steps in (1, 16, 64):
            folder = output/f'psi{dilation:g}-n{steps}'; folder.mkdir()
            cfg['mat'].update(dilation_angle0=f'[{dilation}]', dilation_angle1=f'[{dilation}]')
            cfg['sim'].update(modelname='patch', max_steps=str(steps), is_restarting='true',
                              restarting_from_modelname=str(output/'seed'), restarting_from_frame='0')
            cfg['bc'].update(vbc_val_x1=str(extension/steps), vbc_val_z1=str(-extension/steps))
            run(executable, cfg, folder, 'patch')
            model = cs.ck.Dynearthsol(str(folder/'patch')); frame = model.frames[-1]
            assert len(model.read_field(frame, 'coordinate')) == 4
            actual = model.read_field(frame, 'stress')
            np.testing.assert_allclose(model.read_field(frame, 'strain'),
                                       np.tile([extension, extension, 0.], (len(actual), 1)), atol=1e-18, rtol=1e-12)
            syy = np.frombuffer(cs.ck.fields(folder/f'patch.chkpt.{frame:06d}')['stressyy'], dtype=np.float64)
            error = float(max(np.max(np.abs(actual-[target[0], target[1], 0.])),
                              np.max(np.abs(syy-target[2])))/cs.P)
            case = dict(dilation_degrees=dilation, load_steps=steps, error_over_P=error,
                        expected_xx_zz_yy_Pa=target.tolist(), actual_stress_Pa=actual.tolist(),
                        actual_syy_Pa=syy.tolist(), yield_admissibility=cs.yield_admissibility(actual, syy))
            report['cases'].append(case)
            print(dilation, steps, error, case['yield_admissibility'], flush=True)
    assert cs.ck.digest(executable) == binary_hash
    report['passed'] = all(c['error_over_P'] < 1e-10 and
                          c['yield_admissibility']['outside_yield'] == 0 for c in report['cases'])
    report['scope'] = 'Two equal maximum principals; both dilation regimes; not a general edge/apex qualification.'
    report['files_sha256'] = {str(p.relative_to(output)):cs.ck.digest(p) for p in sorted(output.rglob('*'))
                              if p.is_file() and ('.save.' in p.name or '.chkpt.' in p.name or p.suffix in ('.cfg', '.poly', '.info'))}
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert report['passed'] or show_failure, 'MC edge return violates symmetry, yield and reference; retain evidence'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--binary-manifest', type=Path)
    parser.add_argument('--edges', action='store_true', help='Check the symmetric two-face return')
    parser.add_argument('--show-failure', action='store_true', help='Retain an edge failure without accepting it')
    args = parser.parse_args()
    if args.edges:
        edges(args.executable.resolve(), args.output.resolve(),
              args.binary_manifest.resolve() if args.binary_manifest else None, args.show_failure)
    else:
        study(args.executable.resolve(), args.output.resolve(), args.binary_manifest.resolve() if args.binary_manifest else None)
