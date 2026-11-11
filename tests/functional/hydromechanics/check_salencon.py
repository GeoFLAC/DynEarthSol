#!/usr/bin/env python3
"""Dry plane-strain tunnel unloading against a finite Salencon annulus.

Start with uniform compression held by both walls; restart its equilibrated
frame zero and unload the inner wall in stages. Perfect Mohr-Coulomb
plasticity and small displacement isolate the tunnel's plastic response.
Associated and zero-dilation cases remain separate qualification regimes.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'benchmarks/poromechanics'))
import check_kirsch as ck
from salencon_reference import infinite_radius, plastic_radius, self_check, state

A, R, P, C, F, K, G = 1., 20., 5e5, 1e5, 30., 13e10/6, 1e10


def config(spacing):
    cfg = ck.config(spacing)
    cfg['sim']['modelname'] = 'seed'
    cfg['control']['PT_max_iter'] = '30000'
    cfg['bc'].update(stress_val_z1=str(-P), stress_bc_x1='4', stress_val_x1=str(-P))
    cfg['ic'].update(initial_stress_option='1', initial_stress=f'[-{P},-{P},0]', weakzone_option='0')
    cfg['mat'].update(rheology_type='elasto-plastic', bulk_modulus=f'[{K}]',
                      shear_modulus=f'[{G}]', cohesion0=f'[{C}]', cohesion1=f'[{C}]',
                      friction_angle0=f'[{F}]', friction_angle1=f'[{F}]',
                      dilation_angle0='[0]', dilation_angle1='[0]', max_tension='1e9')
    return cfg


def polar(values, normal):
    nx, nz = normal.T
    rr = nx**2*values[:, 0]+nz**2*values[:, 1]+2*nx*nz*values[:, 2]
    hoop = nz**2*values[:, 0]+nx**2*values[:, 1]-2*nx*nz*values[:, 2]
    shear = nx*nz*(values[:, 1]-values[:, 0])+(nx**2-nz**2)*values[:, 2]
    return rr, hoop, shear


def yield_admissibility(values, syy):
    """Mohr-Coulomb slack of the supplied stress, using sorted principals.

    After mixed-stress projection this is a diagnostic of the averaged
    stress, not a local constitutive-return admissibility test.
    """
    mean = (values[:, 0] + values[:, 1])/2
    radius = np.hypot((values[:, 0] - values[:, 1])/2, values[:, 2])
    principal = np.sort(np.array([mean-radius, mean+radius, syy]), axis=0)
    factor = (1 + np.sin(np.deg2rad(F)))/(1 - np.sin(np.deg2rad(F)))
    slack = principal[0] - factor*principal[2] + 2*C*np.sqrt(factor)
    assert np.isfinite(slack).all()
    return dict(elements=len(slack), outside_yield=int(np.sum(slack < -1e-8*P)),
                min_slack_Pa=float(slack.min()),
                max_violation_over_P=float(max(0., -slack.min())/P))


def compare(folder, dilation=0.):
    path = json.loads((folder/'load-path.json').read_text())
    model = ck.Dynearthsol(str(folder/'salencon'))
    seed = ck.Dynearthsol(str(folder/'seed'))
    initial, final = seed.frames[0], model.frames[-1]
    assert model.ndims == 2 and model.steps[-1] == len(path)
    x = seed.read_field(initial, 'coordinate'); conn = seed.read_field(initial, 'connectivity')
    np.testing.assert_array_equal(model.read_field(final, 'coordinate'), x)
    imposed = np.tile([-P, -P, 0.], (len(conn), 1))
    np.testing.assert_array_equal(seed.read_field(initial, 'stress'), imposed)
    np.testing.assert_array_equal(seed.read_field(initial, 'strain'), np.zeros_like(imposed))
    np.testing.assert_array_equal(seed.read_field(initial, 'plastic strain'), np.zeros(len(conn)))
    centre = x[conn].mean(axis=1); radius = np.linalg.norm(centre, axis=1)
    values = model.read_field(final, 'stress')
    rr, hoop, shear = polar(values, centre/radius[:, None])
    # The innermost straight facets have centroids infinitesimally inside r=a;
    # accuracy observations exclude the wall's geometric approximation.
    band = (radius >= 1.05*A) & (radius <= 5*A)
    exact = state(radius[band], A, R, P, C, F, K, G, dilation=dilation)
    er, et, _ = polar(model.read_field(final, 'strain'), centre/radius[:, None])
    syy = np.frombuffer(ck.fields(folder/f'salencon.chkpt.{final:06d}')['stressyy'], dtype=np.float64)
    errors = {}
    for name, measured, target, scale in (
            ('radial_stress', rr[band], exact[0], P),
            ('hoop_stress', hoop[band], exact[1], P),
            ('out_of_plane', syy[band], exact[2], P),
            ('radial_strain', er[band], exact[4], P/G),
            ('hoop_strain', et[band], exact[5], P/G)):
        errors[name+'_rms'] = float(np.sqrt(np.mean((measured-target)**2))/scale)
    errors['shear_stress_rms'] = float(np.sqrt(np.mean(shear[band]**2))/P)
    # PT accepted physical velocity is total displacement / physical dt; its
    # numerical iteration velocities are not sampled as physical histories.
    node_r = np.linalg.norm(x, axis=1)
    nodes = (node_r >= A-1e-12) & (node_r <= 5*A)
    displacement = np.zeros_like(x)
    for entry in path:
        stage = ck.Dynearthsol(str(folder/entry['model']))
        dt = stage.time[-1]-stage.time[0]
        displacement += stage.read_field(stage.frames[-1], 'velocity')*dt
    radial_u = np.sum(displacement*x/node_r[:, None], axis=1)
    exact_u = state(np.maximum(node_r[nodes], A), A, R, P, C, F, K, G, dilation=dilation)[3]
    wall_u = abs(float(state(A, A, R, P, C, F, K, G, dilation=dilation)[3]))
    errors['displacement_rms_over_wall'] = float(np.sqrt(np.mean((radial_u[nodes]-exact_u)**2))/wall_u)
    errors['displacement_max_over_wall'] = float(np.max(np.abs(radial_u[nodes]-exact_u))/wall_u)
    plastic = model.read_field(final, 'plastic strain')
    vertices = np.concatenate((np.ones((*conn.shape, 1)), x[conn]), axis=2)
    area = np.abs(np.linalg.det(vertices))/2
    plastic_area = float(area[plastic > 1e-8].sum())
    measured_b = float(np.sqrt(A*A+4*plastic_area/np.pi))
    b = plastic_radius(A, R, P, C, F)
    errors['plastic_radius_relative'] = abs(measured_b/b-1)
    assert np.isfinite(list(errors.values())).all()
    accepted = model.read_field(final, 'PT residual force')
    return dict(nodes=len(x), elements=len(conn), errors=errors,
                yield_admissibility=yield_admissibility(values[band], syy[band]),
                max_free_force_over_PR=ck.force_balance(x, conn, values, P, R),
                accepted_max_force_over_PR=float(np.max(np.linalg.norm(accepted, axis=1)))/(P*R),
                plastic_radius_m=measured_b, reference_plastic_radius_m=b,
                reference_wall_displacement_m=-wall_u,
                maximum_abs_strain=float(np.max(np.abs(model.read_field(final, 'strain')))))


def qualify(report):
    """Assess completed spatial results without conflating their stress samples."""
    fine = min(report['cases'], key=lambda c:c['spacing_m'])
    report['execution_completed'] = True
    report['force_qualified'] = bool(fine['max_free_force_over_PR'] < 1e-7)
    report['accuracy_qualified'] = bool(fine['errors']['radial_stress_rms'] < .02 and
                            fine['errors']['hoop_stress_rms'] < .02 and
                            fine['errors']['out_of_plane_rms'] < .02 and
                            fine['errors']['radial_strain_rms'] < .03 and
                            fine['errors']['hoop_strain_rms'] < .03 and
                            fine['errors']['displacement_rms_over_wall'] < .03 and
                            fine['errors']['displacement_max_over_wall'] < .1 and
                            fine['errors']['plastic_radius_relative'] < .05)
    # NMD changes the stress after the material return. Keep its yield slack
    # visible, but test local return admissibility in the affine material cases.
    mixed = report['mixed_stress']
    report['yield_check_applied'] = not mixed
    report['yield_check_scope'] = ('Diagnostic of saved mixed stress; local return tested separately.'
                                   if mixed else 'Saved stress without mixed-stress projection.')
    report['passed'] = bool(report['force_qualified'] and report['accuracy_qualified'] and
                            (mixed or fine['yield_admissibility']['max_violation_over_P'] < 1e-8))


def study(executable, output, spacings, manifest, pt_option=1, dilation=0., mixed=False, load_steps=16, pt_cfl=.25, refined=False):
    if load_steps < 1 or not spacings or min(spacings) <= 0:
        raise ValueError('Require positive load steps and mesh spacings')
    self_check()
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = ck.digest(executable)
    solver = json.loads(manifest.read_text()) if manifest else None
    if solver:
        ck.check_binary(solver, binary_hash)
        shutil.copy2(manifest, output/'solver-qualification.json')
    sources = [Path(__file__), Path(ck.__file__), ROOT/'benchmarks/poromechanics/salencon_reference.py']
    for source in sources: shutil.copy2(source, output/source.name)
    report = dict(executable=str(executable), executable_sha256=binary_hash,
                  solver_commit=solver['implementation_commit'] if solver else None,
                  source_sha256={p.name:ck.digest(p) for p in sources},
                  PT_option=pt_option, PT_CFL=pt_cfl, mixed_stress=mixed, load_steps=load_steps,
                  refined_radius_m=3*A if refined else None,
                  reference=dict(opening_radius_m=A, outer_radius_m=R, pressure_Pa=P,
                                 cohesion_Pa=C, friction_degrees=F, dilation_degrees=dilation,
                                 bulk_modulus_Pa=K, shear_modulus_Pa=G,
                                 finite_plastic_radius_m=plastic_radius(A, R, P, C, F),
                                 infinite_plastic_radius_m=float(infinite_radius(A, P, C, F))), cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for h in spacings:
        folder = output/f'h{h:g}'; folder.mkdir()
        ck.polygon(folder/'annulus.poly', h, refined_radius=3*A if refined else None)
        cfg = config(h)
        cfg['control']['PT_option'] = str(pt_option)
        cfg['control']['PT_CFL'] = str(pt_cfl)
        cfg['control']['PT_max_iter'] = '100000'
        cfg['control']['is_using_mixed_stress'] = str(mixed).lower()
        cfg['mat'].update(dilation_angle0=f'[{dilation}]', dilation_angle1=f'[{dilation}]')
        path = []
        for step in range(load_steps+1):
            phase = 'seed' if step == 0 else ('salencon' if step == load_steps else f'load-{step:03d}')
            if step:
                previous = 'seed' if step == 1 else path[-1]['model']
                cfg['sim'].update(modelname=phase, max_steps=str(step), is_restarting='true',
                                  restarting_from_modelname=str(folder/previous), restarting_from_frame=str(step-1))
                inner = -P*(1-step/load_steps)
                cfg['bc'].update(stress_bc_x1='4', stress_val_x1=str(inner))
            with (folder/f'{phase}.cfg').open('w') as stream: cfg.write(stream)
            with (folder/f'{phase}.log').open('w') as stream:
                try:
                    subprocess.run([str(executable), f'{phase}.cfg'], cwd=folder, env=env, check=True,
                                   stdout=stream, stderr=subprocess.STDOUT, timeout=900)
                except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                    report.update(passed=False, execution_completed=False,
                                  failure=dict(spacing_m=h, phase=phase,
                                  solver_exit_code=getattr(error, 'returncode', None),
                                  reason=type(error).__name__))
                    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
                    raise
            if step:
                path.append(dict(model=phase, step=step, inner_pressure_Pa=-inner))
                (folder/'load-path.json').write_text(json.dumps(path, indent=2)+'\n')
        case = dict(spacing_m=h, folder=folder.name, **compare(folder, dilation),
                    inputs_sha256={p.name:ck.digest(p) for p in sorted(folder.iterdir())
                                   if p.suffix == '.cfg' or p.name in ('annulus.poly', 'load-path.json')},
                    outputs_sha256={p.name:ck.digest(p) for p in sorted(folder.iterdir())
                                    if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.info'})
        report['cases'].append(case)
        (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        print(h, case['errors'], 'force', case['max_free_force_over_PR'], flush=True)
    assert ck.digest(executable) == binary_hash
    qualify(report)
    report['scope'] = 'Dry isotropic unloading, perfect Mohr-Coulomb, plane strain; no fluid injection.'
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert report['passed'], 'Accuracy target not met; retain raw evidence'
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--spacings', type=float, nargs='+', default=[.2, .1, .05])
    parser.add_argument('--binary-manifest', type=Path)
    parser.add_argument('--pt-option', type=int, choices=(0, 1), default=1)
    parser.add_argument('--dilation', type=float, default=0.)
    parser.add_argument('--mixed', action='store_true')
    parser.add_argument('--load-steps', type=int, default=16)
    parser.add_argument('--pt-cfl', type=float, default=.25)
    parser.add_argument('--refined', action='store_true', help='Resolve r < 3a with maximum triangle area h^2/2')
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.spacings,
          args.binary_manifest.resolve() if args.binary_manifest else None,
          args.pt_option, args.dilation, args.mixed, args.load_steps, args.pt_cfl, args.refined)
