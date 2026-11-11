#!/usr/bin/env python3
"""Dry plane-strain circular opening under isotropic far compression.

A quarter annulus has symmetry cuts, a free inner wall and prescribed outer
normal pressure. Compare its stress with the finite-annulus solution and
state the separate infinite-domain Kirsch truncation error. This is the
elastic first step of the tunnel ladder, not the plastic/poroelastic case.
"""

import argparse
import json
import itertools
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'benchmarks/poromechanics'))
sys.path.insert(0, str(ROOT/'tests/functional'))
from Dynearthsol import Dynearthsol
from kirsch_reference import annulus, kirsch
from pt_initial_checkpoint import fields

A, R, P = 1., 20., 1e5


from des_case import digest, read_config, reset_boundaries, run_solver


def check_binary(solver, binary_hash):
    hashes = {value['sha256'] if isinstance(value, dict) else value
              for value in solver['binaries'].values()}
    assert binary_hash in hashes


def polygon(path, spacing, refined_radius=None):
    outer_n = math.ceil(np.pi*R/(8*spacing))
    inner_n = math.ceil(np.pi*A/(2*spacing))
    outer = np.linspace(0, np.pi/2, outer_n+1)
    inner = np.linspace(np.pi/2, 0, inner_n+1)
    points = [(R*np.cos(t), -R*np.sin(t)) for t in outer]
    if refined_radius is not None:
        if not A < refined_radius < R/3:
            raise ValueError('Refined ring must lie inside the outer region seed')
        points.append((0., -refined_radius))
    points += [(A*np.cos(t), -A*np.sin(t)) for t in inner]
    if refined_radius is not None: points.append((refined_radius, 0.))
    markers = [32]*outer_n + [1]*(1+(refined_radius is not None)) + [2]*inner_n + [16]*(1+(refined_radius is not None))
    assert len(points) == len(markers)
    edges = [(i, (i+1)%len(points), marker) for i, marker in enumerate(markers)]
    if refined_radius is not None:
        # A constrained, unlabelled interior arc grades mesh size without a
        # material or boundary-condition interface. Both regions are material 0.
        angles = np.linspace(np.pi/2, 0, math.ceil(np.pi*refined_radius/(2*spacing))+1)
        ids = [outer_n+1]
        for angle in angles[1:-1]:
            ids.append(len(points)); points.append((refined_radius*np.cos(angle), -refined_radius*np.sin(angle)))
        ids.append(len(markers)-1)
        edges += [(i, j, 0) for i, j in zip(ids[:-1], ids[1:])]
    with path.open('w') as stream:
        stream.write(f'{len(points)} 2 0 0\n')
        for i, (x, z) in enumerate(points):
            stream.write(f'{i} {x:.17g} {z:.17g}\n')
        stream.write(f'{len(edges)} 1\n')
        for i, (start, end, marker) in enumerate(edges):
            stream.write(f'{i} {start} {end} {marker}\n')
        if refined_radius is None:
            stream.write(f'0\n1\n0 {R/3} {-R/3} 0 {16*spacing**2}\n')
        else:
            inner_seed = (A+refined_radius)/(2*np.sqrt(2))
            stream.write(f'0\n2\n0 {R/3} {-R/3} 0 {16*spacing**2}\n')
            stream.write(f'1 {inner_seed} {-inner_seed} 0 {spacing**2/2}\n')


def config(spacing):
    cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
    cfg['sim'].update(modelname='kirsch', max_steps='1', output_step_interval='1',
                      has_initial_checkpoint='true', checkpoint_frame_interval='1')
    cfg['mesh'] = dict(meshing_option='90', poly_filename='annulus.poly',
                       xlength=str(R), ylength='1', zlength=str(R), resolution=str(spacing),
                       remeshing_option='0', min_angle='28', max_steiner_factor='0')
    cfg['control'].update(gravity='0', has_hydraulic_diffusion='false', has_thermal_diffusion='false',
                          has_PT='true', PT_option='1', PT_max_iter='200000',
                          PT_relative_tolerance='1e-9', PT_absolute_tolerance='1e-8',
                          fixed_dt='1', has_moving_mesh='false', is_using_mixed_stress='false')
    reset_boundaries(cfg['bc'], vbc='0', vbc_val='0', stress_bc='0', stress_val='0')
    cfg['bc'].update(vbc_x0='1', vbc_z0='1', stress_bc_z1='4', stress_val_z1=str(-P))
    cfg['ic']['has_body_force_adjustment'] = 'false'
    cfg['mat'].update(rheology_type='elastic', bulk_modulus='[1e9]', shear_modulus='[1e9]',
                      alpha='[0]', is_plane_strain='true')
    return cfg


def compare(folder):
    model = Dynearthsol(str(folder/'kirsch'))
    assert model.ndims == 2
    x = model.read_field(0, 'coordinate')
    conn = model.read_field(0, 'connectivity')
    centre = x[conn].mean(axis=1)
    radius = np.linalg.norm(centre, axis=1)
    radial = centre/radius[:, None]
    values = model.read_field(model.frames[-1], 'stress')
    nx, nz = radial.T
    rr = nx**2*values[:, 0] + nz**2*values[:, 1] + 2*nx*nz*values[:, 2]
    hoop = nz**2*values[:, 0] + nx**2*values[:, 1] - 2*nx*nz*values[:, 2]
    shear = nx*nz*(values[:, 1]-values[:, 0]) + (nx**2-nz**2)*values[:, 2]
    band = (radius >= 1.1*A) & (radius <= 5*A)
    exact_rr, exact_hoop = annulus(radius[band], A, R, P)
    errors = dict(radial_rms_over_load=float(np.sqrt(np.mean((rr[band]-exact_rr)**2)))/P,
                  hoop_rms_over_load=float(np.sqrt(np.mean((hoop[band]-exact_hoop)**2)))/P,
                  shear_rms_over_load=float(np.sqrt(np.mean(shear[band]**2)))/P,
                  max_error_over_load=float(max(np.max(np.abs(rr[band]-exact_rr)),
                                                np.max(np.abs(hoop[band]-exact_hoop))))/P)
    assert all(np.isfinite(v) for v in errors.values())
    strain = model.read_field(model.frames[-1], 'strain')
    er = nx**2*strain[:, 0]+nz**2*strain[:, 1]+2*nx*nz*strain[:, 2]
    et = nz**2*strain[:, 0]+nx**2*strain[:, 1]-2*nx*nz*strain[:, 2]
    bulk, modulus = 1e9, 1e9
    lame = bulk-2*modulus/3
    exact_volume = (exact_rr+exact_hoop)/(2*(lame+modulus))
    exact_er = (exact_rr-lame*exact_volume)/(2*modulus)
    exact_et = (exact_hoop-lame*exact_volume)/(2*modulus)
    errors['radial_strain_rms_over_P_G'] = float(np.sqrt(np.mean((er[band]-exact_er)**2)))/(P/modulus)
    errors['hoop_strain_rms_over_P_G'] = float(np.sqrt(np.mean((et[band]-exact_et)**2)))/(P/modulus)
    out_of_plane = np.frombuffer(fields(folder/f'kirsch.chkpt.{model.frames[-1]:06d}')['stressyy'], dtype=np.float64)
    poisson = (3*bulk-2*modulus)/(2*(3*bulk+modulus))
    errors['out_of_plane_rms_over_load'] = float(np.sqrt(np.mean(
        (out_of_plane[band]-poisson*(exact_rr+exact_hoop))**2)))/P
    force_error = force_balance(x, conn, values, P, R)
    return dict(nodes=len(x), elements=len(conn), errors=errors,
                max_free_force_over_PR=force_error,
                maximum_abs_strain=float(np.max(np.abs(model.read_field(model.frames[-1], 'strain')))))


def force_balance(x, conn, values, pressure, outer_radius):
    """Independently assemble the quarter-annulus P1 force residual."""
    vertices = np.concatenate((np.ones((*conn.shape, 1)), x[conn]), axis=2)
    area = np.abs(np.linalg.det(vertices))/2
    gradients = np.linalg.inv(vertices)[:, 1:, :].swapaxes(1, 2)
    tensor = np.zeros((len(conn), 2, 2))
    tensor[:, 0, 0], tensor[:, 1, 1] = values[:, 0], values[:, 1]
    tensor[:, 0, 1] = tensor[:, 1, 0] = values[:, 2]
    force = np.zeros_like(x)
    np.add.at(force, conn, -area[:, None, None]*(gradients@tensor.swapaxes(1, 2)))
    edges = {}
    for cell in conn:
        for edge in itertools.combinations(cell, 2):
            pair = tuple(sorted(int(n) for n in edge))
            edges[pair] = edges.get(pair, 0)+1
    external = np.zeros(2)
    for edge, count in edges.items():
        pair = list(edge)
        if count != 1 or not np.allclose(np.linalg.norm(x[pair], axis=1), outer_radius, rtol=0, atol=1e-10):
            continue
        tangent = x[pair[1]]-x[pair[0]]
        normal_measure = np.array([tangent[1], -tangent[0]])
        if normal_measure@x[pair].mean(axis=0) < 0: normal_measure *= -1
        traction = -pressure*normal_measure
        force[pair] += traction/2
        external += traction
    np.testing.assert_allclose(external, [-pressure*outer_radius, pressure*outer_radius], rtol=2e-14, atol=1e-8)
    force[np.isclose(x[:, 0], 0), 0] = 0
    force[np.isclose(x[:, 1], 0), 1] = 0
    return float(np.max(np.linalg.norm(force, axis=1)))/(pressure*outer_radius)


def study(executable, output, spacings, manifest):
    output.mkdir(parents=True, exist_ok=False)
    binary_hash = digest(executable)
    solver = json.loads(manifest.read_text()) if manifest else None
    if solver:
        check_binary(solver, binary_hash)
        shutil.copy2(manifest, output/'solver-qualification.json')
    for source in (Path(__file__), ROOT/'benchmarks/poromechanics/kirsch_reference.py'):
        shutil.copy2(source, output/source.name)
    summary = dict(executable=str(executable), executable_sha256=binary_hash,
                   solver_commit=solver['implementation_commit'] if solver else None,
                   runner_sha256=digest(Path(__file__)),
                   reference=dict(opening_radius_m=A, outer_radius_m=R, outer_pressure_Pa=P,
                                  finite_boundary_hoop_factor=R**2/(R**2-A**2)), cases=[])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for h in spacings:
        folder = output/f'h{h:g}'; folder.mkdir()
        polygon(folder/'annulus.poly', h)
        with (folder/'input.cfg').open('w') as stream: config(h).write(stream)
        run_solver([str(executable), 'input.cfg'], log=folder / 'run.log', cwd=folder, env=env, check=True, timeout=600)
        case = dict(spacing_m=h, folder=folder.name, **compare(folder), input_sha256=digest(folder/'input.cfg'),
                    mesh_input_sha256=digest(folder/'annulus.poly'),
                    outputs_sha256={p.name:digest(p) for p in sorted(folder.glob('kirsch.*'))
                                    if '.save.' in p.name or '.chkpt.' in p.name or p.suffix == '.info'})
        summary['cases'].append(case)
        (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
        print(h, case['errors'], flush=True)
    assert digest(executable) == binary_hash
    fine = min(summary['cases'], key=lambda c:c['spacing_m'])
    largest = max(summary['cases'], key=lambda c:c['spacing_m'])
    passed = (fine['max_free_force_over_PR'] < 1e-7 and
              fine['errors']['radial_rms_over_load'] < .01 and
              fine['errors']['hoop_rms_over_load'] < .01 and
              fine['errors']['max_error_over_load'] < .03 and
              fine['errors']['radial_strain_rms_over_P_G'] < .01 and
              fine['errors']['hoop_strain_rms_over_P_G'] < .01 and
              fine['errors']['out_of_plane_rms_over_load'] < .005)
    if len(summary['cases']) > 1:
        passed &= fine['errors']['hoop_rms_over_load'] < largest['errors']['hoop_rms_over_load']
    report = dict(summary, passed=bool(passed), analysis_sha256=digest(Path(__file__)),
                  reference_sha256=digest(ROOT/'benchmarks/poromechanics/kirsch_reference.py'),
                  scope='Dry isotropic compression, plane-strain annulus. Plasticity, anisotropic far loading and coupled injection are not qualified.')
    (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    assert passed, 'Selected-mesh accuracy target not met; retain the raw evidence'
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path); parser.add_argument('output', type=Path)
    parser.add_argument('--spacings', type=float, nargs='+', default=[.2, .1, .05, .025])
    parser.add_argument('--binary-manifest', type=Path)
    args = parser.parse_args()
    study(args.executable.resolve(), args.output.resolve(), args.spacings,
          args.binary_manifest.resolve() if args.binary_manifest else None)
