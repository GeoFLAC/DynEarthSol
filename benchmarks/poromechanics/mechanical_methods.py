#!/usr/bin/env python3
"""Explore explicit/GVS/PT Terzaghi accuracy at actual physical output times.

This is a method comparison, not an acceptance test. Inaccurate and failed
cases stay in the report. Use --analyze-only to replay retained, hashed saves.
Run one dimension at a time with a 6 GiB virtual-memory limit.
"""

import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

import consolidation_refinement as cr
import consolidation as pf


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import (check_hashes, check_manifest, digest, file_hashes, material_properties, read_config, run_des,
                      source_diff, source_revision, write_config, write_json)
from des_test_utils import facet_areas, simplex_gradients, simplex_volumes

def cases(settings):
    methods = settings.get('methods', ['pt0', 'pt1', 'explicit', 'gvs-shear', 'gvs-bulk'])
    for option in (0, 1):
        if f'pt{option}' in methods:
            yield f'pt{option}', dict(has_PT='true', PT_option=str(option),
                                      inertial_scaling='1e5', fixed_dt=str(settings['pt_dt_s']))
    for speed, scaling, method in itertools.product(settings['initial_speeds_m_s'],
                                                     settings.get('scalings', [1e4, 1e5, 1e6]),
                                                     [m for m in methods if not m.startswith('pt')]):
        control = dict(has_PT='false', fixed_dt='0', characteristic_speed=str(speed),
                       inertial_scaling=str(scaling), dt_fraction=str(settings['dt_fraction']),
                       use_global_velocity_scaling=str(method != 'explicit').lower())
        if method != 'explicit':
            control['mass_scaling_reference_speed'] = method.split('-')[1]
        yield f'{method}-S{scaling:g}-v{speed:g}', control


def run(output, settings):
    executable = Path(settings.get('executable') or
                      pf.ROOT / f'dynearthsol{settings["dimension"]}d').resolve()
    (output / 'source.diff').write_bytes(source_diff())
    manifest = settings.get('binary_manifest')
    solver = check_manifest(manifest, digest(executable))
    if solver:
        (output/'solver-qualification.json').write_bytes(Path(manifest).read_bytes())
    provenance = dict(source_commit=source_revision(pf.ROOT),
                      solver_commit=solver.get('implementation_commit') if solver else None,
                      binary_manifest_sha256=digest(manifest) if manifest else None,
                      source_diff_sha256=digest(output / 'source.diff'),
                      executable_sha256=digest(executable),
                      analysis_sha256={name: digest(Path(__file__).with_name(name)) for name in
                                       ('mechanical_methods.py', 'consolidation.py',
                                        'consolidation_refinement.py', 'terzaghi_reference.py')},
                      reader_sha256=digest(pf.ROOT / 'Dynearthsol.py'), cases={})
    for name, control in cases(settings):
        folder = output / name / f'terzaghi-{settings["dimension"]}d'
        folder.mkdir(parents=True)
        cfg = cr.configuration(settings['dimension'], settings['spacing_m'], settings['pt_dt_s'],
                               dict(duration_s=settings['duration_s'], sample_interval_s=settings['duration_s']/10))
        # K/G=2 distinguishes the two GVS density floors. All methods use
        # the same physical column, load, material and hydraulic operator.
        cfg['mat']['bulk_modulus'] = '[2e7]'
        cfg['control'].update(control)
        cfg['sim'].pop('output_step_interval', None)
        cfg['sim'].update(max_steps=str(settings['max_steps']),
                          output_time_interval_in_yr=str(settings['duration_s']/10/pf.YEAR2SEC),
                          checkpoint_frame_interval='1', info_display_step_interval='100000')
        write_config(folder / 'input.cfg', cfg)
        record = dict(control=control, input_sha256=digest(folder / 'input.cfg'), status='running')
        provenance['cases'][name] = record
        write_json(output / 'provenance.json', provenance)
        start = time.monotonic()
        try:
            run_des(executable, folder, 'input', None, settings['timeout_s'])
            if digest(executable) != provenance['executable_sha256']:
                raise RuntimeError('Executable changed during comparison')
            record['status'] = 'complete'
        except (subprocess.SubprocessError, RuntimeError) as error:
            record.update(status='failed', error=str(error))
        record.update(elapsed_s=time.monotonic()-start,
                      output_sha256=file_hashes(p for p in folder.glob('terzaghi.*')
                                                if '.save.' in p.name or '.chkpt.' in p.name))
        write_json(output / 'provenance.json', provenance)
        print(name, record['status'], flush=True)


def diagnostics(folder, case, cfg):
    """Independent force assembly and final-state mechanical energy diagnostics."""
    model = pf.Dynearthsol(str(folder / 'terzaghi'))
    frame = model.frames[-1]
    x, conn = case['coordinate'], case['connectivity']
    dim = x.shape[1]
    volume = simplex_volumes(x, conn)
    gradients = simplex_gradients(x, conn)
    values = model.read_field(frame, 'stress')
    stress = np.zeros((len(conn), dim, dim))
    for d in range(dim):
        stress[:, d, d] = values[:, d]
    for component, (a, b) in enumerate(itertools.combinations(range(dim), 2), start=dim):
        stress[:, a, b] = stress[:, b, a] = values[:, component]
    force = np.zeros_like(x)
    np.add.at(force, conn, -volume[:, None, None] * (gradients @ stress.swapaxes(1, 2)))
    facets = {}
    for nodes in conn:
        for face in itertools.combinations(nodes, dim):
            key = tuple(sorted(int(n) for n in face))
            facets[key] = facets.get(key, 0) + 1
    load = case['reference'].load
    for face, count in facets.items():
        face = list(face)
        if count != 1 or not np.all(np.isclose(x[face, -1], 0)):
            continue
        force[face, -1] -= load * facet_areas(x[face]) / dim
    for d in range(dim - 1):
        force[np.isclose(x[:, d], 0) | np.isclose(x[:, d], 1), d] = 0
    force[np.isclose(x[:, -1], -10), -1] = 0
    strain = model.read_field(frame, 'strain')
    trace = strain[:, :dim].sum(axis=1)
    m = material_properties(cfg)
    shear, bulk = m['shear'], m['bulk']
    energy = np.sum(volume * ((bulk-2*shear/3)*trace**2/2 + shear*(
        np.sum(strain[:, :dim]**2, axis=1) + 2*np.sum(strain[:, dim:]**2, axis=1))))
    # Physical reference mass uses the constant initial saturated mixture.
    # This is a diagnostic, not a discrete coupled energy-balance assertion.
    density = .7*2650 + .3*1000
    velocity = model.read_field(frame, 'velocity')
    kinetic = np.sum(volume * density * np.sum(velocity[conn]**2, axis=(1, 2))/(2*(dim+1)))
    return dict(final_saved_stress_force_max_over_total_load=float(np.max(np.linalg.norm(force, axis=1)) /
                                                                  (load*np.sum(volume)/10)),
                final_physical_reference_kinetic_J=float(kinetic),
                final_skeleton_strain_energy_J=float(energy),
                final_physical_reference_kinetic_over_strain_energy=float(kinetic/energy) if energy > 0 else None)


def analyze(output, settings):
    provenance = json.loads((output / 'provenance.json').read_text())
    report = dict(settings=settings, analysis_sha256=digest(Path(__file__)),
                  qualification='Exploratory; stability alone does not establish accuracy.', cases={})
    for name, record in provenance['cases'].items():
        folder = output / name / f'terzaghi-{settings["dimension"]}d'
        if digest(folder / 'input.cfg') != record['input_sha256']:
            raise RuntimeError(f'{name}: changed configuration')
        check_hashes(record['output_sha256'], folder)
        result = dict(status=record['status'], elapsed_s=record['elapsed_s'])
        if record['status'] == 'complete':
            try:
                case = pf.load_case(output / name, settings['dimension'])
                metrics = dict(case['metrics'], **diagnostics(folder, case, read_config(folder / 'input.cfg')))
                if any(value is not None and not np.isfinite(value) for value in metrics.values()):
                    raise RuntimeError('Nonfinite analysis metric')
                result.update(metrics)
                late = case['times'] >= settings['duration_s']/2
                result['late_max_pressure_error_over_load'] = float(case['pressure_error'][late].max())
                result['late_max_settlement_error_over_drained'] = float(case['settlement_error'][late].max())
                result['drainage_time_s'] = case['reference'].height**2/case['reference'].diffusivity
            except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                result.update(status='unqualified', error=str(error))
        else:
            result['error'] = record['error']
        report['cases'][name] = result
    write_json(output / 'summary.json', report)
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--dimension', type=int, choices=(2, 3), default=2)
    parser.add_argument('--executable', type=Path,
                        help='Use a frozen executable instead of the shared-checkout binary')
    parser.add_argument('--binary-manifest', type=Path,
                        help='Bind the frozen executable to its qualified solver commit')
    parser.add_argument('--spacing', type=float, default=1)
    parser.add_argument('--duration', type=float, default=400)
    parser.add_argument('--dt-fraction', type=float, default=.1)
    parser.add_argument('--pt-dt', type=float, default=.125)
    parser.add_argument('--initial-speeds', type=float, nargs='+', default=[1e-12, 1e-6])
    parser.add_argument('--scalings', type=float, nargs='+', default=[1e4, 1e5, 1e6])
    parser.add_argument('--methods', nargs='+', choices=('pt0', 'pt1', 'explicit', 'gvs-shear', 'gvs-bulk'),
                        default=['pt0', 'pt1', 'explicit', 'gvs-shear', 'gvs-bulk'])
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--max-steps', type=int, default=2000000)
    parser.add_argument('--analyze-only', dest='analyze_only', action='store_true')
    args = parser.parse_args()
    if (not all(math.isfinite(v) and v > 0 for v in [args.spacing, args.duration,
            args.dt_fraction, args.pt_dt, args.timeout, *args.initial_speeds, *args.scalings]) or
            args.dt_fraction > 1 or args.max_steps < 1):
        parser.error('Use finite positive settings, dt-fraction <= 1 and max-steps >= 1')
    output = args.output.resolve()
    if args.analyze_only:
        settings = json.loads((output / 'settings.json').read_text())
    else:
        settings = dict(dimension=args.dimension, spacing_m=args.spacing, duration_s=args.duration,
                        executable=str(args.executable.resolve()) if args.executable else None,
                        binary_manifest=str(args.binary_manifest.resolve()) if args.binary_manifest else None,
                        dt_fraction=args.dt_fraction, pt_dt_s=args.pt_dt,
                        initial_speeds_m_s=args.initial_speeds, timeout_s=args.timeout,
                        scalings=args.scalings,
                        methods=args.methods,
                        max_steps=args.max_steps)
        output.mkdir(parents=True, exist_ok=False)
        write_json(output / 'settings.json', settings)
        run(output, settings)
    analyze(output, settings)


if __name__ == '__main__':
    main()
