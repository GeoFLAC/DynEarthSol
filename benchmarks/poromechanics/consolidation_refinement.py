#!/usr/bin/env python3
"""Separate mesh and physical-timestep refinement of DES Terzaghi columns.

Run dimensions sequentially with OMP/BLAS threads set to one and a 6 GiB
virtual-memory ceiling. Existing destinations are refused; --plot-only
reconstructs the report from recorded configurations and raw saves.
"""

import argparse
import configparser
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

import prepare_figures as pf


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def configuration(dim, spacing, dt, settings):
    cfg = pf.make_config(dim)
    cfg['mesh'] = dict(meshing_option='1', meshing_elem_shape='1', xlength='1',
                       ylength='1', zlength='10', resolution=repr(spacing),
                       remeshing_option='0')
    cfg['control'].update(gravity='0', has_thermal_diffusion='false',
                          fixed_dt=repr(dt), PT_absolute_tolerance='1e-8')
    # Step-aligned output gives every refinement the same physical sample times.
    cfg['sim'].pop('output_time_interval_in_yr', None)
    cfg['sim'].update(max_steps=str(round(settings['duration_s'] / dt)),
                      max_time_in_yr=repr(settings['duration_s'] / pf.YEAR2SEC),
                      output_step_interval=str(round(settings['sample_interval_s'] / dt)))
    return cfg


def row_bound(case, cfg, dim):
    """Independent constant-property FE row bound, including held nodes.

    The study requires one hydraulic substep, so that mesh refinement cannot
    also change hydraulic time resolution through automatic subcycling.
    """
    x, conn = case['coordinate'], case['connectivity']
    cells = x[conn]
    delta = cells[:, 1:] - cells[:, :1]
    volume = np.abs(np.linalg.det(delta)) / math.factorial(dim)
    grad = np.linalg.inv(delta).swapaxes(1, 2)
    grad = np.concatenate((-grad.sum(axis=1, keepdims=True), grad), axis=1)
    mat = lambda name: json.loads(cfg['mat'][name])[0]
    bulk = mat('bulk_modulus') + (mat('shear_modulus') / 3 if dim == 2 else 0)
    storage = (mat('porosity') / mat('fluid_bulk_modulus') +
               (mat('biot_coeff') - mat('porosity')) / mat('bulk_modulus_s') +
               mat('biot_coeff')**2 / bulk)
    mobility = mat('hydraulic_perm') / mat('fluid_visc')
    operator = mobility * volume[:, None, None] * (grad @ grad.swapaxes(1, 2))
    mass, rows = np.zeros(len(x)), np.zeros(len(x))
    np.add.at(mass, conn, np.broadcast_to((storage * volume / (dim + 1))[:, None], conn.shape))
    np.add.at(rows, conn, np.abs(operator).sum(axis=2))
    return float(np.min(mass[rows > 0] / rows[rows > 0]))


def plan(settings):
    runs = {}
    for dim in settings['dimensions']:
        for role, values in [('spatial', settings['spacings_m']),
                             ('temporal', settings['timesteps_s'])]:
            for value in values:
                h = value if role == 'spatial' else settings['temporal_spacing_m']
                dt = settings['spatial_dt_s'] if role == 'spatial' else value
                key = f'{dim}d-h{h:.17g}-dt{dt:.17g}'
                run = runs.setdefault(key, dict(dimension=dim, spacing_m=h, dt_s=dt, roles=[]))
                run['roles'].append(role)
    return runs


def simulate(output, settings):
    provenance = dict(source_commit=subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=pf.ROOT, text=True).strip(),
        analysis_sha256={name: digest(Path(__file__).with_name(name)) for name in
                         ('consolidation_refinement.py', 'prepare_figures.py', 'terzaghi_reference.py')},
        reader_sha256=digest(pf.ROOT / 'Dynearthsol.py'), runs=plan(settings))
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for key, run in provenance['runs'].items():
        dim = run['dimension']
        group = output / key
        folder = group / f'terzaghi-{dim}d'
        folder.mkdir(parents=True)
        cfg = configuration(dim, run['spacing_m'], run['dt_s'], settings)
        with (folder / 'input.cfg').open('w') as stream:
            cfg.write(stream)
        exe = pf.ROOT / f'dynearthsol{dim}d'
        run.update(executable=str(exe), executable_sha256=digest(exe),
                   input_sha256=digest(folder / 'input.cfg'), status='running')
        write_json(output / 'provenance.json', provenance)
        print(f'Running {key} ({", ".join(run["roles"])})', flush=True)
        start = time.monotonic()
        try:
            with (folder / 'run.log').open('w') as stream:
                subprocess.run([str(exe), 'input.cfg'], cwd=folder, env=env,
                               stdout=stream, stderr=subprocess.STDOUT,
                               timeout=settings['timeout_s'], check=True)
            if digest(exe) != run['executable_sha256']:
                raise RuntimeError('Executable changed during the run')
        except Exception as error:
            run.update(status='failed', error=str(error))
            write_json(output / 'provenance.json', provenance)
            raise
        run.update(status='complete', elapsed_s=time.monotonic() - start,
                   saves_sha256={p.name: digest(p) for p in sorted(folder.glob('terzaghi.save.*'))})
        write_json(output / 'provenance.json', provenance)


def orders(errors, steps):
    return [math.log(a / b) / math.log(h / k) if a > 0 and b > 0 else None
            for a, b, h, k in zip(errors, errors[1:], steps, steps[1:])]


def analyze(output, settings):
    provenance = json.loads((output / 'provenance.json').read_text())
    samples = settings['sample_interval_s'] * np.arange(
        round(settings['duration_s'] / settings['sample_interval_s']) + 1)
    cases, records = {}, []
    for key, run in provenance['runs'].items():
        dim = run['dimension']
        folder = output / key / f'terzaghi-{dim}d'
        if run['status'] != 'complete' or digest(folder / 'input.cfg') != run['input_sha256']:
            raise RuntimeError(f'{key}: incomplete run or changed configuration')
        for name, sha in run['saves_sha256'].items():
            if digest(folder / name) != sha:
                raise RuntimeError(f'{key}: changed raw output {name}')
        case = pf.load_case(output / key, dim)
        np.testing.assert_allclose(case['times'], samples, rtol=0, atol=1e-8)
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(folder / 'input.cfg')
        bound = row_bound(case, cfg, dim)
        if run['dt_s'] > bound * (1 + 1e-12):
            raise RuntimeError(f'{key}: dt={run["dt_s"]} exceeds the one-substep bound {bound}')
        cases[key] = case
        records.append(dict(key=key, **run, hydraulic_row_bound_s=bound, **case['metrics']))

    summary = dict(settings=settings, provenance=provenance, analysis_sha256=digest(Path(__file__)),
                   runs=records, dimensions={})
    fig, axes = pf.plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for dim in settings['dimensions']:
        spatial = sorted((r for r in records if r['dimension'] == dim and 'spatial' in r['roles']),
                         key=lambda r: -r['spacing_m'])
        temporal = sorted((r for r in records if r['dimension'] == dim and 'temporal' in r['roles']),
                          key=lambda r: -r['dt_s'])
        h = [r['spacing_m'] for r in spatial]
        pe = [r['max_pressure_error_over_load'] for r in spatial]
        se = [r['max_settlement_error_over_drained'] for r in spatial]
        if not all(b['elements'] > a['elements'] for a, b in zip(spatial, spatial[1:])):
            raise RuntimeError(f'{dim}D: nominal refinement did not increase element counts')
        changes = []
        for coarse, fine in zip(temporal, temporal[1:]):
            a, b = cases[coarse['key']], cases[fine['key']]
            np.testing.assert_array_equal(a['coordinate'], b['coordinate'])
            np.testing.assert_array_equal(a['connectivity'], b['connectivity'])
            changes.append(dict(dt_s=coarse['dt_s'], fine_dt_s=fine['dt_s'],
                                pressure=float(np.max(np.abs(a['pressure'][1:] - b['pressure'][1:])) /
                                               a['reference'].load),
                                settlement=float(np.max(np.abs(a['settlement'][1:] - b['settlement'][1:])) /
                                                 a['reference'].drained_settlement)))
        dt = [r['dt_s'] for r in changes]
        pc, sc = [r['pressure'] for r in changes], [r['settlement'] for r in changes]
        checks = dict(spatial_pressure_decreases=all(b < a for a, b in zip(pe, pe[1:])),
                      spatial_settlement_decreases=all(b < a for a, b in zip(se, se[1:])),
                      temporal_pressure_change_decreases=all(b < a for a, b in zip(pc, pc[1:])),
                      temporal_settlement_change_decreases=all(b < a for a, b in zip(sc, sc[1:])),
                      finest_pressure_below_one_percent=pe[-1] < 0.01,
                      finest_settlement_below_half_percent=se[-1] < 0.005)
        summary['dimensions'][str(dim)] = dict(checks=checks, spatial_pressure_orders=orders(pe, h),
            spatial_settlement_orders=orders(se, h), temporal_changes=changes,
            temporal_pressure_orders=orders(pc, dt), temporal_settlement_orders=orders(sc, dt))
        for ax, x, y in [(axes[0, 0], h, pe), (axes[0, 1], h, se),
                          (axes[1, 0], dt, pc), (axes[1, 1], dt, sc)]:
            ax.loglog(x, y, 'o-', label=f'{dim}D')
    for ax in axes.flat:
        ax.grid(alpha=0.25, which='both')
        ax.legend()
    axes[0, 0].set(title='Mesh: pressure error / applied load', xlabel='Spacing (m)')
    axes[0, 1].set(title='Mesh: settlement error / drained value', xlabel='Spacing (m)')
    axes[1, 0].set(title='Time: pressure difference, dt versus dt/2', xlabel='Physical dt (s)')
    axes[1, 1].set(title='Time: settlement difference, dt versus dt/2', xlabel='Physical dt (s)')
    fig.suptitle('DES Terzaghi refinement at common physical sample times')
    pf.save_figure(fig, output, 'terzaghi_refinement')
    summary['passed'] = all(all(d['checks'].values()) for d in summary['dimensions'].values())
    write_json(output / 'summary.json', summary)
    columns = ['dimension', 'spacing_m', 'dt_s', 'nodes', 'elements', 'hydraulic_row_bound_s',
               'max_pressure_error_over_load', 'max_settlement_error_over_drained', 'elapsed_s']
    with (output / 'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    rows = ''.join(f'<tr><td>{html.escape(r["key"])}</td><td>{r["elements"]}</td>'
                   f'<td>{100*r["max_pressure_error_over_load"]:.4g}%</td>'
                   f'<td>{100*r["max_settlement_error_over_drained"]:.4g}%</td></tr>' for r in records)
    (output / 'index.html').write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><title>Terzaghi refinement</title>'
        '<style>body{font:17px/1.5 system-ui;max-width:1100px;margin:32px auto;padding:0 20px}'
        'img{width:100%}td,th{padding:6px 16px;text-align:left}</style>'
        '<h1>Terzaghi mesh and timestep refinement</h1><p>Actual DES 2D/3D output; '
        'fixed mesh, zero gravity, independent Biot reference. Spatial runs share one dt; '
        'temporal runs share one mesh. Each run remains below its independently assembled '
        'one-substep hydraulic bound. Frame zero is pre-load and excluded.</p>'
        '<p>Pressure error is normalized by load, settlement by drained settlement. '
        'Time plots show differences between successive timesteps on the same mesh. '
        'Observed orders apply to these samples; no universal order is asserted.</p>'
        f'<p>Acceptance checks: {"PASS" if summary["passed"] else "FAIL"}. '
        '<a href="summary.json">Full metrics and provenance</a> · <a href="metrics.csv">CSV</a> · '
        '<a href="terzaghi_refinement.svg">SVG</a> · <a href="terzaghi_refinement.pdf">PDF</a></p>'
        '<img src="terzaghi_refinement.png" alt="Mesh and timestep error trends">'
        '<table><tr><th>Run</th><th>Elements</th><th>Pressure/load</th>'
        f'<th>Settlement/drained</th></tr>{rows}</table></html>')
    print(json.dumps(summary['dimensions'], indent=2))
    return summary['passed']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dimensions', nargs='+', type=int, choices=(2, 3), default=[2, 3])
    parser.add_argument('--spacings', nargs='+', type=float, default=[1.0, 0.5, 0.25])
    parser.add_argument('--timesteps', nargs='+', type=float, default=[1.0, 0.5, 0.25, 0.125])
    parser.add_argument('--spatial-dt', type=float, default=0.125)
    parser.add_argument('--temporal-spacing', type=float, default=0.5)
    parser.add_argument('--duration', type=float, default=400.0)
    parser.add_argument('--sample-interval', type=float, default=100.0)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    if args.plot_only:
        settings = json.loads((output / 'settings.json').read_text())
    else:
        values = args.spacings + args.timesteps + [args.spatial_dt, args.temporal_spacing,
                                                  args.duration, args.sample_interval, args.timeout]
        if not all(math.isfinite(x) and x > 0 for x in values):
            parser.error('All spacings, timesteps, intervals and timeout must be finite and positive')
        if len(set(args.spacings)) < 3 or len(set(args.timesteps)) < 3:
            parser.error('Require at least three distinct spacings and timesteps')
        if max(args.spacings + [args.temporal_spacing]) > 1:
            parser.error('Spacing must not exceed the fixed one-metre column width')
        timesteps = sorted(set(args.timesteps), reverse=True)
        if not all(math.isclose(a / b, 2.0) for a, b in zip(timesteps, timesteps[1:])):
            parser.error('Temporal refinement requires successive timestep halvings')
        for a, b in [(args.duration, args.sample_interval)] + [
                (args.sample_interval, dt) for dt in args.timesteps + [args.spatial_dt]]:
            if a < b or not math.isclose(a / b, round(a / b), abs_tol=1e-10, rel_tol=0):
                parser.error('Sample interval must divide duration; every dt must divide sample interval')
        settings = dict(dimensions=sorted(set(args.dimensions)), spacings_m=sorted(set(args.spacings), reverse=True),
                        timesteps_s=sorted(set(args.timesteps), reverse=True), spatial_dt_s=args.spatial_dt,
                        temporal_spacing_m=args.temporal_spacing, duration_s=args.duration,
                        sample_interval_s=args.sample_interval, timeout_s=args.timeout)
        output.mkdir(parents=True, exist_ok=False)
        write_json(output / 'settings.json', settings)
        simulate(output, settings)
    if not analyze(output, settings):
        raise SystemExit('Terzaghi refinement acceptance checks failed; see summary.json')


if __name__ == '__main__':
    main()
