#!/usr/bin/env python3
"""Run existing DES consolidation cases and write a reproducible figure gallery.

Theis plots are analytical references only. Terzaghi plots compare actual DES
output with an independent small-strain Biot reference. No solver source is
modified. Use --plot-only to redraw a previously generated run directory.
"""

import argparse
import configparser
import hashlib
import html
import json
import math
import os
from pathlib import Path
import subprocess
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/des-hydromechanics-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from terzaghi_reference import TerzaghiReference
from theis_reference import drawdown, scheduled_drawdown

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol, YEAR2SEC, scan_frames

COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7']


def save_figure(fig, output, name):
    # Keep dense fields compact in SVG/PDF while preserving vector labels.
    for ax in fig.axes:
        for artist in ax.collections:
            values = artist.get_array()
            if values is not None and values.size > 1000:
                artist.set_rasterized(True)
    outside = [ax.get_legend() for ax in fig.axes
               if ax.get_legend() is not None and not ax.get_legend().get_in_layout()]
    for suffix in ('png', 'svg', 'pdf'):
        fig.savefig(output / f'{name}.{suffix}', dpi=180, bbox_inches='tight',
                    bbox_extra_artists=outside or None)
    plt.close(fig)


def make_config(dim):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(ROOT / f'benchmarks/poromechanics/1d-consolidation-des{dim}d_traction.cfg')
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
        with (case / 'run.log').open('w') as stream:
            subprocess.run(command, cwd=case, env=environment, check=True,
                           stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
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


def theis_figures(output):
    trans, storage, rate = 1e-3, 1e-3, 0.01
    radius = np.geomspace(1, 500, 300)
    time = np.geomspace(1, 14400, 400)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for t, color in zip((60, 600, 3600, 14400), COLORS):
        axes[0, 0].semilogx(radius, -drawdown(radius, t, rate, trans, storage),
                           color=color, label=f'{t / 60:g} min')
    for r, color in zip((10, 25, 50, 100), COLORS):
        axes[0, 1].semilogx(time / 60, -drawdown(r, time, rate, trans, storage),
                           color=color, label=f'r = {r} m')
    schedule = [(0.0, rate), (3600.0, 0.0), (7200.0, -rate), (10800.0, 0.0)]
    axes[1, 0].step([t / 60 for t, _ in schedule] + [240],
                    [q * 3600 for _, q in schedule] + [schedule[-1][1] * 3600],
                    where='post', color=COLORS[0], linewidth=2)
    axes[1, 0].axhline(0, color='0.5', linewidth=0.7)
    axes[1, 0].set(xlabel='Time (min)', ylabel='Extraction rate (m³/h)', title='Pumping → shut-in → injection → shut-in')
    scheduled_time = np.unique(np.r_[np.linspace(0, 14400, 1500), [x[0] for x in schedule]])
    histories = []
    for r, color in zip((10, 25, 50, 100), COLORS):
        head = -scheduled_drawdown(r, scheduled_time, schedule, trans, storage)
        histories.append(head)
        axes[1, 1].plot(scheduled_time / 60, head, color=color, label=f'r = {r} m')
    axes[1, 1].axhline(0, color='0.5', linewidth=0.7)
    axes[0, 0].set(xlabel='Radius from well (m)', ylabel='Head change (m)', title='Constant pumping: radial profiles')
    axes[0, 1].set(xlabel='Time (min)', ylabel='Head change (m)', title='Constant pumping: observation wells')
    axes[1, 1].set(xlabel='Time (min)', ylabel='Head change (m)', title='Response to the scheduled rates')
    for ax in axes.flat:
        ax.grid(alpha=0.22)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 1]):
        ax.legend(fontsize=8)
    fig.suptitle('Theis analytical reference — DES well simulation pending\n'
                 'T = 0.001 m²/s · S = 0.001 · |Q| = 0.01 m³/s · constant properties', fontsize=13)
    save_figure(fig, output, 'theis_reference')
    np.savetxt(output / 'theis-scheduled-head.csv', np.column_stack([scheduled_time] + histories),
               delimiter=',', header='time_s,head_change_r10_m,head_change_r25_m,head_change_r50_m,head_change_r100_m', comments='')
    return dict(transmissivity_m2_s=trans, storativity=storage,
                schedule_seconds_rate_m3_s=schedule, simulation_status='analytical reference only')


def consolidation_figures(output, cases):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ref = cases[2]['reference']
    smooth = np.linspace(0.1, max(c['times'][-1] for c in cases.values()), 600)
    axes[0, 0].plot(smooth, ref.pressure([ref.height], smooth)[:, 0] / 1000, color='0.15', label='Biot analytical')
    axes[0, 1].plot(smooth, ref.settlement(smooth) * 1000, color='0.15', label='Biot analytical')
    for dim, color, marker in ((2, COLORS[0], 'o'), (3, COLORS[1], '^')):
        case = cases[dim]
        times, pressure = case['times'], case['pressure']
        bottom = np.isclose(case['depth'], ref.height)
        active = times > 0
        axes[0, 0].plot(times[active][::4], pressure[active][:, bottom].mean(axis=1)[::4] / 1000,
                       marker, color=color, markersize=4, fillstyle='none', label=f'DES {dim}D')
        axes[0, 1].plot(times[active][::4], case['settlement'][active][::4] * 1000,
                       marker, color=color, markersize=4, fillstyle='none', label=f'DES {dim}D')
        axes[1, 1].semilogy(times[active], case['pressure_error'][active] * 100,
                          color=color, label=f'{dim}D pressure / load')
        axes[1, 1].semilogy(times[active], case['settlement_error'][active] * 100,
                          '--', color=color, label=f'{dim}D settlement / drained')
    case = cases[2]
    for target, color in zip((100, 1000, 5000), COLORS):
        i = np.argmin(np.abs(case['times'] - target))
        depths, inverse = np.unique(case['depth'].round(10), return_inverse=True)
        averaged = np.bincount(inverse, weights=case['pressure'][i]) / np.bincount(inverse)
        axes[1, 0].plot(averaged / 1000, depths, 'o', color=color, markersize=4, fillstyle='none')
        exact_depth = np.linspace(0, ref.height, 200)
        axes[1, 0].plot(ref.pressure(exact_depth, case['times'][i])[0] / 1000, exact_depth,
                       color=color, label=f'{case["times"][i]:.0f} s')
    axes[0, 0].set(xlabel='Time (s)', ylabel='Bottom excess pressure (kPa)', title='Pressure dissipation')
    axes[0, 1].set(xlabel='Time (s)', ylabel='Top settlement (mm)', title='Settlement integrated from axial strain')
    axes[1, 0].set(xlabel='Excess pressure (kPa)', ylabel='Depth below drained top (m)', title='2D profiles: symbols DES, lines analytical')
    axes[1, 0].invert_yaxis()
    axes[1, 1].set(xlabel='Time (s)', ylabel='Normalized absolute error (%)', title='All saved post-load frames; all pressure nodes')
    for ax in axes.flat:
        ax.grid(alpha=0.22)
        ax.legend(fontsize=8)
    fig.suptitle('Actual DES baseline vs classical Terzaghi consolidation\n'
                 'PT option 1 · absolute tolerance 10⁻⁶ · fixed mesh · frame 0 is pre-load', fontsize=13)
    save_figure(fig, output, 'terzaghi_comparison')

    fig, axes = plt.subplots(1, 3, figsize=(9, 6), constrained_layout=True)
    for ax, target in zip(axes, (100, 1000, 5000)):
        i = np.argmin(np.abs(case['times'] - target))
        image = ax.tripcolor(case['coordinate'][:, 0], -case['coordinate'][:, -1],
                             case['connectivity'], case['pressure'][i] / 1000,
                             shading='gouraud', cmap='viridis', vmin=0, vmax=ref.load / 1000)
        ax.set(xlabel='Horizontal position (m)', ylabel='Depth (m)',
               title=f't = {case["times"][i]:.0f} s')
        ax.invert_yaxis()
    fig.colorbar(image, ax=axes, label='Excess pore pressure (kPa)', shrink=0.8)
    fig.suptitle('Actual DES 2D pressure field\nTop drained · bottom and sides sealed · 100 kPa traction load', fontsize=13)
    save_figure(fig, output, 'terzaghi_pressure_fields')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot-only', action='store_true')
    parser.add_argument('--timeout', type=float, default=180)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.plot_only:
        if not output.is_dir():
            parser.error('--plot-only requires an existing run directory')
    else:
        output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    runs = {} if args.plot_only else run_cases(output, args.timeout)
    if args.plot_only and (output / 'summary.json').exists():
        runs = json.loads((output / 'summary.json').read_text()).get('runs', {})
    theis = theis_figures(output)
    cases = {dim: load_case(output, dim) for dim in (2, 3)}
    consolidation_figures(output, cases)
    summary = dict(source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   runs=runs, theis=theis,
                   terzaghi={str(dim): case['metrics'] for dim, case in cases.items()},
                   notes=['Theis figures are analytical references, not DES well simulations.',
                          'Consolidation uses existing DES inputs with PT_option=1, PT_absolute_tolerance=1e-6, gravity=10.',
                          'Frame zero is pre-load and excluded from post-load analytical error metrics.',
                          'These baseline results are not a mesh/time convergence qualification.'])
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    sections = []
    for name, title in [('terzaghi_comparison', 'DES vs analytical consolidation'),
                        ('terzaghi_pressure_fields', 'DES pressure fields'),
                        ('theis_reference', 'Theis analytical reference')]:
        sections.append(f'<section><h2>{title}</h2><a href="{name}.png"><img src="{name}.png" alt="{title}"></a>'
                        f'<p><a href="{name}.png">PNG</a> · <a href="{name}.svg">SVG</a> · <a href="{name}.pdf">PDF</a></p></section>')
    body = ('<!doctype html><html lang="en"><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>Hydromechanics benchmark figures</title>'
            '<style>body{font:17px/1.5 system-ui;max-width:1150px;margin:40px auto;padding:0 24px;background:#f6f7f8;color:#17212b}'
            'section{background:white;padding:24px;margin:24px 0;border-radius:12px}img{width:100%;height:auto}'
            'pre{white-space:pre-wrap;font-size:13px}a{color:#006e9b}</style>'
            '<h1>Hydromechanics benchmark figures</h1>'
            '<p>Terzaghi figures show actual 2D/3D DES simulations. Theis figures show analytical references only; DES well controls are pending.</p>'
            '<p>Baseline comparisons precede mesh/time convergence studies. Frame zero is pre-load and excluded from error metrics.</p>'
            '<p><a href="summary.json">Summary and provenance</a> · '
            '<a href="terzaghi-2d/input.cfg">2D input</a> · <a href="terzaghi-2d/run.log">2D log</a> · '
            '<a href="terzaghi-3d/input.cfg">3D input</a> · <a href="terzaghi-3d/run.log">3D log</a> · '
            '<a href="terzaghi-2d-history.csv">2D data</a> · <a href="terzaghi-3d-history.csv">3D data</a> · '
            '<a href="theis-scheduled-head.csv">Theis data</a></p>' + ''.join(sections) +
            '<section><h2>Measured errors</h2><pre>' + html.escape(json.dumps(summary['terzaghi'], indent=2)) +
            '</pre></section></html>')
    (output / 'index.html').write_text(body)
    print(json.dumps(summary['terzaghi'], indent=2))
    print(f'Figures: {output / "index.html"}')


if __name__ == '__main__':
    main()
