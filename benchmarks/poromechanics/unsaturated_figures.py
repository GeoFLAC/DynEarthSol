#!/usr/bin/env python3
"""Export actual DES U6/U7 profiles for the implementation review.

The input directory contains completed unsaturated-2d/ and unsaturated-3d/
outputs from check_unsaturated_column.py, including their summary.json.
"""

import argparse
import configparser
import json
from pathlib import Path
import sys

import numpy as np

import consolidation_refinement as cr
import prepare_figures as pf
from unsaturated_reference import Retention, srivastava_yeh, steady_pressure


def analyze(source, output):
    provenance = dict(source=str(source), analysis_sha256=cr.digest(Path(__file__)),
                      reference_sha256=cr.digest(Path(__file__).with_name('unsaturated_reference.py')),
                      reader_sha256=cr.digest(pf.ROOT / 'Dynearthsol.py'),
                      helper_sha256={name: cr.digest(Path(__file__).with_name(name)) for name in
                                     ('prepare_figures.py', 'consolidation_refinement.py')},
                      files={}, qualification='Transport stage with chi=1; Bishop mechanics is unfinished.')

    def record(path):
        provenance['files'][str(path.relative_to(source))] = cr.digest(path)

    def load(folder):
        cfg = configparser.ConfigParser()
        cfg.read(folder / 'column.cfg')
        record(folder / 'column.cfg')
        if (folder / 'column.info').exists():
            record(folder / 'column.info')
        model = pf.Dynearthsol(str(folder / 'column'))
        for frame in model.frames:
            record(Path(model.get_fn(frame)))
        times = {r['frame']: r['time'] for r in pf.scan_frames(str(folder / 'column'))}
        return cfg, model, times

    def value(cfg, name):
        return json.loads(cfg['mat'][name])[0]

    histories, steady_rows = [], []
    summaries = {dim: json.loads((source / f'unsaturated-{dim}d/summary.json').read_text()) for dim in (2, 3)}
    for dim, summary in summaries.items():
        if not all(key in summary for key in ('transient', 'modes', 'steady', 'restart')):
            raise RuntimeError(f'{dim}D verification is incomplete')
        record(source / f'unsaturated-{dim}d/summary.json')
    output.mkdir(parents=True, exist_ok=False)
    (output / 'data').mkdir()
    fig, axes = pf.plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for row, dim in enumerate((2, 3)):
        for col, result in enumerate(summaries[dim]['transient']):
            spacing = result['resolution']
            cfg, model, times = load(source / f'unsaturated-{dim}d/transient-{spacing}')
            height = cfg.getfloat('mesh', 'zlength')
            gamma = value(cfg, 'fluid_rho0') * cfg.getfloat('control', 'gravity')
            alpha, ks = value(cfg, 'retention_alpha')*gamma, value(cfg, 'hydraulic_perm')*gamma/value(cfg, 'fluid_visc')
            dtheta = value(cfg, 'porosity')*(1-value(cfg, 'residual_saturation'))
            flux = -json.loads(cfg['bc']['hbc_flux'])[-1]
            coord = model.read_field(0, 'coordinate')
            z = coord[:, -1] + height
            grid = np.linspace(0, height, 401)
            measured = []
            for i, frame in enumerate(model.frames):
                pressure = model.read_field(frame, 'pore pressure')
                kr = np.exp(value(cfg, 'retention_alpha')*pressure)
                exact = srivastava_yeh(z, times[frame], length=height, alpha=alpha,
                                      ks=ks, dtheta=dtheta, qa=0, qb=flux)
                error = float(np.max(np.abs(kr-exact)))
                measured.append(error)
                histories.append((dim, spacing, times[frame], error))
                np.savetxt(output / f'data/transient-{dim}d-h{spacing}-frame{frame}.csv',
                           np.column_stack((z, pressure, kr, exact)), delimiter=',',
                           header='elevation_above_base_m,pressure_Pa,DES_K_over_Ks_from_pressure,reference_K_over_Ks', comments='')
                if frame == 0:
                    continue
                color = pf.COLORS[(i-1) % len(pf.COLORS)]
                axes[row, col].plot(grid, srivastava_yeh(grid, times[frame], length=height,
                    alpha=alpha, ks=ks, dtheta=dtheta, qa=0, qb=flux), color=color, label=f'reference {times[frame]:g} s')
                axes[row, col].scatter(z, kr, s=8, facecolors='none', edgecolors=color)
            np.testing.assert_allclose(measured, result['max_error'], rtol=0, atol=1e-9)
            axes[row, col].set(title=f'{dim}D, h={spacing:g} m (markers: DES)',
                                xlabel='elevation above base (m)', ylabel='K / Ks (from pressure)')
            axes[row, col].legend(fontsize=7)
            axes[row, col].grid(alpha=.2)
    pf.save_figure(fig, output, 'transient-profiles')
    np.savetxt(output / 'data/transient-errors.csv', histories, delimiter=',',
               header='dimension,spacing_m,time_s,max_K_over_Ks_error', comments='')

    fig, axes = pf.plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for col, dim in enumerate((2, 3)):
        for row, law in enumerate(('van_genuchten', 'brooks_corey')):
            for color, result in zip(pf.COLORS, summaries[dim]['steady'][law]):
                spacing = result['resolution']
                cfg, model, times = load(source / f'unsaturated-{dim}d/steady-{law}-{spacing}')
                height = cfg.getfloat('mesh', 'zlength')
                z = model.read_field(0, 'coordinate')[:, -1] + height
                pressure = model.read_field(model.frames[-1], 'pore pressure')
                retention = Retention(law, value(cfg, 'retention_alpha'), value(cfg, 'retention_n'),
                                      value(cfg, 'residual_saturation'))
                options = dict(flux=-json.loads(cfg['bc']['hbc_flux'])[-1],
                               perm=value(cfg, 'hydraulic_perm'), mu=value(cfg, 'fluid_visc'),
                               rho=value(cfg, 'fluid_rho0'), g=cfg.getfloat('control', 'gravity'))
                exact = steady_pressure(retention, z, **options)
                error = float(np.max(np.abs(pressure-exact)))
                np.testing.assert_allclose(error, result['max_error_pa'], rtol=0, atol=1e-7)
                axes[row, col].scatter(z, pressure, s=10, facecolors='none', edgecolors=color,
                                        label=f'DES h={spacing:g}, t={times[model.frames[-1]]:g} s')
                grid = np.linspace(0, height, 401)
                axes[row, col].plot(grid, steady_pressure(retention, grid, **options), color=color,
                                     linewidth=1, label=f'reference h={spacing:g}')
                np.savetxt(output / f'data/steady-{law}-{dim}d-h{spacing}.csv',
                           np.column_stack((z, pressure, exact)), delimiter=',',
                           header='elevation_above_base_m,DES_pressure_Pa,reference_pressure_Pa', comments='')
                steady_rows.append((dim, law, spacing, error))
            axes[row, col].set(title=f'{dim}D {law.replace("_", " ")}',
                                xlabel='elevation above base (m)', ylabel='pressure (Pa)')
            axes[row, col].legend(fontsize=7)
            axes[row, col].grid(alpha=.2)
    pf.save_figure(fig, output, 'steady-profiles')
    cr.write_json(output / 'metrics.json', dict(transient=histories, steady=steady_rows))
    cr.write_json(output / 'provenance.json', provenance)
    (output / 'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Unsaturated verification</title>'
        '<h1>DES unsaturated transport verification</h1><p>Transport stage, chi=1; Bishop mechanics is unfinished.</p>'
        '<p>Lines: independent references. Markers: DES pressure, or conductivity derived from that pressure. '
        'Times come from raw save metadata.</p>'
        '<p><a href="metrics.json">Metrics</a> · <a href="provenance.json">Input/output hashes</a></p>'
        '<img src="transient-profiles.svg" width="1000" alt="Transient conductivity profiles">'
        '<img src="steady-profiles.svg" width="1000" alt="Steady pressure profiles">\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    analyze(args.source.resolve(), args.output.resolve())
