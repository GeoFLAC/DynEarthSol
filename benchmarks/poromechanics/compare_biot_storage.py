#!/usr/bin/env python3
"""Run parent/corrected DES binaries against the same finite-grain Biot column.

Output includes input files, raw DES saves, logs, hashes, CSV histories and
PNG/SVG/PDF figures. --plot-only reuses the recorded inputs and provenance.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import prepare_figures as pf
import numpy as np


def run(output, binaries, biot, timeout, fixed_dt):
    provenance = {}
    for label, bindir in binaries.items():
        base = output / label
        base.mkdir()
        provenance[label] = {}
        for dim in (2, 3):
            case = base / f'terzaghi-{dim}d'
            case.mkdir()
            cfg = pf.make_config(dim)
            bulk = json.loads(cfg['mat']['bulk_modulus'])[0]
            cfg['mat']['biot_coeff'] = f'[{biot}]'
            if biot < 1:
                cfg['mat']['bulk_modulus_s'] = f'[{bulk / (1 - biot):.17g}]'
            if fixed_dt is not None:
                cfg['control']['fixed_dt'] = f'{fixed_dt:.17g}'
            with (case / 'input.cfg').open('w') as stream:
                cfg.write(stream)
            exe = bindir / f'dynearthsol{dim}d'
            provenance[label][str(dim)] = dict(
                executable=str(exe),
                executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                input_sha256=hashlib.sha256((case / 'input.cfg').read_bytes()).hexdigest())
            # Persist before execution so later plotting never attributes an old
            # run to a binary that was rebuilt in the meantime.
            (output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
            print(f'Running {label} {dim}D, alpha={biot}', flush=True)
            with (case / 'run.log').open('w') as stream:
                subprocess.run([str(exe), 'input.cfg'], cwd=case,
                               env=dict(os.environ, OMP_NUM_THREADS='1'),
                               stdout=stream, stderr=subprocess.STDOUT,
                               check=True, timeout=timeout)
            if hashlib.sha256(exe.read_bytes()).hexdigest() != provenance[label][str(dim)]['executable_sha256']:
                raise RuntimeError('Executable changed during the run; rebuild before benchmarking')


def plot(output, biot):
    cases = {label: {dim: pf.load_case(output / label, dim) for dim in (2, 3)}
             for label in ('before', 'after')}
    fig, axes = pf.plt.subplots(3, 2, figsize=(11, 11), constrained_layout=True)
    summary = {}
    for col, dim in enumerate((2, 3)):
        ref = cases['after'][dim]['reference']
        time = np.linspace(0, max(cases[label][dim]['times'][-1] for label in cases), 800)
        axes[0, col].plot(time, ref.pressure([ref.height], time)[:, 0] / 1000,
                          color='0.15', label='Biot analytical')
        axes[1, col].plot(time, ref.settlement(time) * 1000, color='0.15', label='Biot analytical')
        for label, color, marker in (('before', '#D55E00', 'x'), ('after', '#0072B2', 'o')):
            c = cases[label][dim]
            active = c['times'] > 0
            t = c['times'][active]
            p = c['pressure'][active][:, np.isclose(c['depth'], ref.height)].mean(axis=1)
            axes[0, col].plot(t[::3], p[::3] / 1000, marker, color=color,
                              markersize=4, fillstyle='none', label=label.capitalize())
            axes[1, col].plot(t[::3], c['settlement'][active][::3] * 1000, marker,
                              color=color, markersize=4, fillstyle='none', label=label.capitalize())
            axes[2, col].semilogy(t, c['pressure_error'][active] * 100, color=color,
                                  label=f'{label}: pressure / load')
            axes[2, col].semilogy(t, c['settlement_error'][active] * 100, '--', color=color,
                                  label=f'{label}: settlement / drained')
            summary.setdefault(label, {})[str(dim)] = c['metrics']
        axes[0, col].set(title=f'{dim}D: pressure dissipation', ylabel='Bottom excess pressure (kPa)')
        axes[1, col].set(title=f'{dim}D: settlement', ylabel='Settlement (mm)')
        axes[2, col].set(title='All post-load frames and pressure nodes',
                         ylabel='Normalized absolute error (%)')
    for ax in axes.flat:
        ax.set_xlabel('Time (s)')
        ax.grid(alpha=0.22)
        ax.legend(fontsize=8)
    fig.suptitle(f'Actual DES before/after Biot storage correction — α = {biot:g}\n'
                 'Same mesh, inputs and PT tolerances · each frame compared at its actual time', fontsize=13)
    pf.save_figure(fig, output, 'biot_storage_comparison')
    summary['biot'] = biot
    summary['provenance'] = json.loads((output / 'provenance.json').read_text())
    summary['notes'] = ['Frame zero is pre-load and excluded from analytical error metrics.',
                        'This is a coarse-column regression, not a mesh/time convergence claim.',
                        'Pressure error is normalized by applied load; settlement by drained settlement.']
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (output / 'index.html').write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<title>Biot storage correction</title><style>body{font:17px/1.5 system-ui;'
        'max-width:1100px;margin:32px auto;padding:0 20px}img{width:100%}</style>'
        '<h1>Biot storage correction: actual DES results</h1>'
        '<p>Before and after use identical inputs. Each saved pressure and settlement '
        'is compared with the finite-grain Biot analytical solution at its actual time. '
        'These coarse columns precede mesh/time convergence studies.</p>'
        '<p><a href="biot_storage_comparison.png">PNG</a> · '
        '<a href="biot_storage_comparison.svg">SVG</a> · '
        '<a href="biot_storage_comparison.pdf">PDF</a> · '
        '<a href="summary.json">Metrics and provenance</a></p>'
        '<img src="biot_storage_comparison.png" alt="Before and after pressure and settlement">'
        '</html>')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--parent-bindir', type=Path)
    parser.add_argument('--current-bindir', type=Path, default=pf.ROOT)
    parser.add_argument('--biot', type=float, default=0.6)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--fixed-dt-seconds', type=float,
                        help='Use the same explicit physical timestep in both versions')
    parser.add_argument('--plot-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    if args.plot_only:
        biot = json.loads((output / 'settings.json').read_text())['biot']
    else:
        if (args.parent_bindir is None or not 0.3 <= args.biot <= 1 or
                not (args.timeout > 0 and math.isfinite(args.timeout))):
            parser.error('Require --parent-bindir, 0.3 <= --biot <= 1 and positive --timeout')
        if args.fixed_dt_seconds is not None and not (args.fixed_dt_seconds > 0 and math.isfinite(args.fixed_dt_seconds)):
            parser.error('--fixed-dt-seconds must be positive and finite')
        output.mkdir(parents=True, exist_ok=False)
        biot = args.biot
        (output / 'settings.json').write_text(json.dumps(dict(biot=biot, fixed_dt_seconds=args.fixed_dt_seconds)) + '\n')
        run(output, dict(before=args.parent_bindir.resolve(), after=args.current_bindir.resolve()),
            biot, args.timeout, args.fixed_dt_seconds)
    plot(output, biot)


if __name__ == '__main__':
    main()
