#!/usr/bin/env python3
"""Elder's heater problem: selected finest DES mesh and reference grid.

    plot_elder.py STUDY REFERENCE.npz OUTPUT_PREFIX [--times 0.1] [--all-grids]

STUDY is the output folder of tests/functional/hydromechanics/check_elder.py.
The first figure shows the temperature theta = (T - T0)/dT at the selected
time, with the 0.2 and 0.6
isotherms of Voss and Souza (1987). The second shows the heat content of the
layer against tau. --all-grids exposes the refinement controls. Nothing is fitted or smoothed.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(Path(__file__).resolve().parent))
from Dynearthsol import Dynearthsol, scan_frames
from diagram_style import RC, FIGURE_RC, INK, MUTED, panel

T0, DT, KAPPA = 273.0, 10.0, 1e-6
LEVELS = np.linspace(0, 1, 11)


def des_cases(study, all_grids=False):
    summary = json.loads((study/'summary.json').read_text())
    cases = summary['cases']
    if not all_grids:
        cases = [min(cases, key=lambda c: c['spacing_m'])]
    out = []
    for case in cases:
        folder = study/f"h{case['spacing_m']:g}"
        model = Dynearthsol(str(folder/'elder'))
        rows = {round(r['time']*KAPPA, 9): r for r in scan_frames(str(folder/'elder'))}
        coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
        tri = mtri.Triangulation(coord[:, 0], coord[:, 1] + 1, conn)
        out.append(dict(label=f"DES, h = 1/{round(1/case['spacing_m'])} m", spacing_m=case['spacing_m'], model=model, rows=rows, tri=tri,
                        series=case['series'], dt=summary['dt_s']))
    return out


def theta_des(case, tau):
    row = case['rows'][round(tau, 9)]
    return (case['model'].read_field(row['frame'], 'temperature') - T0)/DT


def fields(study, reference, prefix, times, all_grids=False, include_heat=False):
    ref = np.load(reference)
    grids = sorted(int(k.split('_')[1]) for k in ref.files if k.startswith('theta_'))
    des = des_cases(study, all_grids=all_grids)
    if not all_grids:
        grids = [grids[-1]]
    if include_heat:
        assert len(times) == 1 and not all_grids
    columns = [(f'Reference, {n} x {2*n}', ('ref', n)) for n in grids] + [(c['label'], ('des', c)) for c in des]
    with plt.rc_context(FIGURE_RC if include_heat else RC):
        size = (12, 3.9) if include_heat else (2.6*len(columns), 1.45*len(times)+.6)
        fig, axes = plt.subplots(len(times), len(columns)+int(include_heat), figsize=size,
                                 sharex=not include_heat, sharey=not include_heat, squeeze=False,
                                 constrained_layout=include_heat)
        for i, tau in enumerate(times):
            k = int(np.argmin(np.abs(ref['times'] - tau)))
            assert abs(ref['times'][k] - tau) < 1e-12, tau
            for j, (title, (kind, item)) in enumerate(columns):
                ax = axes[i, j]
                if kind == 'ref':
                    X, Z = np.meshgrid(ref[f'x_{item}'], ref[f'z_{item}'])
                    th = ref[f'theta_{item}'][k]
                    im = ax.contourf(X, Z, th, LEVELS, cmap='inferno', extend='both')
                    ax.contour(X, Z, th, [.2, .6], colors='white', linewidths=.6)
                else:
                    th = theta_des(item, tau)
                    im = ax.tricontourf(item['tri'], th, LEVELS, cmap='inferno', extend='both')
                    ax.tricontour(item['tri'], th, [.2, .6], colors='white', linewidths=.6)
                ax.plot([0, 1], [0, 0], color='#b23a24', lw=3, solid_capstyle='butt', clip_on=False)
                ax.set_aspect('equal'); ax.set_xlim(0, 2); ax.set_ylim(0, 1)
                ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1])
                if i == 0: ax.set_title(f'({chr(97+j)}) {title}' if include_heat else title, fontsize=10, color=INK)
                if j == 0: ax.set_ylabel(f'$\\tau$ = {tau:g}\n$z/H$', fontsize=10)
                if i == len(times)-1: ax.set_xlabel('$x/H$', fontsize=10)
        fig.colorbar(im, ax=list(axes[:, :len(columns)].flat), shrink=.6,
                     label=r'$\theta = (T - T_0)/\Delta T$', ticks=[0, .2, .4, .6, .8, 1])
        if include_heat:
            heat(study, reference, prefix, ax=axes[0, -1])
            panel(axes[0, -1], 'c', 'Heat content through plume formation')
        target = str(prefix) if include_heat else f'{prefix}-fields'
        for ext in ('png', 'svg', 'pdf'):
            fig.savefig(f'{target}.{ext}', dpi=200, bbox_inches='tight')
        plt.close(fig)


def heat(study, reference, prefix, all_grids=False, ax=None):
    summary = json.loads((Path(reference).with_suffix('.json')).read_text())
    des = json.loads((study/'summary.json').read_text())
    cases = summary['cases']
    if not all_grids:
        n = max(cases, key=int)
        cases = {n: cases[n]}
        des['cases'] = [min(des['cases'], key=lambda c: c['spacing_m'])]
    own_figure = ax is None
    with plt.rc_context(RC if own_figure else FIGURE_RC):
        if own_figure: fig, ax = plt.subplots(figsize=(5.2, 3.4))
        else: fig = ax.figure
        shades = plt.cm.Greys(np.linspace(.45, .9, len(summary['cases'])))
        for (n, case), color in zip(cases.items(), shades):
            pts = case['at_times']
            ax.plot([p['tau'] for p in pts], [p['heat'] for p in pts], '-', color=color, lw=1.2,
                    label=f'reference, {n} x {2*int(n)}')
        marks = ['o', 's', '^', 'D']
        for case, m in zip(des['cases'], marks):
            s = case['series']
            ax.plot([p['tau'] for p in s], [p['heat_des'] for p in s], m, mfc='none', color='#24618f', ms=6,
                    label=f"DES, h = 1/{round(1/case['spacing_m'])} m")
        ax.set_xlabel(r'$\tau = \kappa t/(\sigma H^2)$'); ax.set_ylabel(r'heat content $\int\theta\,dA\,/H^2$')
        ax.set_xlim(0, .105); ax.legend(frameon=False, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        if own_figure:
            for ext in ('png', 'svg', 'pdf'):
                fig.savefig(f'{prefix}-heat.{ext}', dpi=200, bbox_inches='tight')
            plt.close(fig)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('study', type=Path); p.add_argument('reference', type=Path); p.add_argument('prefix')
    p.add_argument('--times', type=float, nargs='+', default=[0.1])
    p.add_argument('--all-grids', action='store_true', help='show every mesh/reference control')
    a = p.parse_args()
    fields(a.study, a.reference, a.prefix, a.times, all_grids=a.all_grids)
    heat(a.study, a.reference, a.prefix, all_grids=a.all_grids)
