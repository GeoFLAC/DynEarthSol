"""Validate the recorded DES Bishop-cell studies (transient coupling) and plot them.

Interface as qualify_theis. build(study, output, record_source, space, before,
restart) takes the 2D step sequence (with the explicit runs), the 3D step
sequence, the 2D sequence of the build before the fixed-stress correction and
the restart checks.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_bishop_cell as cc

KEY = 'bishop_cell'
NAMES = dict(moderate=r'moderate suction, $\chi\Pi^\prime>0$', strong=r'strong suction, $\Pi^\prime<0$')


def _check(study, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    return summary


def _cases(summary, state, mode='PT'):
    return sorted((c for c in summary['cases'] if c['state'] == state and c['mode'] == mode), key=lambda c: c['steps'])


def build(study, output, record_source, space, before, restart):
    summary, other, old = (_check(s, record_source) for s in (study, space, before))
    checks = {state: json.loads((restart/state/'restart.json').read_text()) for state in cc.STATES}
    for state in cc.STATES:
        record_source(restart/state/'restart.json')
        assert checks[state]['same_step_different_payloads'] == 0
    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'check_bishop_cell.py')
    for name, s in (('summary-3d.json', space), ('summary-before.json', before)):
        shutil.copy2(s/'summary.json', data/name)
    # The page selects the finest moderate-suction PT case from this study.
    selected = _cases(summary, 'moderate')[-1]
    folder = study/f"moderate-pt-{selected['steps']}"
    shutil.copy2(folder/'cell.cfg', data/'input.cfg')
    record_source(folder/'cell.cfg')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        for ax, letter, state in zip(axes[:2], 'ab', cc.STATES):
            p0 = cc.STATES[state][0]
            eps = np.linspace(0, cc.STRAIN, 200)
            ds.analytical(ax, -1e3*eps, cc.exact(state, eps) - p0, color=ds.INK, label='exact rate equation')
            for color, cases, label in ((ds.HOT, _cases(summary, state), 'corrected'), (ds.MUTED, _cases(old, state), 'before correction')):
                c = cases[-1]                                    # the finest step of each build
                ax.plot([-1e3*f['strain'] for f in c['frames']], [f['pressure_mean'] - p0 for f in c['frames']], 'o',
                        ms=4, mfc='white', mec=color, mew=.9, label=f"PT, {c['steps']} steps, {label}")
            x = _cases(summary, state, 'explicit')[0]
            ax.plot([-1e3*f['strain'] for f in x['frames']], [f['pressure_mean'] - p0 for f in x['frames']], 's',
                    ms=3, color=ds.WATER, label='explicit strain form')
            ds.panel(ax, letter, NAMES[state], xlabel=r'compression $-\varepsilon_{zz}$ ($10^{-3}$)', ylabel=r'$p-p_0$ (Pa)')
            ax.legend(fontsize=8)
        ax = axes[2]
        for color, state in ((ds.HOT, 'moderate'), (ds.WATER, 'strong')):
            for s, style, label in ((summary, dict(marker='o', ls='-'), '2D'), (other, dict(marker='s', ls=':'), '3D'),
                                    (old, dict(marker='x', ls='--'), '2D before')):
                cs = _cases(s, state)
                ax.loglog([c['steps'] for c in cs], [c['max_error_over_change'] for c in cs], color=color, mfc='white',
                          ms=4, label=f"{state}, {label}", **style)
        n = np.array([10, 320])
        ax.loglog(n, .5*(10/n), '--', color=ds.MUTED, lw=.8, label='slope $-1$')
        ds.panel(ax, 'c', 'Step convergence', xlabel='PT steps over the loading', ylabel='max error / pressure change')
        ax.legend(fontsize=7, ncol=2, loc='lower left')
        pf.save_figure(fig, output/'figures', KEY)

    fine = {state: _cases(summary, state)[-1] for state in cc.STATES}
    fine3 = {state: _cases(other, state)[-1] for state in cc.STATES}
    stale = {state: _cases(old, state)[-1] for state in cc.STATES}
    pct = lambda x: f'{100*x:.2f}%'
    return dict(key=KEY, theme='unsaturated', title='Bishop coupling in a sealed cell under compression',
        problem=('A sealed, partly saturated soil cell is compressed slowly. Without drainage, the pore pressure follows '
                 'the strain through the Bishop rate equation; the stabilised fixed-stress coupling used with the PT solver must recover it as the '
                 'step shrinks, also where the Bishop pressure Se p decreases with suction.'),
        result=(f"Largest pressure error at {fine['moderate']['steps']} steps: {pct(fine['moderate']['max_error_over_change'])} "
                f"(moderate suction) and {pct(fine['strong']['max_error_over_change'])} (strong suction) in 2D; the error halves with "
                f"each halving of the step. Without the Bishop correction of the fixed-stress source the error stayed at {pct(stale['moderate']['max_error_over_change'])} and "
                f"{pct(stale['strong']['max_error_over_change'])}."),
        scope=(f"Uniaxial compression to a vertical strain of $-{cc.STRAIN*1e3:g}\\times10^{{-3}}$ over {cc.DURATION:g} s, K = G = 10 kPa, Gardner law; 2D and 3D, "
               f"10 to {fine['moderate']['steps']} PT steps; 3D at {fine3['moderate']['steps']} steps: "
               f"{pct(fine3['moderate']['max_error_over_change'])} and {pct(fine3['strong']['max_error_over_change'])}. "
               f"The explicit strain form is exact to {max(_cases(summary, s, 'explicit')[0]['max_error_over_change'] for s in cc.STATES):.0e}. "
               f"A restart halfway reproduces the uninterrupted run exactly; a restart with a quarter of the step stays on the exact path. "
               '(a, b) pressure along the strain path; (c) step convergence.'),
        meaning=('This checks the transient Bishop coupling of the PT scheme: the fixed-stress source corrected by '
                 '$L(\\chi E-D)$, with the increments the mechanics consumed. It does not test flow, which the column '
                 'checks cover, or large strain.'),
        metric=fine['strong']['max_error_over_change'], metric_label='pressure error / change', dimension=2, spacing_m=0.5,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 cases_2d=[{k: v for k, v in c.items() if k != 'frames'} for c in summary['cases']],
                                 cases_3d=[{k: v for k, v in c.items() if k != 'frames'} for c in other['cases']],
                                 before=[{k: v for k, v in c.items() if k != 'frames'} for c in old['cases']],
                                 restart=checks), indent=1))
