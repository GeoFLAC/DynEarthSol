"""Validate the recorded DES rate-and-state velocity-step study and plot it.

Interface as qualify_theis: build(study, output, record_source) replays the
finest step from raw saves, writes figures/rsf.* and data/rsf/, and returns
the review record.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_rsf_step as rs

KEY = 'rsf'


def build(study, output, record_source):
    summary = pf.load_study(study, record_source)
    cases = sorted(summary['cases'], key=lambda c: -c['dt'])
    for case in cases:
        folder = study/f"r{case['refine']:g}"
        des = np.concatenate([rs.read(folder, 'cell'), rs.read(folder, 'step', rs.T1)])
        for column, name in ((0, 'time'), (2, 'mu'), (3, 'stress')):
            np.testing.assert_array_equal(des[:, column], case['series'][name])
        expected = rs.discrete(case['dt'], int(rs.T1/case['dt']), int(rs.T2/case['dt']))
        np.testing.assert_allclose(des[:, 1], expected[:, 1], rtol=1e-14, atol=0)
        for path in folder.iterdir():
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name or '.chkpt.' in path.name: record_source(path)
    finest = cases[-1]
    folder = study/f"r{finest['refine']:g}"
    des = np.concatenate([rs.read(folder, 'cell'), rs.read(folder, 'step', rs.T1)])
    np.testing.assert_allclose(des[:, 2], finest['series']['mu'], rtol=0, atol=1e-15)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'cell.cfg', data/'input.cfg')
    shutil.copy2(folder/'step.cfg', data/'step.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'check_rsf_step.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        tt = np.linspace(rs.T1 - 16, rs.T1 + rs.T2, 2000)
        selected = finest
        for ax, key, scale, label in ((axes[0], 'mu', 1.0, r'friction coefficient $\mu$'),
                                      (axes[1], 'stress', 1e-3, 'principal shear $t$ (kPa)')):
            mu_c = rs.continuum(tt)
            phi = np.arctan(mu_c)
            exact = mu_c if key == 'mu' else (rs.COHESION*np.cos(phi) + (rs.C - rs.BIOT*rs.P)*np.sin(phi))*scale
            ds.analytical(ax, tt - rs.T1, exact, label='continuum')
            ref = rs.discrete(selected['dt'], int(rs.T1/selected['dt']), int(rs.T2/selected['dt']))
            col = 2 if key == 'mu' else 3
            ax.plot(ref[:, 0] - rs.T1, ref[:, col]*scale, '-', color=ds.MUTED, lw=.8, drawstyle='steps-post',
                    label=f"forward Euler, $\\Delta t={selected['dt']:g}$ s")
            t = np.array(selected['series']['time']); v = np.array(selected['series'][key])*scale
            indices = np.flatnonzero(t >= rs.T1-16)
            keep = np.unique(np.r_[indices[::max(1, len(indices)//80)], np.flatnonzero((t >= rs.T1) & (t <= rs.T1+selected['dt']))])
            ds.des(ax, t[keep] - rs.T1, v[keep], label=f"DES, $\\Delta t={selected['dt']:g}$ s", size=3)
            ax.axvline(0, color=ds.MUTED, lw=.6, ls=':')
            ax.set_xlim(-16, rs.T2)
            if key == 'stress':
                ax.set_ylim(np.min(exact) - .6, np.max(exact) + .9)
            ax.legend(fontsize=10, loc='upper right')
        ds.panel(axes[0], 'a', 'Friction across a tenfold rate step', xlabel='$t-t_1$ (s)', ylabel=r'$\mu$')
        ds.panel(axes[1], 'b', 'Shear stress at yield', xlabel='$t-t_1$ (s)', ylabel='principal shear $t$ (kPa)')
        ax = axes[2]
        dt = np.array([c['dt'] for c in cases])
        ax.loglog(dt, [c['mu_vs_continuum'] for c in cases], 'o-', color=ds.HOT, mfc='white', label='DES − continuum')
        ax.loglog(dt, cases[0]['mu_vs_continuum']*dt/dt[0], '--', color=ds.MUTED, lw=.8, label='slope 1')
        ax.set_xticks(dt, [f'{v:g}' for v in dt]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ticks = [1e-4, 2e-4, 5e-4, 1e-3]
        ax.set_yticks(ticks, [f'{v:.0e}'.replace('e-0', 'e-') for v in ticks]); ax.yaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Friction error by time step', xlabel=r'$\Delta t$ (s)', ylabel=r'max $|\mu-\mu_{\rm exact}|$ after the step')
        ax.legend(fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)

    worst = max(max(c['theta_vs_discrete'], c['mu_vs_discrete']) for c in cases)
    return dict(key=KEY, theme='mechanics', title='Rate-and-state friction under a velocity step',
        problem=('A pressurized cell sheared at yield sees its loading rate rise tenfold. Rate-and-state friction jumps by '
                 '$a\\ln10$, then decays over the slip distance $D_c$ to the new steady state ($b>a$: velocity weakening).'),
        result=(f"DES follows its exact discrete recursion to round-off at every step (state and friction identical, stress "
                f"within {max(c['stress_vs_discrete'] for c in cases):.0e}), so the state is updated once per physical time step, "
                f"not once per PT iteration. Against the continuum solution the largest friction error after the step is "
                + ', '.join(f"{c['mu_vs_continuum']:.2e}" for c in cases) + ' for steps of '
                + ', '.join(f"{c['dt']:g}" for c in cases) + f" s: first order, against a jump of {cases[0]['mu_jump']:.4f}."),
        scope=('2D unit cell with every node prescribed, so equilibrium and mesh errors cannot enter; EP-RSF, aging law, '
               f"$a={summary['a']:g}$, $b={summary['b']:g}$, $D_c={summary['dc']:g}$ m, slip rate $V=2w\\dot\\epsilon_{{II}}$ "
               "(total-strain-rate projection); confinement 1 MPa, pore pressure 0.8 MPa, cohesion 100 kPa. "
               'The rate step is a restart from the last checkpoint with the new boundary velocity, which also checks that the '
               'state variable survives restart. Not covered: slip-law evolution, the maximum-shear projection, '
               'stick-slip instability and coupling of friction to flow or dilation. (a) friction; (b) stress; (c) step convergence.'),
        meaning=('This verifies the evolving rate-and-state friction of the effective-stress plasticity, which the strength '
                 'page checks only at constant friction.'),
        metric=cases[-1]['mu_vs_continuum'], metric_label='friction error', dimension=2, spacing_m=1.0,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 discrete_worst=worst,
                                 cases=[{k: v for k, v in c.items() if k != 'series'} for c in cases]), indent=1))
