"""Application page: pore-pressure-driven failure of a soil layer on inclined bedrock.

Interface as qualify_theis: build(study, output, record_source, steps)
re-reads the finest run from raw saves, writes figures/slope.* and
data/slope/, and returns the review record. No exact solution exists;
the comparison is with three estimates computed in slope_application.py:
pointwise yield of the poroelastic infinite slope, Davis's non-associated
slip-plane estimate and the classical infinite-slope factor of safety,
read against the base pressure at which the basal slip reaches fixed
plastic strains (it converges with mesh and step; the last PT equilibrium
does not).
"""
import json
import math
from pathlib import Path
import shutil

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import slope_application as sa

KEY = 'slope'


LEVELS = (1e-3, 2e-3)                    # mean plastic strain of the basal 0.6 m, mid-slope


def _onset(case):
    return next(f['base_pressure'] for f in case['frames'] if f['basal_plastic_fraction'] > 0)


def _curve(case):
    return (np.array([f['base_pressure'] for f in case['frames']]), np.array([f['basal_plastic_strain'] for f in case['frames']]))


def _level(case, strain):
    """Base pressure at which the mean basal plastic strain first reaches strain (it only grows)."""
    p, e = _curve(case)
    return float(np.interp(strain, e, p)) if e.max() >= strain else float('nan')


def build(study, output, record_source, steps):
    summary, fine = pf.load_study(study, record_source), pf.load_study(steps, record_source)
    cases = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    finest = cases[-1]
    folder = study/f"h{finest['spacing_m']:g}"
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.info', '.poly') or '.save.' in path.name or '.chkpt.' in path.name: record_source(path)
    again, _, _ = sa.analyse(folder)
    np.testing.assert_allclose([f['base_pressure'] for f in again], [f['base_pressure'] for f in finest['frames']], rtol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    # input.cfg restarts from the checkpoint that seed.cfg writes (gravity applied as a load); run seed.cfg first.
    shutil.copy2(folder/'slope.cfg', data/'input.cfg')
    for name in ('seed.cfg', 'slope.poly'): shutil.copy2(folder/name, data/name)
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(steps/'summary.json', data/'summary-steps.json')
    shutil.copy2(study/'runner.py', data/'slope_application.py')

    estimates = ((sa.U_POINT, 'pointwise yield', ds.MUTED, ':'), (sa.U_DAVIS, 'Davis, non-associated', ds.INK, '--'),
                 (sa.U_STAR, 'infinite slope, FS = 1', ds.HOT, '-.'))
    refined = fine['cases'][0]
    runs = [(c, summary['steps']) for c in cases] + [(refined, fine['steps'])]
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3, width=13.0)
        ax = axes[0]
        model = sa.Dynearthsol(str(folder/'slope')); rows = sa.scan_frames(str(folder/'slope'))
        x = model.read_field(rows[0]['frame'], 'coordinate'); conn = model.read_field(rows[0]['frame'], 'connectivity')
        pls = model.read_field(rows[-1]['frame'], 'plastic strain')
        tri = pf.matplotlib.tri.Triangulation(x[:, 0], x[:, 1], conn)
        art = ax.tripcolor(tri, facecolors=np.log10(np.maximum(pls, 1e-6)), cmap='magma_r', vmin=-6, vmax=-1)
        ax.set_aspect('equal'); ax.set_xlim(14, 26); ax.set_ylim(-26*math.tan(sa.BETA) - sa.D - .5, -14*math.tan(sa.BETA) + .5)
        fig.colorbar(art, ax=ax, location='bottom', shrink=.7, pad=.18,
                     label=r'$\log_{10}$ plastic strain')
        ds.panel(ax, 'a', f"Plastic strain just before sliding ({finest['spacing_m']:g} m)", xlabel='$x$ (m)', ylabel='$z$ (m)')

        ax = axes[1]
        for color, (case, n) in zip(ds.SERIES, runs):
            p, e = _curve(case)
            keep = e > 0
            ax.semilogy(p[keep]/1e3, e[keep], '--' if n != summary['steps'] else '-', color=color, lw=1.2,
                        label=f"DES, {case['spacing_m']:g} m, {n} steps")
            ax.plot(p[-1]/1e3, e[-1], 'X', color=color, ms=6)
        for u, label, color, ls in estimates:
            ax.axvline(u/1e3, color=color, ls=ls, lw=1.0, label=label)
        ax.plot([], [], 'X', color=ds.INK, ms=6, label='last PT equilibrium')
        ax.set_xlim(12, 22); ax.set_ylim(1e-6, 3e-2)
        ds.panel(ax, 'b', 'Basal slip as the pressure rises', xlabel='base pore pressure $u$ (kPa)',
                 ylabel='mean basal plastic strain')
        ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(0, -.25)).set_in_layout(False)

        ax = axes[2]
        h = np.array([c['spacing_m'] for c in cases])
        series = [(lambda c: _onset(c), 'o', ds.WATER, 'onset of basal yield')]
        series += [((lambda c, t=t: _level(c, t)), m, color, f'plastic strain {t:g}')
                   for t, m, color in zip(LEVELS, ('^', 'v'), (ds.SERIES[1], ds.SERIES[0]))]
        series += [(lambda c: c['frames'][-1]['base_pressure'], 's', ds.HOT, 'last PT equilibrium')]
        for f, m, color, label in series:
            ax.plot(h, [f(c)/1e3 for c in cases], m + '-', color=color, mfc='white', label=label)
            ax.plot([refined['spacing_m']*1.06], [f(refined)/1e3], m, color=color, ms=5)
        ax.plot([], [], 'o', color=ds.MUTED, ms=5, label=f"filled: {fine['steps']} steps")
        for u, label, color, ls in estimates:
            ax.axhline(u/1e3, color=color, ls=ls, lw=1.0)
            ax.text(.75, u/1e3 - .12, label, transform=ax.get_yaxis_transform(), ha='center', va='top', fontsize=10, color=color)
        ax.set_xscale('log'); ax.set_xticks(h, [f'{v:g}' for v in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ax.set_ylim(12, 22)
        ds.panel(ax, 'c', 'Pressure at fixed slip, by mesh', xlabel='mesh spacing (m)', ylabel='base pore pressure (kPa)')
        ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(0, -.25)).set_in_layout(False)
        pf.save_figure(fig, output/'figures', KEY)

    onset = [_onset(c) for c in cases]
    level = {t: [_level(c, t) for c, _ in runs] for t in LEVELS}
    ends = [c['frames'][-1]['base_pressure'] for c, _ in runs]
    spread = (max(level[LEVELS[0]]) - min(level[LEVELS[0]]))/np.mean(level[LEVELS[0]])
    fmt = lambda v: ', '.join(f'{u/1e3:.2f}' for u in v)
    return dict(key=KEY, theme='mechanics', status='preliminary',
        title='Application: pore pressure triggers a translational slide',
        problem=('A 2 m soil layer on 20° bedrock is loaded by its own weight and then by pore pressure rising from the '
                 'bedrock. Plastic flow starts at the base and spreads along it, and the basal slip accelerates as the '
                 'pressure rises.'),
        result=(f"Basal yield starts at {onset[-1]/1e3:.1f} kPa on the {finest['spacing_m']:g} m mesh (pointwise estimate "
                f"{sa.U_POINT/1e3:.1f} kPa). The slip then grows along one curve on every mesh and step count: the basal "
                f"plastic strain reaches {LEVELS[0]:g} at {min(level[LEVELS[0]])/1e3:.2f}–{max(level[LEVELS[0]])/1e3:.2f} kPa, "
                f"at Davis's non-associated estimate ({sa.U_DAVIS/1e3:.1f} kPa), and PT stops finding an equilibrium at "
                f"{min(ends)/1e3:.1f}–{max(ends)/1e3:.1f} kPa, below the classical factor of safety ({sa.U_STAR/1e3:.1f} kPa)."),
        scope=('2D plane strain, 40 m long with 4 m strong ends (they keep the corners elastic), Mohr–Coulomb $c=2$ kPa, '
               r'$\phi=35^\circ$, no dilation, $K=50$ MPa, $G=5.17$ MPa ($\nu=0.45$), saturated; gravity applied as a load by '
               'restart from an unstressed checkpoint; base pressure ramped over 2000 s, quasi-steady seepage. Meshes '
               + ', '.join(f"{c['spacing_m']:g}" for c in cases) + f" m at {summary['steps']} steps and "
               f"{refined['spacing_m']:g} m at {fine['steps']} steps: onset " + fmt(onset)
               + f" kPa (first order toward the estimate); pressure at basal plastic strain {LEVELS[0]:g}: " + fmt(level[LEVELS[0]])
               + f" kPa; at {LEVELS[1]:g}: " + fmt(level[LEVELS[1]]) + ' kPa; last PT equilibrium: ' + fmt(ends)
               + ' kPa. The slip curve converges; the last PT equilibrium does not (it rises with the step count and '
               'falls on the finest mesh), so it marks where the accelerating slip outruns PT, not a collapse pressure. '
               'There is no sharp limit load, probably because the strong ends buttress the layer. (a) plastic strain; (b) basal slip; '
               '(c) pressure at fixed slip.'),
        meaning=('Onset agrees with the poroelastic infinite-slope estimate. Slip develops through Davis’s non-associated '
                 'estimate and accelerates well before the classical factor of safety, as zero-dilation flow predicts: '
                 'limit-equilibrium design values overestimate the stability of such a slope. The rainfall page applies '
                 'the same slope to infiltration.'),
        metric=spread, metric_label=f'spread of u at plastic strain {LEVELS[0]:g}', dimension=2, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 u_point=sa.U_POINT, u_davis=sa.U_DAVIS, u_star=sa.U_STAR, k0=sa.K0, onset=onset,
                                 runs=[dict(spacing_m=c['spacing_m'], steps=n) for c, n in runs],
                                 level={f'{t:g}': v for t, v in level.items()}, last_equilibrium=ends), indent=1))
