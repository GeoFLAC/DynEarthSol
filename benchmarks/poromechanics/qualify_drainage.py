"""Validate the recorded DES coupled-drainage study (unsaturated level U9) and plot it.

Interface as qualify_theis: build(study, output, record_source, steps, sanity)
replays the finest PT run against an independently recomputed reference,
writes figures/drainage.* and data/drainage/, and returns the review record.
steps is the finer-step study on the finest mesh; sanity the 3D check.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_unsaturated_drainage as cd

KEY = 'drainage'


def _folder(root, case):
    return root/(f"h{case['spacing_m']:g}" + ('-explicit' if case['mode'] == 'explicit' else ''))


def build(study, output, record_source, steps, sanity=None, explicit_study=None):
    summary, fine = (pf.load_study(s, record_source) for s in (study, steps))
    s3d = pf.load_study(sanity, record_source) if sanity else None
    reference = cd.REF.solve(summary['reference_grid'], cd.TIMES)
    control_grid = 2*summary['reference_grid']
    control = cd.REF.solve(control_grid, cd.TIMES)
    scale = cd.RHO*cd.G*cd.PROPS['height']
    # Estimate reference uncertainty from its observed second-order spatial
    # convergence and the existing BDF absolute tolerances. This replay
    # allowance is not a rigorous bound on accumulated ODE integration error.
    reference_uncertainty = [
        4/3*float(abs(r['p']-np.interp(r['z'], c['z'], c['p'])).max())/scale
        + 2*cd.REF.pressure_atol_ratio for r, c in zip(reference, control)]
    for r, kept in zip(reference, summary['reference']):
        np.testing.assert_allclose(r['settlement'], kept['settlement'], rtol=1e-6)
    pt = sorted((c for c in summary['cases'] if c['mode'] == 'PT'), key=lambda c: -c['spacing_m'])
    # The mesh study's own explicit case is invalid (characteristic speed left at 0; see its NOTE-explicit.txt).
    explicit = [c for c in pf.load_study(explicit_study, record_source)['cases'] if c['mode'] == 'explicit'] if explicit_study else []
    finest, refined = pt[-1], fine['cases'][0]
    replay_difference = 0.
    for root, case in [(study, finest), (steps, refined)] + [(explicit_study, c) for c in explicit]:
        folder = _folder(root, case)
        for path in folder.glob('drain.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
        again = cd.compare(folder, case['spacing_m'], reference)
        actual = np.array([f['pressure_error'] for f in again])
        saved = np.array([f['pressure_error'] for f in case['frames']])
        for value, kept, uncertainty in zip(actual, saved, reference_uncertainty):
            np.testing.assert_allclose(value, kept, rtol=1e-6, atol=uncertainty)
        replay_difference = max(replay_difference, float(abs(actual-saved).max()))

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(_folder(study, finest)/'drain.cfg', data/'input.cfg')
    for name, root in (('summary.json', study), ('summary-steps.json', steps), ('summary-3d.json', sanity)):
        if root: shutil.copy2(root/'summary.json', data/name)
    shutil.copy2(study/'runner.py', data/'check_unsaturated_drainage.py')
    shutil.copy2(Path(__file__).with_name('drainage_reference.py'), data/'drainage_reference.py')

    height, scale = cd.PROPS['height'], cd.RHO*cd.G*cd.PROPS['height']
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        folder = _folder(study, finest)
        model = cd.Dynearthsol(str(folder/'drain')); rows = cd.scan_frames(str(folder/'drain'))
        coord = model.read_field(0, 'coordinate'); zn = coord[:, -1] + height
        column = np.isclose(coord[:, 0], 0)
        for color, ref in zip(ds.SERIES, reference):
            row = min(rows[1:], key=lambda r: abs(r['time'] - ref['time']))
            p = model.read_field(row['frame'], 'pore pressure')
            ds.analytical(ax, ref['p']/1e3, ref['z'], color=color, label=f"{ref['time']/60:g} min")
            ds.des(ax, p[column][::4]/1e3, zn[column][::4], color=color, label=None)
        ax.plot(cd.RHO*cd.G*(height - np.array([0, height]))/1e3, [0, height], ':', color=ds.MUTED, lw=.8, label='initial')
        ax.axvline(0, color=ds.MUTED, lw=.6)
        ds.panel(ax, 'a', 'Pore pressure as the column drains', xlabel='$p$ (kPa)', ylabel='height $z$ (m)')
        from matplotlib.lines import Line2D
        handles, labels = ax.get_legend_handles_labels()
        handles += [Line2D([], [], color=ds.INK, lw=1.2), Line2D([], [], marker='o', ls='none', ms=4, mfc='white', mec=ds.INK)]
        ax.legend(handles, labels + ['reference (numerical)', 'DES'], loc='upper right', fontsize=10)

        ax = axes[1]
        dense = np.concatenate([np.geomspace(5, 7200, 40)])
        curve = cd.REF.solve(100, dense)
        ds.analytical(ax, dense/60, [-1e3*r['settlement'] for r in curve], label='reference (numerical)')
        t, s = [], []
        for row in rows[1:]:
            e = cd.layer_strain(coord, model.read_field(0, 'connectivity'),
                                model.read_field(row['frame'], 'strain')[:, coord.shape[1]-1], finest['spacing_m'])
            t.append(row['time']); s.append(e)
        e0 = cd.layer_strain(coord, model.read_field(0, 'connectivity'),
                             model.read_field(0, 'strain')[:, coord.shape[1]-1], finest['spacing_m'])
        ds.des(ax, np.array(t)/60, [-1e3*float(np.sum(e - e0)*finest['spacing_m']) for e in s],
               label=f"DES, {finest['spacing_m']:g} m")
        ax.set_xscale('log')
        ds.panel(ax, 'b', 'Settlement of the top', xlabel='$t$ (min)', ylabel='settlement (mm)')

        ax = axes[2]
        h = np.array([c['spacing_m'] for c in pt])
        for color, k in zip(ds.SERIES, range(len(cd.TIMES))):
            ax.loglog(h, [100*c['frames'][k]['pressure_error'] for c in pt], 'o-', color=color, mfc='white',
                      label=f'{cd.TIMES[k]/60:g} min')
            ax.loglog([refined['spacing_m']], [100*refined['frames'][k]['pressure_error']], 'D', color=color, ms=4)
            for c in explicit:
                ax.loglog([c['spacing_m']], [100*c['frames'][k]['pressure_error']], 's', color=color, ms=5, mfc='none')
        ax.plot([], [], 'D', color=ds.MUTED, ms=4, label=f"{fine['steps']} steps")
        if explicit: ax.plot([], [], 's', color=ds.MUTED, ms=5, mfc='none', label='explicit form')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Pressure error by mesh', xlabel='mesh spacing $h$ (m)', ylabel=r'max $|p-p_{\rm ref}|/\rho_fgH$ (%)')
        ax.legend(fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)

    pct = lambda c, k: ' / '.join(f"{100*f[k]:.2f}" for f in c['frames'])
    s3 = s3d['cases'][0] if s3d else None
    return dict(key=KEY, theme='unsaturated', title='Drainage of a saturated column into the unsaturated state',
        problem=('A saturated sand column drains through its base: the top desaturates, pore pressure falls below zero and the '
                 'column settles under the Bishop effective stress and its lighter weight, as in the Liakopoulos experiment. '
                 'No closed form exists; an independent solver of the same continuum supplies the reference.'),
        result=(f"Largest nodal pressure error on the {finest['spacing_m']:g} m mesh with {summary['steps']} steps at 5, 20, 60 and 120 min: "
                f"{pct(finest, 'pressure_error')}% of $\\rho_fgH$; settlement {pct(finest, 'settlement_error')}%. "
                f"With {refined['steps']} steps: {pct(refined, 'pressure_error')}% and {pct(refined, 'settlement_error')}%."),
        scope=('2D, PT, ' + f"{summary['steps']} steps; pressure errors on meshes "
               + '; '.join(f"{c['spacing_m']:g} m: {pct(c, 'pressure_error')}%" for c in pt)
               + '. Late-time pressure converges at about second order with the mesh; the first frame and the settlement '
               'are limited by the step, since the base opens at once and the sequential coupling lags one step. '
               + ''.join(f"The explicit strain form on the {c['spacing_m']:g} m mesh gives {pct(c, 'pressure_error')}% "
                         f"and settlement {pct(c, 'settlement_error')}%. " for c in explicit)
               + (f"3D check ({s3['spacing_m']:g} m, PT): pressure {pct(s3, 'pressure_error')}%, settlement "
                  f"{pct(s3, 'settlement_error')}%. " if s3 else '3D is queued for a larger machine. ')
               + f"Reference: {summary['reference_grid']} cells, BDF; grid differences fall at "
               'second order. (a) pressure profiles; (b) settlement history, DES saved every 5 min; (c) pressure error: circles '
               f"{summary['steps']} steps, diamonds {fine['steps']} steps" + (', open squares the explicit form.' if explicit else '.')),
        meaning=('This checks Richards flow and the Bishop coupling together during drainage: desaturation, suction, the '
                 'effective-stress change and the lighter weight of the drying column.'),
        metric=refined['pressure_error'], metric_label='pressure error', dimension=2, spacing_m=refined['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 reference=summary['reference'], cases=pt, refined=refined, explicit=explicit,
                                 sanity_3d=s3,
                                 reference_pressure_atol_ratio=cd.REF.pressure_atol_ratio,
                                 reference_control_grid=control_grid,
                                 reference_metric_replay_budget=reference_uncertainty,
                                 reference_replay_metric_difference=replay_difference), indent=1))
