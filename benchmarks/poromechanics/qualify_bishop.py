"""Validate the recorded DES Bishop-column studies (U8) and plot them.

Interface as qualify_theis. build(study, output, record_source, space) takes
the 2D study (all laws, every mesh) and the 3D study as the dimension check.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_bishop_column as cb

KEY = 'bishop'
NAMES = dict(van_genuchten='van Genuchten', brooks_corey='Brooks–Corey', gardner='Gardner',
             layered='van Genuchten / Gardner')


def _analytical(ax, x, y, color, law):
    """Analytical curve; the two-law column is dashed so it stands apart from Gardner."""
    ax.plot(x, y, '--' if law == 'layered' else '-', color=color, lw=1.2, label=NAMES[law])


def build(study, output, record_source, space):
    summary, other = pf.load_study(study, record_source), pf.load_study(space, record_source)
    laws = list(dict.fromkeys(c['law'] for c in summary['cases']))
    finest = {law: min((c for c in summary['cases'] if c['law'] == law), key=lambda c: c['spacing_m']) for law in laws}
    for law, case in finest.items():
        folder = study/f"{law}-h{case['spacing_m']:g}"
        for path in folder.glob('bishop.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
        again = cb.compare(folder, case['spacing_m'], law)
        np.testing.assert_allclose(again['u_des'], case['u_des'], rtol=0, atol=1e-15)
        np.testing.assert_allclose(again['u_ref'], case['u_ref'], rtol=1e-12, atol=0)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    law = laws[0]
    shutil.copy2(study/f"{law}-h{finest[law]['spacing_m']:g}"/'bishop.cfg', data/'input.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(space/'summary.json', data/'summary-3d.json')
    shutil.copy2(study/'runner.py', data/'check_bishop_column.py')
    shutil.copy2(Path(__file__).with_name('unsaturated_reference.py'), data/'unsaturated_reference.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        for color, law in zip(ds.SERIES, laws):
            c = finest[law]; z = np.array(c['z'])
            zf = np.linspace(0, cb.H, 400)
            _analytical(ax, 1e3*cb.reference(law, zf), zf, color, law)
            pick = slice(None, None, max(1, len(z)//10))
            ax.plot(1e3*np.array(c['u_des'])[pick], z[pick], 'o', ms=3.5, mfc='white', mec=color, mew=.9)
        for table, name in ((cb.TABLE0, 'start'), (cb.TABLE1, 'end')):
            ax.axhline(table, color=ds.WATER, lw=.8, ls=':')
            ax.text(.98, table + .12, f'water table, {name}', transform=ax.get_yaxis_transform(),
                    ha='right', va='bottom', color=ds.WATER, fontsize=10.5)
        ds.panel(ax, 'a', 'Settlement after drainage', xlabel='$u_z$ (mm)', ylabel='$z$ (m)')
        ds.key(ax)
        ax = axes[1]
        for color, law in zip(ds.SERIES, laws):
            c = finest[law]; z = np.array(c['z'])
            mid = .5*(z[1:] + z[:-1])
            eps = np.diff(c['u_des'])/np.diff(z)
            zf = np.linspace(0, cb.H, 2001)
            ref = np.gradient(cb.reference(law, zf), zf)
            _analytical(ax, 1e3*ref, zf, color, law)
            pick = slice(None, None, max(1, len(mid)//20))
            ax.plot(1e3*eps[pick], mid[pick], 'o', ms=3, mfc='white', mec=color, mew=.8)
        ax.axvline(0, color='#9aa5ab', lw=.8)
        ds.panel(ax, 'b', 'Vertical strain change', xlabel=r'$\Delta\varepsilon_{zz}$ ($10^{-3}$)', ylabel='$z$ (m)')
        ds.key(ax)
        ax = axes[2]
        # Colour marks the law; 2D is a line with dots, 3D open squares, which the 2D line passes through.
        for color, law in zip(ds.SERIES, laws):
            for s, style in ((summary, dict(marker='o', ls='-', ms=3.5)), (other, dict(marker='s', ls='none', ms=7))):
                cs = sorted((c for c in s['cases'] if c['law'] == law), key=lambda c: -c['spacing_m'])
                if not cs: continue
                h = np.array([c['spacing_m'] for c in cs])
                if law == 'layered' and s is summary: style = dict(style, ls='--')
                ax.loglog(h, [c['max_error'] for c in cs], color=color, mfc='none' if s is other else 'white',
                          mew=.9, label=NAMES[law] if s is summary else None, **style)
        cs = sorted((c for c in summary['cases'] if c['law'] == laws[0]), key=lambda c: -c['spacing_m'])
        h = np.array([c['spacing_m'] for c in cs])
        ax.loglog(h, cs[0]['max_error']*(h/h[0])**2, '--', color=ds.MUTED, lw=.8, label='slope 2')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Mesh convergence', xlabel='mesh spacing $h$ (m)', ylabel='max error / max $|u_z|$')
        from matplotlib.lines import Line2D
        handles, labels = ax.get_legend_handles_labels()
        handles += [Line2D([], [], color=ds.INK, marker='o', mfc='white', ms=3.5, ls='-'),
                    Line2D([], [], color=ds.INK, marker='s', mfc='none', ms=7, ls='none')]
        labels += [f"{summary['ndims']}D", f"{other['ndims']}D"]
        ax.legend(handles, labels, fontsize=10, loc='upper left')
        pf.save_figure(fig, output/'figures', KEY)

    finest3 = {law: min((c for c in other['cases'] if c['law'] == law), key=lambda c: c['spacing_m'])
               for law in dict.fromkeys(c['law'] for c in other['cases'])}
    single = [k for k in laws if k != 'layered']
    worst = max(finest[k]['max_error'] for k in single)
    worst3 = max(finest3[k]['max_error'] for k in single if k in finest3)
    drained = max(c['drained_residual'] for s in (summary, other) for c in s['cases'])
    hmin = min(c['spacing_m'] for c in finest.values())
    orders = {law: float(np.polyfit(np.log([c['spacing_m'] for c in summary['cases'] if c['law'] == law]),
                                    np.log([c['max_error'] for c in summary['cases'] if c['law'] == law]), 1)[0])
              for law in laws}
    return dict(key=KEY, theme='unsaturated', title='Bishop effective stress: lowering the water table',
        problem=('The water table under a laterally confined elastic column drops. Above it the soil drains; '
                 'suction weighted by the effective saturation tightens the skeleton, and the drained soil becomes lighter.'),
        result=(f'Largest settlement error on the {hmin:g} m mesh, relative to the top settlement: {100*worst:.3f}% in 2D and '
                f'{100*worst3:.3f}% in 3D for the ' + ', '.join(NAMES[k] for k in single[:-1]) + f' and {NAMES[single[-1]]} laws'
                + (f"; {100*finest['layered']['max_error']:.3f}% and {100*finest3['layered']['max_error']:.3f}% for the column "
                   'with two laws, where saturation jumps at the shared interface nodes.' if 'layered' in finest and 'layered' in finest3 else '.')),
        scope=(f"Column {cb.H:g} m high, water table from {cb.TABLE0:g} m to {cb.TABLE1:g} m; PT mechanics with pressure feedback; "
               f"both states drained to within {drained:.0e} of the hydrostatic pressure change. "
               'Error orders on 2D meshes: ' + ', '.join(f'{NAMES[k]} {v:.2f}' for k, v in orders.items()) + '. '
               "(a) settlement; (b) strain change; (c) mesh convergence, 2D lines and 3D squares."),
        meaning=('This checks the Bishop coupling: the effective-stress increment from each element\'s own retention law, '
                 'including nodes shared by two laws, and gravity on the partly saturated bulk weight. '
                 'Drained end states do not test the transient or the coupling term in the flow equation.'),
        metric=worst, metric_label='settlement error / max |u|', dimension=2, spacing_m=hmin,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 orders=orders, cases_2d=[{k: v for k, v in c.items() if k not in ('z', 'u_des', 'u_ref')} for c in summary['cases']],
                                 cases_3d=[{k: v for k, v in c.items() if k not in ('z', 'u_des', 'u_ref')} for c in other['cases']]), indent=1))
