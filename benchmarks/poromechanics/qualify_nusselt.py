"""Validate the recorded DES steady-convection study and plot it.

Interface as qualify_theis: build(study, output, record_source) recomputes Nu
from the raw fields of the finest runs, writes figures/nusselt.* and
data/nusselt/, and returns the review record.
"""
import json
from pathlib import Path
import shutil
import sys

import numpy as np

import prepare_figures as pf
import diagram_style as ds
from convection_steady import solve
sys.path.insert(0, str(pf.ROOT/'tests/functional/hydromechanics'))
import check_convection_nusselt as cn

KEY = 'nusselt'


def build(study, output, record_source):
    summary = pf.load_study(study, record_source)
    ras = sorted({c['ra'] for c in summary['cases']})
    finest = {ra: min((c for c in summary['cases'] if c['ra'] == ra), key=lambda c: c['spacing_m']) for ra in ras}
    for ra, case in finest.items():
        folder = study/f"ra{ra:g}-h{case['spacing_m']:g}"
        for path in folder.glob('box.*'):
            if path.suffix in ('.cfg', '.info') or '.save.' in path.name: record_source(path)
        again = cn.nusselt_des(folder, ra)[-1]['nu']
        np.testing.assert_allclose(again, case['nu'], rtol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    top = max(ras)
    shutil.copy2(study/f"ra{top:g}-h{finest[top]['spacing_m']:g}"/'box.cfg', data/'input.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'check_convection_nusselt.py')
    shutil.copy2(Path(__file__).with_name('convection_steady.py'), data/'convection_steady.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3)
        ax = axes[0]
        folder = study/f"ra{top:g}-h{finest[top]['spacing_m']:g}"
        model = cn.cc.Dynearthsol(str(folder/'box')); row = cn.cc.scan_frames(str(folder/'box'))[-1]
        coord = model.read_field(0, 'coordinate'); conn = model.read_field(0, 'connectivity')
        theta = (model.read_field(row['frame'], 'temperature') - cn.cc.T_TOP)/cn.cc.DT
        levels = np.linspace(.1, .9, 9)
        tri = pf.matplotlib.tri.Triangulation(coord[:, 0], coord[:, 1] + 1, conn)
        ax.tricontour(tri, theta, levels=levels, colors=ds.HOT, linewidths=1.0)
        _, _, ref, _ = solve(top, 128)
        g = np.linspace(0, 1, ref.shape[0])
        # The reference cell turns either way; match the DES sense of rotation by the x-mirror if needed.
        if np.corrcoef(theta, ref[np.clip(np.round((coord[:, 1]+1)*128).astype(int), 0, 128),
                                     np.clip(np.round(coord[:, 0]*128).astype(int), 0, 128)])[0, 1] < .5:
            ref = ref[:, ::-1]
        ax.contour(g, g, ref, levels=levels, colors=ds.INK, linewidths=.8, linestyles='--')
        ax.set_aspect('equal')
        # Isotherms fill the cell, so the key goes in the title.
        ds.panel(ax, 'a', f'Isotherms at Ra = {top:g}: DES solid, reference dashed', xlabel='$x/H$', ylabel='$z/H$')
        ax = axes[1]
        for color, ra in zip((ds.WATER, ds.HOT), ras):
            series = finest[ra]['series']
            ax.plot([f['time']/1e6 for f in series], [f['nu'] for f in series], 'o-', ms=3, color=color, mfc='white',
                    label=f'Ra = {ra:g}, DES')
            ax.axhline(finest[ra]['nu_ref'], color=color, lw=.8, ls='--')
        ax.plot([], [], '--', color=ds.MUTED, lw=.8, label='steady reference')
        # The volume form measures steady transport; early storage-driven flow sends it far off scale.
        ax.set_ylim(0, 1.3*max(finest[ra]['nu_ref'] for ra in ras))
        ds.panel(ax, 'b', 'Approach to steady transport', xlabel=r'$t$ (Ms)', ylabel='Nu (volume form)')
        ax.legend(fontsize=10)
        ax = axes[2]
        for color, ra in zip((ds.WATER, ds.HOT), ras):
            cs = sorted((c for c in summary['cases'] if c['ra'] == ra), key=lambda c: -c['spacing_m'])
            h = np.array([c['spacing_m'] for c in cs])
            ax.loglog(h, [c['error'] for c in cs], 'o-', color=color, mfc='white', label=f'Ra = {ra:g}, volume')
            ax.loglog(h, [c['wall_error'] for c in cs], 's:', color=color, mfc='white', label=f'Ra = {ra:g}, wall')
        ax.loglog(h, cs[0]['error']*(h/h[0])**2, '--', color=ds.MUTED, lw=.8, label='slope 2')
        ax.set_xticks(h, [f'{x:g}' for x in h]); ax.xaxis.set_minor_locator(pf.matplotlib.ticker.NullLocator())
        ds.panel(ax, 'c', 'Mesh convergence', xlabel='mesh spacing $h/H$', ylabel='|Nu error| / Nu')
        ax.legend(fontsize=10)
        pf.save_figure(fig, output/'figures', KEY)

    h = min(c['spacing_m'] for c in finest.values())
    orders = {}
    for ra in ras:
        cs = sorted((c for c in summary['cases'] if c['ra'] == ra), key=lambda c: -c['spacing_m'])
        orders[ra] = float(np.polyfit(np.log([c['spacing_m'] for c in cs]), np.log([c['error'] for c in cs]), 1)[0])
    ref = summary['reference']
    return dict(key=KEY, theme='heat', title='Steady hydrothermal convection: Nusselt number',
        problem=('Above onset the heated porous cell settles into a steady convection roll that carries several times the '
                 'conductive heat flow. Its Nusselt number has no closed form; an independent fine-grid solution supplies it.'),
        result=('DES Nusselt numbers on the ' + f'{h:g} H mesh: ' + '; '.join(
                f"Ra = {ra:g}: {finest[ra]['nu']:.4f} against {finest[ra]['nu_ref']:.4f} ({100*finest[ra]['error']:.2f}%)" for ra in ras) + '.'),
        scope=('Unit square cell, impermeable walls, insulated sides; the steady state is reached: Nu changes by less than '
               f"{max(c['steady_change'] for c in summary['cases']):.0e} (relative) at the end. Errors on meshes "
               + ', '.join(f'{c:g}' for c in sorted({c['spacing_m'] for c in summary['cases']}, reverse=True)) + ' H: '
               + '; '.join(f"Ra = {ra:g}: " + ', '.join(f"{100*c['error']:.2f}%" for c in sorted(
                   (c for c in summary['cases'] if c['ra'] == ra), key=lambda c: -c['spacing_m'])) for ra in ras)
               + '. At Ra = 100 the volume form flattens near 0.17% above the reference, likely from rebuilding the Darcy flux '
               'from saved fields; the wall-gradient estimate approaches from below ('
               + ', '.join(f"{100*finest[ra]['wall_error']:.2f}% at Ra = {ra:g}" for ra in ras) + ', on the finest mesh). Reference: '
               + '; '.join(f"Ra = {ra}: {r['nu_extrapolated']:.4f} (order {r['observed_order']:.2f} on " + ', '.join(f'${g}^2$' for g in r['grids']) + ' grids)' for ra, r in ref.items())
               + '. (a) isotherms; (b) approach to steady state; (c) mesh convergence.'),
        meaning=('This checks finite-amplitude Darcy convection: buoyancy, heat advection by the Darcy flux and conduction '
                 'together, through the heat they carry. Storage only sets the transient, so a compressible pore fluid '
                 'shortens the runs; the transient itself is not compared.'),
        metric=max(c['error'] for c in finest.values()), metric_label='Nu error', dimension=2, spacing_m=h,
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 orders=orders, reference=ref,
                                 cases=[{k: v for k, v in c.items() if k != 'series'} for c in summary['cases']]), indent=1))
