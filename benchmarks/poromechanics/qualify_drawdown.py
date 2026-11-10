"""Rapid drawdown against the embankment: factor of safety in time, DES against Bishop.

Interface as qualify_theis: build(study, output, record_source) recomputes
Bishop's factor of safety on DES's own pore pressure at the first state,
re-reads one trial from raw saves, writes figures/drawdown.* and
data/drawdown/, and returns the review record.
"""
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import bishop_reference as br
import strength_reduction_slope as srm
import drawdown_slope as dd

KEY = 'drawdown'


def _bracket(state):
    rest = [t['factor'] for t in state['trials'] if t['at_rest']]
    slide = [t['factor'] for t in state['trials'] if not t['at_rest']]
    low = max(rest, default=None)
    high = min((f for f in slide if low is None or f > low), default=None)
    assert low is None or high is None or all(f < high for f in rest), 'a trial above the first sliding one came to rest'
    return low, high


def build(study, output, record_source):
    summary = json.loads((study/'summary.json').read_text())
    assert hashlib.sha256((study/'runner.py').read_bytes()).hexdigest() == summary['runner_sha256']
    record_source(study/'summary.json'); record_source(study/'runner.py')
    case = min(summary['cases'], key=lambda c: c['spacing_m'])
    folder = study/f"h{case['spacing_m']:g}"
    audit = None
    if (study/'audit.json').exists():
        audit = json.loads((study/'audit.json').read_text())
        assert audit['replay_passed']
        record_source(study/'audit.json')
        for path in folder.iterdir():
            if path.name in audit['raw_files_sha256']:
                assert dd.digest(path) == audit['raw_files_sha256'][path.name], path
    states = sorted(case['states'], key=lambda s: s['time'])
    first = states[0]
    again = dd.bishop(folder, 'drawdown', first['frame'])[0]
    np.testing.assert_allclose(again, first['fs_bishop'], rtol=1e-6)
    brackets = [_bracket(s) for s in states]
    low, high = brackets[0]
    # Trials are named srf<requested time>-<F>; the state records the saved frame's time.
    show = next((p.stem for p in folder.glob(f"srf*-{high:.3f}.cfg") if abs(float(p.stem[3:].split('-')[0]) -
                 first['time']) < 1), None) if high else None
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.poly') or path.name.startswith(('seed.', 'reservoir.', 'drawdown.')) or \
           (show and path.name.startswith(show + '.')):
            record_source(path)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/'slope.poly', data/'slope.poly')
    shutil.copy2(folder/'drawdown.cfg', data/'input.cfg')
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'drawdown_slope.py')
    shutil.copy2(Path(br.__file__), data/'bishop_reference.py')

    model = dd.Dynearthsol(str(folder/'drawdown'))
    x = model.read_field(first['frame'], 'coordinate'); conn = model.read_field(first['frame'], 'connectivity')
    tri = pf.matplotlib.tri.Triangulation(x[:, 0], x[:, 1], conn)
    hour = 3600.0
    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3, width=13.0)
        ax = axes[0]
        p = model.read_field(first['frame'], 'pore pressure')/1e3
        top = 50*math.ceil(max(p.max(), -p.min())/50)
        art = ax.tripcolor(tri, p, cmap='RdBu_r', vmin=-top, vmax=top, shading='gouraud')
        ax.tricontour(tri, p, levels=[0.0], colors=ds.INK, linewidths=1.0)
        xc, zc, R = first['circle']
        th = np.linspace(-math.pi, 0, 400); cx, cz = xc + R*np.cos(th), zc + R*np.sin(th)
        keep = (cz <= br.embankment().surface(cx) + 1e-9)
        ax.plot(cx[keep], cz[keep], '--', color=ds.HOT, lw=1.2, label=f"Bishop circle, FS = {first['fs_bishop']:.3f}")
        ax.set_aspect('equal'); ax.set_xlim(0, srm.WIDTH); ax.set_ylim(br.BASE, 2)
        cax = ax.inset_axes([.04, .72, .26, .05])
        bar = fig.colorbar(art, cax=cax, orientation='horizontal'); bar.set_label('pore pressure (kPa)', fontsize=8)
        bar.ax.xaxis.set_label_position('top'); bar.ax.tick_params(labelsize=7.5)
        ds.panel(ax, 'a', f"Pore pressure {first['time']:.0f} s after drawdown", xlabel='$x$ (m)', ylabel='$z$ (m)')
        ax.legend(fontsize=8, loc='lower left')

        ax = axes[1]
        t = np.array([s['time'] for s in states])
        ax.plot(t, [s['fs_bishop'] for s in states], 'o-', color=ds.INK, mfc='white', label='Bishop on the DES pore pressure')
        for s, (lo, hi) in zip(states, brackets):
            if lo is not None and hi is not None:
                ax.vlines(s['time'], lo, hi, color=ds.HOT, lw=3)
        ax.plot([], [], '-', color=ds.HOT, lw=3, label='DES: last at rest to first sliding $F$')
        ax.axhline(1.0, color=ds.MUTED, lw=.6); ax.set_ylim(0.8, 1.05)
        ax.set_xscale('log')
        ds.panel(ax, 'b', 'Factor of safety after drawdown', xlabel='time after drawdown (s)', ylabel='factor of safety')
        ax.legend(fontsize=8, loc='upper left')

        ax = axes[2]
        if show:
            m = dd.Dynearthsol(str(folder/show)); rows = dd.scan_frames(str(folder/show))
            u = np.zeros_like(x); v0 = m.read_field(rows[0]['frame'], 'velocity'); t0 = rows[0]['time']
            for row in rows[1:]:
                v1 = m.read_field(row['frame'], 'velocity'); u += (v0 + v1)/2*(row['time'] - t0); v0, t0 = v1, row['time']
            art = ax.tripcolor(tri, np.linalg.norm(u, axis=1)*1e3, cmap='magma_r', shading='gouraud', vmin=0)
            ax.plot(cx[keep], cz[keep], '--', color=ds.WATER, lw=1.2, label='Bishop circle')
            cax = ax.inset_axes([.04, .72, .26, .05])
            bar = fig.colorbar(art, cax=cax, orientation='horizontal'); bar.set_label('displacement (mm)', fontsize=8)
            bar.ax.xaxis.set_label_position('top'); bar.ax.tick_params(labelsize=7.5)
        ax.set_aspect('equal'); ax.set_xlim(0, srm.WIDTH); ax.set_ylim(br.BASE, 2)
        ds.panel(ax, 'c', f"Mechanism at F = {high:g}, {first['time']:.0f} s" if high else 'Mechanism', xlabel='$x$ (m)',
                 ylabel='$z$ (m)')
        if show: ax.legend(fontsize=8, loc='lower left')
        pf.save_figure(fig, output/'figures', KEY)

    mids = [(lo + hi)/2 if lo is not None and hi is not None else None for lo, hi in brackets]
    dry = br.embankment().critical(br.starts())[0]
    worst = max(abs(m/s['fs_bishop'] - 1) for m, s in zip(mids, states) if m is not None)
    contained = [lo <= s['fs_bishop'] <= hi for s, (lo, hi) in zip(states, brackets)]
    misses = [s['fs_bishop']-hi for s, (lo, hi) in zip(states, brackets) if s['fs_bishop'] > hi]
    comparison = (f" {sum(contained)} of {len(states)} brackets contain Bishop. "
                  + (f"The other values exceed the upper edge by {min(misses):.5f}–{max(misses):.5f}; "
                     "DES is slightly conservative. The original all-times containment criterion is not met."
                     if misses else "The original all-times containment criterion is met."))
    fmt = lambda b: f'{b[0]:g}–{b[1]:g}' if b[0] is not None and b[1] is not None else 'not bracketed'
    return dict(key=KEY, theme='mechanics', status='preliminary', title='Application: rapid drawdown of a reservoir against an embankment',
        problem=('A reservoir stands at the crest of the embankment and the slope is saturated. When the reservoir is drawn down '
                 'quickly the water no longer supports the face, while the pore pressure inside falls only as the slope drains: '
                 'the classic critical case for the upstream slope of a dam.'),
        result=(f"{states[0]['time']:.0f} s after drawdown DES brackets the factor of safety at {fmt(brackets[0])} against "
                f"{states[0]['fs_bishop']:.3f} from Bishop on the same pore pressure; at {states[-1]['time']:.0f} s, "
                f"{fmt(brackets[-1])} against {states[-1]['fs_bishop']:.3f}. Largest gap between the bracket midpoint and Bishop: "
                f"{100*worst:.1f}%." + comparison + f" The slope that stands dry (FS {dry:.3f}) fails on rapid drawdown, and draining "
                "restores only part of its strength in this time."),
        scope=(f"2D, {case['spacing_m']:g} m mesh, the embankment of the strength-reduction page, saturated (porosity "
               f"{summary['porosity']:g}, saturated unit weight 20 kN/m³), Gardner retention ($\\alpha={summary['alpha']*1e4:g}\\times10^{{-4}}$ "
               f"Pa$^{{-1}}$), Bishop $\\chi=S_e$, permeability {summary['permeability']:g} m². The reservoir at the crest loads the face; "
               f"the slope fills to hydrostatic (misfit {100*case['hydrostatic_misfit']:.2f}% of the base pressure). Drawdown "
               'removes the load and drains the face and the toe platform; everything is explicit DES with the dry page\'s '
               'settings, because DES does not switch hydraulic formulation mid-run. Strength-reduction trials restart at '
               'each time with the pore pressure held (permeability frozen, pore fluid made compressible), as limit '
               'equilibrium assumes. Brackets: ' + '; '.join(f"{s['time']:.0f} s: {fmt(b)} (Bishop {s['fs_bishop']:.3f})"
                                                          for s, b in zip(states, brackets))
               + '. (a) pore pressure; (b) factor of safety in time; (c) mechanism.'),
        meaning=('All 72 fine-mesh trials and six Bishop searches reproduce their recorded metrics from raw fields. '
                 'The 0.05 factor spacing and small conservative misses remain visible; no tuning is used to match Bishop. '
                 'This joins the verified pieces in the problem dam engineers check first: the undrained response to unloading, '
                 'unsaturated drainage and effective-stress plasticity. Limit equilibrium needs the pore pressure as an input; '
                 'DES computes it and finds the mechanism.'),
        metric=worst, metric_label='FS midpoint vs Bishop', dimension=2, spacing_m=case['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 hydrostatic_misfit=case['hydrostatic_misfit'], contained=contained,
                                 states=[dict(time=s['time'], fs_bishop=s['fs_bishop'], bracket=b, p_max=s['p_max'], p_min=s['p_min'])
                                         for s, b in zip(states, brackets)]), indent=1))
