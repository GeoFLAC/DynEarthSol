"""Factor of safety of a homogeneous embankment by strength reduction, against Bishop.

Interface as qualify_theis: build(study, output, record_source) re-reads the
largest stable trial on the finest mesh from raw saves, reuses its recorded
Bishop reference, writes figures/embankment.* and
data/embankment/, and returns the review record.
"""
import json
import math
from pathlib import Path
import shutil

import numpy as np

import prepare_figures as pf
import diagram_style as ds
import bishop_reference as br
import strength_reduction_slope as srm

KEY = 'embankment'


def _bracket(case):
    """Largest F at rest and smallest F that keeps sliding; every F in between must be ordered."""
    rest = [t['factor'] for t in case['trials'] if t['at_rest']]
    slide = [t['factor'] for t in case['trials'] if not t['at_rest']]
    low = max(rest); high = min(f for f in slide if f > low)
    assert all(f < high for f in rest), 'a trial above the first sliding one came to rest'
    return low, high


def _displacement(folder, name):
    """Nodal displacement at the end of a trial: saved velocities integrated by the trapezoid rule (fixed mesh)."""
    model = srm.Dynearthsol(str(folder/name)); rows = srm.scan_frames(str(folder/name))
    x0 = model.read_field(rows[0]['frame'], 'coordinate'); conn = model.read_field(rows[0]['frame'], 'connectivity')
    u = np.zeros_like(x0); v0 = model.read_field(rows[0]['frame'], 'velocity'); t0 = rows[0]['time']
    for row in rows[1:]:
        v1 = model.read_field(row['frame'], 'velocity'); u += (v0 + v1)/2*(row['time'] - t0); v0, t0 = v1, row['time']
    return x0, conn, u


def build(study, output, record_source):
    summary = pf.load_study(study, record_source)
    fs_bishop = summary['fs_bishop']
    xc, zc, R = summary['bishop_circle']
    assert fs_bishop > 0 and R > 0 and np.isfinite([fs_bishop, xc, zc, R]).all()
    cases = sorted(summary['cases'], key=lambda c: -c['spacing_m'])
    finest = cases[-1]
    folder = study/f"h{finest['spacing_m']:g}"
    brackets = [_bracket(c) for c in cases]
    low, high = brackets[-1]
    show = f'srf{high:.3f}'
    for path in folder.iterdir():
        if path.suffix in ('.cfg', '.poly') or path.name.startswith(('seed.', 'gravity.', show + '.')):
            record_source(path)
    shown = next(t for t in finest['trials'] if t['factor'] == high)
    again, _, _ = srm.analyse(folder, show)
    np.testing.assert_allclose([f['max_displacement'] for f in again], [f['max_displacement'] for f in shown['frames']], rtol=1e-12)

    data = output/'data'/KEY; data.mkdir(parents=True, exist_ok=True)
    shutil.copy2(folder/f'{show}.cfg', data/'input.cfg')
    for name in ('seed.cfg', 'gravity.cfg', 'slope.poly'): shutil.copy2(folder/name, data/name)
    shutil.copy2(study/'summary.json', data/'summary.json'); shutil.copy2(study/'runner.py', data/'strength_reduction_slope.py')
    shutil.copy2(Path(br.__file__), data/'bishop_reference.py')

    with pf.plt.rc_context(ds.FIGURE_RC):
        fig, axes = ds.result_figure(3, width=13.0)
        ax = axes[0]
        x, conn, u = _displacement(folder, show)
        tri = pf.matplotlib.tri.Triangulation(x[:, 0], x[:, 1], conn)
        art = ax.tripcolor(tri, np.linalg.norm(u, axis=1)*1e3, cmap='magma_r', shading='gouraud', vmin=0)
        th = np.linspace(-math.pi, 0, 400)
        cx, cz = xc + R*np.cos(th), zc + R*np.sin(th)
        keep = (cx >= 0) & (cx <= srm.WIDTH) & (cz <= br.embankment().surface(cx))
        ax.plot(cx[keep], cz[keep], '--', color=ds.WATER, lw=1.4, label=f'Bishop critical circle, FS = {fs_bishop:.3f}')
        ax.set_aspect('equal'); ax.set_xlim(0, srm.WIDTH); ax.set_ylim(br.BASE, 2)
        cax = ax.inset_axes([.04, .72, .26, .05])                 # the air above the toe platform
        bar = fig.colorbar(art, cax=cax, orientation='horizontal'); bar.set_label('displacement (mm)', fontsize=10)
        bar.ax.xaxis.set_label_position('top')
        bar.ax.tick_params(labelsize=7.5)
        ds.panel(ax, 'a', f"Sliding mechanism at F = {high:g} ({finest['spacing_m']:g} m mesh)", xlabel='$x$ (m)', ylabel='$z$ (m)')
        ax.legend(fontsize=10, loc='lower left')

        ax = axes[1]
        trials = sorted(finest['trials'], key=lambda t: t['factor'])
        colors = pf.matplotlib.colormaps['viridis'](np.linspace(0, .9, len(trials)))
        for color, t in zip(colors, trials):
            tt = [f['time'] - t['frames'][0]['time'] for f in t['frames']]
            ax.semilogy(tt, [f['max_displacement']*1e3 for f in t['frames']], '-' if t['at_rest'] else '--', color=color, lw=1.1,
                        label=f"F = {t['factor']:g}")
        ax.plot([], [], '--', color=ds.MUTED, label='still sliding')
        ds.panel(ax, 'b', f"Largest displacement in time ({finest['spacing_m']:g} m)", xlabel='mass-scaled time (s)',
                 ylabel='displacement (mm)')
        ax.legend(fontsize=11, ncol=2, loc='upper center', bbox_to_anchor=(.5, -.25))

        ax = axes[2]
        for color, case, (lo, hi) in zip(ds.SERIES[1:] + [ds.HOT], cases, brackets):
            ts = sorted(case['trials'], key=lambda t: t['factor'])
            ax.semilogy([t['factor'] for t in ts], [t['frames'][-1]['max_displacement']*1e3 for t in ts], '-', color=color, lw=1.0,
                        label=f"DES, {case['spacing_m']:g} m: FS {lo:g}–{hi:g}")
            for t in ts:
                ax.plot(t['factor'], t['frames'][-1]['max_displacement']*1e3, 'o' if t['at_rest'] else 's', color=color,
                        mfc='white' if t['at_rest'] else color, ms=4)
        ax.plot([], [], 'o', color=ds.INK, mfc='white', ms=4, label='at rest'); ax.plot([], [], 's', color=ds.INK, ms=4, label='sliding')
        ax.axvline(fs_bishop, color=ds.INK, ls='--', lw=1.0, label=f'Bishop, {fs_bishop:.3f}')
        ds.panel(ax, 'c', f"Displacement after {summary['duration']:g} s by $F$", xlabel='strength reduction factor $F$',
                 ylabel='displacement (mm)')
        ax.legend(fontsize=9.5, loc='upper left')
        pf.save_figure(fig, output/'figures', KEY)

    mid = (low + high)/2
    mids = [(lo + hi)/2 for lo, hi in brackets]
    d1, d2 = mids[-3] - mids[-2], mids[-2] - mids[-1]
    # Richardson with the observed ratio of successive changes (meshes halve); only when they shrink monotonically.
    limit = mids[-1] - d2/(d1/d2 - 1) if len(mids) >= 3 and 0 < d2 < d1 else None
    return dict(key=KEY, theme='mechanics', title='Factor of safety of an embankment by strength reduction',
        problem=('The standard slope-stability problem of geotechnical practice: how far can the strength of a homogeneous '
                 'embankment be reduced before it slides? Strength reduction answers it with a deforming model that finds its own '
                 'mechanism; limit-equilibrium methods of slices answer it on a trial slip circle.'),
        result=(f"On the {finest['spacing_m']:g} m mesh the slope comes to rest at F = {low:g} and keeps sliding at "
                f"F = {high:g}: FS = {mid:.3f} ± {(high - low)/2:.3f}, against {fs_bishop:.3f} from Bishop’s simplified method "
                f"({100*(mid/fs_bishop - 1):+.1f}%). "
                + (f"The midpoints fall by {d1:.3f} and {d2:.3f} as the mesh halves, toward about {limit:.3f} "
                   f"({100*(limit/fs_bishop - 1):+.1f}%), slightly below a circular mechanism as expected. " if limit else '')
                + 'The mechanism passes through the toe, as Bishop’s critical circle does.'),
        scope=(f"2D plane strain, dry, $c=10$ kPa, $\\phi=20^\\circ$, $\\psi=0$, $E=100$ MPa, $\\nu=0.3$. Gravity by two elastic PT "
               f"steps; each trial then runs explicitly (no PT) for {summary['duration']:g} s of mass-scaled time (inertial scaling "
               f"$10^{{{round(math.log10(summary['scaling']))}}}$, DES's default damping {summary['damping']:g}) on the fixed original mesh; it is coming to rest when its peak "
               "speed over the last frames is below that three quarters through the run, and sliding otherwise. Brackets by mesh: "
               + '; '.join(f"{c['spacing_m']:g} m: {b[0]:g}–{b[1]:g}" for c, b in zip(cases, brackets))
               + '. Displacement at rest on the finest mesh: '
               + ', '.join(f"F = {t['factor']:g}: {t['frames'][-1]['max_displacement']*1e3:.1f} mm"
                           for t in sorted(finest['trials'], key=lambda t: t['factor']) if t['at_rest'])
               + f". Bishop's circle: centre ({xc:.1f}, {zc:.1f}) m, radius {R:.1f} m; converged in the number of slices. "
               'Under elastic gravity ($K_0=0.43$) the deep foundation yields slightly at $F=1$; it does not take part in the '
               'mechanism. (a) mechanism; (b) displacement in time; (c) displacement by $F$.'),
        meaning=('This checks Mohr–Coulomb plasticity and explicit quasi-static equilibrium in the problem slope engineers solve every day: the '
                 'mechanism is found, not assumed, and its factor of safety agrees with limit equilibrium. The hydromechanical '
                 'slope cases build on it.'),
        metric=abs(mid/fs_bishop - 1), metric_label='FS vs Bishop', dimension=2, spacing_m=finest['spacing_m'],
        evidence=json.dumps(dict(executable_sha256=summary['executable_sha256'], source_commit=summary['source_commit'],
                                 fs_bishop=fs_bishop, circle=[xc, zc, R], brackets=brackets, extrapolated=limit,
                                 cases=[dict(spacing_m=c['spacing_m'], trials=[dict(factor=t['factor'], at_rest=t['at_rest'],
                                        max_displacement=t['frames'][-1]['max_displacement'])
                                        for t in c['trials']]) for c in cases]), indent=1))
