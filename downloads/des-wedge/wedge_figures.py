#!/usr/bin/env python3
"""Figures and numbers for the accretionary-wedge tutorial.

usage: python3 wedge_figures.py [rundir] [outdir]
       python3 wedge_figures.py compare out.png model:label[@frame] [model:label[@frame] ...]
       python3 wedge_figures.py sections out.png f1,f2,... model:label [model:label ...]

The first form reads wedge-phi1.* and wedge-phi5.* in rundir (default: .) and
writes wedge-evolution.png, wedge-strainrate.png and wedge-growth.png to outdir;
it prints the wedge width and height of every frame. The second form
compares any runs, e.g. the exercise variants, in one figure: plastic strain
of each at its last frame (or at @frame), then wedge width and height against shortening.
The third form draws plastic strain of each run (rows) at each listed frame
(columns), for models of any geometry.
Output is read with DynEarthSol's Dynearthsol.py, found through DES_DIR or in
../DynEarthSol.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('DES_DIR', os.path.join(HERE, '..', 'DynEarthSol')))
from Dynearthsol import Dynearthsol, YEAR2SEC

MODELS = [('wedge-phi1', 'basal friction 1°'), ('wedge-phi5', 'basal friction 5°')]
COLORS = {'wedge-phi1': '#2F6FB0', 'wedge-phi5': '#C2410C'}
TOP = 32            # bcflag of the top surface (BOUNDZ1)
VE = 2.0            # vertical exaggeration of the cross sections
ZMIN, ZMAX = -4.0, 3.2
SHOW = [5, 10, 15, 20]        # frames at 50, 100, 150, 200 kyr
plt.rcParams.update({'font.size': 9})


def frame(des, fr):
    d = {}
    d['t'] = des.time[des.frames.index(fr)] / YEAR2SEC / 1e3            # kyr
    d['xy'] = des.read_field(fr, 'coordinate') / 1e3                     # km
    d['conn'] = des.read_field(fr, 'connectivity')
    d['pls'] = des.read_field(fr, 'plastic strain')
    sr = des.read_field(fr, 'strain-rate')                               # xx, zz, xz
    d['sr'] = np.sqrt(0.25 * (sr[:, 0] - sr[:, 1])**2 + sr[:, 2]**2)     # second invariant
    d['mat'] = des.read_field(fr, 'material')
    flag = des.read_field(fr, 'bcflag')
    top = np.nonzero(flag & TOP)[0]
    top = top[np.argsort(d['xy'][top, 0])]
    d['surf'] = d['xy'][top]
    return d


def metrics(d):
    """Wall position, thrust front and maximum elevation.

    The thrust front is the farthest brittle-layer element above 2.5 km depth
    whose plastic strain exceeds 0.5: the frontal thrust, not the toe of the
    topography, which diffusion spreads into the foreland.
    """
    sx, sz = d['surf'][:, 0], d['surf'][:, 1]
    wall = sx[0]
    cx = d['xy'][d['conn'], 0].mean(axis=1)
    cz = d['xy'][d['conn'], 1].mean(axis=1)
    faulted = (d['pls'] > 0.5) & (d['mat'] < 0.5) & (cz > -2.5)
    front = cx[faulted].max() if faulted.any() else wall
    return dict(t=d['t'], wall=wall, front=front, width=front - wall, hmax=sz.max())


def interface(ax, d):
    """Draw the edges shared by a brittle-layer and a decollement element."""
    conn, mat = d['conn'], d['mat']
    edges = {}
    for e, (a, b, c) in enumerate(conn):
        for p, q in ((a, b), (b, c), (c, a)):
            k = (min(p, q), max(p, q))
            edges.setdefault(k, []).append(e)
    segs = [k for k, es in edges.items() if len(es) == 2 and mat[es[0]] != mat[es[1]]]
    xy = d['xy']
    for p, q in segs:
        ax.plot(xy[[p, q], 0], xy[[p, q], 1], color='#4FC3F7', lw=0.6)


def section(ax, d, field, **kw):
    tri = mtri.Triangulation(d['xy'][:, 0], d['xy'][:, 1], d['conn'])
    ax.set_facecolor('#F7F6F2')
    im = ax.tripcolor(tri, facecolors=field, rasterized=True, **kw)
    ax.plot(d['surf'][:, 0], d['surf'][:, 1], color='k', lw=0.8)
    ax.set_xlim(0, 40)
    ax.set_ylim(ZMIN - 0.1, ZMAX)
    ax.set_aspect(VE)
    return im


def compare(out, specs):
    pal = ['#2F6FB0', '#7C3AED', '#C2410C', '#059669', '#6B7280']
    n = len(specs)
    fig = plt.figure(figsize=(11, 2.3 * n + 1.2))
    top, bot = fig.subfigures(1, 2, width_ratios=[1.15, 1])
    axes = top.subplots(n, 1, sharex=True, squeeze=False)[:, 0]
    top.subplots_adjust(hspace=0.45)
    series = []
    for i, spec in enumerate(specs):
        spec, last = spec.rsplit('@', 1) if '@' in spec else (spec, None)
        name, label = spec.split(':', 1) if ':' in spec else (spec, spec)
        des = Dynearthsol(name)
        last = int(last) if last else des.frames[-1]
        d = frame(des, last)
        ax = axes[i]
        im = section(ax, d, d['pls'], cmap='magma', vmin=0, vmax=2)
        interface(ax, d)
        m = metrics(d)
        ax.set_title('%s, t = %.0f kyr' % (label, d['t']), loc='left', fontsize=9)
        ax.text(0.99, 0.95, '%.1f km wide, %.1f km high' % (m['width'], m['hmax']),
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
        ax.set_ylabel('z (km)')
        ax.set_xlim(max(0, m['wall'] - 1), 40)
        series.append((label, [metrics(d if fr == last else frame(des, fr))
                               for fr in des.frames[1:] if fr <= last], pal[i % len(pal)]))
    top.colorbar(im, ax=axes, location='bottom', shrink=0.8, aspect=30, extend='max', pad=0.08).set_label('plastic strain')
    for ax, (key, title) in zip(bot.subplots(2, 1, sharex=True), (('width', 'wedge width, wall to thrust front (km)'),
                                                     ('hmax', 'maximum elevation (km)'))):
        for label, ms, c in series:
            ax.plot([m['t'] * 0.05 for m in ms], [m[key] for m in ms], 'o-', ms=3, color=c, label=label)
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax.set_xlabel('shortening (km)')
    fig.savefig(out, dpi=130, bbox_inches='tight')


def sections(out, frames, specs):
    nr, nc = len(specs), len(frames)
    models = {}
    for spec in specs:
        name = spec.split(':', 1)[0]
        if name not in models:
            des = Dynearthsol(name)
            models[name] = (des, {fr: frame(des, fr) for fr in set(frames)})
    zmin = min((data[0]['xy'] if 0 in data else des.read_field(0, 'coordinate') / 1e3)[:, 1].min()
               for des, data in models.values())
    zmax = max(d['xy'][:, 1].max() for _, data in models.values() for d in data.values()) + 0.4
    fig, axes = plt.subplots(nr, nc, figsize=(6 * nc, 1.9 * nr + 0.8), squeeze=False,
                             sharex=True, sharey=True)
    for i, spec in enumerate(specs):
        name, label = spec.split(':', 1) if ':' in spec else (spec, spec)
        _, data = models[name]
        for j, fr in enumerate(frames):
            d = data[fr]
            ax = axes[i, j]
            im = section(ax, d, d['pls'], cmap='magma', vmin=0, vmax=2)
            interface(ax, d)
            ax.set_ylim(zmin - 0.1, max(ZMAX, zmax))
            ax.set_title('%s, t = %.0f kyr' % (label, d['t']), loc='left', fontsize=9)
            if j == 0:
                ax.set_ylabel('z (km)')
    for ax in axes[-1]:
        ax.set_xlabel('x (km)')
    fig.colorbar(im, ax=axes, shrink=0.6, extend='max', pad=0.015).set_label('plastic strain')
    fig.suptitle('Plastic strain, vertical exaggeration %g×; blue lines: decollement boundaries' % VE,
                 fontweight='bold')
    fig.savefig(out, dpi=130, bbox_inches='tight')


if len(sys.argv) > 1 and sys.argv[1] == 'sections':
    sections(sys.argv[2], [int(f) for f in sys.argv[3].split(',')], sys.argv[4:])
    sys.exit(0)

if len(sys.argv) > 1 and sys.argv[1] == 'compare':
    compare(sys.argv[2], sys.argv[3:])
    sys.exit(0)

rundir = sys.argv[1] if len(sys.argv) > 1 else '.'
outdir = sys.argv[2] if len(sys.argv) > 2 else '.'
runs = {}
for name, _ in MODELS:
    des = Dynearthsol(os.path.join(rundir, name))
    history, snapshots = [], {}
    for fr in des.frames:
        d = frame(des, fr)
        history.append(metrics(d))
        if fr in SHOW:
            snapshots[fr] = d
    runs[name] = (des, history, snapshots)

for name, label in MODELS:
    print('\n%s (%s)' % (name, label))
    print('   t(kyr)  wall(km)  front(km)  width(km)  height(km)')
    for m in runs[name][1]:
        print('  %6.0f  %8.2f  %9.2f  %9.2f  %10.2f' % (m['t'], m['wall'], m['front'], m['width'], m['hmax']))

# --- plastic strain ----------------------------------------------------------
fig, axes = plt.subplots(len(SHOW), 2, figsize=(12, 7.8), sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.45)
for j, (name, label) in enumerate(MODELS):
    for i, fr in enumerate(SHOW):
        d = runs[name][2][fr]
        ax = axes[i, j]
        im = section(ax, d, d['pls'], cmap='magma', vmin=0, vmax=2)
        interface(ax, d)
        m = metrics(d)
        ax.text(0.99, 0.95, 'thrust front %.1f km from the wall, %.1f km high' % (m['width'], m['hmax']),
                transform=ax.transAxes, ha='right', va='top', fontsize=8,
                bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))
        if j == 0:
            ax.set_ylabel('z (km)')
        ax.set_title(('%s\n' % label if i == 0 else '') + 't = %.0f kyr, %.1f km shortening' % (d['t'], d['t'] * 0.05),
                     loc='left', fontsize=9, fontweight='bold' if i == 0 else 'normal')
for ax in axes[-1]:
    ax.set_xlabel('x (km)')
cb = fig.colorbar(im, ax=axes, shrink=0.6, extend='max', pad=0.015)
cb.set_label('plastic strain')
fig.suptitle('Plastic strain, vertical exaggeration 2×; blue line: top of the decollement', fontweight='bold')
fig.savefig(os.path.join(outdir, 'wedge-evolution.png'), dpi=130, bbox_inches='tight')
plt.close(fig)

# --- strain rate -------------------------------------------------------------
fig, axes = plt.subplots(len(SHOW), 2, figsize=(12, 7.8), sharex=True, sharey=True)
fig.subplots_adjust(hspace=0.45)
for j, (name, label) in enumerate(MODELS):
    for i, fr in enumerate(SHOW):
        d = runs[name][2][fr]
        ax = axes[i, j]
        im = section(ax, d, np.log10(np.maximum(d['sr'], 1e-20)), cmap='viridis', vmin=-15.5, vmax=-12.5)
        if j == 0:
            ax.set_ylabel('z (km)')
        ax.set_title(('%s\n' % label if i == 0 else '') + 't = %.0f kyr' % d['t'],
                     loc='left', fontsize=9, fontweight='bold' if i == 0 else 'normal')
for ax in axes[-1]:
    ax.set_xlabel('x (km)')
cb = fig.colorbar(im, ax=axes, shrink=0.6, extend='both', pad=0.015)
cb.set_label('log10 strain rate (1/s)')
fig.suptitle('Strain rate: the faults active at each moment, vertical exaggeration 2×', fontweight='bold')
fig.savefig(os.path.join(outdir, 'wedge-strainrate.png'), dpi=130, bbox_inches='tight')
plt.close(fig)

# --- growth ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
ax = axes[0]
for name, label in MODELS:
    for fr, ls in ((10, '--'), (20, '-')):
        d = runs[name][2][fr]
        ax.plot(d['surf'][:, 0], d['surf'][:, 1], ls, color=COLORS[name],
                label='%s, %.0f kyr' % (label.replace('basal friction ', 'φb = '), d['t']))
ax.set_xlim(0, 40)
ax.set_xlabel('x (km)')
ax.set_ylabel('elevation (km)')
ax.set_title('surface topography')
ax.legend(fontsize=8)
for k, (key, ylab, title) in enumerate((('width', 'km', 'wedge width (wall to thrust front)'),
                                        ('hmax', 'km', 'maximum elevation'))):
    ax = axes[k + 1]
    for name, label in MODELS:
        ms = runs[name][1][1:]
        ax.plot([m['t'] * 0.05 for m in ms], [m[key] for m in ms], 'o-', ms=3,
                color=COLORS[name], label=label)
    ax.set_xlabel('shortening (km)')
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.legend(fontsize=8)
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(outdir, 'wedge-growth.png'), dpi=130, bbox_inches='tight')
plt.close(fig)
