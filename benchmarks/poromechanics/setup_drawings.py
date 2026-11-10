#!/usr/bin/env python3
"""Textbook-style setup drawings for every review problem.

The drawings follow the conditions stated in review.SPECS; they explain the
problem and are not meshes or results. The shared review renderer calls this module directly.
Run directly to redraw an existing review in place:

    python3 setup_drawings.py REVIEW
"""
import argparse
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, Rectangle

import diagram_style as ds
import prepare_figures as pf

DRY, WET, SOLID = '#f4f1ea', '#cfe2ee', '#e8edf0'


def fill(ax, x, y, w, h, colors, horizontal=False, z=0):
    """Smooth shading inside a rectangle (saturation, temperature)."""
    ramp = np.linspace(0, 1, 256)
    image = ramp[None, :] if horizontal else ramp[:, None]
    cmap = LinearSegmentedColormap.from_list('fill', colors)
    ax.imshow(image, extent=(x, x+w, y, y+h), origin='lower', cmap=cmap, aspect='auto', zorder=z)


def grains(ax, x, y, w, h, count, seed=1):
    """Faint stipple marking a porous solid."""
    rng = np.random.default_rng(seed)
    p = rng.uniform((x+.04, y+.04), (x+w-.04, y+h-.04), (count, 2))
    ax.scatter(p[:, 0], p[:, 1], s=rng.uniform(2, 6, count), facecolor='white',
               edgecolor='#b9c4cb', lw=.35, zorder=1)


def water_table(ax, x, y, text=None):
    """Hydrogeological water-table mark: an open triangle on the surface."""
    s = .11
    ax.add_patch(Polygon([(x-s, y+1.6*s), (x+s, y+1.6*s), (x, y+.15*s)], fill=False, ec=ds.WATER, lw=ds.BODY, zorder=6))
    for k, half in ((1, .75), (2, .45)):
        ax.plot([x-half*s*1.4, x+half*s*1.4], [y-.55*k*s, y-.55*k*s], color=ds.WATER, lw=ds.THIN+.1, zorder=6)
    if text:
        ds.label(ax, x+.22, y+.2, text, ds.WATER, 10, ha='left')


def side_rollers(ax, x0, x1, y0, y1, n=2, size=.17):
    for yy in np.linspace(y0, y1, n+2)[1:-1]:
        ds.roller(ax, x0, yy, (-1, 0), size); ds.roller(ax, x1, yy, (1, 0), size)


def end_rollers(ax, x0, x1, y, normal, n=2, size=.17):
    for xx in np.linspace(x0, x1, n+2)[1:-1]:
        ds.roller(ax, xx, y, normal, size)


def down_arrows(ax, x0, x1, y, length, color, n=5, bar=False):
    xs = np.linspace(x0+.12, x1-.12, n)
    for xx in xs:
        ds.arrow(ax, (xx, y+length), (xx, y+.02), color, ds.BODY, scale=9)
    if bar:
        ax.plot([xs[0]-.06, xs[-1]+.06], [y+length, y+length], color=color, lw=ds.BODY)


# ----------------------------------------------------------------- columns
COLUMN = dict(
    terzaghi=dict(height=r'$H=10\ \mathrm{m}$', width=r'$W=1\ \mathrm{m}$', body='saturated\nelastic', wet=True),
    gardner=dict(height=r'$H=1\ \mathrm{m}$', width=r'$W=0.05\ \mathrm{m}$', body='Gardner', wet=False),
    van_genuchten=dict(height=r'$H=1\ \mathrm{m}$', width=r'$W=0.05\ \mathrm{m}$', body='van\nGenuchten', wet=False),
    brooks_corey=dict(height=r'$H=1\ \mathrm{m}$', width=r'$W=0.05\ \mathrm{m}$', body='Brooks–\nCorey', wet=False),
    undrained_wave=dict(height=r'$H=1\ \mathrm{m}$', width=r'$W=0.025\ \mathrm{m}$', body='', wet=True))


def column(key, folder):
    c = COLUMN[key]
    with ds.canvas(5.8, 5.4) as (fig, ax):
        x0, y0, w, h = 1.75, 1.0, 1.15, 3.2
        x1, y1 = x0+w, y0+h
        if c['wet']:
            ax.add_patch(Rectangle((x0, y0), w, h, fc=WET, ec='none', zorder=0))
        else:  # saturated at the water table, drier upward
            fill(ax, x0, y0, w, h, [WET, '#e6eef1', DRY])
        grains(ax, x0, y0, w, h, 60)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        if c['body']:
            ds.label(ax, x0+w/2, y0+.62*h, c['body'], ds.MUTED, 9.5, style='italic')
        ds.dimension(ax, (x0-.2, y0), (x0-.2, y1), c['height'], .3, 'left')
        ds.dimension(ax, (x0, y0-.28), (x1, y0-.28), c['width'], .32, 'below')
        ds.label(ax, .2, .15, 'not to scale', ds.MUTED, 7.5, ha='left')
        right = x1+.38
        if key == 'terzaghi':
            down_arrows(ax, x0, x1, y1, .48, ds.LOAD, bar=True)
            ds.label(ax, x0+w/2, y1+.72, r'$\sigma_{zz}=-P$', ds.LOAD)
            ds.leader(ax, (x1, y1), (right, y1+.3), r'drained: $p=0$', ds.WATER, 10)
            ds.leader(ax, (x1+.16, y0+.62*h), (right+.05, y0+.62*h), r'$q_n=0$', ds.WATER, 10)
            ds.leader(ax, (x1-.15, y0), (right+.05, y0+.35), r'sealed: $q_n=0$', ds.WATER, 10)
            ds.arrow(ax, (x1+1.35, y1-.15), (x1+1.35, y1-.85), ds.MUTED, ds.THIN+.2, scale=8)
            ds.label(ax, x1+1.53, y1-.5, '$d$', ds.MUTED, 10.5)
        elif key == 'undrained_wave':
            end_rollers(ax, x0, x1, y1, (0, 1))
            # Initial velocity profile v_z = V sin(pi s), s = (z+H)/H, plotted across the column.
            xc, amp = x0+.2*w, .62*w
            ss = np.linspace(0, 1, 120)
            ax.plot([xc, xc], [y0, y1], color=ds.MUTED, lw=ds.THIN, ls=(0, (2, 2)), zorder=5)
            ax.plot(xc+amp*np.sin(np.pi*ss), y0+ss*h, color=ds.WATER, lw=ds.BODY+.2, zorder=5)
            for s0 in np.linspace(.15, .85, 6):
                xa = xc+amp*np.sin(np.pi*s0)
                ax.plot([xc, xa], [y0+s0*h]*2, color=ds.WATER, lw=ds.THIN, zorder=5)
                ds.arrow(ax, (xa, y0+s0*h-.12), (xa, y0+s0*h+.12), ds.WATER, ds.BODY, scale=8)
            ds.leader(ax, (xc+amp, y0+.5*h), (right+.05, y0+.62*h), r'$v_z=V\sin\pi s$', ds.WATER, 10)
            ds.label(ax, right+.11, y0+.62*h-.34, r'$s=(z+H)/H$', ds.MUTED, 9, ha='left')
            ds.leader(ax, (x1+.16, y0+.3*h), (right+.05, y0+.3*h), r'$q_n=0$ on all faces', ds.WATER, 10)
            ds.label(ax, x0+w/2, y1+.42, r'$u_z=0$')
            ds.label(ax, x1+.42, y0-.12, r'$u_z=0$')
            ds.axes(ax, (4.85, .45))
        else:
            down_arrows(ax, x0, x1, y1, .48, ds.WATER, n=3)
            ds.label(ax, x0+w/2, y1+.72, r'$q_n=-I$', ds.WATER)
            ds.leader(ax, (x1+.16, y0+.7*h), (right+.05, y0+.7*h), r'$q_n=0$', ds.WATER, 10)
            ax.plot([x1, x1+.55], [y0, y0], color=ds.WATER, lw=ds.THIN+.2, zorder=5)
            water_table(ax, x1+.34, y0)
            ds.label(ax, x1+.62, y0+.22, r'$p=0$', ds.WATER, 10, ha='left')
            ds.gravity(ax, 5.2, 4.5)
            ds.axes(ax, (4.85, .45))
            if key == 'brooks_corey':
                # Steady infiltration at I = Ks/2: p = -rho g z/2 meets -1/alpha_p.
                ya = y0+h/(.001*9810*.5)
                ax.plot([x0, x1], [ya, ya], color=ds.HOT, lw=ds.BODY, ls=(0, (4, 3)), zorder=5)
                ds.leader(ax, (x1-.1, ya), (right+.05, ya+.15), 'air entry', ds.HOT, 10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- heat
def heat(folder):
    with ds.canvas(7.2, 3.7) as (fig, ax):
        x0, y0, L, h = 1.35, 1.35, 4.5, .95
        x1, y1 = x0+L, y0+h
        # A thermal front advected to the right (sketch of erfc).
        ramp = np.linspace(0, 1, 400)
        profile = .5*(1-np.tanh((ramp-.42)/.08))
        cmap = LinearSegmentedColormap.from_list('t', ['#f8f6f2', '#f0c3ad'])
        ax.imshow(profile[None, :], extent=(x0, x1, y0, y1), origin='lower', cmap=cmap, aspect='auto', zorder=0)
        grains(ax, x0, y0, L, h, 70)
        ds.wall(ax, x0, y1, L, .13); ds.wall(ax, x0, y0-.13, L, .13)
        ax.add_patch(Rectangle((x0, y0), L, h, fill=False, ec=ds.INK, lw=ds.BODY, zorder=4))
        ax.plot([x0, x0], [y0, y1], color=ds.HOT, lw=3.2, solid_capstyle='butt', zorder=5)
        for yy in (y0+.3*h, y0+.7*h):
            ds.arrow(ax, (x0+2.3, yy), (x0+3.2, yy), ds.WATER, ds.BODY, scale=9)
        ds.arrow(ax, (x1+.08, y0+h/2), (x1+.6, y0+h/2), ds.WATER, ds.BODY, scale=9)
        ds.label(ax, x0+3.75, y0+h/2, r'$q_x>0$', ds.WATER, 10.5)
        ds.label(ax, x0-.1, y1+.4, r'$T=T_{in}$', ds.HOT, 10.5, ha='right')
        ds.label(ax, x0-.1, y0-.4, r'$p=p_0$', ds.WATER, 10.5, ha='right')
        ds.label(ax, x1+.1, y1+.4, r'$\partial T/\partial x=0$', size=10.5, ha='left')
        ds.label(ax, x1+.1, y0-.4, r'$p=p_1$', ds.WATER, 10.5, ha='left')
        ds.label(ax, x0+L/2, y0-.42, r'walls: $q_n=0,\ \partial T/\partial n=0,\ u_n=0$', size=9.5)
        ds.leader(ax, (x0+L*.42, y1-.2), (x0+L*.42+.35, y1+.55), r'front, $v_T$', ds.HOT, 9.5)
        ds.dimension(ax, (x0, y0-.13), (x1, y0-.13), r'$L=1\ \mathrm{m}$', .62, 'below')
        ds.label(ax, x0+L/2, .1, r'initially $T=T_0$;  $g=0$', ds.MUTED, 9)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def injection_temperature(folder):
    """Plan views: radial hot-water injection (Avdonin) and a sealed hot doublet."""
    from matplotlib.patches import Circle, Wedge
    with ds.canvas(8.6, 4.4) as (fig, ax):
        warm = '#f0c3ad'
        # Left: quarter plane, source at the corner, drained far sides.
        x0, y0, s = .9, 1.0, 2.6
        ax.add_patch(Rectangle((x0, y0), s, s, fc=WET, ec='none', zorder=0))
        ax.add_patch(Wedge((x0, y0+s), .85, 270, 360, fc=warm, ec='none', zorder=1))
        ax.add_patch(Rectangle((x0, y0), s, s, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.plot([x0+s, x0+s], [y0, y0+s], color=ds.WATER, lw=2.6, solid_capstyle='butt', zorder=5)
        ax.plot([x0, x0+s], [y0, y0], color=ds.WATER, lw=2.6, solid_capstyle='butt', zorder=5)
        ax.add_patch(Circle((x0, y0+s), .07, fc=ds.HOT, ec=ds.INK, lw=.8, zorder=6))
        for ang in np.radians([-20, -45, -70]):
            ds.arrow(ax, (x0+.95*np.cos(ang), y0+s+.95*np.sin(ang)), (x0+1.45*np.cos(ang), y0+s+1.45*np.sin(ang)), ds.WATER, scale=8)
        ds.label(ax, x0+.1, y0+s+.32, r'$q/4$ at $T_{in}$', ds.HOT, 10, ha='left')
        ds.label(ax, x0+s/2, y0-.32, r'$p=0$, outflow', ds.WATER, 9.5)
        ds.label(ax, x0-.12, y0+s/2, 'symmetry:\nsealed,\ninsulated', ds.MUTED, 8.5, ha='right')
        ds.label(ax, x0+s/2, y0+s+.75, '(a) radial injection, Avdonin', size=10.5)
        ds.dimension(ax, (x0, y0-.55), (x0+s, y0-.55), r'$8\ \mathrm{m}$', .25, 'below')
        # Right: sealed, insulated box with an injector and a producer.
        x1, w, hh = 5.0, 3.2, 2.6
        ax.add_patch(Rectangle((x1, y0), w, hh, fc=WET, ec='none', zorder=0))
        ax.add_patch(Circle((x1+.8, y0+hh/2), .62, fc=warm, ec='none', zorder=1))
        ax.add_patch(Rectangle((x1, y0), w, hh, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.add_patch(Circle((x1+.8, y0+hh/2), .07, fc=ds.HOT, ec=ds.INK, lw=.8, zorder=6))
        ax.add_patch(Circle((x1+2.4, y0+hh/2), .07, fc='white', ec=ds.INK, lw=.8, zorder=6))
        ds.arrow(ax, (x1+1.5, y0+hh/2), (x1+2.2, y0+hh/2), ds.WATER, scale=8)
        ds.label(ax, x1+.8, y0+hh/2-.32, r'$+q$, $T_{in}$', ds.HOT, 9.5)
        ds.label(ax, x1+2.4, y0+hh/2-.32, r'$-q$, $T_{prod}$', size=9.5)
        ds.label(ax, x1+w/2, y0-.32, r'all sides: $q_n=0$, $\partial T/\partial n=0$', size=9.5)
        ds.label(ax, x1+w/2, y0+hh+.75, '(b) sealed hot doublet', size=10.5)
        ds.label(ax, x1+w/2, y0+hh+.3, r'$E-E_0=\rho_fc_fq\int(T_{in}-T_{prod})\,dt$', size=9.5)
        ds.label(ax, x1+w/2, .2, r'initially $T=T_0$, $p=0$; $g=0$; plan views', ds.MUTED, 9)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- point sources
def sources(folder):
    with ds.canvas(5.8, 5.0) as (fig, ax):
        o, s = np.array([1.0, .85]), 2.75
        d = .45*s*np.array([np.cos(np.radians(33)), np.sin(np.radians(33))])
        def P(x, y, z): return o + s*np.array([x, z+1]) + y*d
        corners = {k: P(*k) for k in [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (-1, 0)]}
        def face(keys, color):
            ax.add_patch(Polygon([corners[k] for k in keys], fc=color, ec='none', zorder=0))
        face([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], '#eef3f6')     # top
        face([(1, 0, -1), (1, 1, -1), (1, 1, 0), (1, 0, 0)], '#e3eaee')   # right
        face([(0, 0, -1), (1, 0, -1), (1, 0, 0), (0, 0, 0)], '#f3f6f8')   # front
        hidden = [((0, 1, -1), (1, 1, -1)), ((0, 1, -1), (0, 1, 0)), ((0, 1, -1), (0, 0, -1))]
        for a, b in hidden:
            ax.plot(*np.column_stack((corners[a], corners[b])), color=ds.MUTED, lw=ds.THIN, ls=(0, (3, 2)))
        edges = [((0, 0, -1), (1, 0, -1)), ((1, 0, -1), (1, 0, 0)), ((1, 0, 0), (0, 0, 0)), ((0, 0, 0), (0, 0, -1)),
                 ((1, 0, -1), (1, 1, -1)), ((1, 1, -1), (1, 1, 0)), ((1, 1, 0), (1, 0, 0)),
                 ((0, 0, 0), (0, 1, 0)), ((0, 1, 0), (1, 1, 0))]
        for a, b in edges:
            ax.plot(*np.column_stack((corners[a], corners[b])), color=ds.INK, lw=ds.BODY, zorder=3)
        for n, (x, y, z), sign, color in (('1', (.37, .43, -.53), '+', ds.HOT), ('2', (.61, .57, -.41), '−', ds.WATER)):
            p, floor = P(x, y, z), P(x, y, -1)
            ax.plot(*np.column_stack((p, floor)), color=ds.MUTED, lw=ds.THIN, ls=(0, (2, 2)), zorder=4)
            ax.plot(*floor, 'o', ms=2.5, color=ds.MUTED, zorder=4)
            ax.plot(*p, 'o', ms=13, mfc='white', mec=color, mew=1.3, zorder=6)
            ds.label(ax, *p, sign, color, 11, fontweight='bold', zorder=7)
            ds.label(ax, p[0]-.35 if n == '1' else p[0]+.35, p[1]+.3, n, color, 10.5)
        ds.label(ax, 2.9, 4.65, r'sealed, insulated:  $q_n=0,\ \partial T/\partial n=0,\ u_n=0$', size=9.5)
        ds.label(ax, 2.9, .38, r'1: $+Q,\ +\dot H$ for $0.25\!-\!2.75$ s      2: $-Q,\ -\dot H$ for $0.8\!-\!2.1$ s', size=9)
        ds.dimension(ax, P(0, 0, -1), P(0, 0, 0), r'$1\ \mathrm{m}$', .28, 'left')
        ds.axes(ax, (4.95, 1.05), names=('x', 'z'))
        ds.arrow(ax, (4.95, 1.05), (4.95+.3, 1.05+.2), ds.MUTED, ds.THIN+.2, scale=8)
        ds.label(ax, 5.38, 1.32, '$y$', ds.MUTED, 10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- Theis
def theis(folder):
    """Actual finite quarter-box test, observation band and full well screen."""
    from matplotlib.patches import Wedge
    with ds.canvas(7.0, 4.9) as (fig, ax):
        x0, y0, w = 1.0, 1.3, 2.8
        ax.add_patch(Rectangle((x0, y0), w, w, fc=WET, ec=ds.INK, lw=ds.THICK))
        grains(ax, x0, y0, w, w, 70)
        ax.add_patch(Wedge((x0, y0), w/2, 0, 90, width=w/3,
                          fc='#90bad2', ec='none', alpha=.5, zorder=2))
        theta = np.linspace(0, np.pi/2, 100)
        for radius in (w/6, w/2):
            ax.plot(x0+radius*np.cos(theta), y0+radius*np.sin(theta),
                    color=ds.WATER, lw=ds.THIN, ls='--', zorder=3)
        ax.plot(x0, y0, 'o', color=ds.WATER, ms=6, zorder=6)
        ds.leader(ax, (x0, y0), (.32, .8), r'$Q/4$', ds.WATER, 11)
        ds.leader(ax, (x0+w*.36, y0+w*.36), (2.35, 3.1),
                  r'$10\leq r\leq30\ \mathrm{m}$', ds.WATER, 10, ha='center')
        ds.label(ax, 2.4, 4.65, 'plan view', size=11)
        ds.label(ax, x0+w/2, y0+w+.22, r'$q_n=0$', size=11)
        ds.label(ax, x0+w+.22, y0+w/2, r'$q_n=0$', rotation=90, size=11)
        ds.label(ax, x0+.12, y0+w*.62, 'symmetry', ds.MUTED, 9, rotation=90)
        ds.label(ax, x0+w*.62, y0+.14, 'symmetry', ds.MUTED, 9)
        ds.dimension(ax, (x0, y0), (x0+w, y0), r'$L=60\ \mathrm{m}$', .57, 'below')
        ds.dimension(ax, (x0, y0), (x0, y0+w), '$L$', .48)
        ds.axes(ax, (3.7, .38), length=.35, names=('x', 'y'))
        # Cross-section through the corner well: screen spans the full thickness.
        xl, xr, base, top = 5.0, 6.0, 1.45, 3.75
        ax.add_patch(Rectangle((xl, base), xr-xl, top-base, fc=WET, ec='none'))
        ds.wall(ax, xl, top, xr-xl, .12)
        ds.wall(ax, xl, base-.12, xr-xl, .12)
        for xx in (5.16, 5.34):
            ax.plot([xx, xx], [base, top+.25], color=ds.INK, lw=ds.BODY, zorder=4)
        for yy in np.arange(base+.08, top-.04, .15):
            ax.plot([5.16, 5.34], [yy, yy], color=ds.INK, lw=ds.THIN)
        ds.arrow(ax, (5.25, top+.27), (5.25, top+.7), ds.WATER)
        ds.label(ax, 5.55, top+.5, r'$Q/4$', ds.WATER, 11)
        ds.dimension(ax, (xr, base), (xr, top), r'$B=5\ \mathrm{m}$', .38, 'right')
        ds.label(ax, 5.5, 1.1, 'sealed top and base', size=9)
        ds.label(ax, 5.5, 4.65, 'well screen', size=11)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- averaged aquifer
def aquifer(folder):
    """Cross-section of the actual finite pumping slab and its supports."""
    with ds.canvas(7.4, 4.2) as (fig, ax):
        left, right, base, top = 1.3, 5.6, 1.7, 2.9
        ax.add_patch(Rectangle((left, base), right-left, top-base, fc=WET, ec=ds.INK, lw=ds.BODY))
        grains(ax, left, base, right-left, top-base, 120, seed=6)
        for xx in (2.0, 3.5, 5.0): ds.roller(ax, xx, base, (0, -1), size=.17)
        for yy in (2.05, 2.58): ds.roller(ax, right, yy, (1, 0), size=.17)
        for xx in (left-.2, left):
            ax.plot([xx, xx], [base, top+.24], color=ds.INK, lw=ds.BODY, zorder=5)
        for yy in np.arange(base+.07, top-.03, .14):
            ax.plot([left-.2, left], [yy, yy], color=ds.INK, lw=ds.THIN, zorder=5)
        ds.arrow(ax, (left-.1, top+.27), (left-.1, top+.82), ds.WATER)
        ds.label(ax, left+.21, top+.73, r'$Q/4$', ds.WATER, 11, ha='left')
        for yy in (base+.3, top-.3):
            ds.arrow(ax, (left+1.2, yy), (left+.25, yy), ds.WATER, scale=10)
        ds.label(ax, 3.9, 3.62, r'$q_z=0$', size=10.5)
        ds.label(ax, 3.9, 3.28, r'$\Delta\sigma_{zz}=\Delta\sigma_{xz}=0$', size=10.5)
        ds.arrow(ax, (2.15, 3.5), (2.15, top+.05), ds.INK, scale=9)
        ds.label(ax, 2.4, 3.18, r'$\delta$', size=12)
        ds.label(ax, (left+right)/2, 1.16, r'$q_z=0,\quad U_z=0$', size=10.5)
        ds.label(ax, 6.0, 3.25, r'$p=p_{\rm initial}$', ds.WATER, 10.5)
        ds.label(ax, 6.0, 1.22, r'$U_x=0$', size=10.5)
        ds.dimension(ax, (right+.35, base), (right+.35, top), r'$B=10\ \mathrm{m}$', .5, 'right')
        ds.dimension(ax, (left, base), (right, base), r'$L=400\ \mathrm{m}$', .98, 'below')
        ds.label(ax, left-.12, 1.2, '$x=0$', size=10)
        ds.label(ax, 3.6, .12, 'section through the well at y = 0; not to scale', ds.MUTED, 9, style='italic')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def schematic_convection(folder):
    """Square porous layer heated from below, to scale, with the first-mode cell."""
    with ds.canvas(6.6, 5.0) as (fig, ax):
        x0, y0, w = 1.0, 1.05, 3.0
        shade = pf.matplotlib.colors.LinearSegmentedColormap.from_list('tc', ['#f3d9ce', '#faf7f4', '#dbe8f1'])
        ax.imshow(np.linspace(0, 1, 256)[:, None], extent=(x0, x0+w, y0, y0+w), origin='lower',
                  cmap=shade, aspect='auto', zorder=0)
        # Streamlines of the marginal mode, psi = sin(pi x) sin(pi z).
        gx, gz = np.meshgrid(np.linspace(0, 1, 301), np.linspace(0, 1, 301))
        lines = ax.contour(x0+w*gx, y0+w*gz, np.sin(np.pi*gx)*np.sin(np.pi*gz), levels=[.2, .5, .8],
                           colors=ds.WATER, linewidths=ds.BODY, zorder=2)
        for path in lines.allsegs:
            loop = max(path, key=len)
            left = loop[np.argmin(loop[:, 0])]; right = loop[np.argmax(loop[:, 0])]
            i = np.argmin(loop[:, 0]); up = loop[(i+1) % len(loop)][1] > loop[i-1][1]
            ds.on_curve(ax, loop, left, ds.WATER, reverse=not up)   # warm water rises
            ds.on_curve(ax, loop, right, ds.WATER, reverse=not up)  # cool water sinks
        ds.wall(ax, x0-.14, y0, .14, w); ds.wall(ax, x0+w, y0, .14, w)
        ds.plate(ax, x0-.14, x0+w+.14, y0+w, +1, ds.COLD)
        ds.plate(ax, x0-.14, x0+w+.14, y0, -1, ds.HOT)
        ax.add_patch(Rectangle((x0, y0), w, w, fill=False, ec=ds.INK, lw=ds.BODY, zorder=4))
        ds.label(ax, x0+w/2, y0+w+.45, r'$T=T_0,\quad q_n=0$', ds.COLD)
        ds.label(ax, x0+w/2, y0-.45, r'$T=T_0+\Delta T,\quad q_n=0$', ds.HOT)
        ds.leader(ax, (x0+w+.14, y0+.78*w), (x0+w+.4, y0+.78*w), r'$q_n=0$', size=10.5)
        ds.label(ax, x0+w+.46, y0+.78*w-.36, r'$\partial T/\partial n=0$', size=10.5, ha='left')
        ds.leader(ax, (x0+w*.86, y0+w*.5), (x0+w+.4, y0+.32*w), r'$\psi\propto\sin\pi x\,\sin\pi z$', ds.WATER, 10.5)
        ds.dimension(ax, (x0-.14, y0), (x0-.14, y0+w), r'$H$', .3, 'left')
        ds.dimension(ax, (x0, y0-.14), (x0+w, y0-.14), r'$W=H$', .62, 'below')
        ds.gravity(ax, 6.0, 4.75)
        ds.axes(ax, (5.3, .3))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')

def mandel(folder):
    """Quarter strip under the rigid platen (drawing; not to scale)."""
    with ds.canvas(6.4, 4.9) as (fig, ax):
        x0, y0, w, h = 1.25, 1.0, 3.0, 2.4
        x1, y1 = x0+w, y0+h
        ax.add_patch(Rectangle((x0, y0), w, h, fc=WET, ec='none', zorder=0))
        grains(ax, x0, y0, w, h, 110, seed=8)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.add_patch(Rectangle((x0, y1), w, .16, fc='#5d6b73', ec=ds.INK, lw=ds.BODY, zorder=5))
        down_arrows(ax, x0, x1, y1+.16, .5, ds.LOAD, n=6, bar=True)
        ds.label(ax, x0+w/2, y1+.9, r'total force $F=Pa$ on the half width', ds.LOAD, 10)
        ds.leader(ax, (x1-.3, y1+.1), (x1+.45, y1+.38), 'rigid, frictionless, sealed', size=9.5)
        for yy in np.linspace(y0, y1, 4)[1:-1]:
            ds.roller(ax, x0, yy, (-1, 0), .17)
        end_rollers(ax, x0, x1, y0, (0, -1), n=3)
        for yy in (y0+.25*h, y0+.5*h):
            ds.arrow(ax, (x1-.55, yy), (x1+.1, yy), ds.WATER, ds.BODY, scale=9)
        ds.leader(ax, (x1, y0+.75*h), (x1+.45, y0+.75*h), r'drained: $p=0$', ds.WATER, 10)
        ds.label(ax, x1+.51, y0+.75*h-.36, 'traction-free', size=9.5, ha='left')
        ds.label(ax, x0-.62, y0+.5*h, 'symmetry', ds.MUTED, 8.5, rotation=90, style='italic')
        ds.label(ax, x0+w/2, y0-.42, 'symmetry (mid-plane)', ds.MUTED, 8.5, style='italic')
        ds.dimension(ax, (x0, y0-.55), (x1, y0-.55), '$a$', .25, 'below')
        ds.dimension(ax, (x1+1.6, y0), (x1+1.6, y1), '$b$', .2, 'right')
        ds.axes(ax, (.3, .35))
        assert not ds.check_labels(fig, ax, strict=False)
        pf.save_figure(fig, folder, 'setup')


def cryer(folder):
    """Section of one octant: quarter disk with symmetry axes and the loaded arc."""
    with ds.canvas(5.8, 5.0) as (fig, ax):
        o, R = np.array([1.1, 4.3]), 3.1
        t = np.linspace(0, -np.pi/2, 200)
        arc = o + R*np.column_stack((np.cos(t), np.sin(t)))
        ax.add_patch(Polygon(np.vstack((o, arc)), fc=WET, ec='none', zorder=0))
        rng = np.random.default_rng(9)
        pts = rng.uniform((o[0], o[1]-R), (o[0]+R, o[1]), (220, 2))
        pts = pts[np.hypot(*(pts-o).T) < R-.06]
        ax.scatter(pts[:, 0], pts[:, 1], s=rng.uniform(2, 6, len(pts)), facecolor='white', edgecolor='#b9c4cb', lw=.35, zorder=1)
        ax.plot(arc[:, 0], arc[:, 1], color=ds.INK, lw=ds.THICK, zorder=4)
        ax.plot([o[0], o[0]+R], [o[1], o[1]], color=ds.INK, lw=ds.BODY, zorder=4)
        ax.plot([o[0], o[0]], [o[1], o[1]-R], color=ds.INK, lw=ds.BODY, zorder=4)
        for ang in np.linspace(-.22, -1.35, 6):
            n = np.array([np.cos(ang), np.sin(ang)])
            ds.arrow(ax, o+(R+.55)*n, o+(R+.03)*n, ds.LOAD, ds.BODY, scale=9)
        for xx in (o[0]+.9, o[0]+2.0):
            ds.roller(ax, xx, o[1], (0, 1), .16)
        for yy in (o[1]-.9, o[1]-2.0):
            ds.roller(ax, o[0], yy, (-1, 0), .16)
        ds.label(ax, 4.55, .55, r'$\sigma_n=-P,\ p=0$', ds.LOAD, 10.5)
        ds.dimension(ax, o, o+R*np.array([np.cos(-.95), np.sin(-.95)]), '$R$', 0, 'left')
        ds.label(ax, o[0]+.12, o[1]+.36, 'centre', ds.MUTED, 8.5, ha='left')
        ds.label(ax, o[0]+R/2+.25, o[1]+.4, r'symmetry: $u_n=0,\ q_n=0$', ds.MUTED, 9)
        ds.label(ax, .45, o[1]-R/2, 'symmetry', ds.MUTED, 8.5, rotation=90, style='italic')
        assert not ds.check_labels(fig, ax, strict=False)
        pf.save_figure(fig, folder, 'setup')


def rudnicki(folder):
    """Section through one octant: point source at the corner, diffusion front, far faces."""
    with ds.canvas(5.8, 5.0) as (fig, ax):
        o, L = np.array([1.0, 4.4]), 3.6
        ax.add_patch(Rectangle((o[0], o[1]-L), L, L, fc=SOLID, ec='none', zorder=0))
        t = np.linspace(0, -np.pi/2, 100)
        for radius, alpha in ((2.2, .55), (1.4, .8)):
            ax.add_patch(Polygon(np.vstack((o, o+radius*np.column_stack((np.cos(t), np.sin(t))))), fc=WET, ec='none', alpha=alpha, zorder=0))
        grains(ax, o[0], o[1]-L, L, L, 150, seed=12)
        front = o+2.2*np.column_stack((np.cos(t), np.sin(t)))
        ax.plot(front[:, 0], front[:, 1], color=ds.WATER, lw=ds.THIN+.2, ls=(0, (4, 3)), zorder=3)
        ax.add_patch(Rectangle((o[0], o[1]-L), L, L, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.plot(*o, 'o', ms=9, mfc='white', mec=ds.WATER, mew=1.4, zorder=6)
        ds.label(ax, *o, '+', ds.WATER, 11, fontweight='bold', zorder=7)
        for ang in np.linspace(-.25, -1.32, 4):
            n = np.array([np.cos(ang), np.sin(ang)])
            ds.arrow(ax, o+.35*n, o+1.0*n, ds.WATER, ds.BODY, scale=9)
        for xx in (o[0]+1.2, o[0]+2.6):
            ds.roller(ax, xx, o[1], (0, 1), .16)
        for yy in (o[1]-1.2, o[1]-2.6):
            ds.roller(ax, o[0], yy, (-1, 0), .16)
        ds.label(ax, o[0]+.15, o[1]+.55, r'$Q/8$ at the corner', ds.WATER, 10, ha='left')
        ds.leader(ax, o+2.2*np.array([np.cos(-.5), np.sin(-.5)]), (o[0]+L+.12, o[1]-.6), r'front $r\sim2\sqrt{ct}$', ds.WATER, 9.5)
        ds.label(ax, o[0]+L+.12, o[1]-L+.55, 'far faces:', size=9.5, ha='left')
        ds.label(ax, o[0]+L+.12, o[1]-L+.2, r'$p=0,\ u_n=0$', size=9.5, ha='left')
        ds.label(ax, .5, o[1]-L/2, r'symmetry: $u_n=0,\ q_n=0$', ds.MUTED, 8.5, rotation=90)
        ds.dimension(ax, (o[0], o[1]-L), (o[0]+L, o[1]-L), r'$L=60\ \mathrm{m}$ (not to scale)', .3, 'below')
        assert not ds.check_labels(fig, ax, strict=False)
        pf.save_figure(fig, folder, 'setup')


def _opening(folder, plastic=False):
    """Shared quarter-annulus geometry for the elastic/plastic tunnel ladder."""
    from matplotlib.patches import Wedge
    with ds.canvas(7.0, 4.8) as (fig, ax):
        centre, outer, inner = np.array([1.3, 3.8]), 2.7, .55
        ax.add_patch(Wedge(centre, outer, -90, 0, width=outer-inner,
                          fc=SOLID, ec=ds.INK, lw=ds.THICK))
        if plastic:
            boundary = inner*1.396
            ax.add_patch(Wedge(centre, boundary, -90, 0, width=boundary-inner,
                              fc='#e7cbb5', ec=ds.INK, lw=ds.THIN))
            ds.leader(ax, centre+boundary*np.array([np.cos(-.8), np.sin(-.8)]),
                      (4.6, 3.05), r'plastic zone: $b/a=1.396$', size=9, ha='left')
        for xx in (2.3, 3.4): ds.roller(ax, xx, centre[1], (0, 1), .15)
        for yy in (1.6, 2.7): ds.roller(ax, centre[0], yy, (-1, 0), .15)
        for angle in np.linspace(-1.4, -.15, 5):
            direction = np.array([np.cos(angle), np.sin(angle)])
            ds.arrow(ax, centre+(outer+.35)*direction, centre+outer*direction, ds.LOAD, scale=9)
        angle = -.7
        direction = np.array([np.cos(angle), np.sin(angle)])
        ds.arrow(ax, centre, centre+outer*direction, ds.MUTED, ds.THIN, two=True, scale=7)
        ds.label(ax, 2.7, 2.8, '$R=20a$', size=10, rotation=-40)
        ds.leader(ax, centre+inner*np.array([np.cos(-.25), np.sin(-.25)]), (4.6, 3.75),
                  r'free wall: $\sigma_{rr}=\sigma_{r\theta}=0$', size=9, ha='left')
        ds.label(ax, 4.7, 2.6, r'$\sigma_{rr}=-P$', ds.LOAD, 11, ha='left')
        ds.label(ax, 4.7, 2.1, '$P=500$ kPa' if plastic else '$P=100$ kPa', ds.LOAD, 10, ha='left')
        if plastic:
            ds.label(ax, 4.7, 1.65, r'$P_i:P\longrightarrow0$', ds.LOAD, 10, ha='left')
            ds.label(ax, 4.7, 1.15, r'$\phi=\psi=30^\circ$', size=10, ha='left')
        ds.label(ax, 2.8, 4.45, r'symmetry: $u_z=0$', size=10)
        ds.label(ax, .35, 2.4, r'$u_x=0$', size=10, rotation=90)
        ds.label(ax, 3.7, .4, r'$a=1$ m; plane strain; not to scale', size=9)
        ds.axes(ax, (.3, .35), length=.4)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def kirsch(folder):
    """Elastic quarter annulus with a free opening."""
    _opening(folder)


def salencon(folder):
    """Plastic tunnel after removing the inner support pressure."""
    _opening(folder, plastic=True)


def borehole(folder):
    """Quarter annulus with fluid injected at the wall; pressure diffuses outward."""
    from matplotlib.patches import Wedge
    with ds.canvas(7.0, 4.8) as (fig, ax):
        centre, outer, inner = np.array([1.3, 3.8]), 2.7, .55
        ax.add_patch(Wedge(centre, outer, -90, 0, width=outer-inner, fc=SOLID, ec='none'))
        # Pressurized zone near the wall, fading outward.
        for k, frac in enumerate(np.linspace(1, 0, 6)[:-1]):
            radius = inner + .9*(1 - frac) + .18
            ax.add_patch(Wedge(centre, radius, -90, 0, width=radius-inner, fc=WET, ec='none', alpha=.35, zorder=1))
        ax.add_patch(Wedge(centre, outer, -90, 0, width=outer-inner, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.add_patch(Wedge(centre, inner, -90, 0, fc='#9ec3da', ec=ds.WATER, lw=ds.BODY, zorder=3))
        for xx in (2.3, 3.4): ds.roller(ax, xx, centre[1], (0, 1), .15)
        for yy in (1.6, 2.7): ds.roller(ax, centre[0], yy, (-1, 0), .15)
        for angle in np.linspace(-1.35, -.2, 4):
            direction = np.array([np.cos(angle), np.sin(angle)])
            ds.arrow(ax, centre+(inner+.08)*direction, centre+(inner+.75)*direction, ds.WATER, scale=8)
        angle = -.7
        direction = np.array([np.cos(angle), np.sin(angle)])
        ds.arrow(ax, centre+1.3*direction, centre+outer*direction, ds.MUTED, ds.THIN, two=True, scale=7)
        ds.label(ax, 2.95, 2.55, '$R=20a$', size=10, rotation=-40)
        ds.label(ax, 4.6, 3.75, r'wall $r=a$, from $t=0^+$:', ds.WATER, 9.5, ha='left')
        ds.label(ax, 4.6, 3.3, r'$p=p_w$, $\sigma_{rr}=-p_w$', ds.WATER, 10, ha='left')
        ds.leader(ax, centre+outer*np.array([np.cos(-.3), np.sin(-.3)]), (4.6, 2.35),
                  r'outer: sealed, traction-free', size=9.5, ha='left')
        ds.label(ax, 2.8, 4.45, r'symmetry: $u_z=0$, $q_n=0$', size=10)
        ds.label(ax, .35, 2.4, r'$u_x=0$, $q_n=0$', size=10, rotation=90)
        ds.label(ax, 3.7, .4, r'$a=1$ m; plane strain; initially $p=0$, unstressed; not to scale', size=9)
        ds.axes(ax, (.3, .35), length=.4)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def theis_wells(folder):
    """Full confined slab with two pumps and one injecting line well."""
    with ds.canvas(7.0, 4.8) as (fig, ax):
        x0, y0, w = 1.0, 1.15, 2.8
        ax.add_patch(Rectangle((x0, y0), w, w, fc=WET, ec=ds.INK, lw=ds.THICK))
        grains(ax, x0, y0, w, w, 75, seed=20)
        for x, y, name, sign in ((60, 80, '1', 1), (100, 80, '2', 1), (80, 100, '3', -1)):
            centre = np.array([x0+w*x/160, y0+w*y/160])
            ax.plot(*centre, 'o', color=ds.WATER, ms=5, zorder=6)
            for angle in np.linspace(0, 2*np.pi, 4, endpoint=False):
                direction = np.array([np.cos(angle), np.sin(angle)])
                a, b = centre+.3*direction, centre+.08*direction
                ds.arrow(ax, a if sign > 0 else b, b if sign > 0 else a, ds.WATER, scale=7)
            ds.label(ax, centre[0], centre[1]-.42 if sign > 0 else centre[1]+.42, f'well {name}', size=9)
        ds.label(ax, 2.4, 4.55, 'plan view', size=11)
        ds.label(ax, 2.4, 4.15, r'all faces sealed: $q_n=0$', size=10)
        for y, text in ((3.5, r'$Q_1=Q$'), (3.0, r'$Q_2=0.6Q$'), (2.5, r'$Q_3=-0.25Q$')):
            ds.label(ax, 4.55, y, text, ds.WATER, 11, ha='left')
        ds.label(ax, 4.55, 1.75, r'$B=5\ \mathrm{m}$', size=10, ha='left')
        ds.label(ax, 4.55, 1.25, 'pump/inject: 0–40 s\nrecovery: 40–80 s', size=9, ha='left')
        ds.dimension(ax, (x0, y0), (x0+w, y0), r'$L=160\ \mathrm{m}$', .5, 'below')
        ds.axes(ax, (3.8, .25), length=.35, names=('x', 'y'))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def open_flow(folder):
    """A pressure eigenmode draining through two prescribed-pressure ends."""
    with ds.canvas(7.0, 3.8) as (fig, ax):
        left, right, base, top = 1.25, 5.75, 1.35, 2.35
        ax.add_patch(Rectangle((left, base), right-left, top-base, fc=WET, ec=ds.INK, lw=ds.THICK))
        grains(ax, left, base, right-left, top-base, 70, seed=21)
        ds.arrow(ax, (left+.7, 1.85), (left-.5, 1.85), ds.WATER, scale=10)
        ds.arrow(ax, (right-.7, 1.85), (right+.5, 1.85), ds.WATER, scale=10)
        ds.label(ax, left-.45, 2.5, '$p=0$', ds.WATER, 11)
        ds.label(ax, right+.45, 2.5, '$p=0$', ds.WATER, 11)
        ds.label(ax, 3.5, 3.35, r'$p(x,0)=P\sin(\pi x/L)$', size=12)
        ds.label(ax, 3.5, 2.75, r'sealed sides: $q_n=0$', size=10)
        ds.label(ax, 3.5, .92, r'cross-section $A=1\ \mathrm{m^2}$; transport only', size=9)
        ds.dimension(ax, (left, base), (right, base), r'$L=10\ \mathrm{m}$', .95, 'below')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- Bishop column
def bishop(folder):
    """Confined column under prescribed water tables at 9 m and 5 m."""
    with ds.canvas(6.4, 5.6) as (fig, ax):
        x0, y0, w, h = 1.8, 1.0, 1.15, 3.4
        x1, y1 = x0+w, y0+h
        z = lambda metres: y0 + h*metres/10
        ax.add_patch(Rectangle((x0, y0), w, z(5)-y0, fc=WET, ec='none', zorder=0))
        fill(ax, x0, z(5), w, y1-z(5), [WET, '#e6eef1', DRY])
        grains(ax, x0, y0, w, h, 60)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        for metres, text, dash in ((9, r'$z_0=9\ \mathrm{m}$, start', (0, (4, 3))), (5, r'$z_1=5\ \mathrm{m}$, end', None)):
            ax.plot([x0, x1], [z(metres)]*2, color=ds.WATER, lw=ds.BODY, ls=dash or '-', zorder=5)
            water_table(ax, x0+.3, z(metres))
            ds.leader(ax, (x1, z(metres)), (x1+.42, z(metres)+.12), text, ds.WATER, 10)
        ax.plot([x0, x1], [z(7)]*2, color=ds.MUTED, lw=ds.THIN, ls=(0, (1, 2)), zorder=5)
        ds.leader(ax, (x1, z(7)), (x1+.42, z(7)+.05), 'two-law case: interface', ds.MUTED, 9)
        ds.label(ax, x0+w/2, y1+.3, r'$\sigma_{zz}=0,\ q_n=0$')
        ds.label(ax, x0+w/2, y0-.62, r'$p=\rho_f g\,z_w$, $z_w$: 9 m $\to$ 5 m', ds.WATER, 10)
        ds.leader(ax, (x1, z(2.5)), (x1+.42, z(2.5)), r'$u_x=0,\ q_n=0$', size=10)
        ds.dimension(ax, (x0-.2, y0), (x0-.2, y1), r'$H=10\ \mathrm{m}$', .3, 'left')
        ds.gravity(ax, 5.6, 4.6)
        ds.axes(ax, (5.25, .5))
        ds.label(ax, .2, .15, 'not to scale', ds.MUTED, 7.5, ha='left')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- leaky aquifer
def hantush(folder):
    """Section through the corner well: aquifer, one-element aquitard, held top."""
    with ds.canvas(7.0, 4.6) as (fig, ax):
        x0, x1, base, top, cap = 1.2, 6.2, 1.2, 2.7, 3.3
        ax.add_patch(Rectangle((x0, base), x1-x0, top-base, fc=WET, ec='none', zorder=0))
        grains(ax, x0, base, x1-x0, top-base, 80)
        ax.add_patch(Rectangle((x0, top), x1-x0, cap-top, fc='#e3ddcf', ec='none', zorder=0))
        ax.add_patch(Rectangle((x0, base), x1-x0, cap-base, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.plot([x0, x1], [top, top], color=ds.INK, lw=ds.THIN, zorder=4)
        ds.wall(ax, x0, base-.12, x1-x0, .12)
        # Well screen in the aquifer only.
        for xx in (x0+.06, x0+.24):
            ax.plot([xx, xx], [base, cap+.25], color=ds.INK, lw=ds.BODY, zorder=5)
        for yy in np.arange(base+.08, top-.04, .15):
            ax.plot([x0+.06, x0+.24], [yy, yy], color=ds.INK, lw=ds.THIN, zorder=5)
        ds.arrow(ax, (x0+.15, cap+.27), (x0+.15, cap+.7), ds.WATER)
        ds.label(ax, x0+.42, cap+.5, r'$Q/4$', ds.WATER, 11, ha='left')
        # Leakage through the aquitard, strongest near the well.
        for xx, length in ((1.9, .5), (2.7, .42), (3.6, .34), (4.6, .26)):
            ds.arrow(ax, (xx, cap-.05), (xx, cap-.05-length), ds.WATER, ds.BODY, scale=8)
        ds.label(ax, (x0+x1)/2+.6, cap+.25, r'$p=0$ held above the aquitard', ds.WATER, 10)
        ds.label(ax, 5.3, (top+cap)/2, "aquitard $k_z'$", ds.MUTED, 9.5, style='italic')
        ds.label(ax, 5.3, (base+top)/2, 'aquifer $k$', ds.MUTED, 9.5, style='italic')
        ds.dimension(ax, (x1, base), (x1, top), r'$B=2.5\ \mathrm{m}$', .42, 'right')
        ds.dimension(ax, (x1, top), (x1, cap), "$b'=h$", .42, 'right')
        ds.label(ax, (x0+x1)/2, base-.42, r'sealed base; far side $r=75$ m sealed', size=9.5)
        ax.plot([2.2, 4.6], [base+.35]*2, color=ds.WATER, lw=ds.THIN, ls='--', zorder=5)
        ds.label(ax, 3.4, base+.55, r'compared: $10\leq r\leq30$ m', ds.WATER, 9.5)
        ds.label(ax, 3.6, 4.35, r"leakage factor $\lambda^2=kB\,b'/k_z'=(15\ \mathrm{m})^2$", size=10.5)
        ds.axes(ax, (.35, .45), names=('r', 'z'))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- McTigue heating
def mctigue(folder):
    """Saturated column whose drained, traction-free top is heated at t = 0+."""
    with ds.canvas(6.4, 5.6) as (fig, ax):
        x0, y0, w, h = 1.8, 1.0, 1.15, 3.4
        x1, y1 = x0+w, y0+h
        ax.add_patch(Rectangle((x0, y0), w, h, fc=WET, ec='none', zorder=0))
        fill(ax, x0, y0+.55*h, w, .45*h, [WET, '#f1dcd2'])
        grains(ax, x0, y0, w, h, 60)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ax.plot([x0, x1], [y1, y1], color=ds.HOT, lw=ds.THICK+.6, zorder=5)
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        ds.label(ax, x0+w/2, y1+.62, r'$T=T_0+\Delta T$ for $t>0$', ds.HOT)
        ds.label(ax, x0+w/2, y1+.27, r'$p=0,\ \sigma_{zz}=0$', ds.WATER, 10.5)
        ds.leader(ax, (x1, y0+.45*h), (x1+.42, y0+.45*h), r'sealed, insulated', size=10)
        ds.label(ax, x1+.48, y0+.45*h-.3, r'$u_x=0$', ds.MUTED, 9.5, ha='left')
        ds.label(ax, x0+w/2, y0-.62, r'sealed, insulated, $u_z=0$', size=10)
        ds.dimension(ax, (x0-.2, y0), (x0-.2, y1), r'$L=1\ \mathrm{m}$', .3, 'left')
        ds.arrow(ax, (x1+1.2, y1-.1), (x1+1.2, y1-.95), ds.MUTED, ds.THIN+.2, scale=8)
        ds.label(ax, x1+1.4, y1-.52, 'depth $z$', ds.MUTED, 10, ha='left')
        ds.label(ax, .2, .15, 'not to scale', ds.MUTED, 7.5, ha='left')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- Bishop cell
def bishop_cell(folder):
    """Sealed, partly saturated cell compressed at a prescribed rate."""
    with ds.canvas(6.0, 5.0) as (fig, ax):
        x0, y0, w = 1.6, 1.0, 2.6
        x1, y1 = x0+w, y0+w
        fill(ax, x0, y0, w, w, ['#dcebf3', '#eef1ea'])
        grains(ax, x0, y0, w, w, 70)
        ax.add_patch(Rectangle((x0, y0), w, w, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        ds.plate(ax, x0, x1, y1, +1, ds.LOAD)
        down_arrows(ax, x0+.3, x1-.3, y1+.14, .5, ds.LOAD, n=4)
        ds.label(ax, x0+w/2, y1+.92, r'$\varepsilon_{zz}\to-10^{-3}$ at a constant rate over 1000 s', ds.LOAD, 10.5)
        ds.label(ax, x0+w/2, y0+.62*w, r'$p_0<0$, Gardner $S_e(p)$', ds.MUTED, 10)
        ds.label(ax, x0+w/2, y0+.38*w, r'$k\approx0$: no flow', ds.MUTED, 10)
        ds.leader(ax, (x1, y0+.5*w), (x1+.45, y0+.5*w), r'$q_n=0$ on all faces', ds.WATER, 10)
        ds.label(ax, x1+.51, y0+.5*w-.32, r'$u_n=0$ on sides and base', ds.MUTED, 9.5, ha='left')
        ds.dimension(ax, (x0-.25, y0), (x0-.25, y1), r'$1\ \mathrm{m}$', .3, 'left')
        ds.axes(ax, (5.2, .45))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- coupled drainage (U9)
def drainage(folder):
    """Saturated column opened to drainage at its base; the top desaturates."""
    with ds.canvas(6.4, 5.6) as (fig, ax):
        x0, y0, w, h = 2.0, 1.0, 1.0, 3.4
        x1, y1 = x0+w, y0+h
        fill(ax, x0, y0, w, h, [WET, WET, '#e6eef1', DRY])
        grains(ax, x0, y0, w, h, 60)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        ax.plot([x0, x1], [y1]*2, color=ds.WATER, lw=ds.BODY, ls=(0, (4, 3)), zorder=5)
        water_table(ax, x0+.3, y1)
        ds.leader(ax, (x1, y1), (x1+.42, y1+.2), r'water table, $t=0$', ds.WATER, 10)
        ds.label(ax, x0+w/2, y1+.4, r'$\sigma_{zz}=0,\ q_n=0$')
        ds.label(ax, x0+w/2, y1-.28*h, 'desaturates', ds.MUTED, 9.5)
        ax.plot([x0, x1], [y0, y0], color=ds.WATER, lw=2.6, zorder=5, solid_capstyle='butt')
        ds.label(ax, x0+w/2, y0-.62, r'base: $p=0$ for $t>0$, drains', ds.WATER, 10)
        ds.leader(ax, (x1, y0+.45*h), (x1+.42, y0+.45*h), r'$u_x=0,\ q_n=0$', size=10)
        ds.label(ax, x1+.48, y0+.45*h-.4, r'van Genuchten $S_e(p)$, $\chi=S_e$', ds.MUTED, 9.5, ha='left')
        ds.dimension(ax, (x0-.2, y0), (x0-.2, y1), r'$H=1\ \mathrm{m}$', .3, 'left')
        ds.gravity(ax, 5.7, 3.6)
        ds.axes(ax, (5.25, .5))
        ds.label(ax, .2, .15, 'not to scale', ds.MUTED, 7.5, ha='left')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- pressure-driven yield
def pressure_yield(folder):
    """Confined layer under a surcharge whose pore pressure rises; a yield front moves down."""
    with ds.canvas(6.6, 5.4) as (fig, ax):
        x0, y0, w, h = 1.8, 1.0, 1.4, 3.2
        x1, y1 = x0+w, y0+h
        front = y1 - .6*h
        ax.add_patch(Rectangle((x0, y0), w, h, fc=WET, ec='none', zorder=0))
        ax.add_patch(Rectangle((x0, front), w, y1-front, fc='#f0d9cf', ec='none', zorder=0))
        grains(ax, x0, y0, w, h, 70)
        ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        side_rollers(ax, x0, x1, y0, y1)
        end_rollers(ax, x0, x1, y0, (0, -1))
        down_arrows(ax, x0, x1, y1+.04, .5, ds.LOAD, n=5, bar=True)
        ds.label(ax, x0+w/2, y1+.78, r'surcharge $q=20$ kPa, drained', ds.LOAD, 10.5)
        ax.plot([x0, x1], [front]*2, color=ds.HOT, lw=ds.BODY, ls=(0, (4, 3)), zorder=5)
        ds.arrow(ax, (x0+w/2, front+.05), (x0+w/2, front-.35), ds.HOT)
        ds.leader(ax, (x1, front), (x1+.42, front+.1), r'yield front, $\Delta p=\Delta p^*(d)$', ds.HOT, 10)
        ds.label(ax, x0+w/2, y1-.3*h+.15, 'yielded', ds.HOT, 10)
        ds.label(ax, x0+w/2, y0+.22*h, 'elastic', ds.MUTED, 10)
        ds.leader(ax, (x1, y1-.04), (x1+.42, y1+.3), r'$p=\rho_f g\,(z_w-z)+\Delta p(t)$, ramped', ds.WATER, 10)
        ds.leader(ax, (x1, y0+.18*h), (x1+.42, y0+.18*h), r'$u_n=0,\ q_n=0$', size=10)
        ds.label(ax, x0+w/2, y0-.6, r'Mohr–Coulomb: $c=10$ kPa, $\phi=20^\circ$', ds.MUTED, 10)
        ds.dimension(ax, (x0-.2, y0), (x0-.2, y1), r'$H=2\ \mathrm{m}$', .3, 'left')
        ds.gravity(ax, 5.9, 4.4)
        ds.axes(ax, (5.55, .45))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


# ----------------------------------------------------------------- regional subsidence
def regional(folder):
    """Section through the basin: well field, slab aquifer, bowl and closed outer wall; schedule inset."""
    with ds.canvas(7.4, 4.6) as (fig, ax):
        x0, x1, base, top = .9, 6.4, 1.3, 2.3
        ax.add_patch(Rectangle((x0, base), x1-x0, top-base, fc=WET, ec='none', zorder=0))
        grains(ax, x0, base, x1-x0, top-base, 90)
        ax.add_patch(Rectangle((x0, base), x1-x0, top-base, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        ds.wall(ax, x0, base-.12, x1-x0, .12)
        ds.wall(ax, x1, base, .12, top-base)
        # Exaggerated subsidence bowl of the free top.
        xs = np.linspace(x0, x1, 200)
        ax.plot(xs, top - .45*np.exp(-((xs - x0)/1.4)**2) - .05*np.exp(-((xs - x0)/3.5)**2), color=ds.LOAD,
                lw=ds.BODY, ls=(0, (4, 2)), zorder=5)
        ds.leader(ax, (x0+1.2, top-.25), (x0+2.0, top+.55), 'subsidence (exaggerated)', ds.LOAD, 9.5)
        for xx in (x0+.05, x0+.2):
            ax.plot([xx, xx], [base, top+.3], color=ds.INK, lw=ds.BODY, zorder=6)
        ds.arrow(ax, (x0+.125, top+.32), (x0+.125, top+.75), ds.WATER)
        ds.label(ax, x0+.33, top+.6, r'$Q(t)/4$', ds.WATER, 10.5, ha='left')
        ds.label(ax, x0-.12, (base+top)/2, 'symmetry', ds.MUTED, 8.5, rotation=90)
        ds.label(ax, x1+.42, (base+top)/2, r'$p=0$, $u_n=0$', size=9.5, rotation=90)
        ds.dimension(ax, (x1+.12, base), (x1+.12, top), r'$B=50\ \mathrm{m}$', .62, 'right')
        ds.dimension(ax, (x0, base-.12), (x1, base-.12), r'$L=3\ \mathrm{km}$', .38, 'below')
        ds.label(ax, (x0+x1)/2+.6, (base+top)/2, r'soft aquifer: $K=5$, $G=3$ MPa, $k=10^{-12}\ \mathrm{m^2}$',
                 ds.MUTED, 9)
        # Schedule inset.
        ix, iy, iw, ih = 4.0, 3.1, 2.3, .8
        ax.plot([ix, ix+iw], [iy, iy], color=ds.MUTED, lw=ds.THIN)
        steps_x = [ix, ix, ix+iw/4, ix+iw/4, ix+iw/2, ix+iw/2, ix+3*iw/4, ix+3*iw/4, ix+iw]
        steps_y = [iy, iy+ih, iy+ih, iy, iy, iy+ih, iy+ih, iy, iy]
        ax.plot(steps_x, steps_y, color=ds.WATER, lw=ds.BODY)
        ds.label(ax, ix+iw/2, iy+ih+.22, 'pump 6 months, rest 6 months, twice', size=9)
        ds.label(ax, ix+iw+.05, iy-.05, '2 yr', ds.MUTED, 8.5, ha='left')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def strength(folder):
    """Affine, zero-volume material path at fixed pore pressure."""
    with ds.canvas(7.2, 4.8) as (fig, ax):
        ax.add_patch(Rectangle((2.1, 1.25), 2.3, 2.3, facecolor='#edf4f7', edgecolor=ds.INK, lw=ds.THICK))
        ds.label(ax, 3.25, 2.65, r'$\boldsymbol{\sigma}_0=-C\mathbf{I}$', size=12)
        ds.label(ax, 3.25, 2.15, r'$p=\mathrm{constant}$', color=ds.WATER, size=12)
        ds.label(ax, 3.25, 1.7, r'$\dot\epsilon_v=0$', size=12)
        for x in (2.65, 3.85): ds.arrow(ax, (x, 4.0), (x, 3.6), color=ds.LOAD)
        for y in (1.85, 3.0): ds.arrow(ax, (4.45, y), (4.95, y), color=ds.LOAD)
        ds.label(ax, 3.25, 4.35, r'$v_z(0)=-\dot e$', color=ds.LOAD)
        ds.label(ax, 5.9, 2.5, r'$v_x(1)=\dot e$', color=ds.LOAD)
        ds.label(ax, .95, 2.4, r'$v_x(0)=0$', size=11)
        ds.label(ax, 3.25, .9, r'$v_z(-1)=0$', size=11)
        ds.arrow(ax, (5.35, 1.25), (5.95, 1.25), scale=8)
        ds.arrow(ax, (5.35, 1.25), (5.35, 1.85), scale=8)
        ds.label(ax, 6.1, 1.25, '$x$', size=10)
        ds.label(ax, 5.35, 2.0, '$z$', size=10)
        ds.label(ax, 3.6, .3, 'Unit cell; plane strain in 2D, zero y strain in 3D', size=10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')



def rsf(folder):
    """The strength cell loaded at a rate that steps up tenfold; friction follows rate and state."""
    with ds.canvas(7.6, 4.8) as (fig, ax):
        ax.add_patch(Rectangle((1.6, 1.25), 2.3, 2.3, facecolor='#edf4f7', edgecolor=ds.INK, lw=ds.THICK))
        ds.label(ax, 2.75, 2.75, r'$\boldsymbol{\sigma}_0=-C\mathbf{I}$', size=11.5)
        ds.label(ax, 2.75, 2.3, r'$p=0.8$ MPa, fixed', color=ds.WATER, size=11)
        ds.label(ax, 2.75, 1.85, r'$\mu(V,\theta)$', size=11.5)
        for x in (2.15, 3.35): ds.arrow(ax, (x, 4.0), (x, 3.6), color=ds.LOAD)
        for y in (1.85, 3.0): ds.arrow(ax, (3.95, y), (4.45, y), color=ds.LOAD)
        ds.label(ax, 2.75, 4.35, r'$v_z=-\dot e(t)$', color=ds.LOAD)
        ds.label(ax, 2.75, .9, r'$v_x=\dot e(t)$ on the right; rollers left and below', size=10)
        # Inset: the imposed rate steps up tenfold at t_1.
        x0, y0, w, h = 5.0, 2.2, 2.0, 1.2
        ax.plot([x0, x0+w], [y0, y0], color=ds.INK, lw=ds.THIN); ax.plot([x0, x0], [y0, y0+h], color=ds.INK, lw=ds.THIN)
        ax.plot([x0, x0+w/2, x0+w/2, x0+w], [y0+.15, y0+.15, y0+1.0, y0+1.0], color=ds.LOAD, lw=ds.BODY)
        ds.label(ax, x0+w/2, y0-.22, '$t_1$', size=10)
        ds.label(ax, x0-.15, y0+h+.12, r'$\dot e$', size=10, ha='left')
        ds.label(ax, x0+w+.1, y0+1.0, r'$10\,\dot e_1$', color=ds.LOAD, size=10, ha='left')
        ds.label(ax, x0+.45, y0+.35, r'$\dot e_1$', color=ds.LOAD, size=10)
        ds.label(ax, x0+w/2, 1.45, r'$\dot\theta=1-V\theta/D_c$', size=10.5)
        ds.label(ax, x0+w/2, 1.0, r'$V=2w\,\dot\epsilon_{II}$', size=10.5)
        ds.label(ax, 3.8, .3, 'Unit cell, plane strain; every node prescribed', size=10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def anisotropic_time(folder):
    with ds.canvas(7.2, 4.6) as (fig, ax):
        x0, y0, side = 1.15, 1.05, 2.5
        ax.add_patch(Rectangle((x0, y0), side, side, fc=WET, ec=ds.WATER, lw=ds.THICK))
        grains(ax, x0, y0, side, side, 55)
        ds.label(ax, x0+side/2, y0+side+.35, '$p=0$ on all four faces', ds.WATER, 10.5)
        ds.dimension(ax, (x0, y0), (x0+side, y0), '$1$ m', .3, 'below')
        ds.dimension(ax, (x0, y0), (x0, y0+side), '$1$ m', .3, 'left')
        ds.axes(ax, (x0+.45, y0+.5), .65)
        ds.label(ax, x0+side/2, y0+side-.55, 'two pressure modes', size=10.5)
        ds.label(ax, 5.3, 3.5, 'permeability changes at 0.5 s', size=10.5)
        ds.label(ax, 5.3, 2.85, r'$\mathbf{k}_1=k_0\,\mathrm{diag}(4,1)$', ds.WATER, 10.5)
        ds.arrow(ax, (5.3, 2.5), (5.3, 2.15), color=ds.WATER)
        ds.label(ax, 5.3, 1.8, r'$\mathbf{k}_2=k_0\,\mathrm{diag}(2,18)$', ds.WATER, 10.5)
        ds.label(ax, 5.3, 1.1, 'pressure stays continuous', size=10)
        ds.label(ax, 3.6, .25, 'Fixed geometry; strain feedback into pressure disabled', ds.MUTED, 9.5)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')



def slope(folder):
    """A soil layer on inclined bedrock, strong ends, pore pressure rising from the bedrock."""
    with ds.canvas(7.8, 4.6) as (fig, ax):
        t = np.tan(np.radians(20)); x0, z0, Lw, Dw = .6, 3.7, 5.6, .55
        top = lambda x: z0 - (x - x0)*t
        def band(a, b, color):
            ax.add_patch(Polygon([(a, top(a)), (b, top(b)), (b, top(b) - Dw), (a, top(a) - Dw)],
                                 closed=True, fc=color, ec='none', zorder=1))
        band(x0, x0 + .6, '#d9d6cf'); band(x0 + .6, x0 + Lw - .6, WET); band(x0 + Lw - .6, x0 + Lw, '#d9d6cf')
        outline = [(x0, top(x0)), (x0 + Lw, top(x0 + Lw)), (x0 + Lw, top(x0 + Lw) - Dw), (x0, top(x0) - Dw)]
        ax.add_patch(Polygon(outline, closed=True, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        bed = [(x0, top(x0) - Dw), (x0 + Lw, top(x0 + Lw) - Dw), (x0 + Lw, top(x0 + Lw) - Dw - .35), (x0, top(x0) - Dw - .35)]
        ax.add_patch(Polygon(bed, closed=True, fc='none', ec=ds.INK, lw=ds.THIN, hatch='////', zorder=2))
        for xx in (x0 + 1.0, x0 + 1.6, x0 + Lw - 1.6, x0 + Lw - 1.0):
            ds.arrow(ax, (xx, top(xx) - Dw - .02), (xx, top(xx) - Dw + .32), ds.WATER, scale=8)
        ds.label(ax, 3.4, 4.45, r'bedrock: no slip, water supplied, $p_b(t)$ rising', ds.WATER, 9.5, ha='left')
        ds.leader(ax, (x0 + 4.0, top(x0 + 4.0)), (x0 + 4.4, top(x0 + 4.0) + .5), r'surface: free, $p=0$', size=9.5)
        ds.label(ax, x0 + .3, top(x0) + .32, 'strong', ds.MUTED, 8.5)
        ds.label(ax, x0 + Lw - .3, top(x0 + Lw) + .32, 'strong', ds.MUTED, 8.5)
        ds.label(ax, x0 + Lw/2, top(x0 + Lw/2) - Dw/2 + .02, r'Mohr–Coulomb, $\psi=0$', size=9.5, rotation=-20)
        ds.label(ax, 3.4, 4.05, r'$\beta=20^\circ$, $D=2$ m, $L=40$ m; rollers at both ends', size=9.5, ha='left')
        ds.label(ax, 3.4, 3.68, 'gravity applied as a load from an unstressed start', ds.MUTED, 8.5, ha='left')
        ds.label(ax, 7.6, .35, 'not to scale', ds.MUTED, 7.5, ha='right')
        ds.gravity(ax, .6, 1.6)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')



def rainfall(folder):
    """The slope of the pressure-driven page under rain: sealed bedrock, draining toe, perched table."""
    with ds.canvas(7.8, 4.6) as (fig, ax):
        t = np.tan(np.radians(20)); x0, z0, Lw, Dw = .6, 3.7, 5.6, .55
        top = lambda x: z0 - (x - x0)*t
        def band(a, b, color, fraction=1.0):
            ax.add_patch(Polygon([(a, top(a) - Dw*(1 - fraction)), (b, top(b) - Dw*(1 - fraction)), (b, top(b) - Dw), (a, top(a) - Dw)],
                                 closed=True, fc=color, ec='none', zorder=1))
        band(x0, x0 + Lw, '#f1ede4')                    # unsaturated soil
        band(x0 + .6, x0 + Lw - .6, WET, .55)          # perched water above the bedrock
        band(x0, x0 + .6, '#d9d6cf'); band(x0 + Lw - .6, x0 + Lw, '#d9d6cf')
        ax.plot([x0 + .6, x0 + Lw - .6], [top(x0 + .6) - .45*Dw, top(x0 + Lw - .6) - .45*Dw], color=ds.WATER, lw=ds.THIN, ls=(0, (4, 3)), zorder=3)
        outline = [(x0, top(x0)), (x0 + Lw, top(x0 + Lw)), (x0 + Lw, top(x0 + Lw) - Dw), (x0, top(x0) - Dw)]
        ax.add_patch(Polygon(outline, closed=True, fill=False, ec=ds.INK, lw=ds.THICK, zorder=4))
        bed = [(x0, top(x0) - Dw), (x0 + Lw, top(x0 + Lw) - Dw), (x0 + Lw, top(x0 + Lw) - Dw - .35), (x0, top(x0) - Dw - .35)]
        ax.add_patch(Polygon(bed, closed=True, fc='none', ec=ds.INK, lw=ds.THIN, hatch='////', zorder=2))
        for xx in np.linspace(x0 + .9, x0 + Lw - 1.9, 5):
            ds.arrow(ax, (xx, top(xx) + .55), (xx, top(xx) + .05), ds.WATER, scale=8)
        ax.plot([x0 + Lw, x0 + Lw], [top(x0 + Lw) - Dw, top(x0 + Lw)], color=ds.WATER, lw=2.6, zorder=5, solid_capstyle='butt')
        ds.label(ax, 3.4, 4.45, r'rain $q=10^{-5}$ m/s; the surface never saturates', ds.WATER, 9.5, ha='left')
        ds.label(ax, 3.4, 4.05, r'bedrock sealed and fixed; toe drains, $p=0$', size=9.5, ha='left')
        ds.label(ax, 3.4, 3.68, r'initial suction 5 kPa; Gardner $S_e(p)$, Bishop $\chi=S_e$', ds.MUTED, 8.5, ha='left')
        ds.leader(ax, (x0 + 4.0, top(x0 + 4.0) - .45*Dw), (x0 + 4.4, top(x0 + 4.0) + .15), 'perched water table', ds.WATER, 9.5)
        ds.label(ax, 7.6, .35, r'$\beta=20^\circ$, $D=2$ m, $L=40$ m; strong ends; not to scale', ds.MUTED, 7.5, ha='right')
        ds.gravity(ax, .6, 1.6)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def drawdown(folder):
    """The saturated embankment with the reservoir at the crest, drawn down to the toe platform."""
    with ds.canvas(7.8, 4.6) as (fig, ax):
        k, x0, y0 = .1, .5, .9
        P = lambda x, z: (x0 + k*x, y0 + k*(z + 20))
        ax.add_patch(Polygon([P(0, -10), P(20, -10), P(40, 0), P(0, 0)], closed=True, fc=WET, ec='none', alpha=.55, zorder=1))
        outline = [P(0, -20), P(60, -20), P(60, 0), P(40, 0), P(20, -10), P(0, -10)]
        ax.add_patch(Polygon(outline, closed=True, fc='#dbe8ef', ec=ds.INK, lw=ds.THICK, zorder=2))
        ax.add_patch(Polygon([P(0, -20), P(60, -20), P(60, -21.8), P(0, -21.8)], closed=True, fc='none', ec=ds.INK,
                             lw=ds.THIN, hatch='////', zorder=1))
        for z in (-17, -13):
            ds.roller(ax, *P(0, z), (-1, 0))
        for z in (-16, -10, -4):
            ds.roller(ax, *P(60, z), (1, 0))
        water_table(ax, *P(9, 0))
        ax.plot(*zip(P(0, 0), P(40, 0)), color=ds.WATER, lw=ds.THIN, ls=(0, (4, 3)), zorder=3)
        ds.label(ax, *P(10, 4.0), 'reservoir before drawdown', ds.WATER, 9)
        ds.arrow(ax, P(10, -1.5), P(10, -8.5), ds.WATER, scale=9)
        ds.label(ax, *P(10.8, -5), 'rapid\ndrawdown', ds.WATER, 8.5, ha='left')
        ax.plot(*zip(P(0, -10), P(20, -10)), color=ds.WATER, lw=2.4, zorder=4, solid_capstyle='butt')
        ax.plot(*zip(P(20, -10), P(40, 0)), color=ds.WATER, lw=2.4, zorder=4, solid_capstyle='butt')
        ds.leader(ax, P(31, -5.5), P(36, -10.5), 'after: face and toe drain, $p=0$', ds.WATER, 9)
        ds.label(ax, *P(50, 1.6), 'crest: water level, then sealed', size=8.5)
        ds.label(ax, *P(30, -17), r'saturated, $k=10^{-13}$ m$^2$; Gardner $S_e(p)$, Bishop $\chi=S_e$', size=9)
        ds.label(ax, *P(30, 6.5), r'$c=10$ kPa, $\phi=20^\circ$, $\psi=0$; base fixed, sides on rollers', ds.MUTED, 8.5)
        ds.gravity(ax, *P(54, -12))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def embankment(folder):
    """The homogeneous embankment of the strength-reduction page, drawn to scale (1 unit = 10 m)."""
    with ds.canvas(7.8, 4.6) as (fig, ax):
        k, x0, y0 = .1, .5, .9                            # m -> drawing units; base-left corner
        P = lambda x, z: (x0 + k*x, y0 + k*(z + 20))
        outline = [P(0, -20), P(60, -20), P(60, 0), P(40, 0), P(20, -10), P(0, -10)]
        ax.add_patch(Polygon(outline, closed=True, fc=DRY, ec=ds.INK, lw=ds.THICK, zorder=2))
        ax.add_patch(Polygon([P(0, -20), P(60, -20), P(60, -21.8), P(0, -21.8)], closed=True, fc='none', ec=ds.INK,
                             lw=ds.THIN, hatch='////', zorder=1))
        for z in (-17, -13):
            ds.roller(ax, *P(0, z), (-1, 0))
        for z in (-16, -10, -4):
            ds.roller(ax, *P(60, z), (1, 0))
        # Bishop's critical circle (entry at the toe, exit on the crest) for scale.
        th = np.linspace(np.arccos((20 - 23.41)/22.89), np.arccos((42.5 - 23.41)/22.89), 60)
        cx, cz = 23.41 + 22.89*np.cos(th), 12.64 - 22.89*np.sin(th)
        ax.plot(*zip(*[P(a, b) for a, b in zip(cx, cz)]), color=ds.HOT, lw=ds.BODY, ls=(0, (5, 3)), zorder=3)
        ds.leader(ax, P(35, -7.1), P(39, -12.5), "Bishop's critical circle", ds.HOT, 9)
        ax.annotate('', xy=P(63, 0), xytext=P(63, -10), arrowprops=dict(arrowstyle='<->', color=ds.INK, lw=ds.THIN))
        ds.label(ax, *P(64, -5), r'$H=10$ m', size=10, ha='left')
        ax.annotate('', xy=P(63, -10), xytext=P(63, -20), arrowprops=dict(arrowstyle='<->', color=ds.INK, lw=ds.THIN))
        ds.label(ax, *P(64, -15), r'$H$', size=10, ha='left')
        ds.label(ax, *P(28, -8.1), '2', size=9.5); ds.label(ax, *P(31, -6), '1', size=9.5)
        ax.plot(*zip(P(26, -7), P(30, -7), P(30, -5)), color=ds.INK, lw=ds.THIN)
        ds.label(ax, *P(30, -18), r'$c=10$ kPa, $\phi=20^\circ$, $\psi=0$, $\gamma=20$ kN/m$^3$', size=9.5)
        ds.label(ax, *P(30, 6.5), r'strength reduction: $c/F$, $\tan\phi/F$; the largest $F$ that stands is FS', size=9.5)
        ds.label(ax, *P(30, 3.6), 'base fixed, sides on rollers; gravity applied elastically first', ds.MUTED, 8.5)
        ds.gravity(ax, *P(8, -3))
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')



def content(folder, inflow=1e-6):
    with ds.canvas(7.4, 4.7) as (fig, ax):
        ax.fill_between([2.5, 3.5], 1.2, 3.5, color='#eef3f7')
        ds.wall(ax, 2.28, 1.2, .12, 2.3)
        ds.wall(ax, 3.6, 1.2, .12, 2.3)
        for y in (1.65, 2.35, 3.05):
            ds.roller(ax, 2.5, y, (-1, 0), .09)
            ds.roller(ax, 3.5, y, (1, 0), .09)
        ax.plot([2.5, 3.5], [3.5, 3.5], color=ds.WATER, lw=ds.THICK)
        ax.plot([2.5, 3.5], [1.2, 1.2], color=ds.INK, lw=ds.THICK)
        for x in (2.7, 3.3):
            ds.arrow(ax, (x, .55), (x, 1.15), ds.WATER)
            ds.arrow(ax, (x, 3.58), (x, 4.05), ds.WATER)
        exponent = int(np.floor(np.log10(inflow)))
        coefficient = inflow/10**exponent
        value = ('' if coefficient == 1 else f'{coefficient:g}\\times') + f'10^{{{exponent}}}'
        ds.label(ax, 3., .22, f'$q_{{\\rm in}}={value}$ m/s', ds.WATER)
        ds.dimension(ax, (2.5, 1.2), (2.5, 3.5), '$H=1$ m', .7)
        ds.leader(ax, (3.5, 3.5), (4.5, 4.12), 'Drain: $p=0$', ds.WATER)
        ds.label(ax, 5.65, 3.57, r'Free top: $\sigma_{zz}=0$', size=10.5)
        ds.leader(ax, (3.7, 2.6), (4.5, 2.8), 'Sealed sides', size=10.5)
        ds.label(ax, 5.45, 2.38, '$u_x=0$ (and $u_y=0$ in 3D)', size=10)
        ds.leader(ax, (3.5, 1.2), (4.5, 1.45), 'Porous rigid base', size=10.5)
        ds.label(ax, 5.55, 1.03, '$u_z=0$; prescribed inflow', size=10)
        ds.arrow(ax, (1.02, 3.5), (1.02, 2.8), ds.MUTED)
        ds.label(ax, 1.02, 3.8, '$d$', ds.MUTED)
        ds.label(ax, 3., 2.35, 'Elastic\nporous\nsolid', size=10)
        ds.label(ax, 1.05, .65, '$g=0$', size=10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def thermal_cell(folder):
    with ds.canvas(7.8, 4.6) as (fig, ax):
        ax.add_patch(pf.plt.Rectangle((2., 1.2), 2., 2., fc='#f6ede7', ec=ds.INK, lw=ds.BODY))
        ds.wall(ax, 1.8, 1.2, .12, 2.)
        ds.wall(ax, 4.08, 1.2, .12, 2.)
        ds.wall(ax, 2., .98, 2., .12)
        ds.wall(ax, 2., 3.3, 2., .12)
        for x in (2.4, 3., 3.6):
            ds.arrow(ax, (x, 2.65), (x, 2.17), ds.HOT)
        ds.label(ax, 3., 1.65, r'$T(t),\ p(t)$', size=12)
        ds.dimension(ax, (2., 1.2), (2., 3.2), '1 m', .85)
        ds.leader(ax, (4.18, 3.), (4.8, 3.7), 'Rigid, sealed walls', size=10.5)
        ds.label(ax, 5.85, 3.25, r'$\varepsilon_v=0,\ q_n=0$', size=11)
        ds.leader(ax, (3.55, 2.42), (4.8, 2.5), 'Uniform supplied heat', ds.HOT, size=10.5)
        ds.label(ax, 5.8, 2.05, r'$\dot Q/V=C_0\,(40\ {\rm K/s})$', ds.HOT, size=10)
        ds.label(ax, 5.8, 1.55, '$0<t<0.25$ s; then hold', size=10)
        ds.label(ax, 3., .52, r'$g=0,\quad \alpha_s=0$', size=11)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def plastic_heat(folder):
    with ds.canvas(7.8, 5.5) as (fig, ax):
        x0, x1, z0, z1 = 1.8, 4., 1.6, 3.8
        ax.add_patch(Rectangle((x0, z0), x1-x0, z1-z0, fc=WET, ec=ds.INK, lw=ds.THICK))
        for z in np.linspace(z0+.4, z1-.4, 3):
            ds.roller(ax, x0, z, (-1, 0), .12)
            ds.arrow(ax, (x1+.04, z), (x1+.6, z), ds.LOAD)
        for x in np.linspace(x0+.4, x1-.4, 3):
            ds.roller(ax, x, z1, (0, 1), .12)
            ds.arrow(ax, (x, z0-.6), (x, z0-.04), ds.LOAD)
        ds.label(ax, 2.9, 4.65, r'$v_z=0$', size=11)
        ds.label(ax, 2.9, .55, r'$v_z=+0.5$ m/s', ds.LOAD, size=11)
        ds.dimension(ax, (x0, z0), (x0, z1), '$1$ m', .85)
        ds.label(ax, 5.9, 3.55, r'$v_x=+0.5$ m/s', ds.LOAD, size=11)
        ds.label(ax, 5.9, 2.8, 'Sealed and insulated', ds.WATER, size=10.5)
        ds.label(ax, 5.9, 2.25, 'on every face', ds.WATER, size=10.5)
        ds.label(ax, 5.9, 1.4, r'$\sigma_0=-30$ MPa $\mathbf{I}$', size=10)
        ds.label(ax, 5.9, .9, r'$T_0=300$ K, $p_0=0$', size=10)
        ds.label(ax, 2.9, 2.9, 'Porous EP-RSF', size=10.5)
        ds.label(ax, 2.9, 2.35, r'$\dot\varepsilon_v=0$', size=11)
        ds.label(ax, .6, 4.65, r'$v_x=0$', size=10.5)
        ds.axes(ax, (6.65, 4.3), length=.4)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def plastic_content(folder):
    with ds.canvas(7.8,4.8) as (fig,ax):
        ax.fill_between([2.5,3.5],1.3,3.5,color='#eef3f7')
        ds.wall(ax,2.28,1.3,.12,2.2);ds.wall(ax,3.6,1.3,.12,2.2)
        ds.wall(ax,2.5,1.08,1.,.12)
        ax.plot([2.5,3.5],[3.5,3.5],color=ds.INK,lw=2)
        for y in (1.7, 2.4, 3.1):
            ds.roller(ax, 2.5, y, (-1, 0), .09)
            ds.roller(ax, 3.5, y, (1, 0), .09)
        for x in (2.7, 3.3):
            ds.arrow(ax, (x, 4.2), (x, 3.56), ds.HOT)
        ds.label(ax,3.,4.48,r'$v_z=-10^{-3}$ m/s',ds.HOT,size=11)
        ds.dimension(ax,(2.5,1.3),(2.5,3.5),'$H=1$ m',.8)
        ds.label(ax,3.,2.4,'Porous\nplastic\nskeleton',size=10)
        ds.leader(ax,(3.7,3.15),(4.65,3.65),'Side-normal motion fixed',size=10.5)
        ds.label(ax,5.95,3.18,'All water boundaries sealed',ds.WATER,size=10.5)
        ds.leader(ax,(3.5,1.3),(4.65,1.6),'Fixed base',size=10.5)
        ds.label(ax,5.65,1.15,'$u_z=0$',size=11)
        ds.label(ax,3.,.55,'$g=0$',size=11)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')

def pore_recovery(folder):
    """Drained material cell: tensile opening followed by a fixed-strain hold."""
    with ds.canvas(8, 5.2) as (fig, ax):
        x0, x1, z0, z1 = 2.2, 4.4, 1.4, 3.6
        ax.add_patch(Rectangle((x0, z0), x1-x0, z1-z0, fc=SOLID, ec=ds.INK, lw=ds.THICK))
        ds.wall(ax, x0, z0-.16, x1-x0, .12)
        for z in (1.8, 2.5, 3.2):
            ds.roller(ax, x0, z, (-1, 0), .11)
            ds.roller(ax, x1, z, (1, 0), .11)
        for x in (2.65, 3.3, 3.95):
            ds.arrow(ax, (x, z1+.05), (x, z1+.55), ds.LOAD)
        ds.label(ax, 3.3, 4.5, r'$v_z=10^{-3}$ m/s, then $v_z=0$', ds.LOAD, 11)
        ds.dimension(ax, (x0, z0), (x0, z1), '$H=1$ m', .75)
        ds.dimension(ax, (x0, z0), (x1, z0), '$W=1$ m', .48, side='below')
        ds.label(ax, 3.3, 2.6, 'Plastic\nporous solid', size=11)
        ds.leader(ax, (x1, 2), (5.3, 1.8), 'All vertices: $p=0$', ds.WATER, 10.5)
        ds.label(ax, 6.2, 3.3, '$T=293$ K\n$g=0$', size=11)
        ds.axes(ax, (6.15, 2.2), length=.45)
        ds.label(ax, 4, .35, 'Fixed mesh; plane strain in 2D, zero lateral y strain in 3D', size=10)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def fault_heat(folder):
    """Reference mesh, imposed affine slip and the free normal platen."""
    with ds.canvas(8.2, 5.4) as (fig, ax):
        x0, x1, z0, z1 = 1.7, 4.1, 1.4, 3.8
        ax.add_patch(Rectangle((x0, z0), x1-x0, z1-z0, fc=WET, ec=ds.INK, lw=ds.THICK))
        ax.plot([x0, x1], [z1, z0], color=ds.MUTED, lw=ds.THIN, ls=(0, (4, 3)))
        ds.wall(ax, x0, z0-.14, x1-x0, .14)
        ds.wall(ax, x0, z1, x1-x0, .14)
        for x in (2.15, 2.9, 3.65):
            ds.arrow(ax, (x, 4.4), (x, z1+.17), ds.LOAD)
        ds.arrow(ax, (2.1, 4.78), (3.7, 4.78), ds.LOAD, scale=11)
        ds.label(ax, 2.9, 5.12, '$V=1$ m/s', ds.LOAD, size=11)
        for fraction in (.25, .5, .75):
            z = z0+fraction*(z1-z0)
            ds.arrow(ax, (x1+.06, z), (x1+.06+.9*fraction, z), ds.LOAD)
        ds.leader(ax, (x1, z1+.08), (5.1, 4.35), '$N=50$ MPa', ds.LOAD, 11)
        ds.label(ax, 6.2, 3.65, r'$v_x=V(z+H)/H$', ds.LOAD, size=10.5)
        ds.label(ax, 6.2, 3.15, 'on both side faces', ds.LOAD, size=10)
        ds.label(ax, 3.05, 3.25, '$T_0=300$ K', size=10.5)
        ds.label(ax, 2.55, 2., '$p_0=0$', size=11)
        ds.label(ax, 6.2, 2.45, 'Sealed and insulated', ds.WATER, size=10.5)
        ds.label(ax, 6.2, 1.95, 'on every face', ds.WATER, size=10)
        ds.label(ax, 6.2, 1.35, '$g=0$', size=11)
        ds.dimension(ax, (x0, z0), (x0, z1), '$H=10$ mm', .75)
        ds.label(ax, 2.9, .9, r'base: $v_x=v_z=0$', size=10.5)
        ds.label(ax, 2.9, .35, 'Reference geometry; two triangles', ds.MUTED, size=10)
        ds.axes(ax, (7.35, .45), length=.4)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


def inventory_cell(folder, kind):
    """Shared cell geometry; each problem shows its own motion and exchange."""
    with ds.canvas(8., 5.1) as (fig, ax):
        x0, x1, z0, z1 = 2.05, 4.05, 1.35, 3.35
        ax.add_patch(Rectangle((x0, z0), x1-x0, z1-z0,
                               fc=WET if kind=='mass' else '#f6ede7', ec=ds.INK, lw=ds.THICK))
        ds.wall(ax, x0, z0-.2, x1-x0, .1)
        ds.dimension(ax, (x0, z0), (x0, z1), '1 m', .8)
        ds.label(ax, 3.05, 2.45, 'Saturated cell', size=11)
        ds.label(ax, 3.05, 1.98, r'$M(t),\ T(t),\ p(t)$', size=11)
        ds.label(ax, 3.05, .63, 'unit thickness; $g=0$', ds.MUTED, size=10)
        if kind == 'mass':
            side_rollers(ax, x0, x1, z0, z1, size=.13)
            down_arrows(ax, x0+.1, x1-.1, z1, .55, ds.LOAD, n=3, bar=True)
            ds.label(ax, 3.05, 4.35, r'$v_z=-0.01$ m/s', ds.LOAD, size=11)
            ds.label(ax, 6., 3.45, 'Sealed on every face', ds.WATER, size=10.5)
            ds.label(ax, 6., 2.85, r'$T=273$ K', size=11)
            ds.label(ax, 6., 2.25, r'$0\leq t\leq0.1$ s', size=11)
            ds.label(ax, 6., 1.65, 'Fixed side-normal motion', size=10)
        else:
            ds.wall(ax, x0-.2, z0, .1, z1-z0)
            ds.wall(ax, x0, z1+.1, x1-x0, .1)
            ds.label(ax, 6., 3.95, 'Rigid; $p=0$ on all faces', size=10.5)
            ds.arrow(ax, (x1+.05, 2.85), (4.85, 2.85), ds.WATER)
            ds.label(ax, 6.15, 2.85, 'Outflow at $T(t)$', ds.WATER, size=10.5)
            ds.label(ax, 6., 1.95, 'Insulated thermal walls', size=10.5)
            ds.label(ax, 6., 1.35, r'$T_0=273$ K; $0\leq t\leq5$ s', size=10)
            if kind == 'drained_heat':
                for x in (2.55, 3.05, 3.55):
                    ds.arrow(ax, (x, 2.98), (x, 2.68), ds.HOT)
                ds.label(ax, 3.05, 4.3, r'$\dot Q=30$ MW', ds.HOT, size=12)
            else:
                ds.arrow(ax, (.35, 2.85), (x0-.25, 2.85), ds.HOT)
                ds.label(ax, .95, 3.5, r'$Q_i=0.01$ m$^3$/s', ds.HOT, size=10)
                ds.label(ax, .95, 1.2, '$T_i=323$ K', ds.HOT, size=10.5)
                ds.label(ax, 3.05, 4.3, 'Uniform mixing', size=11)
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


DRAW = dict(content=content, thermal_cell=thermal_cell, plastic_content=plastic_content, drawdown=drawdown, embankment=embankment, rainfall=rainfall, slope=slope, anisotropic_time=anisotropic_time, rsf=rsf, borehole=borehole, injection_temperature=injection_temperature, drainage=drainage, pressure_yield=pressure_yield, strength=strength, regional=regional, nusselt=schematic_convection, bishop_cell=bishop_cell, mctigue=mctigue, bishop=bishop, hantush=hantush, kirsch=kirsch, salencon=salencon, theis_wells=theis_wells, open_flow=open_flow, rudnicki=rudnicki, mandel=mandel, cryer=cryer, convection=schematic_convection,heat=heat, ogata_banks=heat, point_sources=sources, theis=theis, averaged_aquifer=aquifer,
            **{k: (lambda folder, k=k: column(k, folder)) for k in COLUMN})

def convection_cube(folder):
    """Axonometric cube with the actual thermal and water boundaries."""
    with ds.canvas(6.8, 5.5) as (fig, ax):
        origin = np.array([1.5, 1.3]); ex = np.array([2.45, 0.])
        ey = np.array([.85, .7]); ez = np.array([0., 2.4])
        point = lambda i, j, k: origin+i*ex+j*ey+k*ez
        for k, color in ((0, ds.HOT), (1, ds.COLD)):
            face = [point(0, 0, k), point(1, 0, k), point(1, 1, k), point(0, 1, k)]
            ax.add_patch(Polygon(face, fc=color, ec=color, alpha=.15, lw=ds.BODY))
            ax.plot(*np.array(face+[face[0]]).T, color=color, lw=ds.THICK)
        for i, j in ((0, 0), (1, 0), (1, 1), (0, 1)):
            a, b = point(i, j, 0), point(i, j, 1)
            ax.plot([a[0], b[0]], [a[1], b[1]], color=ds.INK, lw=ds.BODY,
                    ls='--' if (i, j) == (0, 1) else '-')
        ds.label(ax, 3.1, 4.85, r'top: $T=T_0$', ds.COLD, 11)
        ds.label(ax, 2.8, .64, r'base: $T=T_0+\Delta T$', ds.HOT, 11)
        ds.leader(ax, point(1, 1, .6), (5.0, 3.45), 'sealed on all faces', ds.WATER, 10)
        ds.label(ax, 5.06, 3.02, 'vertical walls insulated', ds.MUTED, 9.5, ha='left')
        clear = dict(facecolor='white', edgecolor='none', pad=2)
        ds.label(ax, 2.8, 2.65, 'saturated porous cube', ds.MUTED, 10, bbox=clear)
        ds.label(ax, 2.8, 2.16, r'$H=L_x=L_y=1$ m', ds.MUTED, 10, bbox=clear)
        ds.dimension(ax, point(0, 0, 0), point(0, 0, 1), '$H$', .42, 'left')
        ds.gravity(ax, 6.15, 4.8)
        o = np.array([5.3, 1.1])
        for delta, label in (([.65, 0], 'x'), ([.3, .4], 'y'), ([0, .7], 'z')):
            end = o+delta; ds.arrow(ax, o, end, ds.MUTED, scale=8)
            ds.label(ax, end[0]+.12, end[1]+.08, '$'+label+'$', ds.MUTED, 9)
        ds.label(ax, .2, .15, 'schematic; rigid matrix, fixed mesh', ds.MUTED, 8, ha='left')
        ds.check_labels(fig, ax)
        pf.save_figure(fig, folder, 'setup')


DRAW['convection_3d'] = convection_cube

for _key in ('mass', 'drained_heat', 'mixing'):
    DRAW[_key] = lambda folder, key=_key: inventory_cell(folder, key)

DRAW['plastic_heat'] = plastic_heat
DRAW['fault_heat'] = fault_heat


def schematic(key, folder):
    """Render each registered setup through the same drawing primitives."""
    DRAW[key](folder)


if __name__ == '__main__':
    import review
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('review', type=Path)
    p.add_argument('--only', nargs='*')
    a = p.parse_args()
    for key in a.only or [k for k in review.SPECS if k in DRAW]:
        folder = a.review/'fig'/review.SLUGS.get(key, key.replace('_', '-'))
        folder.mkdir(parents=True, exist_ok=True)
        DRAW[key](folder)
        print('drew', key, '->', folder/'setup.svg')
