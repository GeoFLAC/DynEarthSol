"""Textbook-style primitives for the review's setup drawings.

One serif face for text and mathematics (STIX), three line weights, standard
boundary symbols (hatched walls and plates, roller supports), dimension lines
with extension lines, and arrowheads placed on curves along their tangent.
Colour is reserved for physical meaning: hot, cold and water.
"""
from contextlib import contextmanager

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle, Circle

INK, MUTED, HOT, COLD, WATER, LOAD = '#1f2a30', '#5f6b72', '#b23a24', '#24618f', '#1d5a7a', '#8f2f1e'
THICK, BODY, THIN = 1.5, 1.0, .6
# Drawings are shown about 500 px wide; this scale keeps labels at body-text size.
TYPE = 1.45
RC = {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 11,
      'hatch.linewidth': .6, 'hatch.color': INK, 'svg.fonttype': 'path'}


@contextmanager
def canvas(width, height):
    """A figure in drawing units of inches; yields (fig, ax)."""
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(width, height))
        ax.set(xlim=(0, width), ylim=(0, height)); ax.set_aspect('equal'); ax.axis('off')
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        yield fig, ax


def label(ax, x, y, value, color=INK, size=11, ha='center', va='center', **kw):
    return ax.text(x, y, value, color=color, fontsize=TYPE*size, ha=ha, va=va, **kw)


def arrow(ax, start, end, color=INK, lw=BODY, two=False, scale=10):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='<|-|>' if two else '-|>', color=color,
                                 lw=lw, mutation_scale=scale, shrinkA=0, shrinkB=0, zorder=6))


def leader(ax, target, at, text, color=INK, size=10.5, ha='left'):
    """Thin leader line from a feature to its label."""
    ax.plot([target[0], at[0]], [target[1], at[1]], color=MUTED, lw=THIN, zorder=5)
    ax.plot(*target, 'o', ms=2.2, color=MUTED, zorder=5)
    pad = .06 if ha == 'left' else -.06
    return label(ax, at[0]+pad, at[1], text, color, size, ha=ha)


def dimension(ax, a, b, text, offset, side='left'):
    """Engineering dimension: extension lines, open-ended arrows, centred label."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    along = (b-a)/np.linalg.norm(b-a); normal = np.array([-along[1], along[0]])
    if side in ('right', 'below'): normal = -normal
    p, q = a+offset*normal, b+offset*normal
    for base, tip in ((a, p), (b, q)):
        start, end = base+.06*normal, tip+.07*normal
        ax.plot([start[0], end[0]], [start[1], end[1]], color=MUTED, lw=THIN)
    arrow(ax, p, q, lw=THIN+.1, two=True, scale=8)
    mid = (p+q)/2+.16*normal
    angle = np.degrees(np.arctan2(along[1], along[0]))
    if angle > 90 or angle < -90: angle += 180
    return label(ax, *mid, text, rotation=angle if abs(along[1]) > .5 else 0)


def wall(ax, x, y, width, height, hatch='////'):
    """Rigid impermeable wall drawn as a hatched strip."""
    ax.add_patch(Rectangle((x, y), width, height, fc='white', ec=INK, lw=BODY, hatch=hatch, zorder=3))


def plate(ax, x0, x1, y, side, color, depth=.14):
    """Isothermal plate: a heavy coloured face with hatching on its outer side."""
    yb = y if side > 0 else y-depth
    # Hatching takes the edge colour; a zero-width edge keeps only the hatch.
    ax.add_patch(Rectangle((x0, yb), x1-x0, depth, fc='white', ec=INK, lw=0, hatch='////', zorder=3))
    ax.plot([x0, x1], [y, y], color=color, lw=3.2, solid_capstyle='butt', zorder=4)


def roller(ax, x, y, normal, size=.12):
    """Roller support: zero normal displacement, free tangential sliding."""
    n = np.asarray(normal, float); t = np.array([-n[1], n[0]]); tip = np.array([x, y]); base = tip+size*n
    ax.add_patch(Polygon([tip, base+.55*size*t, base-.55*size*t], fill=False, ec=INK, lw=THIN+.2, zorder=4))
    for s in (-1, 1):
        ax.add_patch(Circle(base+.3*size*n+s*.3*size*t, .18*size, fill=False, ec=INK, lw=THIN, zorder=4))
    rail = base+.55*size*n
    ax.plot(*np.column_stack((rail-.8*size*t, rail+.8*size*t)), color=INK, lw=THIN+.2)
    ground = rail+.12*size*n
    ax.add_patch(Polygon([rail-.8*size*t, rail+.8*size*t, ground+.8*size*t, ground-.8*size*t],
                         fc='white', ec=INK, lw=0, hatch='/////', zorder=3))


def axes(ax, origin, length=.45, names=('x', 'z'), color=MUTED):
    o = np.asarray(origin, float)
    arrow(ax, o, o+(length, 0), color, THIN+.2, scale=8)
    arrow(ax, o, o+(0, length), color, THIN+.2, scale=8)
    label(ax, o[0]+length+.1, o[1], f'${names[0]}$', color, 10)
    label(ax, o[0], o[1]+length+.11, f'${names[1]}$', color, 10)


def gravity(ax, x, y, length=.55):
    arrow(ax, (x, y), (x, y-length), INK, BODY, scale=11)
    label(ax, x+.16, y-length/2, '$g$', size=12)


def on_curve(ax, path, at, color, reverse=False, scale=11):
    """Arrowhead on a polyline at the vertex nearest `at`, along its tangent."""
    path = np.asarray(path)
    i = int(np.argmin(np.hypot(*(path-np.asarray(at)).T)))
    i = min(max(i, 1), len(path)-2)
    a, b = (path[i+1], path[i-1]) if reverse else (path[i-1], path[i+1])
    d = (b-a)/np.linalg.norm(b-a)
    c = path[i]
    ax.add_patch(FancyArrowPatch(c-.04*d, c+.04*d, arrowstyle='-|>', color=color, lw=0,
                                 mutation_scale=scale, zorder=6))


def check_labels(fig, ax, strict=True):
    """No two labels may overlap; returns the overlapping pairs."""
    fig.canvas.draw()
    boxes = [(t.get_text(), t.get_window_extent()) for t in ax.texts if t.get_text()]
    clashes = [(ti, tj) for i, (ti, bi) in enumerate(boxes) for tj, bj in boxes[i+1:] if bi.overlaps(bj)]
    assert not (strict and clashes), clashes
    return clashes


# ----------------------------------------------------------------- result figures
# Every result figure uses the same face, colours and marks: the analytical
# solution as a dark solid line, DES as open markers, lettered panel titles.
FIGURE_RC = {**RC, 'font.size': 13, 'axes.edgecolor': INK, 'axes.linewidth': .8, 'axes.labelcolor': INK,
             'xtick.color': INK, 'ytick.color': INK, 'axes.titlesize': 13, 'legend.frameon': False,
             'legend.fontsize': 11, 'axes.grid': True, 'grid.alpha': .15, 'grid.linewidth': .6}
SERIES = ['#1f4e79', '#5b8cbf', '#b8432a', '#7a1f12']   # cool to warm, for times or cases


def result_figure(panels, width=None):
    """Figure and axes in the shared style; call inside plt.rc_context(FIGURE_RC)."""
    width = width or {1: 6.4, 2: 9.0, 3: 11.0}[panels]
    fig, axes = plt.subplots(1, panels, figsize=(width, 3.9), constrained_layout=True)
    return fig, np.atleast_1d(axes)


def analytical(ax, x, y, color=INK, label='analytical', **kw):
    return ax.plot(x, y, '-', color=color, lw=1.2, label=label, **kw)


def des(ax, x, y, color=HOT, label='DES', size=3.5, **kw):
    return ax.plot(x, y, 'o', ms=size, mfc='white', mec=color, mew=.9, ls='none', label=label, **kw)


def panel(ax, letter, title, **labels):
    """Title '(a) ...' when a figure has several panels; plain title otherwise."""
    ax.set(title=f'({letter}) {title}' if letter else title, **labels)
    if ax.get_legend_handles_labels()[0]:
        ax.legend()


def band(ax, a, b, text='compared band', x=None, y=.97, rotation=0, ha='left', va='top'):
    """Shade a range of the x axis and say what it is, inside the panel."""
    ax.axvspan(a, b, color='#eef1f3', zorder=0)
    if text:
        ax.text(a + .02*(b - a) if x is None else x, y, text, transform=ax.get_xaxis_transform(), ha=ha, va=va,
                fontsize=9.5, color=MUTED, rotation=rotation, zorder=1)


def key(ax, dashed=False, loc=None):
    """Legend for panels whose colours mark series: add what is analytical and what is DES."""
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color=INK, lw=1.2, ls='--' if dashed else '-'),
                Line2D([], [], marker='o', ls='none', ms=4, mfc='white', mec=INK, mew=.9)]
    labels += ['analytical' + (' (dashed)' if dashed else ''), 'DES']
    ax.legend(handles, labels, loc=loc)
