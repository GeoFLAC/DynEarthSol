"""Usage: python excess.py <model prefix>.
Excess pore pressure relative to hydrostatic below the local surface: p - rho_w g (z_s(x) - z), z_s from the
mesh's top boundary nodes (side walls excluded). Writes fig_<model>.png and prints the range per frame."""
import sys, os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.tri as mtri
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from Dynearthsol import Dynearthsol
def excess(d, fr):
    xy = d.read_field(fr, "coordinate"); x, z = xy[:, 0], xy[:, 1]; c = d.read_field(fr, "connectivity")
    e = np.sort(np.vstack([c[:, [0, 1]], c[:, [1, 2]], c[:, [2, 0]]]), axis=1)   # boundary edges: in one element
    u, n = np.unique(e, axis=0, return_counts=True); bn = np.unique(u[n == 1])
    top = bn[(z[bn] > z.min() + 0.5 * (z.max() - z.min())) & (x[bn] > x.min() + 1) & (x[bn] < x.max() - 1)]                       # upper half of the boundary
    top = top[np.argsort(x[top])]  # side walls excluded; np.interp extends the end values
    return xy, (d.read_field(fr, "pore pressure") - 1e4 * (np.interp(x, x[top], z[top]) - z)) / 1e6
run = sys.argv[1]; d = Dynearthsol(run)  # model path prefix, e.g. core-complex-poro
for fr in d.frames[::4] + [d.frames[-1]]:
    xy, ex = excess(d, fr); print("%5.1f kyr  excess %6.2f .. %5.2f MPa" % (d.time[d.frames.index(fr)] / 3.15576e10, ex.min(), ex.max()))
frs = [d.frames[i] for i in (4, 10, -1)]
fig, ax = plt.subplots(3, 2, figsize=(13, 7.5), layout="constrained", sharex=True, sharey=True)
for i, fr in enumerate(frs):
    xy, ex = excess(d, fr); c = d.read_field(fr, "connectivity"); tri = mtri.Triangulation(xy[:, 0] / 1e3, xy[:, 1] / 1e3, c)
    t = d.time[d.frames.index(fr)] / 3.15576e10
    im0 = ax[i, 0].tripcolor(tri, ex, shading="gouraud", cmap="RdBu_r", vmin=-10, vmax=10)
    ax[i, 0].set_title("excess pore pressure (local hydrostatic), %.1f kyr" % t, loc="left", fontsize=10)
    im1 = ax[i, 1].tripcolor(tri, d.read_field(fr, "plastic strain"), cmap="viridis", vmin=0, vmax=2)
    ax[i, 1].set_title("plastic strain, %.1f kyr" % t, loc="left", fontsize=10)
for a in ax.flat: a.set_aspect("equal"); a.set_xlim(25, 75); a.axhline(-3, color="k", lw=0.5, ls="--")
fig.colorbar(im0, ax=ax[:, 0], shrink=0.6, label="MPa"); fig.colorbar(im1, ax=ax[:, 1], shrink=0.6)
fig.savefig("fig_%s.png" % os.path.basename(run), dpi=120)
