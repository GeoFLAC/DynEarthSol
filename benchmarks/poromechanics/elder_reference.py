#!/usr/bin/env python3
"""Elder's heater problem, transient, independent of DES.

Elder (1967, J. Fluid Mech. 27, 609-623) heated a porous layer of height H
from below over the central half of a base of width 4H and kept the top
cold. The saline version of Voss and Souza (1987, Water Resour. Res. 23,
1851-1866) is the same problem with salt sinking from the top; there one year
is tau = 0.005. This solver takes the symmetric half, 0 <= x <= 2H, in the
nondimensional Darcy-Boussinesq form with u = psi_z, w = -psi_x:

    lap psi = -Ra theta_x,    theta_t + u theta_x + w theta_z = lap theta,

with tau = kappa t/(sigma H^2). Boundaries are impermeable (psi = 0);
theta = 0 at the top; theta = 1 on the heater, 0 <= x <= 1 at z = 0; the
rest of the base, the symmetry axis x = 0 and the far wall x = 2 are
insulated. The fluid starts at theta = 0. Second-order central differences
on a uniform grid with h = 1/n, a direct Poisson solve, and SSP-RK3 steps
limited by diffusion and advection.

The half domain imposes mirror symmetry. Van Reeuwijk et al. (2009, Water
Resour. Res. 45, W04416) found several steady states at Ra = 400, so late
times depend on the branch; the comparison window ends at tau = 0.1.

    elder_reference.py OUTPUT.npz [--ra 400] [--grids 32 64 128 256]
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import splu

TIMES = (0.005, 0.01, 0.02, 0.05, 0.075, 0.1)   # Voss and Souza's 1, 2, 4, 10, 15 and 20 years


class Grid:
    def __init__(self, n):
        self.n, self.h = n, 1.0/n
        self.x = np.linspace(0, 2, 2*n + 1)
        self.z = np.linspace(0, 1, n + 1)
        self.fixed = np.zeros((n + 1, 2*n + 1), bool)
        self.fixed[-1, :] = True                         # cold top
        self.fixed[0, :n + 1] = True                     # heater, x <= 1
        self.value = np.zeros_like(self.fixed, float)
        self.value[0, :n + 1] = 1.0
        k, m = n - 1, 2*n - 1                            # interior psi unknowns
        lz = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(k, k))
        lx = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(m, m))
        self.poisson = splu(((kron(lz, identity(m)) + kron(identity(k), lx))/self.h**2).tocsc())

    def derivatives(self, theta):
        """Laplacian and gradient with insulated mirrors; fixed nodes are overwritten later."""
        h, p = self.h, np.pad(theta, 1, mode='reflect')
        lap = (p[1:-1, 2:] + p[1:-1, :-2] + p[2:, 1:-1] + p[:-2, 1:-1] - 4*theta)/h**2
        tx = (p[1:-1, 2:] - p[1:-1, :-2])/(2*h)
        tz = (p[2:, 1:-1] - p[:-2, 1:-1])/(2*h)
        return lap, tx, tz

    def velocity(self, tx, ra):
        psi = np.zeros_like(tx)
        rhs = -ra*tx[1:-1, 1:-1]
        psi[1:-1, 1:-1] = self.poisson.solve(rhs.ravel()).reshape(rhs.shape)
        u = np.gradient(psi, self.h, axis=0, edge_order=2)
        w = -np.gradient(psi, self.h, axis=1, edge_order=2)
        return psi, u, w

    def rate(self, theta, ra):
        lap, tx, tz = self.derivatives(theta)
        _, u, w = self.velocity(tx, ra)
        f = lap - u*tx - w*tz
        f[self.fixed] = 0.0
        return f, max(np.abs(u).max(), np.abs(w).max())

    def diagnostics(self, theta):
        h = self.h
        wx = np.full(len(self.x), h); wx[[0, -1]] = h/2
        wz = np.full(len(self.z), h); wz[[0, -1]] = h/2
        top = -(3*theta[-1] - 4*theta[-2] + theta[-3])/(2*h)
        base = -(-3*theta[0] + 4*theta[1] - theta[2])/(2*h)
        heater = base[:self.n + 1]
        wh = wx[:self.n + 1].copy(); wh[-1] = h/2
        return dict(heat=float(wz @ theta @ wx), nu_top=float(top @ wx), nu_heater=float(heater @ wh))


def solve(n, ra=400.0, times=TIMES, cfl=0.4):
    """theta snapshots at the requested tau and the heat-budget history."""
    g = Grid(n)
    theta = np.where(g.fixed, g.value, 0.0)
    t, out, history, steps = 0.0, [], [], 0
    for target in times:
        while t < target - 1e-14:
            k1, speed = g.rate(theta, ra)
            dt = min(cfl*min(g.h**2/4, g.h/max(speed, 1e-30)), target - t)
            a = theta + dt*k1
            k2, _ = g.rate(a, ra)
            b = 0.75*theta + 0.25*(a + dt*k2)
            k3, _ = g.rate(b, ra)
            theta = theta/3 + 2/3*(b + dt*k3)
            theta[g.fixed] = g.value[g.fixed]
            t += dt; steps += 1
            if steps % 200 == 0:
                history.append(dict(tau=t, **g.diagnostics(theta)))
        assert np.isfinite(theta).all()
        out.append(theta.copy())
        history.append(dict(tau=t, **g.diagnostics(theta)))
    return g, np.array(out), history, steps


def study(output, ra, grids):
    data, summary = {}, dict(ra=ra, times=list(TIMES), grids=list(grids), cases={})
    for n in grids:
        g, snaps, history, steps = solve(n, ra)
        data[f'theta_{n}'] = snaps
        data[f'x_{n}'], data[f'z_{n}'] = g.x, g.z
        at = [h for h in history if any(abs(h['tau'] - t) < 1e-12 for t in TIMES)]
        summary['cases'][str(n)] = dict(steps=steps, at_times=at)
        print(json.dumps(dict(n=n, steps=steps, heat=[round(a['heat'], 6) for a in at])), flush=True)
    # Grid convergence: successive differences of each snapshot on the coarse grid's nodes.
    conv = []
    for a, b in zip(grids, grids[1:]):
        r = b//a
        diff = np.abs(data[f'theta_{b}'][:, :, ::r][:, ::r, :] - data[f'theta_{a}'])
        conv.append(dict(coarse=a, fine=b, max_difference=[float(d.max()) for d in diff]))
    summary['successive_max_differences'] = conv
    np.savez_compressed(output, times=np.array(TIMES), ra=ra, **data)
    Path(output).with_suffix('.json').write_text(json.dumps(summary, indent=1))
    return summary


def self_checks():
    # Without buoyancy the base heats the layer by conduction only, with the heat
    # budget d(heat)/dtau = Nu_heater - Nu_top up to the heater-edge singularity.
    _, cond, history, _ = solve(32, ra=0.0, times=(0.04, 0.05))
    h0, h1 = history[-2], history[-1]
    rate = (h1['heat'] - h0['heat'])/(h1['tau'] - h0['tau'])
    mean_flux = 0.5*(h0['nu_heater'] - h0['nu_top'] + h1['nu_heater'] - h1['nu_top'])
    assert abs(rate - mean_flux) < 0.1*abs(mean_flux), (rate, mean_flux)
    # With Ra = 400 plumes carry heat to mid-height above the heater, well
    # beyond conduction, and the layer takes up more heat.
    _, conv, history400, _ = solve(32, times=(0.04, 0.05))
    mid, mid0 = conv[1][16, :33].max(), cond[1][16, :33].max()
    assert mid > 0.2 and mid > 3*mid0, (mid, mid0)
    assert history400[-1]['heat'] > h1['heat']
    print('elder reference: self checks pass', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('--ra', type=float, default=400.0)
    p.add_argument('--grids', type=int, nargs='+', default=[32, 64, 128, 256])
    a = p.parse_args()
    if any(b % a_ for a_, b in zip(a.grids, a.grids[1:])):
        p.error('each grid must refine the previous one by an integer factor')
    self_checks()
    study(a.output, a.ra, a.grids)
