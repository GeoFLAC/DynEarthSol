#!/usr/bin/env python3
"""Bishop's simplified method of slices for a slope with a toe on the left.

The ground surface is a piecewise-linear profile z_s(x); the mass slides
toward -x on a circular surface z_b(x) = zc - sqrt(R^2 - (x - xc)^2) that
stays above the base of the domain. With slices of width b, weight W, base
inclination alpha = atan(dz_b/dx) (positive under the crest) and pore
pressure u on the base (the effective-stress pore pressure: chi p, so suction
counts with chi = S_e as in DES),

    FS = sum[(c b + (W - u b) tan(phi)) / m] / sum(W sin(alpha)),
    m = cos(alpha) (1 + tan(alpha) tan(phi)/FS),

solved by fixed-point iteration (Bishop 1955). The critical circle is found
by minimising FS over entry and exit points on the surface and the depth of
the circle (Nelder-Mead from a grid of starts).

    bishop_reference.py          # prints the self-checks and the dry embankment FS
"""
import math

import numpy as np
from scipy.optimize import minimize


class Slope:
    def __init__(self, profile, base, gamma, cohesion, friction, pore=None, slices=400, unit_weight=None):
        """gamma is a constant unit weight, or unit_weight(x, z) a field (integrated over each slice column)."""
        self.px, self.pz = np.array(profile, dtype=float).T
        self.base, self.gamma, self.c, self.tanphi = base, gamma, cohesion, math.tan(friction)
        self.pore = pore or (lambda x, z: np.zeros_like(x))
        self.unit_weight = unit_weight
        self.n = slices

    def surface(self, x):
        return np.interp(x, self.px, self.pz)

    def circle(self, xa, xb, sag):
        """Circle through the surface points at xa < xb whose chord midpoint sinks by sag (> 0)."""
        za, zb = self.surface(xa), self.surface(xb)
        half = math.hypot(xb - xa, zb - za)/2
        R = (half**2 + sag**2)/(2*sag)
        mx, mz = (xa + xb)/2, (za + zb)/2
        nx, nz = -(zb - za)/(2*half), (xb - xa)/(2*half)        # unit normal pointing up
        return mx + nx*(R - sag), mz + nz*(R - sag), R

    def factor(self, xa, xb, sag):
        if not (self.px[0] <= xa < xb <= self.px[-1]) or sag <= 0: return math.inf
        xc, zc, R = self.circle(xa, xb, sag)
        edges = np.linspace(xa, xb, self.n + 1); x = (edges[:-1] + edges[1:])/2; b = edges[1] - edges[0]
        if np.any(np.abs(x - xc) >= R): return math.inf
        root = np.sqrt(R**2 - (x - xc)**2)
        zb = zc - root
        height = self.surface(x) - zb
        if zb.min() < self.base or np.any(height < -1e-9*R): return math.inf
        if self.unit_weight is None:
            W = self.gamma*np.maximum(height, 0)*b
        else:
            g, w = np.polynomial.legendre.leggauss(16)
            z = zb[:, None] + (np.maximum(height, 0)[:, None])*(g[None, :] + 1)/2
            W = b*np.maximum(height, 0)/2*np.sum(w[None, :]*self.unit_weight(np.repeat(x[:, None], 16, 1), z), axis=1)
        alpha = np.arctan((x - xc)/root)
        u = self.pore(x, zb)
        driving = np.sum(W*np.sin(alpha))
        if driving <= 0: return math.inf
        fs = 1.0
        for _ in range(200):
            m = np.cos(alpha)*(1 + np.tan(alpha)*self.tanphi/fs)
            if np.any(m <= 0.05): return math.inf                  # Bishop's m_alpha must stay positive
            new = np.sum((self.c*b + (W - u*b)*self.tanphi)/m)/driving
            if abs(new - fs) < 1e-12*fs: return new
            fs = new
        return fs

    def critical(self, starts):
        """Smallest FS over circles, searched from each start (entry x, exit x, sag)."""
        best = (math.inf, None)
        for s in starts:
            r = minimize(lambda v: self.factor(*v), s, method='Nelder-Mead',
                         options=dict(xatol=1e-4, fatol=1e-9, maxiter=4000))
            if r.fun < best[0]: best = (r.fun, r.x)
        return best


# Classic homogeneous embankment (as in Griffiths and Lane 1999): H = 10 m, 2:1 face, foundation H deep; toe on the left.
H, RUN, TOE, CREST = 10.0, 20.0, 20.0, 20.0
PROFILE = [(0.0, -H), (TOE, -H), (TOE + RUN, 0.0), (TOE + RUN + CREST, 0.0)]
BASE = -2*H
GAMMA, COHESION, FRICTION = 20e3, 10e3, math.radians(20.0)


def embankment(cohesion=COHESION, friction=FRICTION, pore=None, slices=400):
    return Slope(PROFILE, BASE, GAMMA, cohesion, friction, pore, slices)


def starts():
    return [(xa, xb, sag) for xa in (5.0, 12.0, 18.0, 21.0) for xb in (42.0, 48.0, 55.0) for sag in (2.0, 4.0, 7.0)]


if __name__ == '__main__':
    # Self-check 1: a cohesionless slope fails by shallow sliding parallel to the face, FS -> tan(phi)/tan(beta).
    dry_sand = embankment(cohesion=0.0, friction=math.radians(30.0))
    shallow = dry_sand.factor(TOE + 2, TOE + RUN - 2, 0.05)
    print(f'c = 0, phi = 30: shallow circle FS {shallow:.4f}, infinite slope {math.tan(math.radians(30))/0.5:.4f}')
    # Self-check 2: slice convergence of the critical circle.
    fs, v = embankment().critical(starts())
    for n in (100, 400, 1600):
        print(f'{n} slices: FS {embankment(slices=n).factor(*v):.5f}')
    xc, zc, R = embankment().circle(*v)
    const = embankment(); field = Slope(PROFILE, BASE, None, COHESION, FRICTION, unit_weight=lambda x, z: GAMMA + 0*x)
    print(f'unit-weight field = constant: {field.factor(*v) - const.factor(*v):.2e}')
    print(f'critical circle: entry x {v[0]:.2f}, exit x {v[1]:.2f}, centre ({xc:.2f}, {zc:.2f}), R {R:.2f}, FS {fs:.4f}')
