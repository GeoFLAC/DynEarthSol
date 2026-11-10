#!/usr/bin/env python3
"""Fluid injection into a poroelastic borehole (plane strain, axisymmetric).

A borehole of radius a in an annulus a <= r <= R is loaded at t = 0+ by
fluid at pressure p_w: the wall carries the pore pressure p_w and the
normal traction -p_w. The outer boundary is traction-free and either
drained (p = 0) or sealed (p_r = 0); the material starts unstressed at
p = 0. In an infinite
medium this is the axisymmetric part of Detournay and Cheng (1988).

With radial symmetry the displacement is u = m J/r + C1(t) r + C2(t)/r,
J(r) = INT_a^r p r' dr', m = alpha/(lambda + 2G), so the volumetric strain
is m p + 2 C1 and the total stresses (tension positive) are

    sigma_rr = -2G m J/r^2 + 2(lambda + G) C1 - 2G C2/r^2,
    sigma_tt = 2G m (J/r^2 - p) + 2(lambda + G) C1 + 2G C2/r^2,

with C1, C2 from sigma_rr(a) = -p_w and sigma_rr(R) = 0; C1 is linear in
p_w and J(R). Fluid mass then gives

    S p_t + 2 alpha C1'(t) = (k/mu) (1/r)(r p_r)_r,     S = 1/M + alpha m.

On a finite annulus C1 is not zero (it vanishes as R -> infinity), so the
pressure couples to its own integral. In Laplace space p_bar is the
homogeneous Bessel solution plus a uniform particular part -2 alpha
C1_bar/S, and C1_bar follows from one scalar equation; the result is
inverted by the Stehfest method. An independent Crank-Nicolson solution
of the same coupled equation checks the inversion.
"""
import math

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import i0e, i1e, k0e, k1e


from laplace import stehfest_weights as _stehfest


def _trapezoid_weights(grid):
    w = np.zeros(len(grid)); h = np.diff(grid)
    w[:-1] += h/2; w[1:] += h/2
    return w*grid                                         # weights of INT p r dr


class Borehole:
    def __init__(self, *, radius, outer, bulk, shear, biot, storage_modulus, mobility, wall_pressure,
                 terms=14, grid=4001, sealed=False):
        self.a, self.R, self.pw, self.alpha, self.sealed = radius, outer, wall_pressure, biot, sealed
        self.G, self.lam = shear, bulk - 2*shear/3
        self.m = biot/(self.lam + 2*shear)
        self.S = 1/storage_modulus + biot*self.m
        self.mobility = mobility
        self.c = mobility/self.S
        self.v = _stehfest(terms)
        self.grid = np.linspace(self.a, self.R, grid)
        self.w = _trapezoid_weights(self.grid)
        # C1 = c_p p_w + c_J J(R) from the two traction conditions.
        A = np.array([[2*(self.lam + shear), -2*shear/self.a**2], [2*(self.lam + shear), -2*shear/self.R**2]])
        self.c_p = np.linalg.solve(A, [-1.0, 0.0])[0]
        self.c_J = np.linalg.solve(A, [0.0, 2*shear*self.m/self.R**2])[0]

    def _shapes(self, r, s):
        """Homogeneous solutions phi1 (1 at a, 0 at R) and phi2 (0 at a, 1 at R); sealed: phi1 with
        phi1'(R) = 0 and phi2 = 0."""
        xi = np.sqrt(s/self.c)
        if self.sealed:
            def neumann(x):     # I0(xi x) K1(xi R) + K0(xi x) I1(xi R), scaled consistently
                return i0e(xi*x)*k1e(xi*self.R)*np.exp(xi*(x - self.R)) + k0e(xi*x)*i1e(xi*self.R)*np.exp(xi*(self.R - x))
            return neumann(r)/neumann(self.a), np.zeros_like(np.asarray(r, float))
        def combo(x, y):        # I0(xi x) K0(xi y) - K0(xi x) I0(xi y), scaled consistently
            return i0e(xi*x)*k0e(xi*y)*np.exp(xi*(x - y)) - k0e(xi*x)*i0e(xi*y)*np.exp(xi*(y - x))
        return combo(r, self.R)/combo(self.a, self.R), combo(r, self.a)/combo(self.R, self.a)

    def _pbar(self, r, s):
        g1, g2 = self._shapes(self.grid, s)
        I1, I3 = self.w @ g1, self.w @ (1 - g1 - g2)
        C1 = (self.pw/s)*(self.c_p + self.c_J*I1)/(1 + 2*self.alpha*self.c_J*I3/self.S)
        q = -2*self.alpha*C1/self.S
        f1, f2 = self._shapes(np.asarray(r, float), s)
        return (self.pw/s)*f1 + q*(1 - f1 - f2)

    def pressure(self, r, t):
        ln2 = math.log(2)/t
        return ln2*sum(vk*self._pbar(r, (k + 1)*ln2) for k, vk in enumerate(self.v))

    def fields(self, r, t):
        """Pressure, radial displacement and total radial/hoop stress at radii r."""
        grid = self.grid
        p = self.pressure(grid, t)
        J = np.concatenate([[0.0], np.cumsum(np.diff(grid)*.5*(p[1:]*grid[1:] + p[:-1]*grid[:-1]))])
        G, lam, m = self.G, self.lam, self.m
        C1 = self.c_p*self.pw + self.c_J*J[-1]
        C2 = (2*(lam + G)*C1 + self.pw)*self.a**2/(2*G)
        r = np.asarray(r, float)
        Jr, pr = np.interp(r, grid, J), np.interp(r, grid, p)
        return dict(p=pr, u=m*Jr/r + C1*r + C2/r,
                    srr=-2*G*m*Jr/r**2 + 2*(lam + G)*C1 - 2*G*C2/r**2,
                    stt=2*G*m*(Jr/r**2 - pr) + 2*(lam + G)*C1 + 2*G*C2/r**2)


def crank_nicolson(b, t, n=2000, steps=4000):
    """Independent finite-volume solution of S p_t + 2 alpha C1' = (k/mu)(1/r)(r p_r)_r,
    C1 = c_p p_w + c_J J(R); the last node is unknown when the outer boundary is sealed."""
    r = np.linspace(b.a, b.R, n + 1); h = r[1] - r[0]
    w = _trapezoid_weights(r)
    p = np.zeros(n + 1); p[0] = b.pw
    dt = t/steps
    last = n + 1 if b.sealed else n                    # unknowns are p[1:last]
    k = slice(1, last); m = last - 1
    vol = r[k]*h; vol = vol.copy()
    left, right = b.c*(r[k] - h/2)/h, b.c*(r[k] + h/2)/h
    if b.sealed:
        vol[-1] = r[-1]*h/2; right[-1] = 0.0           # half cell, no flux through R
    lo, up = left/vol, right/vol
    di = -(lo + up)
    ab = np.zeros((3, m))
    ab[0, 1:] = -.5*dt*up[:-1]; ab[1] = 1 - .5*dt*di; ab[2, :-1] = -.5*dt*lo[1:]
    gamma = 2*b.alpha*b.c_J/b.S                        # coupling to J(R): rank-one term
    z = solve_banded((1, 1), ab, np.full(m, gamma))
    wi = w[k]
    C1_old = 0.0                                       # unloaded before t = 0
    for _ in range(steps):
        cur = p[k]
        nxt = np.append(p[2:last], 0.0) if b.sealed else p[2:last + 1]
        rhs = cur + .5*dt*(lo*p[0:last - 1] + di*cur + up*nxt)
        rhs[0] += .5*dt*lo[0]*b.pw
        rhs -= 2*b.alpha*(b.c_p*b.pw + b.c_J*w[0]*b.pw - C1_old)/b.S
        y = solve_banded((1, 1), ab, rhs)
        p[k] = y - z*(wi @ y)/(1 + wi @ z)
        C1_old = b.c_p*b.pw + b.c_J*(w @ p)
    return r, p


if __name__ == '__main__':
    b = Borehole(radius=1.0, outer=20.0, bulk=2e8, shear=1e8, biot=1.0, storage_modulus=1e9,
                 mobility=1e-12/1e-3, wall_pressure=1e6)
    t = 9/b.c                                           # diffusion length 3 a
    r = np.linspace(1.0, 20.0, 200)
    f = b.fields(r, t)
    # Wall and outer conditions.
    assert abs(f['p'][0] - b.pw) < 1e-6*b.pw and abs(f['srr'][0] + b.pw) < 1e-9*b.pw and abs(f['srr'][-1]) < 1e-9*b.pw
    # Radial equilibrium d(sigma_rr)/dr + (sigma_rr - sigma_tt)/r = 0 (one-sided ends excluded).
    rr = np.linspace(1.5, 15, 4001); g = b.fields(rr, t)
    res = (np.gradient(g['srr'], rr) + (g['srr'] - g['stt'])/rr)[1:-1]
    assert np.max(np.abs(res)) < 2e-4*b.pw, np.max(np.abs(res))/b.pw
    # Stehfest inversion of the coupled problem against an independent Crank-Nicolson solution.
    fr, fp = crank_nicolson(b, t)
    err = np.max(np.abs(b.pressure(fr[1:-1:20], t) - fp[1:-1:20]))/b.pw
    assert err < 2e-4, err
    # The finite-annulus coupling matters at this accuracy: the uncoupled pressure differs.
    coupled = b.pressure(r, t); b_unc = Borehole(radius=1.0, outer=20.0, bulk=2e8, shear=1e8, biot=1.0,
                                                 storage_modulus=1e9, mobility=1e-9, wall_pressure=1e6)
    b_unc.c_J = b_unc.c_p = 0.0
    shift = np.max(np.abs(coupled - b_unc.pressure(r, t)))/b.pw
    # Late time: steady logarithmic profile (C1 constant, no coupling left).
    late = b.pressure(r, 1e4*b.R**2/b.c)
    steady = b.pw*np.log(b.R/r)/np.log(b.R/b.a)
    assert np.max(np.abs(late - steady)) < 1e-4*b.pw
    # Sealed outer boundary: the same checks against its own Crank-Nicolson solution, and p -> p_w late.
    bs = Borehole(radius=1.0, outer=20.0, bulk=2e8, shear=1e8, biot=1.0, storage_modulus=1e9,
                  mobility=1e-9, wall_pressure=1e6, sealed=True)
    fr, fp = crank_nicolson(bs, t)
    err_s = np.max(np.abs(bs.pressure(fr[1:-1:20], t) - fp[1:-1:20]))/bs.pw
    assert err_s < 2e-4, err_s
    gs = bs.fields(rr, t)
    res = (np.gradient(gs['srr'], rr) + (gs['srr'] - gs['stt'])/rr)[1:-1]
    assert np.max(np.abs(res)) < 2e-4*bs.pw and abs(bs.fields(np.array([20.0]), t)['srr'][0]) < 1e-9*bs.pw
    assert np.max(np.abs(bs.pressure(r, 1e3*bs.R**2/bs.c) - bs.pw)) < 1e-4*bs.pw
    print(f'borehole reference: 8 checks pass; c = {b.c:.3g} m2/s; Stehfest vs Crank-Nicolson {err:.1e} (drained), '
          f'{err_s:.1e} (sealed) of p_w; finite-annulus coupling shifts p by {shift:.1e} of p_w')
