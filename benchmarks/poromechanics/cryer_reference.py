#!/usr/bin/env python3
"""Cryer's problem: a poroelastic sphere under a sudden external pressure.

A sphere of radius R, drained at its surface (p = 0), receives a uniform
compressive normal traction P at t = 0+. Spherical symmetry makes the strain
irrotational, so (lambda + 2G) eps = alpha p + f(t), and with
S = 1/M + alpha^2/(lambda+2G), c = kappa/(mu S) the pressure obeys

    p_t - c lap p = -beta d<p>/dt,     beta = 4 G alpha^2 / (3 K S (lambda+2G)),

where <p> is the volume average (from sigma_rr(R) = -P). Its undrained value
is p0 = alpha P / (K S (1 + beta)). With q = sqrt(s/c), x = qR and
g(x) = 3 (x coth x - 1)/x^2, the Laplace transform is

    p(r,s)/p0 = (1+beta) (1 - R sinh(q r)/(r sinh x)) / (s (1 + beta (1 - g(x)))).

Its poles s_n = -c xi_n^2/R^2 solve 3 beta (sin xi - xi cos xi) = (1+beta) xi^2 sin xi,
and the residues give the series used here. The centre pressure first rises
above p0 (Mandel-Cryer effect). Cryer (1963); Detournay and Cheng (1993).
The series is checked against a Stehfest inversion of the transform.
"""
import math

import numpy as np
from scipy.optimize import brentq


class Cryer:
    def __init__(self, bulk, shear, biot, porosity, fluid_bulk, grain_bulk, permeability, viscosity,
                 radius, pressure):
        self.R, self.P = radius, pressure
        M = 1/(porosity/fluid_bulk + (biot-porosity)/grain_bulk)
        H = bulk + 4*shear/3                         # lambda + 2G
        self.S = 1/M + biot**2/H
        self.c = permeability/(viscosity*self.S)
        self.beta = 4*shear*biot**2/(3*bulk*self.S*H)
        self.p0 = biot*pressure/(bulk*self.S*(1+self.beta))
        b = self.beta
        f = lambda x: 3*b*(math.sin(x) - x*math.cos(x)) - (1+b)*x*x*math.sin(x)
        grid = np.linspace(1e-3, 400*math.pi, 400*400)
        values = np.array([f(x) for x in grid])
        roots = [brentq(f, grid[i], grid[i+1], xtol=1e-14) for i in np.nonzero(np.sign(values[:-1]) != np.sign(values[1:]))[0]]
        self.xi = np.array(roots)
        x = self.xi
        # s dD/ds at the poles, D(s) = 1 + beta (1 - g(x)).
        sdd = -(b/2)*(3*(x/np.sin(x)**2 - 1/np.tan(x))/x + 6*(x/np.tan(x) - 1)/x**2)
        self.weight = (1+b)/sdd

    def pressure(self, r, t):
        """Pressure p(r,t) in Pa; r = 0 is the centre."""
        r, t = np.broadcast_arrays(np.asarray(r, float), np.asarray(t, float))
        rr, tt = r.ravel()[None, :], t.ravel()[None, :]
        x = self.xi[:, None]
        with np.errstate(invalid='ignore', divide='ignore'):
            shape = np.where(rr > 1e-12*self.R, self.R*np.sin(x*rr/self.R)/(rr*np.sin(x)), x/np.sin(x))
        terms = self.weight[:, None]*(1-shape)*np.exp(-x**2*self.c*tt/self.R**2)
        return (self.p0*terms.sum(0)).reshape(r.shape)

    def transform(self, r, s):
        """Laplace transform of p(r,t)/p0, for the Stehfest check."""
        q = math.sqrt(s/self.c); x = q*self.R
        g = 3*(x/math.tanh(x) - 1)/x**2
        # sinh(qr)/sinh(x) without overflow for large s.
        ratio = lambda y: math.exp(y-x)*(1-math.exp(-2*y))/(1-math.exp(-2*x))
        shape = 2*x*math.exp(-x)/(1-math.exp(-2*x)) if r == 0 else self.R*ratio(q*r)/r
        return (1+self.beta)*(1-shape)/(s*(1+self.beta*(1-g)))


def stehfest(F, t, n=14):
    """Gaver-Stehfest inversion of a smooth, non-oscillating transform."""
    total, half = 0.0, n//2
    for k in range(1, n+1):
        v = sum(j**half*math.factorial(2*j)/(math.factorial(half-j)*math.factorial(j)*math.factorial(j-1)
                * math.factorial(k-j)*math.factorial(2*j-k)) for j in range((k+1)//2, min(k, half)+1))
        total += (-1)**(k+half)*v*F(k*math.log(2)/t)
    return total*math.log(2)/t


if __name__ == '__main__':
    m = Cryer(1e7, 1e7, 1.0, .3, 2.17e9, 37e9, 1.03e-12, 1.002e-3, 1.0, 1e5)
    tc = m.R**2/m.c
    b = m.beta
    root = lambda x: 3*b*(math.sin(x) - x*math.cos(x)) - (1+b)*x*x*math.sin(x)
    assert all(abs(root(x)) < 1e-8*(1+x*x) for x in m.xi[:40])
    # Series against an independent Laplace inversion, at the centre and mid-radius.
    worst = 0
    for t in (.005, .02, .05, .1, .3, .6):
        for r in (0, .5*m.R):
            series = m.pressure(r, t*tc)/m.p0
            numeric = stehfest(lambda s: m.transform(r, s), t*tc)
            worst = max(worst, abs(series-numeric))
    assert worst < 2e-4, worst
    # Undrained start: the interior excess over p0 grows like sqrt(t) as the shell drains.
    early = [m.pressure(.5*m.R, t*tc)/m.p0 - 1 for t in (1e-5, 1e-4)]
    assert 0 < early[0] < 1e-2 and abs(early[1]/early[0]/math.sqrt(10) - 1) < 2e-2, early
    assert abs(m.pressure(m.R, .1*tc)) < 1e-9*m.p0                  # drained surface
    centre = m.pressure(np.zeros(300), np.linspace(1e-4, .4, 300)*tc)/m.p0
    assert centre.max() > 1.01 and centre[-1] < centre.max()       # Mandel-Cryer effect
    print(f'cryer reference: 5 checks pass; c = {m.c:.4g} m2/s, p0 = {m.p0:.4g} Pa, '
          f'peak/p0 = {centre.max():.4f}, series vs Stehfest {worst:.1e}')
