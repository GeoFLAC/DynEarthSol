#!/usr/bin/env python3
"""Rudnicki (1986): a continuous fluid point source in an infinite poroelastic solid.

Water is injected at a constant volume rate Q from t = 0 at the origin. The
displacement is irrotational and vanishes at infinity, so (lambda + 2G) eps =
alpha p, the fluid content is S p with S = 1/M + alpha^2/(lambda+2G), and

    p(r,t) = Q mu/(4 pi kappa r) erfc(r/a),   a = 2 sqrt(c t),   c = kappa/(mu S).

The radial displacement u = (1/r^2) int_0^r eps s^2 ds integrates in closed form:

    u(r,t) = A/r^2 [ r^2/2 erfc(r/a) + a^2/4 erf(r/a) - a r/(2 sqrt(pi)) exp(-r^2/a^2) ],
    A = alpha Q mu/(4 pi kappa (lambda+2G)),

with strains eps_rr = du/dr = eps - 2u/r and eps_tt = u/r, and total stresses
(tension positive) sigma_rr = -4G u/r and sigma_tt = 2G(u/r - eps). These are
the stress changes a nearby fault feels from injection.
Rudnicki, J. W. (1986), Mech. Mater. 5, 383-393.
"""
import math

import numpy as np
from scipy.special import erf, erfc


class Rudnicki:
    def __init__(self, bulk, shear, biot, porosity, fluid_bulk, grain_bulk, permeability, viscosity, rate):
        self.G, self.alpha, self.Q = shear, biot, rate
        self.H = bulk + 4*shear/3                                   # lambda + 2G
        inv_M = porosity/fluid_bulk + (biot-porosity)/grain_bulk
        self.S = inv_M + biot**2/self.H
        self.mobility = permeability/viscosity
        self.c = self.mobility/self.S
        self.A = biot*rate/(4*math.pi*self.mobility*self.H)

    def pressure(self, r, t):
        r, t = np.broadcast_arrays(np.asarray(r, float), np.asarray(t, float))
        a = 2*np.sqrt(self.c*t)
        return self.Q/(4*math.pi*self.mobility*r)*erfc(r/a)

    def displacement(self, r, t):
        r, t = np.broadcast_arrays(np.asarray(r, float), np.asarray(t, float))
        a = 2*np.sqrt(self.c*t); x = r/a
        bracket = r*r/2*erfc(x) + a*a/4*erf(x) - a*r/(2*math.sqrt(math.pi))*np.exp(-x*x)
        return self.A*bracket/(r*r)

    def strains(self, r, t):
        """(volumetric, radial, hoop) strain."""
        eps = self.alpha*self.pressure(r, t)/self.H
        hoop = self.displacement(r, t)/np.asarray(r, float)
        return eps, eps - 2*hoop, hoop

    def stresses(self, r, t):
        """Total (radial, hoop) stress changes, tension positive."""
        eps, _, hoop = self.strains(r, t)
        return -4*self.G*hoop, 2*self.G*(hoop - eps)


if __name__ == '__main__':
    m = Rudnicki(2e8, 1e8, .8, .25, 2.2e9, 1e9, 1e-13, 1e-3, 1e-4)
    r = np.linspace(.5, 20, 400); t = 120.0
    # The closed-form displacement reproduces the defining integral.
    s = np.linspace(1e-9, r[-1], 400001)
    eps_s = m.alpha*m.pressure(s, t)/m.H
    integral = np.concatenate(([0], np.cumsum(.5*(eps_s[1:]*s[1:]**2 + eps_s[:-1]*s[:-1]**2)*np.diff(s))))
    u_num = np.interp(r, s, integral)/r**2
    assert np.allclose(m.displacement(r, t), u_num, rtol=1e-6, atol=0), np.max(np.abs(m.displacement(r, t)/u_num-1))
    # Strains add up to the volumetric strain; the radial strain is du/dr.
    eps, err, ett = m.strains(r, t)
    assert np.allclose(err + 2*ett, eps, rtol=1e-12)
    du = np.gradient(m.displacement(r, t), r, edge_order=2)
    assert np.allclose(du[2:-2], err[2:-2], rtol=1e-4, atol=1e-12)
    # Pressure solves S p_t = (kappa/mu) lap p away from the source.
    dt = 1e-3; rr = np.linspace(2, 10, 801)
    pt = (m.pressure(rr, t+dt) - m.pressure(rr, t-dt))/(2*dt)
    rp = rr*m.pressure(rr, t)
    lap = np.gradient(np.gradient(rp, rr, edge_order=2), rr, edge_order=2)/rr
    assert np.allclose(m.S*pt[5:-5], m.mobility*lap[5:-5], rtol=2e-4)
    # Late time: steady pressure Q/(4 pi kappa r), uniform displacement A/2.
    assert abs(m.displacement(5.0, 1e12)/(m.A/2) - 1) < 1e-4
    print(f'rudnicki reference: 4 checks pass; c = {m.c:.4g} m2/s, S = {m.S:.4g} 1/Pa')
