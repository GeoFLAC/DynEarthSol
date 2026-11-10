#!/usr/bin/env python3
"""Mandel's problem: plane-strain consolidation under rigid frictionless plates.

A strip -a < x < a, -b < z < b is squeezed by rigid, impermeable, frictionless
plates carrying a total compressive force 2F per unit length; the faces
x = +-a are drained and traction-free. Cheng and Detournay (1988) and
Abousleiman et al. (1996) give, with alpha_n the positive roots of
tan(alpha) = alpha (1 - nu)/(nu_u - nu):

    p(x,t) = 2 F B (1 + nu_u)/(3a) sum_n sin a_n/(a_n - sin a_n cos a_n)
             (cos(a_n x/a) - cos a_n) exp(-a_n^2 c t/a^2)
    u_z(b,t)/b = -F(1-nu)/(2Ga) + F(1-nu_u)/(Ga) sum_n sin a_n cos a_n/(a_n - sin a_n cos a_n)
             exp(-a_n^2 c t/a^2)

c = 2 kappa B^2 G (1-nu)(1+nu_u)^2 / (9 mu (1-nu_u)(nu_u-nu)). At the centre
the pressure first rises above its undrained value (Mandel-Cryer effect).
"""
import math

import numpy as np
from scipy.optimize import brentq


class Mandel:
    def __init__(self, bulk, shear, biot, porosity, fluid_bulk, grain_bulk, permeability, viscosity,
                 half_width, force):
        self.a, self.F, self.G = half_width, force, shear
        inv_M = porosity/fluid_bulk + (biot-porosity)/grain_bulk
        self.M = 1/inv_M
        Ku = bulk + biot**2*self.M
        self.B = biot*self.M/Ku
        self.nu = (3*bulk-2*shear)/(2*(3*bulk+shear))
        self.nu_u = (3*Ku-2*shear)/(2*(3*Ku+shear))
        nu, nuu, B = self.nu, self.nu_u, self.B
        self.c = 2*permeability*B**2*shear*(1-nu)*(1+nuu)**2/(9*viscosity*(1-nuu)*(nuu-nu))
        self.p0 = force*B*(1+nuu)/(3*half_width)   # undrained pressure at t = 0+
        k = (1-nu)/(nuu-nu)
        f = lambda x: math.sin(x) - k*x*math.cos(x)   # tan x = k x, without poles
        roots = []
        for n in range(1, 400):
            lo, hi = (n-1)*math.pi + 1e-12, (n-.5)*math.pi - 1e-12
            if f(lo)*f(hi) > 0:
                lo, hi = (n-.5)*math.pi + 1e-12, n*math.pi - 1e-12
            roots.append(brentq(f, lo, hi, xtol=1e-15, maxiter=200))
        self.roots = np.array(roots)

    def pressure(self, x, t):
        x, t = np.broadcast_arrays(np.asarray(x, float), np.asarray(t, float))
        r = self.roots[:, None]
        s, c = np.sin(r), np.cos(r)
        terms = s/(r-s*c)*(np.cos(r*x.ravel()/self.a)-c)*np.exp(-r**2*self.c*t.ravel()/self.a**2)
        return (2*self.F*self.B*(1+self.nu_u)/(3*self.a)*terms.sum(0)).reshape(x.shape)

    def settlement(self, height, t):
        """Downward plate displacement for a half height b = height."""
        t = np.asarray(t, float)
        r = self.roots[:, None]
        s, c = np.sin(r), np.cos(r)
        series = (s*c/(r-s*c)*np.exp(-r**2*self.c*t.ravel()/self.a**2)).sum(0)
        uz = -self.F*(1-self.nu)/(2*self.G*self.a) + self.F*(1-self.nu_u)/(self.G*self.a)*series
        return (-uz*height).reshape(t.shape)


if __name__ == '__main__':
    m = Mandel(1e7, 1e7, 1.0, .3, 2.17e9, 37e9, 1.03e-12, 1.002e-3, 1.0, 1e5)
    # Root equation and undrained limit.
    k = (1-m.nu)/(m.nu_u-m.nu)
    assert np.allclose(np.tan(m.roots[:50]), k*m.roots[:50], rtol=1e-8)
    assert abs(m.pressure(0.3, 0.0)/m.p0 - 1) < 2e-3          # uniform undrained pressure
    assert abs(m.pressure(m.a, 5.0)) < 1e-9*m.p0              # drained edge
    # Mandel-Cryer effect: the centre pressure first rises, then decays.
    tc = m.a**2/m.c
    centre = m.pressure(np.zeros(200), np.linspace(1e-4, .3, 200)*tc)
    assert centre.max() > 1.01*m.p0 and centre[-1] < centre.max()
    # Settlement: undrained F(1-nu_u)/(2Ga)*b, drained F(1-nu)/(2Ga)*b.
    assert abs(m.settlement(1, 0)/(m.F*(1-m.nu_u)/(2*m.G*m.a)) - 1) < 2e-3
    assert abs(m.settlement(1, 50*tc)/(m.F*(1-m.nu)/(2*m.G*m.a)) - 1) < 1e-9
    print(f'mandel reference: 5 checks pass; c = {m.c:.4g} m2/s, p0 = {m.p0:.4g} Pa, '
          f'peak/p0 = {centre.max()/m.p0:.4f}')
