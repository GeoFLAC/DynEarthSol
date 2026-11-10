#!/usr/bin/env python3
"""Pressure-driven Mohr-Coulomb yielding of a laterally confined layer.

A saturated layer 0 <= d <= H (d depth below the top) carries the DES
legacy initial stress, isotropic and lithostatic, with hydrostatic pore
pressure, and a drained surcharge q on its top. Under uniaxial strain both
horizontal stresses stay equal, kappa = nu/(1 - nu) times the vertical
increment, so the effective stresses (tension positive) are

    sigma'_v = -(gamma' d + q) + alpha dp,
    sigma'_h = -(gamma' d + kappa q) + kappa alpha dp,     gamma' = (rho_b - alpha rho_f) g,

when the pore pressure rises by dp everywhere. The surcharge keeps the top
in compression: without it the top reaches the Mohr-Coulomb apex (mean
effective tension c cot(phi)) just after yielding and loses all strength.
Yield, (s1 - s3)/2 + (s1 + s3)/2 sin(phi) = c cos(phi) on the principal
stresses, starts at the first root dp*(d) of this trial function; dp*
increases with depth, so a yield front moves down as the pressure rises.
Each depth is independent under uniaxial strain, so the front stays exact
after the layer above has yielded. The yielded top reaches the apex when
sigma'_v = c cot(phi), which bounds the useful pressure range.
"""
import math

import numpy as np
from scipy.optimize import brentq


class PressureYield:
    def __init__(self, *, bulk, shear, biot, cohesion, friction_deg, rho_bulk, rho_fluid, surcharge=0.0, g=9.81):
        self.alpha, self.c, self.q = biot, cohesion, surcharge
        self.phi = math.radians(friction_deg)
        nu = (3*bulk - 2*shear)/(2*(3*bulk + shear))
        self.kappa = nu/(1 - nu)
        self.gamma_eff = (rho_bulk - biot*rho_fluid)*g

    def stresses(self, d, dp):
        d, dp = np.broadcast_arrays(np.asarray(d, float), np.asarray(dp, float))
        sv = -(self.gamma_eff*d + self.q) + self.alpha*dp
        sh = -(self.gamma_eff*d + self.kappa*self.q) + self.kappa*self.alpha*dp
        return sv, sh

    def yield_function(self, d, dp):
        """Mohr-Coulomb function of the elastic trial state (negative inside)."""
        sv, sh = self.stresses(d, dp)
        s1, s3 = np.maximum(sv, sh), np.minimum(sv, sh)
        return (s1 - s3)/2 + (s1 + s3)/2*math.sin(self.phi) - self.c*math.cos(self.phi)

    def onset(self, d, top=1e7):
        """First pressure rise at which depth d yields (vectorised over d)."""
        out = []
        for dd in np.atleast_1d(np.asarray(d, float)):
            if self.yield_function(dd, 0.0) >= 0:
                raise ValueError(f'depth {dd} is already yielding under the surcharge')
            grid = np.linspace(0, top, 20001)
            f = self.yield_function(dd, grid)
            k = int(np.argmax(f > 0))
            out.append(brentq(lambda x: float(self.yield_function(dd, x)), grid[k-1], grid[k], xtol=1e-9, rtol=1e-14))
        return np.array(out) if np.ndim(d) else out[0]

    def apex(self, d):
        """Pressure rise at which the yielded state at depth d reaches the apex."""
        return (self.c/math.tan(self.phi) + self.gamma_eff*np.asarray(d, float) + self.q)/self.alpha


if __name__ == '__main__':
    ref = PressureYield(bulk=1e8, shear=6e7, biot=1.0, cohesion=1e4, friction_deg=30.0, rho_bulk=2155.0,
                        rho_fluid=1000.0, surcharge=5e4)
    d = np.linspace(0, 10, 41)
    on = ref.onset(d)
    # The onset lies on the yield surface and the trial state is elastic just before it.
    assert np.allclose(ref.yield_function(d, on), 0, atol=1e-6*ref.c)
    assert np.all(ref.yield_function(d, .999*on) < 0)
    # Without a surcharge the closed form of the vertical-major branch is recovered.
    free = PressureYield(bulk=1e8, shear=6e7, biot=1.0, cohesion=1e4, friction_deg=30.0, rho_bulk=2155.0, rho_fluid=1000.0)
    k, sp = free.kappa, math.sin(free.phi)
    closed = 2*(free.c*math.cos(free.phi) + free.gamma_eff*d*sp)/((1 - k) + (1 + k)*sp)
    assert np.allclose(free.onset(d), closed, rtol=1e-10)
    # The front moves down and the apex stays ahead of the onset.
    assert np.all(np.diff(on) > 0) and np.all(ref.apex(d) > on)
    print(f'pressure-yield reference: 4 checks pass; onset {on[0]:.0f} Pa at the top, {on[-1]:.0f} Pa at 10 m; '
          f'apex at the top {ref.apex(0.0):.0f} Pa')
