#!/usr/bin/env python3
"""McTigue (1986): a drained half-space heated suddenly at its surface.

The surface z = 0 of a poroelastic half-space z >= 0 (depth) is raised by
dT at t = 0, held drained (p = 0) and traction-free. Heat conducts with
diffusivity kappa, without advection. The strain is uniaxial, so the
vertical total stress vanishes and, with tension-positive stress,

    eps = (alpha p + K alpha_s T)/H,    H = K + 4G/3,

where alpha_s is the volumetric thermal expansion of the solid. With fluid
content zeta = alpha eps + S0 p - beta T and beta = phi alpha_f + (alpha -
phi) alpha_s, Darcy flow gives

    S p_t - (k/mu) p_zz = S gamma T_t,   S = S0 + alpha^2/H,
    gamma = (beta - alpha K alpha_s/H)/S,   c = k/(mu S).

For T = dT erfc(z/(2 sqrt(kappa t))) the solution with p(0, t) = 0 is

    p = gamma dT/(1 - c/kappa) [erfc(z/(2 sqrt(kappa t))) - erfc(z/(2 sqrt(c t)))].

The surface heave is the integral of eps over depth.
McTigue, D. F. (1986), J. Geophys. Res. 91(B9), 9533-9542.
"""
import math

import numpy as np
from scipy.special import erfc


class McTigue:
    def __init__(self, bulk, shear, biot, porosity, fluid_bulk, grain_bulk, permeability, viscosity,
                 alpha_solid, alpha_fluid, kappa, dT):
        self.H = bulk + 4*shear/3
        S0 = porosity/fluid_bulk + (biot - porosity)/grain_bulk
        self.S = S0 + biot**2/self.H
        beta = porosity*alpha_fluid + (biot - porosity)*alpha_solid
        self.gamma = (beta - biot*bulk*alpha_solid/self.H)/self.S
        self.c = permeability/(viscosity*self.S)
        self.kappa, self.dT, self.biot = kappa, dT, biot
        self.thermal = bulk*alpha_solid/self.H        # strain per kelvin at p = 0
        if abs(self.c/kappa - 1) < 1e-6:
            raise ValueError('c and kappa must differ (the c = kappa limit is a separate form)')

    def temperature(self, z, t):
        return self.dT*erfc(np.asarray(z, float)/(2*np.sqrt(self.kappa*t)))

    def pressure(self, z, t):
        z = np.asarray(z, float)
        a, b = erfc(z/(2*np.sqrt(self.kappa*t))), erfc(z/(2*np.sqrt(self.c*t)))
        return self.gamma*self.dT/(1 - self.c/self.kappa)*(a - b)

    def strain(self, z, t):
        """Vertical strain, tension positive."""
        return self.biot*self.pressure(z, t)/self.H + self.thermal*self.temperature(z, t)

    def heave(self, t):
        """Upward surface displacement of the half-space."""
        root = 2/math.sqrt(math.pi)
        thermal = self.thermal*self.dT*root*math.sqrt(self.kappa*t)
        drained = (self.biot/self.H)*self.gamma*self.dT/(1 - self.c/self.kappa) \
            * root*(math.sqrt(self.kappa*t) - math.sqrt(self.c*t))
        return thermal + drained


if __name__ == '__main__':
    m = McTigue(1e9, 1e9, .8, .2, 2.2e9, 5e9, 1e-19, 1e-3, 3e-5, 2.07e-4, 1e-6, 10.0)
    t = 4e4
    z = np.linspace(0, 2, 40001)
    # Boundary and initial conditions.
    assert abs(m.pressure(0.0, t)) < 1e-12*abs(m.gamma*m.dT)
    assert abs(m.pressure(.5, 1e-6)) < 1e-12*abs(m.gamma*m.dT)
    # The pressure solves S p_t - (k/mu) p_zz = S gamma T_t.
    dt = 1.0; zz = z[100:-100]
    pt = (m.pressure(zz, t+dt) - m.pressure(zz, t-dt))/(2*dt)
    Tt = (m.temperature(zz, t+dt) - m.temperature(zz, t-dt))/(2*dt)
    h = z[1]-z[0]
    p = m.pressure(z, t)
    pzz = (p[2:] - 2*p[1:-1] + p[:-2])[99:-99]/h**2
    lhs, rhs = pt - m.c*pzz, m.gamma*Tt
    assert np.allclose(lhs, rhs, rtol=0, atol=1e-5*np.abs(rhs).max()), np.max(np.abs(lhs-rhs))/np.abs(rhs).max()
    # Independent solution: Crank-Nicolson in time, central differences in z, forced by
    # the exact temperature rate, from the closed form at a small start time.
    from scipy.sparse import diags, identity
    from scipy.sparse.linalg import splu
    n, length = 1201, 1.5
    zf = np.linspace(0, length, n); dz = zf[1] - zf[0]
    lap = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n-2, n-2))/dz**2
    k = 0.5
    t0, steps = 200.0, 4000
    pf = m.pressure(zf, t0)[1:-1]; tau = (t - t0)/steps
    A = splu((identity(n-2) - k*tau*m.c*lap).tocsc()); B = (identity(n-2) + k*tau*m.c*lap).tocsr()
    for i in range(steps):
        ta, tb = t0 + i*tau, t0 + (i+1)*tau
        source = m.gamma*(m.temperature(zf[1:-1], tb) - m.temperature(zf[1:-1], ta))
        pf = A.solve(B @ pf + source)
    err = np.max(np.abs(pf - m.pressure(zf[1:-1], t)))/np.abs(m.pressure(zf, t)).max()
    assert err < 2e-4, err
    # The closed-form heave is the integral of the strain.
    assert abs(np.trapezoid(m.strain(z, t), z)/m.heave(t) - 1) < 1e-6
    print(f'mctigue reference: 5 checks pass; gamma = {m.gamma:.4g} Pa/K, c/kappa = {m.c/m.kappa:.3f}')
