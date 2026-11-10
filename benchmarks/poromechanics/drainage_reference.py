#!/usr/bin/env python3
"""Coupled gravity drainage of a poroelastic column (unsaturated level U9).

A laterally confined column 0 <= z <= H (z up) is saturated and hydrostatic
with its water table at the top, p = rho_f g (H - z). At t = 0 its base is
opened to p = 0; the top is sealed and traction-free, the base fixed. Water
drains through the base, the upper part desaturates, and the column settles
under the Bishop effective stress and its changing weight (Liakopoulos-type
drainage). With z up, compression-positive p, tension-positive strain and
chi = Se, the equations solved are the DES strain-form rate equations:

    c(p) p_t + alpha Se eps_t = -q_z,          c = S0 + phi (1 - Sr) dSe/dp,
    q = -(k kr/mu) (p_z + rho_f g),
    eps_t = (sigma_t + alpha Pi'(p) p_t)/H_m,   Pi = Se p,   H_m = K + 4G/3,
    sigma_t(z) = -g rho_f phi (1 - Sr) int_z^H dSe/dp p_t dz'.

The last two make the storage nonlocal: drying above a point unloads it.
Vertex-centred finite volumes in z (half cells at the ends), the dense
coupled storage matrix inverted at each evaluation, and an implicit BDF
integrator with tight tolerances; three grids give an observed order and a
Richardson-extrapolated reference. This is an independent solution of the
same continuum, not of DES's discretization.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp


class Drainage:
    def __init__(self, *, height, bulk, shear, biot, porosity, fluid_bulk, grain_bulk, permeability, viscosity,
                 vg_alpha, vg_n, sr, rho=1000.0, g=9.81):
        self.H, self.alpha, self.phi, self.sr = height, biot, porosity, sr
        self.Hm = bulk + 4*shear/3
        self.S0 = porosity/fluid_bulk + (biot - porosity)/grain_bulk
        self.k, self.mu, self.rho, self.g = permeability, viscosity, rho, g
        self.a, self.n = vg_alpha, vg_n
        self.m = 1 - 1/vg_n

    # van Genuchten-Mualem, saturated for p >= 0
    def se(self, p):
        s = np.maximum(-p, 0.0)
        return (1 + (self.a*s)**self.n)**(-self.m)

    def dse(self, p):
        s = np.maximum(-p, 0.0)
        x = (self.a*s)**self.n
        d = self.m*self.n*self.a*(self.a*s)**(self.n - 1)*(1 + x)**(-self.m - 1)
        return np.where(p < 0, d, 0.0)

    def kr(self, p):
        se = self.se(p)
        return np.sqrt(se)*(1 - (1 - se**(1/self.m))**self.m)**2

    def _grid(self, n):
        z = np.linspace(0, self.H, n + 1)
        h = z[1] - z[0]
        w = np.full(n + 1, h); w[0] = w[-1] = h/2            # control-volume widths
        return z, h, w

    def solve(self, n, times, rtol=1e-8):
        """Pressure at the nodes, settlement (top displacement, up positive) and
        cumulative outflow per unit area at the given times."""
        z, h, w = self._grid(n)
        p0 = self.rho*self.g*(self.H - z)
        free = slice(1, None)                                # node 0 holds p = 0 for t > 0
        # Strictly-above weights for int_z^H f dz' on the vertex grid (trapezoid).
        upper = np.zeros((n + 1, n + 1))                     # row i: trapezoid weights of int_{z_i}^{H}
        for i in range(n):
            upper[i, i] = upper[i, n] = h/2
            upper[i, i + 1:n] = h
        coef = self.g*self.rho*self.phi*(1 - self.sr)

        def storage(p):
            se, dse = self.se(p), self.dse(p)
            pip = se + p*dse                                  # Pi'
            c = self.S0 + self.phi*(1 - self.sr)*dse
            M = np.diag(c + self.alpha**2*se*pip/self.Hm)
            M -= (self.alpha*se/self.Hm)[:, None]*coef*upper*dse[None, :]
            return M

        def flux(p):
            kr = self.kr(p)
            kf = 0.5*(kr[1:] + kr[:-1])                       # arithmetic face mobility, as DES
            return -(self.k*kf/self.mu)*((p[1:] - p[:-1])/h + self.rho*self.g)   # up-positive, faces

        def rate(t, y):
            p = np.concatenate(([0.0], y))
            q = flux(p)
            div = np.zeros(n + 1)
            div[:-1] -= q; div[1:] += q                       # net inflow into each control volume
            M = storage(p)
            # Node 0 is fixed: solve for the free nodes with the coupling to node 0 dropped (its rate is 0).
            return np.linalg.solve(M[free, free]*w[free, None], div[free])

        y0 = p0[free].copy()
        sol = solve_ivp(rate, (0.0, max(times)), y0, method='BDF', t_eval=sorted(times), rtol=rtol, atol=1e-6*p0.max(),
                        dense_output=False)
        if not sol.success:
            raise RuntimeError(sol.message)
        out = []
        se0, pi0 = self.se(p0), self.se(p0)*p0
        for k, t in enumerate(sol.t):
            p = np.concatenate(([0.0], sol.y[:, k]))
            dse = self.se(p) - se0
            dsig = -coef*(upper @ dse)
            deps = (dsig + self.alpha*(self.se(p)*p - pi0))/self.Hm
            settle = float(w @ deps)                          # int_0^H d eps dz: top displacement
            out.append(dict(time=float(t), z=z, p=p, settlement=settle,
                            saturation=self.sr + (1 - self.sr)*self.se(p)))
        return out


def liakopoulos_like():
    """Parameters of the U9 check: a 1 m sand column after Liakopoulos (1965),
    with a van Genuchten law of similar range."""
    E, nu = 1.3e6, 0.4
    return dict(height=1.0, bulk=E/(3*(1 - 2*nu)), shear=E/(2*(1 + nu)), biot=1.0, porosity=0.2975,
                fluid_bulk=2.2e9, grain_bulk=1e20, permeability=4.5e-13, viscosity=1e-3,
                vg_alpha=1/3000.0, vg_n=1.8, sr=0.1)


if __name__ == '__main__':
    d = Drainage(**liakopoulos_like())
    times = [300.0, 1200.0, 3600.0, 7200.0]
    runs = {n: d.solve(n, times) for n in (50, 100, 200)}
    # The trapezoid weights integrate exactly what they should.
    z, h, w = d._grid(10)
    for k, t in enumerate(times):
        s = [runs[n][k]['settlement'] for n in (50, 100, 200)]
        p = [runs[n][k]['p'][::{50: 1, 100: 2, 200: 4}[n]] for n in (50, 100, 200)]
        e1, e2 = np.max(np.abs(p[0] - p[1])), np.max(np.abs(p[1] - p[2]))
        print(f't = {t:6.0f} s: settlement {s} m, pressure change 50->100 {e1:.3g} Pa, 100->200 {e2:.3g} Pa,'
              f' order {math.log2(e1/e2):.2f}, min saturation {runs[200][k]["saturation"].min():.3f}')
