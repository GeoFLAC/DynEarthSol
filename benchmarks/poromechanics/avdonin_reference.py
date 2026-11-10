#!/usr/bin/env python3
"""Avdonin (1964): hot water injected from a line source into a porous layer.

Fluid at temperature T_in enters at a rate q (per unit thickness) from a
line at r = 0 into a layer at T0, with no heat loss through its top and
base. The flow is radial, v = q/(2 pi r), and the temperature obeys

    C_m T_t + rho_f c_f v T_r = lambda (1/r) (r T_r)_r,

with C_m the bulk heat capacity. Its similarity solution is

    theta = (T - T0)/(T_in - T0) = Gamma(nu, u)/Gamma(nu),
    u = r^2 C_m/(4 lambda t),    nu = rho_f c_f q/(4 pi lambda),

the regularized upper incomplete gamma function. The injected heat
rho_f c_f q (T_in - T0) t is all stored in the layer, since no heat leaves.
"""
import math

import numpy as np
from scipy.special import gammaincc


class Avdonin:
    def __init__(self, *, rate, conductivity, capacity, fluid_capacity, t0, t_in):
        self.q, self.lam, self.C, self.cf = rate, conductivity, capacity, fluid_capacity
        self.nu = fluid_capacity*rate/(4*math.pi*conductivity)
        self.t0, self.t_in = t0, t_in

    def theta(self, r, t):
        u = np.asarray(r, float)**2*self.C/(4*self.lam*np.asarray(t, float))
        return gammaincc(self.nu, u)

    def temperature(self, r, t):
        return self.t0 + (self.t_in - self.t0)*self.theta(r, t)

    def front(self, t):
        """Radius of the advective front, where the injected volume fills the pores' heat capacity."""
        return math.sqrt(self.cf*self.q*t/(math.pi*self.C))


if __name__ == '__main__':
    a = Avdonin(rate=2e-5, conductivity=2.0, capacity=2.5e6, fluid_capacity=4.18e6, t0=283.0, t_in=333.0)
    t = 5e5
    r = np.linspace(.05, 20, 4001)
    # The similarity solution satisfies the advection-conduction equation.
    dt = 1.0
    Tt = (a.temperature(r, t+dt) - a.temperature(r, t-dt))/(2*dt)
    T = a.temperature(r, t)
    Tr = np.gradient(T, r, edge_order=2)
    lap = np.gradient(r*Tr, r, edge_order=2)/r
    lhs = a.C*Tt + a.cf*a.q/(2*math.pi*r)*Tr
    res = (lhs - a.lam*lap)[20:-20]
    assert np.max(np.abs(res)) < 1e-3*np.max(np.abs(a.C*Tt)), np.max(np.abs(res))/np.max(np.abs(a.C*Tt))
    # All injected heat is stored: int C (T - T0) 2 pi r dr = rho_f c_f q (T_in - T0) t.
    rr = np.linspace(0, 60, 600001)
    stored = np.trapezoid(a.C*(a.temperature(rr, t) - a.t0)*2*math.pi*rr, rr)
    assert abs(stored/(a.cf*a.q*(a.t_in - a.t0)*t) - 1) < 1e-6, stored/(a.cf*a.q*(a.t_in - a.t0)*t)
    # Near the line the water arrives at T_in; far away the layer stays at T0.
    assert abs(a.theta(1e-6, t) - 1) < 1e-6 and a.theta(60.0, t) < 1e-12
    print(f'avdonin reference: 3 checks pass; nu = {a.nu:.3f}, front at {a.front(t):.2f} m after {t:g} s')
