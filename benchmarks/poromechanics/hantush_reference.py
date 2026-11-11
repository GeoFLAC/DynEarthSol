#!/usr/bin/env python3
"""Hantush-Jacob (1955): pumping from a leaky confined aquifer.

A fully penetrating well pumps at rate Q from an aquifer of transmissivity T
and storativity S. Above it lies an aquitard of thickness b' and permeability
k' whose top is held at the initial head. The aquitard stores no water and
its flow is vertical, so it leaks a flux proportional to the local drawdown,
and

    s(r,t) = Q/(4 pi T) W(u, r/lambda),   u = r^2 S/(4 T t),
    W(u, beta) = int_u^inf exp(-y - beta^2/(4y)) dy/y,   lambda^2 = T b'/K'.

In pressure units, as in theis_reference.py, T = k B/mu, S = S_s B and the
aquitard leakance is K'/b' = k'/(mu b'), so lambda^2 = k B b'/k'.
W tends to the Theis E1(u) as lambda grows and to the steady 2 K0(r/lambda)
at late time.
Hantush, M. S. and Jacob, C. E. (1955), Trans. AGU 36(1), 95-100.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import exp1, k0


def well_function(u, beta):
    """Leaky well function W(u, beta), elementwise; u > 0, beta >= 0."""
    u, beta = np.broadcast_arrays(np.asarray(u, float), np.asarray(beta, float))
    out = np.empty(u.shape)
    for i, (a, b) in enumerate(zip(u.ravel(), beta.ravel())):
        if b == 0:
            out.flat[i] = exp1(a)
            continue
        # y = exp(s): the integrand exp(-y - b^2/(4y)) peaks at y = b/2.
        f = lambda s: math.exp(-math.exp(s) - b*b/(4*math.exp(s)))
        lo, peak = math.log(a), math.log(max(b/2, a))
        parts = [(lo, peak), (peak, peak + 60)] if peak > lo else [(lo, lo + 60)]
        out.flat[i] = sum(quad(f, x0, x1, epsabs=0, epsrel=1e-12, limit=200)[0] for x0, x1 in parts)
    return out


def drawdown(radius, time, rate, transmissivity, storativity, leakage_factor):
    """Drawdown for a constant rate from t = 0 (zero for t <= 0)."""
    radius, time = np.broadcast_arrays(np.asarray(radius, float), np.asarray(time, float))
    out = np.zeros(radius.shape)
    on = time > 0
    if rate == 0 or not on.any():
        return out
    u = radius[on]**2*storativity/(4*transmissivity*time[on])
    out[on] = rate/(4*math.pi*transmissivity)*well_function(u, radius[on]/leakage_factor)
    return out


def scheduled_drawdown(radius, time, schedule, transmissivity, storativity, leakage_factor):
    """Superpose absolute rates given by ordered (start_seconds, rate) pairs."""
    out = np.zeros(np.broadcast(np.asarray(radius), np.asarray(time)).shape)
    previous = 0.0
    for start, rate in schedule:
        lag = np.asarray(time, float) - start
        out += drawdown(radius, lag, rate - previous, transmissivity, storativity, leakage_factor)
        previous = rate
    return out


if __name__ == '__main__':
    beta = np.array([0.01, 0.1, 0.5, 1.0, 2.0])
    # Large lambda recovers Theis.
    assert np.allclose(well_function(1e-3, 1e-9), exp1(1e-3), rtol=1e-9)
    # Late time is the steady 2 K0(beta).
    assert np.allclose(well_function(1e-12, beta), 2*k0(beta), rtol=1e-6)
    # Hantush's identity W(u, b) + W(b^2/(4u), b) = 2 K0(b), an independent check
    # of the quadrature on both sides of the peak.
    for u in (1e-3, 0.05, 0.4, 3.0):
        assert np.allclose(well_function(u, beta) + well_function(beta**2/(4*u), beta), 2*k0(beta), rtol=1e-10)
    # The drawdown solves S s_t = T (s_rr + s_r/r) - T s/lambda^2.
    T, S, lam, Q = 2.0, 1.0, 15.0, 1e-2          # u from 0.1 to 3.8: transient, not steady
    r = np.linspace(5, 30, 501); t, dt = 30.0, 1e-3
    s = drawdown(r, t, Q, T, S, lam)
    st = (drawdown(r, t+dt, Q, T, S, lam) - drawdown(r, t-dt, Q, T, S, lam))/(2*dt)
    sr = np.gradient(s, r, edge_order=2); srr = np.gradient(sr, r, edge_order=2)
    lhs, rhs = S*st[5:-5], (T*(srr + sr/r) - T*s/lam**2)[5:-5]
    assert np.allclose(lhs, rhs, rtol=1e-3, atol=1e-6*np.abs(lhs).max()), np.max(np.abs(lhs-rhs)/np.abs(lhs).max())
    print('hantush reference: 4 checks pass')
