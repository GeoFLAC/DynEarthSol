#!/usr/bin/env python3
"""Steady finite-amplitude convection in a square porous cell, independent of DES.

The Horton-Rogers-Lapwood cell of convection_reference.py beyond onset: unit
square, impermeable walls, theta = 1 at the base and 0 at the top, insulated
sides. In nondimensional Darcy-Boussinesq form, with u = psi_z, w = -psi_x,

    lap psi = -Ra theta_x,    u theta_x + w theta_z = lap theta,

and the Nusselt number Nu is the mean of -theta_z on the top (and of
-theta_z on the base; both are reported). Second-order central differences on
an (n+1)^2 node grid, solved by relaxed Picard iteration from a one-cell
seed, then Richardson extrapolation over three grids.
"""
import math

import numpy as np
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import splu, spsolve


def _operators(n):
    h = 1.0/n
    m = n + 1
    # Second derivative with mirror (insulated) ends, and with Dirichlet ends removed.
    d2n = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(m, m)).tolil()
    d2n[0, 1] = 2.0; d2n[-1, -2] = 2.0
    d1n = diags([-1.0, 1.0], [-1, 1], shape=(m, m)).tolil()
    d1n[0, 1] = 0.0; d1n[-1, -2] = 0.0             # theta_x = 0 at the sides
    return h, m, d2n.tocsr()/h**2, d1n.tocsr()/(2*h)


def solve(Ra, n, tol=1e-11, relax=0.7, iterations=5000):
    """Return (Nu_top, Nu_base, theta, psi) on the (n+1)^2 grid; index [k, i] = (z_k, x_i)."""
    h, m, d2x, d1x = _operators(n)
    x = np.linspace(0, 1, m); z = x.copy()
    X, Z = np.meshgrid(x, z)
    theta = 1 - Z + 0.05*np.cos(np.pi*X)*np.sin(np.pi*Z)
    # psi: Dirichlet zero on all walls -> interior unknowns only.
    k = m - 2
    lap1 = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(k, k))/h**2
    lap_psi = splu((kron(identity(k), lap1) + kron(lap1, identity(k))).tocsc())
    # theta: unknown in rows 1..m-2 (fixed top/base), all columns (insulated sides).
    rows = m - 2
    d2z = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(rows, rows))/h**2
    lap_t = kron(identity(rows), d2x) + kron(d2z, identity(m))
    dz_t = kron(diags([-1.0, 1.0], [-1, 1], shape=(rows, rows))/(2*h), identity(m))
    dx_t = kron(identity(rows), d1x)
    # Boundary contributions of theta = 1 at z = 0 (row k = 0) to the first interior row.
    bc_lap = np.zeros((rows, m)); bc_lap[0, :] = 1.0/h**2
    bc_dz = np.zeros((rows, m)); bc_dz[0, :] = -1.0/(2*h)
    psi = np.zeros((m, m))
    for it in range(iterations):
        tx = (d1x @ theta.T).T                                    # theta_x on every node
        psi_in = lap_psi.solve((-Ra*tx[1:-1, 1:-1]).ravel())
        psi[1:-1, 1:-1] = psi_in.reshape(k, k)
        u = np.gradient(psi, h, axis=0, edge_order=2)            # psi_z
        w = -np.gradient(psi, h, axis=1, edge_order=2)           # -psi_x
        U, W = u[1:-1, :].ravel(), w[1:-1, :].ravel()
        A = (lap_t - diags(U) @ dx_t - diags(W) @ dz_t).tocsc()
        rhs = -(bc_lap - W.reshape(rows, m)*bc_dz).ravel()
        new = theta.copy()
        new[1:-1, :] = spsolve(A, rhs).reshape(rows, m)
        change = np.max(np.abs(new - theta))
        theta = theta + relax*(new - theta)
        if change < tol:
            break
    else:
        raise RuntimeError(f'Picard iteration did not converge for Ra = {Ra}, n = {n}')
    top = -(3*theta[-1] - 4*theta[-2] + theta[-3])/(2*h)
    base = -(-3*theta[0] + 4*theta[1] - theta[2])/(2*h)
    trap = lambda f: h*(f.sum() - .5*(f[0] + f[-1]))
    return trap(top), trap(base), theta, psi


def nusselt(Ra, grids=(32, 64, 128)):
    """Richardson-extrapolated Nu and the finest-grid estimate of its error."""
    nu = [solve(Ra, n)[0] for n in grids]
    order = math.log(abs(nu[0] - nu[1])/abs(nu[1] - nu[2]), 2)
    extrapolated = nu[2] + (nu[2] - nu[1])/(2**2 - 1)
    return dict(Ra=Ra, grids=list(grids), nu=nu, observed_order=order, nu_extrapolated=extrapolated,
                error_estimate=abs(extrapolated - nu[2]))


if __name__ == '__main__':
    # Below onset the conductive state is the only steady state: Nu = 1.
    assert abs(solve(30.0, 32)[0] - 1) < 1e-8
    # Top and base heat flows balance at steady state.
    t, b, _, _ = solve(100.0, 64)
    assert abs(t - b) < 2e-3*t, (t, b)
    for Ra in (100.0, 200.0):
        r = nusselt(Ra)
        print(f"Ra = {Ra:g}: Nu = {r['nu']} -> {r['nu_extrapolated']:.5f} (order {r['observed_order']:.2f})")
