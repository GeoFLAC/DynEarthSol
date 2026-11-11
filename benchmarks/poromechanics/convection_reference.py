#!/usr/bin/env python3
"""Horton-Rogers-Lapwood linear stability of a porous layer heated from below.

Boussinesq Darcy flow, q = -(k/mu)(grad p + rho0*(1 - alpha*(T - T0))*g*e_z),
incompressible, with conduction and advection
C dT/dt + rho0*cf*q.grad T = lambda*lap T, C the bulk heat capacity.
Top and bottom: impermeable, fixed T; sides: impermeable, insulated.
A box of width L and height H admits theta ~ cos(m*pi*x/L)*sin(n*pi*z/H).
With a = m*pi*H/L and K2 = a^2 + (n*pi)^2, the growth rate is

    sigma = lambda/(C*H^2) * (Ra*a^2/K2 - K2),
    Ra = rho0*g*alpha*dT*k*H*rho0*cf/(mu*lambda).

Onset (sigma = 0) is Ra = K2^2/a^2, minimal 4*pi^2 at a = pi (L = H, m = n = 1).
Reference: Nield and Bejan, Convection in Porous Media, section 6.2.
"""
import math


def rayleigh(rho0, g, alpha, dT, k, H, cf, mu, conductivity):
    return rho0*g*alpha*dT*k*H*rho0*cf/(mu*conductivity)


def mode_numbers(H, L, m, n, transverse=0, width_y=None):
    a = m*math.pi*H/L
    a2 = a*a
    if transverse:
        a2 += (transverse*math.pi*H/(L if width_y is None else width_y))**2
    return a2, a2+(n*math.pi)**2


def growth_rate(Ra, conductivity, capacity, H, L, m=1, n=1, transverse=0, width_y=None):
    a2, K2 = mode_numbers(H,L,m,n,transverse,width_y)
    return conductivity/(capacity*H*H)*(Ra*a2/K2 - K2)


def critical_rayleigh(H, L, m=1, n=1, transverse=0, width_y=None):
    a2, K2 = mode_numbers(H,L,m,n,transverse,width_y)
    return K2*K2/a2


def storage_growth_rate(Ra, hydraulic_ratio, conductivity, capacity, H, L,
                        m=1, n=1, transverse=0, width_y=None):
    """Dominant real mode with finite pressure storage, D=Dp/(lambda/C).

    The linear pressure/heat equations give lambda_mode^2 + A*lambda_mode
    + B=0, A=(D+1)*K2-Ra, B=D*(K2^2-Ra*a2). The incompressible growth
    rate follows as D tends to infinity. This exact finite-D reference lets
    a scaled-storage test retain explicit pressure subcycling economically.
    """
    if hydraulic_ratio<=0:raise ValueError('Pressure diffusivity must be positive')
    a2,K2=mode_numbers(H,L,m,n,transverse,width_y)
    A=(hydraulic_ratio+1)*K2-Ra
    B=hydraulic_ratio*(K2*K2-Ra*a2)
    discriminant=A*A-4*B
    if discriminant<0:raise ValueError('This mode is oscillatory, not a real exponential')
    root=math.sqrt(discriminant)
    value=-2*B/(A+root) if A>0 else .5*(-A+root)
    return conductivity/(capacity*H*H)*value


if __name__ == '__main__':
    assert abs(critical_rayleigh(1, 1) - 4*math.pi**2) < 1e-12
    assert abs(growth_rate(4*math.pi**2, 1, 1, 1, 1)) < 1e-12
    # The square cell is the most unstable single-cell box.
    assert all(critical_rayleigh(1, L) >= critical_rayleigh(1, 1) for L in (0.5, 0.9, 1.1, 2))
    # A doubled Ra above onset: sigma*C*H^2/lambda = Ra/2 - 2*pi^2 for a = pi.
    assert abs(growth_rate(8*math.pi**2, 1, 1, 1, 1) - 2*math.pi**2) < 1e-12
    print('convection reference: 4 checks pass')
