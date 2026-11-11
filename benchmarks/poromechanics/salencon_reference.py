"""Perfectly plastic circular opening in plane strain.

The infinite-domain stress solution is Salencon's cylindrical-cavity solution:
https://docs.itascacg.com/itasca940/3dec/docproject/source/examples/CylindricalHoleInAnIMCM.html
The finite-annulus extension follows radial equilibrium and matching at the
elastic/plastic interface. Displacement is relative to a uniformly compressed
initial state. Stress is tension positive; displacement is positive outward.
The radial and hoop stresses must remain the extreme principal stresses.
"""

import numpy as np
from scipy.optimize import brentq


def constants(friction, cohesion):
    sine = np.sin(np.deg2rad(friction))
    if not 0 < friction < 90 or cohesion <= 0:
        raise ValueError('Require 0 < friction < 90 degrees and positive cohesion')
    k = (1+sine)/(1-sine)
    q = 2*cohesion*np.cos(np.deg2rad(friction))/(1-sine)
    return k, q


def plastic_radius(opening, outer, pressure, cohesion, friction, inner=0.):
    """Finite outer boundary at r=outer; require a partially plastic annulus."""
    k, q = constants(friction, cohesion)
    if not (0 < opening < outer and 0 <= inner < pressure):
        raise ValueError('Invalid radii or compressive pressures')

    def mismatch(b):
        radial = (inner+q/(k-1))*(b/opening)**(k-1)-q/(k-1)
        interface = ((2*pressure-q)*outer**2+q*b**2)/((k+1)*outer**2-(k-1)*b**2)
        return radial-interface

    if not mismatch(opening) < 0 < mismatch(outer):
        raise ValueError('Reference requires a nonempty plastic zone and elastic exterior')
    return brentq(mismatch, opening, outer, xtol=1e-13)


def state(radius, opening, outer, pressure, cohesion, friction, bulk, shear, inner=0., dilation=0.):
    """Return rr, hoop, out-of-plane stress, radial displacement, er and et.

    On the intermediate-stress branch, e_r^p = -N_psi e_theta^p and e_y^p=0.
    Integrating u' + N_psi*u/r = e_r^e + N_psi*e_theta^e and matching the
    elastic displacement fixes u, including zero and associated dilation.
    """
    r = np.asarray(radius, float)
    if np.any((r < opening) | (r > outer)) or bulk <= 0 or shear <= 0:
        raise ValueError('Invalid radii or elastic moduli')
    k, q = constants(friction, cohesion)
    if not 0 <= dilation <= friction:
        raise ValueError('Require 0 <= dilation <= friction')
    sine = np.sin(np.deg2rad(dilation))
    npsi = (1+sine)/(1-sine)
    b = plastic_radius(opening, outer, pressure, cohesion, friction, inner)
    c = (inner+q/(k-1))/opening**(k-1)
    sb = c*b**(k-1)-q/(k-1)
    ae = (pressure*outer**2-sb*b**2)/(outer**2-b**2)
    be = b**2*outer**2*(pressure-sb)/(outer**2-b**2)
    plastic = r < b
    rr = -np.where(plastic, c*r**(k-1)-q/(k-1), ae-be/r**2)
    hoop = -np.where(plastic, k*c*r**(k-1)-q/(k-1), ae+be/r**2)
    lame = bulk-2*shear/3
    elastic_volume = (rr+hoop+2*pressure)/(2*(lame+shear))
    syy = -pressure+lame*elastic_volume
    ub = (pressure-ae)*b/(2*(lame+shear))-be/(2*shear*b)

    nu = lame/(2*(lame+shear))
    constant = (1+npsi)*(1-2*nu)*(pressure+q/(k-1))/(2*shear)
    power = c*((1-nu)-npsi*nu+k*(npsi*(1-nu)-nu))/(2*shear)

    def integral(x):
        return constant*x**(npsi+1)/(npsi+1)-power*x**(npsi+k)/(npsi+k)

    ue = (pressure-ae)*r/(2*(lame+shear))-be/(2*shear*r)
    up = (b**npsi*ub+integral(r)-integral(b))/r**npsi
    u = np.where(plastic, up, ue)
    et = u/r
    er = np.where(plastic, constant-power*r**(k-1)-npsi*et, elastic_volume-et)
    # Yielding in this reference uses hoop and radial principal stresses only.
    if np.any(syy < hoop-1e-8*pressure) or np.any(syy > rr+1e-8*pressure):
        raise ValueError('Out-of-plane stress leaves the intermediate principal branch')
    return rr, hoop, syy, u, er, et


def infinite_radius(opening, pressure, cohesion, friction, inner=0.):
    k, q = constants(friction, cohesion)
    return opening*(2*(pressure+q/(k-1))/((k+1)*(inner+q/(k-1))))**(1/(k-1))


def self_check():
    a, outer, pressure, cohesion, friction, bulk, shear = 1., 20., 5e6, 1e6, 30., 13e9/6, 1e9
    b = plastic_radius(a, outer, pressure, cohesion, friction)
    k, q = constants(friction, cohesion)
    r = np.linspace(a+1e-4, 5., 301)
    rr, hoop, syy, u, er, et = state(r, a, outer, pressure, cohesion, friction, bulk, shear)
    boundary = state(np.array([a, outer]), a, outer, pressure, cohesion, friction, bulk, shear)
    np.testing.assert_allclose(boundary[0], [0., -pressure], atol=1e-8)
    plastic = r < b
    np.testing.assert_allclose(-hoop[plastic]+k*rr[plastic], q, rtol=1e-14)
    assert np.all(-hoop[~plastic]+k*rr[~plastic] <= q+1e-7)
    delta = 1e-6
    left = state(r-delta, a, outer, pressure, cohesion, friction, bulk, shear)
    right = state(r+delta, a, outer, pressure, cohesion, friction, bulk, shear)
    np.testing.assert_allclose((right[0]-left[0])/(2*delta)+(rr-hoop)/r, 0., atol=.01)
    np.testing.assert_allclose((right[3]-left[3])/(2*delta), er, atol=1e-11)
    interface = state(np.array([b-1e-8, b+1e-8]), a, outer, pressure, cohesion, friction, bulk, shear)
    for values in interface:
        np.testing.assert_allclose(values[0], values[1], rtol=2e-7, atol=1e-10)
    # In-plane elastic reconstruction plus zero-volume plastic strain.
    lame = bulk-2*shear/3
    ev = er+et
    epr = er-(rr+pressure-lame*ev)/(2*shear)
    ept = et-(hoop+pressure-lame*ev)/(2*shear)
    np.testing.assert_allclose(epr+ept, 0., atol=1e-17)
    np.testing.assert_allclose(epr[~plastic], 0., atol=1e-17)
    # Limit of our finite-boundary construction versus the published infinite solution.
    bi = infinite_radius(a, pressure, cohesion, friction)
    np.testing.assert_allclose(plastic_radius(a, 1e6, pressure, cohesion, friction), bi, rtol=1e-11)
    nu = (3*bulk-2*shear)/(2*(3*bulk+shear))
    rp = r[plastic]
    for dilation in (0., 15., 30.):
        sine = np.sin(np.deg2rad(dilation)); npsi = (1+sine)/(1-sine)
        chi = ((2*nu-1)*(pressure+q/(k-1))+
               (1-nu)*(k*k-1)/(k+npsi)*q/(k-1)*(bi/a)**(k-1)*(bi/rp)**(npsi+1)+
               ((1-nu)*(k*npsi+1)/(k+npsi)-nu)*q/(k-1)*(rp/a)**(k-1))
        ui = state(rp, a, 1e6, pressure, cohesion, friction, bulk, shear, dilation=dilation)[3]
        np.testing.assert_allclose(-ui, rp*chi/(2*shear), rtol=1e-10)
        left = state(r-delta, a, outer, pressure, cohesion, friction, bulk, shear, dilation=dilation)
        right = state(r+delta, a, outer, pressure, cohesion, friction, bulk, shear, dilation=dilation)
        actual = state(r, a, outer, pressure, cohesion, friction, bulk, shear, dilation=dilation)
        np.testing.assert_allclose((right[3]-left[3])/(2*delta), actual[4], atol=1e-11)
    return dict(plastic_radius_m=b, infinite_plastic_radius_m=bi, checks=10)


if __name__ == '__main__':
    print(self_check())
