"""Elastic circular opening: Kirsch and the finite pressure annulus.

Stress is tension positive. The full plane-strain Kirsch stress solution is
given, for example, in the Itasca verification problem:
https://docs.itascacg.com/itasca950/3dec/docproject/source/examples/CylindricalHoleInAnIEM.html
The finite annulus follows directly from radial equilibrium and Hooke's law;
it separates truncation error from the infinite-domain Kirsch limit.
"""

import numpy as np


def kirsch(radius, angle, opening, stress_x, stress_z):
    """Return radial, hoop and radial/tangential shear stresses."""
    r, theta = np.broadcast_arrays(np.asarray(radius, float), np.asarray(angle, float))
    if opening <= 0 or np.any(r < opening):
        raise ValueError('Require positive opening radius and observations outside it')
    ratio = (opening/r)**2
    mean, difference = .5*(stress_x+stress_z), .5*(stress_x-stress_z)
    rr = mean*(1-ratio) + difference*(1-4*ratio+3*ratio**2)*np.cos(2*theta)
    hoop = mean*(1+ratio) - difference*(1+3*ratio**2)*np.cos(2*theta)
    shear = -difference*(1+2*ratio-3*ratio**2)*np.sin(2*theta)
    return rr, hoop, shear


def annulus(radius, opening, outer_radius, pressure, inner_pressure=0.):
    """Finite annulus with compressive pressures; return rr and hoop stress."""
    r = np.asarray(radius, float)
    if not (0 < opening < outer_radius) or np.any((r < opening) | (r > outer_radius)):
        raise ValueError('Require opening <= radius <= outer_radius')
    a2, b2 = opening**2, outer_radius**2
    constant = (inner_pressure*a2-pressure*b2)/(b2-a2)
    inverse_square = a2*b2*(inner_pressure-pressure)/(b2-a2)
    return constant-inverse_square/r**2, constant+inverse_square/r**2


def displacement(radius, opening, outer_radius, pressure, bulk, shear, inner_pressure=0.):
    """Plane-strain annular displacement, positive outward."""
    r = np.asarray(radius, float)
    a2, b2 = opening**2, outer_radius**2
    constant = (inner_pressure*a2-pressure*b2)/(b2-a2)
    inverse_square = a2*b2*(inner_pressure-pressure)/(b2-a2)
    lame = bulk-2*shear/3
    return constant*r/(2*(lame+shear)) + inverse_square/(2*shear*r)


if __name__ == '__main__':
    angle = np.linspace(-np.pi, np.pi, 101)
    rr, hoop, shear = kirsch(np.ones_like(angle), angle, 1., -2., -1.)
    np.testing.assert_allclose(rr, 0., atol=1e-15)
    np.testing.assert_allclose(shear, 0., atol=1e-15)
    np.testing.assert_allclose(kirsch(1., 0., 1., -1., -1.)[1], -2.)
    far = kirsch(1e6, angle, 1., -2., -1.)
    np.testing.assert_allclose(far[0], -2*np.cos(angle)**2-np.sin(angle)**2, atol=1e-11)
    np.testing.assert_allclose(annulus(np.array([1., 20.]), 1., 20., 3., 1.)[0], [-1., -3.])
    r = np.linspace(1.1, 5., 100)
    rr, hoop = annulus(r, 1., 20., 1.)
    delta = 1e-5
    derivative = (annulus(r+delta, 1., 20., 1.)[0]-annulus(r-delta, 1., 20., 1.)[0])/(2*delta)
    np.testing.assert_allclose(derivative+(rr-hoop)/r, 0, atol=1e-9)
    np.testing.assert_allclose(annulus(r, 1., 1e6, 1.)[0], kirsch(r, 0, 1., -1., -1.)[0], atol=1e-11)
    bulk, modulus = 8., 2.
    u = displacement(r, 1., 20., 1., bulk, modulus)
    er = (displacement(r+delta, 1., 20., 1., bulk, modulus)-
          displacement(r-delta, 1., 20., 1., bulk, modulus))/(2*delta)
    et = u/r
    lame = bulk-2*modulus/3
    np.testing.assert_allclose(2*modulus*er+lame*(er+et), rr, atol=1e-9)
    np.testing.assert_allclose(2*modulus*et+lame*(er+et), hoop, atol=1e-9)
    print('Kirsch/annulus reference: boundary, far-field, concentration, equilibrium and limit checks pass')
