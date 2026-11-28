#!/usr/bin/env python3
"""Check the thickness-averaged aquifer reference against its limits,
the printed Bear-Corapcioglu (1981b) results, balance laws, boundary
conditions and an independent finite-volume solution."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import exp1

sys.path.insert(0, str(Path(__file__).resolve().parents[3] /
                       'benchmarks' / 'poromechanics'))
from averaged_aquifer_reference import AveragedAquifer, bear_corapcioglu_1981b

YEAR = 365.25 * 86400
PROPS = dict(thickness=50.0, bulk=2e8, shear=1e8, biot=0.8, porosity=0.25,
             fluid_bulk=2.2e9, grain_bulk=3.7e10, permeability=1e-12, viscosity=1e-3)
RATE, RW, R = 0.01, 0.1, 2000.0


def finite_volume(aq, radius, times, outer_hydraulic, inner_mechanics, outer_mechanics,
                  cells=400):
    """Independent radial finite-volume solution of the same averaged equations.

    The uniform areal strain h is eliminated through the outer mechanical
    condition, which makes it a linear functional of the cell pressures.
    """
    m, storage, coupling = aq.mobility, aq.specific_storage, aq.coupling
    lam_ps = 2 * aq.shear * aq.lame / (aq.lame + 2 * aq.shear)
    shear, a2 = aq.shear, aq.biot / (2 * (aq.lame + aq.shear))
    faces = np.geomspace(RW, R, cells + 1)
    centres = np.sqrt(faces[:-1] * faces[1:])
    weight = (faces[1:]**2 - faces[:-1]**2) / 2
    matrix = np.zeros((cells, cells))
    for i in range(cells - 1):
        g = m * faces[i + 1] / (centres[i + 1] - centres[i])
        matrix[i, i] -= g
        matrix[i + 1, i + 1] -= g
        matrix[i, i + 1] += g
        matrix[i + 1, i] += g
    if outer_hydraulic == 'fixed_pressure':
        matrix[-1, -1] -= m * faces[-1] / (faces[-1] - centres[-1])
    source = np.zeros(cells)
    source[0] = -RATE / (2 * math.pi * aq.thickness)
    inner = 0.0 if inner_mechanics == 'fixed' else (lam_ps + 2 * shear) * RW**2 / (2 * shear)
    area = (R**2 - RW**2) / 2
    if coupling == 0:
        functional = np.zeros(cells)
    elif outer_mechanics == 'fixed':
        functional = -a2 * weight / (inner + area)
    else:
        functional = 2 * shear * a2 * weight / ((lam_ps + 2 * shear) * R**2 - 2 * shear * (inner + area))
    mass = np.diag(weight * storage) + coupling * np.outer(weight, functional)
    jac = np.linalg.solve(mass, matrix)
    forcing = np.linalg.solve(mass, source)
    sol = solve_ivp(lambda t, p: jac @ p + forcing, (0, times[-1]), np.zeros(cells),
                    method='BDF', jac=jac, t_eval=times, rtol=1e-9, atol=1e-6)
    pressure = np.array([np.interp(np.log(radius), np.log(centres), row) for row in sol.y.T])
    return pressure, sol.y.T @ functional


def radial_stress(aq, fields, radius):
    eh, ez, ur, p = (fields[k] for k in ('areal_strain', 'vertical_strain',
                                         'radial_displacement', 'pressure'))
    volumetric = eh + ez
    srr = aq.lame * volumetric + 2 * aq.shear * (eh - ur / radius) - aq.biot * p
    stt = aq.lame * volumetric + 2 * aq.shear * ur / radius - aq.biot * p
    szz = aq.lame * volumetric + 2 * aq.shear * ez - aq.biot * p
    return srr, stt, szz


class AveragedAquiferChecks(unittest.TestCase):
    def setUp(self):
        self.aq = AveragedAquifer(**PROPS)
        self.radius = np.array([1.0, 10.0, 100.0, 500.0, 1000.0])
        self.times = np.array([1e3, 1e4, 1e5, 1e6, 1e7])

    def test_paper_closed_forms_are_the_terzaghi_limit(self):
        aq = AveragedAquifer(thickness=142.0, bulk=1e8, shear=5e7, biot=1.0, porosity=0.3,
                             fluid_bulk=math.inf, grain_bulk=math.inf,
                             permeability=1e-12, viscosity=1e-3)
        radius = np.array([10., 100., 1000., 5000.])
        for well_radius in (None, 0.2):
            ours = aq.line_well(radius, 1e7, RATE, well_radius=well_radius)
            paper = bear_corapcioglu_1981b(radius, 1e7, RATE, aq.transmissivity,
                                           aq.diffusivity, aq.thickness, well_radius)
            for key, value in paper.items():
                np.testing.assert_allclose(ours[key], value, rtol=1e-13, err_msg=key)

    def test_part2_example_radial_displacement_peak(self):
        # T=95 cm2/s, B=142 m, Q=50 L/s, Cv=6e3 cm2/s, 3 years: Fig. 3 peak
        # 17.1 cm near 8.56 km; the peak satisfies W(u) = (1-exp(-u))/u.
        cv, time = 0.6, 3 * YEAR
        radius = np.linspace(1000, 20000, 19001)
        ur = bear_corapcioglu_1981b(radius, time, 0.05, 95e-4, cv, 142.0)['radial_displacement']
        peak = np.argmax(-ur)
        u_peak = brentq(lambda u: exp1(u) - (1 - math.exp(-u)) / u, 0.1, 1)
        self.assertAlmostEqual(u_peak, 0.3214, places=4)
        self.assertAlmostEqual(radius[peak], 2 * math.sqrt(u_peak * cv * time), delta=1.0)
        self.assertAlmostEqual(-ur[peak], 0.1706, delta=5e-4)

    def test_infinite_aquifer_splits_strain_equally(self):
        fields = self.aq.line_well(self.radius, 1e5, RATE, well_radius=RW)
        np.testing.assert_allclose(fields['vertical_strain'], fields['areal_strain'], rtol=1e-14)
        _, _, szz = radial_stress(self.aq, fields, self.radius)
        np.testing.assert_allclose(szz, 0, atol=1e-12 * np.max(np.abs(fields['pressure'])))

    def test_inversion_converges(self):
        args = (self.radius, self.times, RATE)
        kwargs = dict(well_radius=RW, outer_radius=R)
        coarse = self.aq.finite_well(*args, talbot_terms=24, **kwargs)
        fine = self.aq.finite_well(*args, talbot_terms=40, **kwargs)
        for key in ('pressure', 'radial_displacement', 'subsidence', 'h'):
            scale = np.max(np.abs(fine[key]))
            self.assertLess(np.max(np.abs(coarse[key] - fine[key])) / scale, 1e-6, key)

    def test_matches_independent_finite_volume_solver(self):
        for outer_hydraulic, inner, outer in (('fixed_pressure', 'fixed', 'fixed'),
                                              ('no_flow', 'free', 'free'),
                                              ('fixed_pressure', 'fixed', 'free')):
            exact = self.aq.finite_well(self.radius, self.times, RATE, well_radius=RW,
                                        outer_radius=R, outer_hydraulic=outer_hydraulic,
                                        inner_mechanics=inner, outer_mechanics=outer)
            pressure, h = finite_volume(self.aq, self.radius, self.times, outer_hydraulic,
                                        inner, outer)
            label = f'{outer_hydraulic} {inner}/{outer}'
            self.assertLess(np.max(np.abs(exact['pressure'] - pressure)) /
                            np.max(np.abs(exact['pressure'])), 1e-4, label)
            self.assertLess(np.max(np.abs(exact['h'] - h)) / np.max(np.abs(exact['h'])), 1e-3, label)

    def test_sealed_aquifer_fluid_balance(self):
        radius = np.geomspace(RW, R, 4001)
        for inner, outer in (('fixed', 'fixed'), ('free', 'free')):
            fields = self.aq.finite_well(radius, self.times, RATE, well_radius=RW, outer_radius=R,
                                         outer_hydraulic='no_flow', inner_mechanics=inner,
                                         outer_mechanics=outer)
            stored = np.trapezoid(2 * math.pi * radius * self.aq.thickness *
                                  fields['fluid_content'], radius, axis=1)
            np.testing.assert_allclose(stored, -RATE * self.times, rtol=1e-5)

    def test_mechanical_boundaries_and_equilibrium(self):
        # Central differences: truncation ~step^2, inversion noise ~1e-9/step.
        step = 1e-3
        radius = np.array([RW, 3.0, 300.0, R])
        probe = np.concatenate([radius, radius[1:3] * (1 - step), radius[1:3] * (1 + step)])
        for inner, outer in (('fixed', 'fixed'), ('free', 'free'), ('fixed', 'free'),
                             ('free', 'fixed')):
            fields = self.aq.finite_well(probe, self.times, RATE, well_radius=RW, outer_radius=R,
                                         inner_mechanics=inner, outer_mechanics=outer)
            srr, stt, szz = radial_stress(self.aq, fields, probe)
            scale = np.max(np.abs(self.aq.biot * fields['pressure']))
            ur = fields['radial_displacement']
            for index, option in ((0, inner), (3, outer)):
                value = ur[:, index] / probe[index] if option == 'fixed' else srr[:, index] / scale
                self.assertLess(np.max(np.abs(value)), 1e-10, f'{inner}/{outer} boundary {index}')
            self.assertLess(np.max(np.abs(szz)) / scale, 1e-10)
            dsrr = (srr[:, 6:8] - srr[:, 4:6]) / (2 * step * radius[1:3])
            residual = dsrr + (srr[:, 1:3] - stt[:, 1:3]) / radius[1:3]
            self.assertLess(np.max(np.abs(residual * radius[1:3])) / scale, 1e-5,
                            f'{inner}/{outer} equilibrium')
            # Areal strain is the divergence of the radial displacement.
            rur = probe * ur
            divergence = (rur[:, 6:8] - rur[:, 4:6]) / (2 * step * radius[1:3]) / radius[1:3]
            eh = fields['areal_strain'][:, 1:3]
            self.assertLess(np.max(np.abs(divergence - eh)) / np.max(np.abs(eh)), 1e-5,
                            f'{inner}/{outer} kinematics')

    def test_large_domain_recovers_line_well(self):
        radius, times = self.radius, self.times[:4]
        for mechanics in ('plane_stress', 'uniaxial'):
            aq = AveragedAquifer(**PROPS, mechanics=mechanics)
            finite = aq.finite_well(radius, times, RATE, well_radius=RW, outer_radius=2e5)
            line = aq.line_well(radius[None, :], times[:, None], RATE, well_radius=RW)
            for key in ('pressure', 'radial_displacement', 'subsidence'):
                scale = np.max(np.abs(line[key])) or 1.0
                self.assertLess(np.max(np.abs(finite[key] - line[key])) / scale, 1e-4,
                                f'{mechanics} {key}')

    def test_steady_state_is_thiem(self):
        late = np.array([100 * R**2 / self.aq.diffusivity])
        fields = self.aq.finite_well(self.radius, late, RATE, well_radius=RW, outer_radius=R)
        thiem = RATE / (2 * math.pi * self.aq.mobility * self.aq.thickness) * np.log(self.radius / R)
        np.testing.assert_allclose(fields['pressure'][0], thiem, rtol=1e-8)

    def test_rejects_invalid_inputs(self):
        with self.assertRaises(ValueError):
            AveragedAquifer(**dict(PROPS, biot=0.2))  # porosity > Biot coefficient
        with self.assertRaises(ValueError):
            AveragedAquifer(**PROPS, mechanics='plane_strain')
        with self.assertRaises(ValueError):
            self.aq.finite_well(self.radius, self.times, RATE, well_radius=2.0, outer_radius=R)
        with self.assertRaises(ValueError):
            self.aq.finite_well(self.radius, [0.0], RATE, well_radius=RW, outer_radius=R)
        with self.assertRaises(ValueError):
            self.aq.finite_well([np.nan], self.times, RATE, well_radius=RW, outer_radius=R)
        with self.assertRaises(ValueError):
            self.aq.finite_well(self.radius, self.times, math.nan, well_radius=RW, outer_radius=R)
        with self.assertRaises(ValueError):
            AveragedAquifer(**PROPS, mechanics='uniaxial').finite_well(
                self.radius, self.times, RATE, well_radius=RW, outer_radius=R,
                inner_mechanics='free')


if __name__ == '__main__':
    unittest.main()
