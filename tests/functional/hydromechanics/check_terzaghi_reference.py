#!/usr/bin/env python3
"""Check the Terzaghi reference against boundary conditions and balance laws."""

from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[3] /
                       'benchmarks' / 'poromechanics'))
from terzaghi_reference import TerzaghiReference


class TerzaghiReferenceChecks(unittest.TestCase):
    def setUp(self):
        self.ref = TerzaghiReference(height=10, load=1e5, bulk=1e7, shear=1e7,
                                    porosity=0.3, fluid_bulk=2.17e9,
                                    grain_bulk=37e9, biot=1.0,
                                    permeability=1.03e-12, viscosity=1.002e-3)

    def test_undrained_constitutive_balance(self):
        ref = self.ref
        strain = -ref.settlement(0)[0] / ref.height
        pressure = ref.pressure([1, 5, 10], 0)[0]
        np.testing.assert_allclose(ref.constrained_modulus * strain - ref.biot * pressure,
                                   -ref.load, atol=1e-9, rtol=1e-13)
        np.testing.assert_allclose(ref.storage * pressure + ref.biot * strain,
                                   0, atol=1e-16)

    def test_drained_limit(self):
        ref = self.ref
        late = 100 * ref.height**2 / ref.diffusivity
        np.testing.assert_allclose(ref.pressure([0, 1, 5, 10], late), 0, atol=1e-90)
        self.assertAlmostEqual(ref.settlement(late)[0], ref.load * ref.height /
                               ref.constrained_modulus, places=14)

    def test_drained_top_and_sealed_bottom(self):
        ref = self.ref
        times = np.array([1, 100, 1000, 10000])
        np.testing.assert_array_equal(ref.pressure([0], times), 0)
        step = 1e-3
        # Even reflection about the sealed bottom gives equal pressures.
        # Test the one-sided boundary slope, normalized by p_initial/H.
        values = ref.pressure([ref.height - step, ref.height], times)
        normalized_slope = (values[:, 1] - values[:, 0]) / step * ref.height / ref.initial_pressure
        self.assertLess(np.max(np.abs(normalized_slope)), 1e-3)

    def test_diffusion_equation(self):
        ref = self.ref
        depth = np.array([1., 5., 9.])
        dx, dt = 1e-3, 1e-2
        for time in (100., 1000.):
            rate = (ref.pressure(depth, time + dt) - ref.pressure(depth, time - dt)) / (2 * dt)
            curvature = (ref.pressure(depth + dx, time) - 2 * ref.pressure(depth, time) +
                         ref.pressure(depth - dx, time)) / dx**2
            scale = ref.initial_pressure * ref.diffusivity / ref.height**2
            self.assertLess(np.max(np.abs(rate - ref.diffusivity * curvature)) / scale, 1e-5)

    def test_settlement_from_independent_pressure_integral(self):
        ref = self.ref
        for time in (1., 100., 1000., 10000.):
            integrated, _ = quad(lambda d: ref.pressure(d, time)[0, 0], 0, ref.height)
            expected = (ref.load * ref.height - ref.biot * integrated) / ref.constrained_modulus
            self.assertAlmostEqual(ref.settlement(time)[0], expected, places=13)

    def test_load_pulse_has_the_correct_undrained_jump(self):
        ref = self.ref
        depth, end = [0., 1., 5., 10.], 1000.
        before = ref.pressure(depth, end)[0]
        p, settlement = ref.load_history(depth, [0., end], [(0., ref.load), (end, -ref.load)])
        np.testing.assert_allclose(p[0], ref.pressure(depth, 0)[0], atol=1e-10)
        np.testing.assert_allclose(p[1]-before, -ref.pressure(depth, 0)[0], atol=1e-10)
        self.assertAlmostEqual(settlement[1], ref.settlement(end)[0]-ref.settlement(0)[0], places=14)

    def test_load_history_settlement_from_pressure_integral(self):
        ref = self.ref
        events = [(0., ref.load), (1000., -.5*ref.load)]
        for time in (100., 1000., 2000., 10000.):
            pressure = lambda d: ref.load_history(d, time, events)[0][0, 0]
            integral = quad(pressure, 0, ref.height, epsabs=1e-7)[0]
            load = sum(increment for start, increment in events if time >= start)
            expected = (load*ref.height-ref.biot*integral)/ref.constrained_modulus
            actual = ref.load_history([0.], time, events)[1][0]
            self.assertAlmostEqual(actual, expected, places=13)

    def test_load_history_starts_only_after_its_first_event(self):
        ref = self.ref
        p, settlement = ref.load_history([0., 5., 10.], [0., 10.], [(100., ref.load)])
        np.testing.assert_array_equal(p, 0)
        np.testing.assert_array_equal(settlement, 0)

    def test_load_history_rejects_invalid_events(self):
        for events in ([(1., 1.), (0., 1.)], [(-1., 1.)], [(0., np.nan)], [[0., 1., 2.]]):
            with self.subTest(events=events), self.assertRaises(ValueError):
                self.ref.load_history([5.], [1.], events)


if __name__ == '__main__':
    unittest.main()
