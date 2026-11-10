#!/usr/bin/env python3
"""Check the analytical oracle; this does not run or qualify the DES solver."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[3] /
                       'benchmarks' / 'poromechanics'))
from theis_reference import drawdown, scheduled_drawdown


class TheisReferenceChecks(unittest.TestCase):
    def test_integral_reference(self):
        # Integrate the defining well-function integral independently of exp1.
        # With r=2, S=T=1, t=1/u, the prefactor Q/(4*pi*T) equals one.
        for u in (1e-6, 1e-3, 0.1, 1.0, 10.0):
            expected, _ = quad(lambda x: math.exp(-x) / x, u, np.inf,
                               epsabs=1e-14, epsrel=1e-12)
            actual = float(drawdown(2.0, 1 / u, 4 * math.pi, 1.0, 1.0))
            self.assertAlmostEqual(actual / expected, 1.0, places=10)

    def test_independent_dimensional_fixture(self):
        # E1(1)=0.21938393439552027; choose u=1 without unit-rate parameters.
        result = drawdown(20.0, 200.0, 0.01, 1e-3, 0.002)
        expected = 0.17458018796997565
        self.assertAlmostEqual(float(result), expected, places=13)

    def test_initial_state_and_injection(self):
        times = np.array([0.0, 100.0, 1000.0])
        pumping = drawdown(10.0, times, 0.01, 1e-3, 1e-3)
        injection = drawdown(10.0, times, -0.01, 1e-3, 1e-3)
        self.assertEqual(pumping[0], 0.0)
        np.testing.assert_array_equal(injection, -pumping)
        self.assertTrue((np.diff(pumping) > 0).all())
        np.testing.assert_array_equal(
            scheduled_drawdown(10.0, times, [], 1e-3, 1e-3), np.zeros(3))

    def test_spatial_broadcast(self):
        result = drawdown(np.array([1.0, 10.0, 100.0])[:, None],
                          np.array([1.0, 100.0, 10000.0])[None, :],
                          0.01, 1e-3, 1e-3)
        self.assertEqual(result.shape, (3, 3))
        self.assertTrue((np.diff(result, axis=0) < 0).all())

    def test_shut_in_against_finite_pulse_integral(self):
        # For t>tau, a finite pumping pulse integrates exp(-a/v)/v over
        # v in [t-tau,t], independently of the rate-change superposition.
        radius, trans, storage, rate, duration = 10.0, 1e-3, 1e-3, 0.01, 50.0
        a = radius**2 * storage / (4 * trans)
        schedule = [(0.0, rate), (duration, 0.0)]
        for time in (60.0, 100.0, 1000.0):
            integral, _ = quad(lambda v: math.exp(-a / v) / v,
                               time - duration, time, epsabs=1e-13)
            expected = rate / (4 * math.pi * trans) * integral
            actual = float(scheduled_drawdown(radius, time, schedule, trans, storage))
            self.assertAlmostEqual(actual / expected, 1.0, places=11)

    def test_delayed_start_and_rate_change(self):
        schedule = [(10.0, 0.01), (20.0, -0.01)]
        times = np.array([0.0, 10.0, 15.0, 20.0, 30.0])
        actual = scheduled_drawdown(10.0, times, schedule, 1e-3, 1e-3)
        np.testing.assert_array_equal(actual[:2], 0.0)
        # Integrate the two physical rate intervals directly, without calling
        # the constant-rate reference to construct the expected response.
        def impulse(age):
            return math.exp(-25.0 / age) / age if age > 0 else 0.0

        prefactor = 0.01 / (4 * math.pi * 1e-3)
        expected = np.array([0.0, 0.0,
                             prefactor * quad(impulse, 0, 5)[0],
                             prefactor * quad(impulse, 0, 10)[0],
                             prefactor * (quad(impulse, 10, 20)[0] -
                                          quad(impulse, 0, 10)[0])])
        np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=0)

    def test_invalid_inputs(self):
        for args in ((0, 1, 1, 1, 1), (1, -1, 1, 1, 1),
                     (1, 1, 1, 0, 1), (1, 1, 1, 1, -1),
                     (np.nan, 1, 1, 1, 1), (1, np.inf, 1, 1, 1),
                     (1, 1, np.inf, 1, 1)):
            with self.assertRaises(ValueError):
                drawdown(*args)
        for schedule in ([(1, 1), (1, 0)], [(2, 1), (1, 0)],
                         [(-1, 1)], [(0, np.nan)], [(np.inf, 1)]):
            with self.assertRaises(ValueError):
                scheduled_drawdown(1, 1, schedule, 1, 1)


if __name__ == '__main__':
    unittest.main()
