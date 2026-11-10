#!/usr/bin/env python3
"""Check the unsaturated-flow references and the explicit prototype update:
constitutive limits and derivatives, the pointwise content inversion, the
Srivastava-Yeh series against an independent solution, the explicit column
against exact transient and steady infiltration, with its fluid balance and
stability bound, and the Bishop column settlement."""

import math
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parents[3] /
                       'benchmarks' / 'poromechanics'))
from unsaturated_reference import (Column, Retention, advance_content, bishop_column,
                                   bishop_pressure,
                                   content, content_slope, eigenvalues, srivastava_yeh,
                                   steady_pressure)

GAMMA = 1000.0 * 9.81
PHI = 0.4
LAWS = (Retention('van_genuchten', 1 / 5e3, 1.6, 0.1),
        Retention('brooks_corey', 1 / 3e3, 0.8, 0.05),
        Retention('gardner', 2 / GAMMA))
# Srivastava and Yeh (1991) example: alpha = 0.1/cm, L = 100 cm, Ks = 1 cm/h.
SY = dict(length=1.0, alpha=10.0, ks=1 / 3.6e5, dtheta=0.4)
SY_Q = dict(qa=0.1 * SY['ks'], qb=0.9 * SY['ks'])


def sy_matrix_exponential(z, T):
    """Independent solution of K_T = K_ZZ + K_Z, K(0) = 1, K_Z + K = qb at L,
    by central differences in Z and an exact matrix exponential in T."""
    Z = SY['alpha'] * z
    h, n = Z[1] - Z[0], len(Z)
    b = SY_Q['qb'] / SY['ks']
    a = np.zeros((n + 1, n + 1))  # last unknown is the constant 1
    for i in range(1, n - 1):
        a[i, i - 1:i + 2] = [1 / h**2 - 0.5 / h, -2 / h**2, 1 / h**2 + 0.5 / h]
    ghost = 1 / h**2 + 0.5 / h  # K_n = K_{n-2} + 2h(b - K_{n-1})
    a[n - 1, n - 2] = 2 / h**2
    a[n - 1, n - 1] = -2 / h**2 - 2 * h * ghost
    a[n - 1, n] = 2 * h * ghost * b
    k0 = srivastava_yeh(z, 0, **SY, **SY_Q)
    return (expm(a * T) @ np.append(k0, 1.0))[:n]


class Constitutive(unittest.TestCase):
    def test_saturated_limits(self):
        for ret in LAWS:
            with self.subTest(ret.model):
                p = np.array([ret.p_sat, 0.0, 1e5])
                np.testing.assert_array_equal(ret.se(p), 1.0)
                np.testing.assert_array_equal(ret.kr(p), 1.0)
                np.testing.assert_array_equal(ret.dse(p), 0.0)
                np.testing.assert_array_equal(bishop_pressure(ret, p), p)
                self.assertLess(float(ret.se(-1e9)), 1e-3)
                self.assertLess(float(ret.kr(-1e9)), 1e-9)

    def test_monotone_and_derivative(self):
        p = -np.geomspace(1.0, 1e7, 2000)
        for ret in LAWS:
            with self.subTest(ret.model):
                se, kr = ret.se(p), ret.kr(p)
                self.assertTrue(np.all(np.diff(se) <= 0) and np.all(np.diff(kr) <= 0))
                h = 1e-5 * np.abs(p)
                fd = (ret.se(p + h) - ret.se(p - h)) / (2 * h)
                exact = ret.dse(p)
                use = (exact > 1e-12 * ret.alpha) & (np.abs(p + 1 / ret.alpha) > 1e-3 / ret.alpha)
                np.testing.assert_allclose(fd[use], exact[use], rtol=1e-5, atol=1e-9 * ret.alpha)
                slope = content_slope(ret, PHI, 1e-9, p)
                np.testing.assert_allclose(slope, PHI * (1 - ret.sr) * exact + 1e-9)

    def test_closed_forms(self):
        vg, bc, ga = LAWS
        s = 7e3
        se = (1 + (vg.alpha * s)**vg.n) ** -vg.m
        self.assertAlmostEqual(float(vg.se(-s)), se, places=15)
        self.assertAlmostEqual(float(vg.kr(-s)), math.sqrt(se) *
                               (1 - (1 - se**(1 / vg.m))**vg.m)**2, places=14)
        self.assertAlmostEqual(float(bc.se(-s)), (bc.alpha * s)**-bc.n, places=15)
        self.assertAlmostEqual(float(bc.kr(-s)), (bc.alpha * s)**(-bc.n * (3 + 2 / bc.n)),
                               places=15)
        self.assertEqual(float(bc.se(-0.999 / bc.alpha)), 1.0)
        self.assertAlmostEqual(float(ga.kr(-s)), math.exp(-ga.alpha * s), places=15)
        for ret in LAWS:
            target = np.array([1e-6, 0.3, 0.9, 1 - 1e-9])
            np.testing.assert_allclose(ret.se(ret.pressure_of_se(target)), target, rtol=1e-9)

    def test_rejects_invalid(self):
        for args in (('nope', 1.0), ('gardner', 0.0), ('van_genuchten', 1.0, 1.0),
                     ('brooks_corey', 1.0, 0.0), ('gardner', 1.0, None, 1.0)):
            with self.assertRaises(ValueError):
                Retention(*args)


class Inversion(unittest.TestCase):
    def test_saturated_branch_is_the_existing_update(self):
        """Saturated nodes take p + dm/s exactly, as the current DES update."""
        rng = np.random.default_rng(3)
        for ret in LAWS:
            p = rng.uniform(1e5, 1e6, 500)
            dm = rng.uniform(-1e-4, 1e-4, 500)  # trial p + dm/s stays saturated
            np.testing.assert_array_equal(advance_content(ret, PHI, 1e-9, p, dm), p + dm / 1e-9)

    def test_residual_round_off(self):
        rng = np.random.default_rng(5)
        for ret in LAWS:
            for s in (1e-7, 1e-9, 0.0):
                with self.subTest(ret.model, storage=s):
                    p = np.concatenate([-np.geomspace(1, 1e6, 400), rng.uniform(-2e4, 2e5, 400)])
                    dm = rng.normal(0, 1e-3, p.size)
                    if s == 0:
                        target = ret.sat(p) + dm / PHI
                        keep = (ret.se(p) < 1) & (target > ret.sr + 1e-6) & (target < 1)
                        p, dm = p[keep], dm[keep]
                    x = advance_content(ret, PHI, s, p, dm)
                    g = PHI * (ret.sat(x) - ret.sat(p)) + s * (x - p) - dm
                    self.assertLess(np.max(np.abs(g)), 1e-14)

    def test_wetting_to_saturation_and_back(self):
        for ret in LAWS:
            p0 = np.array([-3e4])
            dm = (content(ret, PHI, 1e-9, 5e3) - content(ret, PHI, 1e-9, p0)).item()
            p1 = advance_content(ret, PHI, 1e-9, p0, dm)
            self.assertAlmostEqual(p1.item(), 5e3, delta=1e-6)
            np.testing.assert_allclose(advance_content(ret, PHI, 1e-9, p1, -dm), p0, rtol=1e-9)


class SrivastavaYeh(unittest.TestCase):
    def test_eigenvalues(self):
        L = SY['alpha'] * SY['length']
        lam = np.array(eigenvalues(L, 50))
        np.testing.assert_allclose(np.tan(lam * L) + 2 * lam, 0, atol=1e-9)
        self.assertTrue(np.all(np.diff(lam) > 0))

    def test_initial_and_steady_limits(self):
        z = np.linspace(0, SY['length'], 101)
        Z = SY['alpha'] * z
        a, b = SY_Q['qa'] / SY['ks'], SY_Q['qb'] / SY['ks']
        np.testing.assert_allclose(srivastava_yeh(z, 0, **SY, **SY_Q), a - (a - 1) * np.exp(-Z))
        t = 200 * SY['dtheta'] / (SY['alpha'] * SY['ks'])
        np.testing.assert_allclose(srivastava_yeh(z, t, **SY, **SY_Q), b - (b - 1) * np.exp(-Z),
                                   atol=1e-12)

    def test_against_matrix_exponential(self):
        z = np.linspace(0, SY['length'], 801)
        for T in (0.05, 0.5, 2.0):
            t = T * SY['dtheta'] / (SY['alpha'] * SY['ks'])
            exact = srivastava_yeh(z, t, **SY, **SY_Q)
            self.assertLess(np.max(np.abs(sy_matrix_exponential(z, T) - exact)), 3e-5)


class ExplicitColumn(unittest.TestCase):
    def test_gardner_transient_second_order(self):
        """Gardner column without storage, local bound, against the exact series."""
        ret = Retention('gardner', SY['alpha'] / GAMMA)
        perm = SY['ks'] * 1e-3 / GAMMA
        t = 0.5 * SY['dtheta'] / (SY['alpha'] * SY['ks'])
        errors = []
        for n in (21, 41, 81):
            z = np.linspace(0, SY['length'], n)
            col = Column(z, ret, phi=SY['dtheta'], storage=0.0, perm=perm)
            p0 = GAMMA * np.log(srivastava_yeh(z, 0, **SY, **SY_Q)) / SY['alpha']
            p, info = col.run(p0, t, SY_Q['qb'], 0.0, mode='local')
            self.assertLess(abs(info['balance']), 1e-15)
            self.assertLessEqual(info['courant'], 1 + 1e-12)
            k = np.exp(SY['alpha'] * p / GAMMA)
            errors.append(np.max(np.abs(k - srivastava_yeh(z, t, **SY, **SY_Q))))
        self.assertLess(errors[-1], 5e-4)
        for coarse, fine in zip(errors, errors[1:]):
            self.assertGreater(coarse / fine, 3.5)

    def test_steady_infiltration(self):
        """Any storage gives the same steady profile; the error falls with h."""
        perm = 1e-12
        q = 0.5 * perm * GAMMA / 1e-3
        for ret, ratio in ((Retention('van_genuchten', 1 / 2e3, 2.0, 0.1), 3.5),
                           (Retention('brooks_corey', 1 / 1e3, 1.5, 0.05), 1.5)):
            errors = []
            for n in (11, 21):
                z = np.linspace(0, 2.0, n)
                col = Column(z, ret, phi=PHI, storage=2e-5, perm=perm)
                p, info = col.run(-GAMMA * z, 4e5, q, 0.0)
                self.assertLessEqual(info['courant'], 1 + 1e-12)
                self.assertLess(abs(info['balance']), 1e-13)
                errors.append(np.max(np.abs(p - steady_pressure(ret, z, flux=q, perm=perm))))
            with self.subTest(ret.model):
                self.assertGreater(errors[0] / errors[1], ratio)

    def test_saturated_column_is_linear_diffusion(self):
        """All-saturated nodes reproduce the existing linear explicit update."""
        ret = LAWS[0]
        z = np.linspace(0, 10.0, 21)
        col = Column(z, ret, phi=PHI, storage=1e-9, perm=1e-13)
        p = GAMMA * (12.0 - z) + 1e3 * np.sin(z)
        dt = col.stable_dt(p)
        new, courant, _ = col.step(p, dt, 0.0, p[0])
        net = col.fluxes(p)[0]
        expected = p + dt * net / (col.vol * 1e-9)
        expected[0] = p[0]
        np.testing.assert_array_equal(new, expected)
        self.assertAlmostEqual(courant, 1.0, places=12)
        h = z[1] - z[0]
        self.assertAlmostEqual(dt, 0.5 * h * h * 1e-3 * 1e-9 / 1e-13, delta=1e-12 * dt)

    def test_saturated_bound_holds_for_dry_infiltration(self):
        """The existing saturated timestep needs no change when nodes desaturate."""
        ret = LAWS[0]
        perm = 1e-12
        z = np.linspace(0, 2.0, 21)
        col = Column(z, ret, phi=PHI, storage=1e-7, perm=perm)
        q = 0.5 * perm * GAMMA / 1e-3
        p, info = col.run(-GAMMA * z, 200.0, q, 0.0)
        self.assertLessEqual(info['courant'], 1 + 1e-12)
        self.assertLess(abs(info['balance']), 1e-15)
        self.assertGreater(info['stored'], 0.99 * q * 200.0)
        ref, _ = col.run(-GAMMA * z, 200.0, q, 0.0, mode='local', safety=0.05)
        self.assertLess(np.max(np.abs(p - ref)), 0.02 * np.max(np.abs(p + GAMMA * z)))

    def test_local_bound_retries_keep_monotone_steps(self):
        ret = Retention('gardner', SY['alpha'] / GAMMA)
        perm = SY['ks'] * 1e-3 / GAMMA
        z = np.linspace(0, SY['length'], 41)
        col = Column(z, ret, phi=SY['dtheta'], storage=0.0, perm=perm)
        p0 = GAMMA * np.log(srivastava_yeh(z, 0, **SY, **SY_Q)) / SY['alpha']
        _, info = col.run(p0, 2e4, SY_Q['qb'], 0.0, mode='local', safety=4.0)
        self.assertGreater(info['retries'], 0)
        self.assertLessEqual(info['courant'], 1 + 1e-12)


COLUMN = dict(height=10.0, phi=0.4, biot=0.9, bulk=2e7, shear=1e7)
MOD = COLUMN['bulk'] + 4 * COLUMN['shear'] / 3


def pi_antiderivative(ret, p):
    """G(p) with G' = Se*p and G(0) = 0, closed form for Gardner and Brooks-Corey."""
    if p >= ret.p_sat:
        return 0.5 * p * p
    if ret.model == 'gardner':
        a = ret.alpha
        return math.exp(a * p) * (p / a - 1 / a**2) + 1 / a**2
    a, n, s0 = ret.alpha, ret.n, 1 / ret.alpha
    return 0.5 * s0 * s0 + a**-n * ((-p)**(2 - n) - s0**(2 - n)) / (2 - n)


class BishopColumn(unittest.TestCase):
    def test_saturated_limit(self):
        """Both tables above the top: uniform Delta p, no density change."""
        z = np.linspace(0, 10, 11)
        for ret in LAWS:
            u = bishop_column(ret, z, table0=15.0, table1=12.0, **COLUMN)
            np.testing.assert_allclose(u, COLUMN['biot'] * (-3 * GAMMA) * z / MOD, rtol=1e-12)

    def test_closed_forms_without_density(self):
        z = np.array([1.0, 3.5, 6.0, 10.0])
        for ret in LAWS[1:]:
            for t0, t1 in ((10.0, 4.0), (7.0, 2.5), (2.0, 8.0)):
                u = bishop_column(ret, z, table0=t0, table1=t1, density=False, **COLUMN)
                # int_0^z pi dz = -(G(p(z)) - G(p(0)))/gamma, p = -gamma*(z - t)
                exact = [COLUMN['biot'] / MOD / GAMMA * sum(
                    sign * -(pi_antiderivative(ret, -GAMMA * (x - t)) - pi_antiderivative(ret, GAMMA * t))
                    for sign, t in ((1, t1), (-1, t0))) for x in z]
                np.testing.assert_allclose(u, exact, rtol=1e-9)

    def test_density_term_against_direct_integration(self):
        """Integrate stress from the top, then strain from the base, on a fine grid."""
        z = np.linspace(0, 10, 200001)
        for ret in LAWS:
            p0, p1 = -GAMMA * (z - 9.0), -GAMMA * (z - 3.0)
            dse = ret.se(p1) - ret.se(p0)
            above = cumulative_trapezoid(dse[::-1], z, initial=0)[::-1]  # int_z^H dSe
            dsig = -GAMMA * COLUMN['phi'] * (1 - ret.sr) * above
            eps = (dsig + COLUMN['biot'] * (bishop_pressure(ret, p1) - bishop_pressure(ret, p0))) / MOD
            u = cumulative_trapezoid(eps, z, initial=0)
            pick = np.arange(0, z.size, 20000)
            exact = bishop_column(ret, z[pick], table0=9.0, table1=3.0, **COLUMN)
            np.testing.assert_allclose(u[pick], exact, rtol=1e-7, atol=1e-12)

    def test_lowering_settles_and_drainage_unloads(self):
        for ret in LAWS:
            wet = bishop_column(ret, 10.0, table0=10.0, table1=4.0, density=False, **COLUMN)[0]
            dry = bishop_column(ret, 10.0, table0=10.0, table1=4.0, **COLUMN)[0]
            self.assertLess(wet, dry)
            self.assertLess(dry, 0)
            back = bishop_column(ret, 10.0, table0=4.0, table1=10.0, **COLUMN)[0]
            self.assertAlmostEqual(back, -dry, delta=1e-12 * abs(dry))


if __name__ == '__main__':
    unittest.main()
