#!/usr/bin/env python3
"""Check the C++ retention header (unsaturated.hpp) against the Python
reference: constitutive values, inverse, Bishop pressure, the bit-exact
saturated fast path, and the nodal content update for single- and
mixed-material patches. Compiles a small driver with g++, warnings as errors.
"""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'benchmarks' / 'poromechanics'))
from unsaturated_reference import Retention, advance_content, bishop_pressure

MODEL = {'van_genuchten': 1, 'brooks_corey': 2, 'gardner': 3}
LAWS = (Retention('van_genuchten', 1 / 5e3, 1.6, 0.1),
        Retention('van_genuchten', 1 / 2e4, 3.0, 0.0),
        Retention('brooks_corey', 1 / 3e3, 0.8, 0.05),
        Retention('gardner', 2 / 9810.0))
PHI = 0.4


def spec(r):
    return f'{MODEL[r.model]} {r.alpha:.17g} {(r.n or 0.0):.17g} {r.sr:.17g}'


class Header(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.exe = Path(cls.tmp.name) / 'unsaturated_driver'
        subprocess.run(['g++', '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wno-unknown-pragmas',
                        '-Werror', f'-I{ROOT}', str(Path(__file__).with_name('unsaturated_driver.cxx')),
                        '-o', str(cls.exe)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_driver(self, lines):
        out = subprocess.run([str(self.exe)], input='\n'.join(lines) + '\n', text=True,
                             capture_output=True, check=True).stdout
        return np.array([[float(v) for v in row.split()] for row in out.splitlines()])

    def test_constitutive_values(self):
        p = np.concatenate([-np.geomspace(1e-3, 1e8, 300), [0.0, 1e3, 1e6]])
        for r in LAWS:
            with self.subTest(r.model, n=r.n):
                got = self.run_driver([f'eval {spec(r)} {x:.17g}' for x in p])
                want = np.column_stack([r.se(p), r.dse(p), r.sat(p), r.kr(p), bishop_pressure(r, p)])
                np.testing.assert_allclose(got, want, rtol=2e-13, atol=1e-300)

    def test_saturated_model_is_inert(self):
        got = self.run_driver([f'eval 0 1e-4 2 0.1 {x:.17g}' for x in (-1e6, -1.0, 0.0, 5e4)])
        np.testing.assert_array_equal(got[:, :4], [[1, 0, 1, 1]] * 4)
        np.testing.assert_array_equal(got[:, 4], [-1e6, -1.0, 0.0, 5e4])

    def test_inverse(self):
        se = np.array([1e-8, 1e-3, 0.2, 0.7, 0.999, 1 - 1e-10])
        for r in LAWS:
            got = self.run_driver([f'pinv {spec(r)} {x:.17g}' for x in se])[:, 0]
            np.testing.assert_allclose(got, r.pressure_of_se(se), rtol=1e-12)
            np.testing.assert_allclose(r.se(got), se, rtol=1e-9)

    def test_saturated_fast_path_matches_existing_update(self):
        """dm = -(tdot*dt) gives exactly p - tdot*dt/hmass, as hydraulics.cxx."""
        rng = np.random.default_rng(7)
        hmass = 3e-10
        p = rng.uniform(1e4, 1e6, 300)
        tdot = rng.uniform(-1e-6, 1e-6, 300)
        dt = 0.37
        lines = [f'advance {hmass:.17g} {a:.17g} {-(b * dt):.17g} 2 {spec(LAWS[0])} 0.1 {spec(LAWS[2])} 0.05'
                 for a, b in zip(p, tdot)]
        got = self.run_driver(lines)
        np.testing.assert_array_equal(got[:, 0], p - tdot * dt / hmass)
        np.testing.assert_array_equal(got[:, 2], 0)

    def test_single_law_matches_reference(self):
        rng = np.random.default_rng(11)
        w, hmass = 0.25 * PHI, 0.25 * 1e-9
        for r in LAWS:
            with self.subTest(r.model, n=r.n):
                p = np.concatenate([-np.geomspace(1, 1e6, 200), rng.uniform(-2e4, 2e5, 200)])
                dm = rng.normal(0, 2.5e-4, p.size)
                got = self.run_driver([f'advance {hmass:.17g} {a:.17g} {b:.17g} 1 {spec(r)} {w:.17g}'
                                       for a, b in zip(p, dm)])
                self.assertLess(np.max(np.abs(got[:, 1] - dm)), 1e-15)
                np.testing.assert_array_equal(got[:, 2], 0)
                # Python solves the same root per unit volume.
                want = advance_content(r, PHI, 1e-9, p, dm / 0.25)
                np.testing.assert_allclose(got[:, 0], want, rtol=1e-9, atol=1e-6)

    def test_mixed_patch(self):
        """Different laws, weights and a saturated-model element at one node."""
        rng = np.random.default_rng(13)
        laws = [(LAWS[0], 0.03), (LAWS[2], 0.05), (LAWS[3], 0.02)]
        hmass = 2e-10
        tail = ' '.join(f'{spec(r)} {w:.17g}' for r, w in laws) + ' 0 1 2 0 0.04'
        def change(p, x):
            return hmass * (x - p) + sum(w * (r.sat(x) - r.sat(p)) for r, w in laws)
        for p in (-5e4, -2e3, -10.0, 0.0, 3e3):
            for dm in rng.normal(0, 5e-4, 20):
                x, c, bad = self.run_driver([f'advance {hmass:.17g} {p:.17g} {dm:.17g} 4 {tail}'])[0]
                self.assertLess(abs(c - dm), 1e-15)
                self.assertEqual(bad, 0)
                if p >= 0 and p + dm / hmass >= 0:  # whole patch stays saturated
                    self.assertEqual(x, p + dm / hmass)
                    continue
                lo, hi = sorted((p, p + dm / hmass))
                root = brentq(lambda y: change(p, y) - dm, lo, hi, xtol=1e-9, rtol=1e-15)
                self.assertAlmostEqual(x, root, delta=1e-6 * max(1.0, abs(root)))

    def test_absent_material_is_ignored(self):
        """A zero-weight entry (an absent material) changes neither result nor
        saturation threshold: Brooks-Corey here would mark p in (-1/alpha, 0)
        as saturated if it counted."""
        vg, bc = LAWS[0], LAWS[2]
        for p, dm in ((-1e3, 2e-4), (-200.0, -3e-4), (-5e4, 1e-3), (500.0, -1e-3)):
            one = self.run_driver([f'advance 1e-10 {p:.17g} {dm:.17g} 1 {spec(vg)} 0.1'])[0]
            two = self.run_driver([f'advance 1e-10 {p:.17g} {dm:.17g} 2 {spec(bc)} 0 {spec(vg)} 0.1'])[0]
            np.testing.assert_array_equal(one, two)

    def test_wetting_through_saturation_and_back(self):
        r = LAWS[0]
        w, hmass = 0.1, 1e-10
        dm = hmass * (5e3 + 3e4) + w * (r.sat(5e3) - r.sat(-3e4))
        x, _, _ = self.run_driver([f'advance {hmass:.17g} -30000 {dm:.17g} 1 {spec(r)} {w:.17g}'])[0]
        self.assertAlmostEqual(x, 5e3, delta=1e-6)
        back, _, _ = self.run_driver([f'advance {hmass:.17g} {x:.17g} {-dm:.17g} 1 {spec(r)} {w:.17g}'])[0]
        self.assertAlmostEqual(back, -3e4, delta=1e-5)

    def test_non_finite_input_is_flagged(self):
        r = LAWS[0]
        _, _, bad = self.run_driver([f'advance 1e-10 nan 1e-4 1 {spec(r)} 0.1'])[0]
        self.assertEqual(bad, 1)


if __name__ == '__main__':
    unittest.main()
