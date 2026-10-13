"""The independent reference evaluators of the potential validation (pure Python / numpy, no FAC)."""

import math
import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "shared"))
import potential_validation as pv      # noqa: E402

R = np.array([1e-3, 0.05, 0.4, 1.0, 3.0, 10.0, 40.0, 150.0])


class ReferenceEvaluators(unittest.TestCase):
    def test_density_normalisation(self):
        for name, shells in pv.DENSITIES.items():
            # integral of D dr by closed form: sum of N
            total = sum(n for n, _, _ in shells)
            r = np.geomspace(1e-9, 400.0, 400001)
            d = pv.density_values(shells, r)
            self.assertAlmostEqual(np.trapz(d, r), total, delta=1e-3 * total, msg=name)

    def test_closed_form_matches_direct_kernel_quadrature(self):
        for name in ("compact_1s", "valence", "diffuse"):
            for mu in (0.01, 0.25, 1.0):
                cf, _ = pv.closed_form_rV(pv.DENSITIES[name], R, mu)
                gl = pv.direct_rV(pv.DENSITIES[name], R, mu)
                scale = np.max(np.abs(cf))
                self.assertLess(np.max(np.abs(cf - gl)) / scale, 1e-10, (name, mu))

    def test_coulomb_limit_of_the_screened_closed_form(self):
        cf0, _ = pv.closed_form_rV(pv.DENSITIES["valence"], R, 0.0)
        cf_small, _ = pv.closed_form_rV(pv.DENSITIES["valence"], R, 1e-9)
        self.assertLess(np.max(np.abs(cf_small - cf0) / np.max(cf0)), 1e-6)

    def test_coulomb_asymptote_is_N_and_screened_decays(self):
        cf0, _ = pv.closed_form_rV(pv.DENSITIES["multishell"], np.array([300.0]), 0.0)
        self.assertAlmostEqual(cf0[0], 24.0, places=6)
        cf, _ = pv.closed_form_rV(pv.DENSITIES["multishell"], np.array([100.0]), 0.25)
        self.assertLess(cf[0], 1e-9)

    def test_jps_is_outer_integral_at_r0(self):
        r = np.array([1e-6])
        _, jps = pv.closed_form_rV(pv.DENSITIES["compact_1s"], r, 0.0)
        # int_r0^inf D/t dt for D = c t^2 e^{-2at}  ->  c/(2a)^2 = N (2a)^3/2 /(2a)^2 * ... = N*a*2 ... closed: N (2a)/2
        self.assertAlmostEqual(jps, 2.0 * 26.0 / 2.0, places=5)


if __name__ == "__main__":
    unittest.main()
