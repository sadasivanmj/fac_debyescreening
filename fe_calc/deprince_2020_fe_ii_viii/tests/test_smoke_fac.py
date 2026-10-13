"""Paper A: smoke tests that run FAC (skipped when the extension is not built).

Order follows the validation plan: Fe II unscreened with both orbital models, Fe II with
screening, Fe III with the same settings and the Fe II IP by independent subtraction.
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import fac_run_support as sup                        # noqa: E402
import testing_support as ts                         # noqa: E402

CALC = PAPER_DIR / "calculate_energies.py"


@unittest.skipIf(sup.find_extension() is None, "pfac extension not built (fe_calc/shared/build_pfac.sh)")
class SmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.res = Path(cls.tmp.name)
        for args in (("--ions", "fe_ii", "--models", "gref", "ref", "--mu", "0"),
                     ("--ions", "fe_ii", "--models", "gref", "--mu", "0.1"),
                     ("--ions", "fe_iii", "--models", "gref", "--mu", "0", "0.1")):
            proc = ts.run_script(CALC, "--skip-preflight", "--results-dir", str(cls.res), *args)
            assert proc.returncode == 0, proc.stdout + proc.stderr

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def manifest(self, ion, model, case):
        return json.loads((self.res / "runs" / ion / model / case / "manifest.json").read_text())

    def test_raw_log_and_outputs_exist_and_are_inspected(self):
        d = self.res / "runs" / "fe_ii" / "gref" / "mu0p00"
        for name in ("input.json", "run.log", "levels.txt", "final.pot", "manifest.json", "worker_result.json"):
            self.assertTrue((d / name).exists(), name)
        m = self.manifest("fe_ii", "gref", "mu0p00")
        self.assertEqual((m["status"], m["convergence_status"], m["exit_code"]), ("ok", "converged", 0))
        self.assertEqual(m["electron_count"], 25)
        self.assertEqual(m["n_levels"], 986)               # all J and both parities of the 7 configurations
        self.assertTrue(all(m["checks"].values()), m["checks"])

    def test_mu_zero_path_is_truly_unscreened(self):
        m = self.manifest("fe_ii", "gref", "mu0p00")
        self.assertTrue(m["checks"]["unscreened_path"])    # mps = -1: PlasmaScreen never called
        self.assertEqual(m["options"]["debye_mode"], 0)
        self.assertEqual(m["options"]["ee_screen"], 0)
        self.assertEqual(m["options"]["temperature_eV"], None)
        self.assertEqual(m["mu_steps"], [])
        self.assertAlmostEqual(m["asymptotic_charge"], -2.0, places=3)   # free ion: -(zion + 1)

    def test_nonzero_screening_applied_to_both_interactions(self):
        m = self.manifest("fe_ii", "gref", "mu0p10")
        self.assertAlmostEqual(m["mu_used"], 0.1, places=6)
        self.assertEqual((m["options"]["debye_mode"], m["options"]["ee_screen"], m["options"]["hxs_iter0"]), (1, 1, 8))
        self.assertTrue(m["checks"]["ee_direct_potential_screened"])     # asymptotic r*V -> 0
        self.assertAlmostEqual(m["asymptotic_charge"], 0.0, places=3)
        self.assertLess(abs(m["total_energy"] - self.manifest("fe_ii", "gref", "mu0p00")["total_energy"]), 1200)
        self.assertGreater(m["total_energy"], self.manifest("fe_ii", "gref", "mu0p00")["total_energy"])

    def test_ref_and_gref_distinct_with_same_configuration_space(self):
        a, b = self.manifest("fe_ii", "gref", "mu0p00"), self.manifest("fe_ii", "ref", "mu0p00")
        self.assertNotAlmostEqual(a["total_energy"], b["total_energy"], places=2)
        self.assertEqual(a["configurations"], b["configurations"])
        self.assertEqual(a["n_levels"], b["n_levels"])

    def test_fe_ii_ip_by_independent_subtraction(self):
        for case, expected in (("mu0p00", 15.5203), ("mu0p10", 10.2258)):     # validated gref IPs
            e2 = self.manifest("fe_ii", "gref", case)["total_energy"]
            e3 = self.manifest("fe_iii", "gref", case)["total_energy"]
            self.assertAlmostEqual(e3 - e2, expected, places=3, msg=case)
            ip, _ = __import__("analysis_support").ionization_potential(
                {**self.manifest("fe_ii", "gref", case)}, {**self.manifest("fe_iii", "gref", case)})
            self.assertAlmostEqual(ip, e3 - e2, places=6)


if __name__ == "__main__":
    unittest.main()
