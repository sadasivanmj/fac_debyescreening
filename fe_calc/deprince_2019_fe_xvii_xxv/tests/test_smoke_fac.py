"""Paper B: smoke tests that run FAC (skipped when the extension is not built).

Order follows the validation plan: Fe XVII unscreened, A and B; Fe XXVI unscreened and with
nuclear screening (with and without electron-electron screening); the Fe XXV IP from the
matching Fe XXVI energy.
"""

import csv
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
ANALYSIS = PAPER_DIR / "analyze_ionization_potentials.py"

# Independent check (point-nucleus Dirac 1s in V = -Z exp(-mu r)/r, solved separately from FAC):
# shift of the Fe XXVI ground-state energy relative to mu = 0, in eV.
INDEPENDENT_FE_XXVI_SHIFT = {"A": 70.5485, "C": 175.6218}


@unittest.skipIf(sup.find_extension() is None, "pfac extension not built (fe_calc/shared/build_pfac.sh)")
class SmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.res = Path(cls.tmp.name)
        for args in (("--ions", "fe_xvii", "--models", "ground_only", "reference_only", "--cases", "unscreened", "A", "B"),
                     ("--ions", "fe_xxv", "fe_xxvi", "--models", "ground_only")):
            proc = ts.run_script(CALC, "--skip-preflight", "--results-dir", str(cls.res), *args)
            assert proc.returncode == 0, proc.stdout + proc.stderr

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def m(self, ion, model, case):
        return json.loads((self.res / "runs" / ion / model / case / "manifest.json").read_text())

    def test_fe_xvii_runs_ok_with_raw_outputs(self):
        for case in ("unscreened", "A", "B"):
            d = self.res / "runs" / "fe_xvii" / "ground_only" / case
            for name in ("input.json", "run.log", "levels.txt", "final.pot", "manifest.json"):
                self.assertTrue((d / name).exists(), (case, name))
            m = self.m("fe_xvii", "ground_only", case)
            self.assertEqual((m["status"], m["convergence_status"], m["electron_count"], m["charge_state"]), ("ok", "converged", 10, 16))
            self.assertEqual(m["ground_level_identifier"].split(";")[:2], ["ILEV=0", "2J=0"])      # 1S0

    def test_unscreened_path_is_plain_fac(self):
        m = self.m("fe_xvii", "ground_only", "unscreened")
        self.assertTrue(m["checks"]["unscreened_path"])
        self.assertEqual((m["options"]["debye_mode"], m["options"]["ee_screen"], m["mu_steps"]), (0, 0, []))

    def test_nuclear_only_versus_fully_screened_physics(self):
        e0, ea, eb = (self.m("fe_xvii", "ground_only", c)["total_energy"] for c in ("unscreened", "A", "B"))
        self.assertLess(e0, eb)                      # nuclear screening raises the energy
        self.assertLess(eb, ea)                      # e-e screening weakens the repulsion: lower than nuclear-only
        self.assertLess(eb - ea, -10.0)
        a, b = self.m("fe_xvii", "ground_only", "A"), self.m("fe_xvii", "ground_only", "B")
        # e-e direct potential really screened: r*V -> 0 when both are screened, N-1 = 9 when only the nucleus is
        self.assertAlmostEqual(a["asymptotic_charge"], 9.0, places=2)
        self.assertAlmostEqual(b["asymptotic_charge"], 0.0, places=2)
        self.assertTrue(a["checks"]["ee_direct_potential_unscreened"])
        self.assertTrue(b["checks"]["ee_direct_potential_screened"])
        self.assertEqual((a["options"]["ee_screen"], b["options"]["ee_screen"]), (0, 1))
        self.assertAlmostEqual(a["mu_ne_used"], 0.1, places=6)
        self.assertAlmostEqual(b["mu_ne_used"], 0.1, places=6)

    def test_ground_only_and_reference_only_are_distinct_models(self):
        g, r = self.m("fe_xvii", "ground_only", "unscreened"), self.m("fe_xvii", "reference_only", "unscreened")
        self.assertEqual((len(g["configurations"]), len(r["configurations"])), (1, 5))
        self.assertLess(g["n_levels"], r["n_levels"])
        self.assertNotAlmostEqual(g["total_energy"], r["total_energy"], places=3)

    def test_fe_xxvi_responds_to_nuclear_screening_only(self):
        e0 = self.m("fe_xxvi", "ground_only", "unscreened")["total_energy"]
        for nuc, full, key in (("A", "B", "A"), ("C", "D", "C")):
            a, b = self.m("fe_xxvi", "ground_only", nuc), self.m("fe_xxvi", "ground_only", full)
            self.assertEqual((a["status"], b["status"]), ("ok", "ok"))
            self.assertTrue(a["checks"]["nuclear_screening_term_present"])
            self.assertEqual(a["total_energy"], b["total_energy"])            # mu_ee has no effect on one electron
            self.assertAlmostEqual(a["total_energy"] - e0, INDEPENDENT_FE_XXVI_SHIFT[key], delta=5e-3)
        self.assertNotEqual(self.m("fe_xxvi", "ground_only", "A")["total_energy"], e0)

    def test_fe_xxv_ip_from_matching_fe_xxvi_energy(self):
        with ts.poisoned_pfac_dir() as poison:
            out = self.res / "analysis"
            proc = ts.run_script(ANALYSIS, "--results-dir", str(self.res), "--output-dir", str(out), "--models", "ground_only",
                                 "--allow-incomplete", env_extra={"PYTHONPATH": poison})
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        with open(out / "ionization_potentials.csv") as f:
            rows = {(r["ion"], r["case"]): r for r in csv.DictReader(f)}
        for case in ("unscreened", "A", "B", "C", "D"):
            e26, e25 = (self.m(i, "ground_only", case)["total_energy"] for i in ("fe_xxvi", "fe_xxv"))
            self.assertAlmostEqual(float(rows[("Fe XXV", case)]["ip_calc_eV"]), e26 - e25, places=3, msg=case)
        self.assertEqual(rows[("Fe XXV", "A")]["status"], "ok")


if __name__ == "__main__":
    unittest.main()
