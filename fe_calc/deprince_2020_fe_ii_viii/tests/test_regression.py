"""Paper A: comparison of the saved energies with earlier validated results (needs results/energies.csv)."""

import csv
import sys
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "linux_validated_fe_ii_iv_energies.csv"
ENERGIES = PAPER_DIR / "results" / "energies.csv"
CASE_OF_MU = {"0.0": "mu0p00", "0.1": "mu0p10", "0.25": "mu0p25"}


def load_energies():
    with open(ENERGIES, newline="") as f:
        return {(r["ion"], r["orbital_optimization_model"], r["case"]): r for r in csv.DictReader(f)}


class RegressionAgainstValidatedResults(unittest.TestCase):
    def setUp(self):
        if not ENERGIES.exists():
            self.skipTest("results/energies.csv not present; run calculate_energies.py")
        self.new = load_energies()

    def test_fe_ii_to_iv_match_linux_validated_energies(self):
        with open(FIXTURE, newline="") as f:
            old = list(csv.DictReader(f))
        self.assertEqual(len(old), 18)
        worst = 0.0
        for r in old:
            key = (r["ion"], r["model"], CASE_OF_MU[r["mu_requested"]])
            self.assertIn(key, self.new, key)
            self.assertEqual(self.new[key]["status"], "ok", key)
            worst = max(worst, abs(float(self.new[key]["total_energy"]) - float(r["E0_eV"])))
        # same C patches and sequence; platform/flags differ, so a tolerance, not bit-identity
        self.assertLess(worst, 1e-3, f"max |dE| = {worst} eV")

    def test_ref_and_gref_energies_are_distinct_with_same_ci_space(self):
        for ion in cfg.ION_ORDER:
            for case in cfg.CASES:
                a, b = self.new[(ion, "gref", case)], self.new[(ion, "ref", case)]
                self.assertNotEqual(a["total_energy"], b["total_energy"], (ion, case))

    def test_all_48_energies_ok_and_finite(self):
        self.assertEqual(len(self.new), 48)
        for key, r in self.new.items():
            self.assertEqual(r["status"], "ok", key)
            self.assertEqual(r["convergence_status"], "converged", key)
            self.assertLess(float(r["total_energy"]), 0.0, key)

    def test_single_build_provenance_across_all_runs(self):
        self.assertEqual(len({r["build_provenance"] for r in self.new.values()}), 1)
        self.assertEqual(len({r["fac_source_commit"] for r in self.new.values()}), 1)


if __name__ == "__main__":
    unittest.main()
