"""Paper B: comparison of the saved energies with earlier validated results and physical sanity
checks over the full dataset (needs results/energies.csv)."""

import csv
import sys
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "validated_ground_only_energies_cases_A_D.csv"
ENERGIES = PAPER_DIR / "results" / "energies.csv"


class RegressionAndSanity(unittest.TestCase):
    def setUp(self):
        if not ENERGIES.exists():
            self.skipTest("results/energies.csv not present; run calculate_energies.py")
        with open(ENERGIES, newline="") as f:
            self.rows = list(csv.DictReader(f))
        self.e = {(r["ion"], r["configuration_model"], r["case"]): float(r["total_energy"]) for r in self.rows if r["status"] == "ok"}

    def test_ground_only_matches_earlier_validated_energies(self):
        with open(FIXTURE, newline="") as f:
            old = list(csv.DictReader(f))
        self.assertEqual(len(old), 40)
        worst = 0.0
        for r in old:
            key = (r["ion"], "ground_only", r["case"])
            self.assertIn(key, self.e, key)
            worst = max(worst, abs(self.e[key] - float(r["total_energy_eV"])))
        self.assertLess(worst, 1e-3, f"max |dE| = {worst} eV")      # same C code/sequence, different compiler flags

    def test_every_run_ok_converged_and_single_build(self):
        ok = [r for r in self.rows if r["status"] == "ok" and r["convergence_status"] == "converged"]
        self.assertEqual(len(ok), len(self.rows))
        self.assertEqual(len({r["build_provenance"] for r in self.rows}), 1)
        self.assertEqual(len({r["fac_source_commit"] for r in self.rows}), 1)
        self.assertGreaterEqual(len(self.rows), 50)

    def test_screening_ordering_for_every_ion(self):
        for model in {m for _, m, _ in self.e}:
            for ion in cfg.ION_ORDER:
                u, a, b, c, d = (self.e[(ion, model, k)] for k in ("unscreened", "A", "B", "C", "D"))
                self.assertLess(u, a, (ion, model))        # nuclear screening raises the energy
                self.assertLess(a, c, (ion, model))        # stronger screening raises it more
                if cfg.IONS[ion][0] > 1:
                    self.assertLess(b, a, (ion, model))    # e-e screening lowers it
                    self.assertLess(d, c, (ion, model))
                else:                                       # one electron: no e-e interaction
                    self.assertEqual((a, c), (b, d))

    def test_ionization_potentials_positive_and_decrease_with_screening(self):
        for model in {m for _, m, _ in self.e}:
            for k, ion in enumerate(cfg.ION_ORDER[:-1]):
                nxt = cfg.ION_ORDER[k + 1]
                ip = {c: self.e[(nxt, model, c)] - self.e[(ion, model, c)] for c in cfg.CASES}
                self.assertTrue(all(v > 0 for v in ip.values()), (ion, model))
                self.assertGreater(ip["unscreened"], ip["A"])
                self.assertGreater(ip["B"], ip["A"])           # e-e screening increases the IP relative to nuclear-only
                self.assertGreater(ip["A"], ip["C"])

    def test_reference_only_lowers_energy_and_is_unchanged_only_where_symmetry_forbids_mixing(self):
        models = {m for _, m, _ in self.e}
        if models != {"ground_only", "reference_only"}:
            self.skipTest("reference_only not in the dataset")
        # Same orbitals (optimised on the ground configuration); adding configurations to the CI can only
        # lower the lowest level.  Identical energies are expected where the extra references cannot mix
        # with the ground level for J / parity reasons (Fe XVIII 2P3/2, Fe XIX 3P2) and for one electron.
        identical = set()
        for ion in cfg.ION_ORDER:
            for case in cfg.CASES:
                d = self.e[(ion, "reference_only", case)] - self.e[(ion, "ground_only", case)]
                self.assertLessEqual(d, 1e-6, (ion, case))
                if abs(d) < 1e-6:
                    identical.add(ion)
        self.assertEqual(identical, {"Fe XVIII", "Fe XIX", "Fe XXVI"})
        self.assertLess(self.e[("Fe XXIII", "reference_only", "unscreened")] - self.e[("Fe XXIII", "ground_only", "unscreened")], -5.0)   # 2s2-2p2 mixing

if __name__ == "__main__":
    unittest.main()
