"""Paper A: the IP-analysis program on synthetic saved energies (no FAC, no calculations)."""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402
import testing_support as ts                         # noqa: E402

ANALYSIS = PAPER_DIR / "analyze_ionization_potentials.py"


def synthetic_energy(ion, model, case):
    """Deterministic made-up total energies (eV); the IP is their difference."""
    n = cfg.IONS[ion][0]
    return -1000.0 * n + 3.25 * (cfg.ION_ORDER.index(ion) ** 2) - 7.5 * cfg.CASES[case] * (26 - n) + (0.5 if model == "ref" else 0.0)


def records(skip=(), override=None):
    out = []
    for ion in cfg.ION_ORDER:
        for model in cfg.ORBITAL_MODELS:
            for ci, (case, mu) in enumerate(cfg.CASES.items()):
                if (ion, model, case) in skip:
                    continue
                r = {"ion": ion, "ion_order": cfg.ION_ORDER.index(ion), "electron_count": cfg.IONS[ion][0],
                     "configuration_model": cfg.CONFIGURATION_MODEL, "orbital_optimization_model": model,
                     "case": case, "case_order": ci, "mu_ne": mu, "mu_ee": mu,
                     "total_energy": synthetic_energy(ion, model, case)}
                r.update((override or {}).get((ion, model, case), {}))
                out.append(r)
    return out


class AnalysisOnSavedEnergies(unittest.TestCase):
    def run_analysis(self, recs, *extra):
        self.tmp = tempfile.TemporaryDirectory()
        res = Path(self.tmp.name) / "results"
        ts.write_synthetic_results(res, recs)
        with ts.poisoned_pfac_dir() as poison:       # any import of pfac would raise
            proc = ts.run_script(ANALYSIS, "--results-dir", str(res), *extra, env_extra={"PYTHONPATH": poison})
        rows = []
        if (res / "ionization_potentials.csv").exists():
            with open(res / "ionization_potentials.csv") as f:
                rows = list(csv.DictReader(f))
        return proc, rows, res

    def tearDown(self):
        getattr(self, "tmp", None) and self.tmp.cleanup()

    def test_complete_dataset_independent_recomputation_and_no_fac_import(self):
        proc, rows, res = self.run_analysis(records())
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual(len(rows), 42)
        for r in rows:
            lo, up = r["ion"], cfg.ION_ORDER[cfg.ION_ORDER.index(r["ion"]) + 1]
            expected = synthetic_energy(up, r["orbital_optimization_model"], r["case"]) - synthetic_energy(lo, r["orbital_optimization_model"], r["case"])
            self.assertAlmostEqual(float(r["ip_calc_eV"]), expected, places=3)
        self.assertTrue((res / "comparison_summary.md").exists())

    def test_ip_pairs_adjacent_ions_only(self):
        _, rows, _ = self.run_analysis(records())
        self.assertEqual([r["ionization"] for r in rows[:1]], ["Fe II -> Fe III"])
        self.assertEqual(rows[-1]["ionization"], "Fe VIII -> Fe IX")

    def test_missing_energy_fails_loudly_and_is_not_replaced(self):
        # remove Fe III gref mu=0.1: the gref IPs involving it must NOT be filled from the ref model
        proc, rows, _ = self.run_analysis(records(skip={("Fe III", "gref", "mu0p10")}))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("missing energy", proc.stdout)
        bad = [r for r in rows if r["status"] != "ok"]
        self.assertEqual({(r["ionization"], r["orbital_optimization_model"], r["case"]) for r in bad},
                         {("Fe II -> Fe III", "gref", "mu0p10"), ("Fe III -> Fe IV", "gref", "mu0p10")})
        self.assertTrue(all(r["ip_calc_eV"] == "" for r in bad))

    def test_incompatible_source_commit_blocks_the_ip(self):
        proc, rows, _ = self.run_analysis(records(override={("Fe III", "gref", "mu0p25"): {"fac_source_commit": "b" * 40}}))
        self.assertEqual(proc.returncode, 2)
        bad = [r for r in rows if r["status"] != "ok"]
        self.assertEqual(len(bad), 2)
        self.assertTrue(all("incompatible fac_source_commit" in r["status"] for r in bad))

    def test_failed_and_duplicated_records_reported(self):
        recs = records(override={("Fe V", "ref", "mu0p00"): {"status": "not_converged", "convergence_status": "not_converged"}})
        recs.append(dict(recs[0]))
        proc, rows, _ = self.run_analysis(recs)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("invalid run", proc.stdout)
        self.assertIn("duplicated record", proc.stdout)

    def test_allow_incomplete_gives_exit_zero_but_still_reports(self):
        proc, _, _ = self.run_analysis(records(skip={("Fe II", "ref", "mu0p00")}), "--allow-incomplete")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("PROBLEM", proc.stdout)

    def test_models_never_mixed(self):
        _, rows, _ = self.run_analysis(records())
        gref = {(r["ion"], r["case"]): float(r["ip_calc_eV"]) for r in rows if r["orbital_optimization_model"] == "gref"}
        ref = {(r["ion"], r["case"]): float(r["ip_calc_eV"]) for r in rows if r["orbital_optimization_model"] == "ref"}
        self.assertEqual(set(gref), set(ref))        # separate rows for each model


if __name__ == "__main__":
    unittest.main()
