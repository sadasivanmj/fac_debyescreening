"""Paper B: the IP-analysis program on synthetic saved energies (no FAC, no calculations)."""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402
import reference_data as ref                         # noqa: E402
import testing_support as ts                         # noqa: E402

ANALYSIS = PAPER_DIR / "analyze_ionization_potentials.py"


def synthetic_energy(ion, model, case):
    n = cfg.IONS[ion][0]
    mu_ne, mu_ee = cfg.CASES[case]
    return (-1500.0 * n + 2.5 * cfg.ION_ORDER.index(ion) ** 2 + 90.0 * mu_ne * (26 - n + 1) * (1 if n else 0)
            - 40.0 * mu_ee * n + (0.7 if model == "reference_only" else 0.0))


def records(models=("ground_only",), skip=(), override=None):
    out = []
    for ion in cfg.ION_ORDER:
        for model in models:
            for ci, (case, (mu_ne, mu_ee)) in enumerate(cfg.CASES.items()):
                if (ion, model, case) in skip:
                    continue
                r = {"ion": ion, "ion_order": cfg.ION_ORDER.index(ion), "electron_count": cfg.IONS[ion][0],
                     "configuration_model": model, "orbital_optimization_model": cfg.ORBITAL_OPTIMIZATION,
                     "case": case, "case_order": ci, "mu_ne": mu_ne, "mu_ee": mu_ee,
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
        self.assertEqual(len(rows), 45)                                 # 9 ionizations x 5 cases
        for r in rows:
            up = cfg.ION_ORDER[cfg.ION_ORDER.index(r["ion"]) + 1]
            expected = synthetic_energy(up, "ground_only", r["case"]) - synthetic_energy(r["ion"], "ground_only", r["case"])
            self.assertAlmostEqual(float(r["ip_calc_eV"]), expected, places=3)
        self.assertIn("Fe XXV -> Fe XXVI", {r["ionization"] for r in rows})
        self.assertNotIn("Fe XXVI -> ", " ".join(r["ionization"] for r in rows))   # Fe XXVI is not a Table row
        self.assertTrue((res / "comparison_summary.md").exists())

    def test_fe_xxv_ip_uses_matching_fe_xxvi_case(self):
        _, rows, _ = self.run_analysis(records())
        for case in cfg.CASES:
            row = next(r for r in rows if r["ion"] == "Fe XXV" and r["case"] == case)
            self.assertAlmostEqual(float(row["energy_upper_eV"]), synthetic_energy("Fe XXVI", "ground_only", case), places=6)
            self.assertAlmostEqual(float(row["energy_lower_eV"]), synthetic_energy("Fe XXV", "ground_only", case), places=6)

    def test_ee_effect_and_shifts_and_published_columns(self):
        _, rows, _ = self.run_analysis(records())
        ip = {(r["ion"], r["case"]): float(r["ip_calc_eV"]) for r in rows}
        for r in rows:
            if r["case"] in ("B", "D"):
                partner = {"B": "A", "D": "C"}[r["case"]]
                self.assertAlmostEqual(float(r["ee_effect_calc_eV"]), ip[(r["ion"], r["case"])] - ip[(r["ion"], partner)], places=3)
                self.assertAlmostEqual(float(r["ee_effect_published_eV"]),
                                       ref.TABLE4[r["ion"]][r["case"]] - ref.TABLE4[r["ion"]][partner], places=3)
            if r["case"] == "unscreened":
                self.assertEqual(float(r["ip_published_eV"]), ref.TABLE2_MU0[r["ion"]])
                self.assertEqual(float(r["ip_nist_eV"]), ref.NIST_IP[r["ion"]])
            else:
                self.assertEqual(float(r["ip_published_eV"]), ref.TABLE4[r["ion"]][r["case"]])
                self.assertAlmostEqual(float(r["shift_vs_unscreened_calc_eV"]), ip[(r["ion"], r["case"])] - ip[(r["ion"], "unscreened")], places=3)

    def test_missing_energy_fails_loudly_and_models_not_substituted(self):
        recs = records(models=("ground_only", "reference_only"), skip={("Fe XXVI", "ground_only", "B")})
        proc, rows, _ = self.run_analysis(recs)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("missing energy", proc.stdout)
        bad = [r for r in rows if r["status"] != "ok"]
        self.assertEqual([(r["configuration_model"], r["case"], r["ion"]) for r in bad], [("ground_only", "B", "Fe XXV")])
        self.assertEqual(bad[0]["ip_calc_eV"], "")                       # not filled from reference_only

    def test_incompatible_models_or_builds_block_the_ip(self):
        recs = records(override={("Fe XXIV", "ground_only", "C"): {"build_provenance": "another-build"}})
        proc, rows, _ = self.run_analysis(recs)
        self.assertEqual(proc.returncode, 2)
        bad = [r for r in rows if r["status"] != "ok"]
        self.assertEqual(len(bad), 2)
        self.assertTrue(all("incompatible build_provenance" in r["status"] for r in bad))

    def test_failed_and_duplicated_records_reported(self):
        recs = records(override={("Fe XX", "ground_only", "A"): {"status": "timeout", "convergence_status": "unknown"}})
        recs.append(dict(recs[3]))
        proc, _, _ = self.run_analysis(recs)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("invalid run", proc.stdout)
        self.assertIn("duplicated record", proc.stdout)

    def test_models_are_reported_separately(self):
        _, rows, _ = self.run_analysis(records(models=("ground_only", "reference_only")))
        self.assertEqual(len(rows), 90)
        a = {(r["ion"], r["case"]): r["ip_calc_eV"] for r in rows if r["configuration_model"] == "ground_only"}
        b = {(r["ion"], r["case"]): r["ip_calc_eV"] for r in rows if r["configuration_model"] == "reference_only"}
        self.assertEqual(set(a), set(b))


if __name__ == "__main__":
    unittest.main()
