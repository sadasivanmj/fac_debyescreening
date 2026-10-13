"""The A/B analysis program on synthetic saved energies (no FAC, no calculations)."""

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parents[0] / "shared"))
import testing_support as ts      # noqa: E402

ANALYSIS = HERE / "analyze_densitytosz_experiment.py"
IONS = [("Fe II", 25), ("Fe III", 24), ("Fe IV", 23)]
CASES = [("mu0p00", 0.0), ("mu0p10", 0.1), ("mu0p25", 0.25)]


def energy(ion_i, mu, variant):
    # made-up energies: B differs from A by dE = -0.01 (Fe II), -0.03 (Fe III), -0.04 (Fe IV) times mu (so dE = 0 at mu = 0)
    base = -34000.0 - 20.0 * ion_i + 900.0 * mu * (1 + 0.1 * ion_i)
    return base + (0.0 if variant == "A" else -(0.01 + 0.01 * ion_i * ion_i) * mu * 10)


def records(variant, skip=(), control_error=0.0):
    out = []
    for i, (ion, n) in enumerate(IONS):
        for k, (case, mu) in enumerate(CASES):
            if (ion, case) in skip:
                continue
            e = energy(i, mu, variant) + (control_error if (mu == 0 and variant == "A") else 0.0)
            out.append({"ion": ion, "ion_order": i, "electron_count": n, "configuration_model": "reference_configurations",
                        "orbital_optimization_model": "gref", "case": case, "case_order": k, "mu_ne": mu, "mu_ee": mu,
                        "total_energy": e, "direct_potential": "coulomb" if variant == "A" else "screened",
                        "build_provenance": f"build-{variant}", "fac_source_commit": f"commit-{variant}"})
    return out


def run(skip_A=(), control_error=0.0, *extra):
    tmp = tempfile.TemporaryDirectory()
    res = Path(tmp.name)
    ts.write_synthetic_results(res / "variant_A", records("A", skip_A, control_error))
    ts.write_synthetic_results(res / "variant_B", records("B"))
    with ts.poisoned_pfac_dir() as poison:
        p = ts.run_script(ANALYSIS, "--results-dir", str(res), *extra, env_extra={"PYTHONPATH": poison})
    return tmp, res, p


class AnalysisTests(unittest.TestCase):
    def test_complete_dataset_identities_and_no_fac_import(self):
        tmp, res, p = run()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rows = list(csv.DictReader(open(res / "ip_comparison.csv")))
        self.assertEqual(len(rows), 6)
        for r in rows:
            self.assertAlmostEqual(float(r["delta_I_BminusA_eV"]), float(r["delta_I_from_dE_eV"]), places=9)
            self.assertLess(float(r["identity_residual_eV"]), 1e-8)
        lo = next(r for r in rows if r["ionization"] == "Fe II -> Fe III" and r["case"] == "mu0p10")
        eB = {i: energy(i, 0.1, "B") for i in range(3)}
        eA = {i: energy(i, 0.1, "A") for i in range(3)}
        self.assertAlmostEqual(float(lo["IP_B_eV"]), eB[1] - eB[0], places=5)
        self.assertAlmostEqual(float(lo["delta_I_BminusA_eV"]), (eB[1] - eB[0]) - (eA[1] - eA[0]), places=5)
        tmp.cleanup()

    def test_variant_a_failures_reported_not_filled(self):
        tmp, res, p = run(skip_A={("Fe III", "mu0p25"), ("Fe II", "mu0p25")})
        self.assertEqual(p.returncode, 2)
        self.assertIn("missing energy", p.stdout)
        rows = {(r["ionization"], r["case"]): r for r in csv.DictReader(open(res / "ip_comparison.csv"))}
        self.assertEqual(rows[("Fe II -> Fe III", "mu0p25")]["IP_A_eV"], "")
        self.assertEqual(rows[("Fe II -> Fe III", "mu0p25")]["delta_I_BminusA_eV"], "")
        self.assertEqual(rows[("Fe II -> Fe III", "mu0p25")]["status"], "B only (A did not converge)")
        tmp.cleanup()

    def test_allow_incomplete_exit_zero(self):
        tmp, res, p = run(skip_A={("Fe IV", "mu0p10")}, extra=None) if False else run(skip_A={("Fe IV", "mu0p10")})
        self.assertEqual(p.returncode, 2)
        tmp.cleanup()
        tmp, res, p = run({("Fe IV", "mu0p10")}, 0.0, "--allow-incomplete")
        self.assertEqual(p.returncode, 0)
        tmp.cleanup()

    def test_unscreened_control_must_agree(self):
        tmp, res, p = run(control_error=0.01)
        self.assertEqual(p.returncode, 2)
        self.assertIn("unscreened control differs", p.stdout)
        tmp.cleanup()


if __name__ == "__main__":
    unittest.main()
