"""Paper B: resumability, no silent overwrite, timeout handling (runs FAC; skipped without a build)."""

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
ARGS = ("--skip-preflight", "--ions", "fe_xxvi", "--models", "ground_only")


@unittest.skipIf(sup.find_extension() is None, "pfac extension not built")
class ResumeAndFailure(unittest.TestCase):
    def test_resume_overwrite_timeout_and_retry(self):
        with tempfile.TemporaryDirectory() as d:
            res = Path(d)
            first = ts.run_script(CALC, *ARGS, "--cases", "unscreened", "--results-dir", str(res))
            self.assertEqual(first.returncode, 0, first.stdout)
            run = res / "runs" / "fe_xxvi" / "ground_only" / "unscreened"
            stamp = (run / "manifest.json").stat().st_mtime_ns
            e0 = json.loads((run / "manifest.json").read_text())["total_energy"]

            again = ts.run_script(CALC, *ARGS, "--cases", "unscreened", "--results-dir", str(res))
            self.assertIn("skipped (already ok)", again.stdout)
            self.assertEqual((run / "manifest.json").stat().st_mtime_ns, stamp)

            over = ts.run_script(CALC, *ARGS, "--cases", "unscreened", "--results-dir", str(res), "--overwrite")
            self.assertEqual(over.returncode, 0)
            moved = list(run.parent.glob("unscreened.superseded-*"))
            self.assertEqual(len(moved), 1)
            self.assertAlmostEqual(json.loads((moved[0] / "manifest.json").read_text())["total_energy"], e0, places=6)

            fail = ts.run_script(CALC, *ARGS, "--cases", "A", "--results-dir", str(res), "--timeout", "0.05")
            self.assertEqual(fail.returncode, 1)
            bad = res / "runs" / "fe_xxvi" / "ground_only" / "A"
            m = json.loads((bad / "manifest.json").read_text())
            self.assertEqual(m["status"], "timeout")
            self.assertIsNone(m["total_energy"])
            self.assertTrue((bad / "run.FAILED.log").exists())

            retry = ts.run_script(CALC, *ARGS, "--cases", "A", "--results-dir", str(res))
            self.assertEqual(retry.returncode, 0, retry.stdout)
            self.assertEqual(len(list(bad.parent.glob("A.failed-attempt-*"))), 1)
            self.assertEqual(json.loads((bad / "manifest.json").read_text())["status"], "ok")


class PlanningAndPreflight(unittest.TestCase):
    def test_check_only_passes_without_fac(self):
        p = ts.run_script(CALC, "--check-only")
        self.assertEqual(p.returncode, 0, p.stdout)
        self.assertIn("configuration check: passed", p.stdout)

    def test_selection_of_ions_cases_models(self):
        p = ts.run_script(CALC, "--skip-preflight", "--dry-run", "--ions", "fe_xvii", "fe_xxvi",
                          "--cases", "unscreened", "B", "--models", "ground_only", "reference_only")
        self.assertEqual(p.returncode, 0, p.stdout)
        self.assertIn("8 planned run(s)", p.stdout)

    @unittest.skipIf(sup.find_extension() is None, "pfac extension not built")
    def test_preflight_reports_provenance_and_workaround(self):
        p = ts.run_script(CALC, "--preflight")
        self.assertEqual(p.returncode, 0, p.stdout)
        for text in ("stable Debye DensityToSZ (v2) compiled in", "Fe XXVI one-electron nuclear-screening fix",
                     "jmaxk uninitialised-variable defect: source defect PRESENT (not fixed)"):
            self.assertIn(text, p.stdout)


if __name__ == "__main__":
    unittest.main()
