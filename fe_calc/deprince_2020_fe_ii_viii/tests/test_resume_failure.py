"""Paper A: resumability, no silent overwrite, timeout handling (runs FAC; skipped without a build)."""

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
ARGS = ("--skip-preflight", "--ions", "fe_ix", "--models", "gref")


@unittest.skipIf(sup.find_extension() is None, "pfac extension not built")
class ResumeAndFailure(unittest.TestCase):
    def test_resume_overwrite_timeout_and_retry(self):
        with tempfile.TemporaryDirectory() as d:
            res = Path(d)
            first = ts.run_script(CALC, *ARGS, "--mu", "0", "--results-dir", str(res))
            self.assertEqual(first.returncode, 0, first.stdout)
            run = res / "runs" / "fe_ix" / "gref" / "mu0p00"
            stamp = (run / "manifest.json").stat().st_mtime_ns
            e0 = json.loads((run / "manifest.json").read_text())["total_energy"]

            again = ts.run_script(CALC, *ARGS, "--mu", "0", "--results-dir", str(res))      # resume: nothing is redone
            self.assertIn("skipped (already ok)", again.stdout)
            self.assertEqual((run / "manifest.json").stat().st_mtime_ns, stamp)

            over = ts.run_script(CALC, *ARGS, "--mu", "0", "--results-dir", str(res), "--overwrite")
            self.assertEqual(over.returncode, 0)
            moved = list(run.parent.glob("mu0p00.superseded-*"))                             # old result moved aside, not deleted
            self.assertEqual(len(moved), 1)
            self.assertAlmostEqual(json.loads((moved[0] / "manifest.json").read_text())["total_energy"], e0, places=6)

            fail = ts.run_script(CALC, *ARGS, "--mu", "0.1", "--results-dir", str(res), "--timeout", "0.05")
            self.assertEqual(fail.returncode, 1)
            bad = res / "runs" / "fe_ix" / "gref" / "mu0p10"
            m = json.loads((bad / "manifest.json").read_text())
            self.assertEqual(m["status"], "timeout")
            self.assertIsNone(m["total_energy"])
            self.assertTrue((bad / "run.FAILED.log").exists())                               # diagnostics kept under a trackable name
            self.assertIn("timeout", (res / "energies.csv").read_text())

            retry = ts.run_script(CALC, *ARGS, "--mu", "0.1", "--results-dir", str(res))
            self.assertEqual(retry.returncode, 0, retry.stdout)
            self.assertEqual(len(list(bad.parent.glob("mu0p10.failed-attempt-*"))), 1)       # failed attempt preserved
            self.assertEqual(json.loads((bad / "manifest.json").read_text())["status"], "ok")


class PreflightAndCheckOnly(unittest.TestCase):
    def test_check_only_passes_without_fac(self):
        p = ts.run_script(CALC, "--check-only")
        self.assertEqual(p.returncode, 0, p.stdout)
        self.assertIn("configuration check: passed", p.stdout)

    @unittest.skipIf(sup.find_extension() is None, "pfac extension not built")
    def test_preflight_reports_provenance_and_workaround(self):
        p = ts.run_script(CALC, "--preflight")
        self.assertEqual(p.returncode, 0, p.stdout)
        for text in ("stable Debye DensityToSZ (v2) compiled in", "Fe XXVI one-electron nuclear-screening fix",
                     "jmaxk uninitialised-variable defect: source defect PRESENT (not fixed)",
                     "BUILD_INFO matches extension"):
            self.assertIn(text, p.stdout)


if __name__ == "__main__":
    unittest.main()
