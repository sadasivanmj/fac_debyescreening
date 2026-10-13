"""Tests of the shared modules (no FAC needed)."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import analysis_support as an                       # noqa: E402
import config_validation as cv                      # noqa: E402
import fac_run_support as sup                       # noqa: E402

LEVEL_TABLE = """FAC 1.1.5
E0\t= 0, -2.55833916E+04

NELE\t= 6
NLEV\t= 2
  ILEV  IBASE    ENERGY       P   VNL         2J
     0     -1  0.00000000E+00 0   201          0 1*2.2*4                          2p2
     1     -1  8.78461239E+00 0   201          2 1*2.2*4                          2p2
"""


class ConfigValidation(unittest.TestCase):
    def test_valid_configuration(self):
        self.assertEqual(cv.electron_count("1s2 2s2 2p6 3d1"), 11)

    def test_bad_syntax_occupancy_and_shell_rejected(self):
        for bad in ("1s3", "2p7", "1p1", "1s2 1s1", "xx", "", "2s0"):
            with self.assertRaises(ValueError, msg=bad):
                cv.parse_configuration(bad)

    def test_wrong_electron_count_and_duplicates_reported(self):
        problems = cv.validate_configuration_set("X", 3, ["1s2 2s1", "2s1 1s2", "1s2 2p2"])
        self.assertTrue(any("duplicate" in p for p in problems))
        self.assertTrue(any("has 4 electrons" in p for p in problems))


class FacOutputParsing(unittest.TestCase):
    def test_level_table(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "levels.txt"
            p.write_text(LEVEL_TABLE)
            t = sup.parse_level_table(p)
        self.assertEqual((t["nele"], t["nlev_header"], len(t["levels"])), (6, 2, 2))
        self.assertAlmostEqual(t["e0_ev"], -25583.3916)
        g = sup.ground_level(t)
        self.assertEqual((g["ilev"], g["two_j"]), (0, 0))

    def test_log_scan(self):
        s = sup.scan_log("fine\nMaximum iteration reached in OptimizeRadial\n")
        self.assertTrue(s["max_iteration_reached"])
        self.assertFalse(sup.scan_log("all good\n")["max_iteration_reached"])


class IpDefinition(unittest.TestCase):
    def rec(self, ion, n, e, **kw):
        r = {"ion": ion, "electron_count": n, "total_energy": e, "configuration_model": "m",
             "orbital_optimization_model": "o", "case": "c", "mu_ne": 0.0, "mu_ee": 0.0,
             "energy_unit": "eV", "fac_source_commit": "x", "build_provenance": "b"}
        r.update(kw)
        return r

    def test_ip_is_energy_difference_in_ev(self):
        ip, why = an.ionization_potential(self.rec("A", 3, -100.0), self.rec("B", 2, -60.5))
        self.assertAlmostEqual(ip, 39.5, places=9)

    def test_incompatible_records_refused(self):
        for field, value in (("orbital_optimization_model", "other"), ("case", "d"), ("mu_ee", 0.1),
                             ("fac_source_commit", "y"), ("build_provenance", "other-build"),
                             ("configuration_model", "n")):
            ip, why = an.ionization_potential(self.rec("A", 3, -100.0), self.rec("B", 2, -60.0, **{field: value}))
            self.assertIsNone(ip, field)
            self.assertIn(field, why)

    def test_non_adjacent_charge_states_refused(self):
        ip, why = an.ionization_potential(self.rec("A", 4, -100.0), self.rec("C", 2, -60.0))
        self.assertIsNone(ip)


class RunDirectories(unittest.TestCase):
    def test_archive_moves_and_keeps_content(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d) / "run"
            run.mkdir()
            (run / "keep.txt").write_text("x")
            moved = sup.archive_previous(run, "superseded")
            self.assertFalse(run.exists())
            self.assertEqual((moved / "keep.txt").read_text(), "x")


if __name__ == "__main__":
    unittest.main()
