"""Paper A: configuration data (no FAC needed)."""

import re
import sys
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
REPO = PAPER_DIR.parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402
import validate_inputs                               # noqa: E402
from config_validation import canonical, electron_count   # noqa: E402


class Configurations(unittest.TestCase):
    def test_input_validation_passes(self):
        self.assertEqual(validate_inputs.validate(), [])

    def test_electron_counts_and_charge_states(self):
        for ion, n in zip(cfg.ION_ORDER, range(25, 17, -1)):
            electrons, configs, _ = cfg.IONS[ion]
            self.assertEqual(electrons, n)
            self.assertTrue(all(electron_count(c) == n for c in configs), ion)
            self.assertEqual(cfg.charge_state(ion), 26 - n)

    def test_ground_configuration_is_first_and_has_no_hole(self):
        for ion, (n, configs, _) in cfg.IONS.items():
            self.assertTrue(configs[0].startswith("1s2 2s2 2p6 3s2 3p6"), ion)

    def test_hole_notation_translated(self):
        fe3 = cfg.IONS["Fe III"][1]
        self.assertIn("1s2 2s2 2p6 3s2 3p5 3d7", fe3)          # [3p]3d^7
        self.assertIn("1s2 2s2 2p5 3s2 3p6 3d5 4s2", fe3)       # [2p]3d^5 4s^2
        self.assertIn("1s1 2s2 2p6 3s2 3p6 3d6 4s1", fe3)       # [1s]3d^6 4s
        self.assertEqual(len(fe3), 10)
        self.assertEqual([len(cfg.IONS[i][1]) for i in cfg.ION_ORDER], [7, 10, 4, 4, 4, 4, 4, 4])

    def test_fe_ix_is_paper_iii_reference_set(self):
        self.assertEqual(cfg.IONS["Fe IX"][1],
                         ["1s2 2s2 2p6 3s2 3p6", "1s2 2s2 2p6 3s2 3p5 3d1",
                          "1s2 2s2 2p5 3s2 3p6 3d1", "1s1 2s2 2p6 3s2 3p6 3d1"])

    def test_no_generated_excitations(self):
        # reduced model: exactly the listed reference configurations, nothing more
        self.assertLessEqual(max(len(c) for _, c, _ in cfg.IONS.values()), 10)

    def test_fe_ii_to_iv_equal_validated_lists(self):
        base = REPO / "research" / "debye-screening" / "cases" / "fe"
        if not base.exists():
            self.skipTest("validated Fe II-IV scripts not present")
        for ion, name in (("Fe II", "fe_ii"), ("Fe III", "fe_iii"), ("Fe IV", "fe_iv")):
            text = (base / f"{name}_reference.py").read_text()
            block = text[text.index("REFERENCE_CONFIGS = ["):]
            old = re.findall(r'"([^"]*)"', block[:block.index("]")])
            self.assertEqual({canonical(c) for c in old}, {canonical(c) for c in cfg.IONS[ion][1]}, ion)

    def test_ion_name_normalisation(self):
        for text in ("fe_ii", "Fe II", "FeII", "FE_II"):
            self.assertEqual(cfg.normalize_ion(text), "Fe II")
        with self.assertRaises(ValueError):
            cfg.normalize_ion("fe_x")


if __name__ == "__main__":
    unittest.main()
