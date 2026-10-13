"""Paper B: configuration data (no FAC needed)."""

import sys
import unittest
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PAPER_DIR))
sys.path.insert(0, str(PAPER_DIR.parent / "shared"))
import configurations as cfg                         # noqa: E402
import validate_inputs                               # noqa: E402
from config_validation import canonical, electron_count   # noqa: E402


class Configurations(unittest.TestCase):
    def test_input_validation_passes(self):
        self.assertEqual(validate_inputs.validate(), [])

    def test_electron_counts_and_charge_states(self):
        for ion, n in zip(cfg.ION_ORDER, range(10, 0, -1)):
            electrons, refs = cfg.IONS[ion]
            self.assertEqual(electrons, n)
            self.assertTrue(all(electron_count(c) == n for c in refs), ion)
            self.assertEqual(cfg.charge_state(ion), 26 - n)

    def test_published_reference_lists(self):
        # Sec. 2.1 of Deprince et al. (2019): number of listed reference configurations
        self.assertEqual([len(cfg.IONS[i][1]) for i in cfg.ION_ORDER[:-1]], [5, 2, 5, 6, 6, 6, 6, 8, 3])
        self.assertEqual(cfg.IONS["Fe XVII"][1][0], "1s2 2s2 2p6")
        self.assertIn("1s1 2s2 2p6 3p1", cfg.IONS["Fe XVII"][1])         # 1s 2s2 2p6 3p
        self.assertIn("1s1 2s1 2p6", cfg.IONS["Fe XIX"][1])
        self.assertEqual(cfg.IONS["Fe XXV"][1], ["1s2", "1s1 2s1", "1s1 2p1"])
        self.assertEqual(cfg.IONS["Fe XXVI"][1], ["1s1"])               # daughter ion, 1s only

    def test_fe_xxiv_is_all_three_electron_occupations_of_1s_2s_2p(self):
        configs = cfg.IONS["Fe XXIV"][1]
        self.assertEqual(configs[0], "1s2 2s1")
        self.assertEqual(len({canonical(c) for c in configs}), 8)
        self.assertIn("2p3", configs)
        self.assertIn("1s1 2s1 2p1", configs)

    def test_models_are_distinct_and_reduced(self):
        for ion in cfg.ION_ORDER:
            g, r = cfg.configurations(ion, "ground_only"), cfg.configurations(ion, "reference_only")
            self.assertEqual(g, [cfg.IONS[ion][1][0]])
            self.assertEqual(r, cfg.IONS[ion][1])
        self.assertEqual(cfg.configurations("Fe XXVI", "ground_only"), cfg.configurations("Fe XXVI", "reference_only"))

    def test_screening_cases_keep_nuclear_and_ee_separate(self):
        self.assertEqual(cfg.CASES, {"unscreened": (0.0, 0.0), "A": (0.1, 0.0), "B": (0.1, 0.1),
                                     "C": (0.25, 0.0), "D": (0.25, 0.25)})

    def test_ion_name_normalisation(self):
        for text in ("fe_xvii", "Fe XVII", "FeXVII"):
            self.assertEqual(cfg.normalize_ion(text), "Fe XVII")


if __name__ == "__main__":
    unittest.main()
