"""Configurations and run conditions for Deprince et al. (2019), Fe XVII - Fe XXV (+ Fe XXVI).

  Paper B: J. Deprince et al., A&A 626, A83 (2019).  The list of reference (spectroscopic)
  configurations is in Sec. 2.1 of the paper (its Table 1 is the table of mu values).
  Notation of the paper: Fe XVII 1s2 2s2 2p6, 1s2 2s2 2p5 3s, ...; Fe XXIV is "(1s+2s+2p)^3",
  i.e. every 3-electron occupation of 1s, 2s and 2p.

  Fe XXVI (hydrogen-like) is not in the paper: it is the daughter ion needed for IP(Fe XXV).
  Its only configuration is 1s1.

Two reduced models, always kept apart:
  ground_only     : the ground configuration alone (the existing validated benchmark)
  reference_only  : the listed reference configurations, with orbitals optimised on the ground
                    configuration only and the references added in Structure.
NEITHER generates the single/double excitations of the paper's restricted active space, and
neither is the published MCDF calculation (no Breit, QED or extended-average-level optimisation).
"""

import itertools

PAPER = "deprince_2019_fe_xvii_xxv"
Z_IRON = 26

CONFIGURATION_MODELS = ("ground_only", "reference_only")
ORBITAL_OPTIMIZATION = "ground"     # both models optimise orbitals on the ground configuration

# Screening conditions (mu_ne, mu_ee) in atomic units; the two are controlled separately.
#   unscreened : plain FAC, PlasmaScreen is never called
#   A (0.1, 0)   nuclear screening only; electron-electron repulsion stays plain Coulomb
#   B (0.1, 0.1) nuclear and electron-electron screening
#   C (0.25, 0)  D (0.25, 0.25) likewise at mu = 0.25
CASES = {"unscreened": (0.0, 0.0), "A": (0.1, 0.0), "B": (0.1, 0.1), "C": (0.25, 0.0), "D": (0.25, 0.25)}


def _fe_xxiv_configurations():
    """(1s+2s+2p)^3: all ways to put 3 electrons into 1s (2), 2s (2), 2p (6); ground first."""
    out = []
    for n1s, n2s, n2p in itertools.product(range(3), range(3), range(7)):
        if n1s + n2s + n2p == 3:
            out.append(" ".join(f"{name}{n}" for name, n in (("1s", n1s), ("2s", n2s), ("2p", n2p)) if n))
    ground = "1s2 2s1"
    return [ground] + sorted(c for c in out if c != ground)


# ion -> (electrons, reference configurations, ground configuration first)
IONS = {
    "Fe XVII":  (10, ["1s2 2s2 2p6", "1s2 2s2 2p5 3s1", "1s2 2s2 2p5 3p1", "1s1 2s2 2p6 3s1", "1s1 2s2 2p6 3p1"]),
    "Fe XVIII": (9,  ["1s2 2s2 2p5", "1s1 2s2 2p6"]),
    "Fe XIX":   (8,  ["1s2 2s2 2p4", "1s2 2s1 2p5", "1s2 2p6", "1s1 2s2 2p5", "1s1 2s1 2p6"]),
    "Fe XX":    (7,  ["1s2 2s2 2p3", "1s2 2s1 2p4", "1s2 2p5", "1s1 2s2 2p4", "1s1 2s1 2p5", "1s1 2p6"]),
    "Fe XXI":   (6,  ["1s2 2s2 2p2", "1s2 2s1 2p3", "1s2 2p4", "1s1 2s2 2p3", "1s1 2s1 2p4", "1s1 2p5"]),
    "Fe XXII":  (5,  ["1s2 2s2 2p1", "1s2 2s1 2p2", "1s2 2p3", "1s1 2s2 2p2", "1s1 2s1 2p3", "1s1 2p4"]),
    "Fe XXIII": (4,  ["1s2 2s2", "1s2 2s1 2p1", "1s2 2p2", "1s1 2s2 2p1", "1s1 2s1 2p2", "1s1 2p3"]),
    "Fe XXIV":  (3,  _fe_xxiv_configurations()),
    "Fe XXV":   (2,  ["1s2", "1s1 2s1", "1s1 2p1"]),
    "Fe XXVI":  (1,  ["1s1"]),
}
ION_ORDER = list(IONS)


def configurations(ion: str, model: str) -> list[str]:
    refs = IONS[ion][1]
    return refs[:1] if model == "ground_only" else refs


def slug(ion: str) -> str:
    return ion.lower().replace(" ", "_")


def normalize_ion(text: str) -> str:
    """'fe_xvii', 'Fe XVII', 'FeXVII' -> 'Fe XVII'."""
    t = text.replace("_", " ").strip().upper()
    if t.startswith("FE") and " " not in t:
        t = "FE " + t[2:]
    name = "Fe " + t.split()[-1]
    if name not in IONS:
        raise ValueError(f"unknown ion {text!r}; choose from {ION_ORDER}")
    return name


def charge_state(ion: str) -> int:
    return Z_IRON - IONS[ion][0]
