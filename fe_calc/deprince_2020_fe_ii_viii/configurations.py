"""Reference configurations and run conditions for Deprince et al. (2020), Fe II - Fe VIII.

  Paper A: J. Deprince et al., A&A 643, A57 (2020), Table 1 (Fe II - Fe VIII).
  Fe IX  : needed for IP(Fe VIII); its reference configurations are in Paper III,
           J. Deprince et al., A&A 635, A70 (2020), Table 1:  3p6, [3p]3d, [2p]3d, [1s]3d.

Notation: [nl] is a hole in subshell nl; the electron goes to 3d (or 4s) instead.
All ions share the closed core 1s2 2s2 2p6 3s2 3p6 (18 electrons).  The first
configuration of every list is the ground configuration.

REDUCED MODEL: only these configurations are used.  No single or double excitations are
generated; the paper's active spaces (12 268 - 85 798 CSFs) are NOT reproduced.
"""

PAPER = "deprince_2020_fe_ii_viii"
Z_IRON = 26
CONFIGURATION_MODEL = "reference_configurations"

# Orbital-optimisation models (same CI space for both):
#   ref  : OptimizeRadial on the average of ALL reference configurations
#   gref : OptimizeRadial on the ground configuration only; the other references are
#          added only in Structure.  gref is the primary reduced-model result.
ORBITAL_MODELS = ("gref", "ref")

# Run conditions: electron-nucleus and electron-electron screening equal (Paper A).
CASES = {"mu0p00": 0.0, "mu0p10": 0.1, "mu0p25": 0.25}

CORE = {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6}


def config(hole, n3d, n4s):
    """FAC string: closed core (optionally with one hole) + 3d^n3d 4s^n4s."""
    shells = dict(CORE)
    if hole:
        shells[hole] -= 1
    parts = [f"{name}{occ}" for name, occ in shells.items()]
    if n3d:
        parts.append(f"3d{n3d}")
    if n4s:
        parts.append(f"4s{n4s}")
    return " ".join(parts)


def with_holes(valence):
    """`valence` = [(3d, 4s) of the ground configuration, then (3d, 4s) of each hole
    configuration].  Returns the ground configuration followed by the same valence
    occupations with a hole in 3p, 2p and 1s (the [3p], [2p], [1s] rows of Table 1)."""
    ground, extra = valence[0], valence[1:]
    out = [config(None, *ground)]
    for hole in ("3p", "2p", "1s"):
        out += [config(hole, *v) for v in extra]
    return out


# ion -> (number of electrons, configurations, source table)
IONS = {
    "Fe II":   (25, with_holes([(6, 1), (7, 1), (6, 2)]), "Paper IV Table 1"),
    "Fe III":  (24, with_holes([(6, 0), (7, 0), (6, 1), (5, 2)]), "Paper IV Table 1"),
    "Fe IV":   (23, with_holes([(5, 0), (6, 0)]), "Paper IV Table 1"),
    "Fe V":    (22, with_holes([(4, 0), (5, 0)]), "Paper IV Table 1"),
    "Fe VI":   (21, with_holes([(3, 0), (4, 0)]), "Paper IV Table 1"),
    "Fe VII":  (20, with_holes([(2, 0), (3, 0)]), "Paper IV Table 1"),
    "Fe VIII": (19, with_holes([(1, 0), (2, 0)]), "Paper IV Table 1"),
    "Fe IX":   (18, with_holes([(0, 0), (1, 0)]), "Paper III Table 1"),
}
ION_ORDER = list(IONS)


def slug(ion: str) -> str:
    return ion.lower().replace(" ", "_")


def normalize_ion(text: str) -> str:
    """'fe_ii', 'Fe II', 'FeII' -> 'Fe II'."""
    t = text.replace("_", " ").strip().upper()
    if t.startswith("FE") and " " not in t:
        t = "FE " + t[2:]
    name = "Fe " + t.split()[-1]
    if name not in IONS:
        raise ValueError(f"unknown ion {text!r}; choose from {ION_ORDER}")
    return name


def charge_state(ion: str) -> int:
    return Z_IRON - IONS[ion][0]


def case_name(mu: float) -> str:
    return f"mu{mu:.2f}".replace(".", "p")
