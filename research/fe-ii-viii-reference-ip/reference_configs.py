# Reference configurations of Deprince et al. for Fe II - Fe IX, written out for FAC.
#
#   Fe II - Fe VIII : A&A 643, A57 (2020), Table 1   (Paper IV)
#   Fe IX           : A&A 635, A70 (2020), Table 1   (Paper III), needed for IP(Fe VIII)
#
# Notation: [nl] means one electron removed from the nl subshell, and the
# electron is put into the 3d (or 4s) subshell instead.  All ions share the
# closed core 1s2 2s2 2p6 3s2 3p6 (18 electrons).  The first entry of each list
# is the ground configuration.

CORE = {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6}


def config(hole, n3d, n4s):
    """FAC configuration string: closed core (with an optional hole) + 3d^n3d 4s^n4s."""
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
    """Ground valence occupation + the same valence with one more 3d/4s electron
    in each of [3p], [2p], [1s].  `valence` lists (n3d, n4s) of the ground
    configuration followed by the (n3d, n4s) of the hole configurations."""
    ground, extra = valence[0], valence[1:]
    out = [config(None, *ground)]
    for hole in ("3p", "2p", "1s"):
        out += [config(hole, *v) for v in extra]
    return out


# ion -> (electrons, list of configurations)
IONS = {
    "Fe II":   (25, with_holes([(6, 1), (7, 1), (6, 2)])),
    "Fe III":  (24, with_holes([(6, 0), (7, 0), (6, 1), (5, 2)])),
    "Fe IV":   (23, with_holes([(5, 0), (6, 0)])),
    "Fe V":    (22, with_holes([(4, 0), (5, 0)])),
    "Fe VI":   (21, with_holes([(3, 0), (4, 0)])),
    "Fe VII":  (20, with_holes([(2, 0), (3, 0)])),
    "Fe VIII": (19, with_holes([(1, 0), (2, 0)])),
    "Fe IX":   (18, with_holes([(0, 0), (1, 0)])),
}
# Fe II and Fe III have the [3p] configurations first, then [2p], then [1s], as in
# Table 1 (3p holes with 3d^7 4s, 3d^6 4s^2, ...); the order does not matter to FAC.
# The ground configuration listed first is the one used by the "gref" orbital model.

# (ion with one more electron removed): used to pair ions in the IP calculation
ION_ORDER = list(IONS)


def count_electrons(cfg):
    return sum(int(shell[2:]) for shell in cfg.split())


if __name__ == "__main__":
    for ion, (n, cfgs) in IONS.items():
        assert all(count_electrons(c) == n for c in cfgs), ion
        print(f"{ion}: {n} electrons, {len(cfgs)} configurations")
        for c in cfgs:
            print("   ", c)
