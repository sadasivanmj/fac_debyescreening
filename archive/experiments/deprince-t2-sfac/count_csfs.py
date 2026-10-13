#!/usr/bin/env python3
"""Count jj-coupled CSFs of active spaces, to compare with Deprince et al.

The paper's active space: all single and double excitations from the Table 1
reference configurations into the orbital set {1s, 2s, 2p, 3s, 3p, 3d, 4s}
("configurations including n = 3 and 4s orbitals").  For every
nonrelativistic configuration all relativistic subshell distributions are
enumerated, the M-generating polynomial of each subshell is multiplied, and
N(J) = c(M=J) - c(M=J+1) gives the number of CSFs of each J.  This is the
number GRASP generates for given (J, parity); FAC builds the same jj-coupled
basis, so its level count over the same blocks is identical.

Writes results/csf_counts.csv and prints a summary.
"""

import csv
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import deprince_t2_sfac as drv

# relativistic subshells: (nonrel orbital, 2j)
RSUB = {"s": [1], "p": [1, 3], "d": [3, 5]}


@lru_cache(None)
def subshell_m(j2: int, q: int) -> Counter:
    """Counter of 2M over all q-electron determinants in a j subshell."""
    ms = range(-j2, j2 + 1, 2)
    return Counter(sum(c) for c in combinations(ms, q))


def mult(a: Counter, b: Counter) -> Counter:
    out = Counter()
    for x, cx in a.items():
        for y, cy in b.items():
            out[x + y] += cx * cy
    return out


@lru_cache(None)
def config_j_counts(occ_key: tuple) -> Counter:
    """{2J: number of CSFs} for a nonrelativistic configuration."""
    total = Counter()
    # distribute q electrons of each nl over its j-subshells
    dists = [[]]
    for orb, q in zip(drv.ORBITALS, occ_key):
        if q == 0:
            continue
        subs = RSUB[orb[1]]
        opts = []
        if len(subs) == 1:
            opts = [((subs[0], q),)]
        else:
            lo, hi = subs
            for qlo in range(0, min(q, lo + 1) + 1):
                qhi = q - qlo
                if qhi <= hi + 1:
                    opts.append(((lo, qlo), (hi, qhi)))
        dists = [d + list(o) for d in dists for o in opts]
    for d in dists:
        m = Counter({0: 1})
        for j2, q in d:
            m = mult(m, subshell_m(j2, q))
        for m2, c in m.items():
            if m2 >= 0:
                total[m2] += c - m.get(m2 + 2, 0)
    return Counter({k: v for k, v in total.items() if v})


def space(ion: str, active: list[str]) -> list[dict]:
    key = lambda o: tuple(o[x] for x in drv.ORBITALS)
    seen = {}
    for spec in drv.REFS[ion]:
        for occ in drv.sd_excitations(drv.parse_occ(spec), active):
            seen[key(occ)] = occ
    return list(seen.values())


def summarize(configs, p0=None, j2=None):
    n = 0
    for occ in configs:
        if p0 is not None and drv.parity(occ) != p0:
            continue
        c = config_j_counts(tuple(occ[x] for x in drv.ORBITALS))
        n += sum(c.values()) if j2 is None else c.get(j2, 0)
    return n


def main():
    rows = []
    for ion in drv.IONS:
        ref0 = drv.parse_occ(drv.REFS[ion][0])
        p0, j0 = drv.parity(ref0), drv.GROUND_2J[ion]
        full = space(ion, drv.ORBITALS)
        sdv = drv.sd_excitations(ref0, drv.VALENCE)
        row = {"ion": ion, "published_NCSF": drv.PUBLISHED_CSF.get(ion, ""),
               "paperAS_nconf": len(full),
               "paperAS_CSF_allJ_bothP": summarize(full),
               "paperAS_CSF_allJ_groundP": summarize(full, p0),
               "paperAS_CSF_groundJP": summarize(full, p0, j0),
               "sdv_nconf": len([o for o in sdv if drv.parity(o) == p0]),
               "sdv_CSF_groundJP": summarize(sdv, p0, j0),
               "ref_CSF_allJ": summarize([drv.parse_occ(s) for s in drv.REFS[ion]])}
        rows.append(row)
        print(row, flush=True)
    out = drv.HERE / "results" / "csf_counts.csv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
