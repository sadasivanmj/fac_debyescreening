#!/usr/bin/env python3
"""Exact first-order Debye-Hueckel check of the screened total energy (both builds).

To first order in mu, exp(-mu r)/r = 1/r - mu, so for ANY normalised N-electron wavefunction the total energy shifts by
    dE = mu * (N Z - N(N-1)/2)    (Hartree)                      [e-e and e-nucleus both screened]
This must hold for the energies FAC reports, whatever the orbitals.  The check compares E(mu) - E(0) of both builds
(A: Coulomb direct potential, B: screened direct potential) with it at small mu, for H-like to open-shell ions.

Usage: python first_order_check.py [--results-dir DIR]      (worker mode is internal)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
import fac_run_support as sup      # noqa: E402

HARTREE_EV = 27.211386
Z = 26
IONS = {"Fe XXVI (1s)": ["1s1"], "Fe XXV (1s2)": ["1s2"], "Fe XVII (2p6, N=10)": ["1s2 2s2 2p6"],
        "Fe IX (3p6, N=18)": ["1s2 2s2 2p6 3s2 3p6"], "Fe III (3d6, N=24)": ["1s2 2s2 2p6 3s2 3p6 3d6"],
        "Fe II (3d6 4s, N=25)": ["1s2 2s2 2p6 3s2 3p6 3d6 4s1"]}
MUS = [0.0, 0.001, 0.003, 0.01]


def worker(config: str, mu: float, outdir: Path) -> int:
    outdir.mkdir(parents=True)
    os.chdir(outdir)
    sys.path.insert(0, str(sup.find_build_dir()))
    from pfac import fac
    n = sum(int(t[2:]) for t in config.split())
    fac.SetAtom("Fe")
    fac.Config(config, group="g")
    fac.OptimizeRadial(["g"])
    if mu > 0:
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        fac.SetOption("radial:hxs_iter0", 8)
        t = 4 * math.pi * 1e-6 * 0.52917721067**3 / mu**2 * 27.211386018
        fac.PlasmaScreen(Z - n, 1e-6, t, 1, 0.0, 1)
        fac.OptimizeRadial(["g"])
    fac.Structure("l.b", ["g"])
    fac.MemENTable("l.b")
    fac.PrintTable("l.b", "l.txt", 1)
    table = sup.parse_level_table(Path("l.txt"))
    Path("result.json").write_text(json.dumps({"E0_eV": table["e0_ev"], "N": n}))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worker", nargs=3)
    ap.add_argument("--results-dir", type=Path, default=HERE / "results" / "first_order_check")
    ap.add_argument("--repo-a", type=Path, default=sup.REPO_ROOT.parent / "fac_densitytosz_variant_a")
    ap.add_argument("--repo-b", type=Path, default=sup.REPO_ROOT)
    args = ap.parse_args()
    if args.worker:
        return worker(args.worker[0], float(args.worker[1]), Path(args.worker[2]))
    out = args.results_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, cfgs in IONS.items():
        n = sum(int(t[2:]) for t in cfgs[0].split())
        e = {}
        for variant, repo in (("A", args.repo_a), ("B", args.repo_b)):
            for mu in MUS:
                if mu == 0.0 and variant == "A":
                    e[("A", 0.0)] = None          # mu = 0 is plain FAC and identical in both builds
                    continue
                d = out / "runs" / variant / name.split(" (")[0].replace(" ", "_") / f"mu{mu:g}"
                p = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", cfgs[0], repr(mu), str(d)],
                                   capture_output=True, text=True, env=dict(os.environ, FE_CALC_REPO=str(repo.resolve())), timeout=900)
                e[(variant, mu)] = json.loads((d / "result.json").read_text())["E0_eV"] if p.returncode == 0 else None
        e0 = e[("B", 0.0)]
        for variant in "AB":
            for mu in MUS[1:]:
                x = e[(variant, mu)]
                first = mu * (n * Z - n * (n - 1) / 2) * HARTREE_EV
                rows.append({"ion": name, "N": n, "variant": variant, "mu": mu, "dE_FAC_eV": None if x is None else x - e0,
                             "dE_first_order_exact_eV": first, "deviation_eV": None if x is None else x - e0 - first,
                             "deviation_percent": None if x is None else 100 * (x - e0 - first) / first})
        print(name, {(r["variant"], r["mu"]): round(r["deviation_percent"], 3) for r in rows if r["ion"] == name and r["deviation_percent"] is not None}, flush=True)
    with open(out / "first_order_check.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
