#!/usr/bin/env python3
"""Why does variant A (Coulomb direct potential + screened nucleus and Slater integrals) fail at mu >= 0.05?

For Fe II, Fe III, Fe IV (gref orbitals, reference configurations of Deprince et al. 2020) and a ladder of mu, each run in
a fresh process on BOTH builds, this records whether the SCF completes, FAC's own error message, the ground-state energy
and the shape of the total electron potential energy r*V_tot(r) = r*(Vc + U) + ZPS read from the converged potential
(negative = attractive).  A build is selected with FE_CALC_REPO (repository root of the worktree); the same single SCF
step (free ion, then PlasmaScreen at mu with radial:hxs_iter0 = 8) is used for both.

Usage: python variant_a_diagnosis.py [--results-dir DIR] [--repo-a PATH] [--repo-b PATH]
Option --csi VALUE sets radial:csi (e.g. 0) in both builds, to separate the effect of the direct potential from FAC's asymptotic-charge
adjustment (csi/ihx) of the mean field.
Worker mode (internal): python variant_a_diagnosis.py --worker ION MU OUTDIR [CSI]
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
sys.path.insert(0, str(HERE.parent / "deprince_2020_fe_ii_viii"))
import configurations as cfg       # noqa: E402
import fac_run_support as sup      # noqa: E402

MU_LADDER = [0.001, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05]
DENSITY_A3, BOHR_A, HARTREE_EV_FAC = 1e-6, 0.52917721067, 27.211386018


def worker(ion: str, mu: float, outdir: Path, csi: float | None = None) -> int:
    outdir.mkdir(parents=True)
    os.chdir(outdir)
    build = sup.find_build_dir()
    sys.path.insert(0, str(build))
    from pfac import fac
    electrons, configs, _ = cfg.IONS[ion]
    fac.SetAtom("Fe")
    if csi is not None:
        fac.SetOption("radial:csi", csi)
    fac.Config(configs[0], group="g0")
    fac.Config(*configs[1:], group="ref")
    fac.OptimizeRadial(["g0"])
    fac.SetOption("orbital:debye_mode", 1)
    fac.SetOption("radial:ee_screen", 1)
    fac.SetOption("radial:hxs_iter0", 8)
    t = 4.0 * math.pi * DENSITY_A3 * BOHR_A**3 / mu**2 * HARTREE_EV_FAC
    fac.PlasmaScreen(cfg.Z_IRON - electrons, DENSITY_A3, t, 1, 0.0, 1)
    fac.OptimizeRadial(["g0"])
    fac.GetPotential("final.pot")
    fac.Structure("levels.b", ["g0", "ref"])
    fac.MemENTable("levels.b")
    fac.PrintTable("levels.b", "levels.txt", 1)
    table = sup.parse_level_table(Path("levels.txt"))
    rows = [l.split() for l in open("final.pot") if l.strip() and not l.startswith("#")]
    r = [float(x[1]) for x in rows]
    rv = [float(x[4]) + float(x[5]) + float(x[17]) for x in rows]
    k = min(range(len(rv)), key=lambda i: rv[i])
    cross = next((r[i] for i in range(k, len(rv)) if rv[i] > 0), None)
    Path("result.json").write_text(json.dumps({"E0_eV": table["e0_ev"], "rVtot_min": rv[k], "r_at_min": r[k],
                                               "r_where_rVtot_turns_repulsive": cross, "rVtot_last": rv[-1]}))
    return 0


def run(repo: Path, ion: str, mu: float, outdir: Path, csi: float | None = None) -> dict:
    env = dict(os.environ, FE_CALC_REPO=str(repo))
    p = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", ion, repr(mu), str(outdir)] + ([repr(csi)] if csi is not None else []),
                       capture_output=True, text=True, env=env, timeout=600)
    msgs = [l.strip() for l in (p.stdout + p.stderr).splitlines() if "RadialBound" in l or "RadialSolver" in l or "Error" in l or "Debye bracket" in l]
    rec = {"ion": ion, "mu": mu, "exit_code": p.returncode, "message": " | ".join(msgs)[:300],
           "status": "ok" if p.returncode == 0 else "failed"}
    res = outdir / "result.json"
    if p.returncode == 0 and res.exists():
        rec.update(json.loads(res.read_text()))
    else:
        (outdir / "run.FAILED.log").write_text(p.stdout + p.stderr)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worker", nargs="+", metavar="ARG", help=argparse.SUPPRESS)
    ap.add_argument("--csi", type=float, default=None, help="set radial:csi (default: FAC default 0.55)")
    ap.add_argument("--results-dir", type=Path, default=HERE / "results" / "variant_A_diagnosis")
    ap.add_argument("--repo-a", type=Path, default=sup.REPO_ROOT.parent / "fac_densitytosz_variant_a")
    ap.add_argument("--repo-b", type=Path, default=sup.REPO_ROOT)
    args = ap.parse_args()
    if args.worker:
        return worker(args.worker[0], float(args.worker[1]), Path(args.worker[2]), float(args.worker[3]) if len(args.worker) > 3 else None)
    out = args.results_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for variant, repo in (("A", args.repo_a), ("B", args.repo_b)):
        for ion in ("Fe II", "Fe III", "Fe IV"):
            for mu in MU_LADDER:
                d = out / variant / cfg.slug(ion) / f"mu{mu:g}".replace(".", "p")
                rec = run(repo.resolve(), ion, mu, d, args.csi)
                rec["variant"] = variant
                rows.append(rec)
                print(variant, ion, mu, rec["status"], rec.get("E0_eV", ""), rec["message"][:90], flush=True)
    fields = ["variant", "ion", "mu", "status", "exit_code", "E0_eV", "rVtot_min", "r_at_min",
              "r_where_rVtot_turns_repulsive", "rVtot_last", "message"]
    with open(out / "variant_a_diagnosis.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
