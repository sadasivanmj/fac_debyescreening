#!/usr/bin/env python3
"""Focused native-SFAC diagnostics for the screened radial-solver failure.

Original symptom (Fe II, mu = 0.1, Table 1 reference configurations):

    Max iteration before finding correct nodes in RadialBound a:
    0 0 279 3 2 -8.12838e-154
    Error occured in RadialSolver: -2 3 2 0.0000000000E+00

Each case below is a separate sfac process in diagnostics/<case>/ with its
input, log, exit status and (on failure) FAC's error.pot preserved.  Results
go to diagnostics/diagnostics.json and diagnostics/diagnostics.md.

Binaries (see REPORT.md for provenance):
    installed      ~/facinst/bin/sfac (bit-identical to ee-screening-scf ff9a5c9)
    scf, scf_fix   ee-screening-scf build without / with debye_sfac_fixes.patch
    ee, ee_fix     fix-debye-ee-screening build without / with the patch
"""

from __future__ import annotations

import json
import math
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import deprince_t2_sfac as drv

HERE = Path(__file__).resolve().parent
OUT = HERE / "diagnostics"
WT = Path("/data/project/fac_worktrees")
BIN = {
    "installed": Path.home() / "facinst/bin/sfac",
    "scf": WT / "scf/inst_unpatched/bin/sfac",
    "scf_fix": WT / "scf/inst/bin/sfac",
    "ee": WT / "ee/inst_unpatched/bin/sfac",
    "ee_fix": WT / "ee/inst/bin/sfac",
}
FEII_REFS = [drv.occ_to_fac(drv.parse_occ(s)) for s in drv.REFS["FeII"]]


def tps(mu, nps):
    n_au = nps * drv.RBOHR_ANG**3
    return (1 / mu) ** 2 * 4 * math.pi * n_au * drv.HARTREE_EV


def deck(configs, mu=0.0, nps=1e-6, ee=1, en=1, ladder=False, zps=1.0,
         config_energy=False, hxs_iter0=8, structure=True):
    """SFAC input.  ladder=False: one cold OptimizeRadial at mu (the original
    driver's scheme); ladder=True: free-ion start, then mu in 0.05 steps."""
    L = ["SetAtom('Fe')",
         "Config('g1', " + ", ".join(repr(c) for c in configs) + ")"]
    screen = lambda m: (f"PlasmaScreen({zps:.1f}, {nps:.10g}, "
                        f"{tps(m, nps)!r}, 1, 0.0, 1)")
    opts = []
    if mu > 0:
        opts = [f"SetOption('orbital:debye_mode', {en})",
                f"SetOption('radial:ee_screen', {ee})"]
    if not ladder:
        L += opts + ([screen(mu)] if mu > 0 else [])
        L += ["ConfigEnergy(0)"] if config_energy else []
        L += ["OptimizeRadial(['g1'])"]
        L += ["ConfigEnergy(1)"] if config_energy else []
    else:
        L += ["OptimizeRadial(['g1'])"]
        if mu > 0:
            L += opts + [f"SetOption('radial:hxs_iter0', {hxs_iter0})"]
            for m in drv.ladder(mu) if mu >= 0.05 else [mu]:
                L += [screen(m), "OptimizeRadial(['g1'])"]
    L += ["GetPotential('final.pot')"]
    if structure:
        L += ["Structure('d.lev.b', ['g1'])", "MemENTable('d.lev.b')",
              "PrintTable('d.lev.b', 'd.lev', 1)"]
    return "\n".join(L) + "\n"


def pot_tail(path: Path) -> dict:
    """Asymptotic effective charge r*V = Vc*r + U*r + ZPS seen by an electron."""
    rows = [l.split() for l in path.read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    rows = [[float(x) for x in r] for r in rows if len(r) > 17]
    hdr = drv.parse_pot(path)
    rv = [(r[1], r[4] + r[5] + r[17]) for r in rows]
    big = [x for x in rv if x[0] > 50]
    pos = [r for r, v in rv if v > 0 and r > 0.5]
    return {"rV_at_r>50_max": max(v for _, v in big) if big else None,
            "V_positive_from_r": min(pos) if pos else None,
            "ips": hdr.get("ips"), "rps": hdr.get("rps"), "dps": hdr.get("dps"),
            "HXS": hdr.get("HXS"), "orb_3d3/2_Ha": hdr["orb"].get("3,2"),
            "orb_4s_Ha": hdr["orb"].get("4,-1")}


def run(case):
    name, binary, sf, note = case
    wd = OUT / name
    wd.mkdir(parents=True, exist_ok=True)
    for f in wd.iterdir():
        f.unlink()
    (wd / "in.sf").write_text(sf)
    exe = BIN[binary]
    with open(wd / "sfac.log", "w") as log:
        rc = subprocess.run([str(exe), "in.sf"], cwd=wd, stdout=log,
                            stderr=subprocess.STDOUT, timeout=7200).returncode
    log = (wd / "sfac.log").read_text().strip().splitlines()
    rec = {"case": name, "binary": binary, "sfac": str(exe),
           "sfac_sha256": drv.sha256(exe), "returncode": rc, "note": note,
           "log_tail": [l for l in log if "optimize loop" not in l][-2:]}
    if (wd / "d.lev").exists():
        rec.update(drv.parse_lev(wd / "d.lev"))
    for p in ("final.pot", "error.pot"):
        if (wd / p).exists():
            rec[p] = pot_tail(wd / p)
    rec["status"] = "ok" if rc == 0 and rec.get("E0_eV") is not None else "FAILED"
    return rec


def cases():
    g = [FEII_REFS[0]]
    C = []
    add = lambda *a: C.append(a)
    # 1. original failure and its reproduction with a single configuration
    add("A1_orig_refs_mu0.1", "installed", deck(FEII_REFS, 0.1, nps=0.01),
        "original driver input (nps=0.01 A^-3, cold start)")
    add("A2_ground_mu0.1", "installed", deck(g, 0.1, nps=0.01),
        "ground configuration only")
    # 2. incremental reference configurations (original binary and input style)
    for k in range(2, len(FEII_REFS) + 1):
        add(f"A3_refs1-{k}_mu0.1", "installed", deck(FEII_REFS[:k], 0.1, nps=0.01),
            f"first {k} Table 1 references")
    # 3. screening switches
    add("B1_en_only_mu0.1", "installed", deck(g, 0.1, nps=0.01, ee=0),
        "e-n screening only (radial:ee_screen 0)")
    add("B2_en_ee_mu0.1", "installed", deck(g, 0.1, nps=0.01, ee=1),
        "e-n and e-e screening")
    # 4. cause 1: ion-sphere truncation of the screened Hartree potential
    add("C1_nps1e-2", "installed", deck(g, 0.1, nps=1e-2),
        "rps=5.44 a0: screened u(r) frozen beyond rps")
    add("C2_nps1e-6", "installed", deck(g, 0.1, nps=1e-6),
        "rps=117 a0: potential decays correctly, still fails")
    # 5. cause 2: RadialBound bracketing with finitely many bound states
    add("D1_fix_cold_mu0.1", "scf_fix", deck(g, 0.1),
        "patched solver, cold start: genuine unbound transient at SCF iter 2")
    add("D2_fix_ladder_mu0.1", "scf_fix", deck(g, 0.1, ladder=True),
        "patched solver, free start + warm ladder")
    add("D3_fix_ladder_mu0.25", "scf_fix", deck(g, 0.25, ladder=True),
        "patched solver, ladder to 0.25")
    add("D4_fix_ladder_refs_mu0.25", "scf_fix", deck(FEII_REFS, 0.25, ladder=True),
        "patched solver, all Table 1 refs, ladder to 0.25")
    add("D5_unpatched_ladder_mu0.1", "scf", deck(g, 0.1, ladder=True),
        "unpatched solver, ladder: hxs_iter0 unknown option, iteration-0 "
        "exchange-free potential")
    add("D6_fix_ladder_hxs0_mu0.1", "scf_fix", deck(g, 0.1, ladder=True, hxs_iter0=0),
        "patched, ladder without ramp offset")
    # 6. invariance checks
    add("E1_free_unpatched", "scf", deck(FEII_REFS), "free ion, unpatched")
    add("E2_free_patched", "scf_fix", deck(FEII_REFS), "free ion, patched")
    add("E3_free_installed_CE", "installed", deck(FEII_REFS, config_energy=True),
        "original baseline with ConfigEnergy (expect -34561.2841 eV)")
    add("E4_zps3_mu0.1", "scf_fix", deck(g, 0.1, ladder=True, zps=3.0),
        "zps=3 instead of 1: must not change Debye-mode results")
    add("E5_nps1e-4_mu0.1", "scf_fix", deck(g, 0.1, ladder=True, nps=1e-4),
        "nps=1e-4 (T scaled): lambda_D unchanged")
    add("E6_mu1e-4_ladder", "scf_fix", deck(g, 1e-4, ladder=True),
        "large but finite lambda_D = 1e4 a0")
    # 7. fix-debye-ee-screening branch (unscreened SCF direct potential)
    add("F1_ee_branch_ladder_mu0.1", "ee_fix", deck(g, 0.1, ladder=True),
        "fix-debye-ee-screening + patch: DensityToSZ unscreened")
    add("F2_ee_branch_free", "ee_fix", deck(FEII_REFS), "fix-debye-ee-screening free")
    return C


def main():
    OUT.mkdir(exist_ok=True)
    only = sys.argv[1:]
    todo = [c for c in cases() if not only or any(c[0].startswith(o) for o in only)]
    with ThreadPoolExecutor(8) as ex:
        recs = list(ex.map(run, todo))
    old = {}
    if (OUT / "diagnostics.json").exists():
        old = {r["case"]: r for r in json.loads((OUT / "diagnostics.json").read_text())}
    old.update({r["case"]: r for r in recs})
    allr = [old[k] for k in sorted(old)]
    (OUT / "diagnostics.json").write_text(json.dumps(allr, indent=1))
    md = ["# Radial-solver diagnostics (native SFAC)", "",
          "| case | binary | rc | status | E0 (eV) | r·V(r>50) max | "
          "ips | 3d3/2 (Ha) | last log line | note |",
          "|---|---|---:|---|---:|---:|---:|---:|---|---|"]
    for r in allr:
        p = r.get("final.pot") or r.get("error.pot") or {}
        rv = p.get("rV_at_r>50_max")
        md.append(f"| {r['case']} | {r['binary']} | {r['returncode']} | {r['status']} | "
                  f"{drv.fmt(r.get('E0_eV'))} | {'' if rv is None else f'{rv:.3f}'} | "
                  f"{p.get('ips', '')} | {p.get('orb_3d3/2_Ha', '')} | "
                  f"{(r['log_tail'] or [''])[-1][:70]} | {r['note']} |")
    (OUT / "diagnostics.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
