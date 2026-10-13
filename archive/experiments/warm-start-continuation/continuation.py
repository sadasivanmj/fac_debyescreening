#!/usr/bin/env python3
"""Warm-start (continuation) experiment for the screened Fe II / Fe III SCF.

usage: python3 continuation.py MODE ION MODEL OUTDIR MU1 [MU2 ...]
  MODE  ladder : free-ion OptimizeRadial, then in the SAME process
                 PlasmaScreen(mu_k) + OptimizeRadial for each listed mu_k.
                 The converged orbitals of step k-1 stay in FAC's orbital
                 table and are the starting point of step k (no file I/O).
        repeat : as ladder, then one more OptimizeRadial at the LAST mu
                 with nothing changed.
        cold   : debye_mode/ee_screen/PlasmaScreen(MU1) set BEFORE the first
                 OptimizeRadial (no free-ion step).
Physical inputs are those of verify_ip_deprince.py (feature branch ee326b02):
Table 1 reference configurations; nps = 1e-6 A^-3, ups = 0, m = 1, vxf = 1,
zps = ion charge, T chosen so that dps = 1/mu; orbital:debye_mode 1,
radial:ee_screen 1; model ref = orbitals on the average of all references,
gref = on the ground reference.  Nothing is tuned.

After every OptimizeRadial the potential is dumped (GetPotential) and its
header is checked: mps, dps (effective 1/mu), HXS/AHX (exchange scale) and
the orbital energies of the 3d and 4s orbitals.  With SetOptimizePrint(1)
FAC prints one 'optimize loop' line per SCF iteration; at iteration 0,
tol = 1, atol = 10 means no previous orbital existed (cold), anything else
means the previous orbitals were reused (warm).
"""
import json
import math
import re
import sys
from pathlib import Path

RBOHR_ANG, FAC_HARTREE_EV, NPS_A3, UPS = 0.52917721067, 27.211386018, 1e-6, 0.0
REFS = {
    "Fe II": ["1s2 2s2 2p6 3s2 3p6 3d6 4s1", "1s2 2s2 2p6 3s2 3p5 3d7 4s1",
              "1s2 2s2 2p6 3s2 3p5 3d6 4s2", "1s2 2s2 2p5 3s2 3p6 3d7 4s1",
              "1s2 2s2 2p5 3s2 3p6 3d6 4s2", "1s1 2s2 2p6 3s2 3p6 3d7 4s1",
              "1s1 2s2 2p6 3s2 3p6 3d6 4s2"],
    "Fe III": ["1s2 2s2 2p6 3s2 3p6 3d6", "1s2 2s2 2p6 3s2 3p5 3d7",
               "1s2 2s2 2p6 3s2 3p5 3d6 4s1", "1s2 2s2 2p6 3s2 3p5 3d5 4s2",
               "1s2 2s2 2p5 3s2 3p6 3d7", "1s2 2s2 2p5 3s2 3p6 3d6 4s1",
               "1s2 2s2 2p5 3s2 3p6 3d5 4s2", "1s1 2s2 2p6 3s2 3p6 3d7",
               "1s1 2s2 2p6 3s2 3p6 3d6 4s1", "1s1 2s2 2p6 3s2 3p6 3d5 4s2"],
}
NELE = {"Fe II": 25, "Fe III": 24}
WATCH = {"3,2": "3d3/2", "3,-3": "3d5/2", "4,-1": "4s1/2", "3,-1": "3s1/2"}


def tps_for_mu(mu):
    return (1.0 / mu) ** 2 * 4.0 * math.pi * NPS_A3 * RBOHR_ANG ** 3 * (1.0 + UPS) * FAC_HARTREE_EV


def pot_info(p):
    h, orb = {}, {}
    for line in open(p):
        if not line.startswith("#"):
            break
        m = re.match(r"#\s+(\w+) =\s+(\S+)$", line)
        if m:
            h[m.group(1)] = m.group(2)
        f = line[1:].split()
        if len(f) == 9 and f[0].isdigit():
            orb[f"{f[1]},{f[2]}"] = float(f[8])
    return {"mps": int(h["mps"]), "dps": float(h["dps"]), "ips": int(h["ips"]),
            "HXS": float(h["HXS"]), "AHX": float(h["AHX"]),
            "e_Ha": {WATCH[k]: v for k, v in orb.items() if k in WATCH}}


def main():
    mode, ion, model, out = sys.argv[1:5]
    mus = [float(x) for x in sys.argv[5:]]
    from pfac import fac
    print("IMPORTED", fac.__file__, flush=True)
    out = Path(out)
    out.mkdir(parents=True)  # refuses to reuse an existing directory
    import os
    os.chdir(out)
    fac.SetAtom("Fe")
    fac.SetOptimizePrint(1)
    refs = REFS[ion]
    if model == "ref":
        fac.Config(*refs, group="ref")
        opt = ["ref"]
    else:
        fac.Config(refs[0], group="g0")
        fac.Config(*refs[1:], group="ref")
        opt = ["g0"]
    zps = 26.0 - NELE[ion]

    def screen(mu):
        t = tps_for_mu(mu)
        fac.PlasmaScreen(zps, NPS_A3, t, 1, UPS, 1)
        print(f"STEP mu={mu} PlasmaScreen({zps}, {NPS_A3}, {t!r}, 1, {UPS}, 1)", flush=True)

    def optimize(tag, mu):
        print(f"STEP {tag} OptimizeRadial start", flush=True)
        fac.OptimizeRadial(opt)
        fac.GetPotential(f"{tag}.pot")
        info = pot_info(f"{tag}.pot")
        ok = (mu == 0 and info["mps"] != 1) or (mu > 0 and info["mps"] == 1
                                                 and abs(info["dps"] * mu - 1) < 1e-6)
        print("RESULT", json.dumps({"tag": tag, "mu_requested": mu, "screening_ok": ok, **info}),
              flush=True)

    if mode == "cold":
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        screen(mus[0])
        optimize(f"cold_mu{mus[0]:g}", mus[0])
    else:
        optimize("free", 0.0)
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        for mu in mus:
            screen(mu)
            optimize(f"step_mu{mu:g}", mu)
        if mode == "repeat":
            optimize(f"repeat_mu{mus[-1]:g}", mus[-1])
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
