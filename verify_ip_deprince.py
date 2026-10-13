#!/usr/bin/env python3
"""Fe II and Fe III ionization potentials vs Deprince et al., A&A 643, A57 (2020), Table 2.

Stage A (computed): reference-configuration space of Table 1, mu = 0, 0.1, 0.25 a0^-1.
Stage B (counted, not computed): the paper's SD active space (single + double
excitations from every reference into {1s, 2s, 2p, 3s, 3p, 3d, 4s}).

usage:
  python3 verify_ip_deprince.py [--out DIR] [--jobs N] [--hxs-iter0 X]   # driver
  python3 verify_ip_deprince.py --case ION MU MODEL DIR     # one FAC case (internal)

The extension imported is whatever `pfac` resolves to; the driver prints it and
refuses to run unless it lies inside this worktree's build/ directory.

Models (both Stage A, same CI space = all reference configurations, all J):
  ref   orbitals optimized on the average of all reference configurations
        (FAC's analogue of the paper's EAL over the references)
  gref  orbitals optimized on the ground reference configuration only
Protocol per case (one FAC process): free-ion OptimizeRadial (fixes the grid);
for mu > 0: orbital:debye_mode 1, radial:ee_screen 1, then
PlasmaScreen(zps, 1e-6 A^-3, T(mu_k), 1, 0.0, 1) and OptimizeRadial for
mu_k = 0.05, 0.10, ..., mu (warm-started ladder); Structure; PrintTable.
--hxs-iter0 X additionally sets radial:hxs_iter0 = X together with the Debye
options (screened steps only; needs the fix/hxs-iter0 build).  Default: unset.
E0 = lowest level (FAC prints level energies in eV).
"""

import csv
import json
import math
import os
import re
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HARTREE_EV = 27.211386            # requested conversion; FAC itself uses 27.211386018
RBOHR_ANG = 0.52917721067         # FAC consts.h _RBOHR
FAC_HARTREE_EV = 27.211386018     # FAC consts.h _HARTREE_EV (used for T in eV)
NPS_A3 = 1e-6                     # electron density [A^-3]; only sets T for a given mu
UPS = 0.0
LADDER_STEP = 0.05
HXS_ITER0 = None                  # radial:hxs_iter0 for screened steps (None = unset)
MU_VALUES = [0.0, 0.1, 0.25]
MODELS = ["ref", "gref"]
NELE = {"Fe II": 25, "Fe III": 24, "Fe IV": 23}

REFS = {
    "Fe II": [
        "1s2 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d7 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s2",
    ],
    "Fe III": [
        "1s2 2s2 2p6 3s2 3p6 3d6",
        "1s2 2s2 2p6 3s2 3p5 3d7",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d5 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d5 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d5 4s2",
    ],
    "Fe IV": [
        "1s2 2s2 2p6 3s2 3p6 3d5",
        "1s2 2s2 2p6 3s2 3p5 3d6",
        "1s2 2s2 2p5 3s2 3p6 3d6",
        "1s1 2s2 2p6 3s2 3p6 3d6",
    ],
}
PUBLISHED_CSF = {"Fe II": 12268, "Fe III": 31501, "Fe IV": 57189}
TABLE2 = {  # Deprince et al. (2020) Table 2, eV
    ("Fe II", 0.0): 18.89, ("Fe II", 0.1): 14.20, ("Fe II", 0.25): 7.78,
    ("Fe III", 0.0): 32.93, ("Fe III", 0.1): 25.60, ("Fe III", 0.25): 16.42,
}
IPS = [("Fe II", "Fe III"), ("Fe III", "Fe IV")]
ORBITALS = ["1s", "2s", "2p", "3s", "3p", "3d", "4s"]
CAP = {"s": 2, "p": 6, "d": 10}
LVAL = {"s": 0, "p": 1, "d": 2}


# ---------------------------------------------------------------- configs --
def parse_cfg(s):
    occ = dict.fromkeys(ORBITALS, 0)
    for tok in s.split():
        m = re.fullmatch(r"(\d[spd])(\d+)", tok)
        if not m or m.group(1) not in occ:
            raise ValueError(f"bad token {tok!r} in {s!r}")
        occ[m.group(1)] = int(m.group(2))
        if occ[m.group(1)] > CAP[tok[1]]:
            raise ValueError(f"over-full subshell {tok!r} in {s!r}")
    return occ


def check_refs():
    for ion, cfgs in REFS.items():
        for c in cfgs:
            n = sum(parse_cfg(c).values())
            if n != NELE[ion]:
                raise SystemExit(f"{ion}: {c!r} has {n} electrons, expected {NELE[ion]}")


def parity(occ):
    return sum(LVAL[o[1]] * q for o, q in occ.items()) % 2


def sd_space(refs):
    key = lambda o: tuple(o[x] for x in ORBITALS)
    space = {}
    for r in refs:  # SD of each reference separately (own `seen`), then the union
        front = [r]
        seen = {key(r): r}
        for _ in range(2):
            nxt = []
            for occ in front:
                for a in ORBITALS:
                    if not occ[a]:
                        continue
                    for b in ORBITALS:
                        if b != a and occ[b] < CAP[b[1]]:
                            new = dict(occ)
                            new[a] -= 1
                            new[b] += 1
                            if key(new) not in seen:
                                seen[key(new)] = new
                                nxt.append(new)
            front = nxt
        space.update(seen)
    return list(space.values())


# jj-coupled CSF count per J: N(J) = c(M=J) - c(M=J+1)
RSUB = {"s": [1], "p": [1, 3], "d": [3, 5]}


@lru_cache(None)
def _sub_m(j2, q):
    return Counter(sum(c) for c in combinations(range(-j2, j2 + 1, 2), q))


def _mult(a, b):
    out = Counter()
    for x, cx in a.items():
        for y, cy in b.items():
            out[x + y] += cx * cy
    return out


@lru_cache(None)
def csf_by_2j(key):
    dists = [[]]
    for orb, q in zip(ORBITALS, key):
        if not q:
            continue
        subs = RSUB[orb[1]]
        if len(subs) == 1:
            opts = [((subs[0], q),)]
        else:
            lo, hi = subs
            opts = [((lo, a), (hi, q - a)) for a in range(0, min(q, lo + 1) + 1) if q - a <= hi + 1]
        dists = [d + list(o) for d in dists for o in opts]
    tot = Counter()
    for d in dists:
        m = Counter({0: 1})
        for j2, q in d:
            m = _mult(m, _sub_m(j2, q))
        for m2, c in m.items():
            if m2 >= 0:
                tot[m2] += c - m.get(m2 + 2, 0)
    return Counter({k: v for k, v in tot.items() if v})


def count_csfs(occs, p=None, j2=None):
    n = 0
    for o in occs:
        if p is not None and parity(o) != p:
            continue
        c = csf_by_2j(tuple(o[x] for x in ORBITALS))
        n += sum(c.values()) if j2 is None else c.get(j2, 0)
    return n


# ------------------------------------------------------------------ worker --
def tps_for_mu(mu):
    n_au = NPS_A3 * RBOHR_ANG ** 3
    return (1.0 / mu) ** 2 * 4.0 * math.pi * n_au * (1.0 + UPS) * FAC_HARTREE_EV


def ladder(mu):
    k = math.ceil(mu / LADDER_STEP - 1e-9)
    return [round(LADDER_STEP * i, 10) for i in range(1, k)] + [mu]


def run_case(ion, mu, model, wd, hxs_iter0=None):
    from pfac import fac
    print("IMPORTED", fac.__file__, flush=True)
    os.chdir(wd)
    fac.SetAtom("Fe")
    refs = REFS[ion]
    if model == "ref":
        fac.Config(*refs, group="ref")
        opt = ["ref"]
    else:
        fac.Config(refs[0], group="g0")
        fac.Config(*refs[1:], group="ref")
        opt = ["g0"]
    fac.OptimizeRadial(opt)
    print("STAGE free-ion SCF done", flush=True)
    if mu > 0:
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        if hxs_iter0 is not None:
            fac.SetOption("radial:hxs_iter0", hxs_iter0)
            print(f"OPTION radial:hxs_iter0 = {hxs_iter0}", flush=True)
        zps = 26.0 - NELE[ion]
        for m in ladder(mu):
            t = tps_for_mu(m)
            fac.PlasmaScreen(zps, NPS_A3, t, 1, UPS, 1)
            print(f"STAGE ladder mu={m} PlasmaScreen({zps}, {NPS_A3}, {t!r}, 1, {UPS}, 1)", flush=True)
            fac.OptimizeRadial(opt)
            print(f"STAGE ladder mu={m} SCF done", flush=True)
    fac.GetPotential("final.pot")
    groups = ["ref"] if model == "ref" else ["g0", "ref"]
    fac.Structure("e.lev.b", groups)
    fac.MemENTable("e.lev.b")
    fac.PrintTable("e.lev.b", "e.lev", 1)
    print("DONE", flush=True)


# ------------------------------------------------------------------ driver --
def parse_lev(p):
    t = Path(p).read_text().splitlines()
    out = {}
    for l in t:
        if l.startswith("E0"):
            out["E0_eV"] = float(l.split(",")[1])
        elif l.startswith("NELE"):
            out["nele"] = int(l.split("=")[1])
        elif l.startswith("NLEV"):
            out["nlev"] = int(l.split("=")[1])
    rows = [l.split() for l in t if l.strip() and l.split()[0].isdigit()]
    g = min(rows, key=lambda r: float(r[2]))
    out["ground_2j"], out["ground_p"], out["ground_dE"] = int(g[5]), int(g[3]), float(g[2])
    out["ground_label"] = " ".join(g[6:])
    return out


def parse_pot(p):
    hdr, rows = {}, []
    for l in open(p):
        if l.startswith("#"):
            m = re.match(r"#\s+(\w+) =\s+(\S+)$", l)
            if m:
                hdr[m.group(1)] = m.group(2)
        elif l.strip():
            rows.append(l.split())
    last = rows[-1]
    # columns: i r Z Z-Zeff r*Vc r*U ... ZPS(17); q_inf = r(Vc+U) + ZPS at r_max
    return {"mps": int(hdr["mps"]), "dps": float(hdr["dps"]), "ips": int(hdr["ips"]),
            "maxrp": int(hdr["maxrp"]), "rmax": float(last[1]),
            "q_inf": float(last[4]) + float(last[5]) + float(last[17])}


def one(ion, mu, model, out):
    wd = out / "runs" / model / ion.replace(" ", "") / f"mu{mu:g}"
    rec = {"model": model, "ion": ion, "mu": mu, "workdir": str(wd)}
    if wd.exists():
        rec.update(status="FAILED", problem="workdir exists; refusing to overwrite")
        return rec
    wd.mkdir(parents=True)
    t0 = time.time()
    with open(wd / "out.log", "w") as log:
        rc = subprocess.run([sys.executable, "-u", __file__, "--case", ion, repr(mu), model, str(wd)]
                            + ([repr(HXS_ITER0)] if HXS_ITER0 is not None else []),
                            stdout=log, stderr=subprocess.STDOUT, env=os.environ.copy()).returncode
    rec["runtime_s"] = round(time.time() - t0, 1)
    text = (wd / "out.log").read_text()
    rec["imported"] = next((l.split()[1] for l in text.splitlines() if l.startswith("IMPORTED")), None)
    probs = [l for l in text.splitlines()
             if re.search(r"Maximum iteration|Max iteration|rror|nan|inf\b", l)]
    rec["solver_messages"] = probs[:4]
    stages = [l for l in text.splitlines() if l.startswith("STAGE")]
    rec["last_stage"] = stages[-1] if stages else "none"
    if rc != 0 or "DONE" not in text or not (wd / "e.lev").exists():
        rec.update(status="FAILED", problem=f"rc={rc}; last stage: {rec['last_stage']}")
        return rec
    if any("Maximum iteration reached in OptimizeRadial" in l for l in probs):
        rec.update(status="FAILED", problem="SCF not converged")
        return rec
    rec.update(parse_lev(wd / "e.lev"))
    rec.update(parse_pot(wd / "final.pot"))
    errs = []
    if rec.get("nele") != NELE[ion]:
        errs.append(f"NELE {rec.get('nele')} != {NELE[ion]}")
    if not math.isfinite(rec.get("E0_eV", float("nan"))):
        errs.append("E0 not finite")
    if mu > 0 and HXS_ITER0 is not None and f"OPTION radial:hxs_iter0 = {HXS_ITER0}" not in text:
        errs.append("hxs_iter0 requested but not applied by the worker")
    if mu > 0:
        if rec["mps"] != 1 or abs(rec["dps"] * mu - 1) > 1e-6:
            errs.append(f"screening mismatch mps={rec['mps']} dps={rec['dps']}")
        if abs(rec["q_inf"]) > 1e-3:
            errs.append(f"r*V(r_max) = {rec['q_inf']:.4g}: e-e screening not active in SCF")
        if rec["ips"] != 0:
            errs.append(f"ion-sphere cut at ips={rec['ips']}")
    else:
        zion = 26 - NELE[ion]
        if abs(rec["q_inf"] + (zion + 1)) > 1e-3:
            errs.append(f"unscreened r*V(r_max) = {rec['q_inf']:.4g}, expected {-(zion + 1)}")
    rec["status"] = "FAILED" if errs else "ok"
    rec["problem"] = "; ".join(errs)
    return rec


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--case":
        ion, mu, model, wd = sys.argv[2:6]
        run_case(ion, float(mu), model, wd,
                 float(sys.argv[6]) if len(sys.argv) > 6 else None)
        return
    args = sys.argv[1:]
    out = Path(args[args.index("--out") + 1] if "--out" in args else ROOT / "verify_ip_deprince_out")
    jobs = int(args[args.index("--jobs") + 1]) if "--jobs" in args else 4
    global HXS_ITER0
    if "--hxs-iter0" in args:
        HXS_ITER0 = float(args[args.index("--hxs-iter0") + 1])
    print("OPTION radial:hxs_iter0 =", HXS_ITER0, "(None = not set)", flush=True)
    from pfac import fac
    ext = Path(fac.__file__).resolve()
    print("IMPORTED", ext, flush=True)
    if (ROOT / "build") not in ext.parents:
        sys.exit(f"ABORT: pfac extension {ext} is not this worktree's build ({ROOT / 'build'})")
    if out.exists():
        sys.exit(f"ABORT: {out} exists; refusing to overwrite")
    check_refs()
    out.mkdir(parents=True)

    # Stage B bookkeeping: CSF counts
    csf_rows = []
    for ion in REFS:
        refs = [parse_cfg(c) for c in REFS[ion]]
        sd = sd_space(refs)
        csf_rows.append({"ion": ion, "n_ref_configs": len(refs),
                         "ref_csfs_allJ_bothP": count_csfs(refs),
                         "sd_configs": len(sd), "sd_csfs_allJ_bothP": count_csfs(sd),
                         "sd_csfs_allJ_groundP": count_csfs(sd, parity(refs[0])),
                         "published_csfs": PUBLISHED_CSF[ion]})

    cases = [(ion, mu, model) for model in MODELS for ion in REFS for mu in MU_VALUES]
    with ThreadPoolExecutor(jobs) as ex:
        recs = list(ex.map(lambda c: one(*c, out), cases))
    for r in recs:
        print(json.dumps({k: r.get(k) for k in ("model", "ion", "mu", "status", "E0_eV", "nlev",
                                                "ground_2j", "ground_p", "dps", "q_inf",
                                                "problem")}), flush=True)
    for row in csf_rows:  # ground (P, 2J) block of the SD space, from the computed ground level
        g = next((r for r in recs if r["ion"] == row["ion"] and r["model"] == "ref"
                  and r["mu"] == 0 and r["status"] == "ok"), None)
        if g:
            refs = [parse_cfg(c) for c in REFS[row["ion"]]]
            row["ground_block"] = f"P={g['ground_p']} 2J={g['ground_2j']}"
            row["sd_csfs_ground_block"] = count_csfs(sd_space(refs), g["ground_p"], g["ground_2j"])
            row["fac_nlev_ref"] = g["nlev"]

    E = {(r["model"], r["ion"], r["mu"]): r for r in recs}
    ip_rows = []
    for model in MODELS:
        for lo, hi in IPS:
            ip0 = None
            for mu in MU_VALUES:
                a, b = E[(model, lo, mu)], E[(model, hi, mu)]
                row = {"model": model, "ip": f"{lo} -> {hi}", "mu": mu,
                       "deprince_eV": TABLE2[(lo, mu)]}
                if a["status"] != "ok" or b["status"] != "ok":
                    row.update(status="NOT COMPUTED",
                               note=f"{lo}: {a['status']} {a.get('problem', '')}; "
                                    f"{hi}: {b['status']} {b.get('problem', '')}")
                elif a["dps"] != b["dps"] or a["mps"] != b["mps"] or a["imported"] != b["imported"]:
                    row.update(status="NOT COMPUTED", note="inconsistent screening or binary")
                else:
                    ip = b["E0_eV"] - a["E0_eV"]  # FAC E0 is in eV (dbase.c: x HARTREE_EV)
                    row.update(status="ok", fac_eV=round(ip, 4), fac_Ha=ip / FAC_HARTREE_EV,
                               diff_eV=round(ip - TABLE2[(lo, mu)], 4))
                    if mu == 0:
                        ip0 = ip
                    if ip0 is not None:
                        row["dIP_eV"] = round(ip - ip0, 4)
                        row["deprince_dIP_eV"] = round(TABLE2[(lo, mu)] - TABLE2[(lo, 0.0)], 2)
                ip_rows.append(row)

    with open(out / "energies.csv", "w", newline="") as f:
        keys = ["model", "ion", "mu", "status", "E0_eV", "nele", "nlev", "ground_p", "ground_2j",
                "ground_label", "mps", "dps", "ips", "q_inf", "rmax",
                "runtime_s", "problem", "solver_messages", "imported", "workdir"]
        w = csv.DictWriter(f, keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(recs)
    for name, rows in (("verify_ip_deprince_results.csv", ip_rows), ("csf_counts.csv", csf_rows)):
        with open(out / name, "w", newline="") as f:
            keys = list(dict.fromkeys(k for r in rows for k in r))
            w = csv.DictWriter(f, keys)
            w.writeheader()
            w.writerows(rows)
    (out / "provenance.json").write_text(json.dumps({
        "extension": str(ext), "python": sys.version, "git_head": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
        "nps_A3": NPS_A3, "ups": UPS, "ladder_step": LADDER_STEP, "models": MODELS,
        "hxs_iter0": HXS_ITER0}, indent=1))
    print("\nIP table:")
    for r in ip_rows:
        print(json.dumps(r))
    print("\nCSF counts:")
    for r in csf_rows:
        print(json.dumps(r))
    bad = [r for r in ip_rows if r["status"] != "ok"]
    print(f"\n{len(ip_rows) - len(bad)}/{len(ip_rows)} IPs computed; results in {out}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
