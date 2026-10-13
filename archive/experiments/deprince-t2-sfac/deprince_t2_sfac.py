#!/usr/bin/env python3
"""Reproduce Deprince et al. (2020, A&A 643, A57) Table 2 with native SFAC.

Ionization potentials of Fe II-Fe VIII for Debye screening parameters
mu = 0, 0.1 and 0.25 a0^-1, with both electron-nucleus and electron-electron
Debye screening.  Every (model, ion, mu) case runs as a separate native
``sfac`` process on a generated ``.sf`` file in its own directory.  No PFAC.

Computational protocol (identical for every ion and every mu, see REPORT.md):

1. Radial orbitals are optimized for the *free* ion, on the average of the
   Table 1 reference configurations (models ref, sdv, sdp) or on the ground
   reference alone (gref, gsdv, gsdp).  This also fixes the radial grid, so
   all mu values of one ion share the same grid.
2. For mu > 0, Debye screening of the e-n and e-e interactions is switched on
   (``orbital:debye_mode`` 1, ``radial:ee_screen`` 1) and the orbitals are
   re-optimized along a ladder mu = 0.05, 0.10, ... up to the target, each
   step warm-started from the previous one (``radial:hxs_iter0``).
3. ``Structure`` diagonalizes the Hamiltonian for the model's configuration
   list, and E0 is the lowest level, printed by ``PrintTable`` in eV.

Models (orbital optimization / CI space):
    ref   orbitals optimized on the statistical-weight average of all Table 1
          reference configurations (FAC's analogue of an EAL over the
          references); CI over the references (all J, both parities).
    sdv   ref orbitals; CI over references + all single and double
          excitations of the ground reference within {3s, 3p, 3d, 4s} with the
          ground parity, restricted to the ground (parity, 2J) block.
    gref  orbitals optimized on the ground reference configuration only;
          CI as in ref.
    gsdv  gref orbitals; CI as in sdv.
    sdp   ref orbitals; CI over the paper's active space: references + all
          single and double excitations of every Table 1 reference into
          {1s, 2s, 2p, 3s, 3p, 3d, 4s}, restricted to the ground (parity, 2J).
    gsdp  gref orbitals; CI as in sdp.

Usage:
    python deprince_t2_sfac.py [--models ref,sdv] [--ions FeII,FeIII]
        [--mus 0,0.1,0.25] [--jobs 8] [--force] [--sfac PATH]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SFAC = "/data/project/fac_worktrees/scf/inst/bin/sfac"

RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018

# PlasmaScreen(zps, nps[A^-3], tps[eV], m=1, ups=0, vxf=1).  In Debye mode
# (m=1) FAC uses lambda_D = sqrt(T/(4 pi n (1+ups))) [a.u.], so only T/n
# matters.  A small density keeps the ion-sphere radius rps = (3 zps /
# (4 pi n))^(1/3) far outside the ion: DensityToSZ truncates the screened
# Hartree integral at rps when ups == 0 (see REPORT.md, failure cause 1).
NPS_A3 = 1.0e-6
UPS = 0.0
VXF = 1
LADDER_STEP = 0.05
HXS_ITER0 = 8

ORBITALS = ["1s", "2s", "2p", "3s", "3p", "3d", "4s"]
CAPACITY = {"s": 2, "p": 6, "d": 10}
VALENCE = ["3s", "3p", "3d", "4s"]

# Table 1 of Deprince et al.; [nl] is a hole, written out explicitly.
REFS = {
    "FeII": ["3d6 4s1", "[3p] 3d7 4s1", "[3p] 3d6 4s2", "[2p] 3d7 4s1",
             "[2p] 3d6 4s2", "[1s] 3d7 4s1", "[1s] 3d6 4s2"],
    "FeIII": ["3d6", "[3p] 3d7", "[3p] 3d6 4s1", "[3p] 3d5 4s2", "[2p] 3d7",
              "[2p] 3d6 4s1", "[2p] 3d5 4s2", "[1s] 3d7", "[1s] 3d6 4s1",
              "[1s] 3d5 4s2"],
    "FeIV": ["3d5", "[3p] 3d6", "[2p] 3d6", "[1s] 3d6"],
    "FeV": ["3d4", "[3p] 3d5", "[2p] 3d5", "[1s] 3d5"],
    "FeVI": ["3d3", "[3p] 3d4", "[2p] 3d4", "[1s] 3d4"],
    "FeVII": ["3d2", "[3p] 3d3", "[2p] 3d3", "[1s] 3d3"],
    "FeVIII": ["3d1", "[3p] 3d2", "[2p] 3d2", "[1s] 3d2"],
    # Not in Table 1 (needed for IP(Fe VIII)).  Same pattern as Fe IV-VIII:
    # ground [Ar] plus the 3p, 2p and 1s holes with the electron in 3d.
    "FeIX": ["3d0", "[3p] 3d1", "[2p] 3d1", "[1s] 3d1"],
}
NELE = {"FeII": 25, "FeIII": 24, "FeIV": 23, "FeV": 22, "FeVI": 21,
        "FeVII": 20, "FeVIII": 19, "FeIX": 18}
# Expected ground terms: a6D9/2, 5D4, 6S5/2, 5D0, 4F3/2, 3F2, 2D3/2, 1S0.
GROUND_2J = {"FeII": 9, "FeIII": 8, "FeIV": 5, "FeV": 0, "FeVI": 3,
             "FeVII": 4, "FeVIII": 3, "FeIX": 0}
IONS = list(REFS)
TARGET_IONS = IONS[:-1]
NEXT_ION = dict(zip(IONS[:-1], IONS[1:]))
ROMAN = {"FeII": "Fe II", "FeIII": "Fe III", "FeIV": "Fe IV", "FeV": "Fe V",
         "FeVI": "Fe VI", "FeVII": "Fe VII", "FeVIII": "Fe VIII",
         "FeIX": "Fe IX"}

# Deprince et al. Table 2: NIST central value (uncertainty as printed) and
# the published MCDF IPs (eV) for mu = 0, 0.1, 0.25.
REFERENCE = {
    "FeII": ("16.19920", "0.00005", "18.89", "14.20", "7.78"),
    "FeIII": ("30.651", "0.012", "32.93", "25.60", "16.42"),
    "FeIV": ("54.91", "0.04", "55.37", "44.65", "29.64"),
    "FeV": ("75.0", "0.2", "74.92", "61.56", "42.83"),
    "FeVI": ("98.985", "0.015", "98.20", "82.20", "59.79"),
    "FeVII": ("124.98", "0.01", "123.60", "104.96", "78.86"),
    "FeVIII": ("151.060", "0.012", "151.14", "129.85", "99.98"),
}
PUBLISHED_CSF = {"FeII": 12268, "FeIII": 31501, "FeIV": 57189, "FeV": 81237,
                 "FeVI": 85798, "FeVII": 71135, "FeVIII": 43271}
MU_VALUES = [0.0, 0.1, 0.25]
MODELS = ["ref", "sdv", "gref", "gsdv", "sdp", "gsdp"]


# ---------------------------------------------------------------- configs --
def parse_occ(spec: str) -> dict[str, int]:
    """'[2p] 3d7 4s1' -> full occupation dict (core shells filled)."""
    occ = {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 0, "4s": 0}
    for tok in spec.split():
        m = re.fullmatch(r"\[(\d[spd])\]", tok)
        if m:
            occ[m.group(1)] -= 1
            continue
        m = re.fullmatch(r"(\d[spd])(\d+)", tok)
        if not m:
            raise ValueError(f"bad config token {tok!r} in {spec!r}")
        occ[m.group(1)] = int(m.group(2))
    return occ


def occ_to_fac(occ: dict[str, int]) -> str:
    return " ".join(f"{o}{occ[o]}" for o in ORBITALS if occ[o] > 0)


def parity(occ: dict[str, int]) -> int:
    lval = {"s": 0, "p": 1, "d": 2}
    return sum(lval[o[1]] * occ[o] for o in occ) % 2


def sd_excitations(ref: dict[str, int], active: list[str]) -> list[dict]:
    """All configurations within two one-electron moves of ref in active."""
    def moves(occ):
        for a in active:
            if occ[a] == 0:
                continue
            for b in active:
                if b != a and occ[b] < CAPACITY[b[1]]:
                    new = dict(occ)
                    new[a] -= 1
                    new[b] += 1
                    yield new
    key = lambda o: tuple(o[x] for x in ORBITALS)
    seen = {key(ref): ref}
    frontier = [ref]
    for _ in range(2):
        nxt = []
        for occ in frontier:
            for new in moves(occ):
                if key(new) not in seen:
                    seen[key(new)] = new
                    nxt.append(new)
        frontier = nxt
    return list(seen.values())


def model_configs(model: str, ion: str) -> dict[str, list[str]]:
    refs = [parse_occ(s) for s in REFS[ion]]
    for occ in refs:
        if sum(occ.values()) != NELE[ion]:
            raise ValueError(f"{ion}: {occ_to_fac(occ)} has wrong electron count")
    if model in ("ref", "sdv", "sdp"):
        groups = {"ref": [occ_to_fac(o) for o in refs]}
    elif model in ("gref", "gsdv", "gsdp"):
        groups = {"g0": [occ_to_fac(refs[0])],
                  "ref": [occ_to_fac(o) for o in refs[1:]]}
    else:
        raise ValueError(f"unknown model {model}")
    if model.endswith(("sdv", "sdp")):
        p0 = parity(refs[0])
        have = {c for v in groups.values() for c in v}
        if model.endswith("sdv"):
            occs = sd_excitations(refs[0], VALENCE)
        else:
            occs = [o for r in refs for o in sd_excitations(r, ORBITALS)]
        sd = []
        for o in occs:
            c = occ_to_fac(o)
            if parity(o) == p0 and c not in have:
                have.add(c)
                sd.append(c)
        groups["sd"] = sd
    return groups


def opt_group(model: str) -> str:
    return "g0" if model.startswith("g") else "ref"


# --------------------------------------------------------------- sf input --
def mu_tag(mu: float) -> str:
    return "mu" + (f"{mu:.4f}".rstrip("0").rstrip(".")).replace(".", "p")


def tps_for_mu(mu: float) -> float:
    """Temperature [eV] giving lambda_D = 1/mu [a0] at NPS_A3 and UPS."""
    n_au = NPS_A3 * RBOHR_ANG**3
    return (1.0 / mu) ** 2 * 4.0 * math.pi * n_au * (1.0 + UPS) * HARTREE_EV


def ladder(mu: float) -> list[float]:
    k = math.ceil(mu / LADDER_STEP - 1e-9)
    steps = [round(LADDER_STEP * i, 10) for i in range(1, k)]
    return steps + [mu]


def sf_input(model: str, ion: str, mu: float) -> str:
    groups = model_configs(model, ion)
    znet = 26 - NELE[ion]
    lines = ["SetAtom('Fe')"]
    for g, cfgs in groups.items():
        for i in range(0, len(cfgs), 10):
            lines.append(f"Config('{g}', "
                         + ", ".join(repr(c) for c in cfgs[i:i + 10]) + ")")
    opt = f"OptimizeRadial(['{opt_group(model)}'])"
    lines.append(opt)
    if mu > 0:
        lines += ["SetOption('orbital:debye_mode', 1)",
                  "SetOption('radial:ee_screen', 1)",
                  f"SetOption('radial:hxs_iter0', {HXS_ITER0})"]
        for m in ladder(mu):
            lines.append(f"PlasmaScreen({znet:.1f}, {NPS_A3:.10g}, "
                         f"{tps_for_mu(m)!r}, 1, {UPS:.1f}, {VXF})")
            lines.append(opt)
    lines.append("GetPotential('final.pot')")
    names = ", ".join(f"'{g}'" for g in groups)
    if model.endswith(("sdv", "sdp")):
        p0 = parity(parse_occ(REFS[ion][0]))
        lines.append(f"Structure({p0}, [{GROUND_2J[ion]}])")
    lines += [f"Structure('e.lev.b', [{names}])",
              "MemENTable('e.lev.b')",
              "PrintTable('e.lev.b', 'e.lev', 1)", ""]
    return "\n".join(lines)


# ----------------------------------------------------------------- parse --
def parse_lev(path: Path) -> dict:
    text = path.read_text().splitlines()
    out = {"E0_eV": None, "nele": None, "nlev": None}
    rows = []
    for line in text:
        if line.startswith("E0"):
            out["E0_eV"] = float(line.split(",", 1)[1])
        elif line.startswith("NELE"):
            out["nele"] = int(line.split("=")[1])
        elif line.startswith("NLEV"):
            out["nlev"] = int(line.split("=")[1])
        else:
            f = line.split()
            if len(f) >= 6 and f[0].isdigit():
                rows.append((float(f[2]), int(f[3]), int(f[5])))
    if rows:
        e, p, j2 = min(rows)
        out.update(ground_rel_eV=e, ground_parity=p, ground_2j=j2)
    return out


def parse_pot(path: Path) -> dict:
    out = {"orb": {}}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        m = re.match(r"#\s+(\w+) =\s+(\S+)$", line)
        if m:
            out[m.group(1)] = m.group(2)
            continue
        f = line[1:].split()
        if len(f) == 9 and f[0].isdigit():
            out["orb"][f"{f[1]},{f[2]}"] = float(f[8])
    return out


# ------------------------------------------------------------ provenance --
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(sfac: Path) -> dict:
    info = {"sfac": str(sfac), "sfac_sha256": sha256(sfac)}
    root = sfac.resolve().parent
    while root != root.parent and not (root / ".git").exists():
        root = root.parent
    if (root / ".git").exists():
        git = lambda *a: subprocess.run(["git", "-C", str(root), *a],
                                        capture_output=True, text=True).stdout
        diff = git("diff", "--", "faclib", "sfac")
        info.update(source_tree=str(root),
                    git_commit=git("rev-parse", "HEAD").strip(),
                    git_describe=git("log", "-1", "--format=%h %s").strip(),
                    source_diff_sha256=hashlib.sha256(diff.encode()).hexdigest(),
                    source_diff_files=git("diff", "--stat", "--", "faclib",
                                          "sfac").strip())
        lib = root / "libfac.a"
        if lib.exists():
            newer = [str(p.relative_to(root)) for p in
                     list((root / "faclib").glob("*.[ch]")) +
                     list((root / "sfac").glob("*.[ch]"))
                     if p.stat().st_mtime > sfac.stat().st_mtime]
            info["sources_newer_than_binary"] = newer
    return info


# ------------------------------------------------------------------- run --
OUTPUTS = ["e.lev", "e.lev.b", "final.pot", "sfac.log", "run.json",
           "error.pot", "error.bpt", "error.wav"]


def run_case(model, ion, mu, root, sfac, prov, force, timeout):
    wd = root / "runs" / model / ion / mu_tag(mu)
    wd.mkdir(parents=True, exist_ok=True)
    sf = sf_input(model, ion, mu)
    sf_hash = hashlib.sha256(sf.encode()).hexdigest()
    status_file = wd / "run.json"
    if not force and status_file.exists():
        old = json.loads(status_file.read_text())
        if (old.get("status") == "ok" and old.get("sf_sha256") == sf_hash
                and old.get("sfac_sha256") == prov["sfac_sha256"]
                and (wd / "e.lev").exists()
                and (wd / "e.lev").stat().st_mtime >= old["t_start"]):
            old["reused"] = True
            return old
    for name in OUTPUTS:
        (wd / name).unlink(missing_ok=True)
    (wd / "in.sf").write_text(sf)
    groups = model_configs(model, ion)
    rec = {"model": model, "ion": ion, "mu": mu, "workdir": str(wd),
           "sf_sha256": sf_hash, "n_config": sum(map(len, groups.values())),
           "configs": groups, "nps_A3": NPS_A3, "ups": UPS,
           "zps": 26 - NELE[ion], "ladder": ladder(mu) if mu > 0 else [],
           "tps_eV": tps_for_mu(mu) if mu > 0 else None,
           "lambda_target_a0": 1 / mu if mu > 0 else None,
           "debye_mode": 1 if mu > 0 else 0, "ee_screen": 1 if mu > 0 else 0,
           "hxs_iter0": HXS_ITER0 if mu > 0 else 0, **prov,
           "t_start": time.time(), "reused": False}
    try:
        with open(wd / "sfac.log", "w") as log:
            proc = subprocess.run([str(sfac), "in.sf"], cwd=wd, stdout=log,
                                  stderr=subprocess.STDOUT, timeout=timeout)
        rec["returncode"] = proc.returncode
    except subprocess.TimeoutExpired:
        rec["returncode"] = "timeout"
    rec["t_end"] = time.time()
    rec["runtime_s"] = round(rec["t_end"] - rec["t_start"], 1)
    rec["status"], rec["problems"] = "failed", []
    lev, pot = wd / "e.lev", wd / "final.pot"
    if rec["returncode"] != 0:
        rec["problems"].append(f"sfac exit {rec['returncode']}")
        tail = (wd / "sfac.log").read_text().strip().splitlines()[-3:]
        rec["log_tail"] = tail
    elif not lev.exists() or lev.stat().st_mtime < rec["t_start"]:
        rec["problems"].append("e.lev missing or stale")
    else:
        rec.update(parse_lev(lev))
        p = parse_pot(pot) if pot.exists() else {}
        rec["orbital_energies_Ha"] = p.get("orb", {})
        rec["mps"] = int(p["mps"]) if "mps" in p else None
        rec["dps_a0"] = float(p["dps"]) if "dps" in p else None
        rec["ips"] = int(p["ips"]) if "ips" in p else None
        if rec["E0_eV"] is None or not math.isfinite(rec["E0_eV"]):
            rec["problems"].append("E0 missing or non-finite")
        if rec["nele"] != NELE[ion]:
            rec["problems"].append(f"NELE {rec['nele']} != {NELE[ion]}")
        if mu > 0:
            if rec["mps"] != 1:
                rec["problems"].append(f"mps={rec['mps']}, Debye mode inactive")
            elif abs(rec["dps_a0"] * mu - 1) > 1e-6:
                rec["problems"].append(f"dps={rec['dps_a0']} != 1/mu")
        if rec.get("ground_2j") != GROUND_2J[ion]:
            rec.setdefault("warnings", []).append(
                f"ground 2J={rec.get('ground_2j')} (expected {GROUND_2J[ion]})")
        if not rec["problems"]:
            rec["status"] = "ok"
    status_file.write_text(json.dumps(rec, indent=1))
    return rec


# ---------------------------------------------------------------- tables --
def fmt(x, nd=4):
    return "" if x is None else f"{x:.{nd}f}"


def build_tables(recs, root):
    by = {(r["model"], r["ion"], r["mu"]): r for r in recs}
    models = sorted({r["model"] for r in recs}, key=MODELS.index)
    mus = sorted({r["mu"] for r in recs})
    res = root / "results"
    res.mkdir(exist_ok=True)

    with open(res / "energies.csv", "w", newline="") as fh:
        cols = ["model", "ion", "mu", "status", "E0_eV", "nele", "nlev",
                "n_config", "ground_parity", "ground_2j", "dps_a0", "mps",
                "runtime_s", "problems", "warnings", "sfac", "sfac_sha256",
                "git_commit", "source_diff_sha256", "workdir"]
        w = csv.DictWriter(fh, cols, extrasaction="ignore")
        w.writeheader()
        for r in sorted(recs, key=lambda r: (MODELS.index(r["model"]),
                                             IONS.index(r["ion"]), r["mu"])):
            row = dict(r)
            row["problems"] = "; ".join(r.get("problems", []))
            row["warnings"] = "; ".join(r.get("warnings", []))
            w.writerow(row)

    rows = []
    for model in models:
        for ion in TARGET_IONS:
            nist, nist_unc, *pub = REFERENCE[ion]
            for mu in mus:
                a, b = by.get((model, ion, mu)), by.get((model, NEXT_ION[ion], mu))
                pub_ip = float(pub[MU_VALUES.index(mu)]) if mu in MU_VALUES else None
                row = {"model": model, "ion": ROMAN[ion], "mu": mu,
                       "NIST_IP_eV": nist, "NIST_unc_eV": nist_unc,
                       "published_IP_eV": pub[MU_VALUES.index(mu)]
                       if mu in MU_VALUES else ""}
                if a is None or b is None:
                    row["status"] = "NOT RUN"
                elif a["status"] != "ok" or b["status"] != "ok":
                    bad = [x["ion"] for x in (a, b) if x["status"] != "ok"]
                    row["status"] = "FAILED (" + ",".join(bad) + ")"
                else:
                    same = all(a[k] == b[k] for k in
                               ("sfac_sha256", "debye_mode", "ee_screen",
                                "nps_A3", "ups", "dps_a0", "hxs_iter0"))
                    if not same:
                        row["status"] = "FAILED (incompatible settings)"
                    else:
                        ip = b["E0_eV"] - a["E0_eV"]
                        row.update(status="ok", FAC_IP_eV=ip,
                                   E0_ion_eV=a["E0_eV"], E0_next_eV=b["E0_eV"])
                        if pub_ip is not None:
                            d = ip - pub_ip
                            row.update(err_vs_pub_eV=d,
                                       signed_rel_vs_pub_pct=100 * d / pub_ip,
                                       abs_rel_vs_pub_pct=100 * abs(d) / abs(pub_ip))
                        if mu == 0:
                            d = ip - float(nist)
                            row.update(err_vs_NIST_eV=d,
                                       signed_rel_vs_NIST_pct=100 * d / float(nist),
                                       abs_rel_vs_NIST_pct=100 * abs(d) / float(nist))
                        warn = a.get("warnings", []) + b.get("warnings", [])
                        if warn:
                            row["warnings"] = "; ".join(warn)
                rows.append(row)
    cols = ["model", "ion", "mu", "status", "NIST_IP_eV", "NIST_unc_eV",
            "published_IP_eV", "FAC_IP_eV", "err_vs_pub_eV",
            "signed_rel_vs_pub_pct", "abs_rel_vs_pub_pct", "err_vs_NIST_eV",
            "signed_rel_vs_NIST_pct", "abs_rel_vs_NIST_pct", "E0_ion_eV",
            "E0_next_eV", "warnings"]
    with open(res / "ip_comparison.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    md = [f"# Deprince et al. Table 2 with native SFAC",
          f"Generated {dt.datetime.now().isoformat(timespec='seconds')}.",
          "FAC values are computed here; NIST and published values are quoted "
          "from Deprince et al. (2020), Table 2.", ""]
    ipr = {(r["model"], r["ion"], r["mu"]): r for r in rows}
    cell = lambda r, k, nd: fmt(r.get(k), nd) if r["status"] == "ok" else r["status"]
    for model in models:
        md += [f"## Model `{model}`", "", "### Table A: reproduction of Table 2", "",
               "| Ion | NIST | " + " | ".join(
                   f"Pub μ={mu:g} | FAC μ={mu:g} | δ% | ε%" for mu in mus)
               + " |",
               "|---|---:|" + "---:|---:|---:|---:|" * len(mus)]
        for ion in TARGET_IONS:
            cells = [ROMAN[ion], REFERENCE[ion][0]]
            for mu in mus:
                r = ipr[(model, ROMAN[ion], mu)]
                cells += [r["published_IP_eV"], cell(r, "FAC_IP_eV", 3),
                          cell(r, "signed_rel_vs_pub_pct", 2),
                          cell(r, "abs_rel_vs_pub_pct", 2)]
            md.append("| " + " | ".join(cells) + " |")
        if 0.0 in mus:
            md += ["", "### Table B: free-ion IP vs NIST", "",
                   "| Ion | NIST (eV) | FAC μ=0 (eV) | FAC−NIST (eV) | δ% | ε% |",
                   "|---|---:|---:|---:|---:|---:|"]
            for ion in TARGET_IONS:
                r = ipr[(model, ROMAN[ion], 0.0)]
                md.append(f"| {ROMAN[ion]} | {REFERENCE[ion][0]} | "
                          f"{cell(r, 'FAC_IP_eV', 3)} | {cell(r, 'err_vs_NIST_eV', 3)} | "
                          f"{cell(r, 'signed_rel_vs_NIST_pct', 2)} | "
                          f"{cell(r, 'abs_rel_vs_NIST_pct', 2)} |")
        mono = []
        for ion in TARGET_IONS:
            ips = [(mu, ipr[(model, ROMAN[ion], mu)].get("FAC_IP_eV")) for mu in mus]
            ips = [(mu, ip) for mu, ip in ips if ip is not None]
            ok = all(a[1] > b[1] for a, b in zip(ips, ips[1:]))
            mono.append(f"{ROMAN[ion]}: {'decreasing' if ok else 'NOT monotonic'}"
                        f" ({len(ips)} values)")
        md += ["", "Monotonic decrease of IP with μ: " + "; ".join(mono), ""]
    md += ["## Run status", "", "| model | ion | μ | status | E0 (eV) | NLEV | "
           "ground 2J | runtime (s) | notes |", "|---|---|---:|---|---:|---:|---:|---:|---|"]
    for r in sorted(recs, key=lambda r: (MODELS.index(r["model"]),
                                         IONS.index(r["ion"]), r["mu"])):
        notes = "; ".join(r.get("problems", []) + r.get("warnings", [])
                          + r.get("log_tail", [])[-1:])
        md.append(f"| {r['model']} | {ROMAN[r['ion']]} | {r['mu']:g} | "
                  f"{r['status']}{' (reused)' if r.get('reused') else ''} | "
                  f"{fmt(r.get('E0_eV'), 4)} | {r.get('nlev', '')} | "
                  f"{r.get('ground_2j', '')} | {r.get('runtime_s', '')} | {notes} |")
    (res / "summary.md").write_text("\n".join(md) + "\n")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sfac", default=os.environ.get("SFAC", DEFAULT_SFAC))
    ap.add_argument("--models", default="ref")
    ap.add_argument("--ions", default=",".join(IONS))
    ap.add_argument("--mus", default="0,0.1,0.25")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=6 * 3600)
    ap.add_argument("--root", default=str(HERE))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="write .sf inputs and print config counts only")
    a = ap.parse_args()

    sfac = Path(a.sfac).expanduser().resolve()
    if not (sfac.is_file() and os.access(sfac, os.X_OK)):
        sys.exit(f"sfac not found or not executable: {sfac}")
    root = Path(a.root).resolve()
    models = a.models.split(",")
    ions = a.ions.split(",")
    mus = [float(x) for x in a.mus.split(",")]
    for m in models:
        if m not in MODELS:
            sys.exit(f"unknown model {m}")
    for i in ions:
        if i not in IONS:
            sys.exit(f"unknown ion {i}")
    cases = [(m, i, mu) for m in models for i in ions for mu in mus]
    if a.dry_run:
        for m, i, mu in cases:
            g = model_configs(m, i)
            print(m, i, mu, {k: len(v) for k, v in g.items()})
        return

    prov = provenance(sfac)
    print(json.dumps(prov, indent=1))
    if prov.get("sources_newer_than_binary"):
        print("WARNING: sources newer than binary:", prov["sources_newer_than_binary"])
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(run_case, m, i, mu, root, sfac, prov, a.force,
                          a.timeout): (m, i, mu) for m, i, mu in cases}
        recs = []
        for f in futs:
            r = f.result()
            recs.append(r)
            print(f"[{r['status']:6s}] {r['model']:4s} {r['ion']:7s} mu={r['mu']:<5g}"
                  f" E0={fmt(r.get('E0_eV'))} t={r.get('runtime_s')}s "
                  f"{'; '.join(r.get('problems', []))}", flush=True)
    # merge with earlier records so tables cover everything run so far
    allrecs = {}
    for f in (root / "runs").glob("*/*/*/run.json"):
        r = json.loads(f.read_text())
        allrecs[(r["model"], r["ion"], r["mu"])] = r
    for r in recs:
        allrecs[(r["model"], r["ion"], r["mu"])] = r
    rows = build_tables(list(allrecs.values()), root)
    n = {"ok": 0, "failed": 0}
    for r in allrecs.values():
        n["ok" if r["status"] == "ok" else "failed"] += 1
    expected = {(m, i, mu) for m in MODELS for i in IONS for mu in MU_VALUES}
    print(f"\nruns: {n['ok']} ok, {n['failed']} failed, "
          f"{len(expected - set(allrecs))} of the standard grid not yet run")
    print((root / "results" / "summary.md").read_text())


if __name__ == "__main__":
    main()
