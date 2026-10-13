#!/usr/bin/env python3
"""Fe II-IV ground energies and IPs in a single+double (SD) active space, Debye screening.

Deprince et al., A&A 643, A57 (2020): Table 1 references, Table 2 IPs.

usage:
  python3 verify_ip_deprince_full_as.py [--out DIR] [--jobs N] [--threads T]
          [--space sd|ref] [--only ION:MU:MODEL,...] [--hxs-iter0 X] [--audit-only]
          [--mpi N]   (OpenMP build only: fac.InitializeMPI(N) in every worker)
          [--mem-budget-gb G]  (jobs start only while the sum of their memory
                                estimates MEM_EST_GB stays <= G; default 11)
  python3 verify_ip_deprince_full_as.py --case JOB.json        (internal worker)

Active space (--space sd, default): for EACH Table 1 reference separately, all
configurations reachable by one or two one-electron moves between the orbitals
{1s,2s,2p,3s,3p,3d,4s} (Pauli-allowed), with a per-reference `seen` set; the
global space is the union of the per-reference sets.  This is the paper's
recipe as stated; it does NOT reproduce the published CSF counts (see the
audit).  --space ref uses the reference configurations only.

CI: FAC Structure over all configuration groups, restricted to EVEN parity and
all J (Structure(0, [-1])); the ground energy E0 is the lowest level of all even
J blocks.  Odd-parity configurations are omitted from the CI groups (they
need a 2p or 3p hole); odd reference configurations stay in the orbital-
optimization group, so the orbital models are unchanged.

Orbital models (identical to verify_ip_deprince.py):
  ref   OptimizeRadial on the average of ALL reference configurations
  gref  OptimizeRadial on the ground reference configuration only
Protocol per job (one FAC process): free-ion OptimizeRadial (fixes the grid);
for mu > 0: orbital:debye_mode 1, radial:ee_screen 1, radial:hxs_iter0 X, and
PlasmaScreen(zps, 1e-6 A^-3, T(mu_k), 1, 0.0, 1) + OptimizeRadial for
mu_k = 0.05, 0.10, ..., mu; GetPotential; Structure; PrintTable.

Restart: a job is reused only if its result.json has status ok and the same
input hash (job definition + configuration lists + extension sha256 + git
HEAD).  Any other existing job directory is renamed to <dir>.attempt<k>
(kept) before the job is rerun.
"""

import csv
import hashlib
import json
import math
import os
import re
import resource
import subprocess
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FAC_HARTREE_EV = 27.211386018     # consts.h; FAC prints level energies in eV
HARTREE_EV = 27.211386            # requested for unit conversion of differences
RBOHR_ANG = 0.52917721067         # consts.h _RBOHR
NPS_A3, UPS, LADDER_STEP = 1e-6, 0.0, 0.05
MU_VALUES = [0.0, 0.1, 0.25]
MODELS = ["ref", "gref"]
IONS = ["Fe II", "Fe III", "Fe IV"]
NELE = {"Fe II": 25, "Fe III": 24, "Fe IV": 23}
ORBITALS = ["1s", "2s", "2p", "3s", "3p", "3d", "4s"]
CAP = {"s": 2, "p": 6, "d": 10}
LVAL = {"s": 0, "p": 1, "d": 2}
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
    "Fe IV": ["1s2 2s2 2p6 3s2 3p6 3d5", "1s2 2s2 2p6 3s2 3p5 3d6",
              "1s2 2s2 2p5 3s2 3p6 3d6", "1s1 2s2 2p6 3s2 3p6 3d6"],
}
PUBLISHED_CSF = {"Fe II": 12268, "Fe III": 31501, "Fe IV": 57189}
TABLE2 = {("Fe II", 0.0): 18.89, ("Fe II", 0.1): 14.20, ("Fe II", 0.25): 7.78,
          ("Fe III", 0.0): 32.93, ("Fe III", 0.1): 25.60, ("Fe III", 0.25): 16.42}
NIST = {"Fe II": 16.19920, "Fe III": 30.651}
IPS = [("Fe II", "Fe III"), ("Fe III", "Fe IV")]
# Upper-bound peak RSS per job (even parity, all J, SD space).  Measured
# 2026-10-13: serial build Fe II 1.19 GB; OpenMP build (InitializeMPI 4) Fe II
# 1.72, Fe III 3.65, Fe IV 9.01 GB.  The serial build needed ~30% less for
# Fe II, so the OpenMP values are used as upper bounds for Fe III / Fe IV.
# Production uses the serial build: the OpenMP build segfaults in PotentialHX
# at the first screened OptimizeRadial (see the report).
MEM_EST_GB = {"Fe II": 2.0, "Fe III": 4.0, "Fe IV": 9.5}


# --------------------------------------------------------- active space --
def parse_cfg(s):
    occ = dict.fromkeys(ORBITALS, 0)
    for tok in s.split():
        m = re.fullmatch(r"(\d[spd])(\d+)", tok)
        if not m or m.group(1) not in occ or occ[m.group(1)]:
            raise ValueError(f"bad or repeated token {tok!r} in {s!r}")
        occ[m.group(1)] = int(m.group(2))
    return occ


def fac_str(occ):
    return " ".join(f"{o}{occ[o]}" for o in ORBITALS if occ[o] > 0)


def key(occ):
    return tuple(occ[o] for o in ORBITALS)


def parity(occ):
    return sum(LVAL[o[1]] * q for o, q in occ.items()) % 2


def sd_of(ref):
    """{key: (occ, excitation order)} for one reference, own seen set."""
    seen = {key(ref): (ref, 0)}
    front = [ref]
    for order in (1, 2):
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
                            seen[key(new)] = (new, order)
                            nxt.append(new)
        front = nxt
    return seen


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
def csf_by_2j(k):
    """jj-coupled CSFs per 2J of one configuration: N(J) = c(M=J) - c(M=J+1)."""
    dists = [[]]
    for orb, q in zip(ORBITALS, k):
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
    return Counter({x: v for x, v in tot.items() if v})


def build_space(ion, space):
    """Validated configuration manifest for one ion."""
    refs = [parse_cfg(c) for c in REFS[ion]]
    entries = {}
    for i, r in enumerate(refs):
        if sum(r.values()) != NELE[ion]:
            raise SystemExit(f"{ion} reference {fac_str(r)}: wrong electron count")
        gen = sd_of(r) if space == "sd" else {key(r): (r, 0)}
        for k, (occ, order) in gen.items():
            e = entries.setdefault(k, {"config": fac_str(occ), "parity": parity(occ),
                                       "is_reference": False, "generated_by": []})
            e["generated_by"].append({"ref": i, "order": order})
            if order == 0:
                e["is_reference"] = True
    for k, e in entries.items():
        occ = dict(zip(ORBITALS, k))
        if sum(k) != NELE[ion] or any(q < 0 or q > CAP[o[1]] for o, q in occ.items()):
            raise SystemExit(f"{ion}: invalid configuration {e['config']}")
        if parse_cfg(e["config"]) != occ:
            raise SystemExit(f"{ion}: round-trip mismatch {e['config']}")
        e["csf_by_2j"] = {str(j): n for j, n in sorted(csf_by_2j(k).items())}
    if len({e["config"] for e in entries.values()}) != len(entries):
        raise SystemExit(f"{ion}: duplicate configuration strings")
    for r in refs:
        if key(r) not in entries:
            raise SystemExit(f"{ion}: reference {fac_str(r)} missing")
    return [entries[k] for k in sorted(entries, reverse=True)]


def dims(entries, p=None):
    out = Counter()
    for e in entries:
        if p is None or e["parity"] == p:
            for j, n in e["csf_by_2j"].items():
                out[int(j)] += n
    return dict(sorted(out.items()))


# -------------------------------------------------------------- worker --
def tps_for_mu(mu):
    return (1.0 / mu) ** 2 * 4.0 * math.pi * NPS_A3 * RBOHR_ANG ** 3 * (1.0 + UPS) * FAC_HARTREE_EV


def ladder(mu):
    k = math.ceil(mu / LADDER_STEP - 1e-9)
    return [round(LADDER_STEP * i, 10) for i in range(1, k)] + [mu]


def stamp(msg):
    print(f"STAGE {time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}", flush=True)


def worker(job_file):
    job = json.loads(Path(job_file).read_text())
    from pfac import fac
    print("IMPORTED", fac.__file__, flush=True)
    os.chdir(job["workdir"])
    ion, mu, model = job["ion"], job["mu"], job["model"]
    refs, ci = job["ref_configs"], job["ci_configs"]
    if job.get("mpi"):
        fac.InitializeMPI(job["mpi"])
        print("OPTION InitializeMPI", job["mpi"], flush=True)
    fac.SetAtom("Fe")
    fac.SetOptimizePrint(1)
    if model == "ref":
        fac.Config(*refs, group="ref")
        opt, groups = ["ref"], ["ref"]
    else:
        fac.Config(refs[0], group="g0")
        fac.Config(*refs[1:], group="ref")
        opt, groups = ["g0"], ["g0", "ref"]
    extra = [c for c in ci if c not in refs]
    for i in range(0, len(extra), 50):
        fac.Config(*extra[i:i + 50], group="sd")
    if extra:
        groups.append("sd")
    stamp(f"configs ref={len(refs)} extra={len(extra)} groups={groups}")
    fac.OptimizeRadial(opt)
    stamp("free-ion SCF done")
    if mu > 0:
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        fac.SetOption("radial:hxs_iter0", job["hxs_iter0"])
        print(f"OPTION radial:hxs_iter0 = {job['hxs_iter0']}", flush=True)
        for m in ladder(mu):
            t = tps_for_mu(m)
            fac.PlasmaScreen(26.0 - NELE[ion], NPS_A3, t, 1, UPS, 1)
            stamp(f"ladder mu={m} T={t!r}")
            fac.OptimizeRadial(opt)
            stamp(f"ladder mu={m} SCF done")
    fac.GetPotential("final.pot")
    # even parity, all J.  [-1]: SetSymmetry treats j[0] < 0 as "all J"; an empty
    # list must not be used (PFAC IntFromList leaves the pointer uninitialised and
    # PStructure frees it -> munmap_chunk() abort).
    fac.Structure(0, [-1])
    stamp("Structure start")
    fac.Structure("e.lev.b", groups)
    stamp("Structure done")
    fac.MemENTable("e.lev.b")
    fac.PrintTable("e.lev.b", "e.lev", 1)
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print("PEAK_RSS_KB", ru.ru_maxrss, flush=True)
    print("CPU_S", ru.ru_utime + ru.ru_stime, flush=True)
    print("DONE", flush=True)


# -------------------------------------------------------------- driver --
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


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
    out.update(ground_2j=int(g[5]), ground_p=int(g[3]), ground_label=" ".join(g[6:]))
    out["fac_block_dims"] = dict(sorted(Counter(int(r[5]) for r in rows).items()))
    out["fac_parities"] = sorted({int(r[3]) for r in rows})
    out["lowest_per_2j_eV"] = {}
    for r in rows:
        j, e = int(r[5]), float(r[2])
        if j not in out["lowest_per_2j_eV"] or e < out["lowest_per_2j_eV"][j]:
            out["lowest_per_2j_eV"][j] = e
    return out


def parse_pot(p):
    hdr, last = {}, None
    for l in open(p):
        if l.startswith("#"):
            m = re.match(r"#\s+(\w+) =\s+(\S+)$", l)
            if m:
                hdr[m.group(1)] = m.group(2)
        elif l.strip():
            last = l.split()
    return {"mps": int(hdr["mps"]), "dps": float(hdr["dps"]), "ips": int(hdr["ips"]),
            "HXS": float(hdr["HXS"]), "maxrp": int(hdr["maxrp"]), "rmax": float(last[1]),
            "q_inf": float(last[4]) + float(last[5]) + float(last[17])}


def run_job(job, out, env_extra):
    wd = Path(job["workdir"])
    res_file = wd / "result.json"
    if res_file.exists():
        old = json.loads(res_file.read_text())
        if old.get("status") == "ok" and old.get("input_hash") == job["input_hash"]:
            old["reused"] = True
            return old
    if wd.exists():
        k = 1
        while wd.with_name(wd.name + f".attempt{k}").exists():
            k += 1
        wd.rename(wd.with_name(wd.name + f".attempt{k}"))
    wd.mkdir(parents=True)
    (wd / "job.json").write_text(json.dumps(job, indent=1))
    rec = {k: job[k] for k in ("ion", "mu", "model", "space", "hxs_iter0", "input_hash",
                               "n_ci_configs", "counted_even_dims", "extension_sha256",
                               "git_head")}
    rec.update(status="running", start=time.strftime("%Y-%m-%dT%H:%M:%S"), workdir=str(wd))
    (wd / "manifest.json").write_text(json.dumps(rec, indent=1))
    env = os.environ.copy()
    env.update(env_extra)
    t0 = time.time()
    with open(wd / "out.log", "w") as log:
        rc = subprocess.run([sys.executable, "-u", __file__, "--case", str(wd / "job.json")],
                            stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    rec.update(end=time.strftime("%Y-%m-%dT%H:%M:%S"), wall_s=round(time.time() - t0, 1), rc=rc,
               omp_num_threads=env.get("OMP_NUM_THREADS"))
    text = (wd / "out.log").read_text()
    lines = text.splitlines()
    rec["imported"] = next((l.split()[1] for l in lines if l.startswith("IMPORTED")), None)
    rec["peak_rss_MB"] = next((round(int(l.split()[1]) / 1024, 1) for l in lines
                               if l.startswith("PEAK_RSS_KB")), None)
    rec["cpu_s"] = next((round(float(l.split()[1]), 1) for l in lines if l.startswith("CPU_S")), None)
    rec["stages"] = [l[6:] for l in lines if l.startswith("STAGE")]
    rec["messages"] = [l for l in lines if re.search(r"Maximum iteration|Max iteration|rror|"
                                                     r"No Debye|nan\b", l)][:6]
    errs = []
    if rc != 0 or "DONE" not in text or not (wd / "e.lev").exists():
        errs.append(f"rc={rc}; last stage: {rec['stages'][-1] if rec['stages'] else 'none'}")
    elif any("Maximum iteration reached in OptimizeRadial" in m for m in rec["messages"]):
        errs.append("SCF not converged")
    else:
        rec.update(parse_lev(wd / "e.lev"))
        rec.update(parse_pot(wd / "final.pot"))
        rec["e_lev_sha256"] = sha256(wd / "e.lev")
        if rec["nele"] != NELE[job["ion"]]:
            errs.append(f"NELE {rec['nele']}")
        if not math.isfinite(rec["E0_eV"]):
            errs.append("E0 not finite")
        if rec["fac_parities"] != [0]:
            errs.append(f"parities {rec['fac_parities']}")
        counted = {int(j): n for j, n in job["counted_even_dims"].items()}
        if rec["fac_block_dims"] != counted:
            errs.append("FAC block dimensions differ from counted CSFs")
        if rec["HXS"] != 1.0:
            errs.append(f"final HXS {rec['HXS']}")
        mu = job["mu"]
        if mu > 0:
            if rec["mps"] != 1 or abs(rec["dps"] * mu - 1) > 1e-6:
                errs.append(f"screening mismatch mps={rec['mps']} dps={rec['dps']}")
            if abs(rec["q_inf"]) > 1e-3:
                errs.append(f"r*V(rmax)={rec['q_inf']}: e-e screening not active")
            if rec["ips"] != 0:
                errs.append(f"ion-sphere cut ips={rec['ips']}")
            if f"OPTION radial:hxs_iter0 = {job['hxs_iter0']}" not in text:
                errs.append("hxs_iter0 not applied")
        elif abs(rec["q_inf"] + (27 - NELE[job["ion"]])) > 1e-3:
            errs.append(f"unscreened r*V(rmax)={rec['q_inf']}")
    rec["status"] = "FAILED" if errs else "ok"
    rec["problem"] = "; ".join(errs)
    (wd / "manifest.json").write_text(json.dumps(rec, indent=1))
    if not errs:
        tmp = wd / "result.json.tmp"
        tmp.write_text(json.dumps(rec, indent=1))
        tmp.rename(res_file)          # atomic checkpoint
    print(f"[{rec['end']}] {job['ion']} mu={job['mu']} {job['model']}: {rec['status']} "
          f"E0={rec.get('E0_eV')} wall={rec['wall_s']}s rss={rec['peak_rss_MB']}MB {rec['problem']}",
          flush=True)
    return rec


def main():
    a = sys.argv[1:]
    if a and a[0] == "--case":
        worker(a[1])
        return

    def opt(name, default):
        return a[a.index(name) + 1] if name in a else default

    out = Path(opt("--out", ROOT / "verify_ip_deprince_full_as_out")).resolve()
    jobs_n = int(opt("--jobs", 1))
    threads = opt("--threads", "1")
    space = opt("--space", "sd")
    hxs_iter0 = float(opt("--hxs-iter0", 8))
    only = opt("--only", None)
    mpi = int(opt("--mpi", 0))
    budget = float(opt("--mem-budget-gb", 11))
    script_sha = sha256(Path(__file__).resolve())
    from pfac import fac
    ext = Path(fac.__file__).resolve()
    ext_sha = sha256(ext)
    head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                          text=True).stdout.strip()
    print("IMPORTED", ext, ext_sha[:16], "HEAD", head, flush=True)
    out.mkdir(parents=True, exist_ok=True)

    # ---- active-space construction and audit
    spaces, summary = {}, []
    for ion in IONS:
        ent = build_space(ion, space)
        spaces[ion] = ent
        p0 = parse_cfg(REFS[ion][0])
        even = [e for e in ent if e["parity"] == parity(p0)]
        d_all, d_even = dims(ent), dims(ent, parity(p0))
        f = out / "active_space" / f"{space}_{ion.replace(' ', '')}.json"
        f.parent.mkdir(exist_ok=True)
        f.write_text(json.dumps(ent, indent=1))
        summary.append({"ion": ion, "space": space, "n_configs": len(ent),
                        "n_configs_even": len(even), "csf_all": sum(d_all.values()),
                        "csf_even": sum(d_even.values()), "largest_even_block": max(d_even.values()),
                        "even_dims": d_even, "published_csf": PUBLISHED_CSF[ion]})
        print("SPACE", json.dumps(summary[-1]), flush=True)
    with open(out / f"active_space_{space}_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    if "--audit-only" in a:
        return

    # ---- jobs
    want = None if only is None else {tuple(x.split(":")) for x in only.split(",")}
    jobs = []
    for mu in MU_VALUES:
        for ion in IONS:
            for model in MODELS:
                if want and (ion.replace(" ", ""), f"{mu:g}", model) not in want:
                    continue
                ent = spaces[ion]
                p0 = parity(parse_cfg(REFS[ion][0]))
                ci = [e["config"] for e in ent if e["parity"] == p0 or e["is_reference"]]
                job = {"ion": ion, "mu": mu, "model": model, "space": space,
                       "hxs_iter0": hxs_iter0 if mu > 0 else None,
                       "ref_configs": REFS[ion], "ci_configs": ci, "n_ci_configs": len(ci),
                       "counted_even_dims": {str(j): n for j, n in dims(
                           [e for e in ent if e["config"] in set(ci)], p0).items()},
                       "mpi": mpi, "script_sha256": script_sha, "nps_A3": NPS_A3, "ups": UPS, "ladder": ladder(mu) if mu > 0 else [],
                       "extension": str(ext), "extension_sha256": ext_sha, "git_head": head}
                job["input_hash"] = hashlib.sha256(json.dumps(job, sort_keys=True).encode()).hexdigest()
                job["workdir"] = str(out / "runs" / space / model / ion.replace(" ", "") / f"mu{mu:g}")
                jobs.append(job)
    env_extra = {"OMP_NUM_THREADS": str(threads)}
    cond, used = threading.Condition(), [0.0]

    def budgeted(job):
        need = min(MEM_EST_GB[job["ion"]], budget)
        with cond:
            cond.wait_for(lambda: used[0] + need <= budget + 1e-9)
            used[0] += need
        try:
            return run_job(job, out, env_extra)
        finally:
            with cond:
                used[0] -= need
                cond.notify_all()

    print(f"SCHEDULE jobs={jobs_n} mem_budget_GB={budget} est={MEM_EST_GB} "
          f"OMP_NUM_THREADS={threads} InitializeMPI={mpi}", flush=True)
    with ThreadPoolExecutor(jobs_n) as ex:   # threads only wait on separate FAC processes
        recs = list(ex.map(budgeted, jobs))

    # ---- collect all checkpointed results present for this space
    results = {}
    for f in sorted((out / "runs" / space).glob("*/*/*/result.json")):
        r = json.loads(f.read_text())
        results[(r["model"], r["ion"], r["mu"])] = r
    keys = ["model", "ion", "mu", "status", "E0_eV", "ground_p", "ground_2j", "ground_label",
            "nlev", "n_ci_configs", "fac_block_dims", "dps", "q_inf", "HXS", "hxs_iter0",
            "wall_s", "cpu_s", "peak_rss_MB", "omp_num_threads", "reused", "input_hash", "workdir"]
    with open(out / f"energies_{space}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, keys, extrasaction="ignore")
        w.writeheader()
        for k in sorted(results):
            w.writerow(results[k])
    ip_rows = []
    for model in MODELS:
        for lo, hi in IPS:
            ip0 = None
            for mu in MU_VALUES:
                ra, rb = results.get((model, lo, mu)), results.get((model, hi, mu))
                row = {"model": model, "ip": f"{lo} -> {hi}", "mu": mu, "space": space,
                       "deprince_eV": TABLE2[(lo, mu)]}
                if not (ra and rb and ra["status"] == rb["status"] == "ok"):
                    row["status"] = "UNAVAILABLE"
                elif (ra["extension_sha256"] != rb["extension_sha256"] or ra["dps"] != rb["dps"]
                      or ra["hxs_iter0"] != rb["hxs_iter0"]):
                    row["status"] = "INCONSISTENT"
                else:
                    ip = rb["E0_eV"] - ra["E0_eV"]
                    row.update(status="ok", fac_eV=round(ip, 4),
                               diff_eV=round(ip - TABLE2[(lo, mu)], 4),
                               diff_pct=round(100 * (ip - TABLE2[(lo, mu)]) / TABLE2[(lo, mu)], 2))
                    if mu == 0:
                        ip0 = ip
                        row["nist_eV"] = NIST[lo]
                        row["diff_nist_eV"] = round(ip - NIST[lo], 4)
                    if ip0 is not None:
                        row["dIP_eV"] = round(ip - ip0, 4)
                        row["deprince_dIP_eV"] = round(TABLE2[(lo, mu)] - TABLE2[(lo, 0.0)], 2)
                ip_rows.append(row)
    with open(out / f"ip_{space}.csv", "w", newline="") as f:
        ks = list(dict.fromkeys(k for r in ip_rows for k in r))
        w = csv.DictWriter(f, ks)
        w.writeheader()
        w.writerows(ip_rows)
    for r in ip_rows:
        print("IP", json.dumps(r), flush=True)
    bad = [r for r in recs if r["status"] != "ok"]
    print(f"{len(recs) - len(bad)}/{len(recs)} jobs ok this invocation", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
