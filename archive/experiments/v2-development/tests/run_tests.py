#!/usr/bin/env python3
"""Controlled tests of the Debye-screened DensityToSZ (isolated worktree only).

usage: python3 run_tests.py <run_id> <case_prefix ...>
Outputs go to tests/out/<run_id>/<build>/<case>/ (refuses existing dirs).
Builds: base = _baseline_pfac (unpatched c4a016a), patch = _patched_pfac.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path("/tmp/fac-debye-scf-test")
TESTS = ROOT / "tests"
BUILDS = {"base": ROOT / "_baseline_pfac", "patch": ROOT / "_patched_pfac",
          "diag": ROOT / "_diag_pfac", "v2": ROOT / "_patched_v2_pfac"}
HE = ["1s2"]
FEIX = ["1s2 2s2 2p6 3s2 3p6"]
FEIII = ["1s2 2s2 2p6 3s2 3p6 3d6"]

CASES = {
    # Test 1: unscreened regression (no PlasmaScreen)
    "T1_He_free": dict(configs=HE),
    "T1_FeIX_free": dict(configs=FEIX),
    "T1_FeIII_free": dict(configs=FEIII),
    "T1_FeIX_free_ee1": dict(configs=FEIX, debye_mode=1, ee_screen=1),
    # Test 2: minimal screened controls (both screenings on)
    "T2_He_mu0.1": dict(configs=HE, mu=0.1, nps=1e-6, zps=24.0, debye_mode=1, ee_screen=1),
    "T2_He_mu0.25": dict(configs=HE, mu=0.25, nps=1e-6, zps=24.0, debye_mode=1, ee_screen=1),
    "T2_FeIX_mu0.1": dict(configs=FEIX, mu=0.1, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1),
    "T2_FeIX_mu0.25": dict(configs=FEIX, mu=0.25, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1),
    "T2_FeIX_mu0.1_nps0.01": dict(configs=FEIX, mu=0.1, nps=1e-2, zps=8.0, debye_mode=1, ee_screen=1),
    # Test 4: existing Fe III diagnostic settings, only ee_screen changed
    "T4_FeIII_diag_ee0": dict(configs=FEIII, mu=0.1, nps=1e-2, zps=2.0, debye_mode=1, ee_screen=0, config_energy=True),
    "T4_FeIII_diag_ee1": dict(configs=FEIII, mu=0.1, nps=1e-2, zps=2.0, debye_mode=1, ee_screen=1, config_energy=True),
    "T4_FeIII_free_CE": dict(configs=FEIII, config_energy=True),
    # Test 5: edge cases
    "T5_FeIX_mu1e-8": dict(configs=FEIX, mu=1e-8, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1),
    "T5_FeIX_mu1e-4": dict(configs=FEIX, mu=1e-4, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1),
    "T5_FeIX_mu1.0": dict(configs=FEIX, mu=1.0, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1),
    "T5_He_mu1e-4": dict(configs=HE, mu=1e-4, nps=1e-6, zps=24.0, debye_mode=1, ee_screen=1),
    "T5_He_mu1e-3": dict(configs=HE, mu=1e-3, nps=1e-6, zps=24.0, debye_mode=1, ee_screen=1),
}


def pot_header(p):
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
    return h, orb


def main():
    run_id, prefixes = sys.argv[1], sys.argv[2:]
    builds = [b for b in BUILDS if b in prefixes] or ["base", "patch"]
    prefixes = [p for p in prefixes if p not in BUILDS]
    res = []
    for build in builds:
        for name, spec in CASES.items():
            if prefixes and not any(name.startswith(p) for p in prefixes):
                continue
            wd = TESTS / "out" / run_id / build / name
            if wd.exists():
                sys.exit(f"refusing to overwrite {wd}")
            wd.mkdir(parents=True)
            spec = dict(spec, expect_ext=str(BUILDS[build]))
            env = {"PYTHONPATH": str(BUILDS[build]), "PATH": "/usr/bin:/bin",
                   "HOME": str(Path.home())}
            with open(wd / "out.log", "w") as log:
                rc = subprocess.run([sys.executable, str(TESTS / "run_case.py"),
                                     json.dumps(spec)], cwd=wd, env=env,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    timeout=3600).returncode
            log = (wd / "out.log").read_text().splitlines()
            checks = [l for l in log if l.startswith("CHECK") or l == "DONE"]
            rec = {"build": build, "case": name, "rc": rc,
                   "last_stage": checks[-1] if checks else "none",
                   "imported": next((l.split()[1] for l in log if l.startswith("IMPORTED")), None),
                   "scf_iters": sum(1 for l in log if l.startswith("optimize loop")),
                   "errors": [l for l in log if re.search(r"rror|Max iteration|unbound|nan", l)][:3]}
            if (wd / "e.lev").exists():
                t = (wd / "e.lev").read_text()
                rec["E0_eV"] = float(re.search(r"E0\s*=\s*\d+,\s*(\S+)", t).group(1))
            for pf in ("final.pot", "error.pot"):
                if (wd / pf).exists():
                    h, orb = pot_header(wd / pf)
                    rec[pf] = {k: h.get(k) for k in ("mps", "dps", "ips", "rps", "HXS", "maxrp", "rmax")}
                    rec[pf]["orb"] = orb
            res.append(rec)
            print(json.dumps({k: v for k, v in rec.items() if k not in ("final.pot", "error.pot")}), flush=True)
    (TESTS / "out" / run_id / "results.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
