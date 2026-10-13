#!/usr/bin/env python3
"""Run kernel-dump cases with the TEMPORARY diagnostic build (_diag_pfac)."""
import json, subprocess, sys
from pathlib import Path
T = Path("/tmp/fac-debye-scf-test/tests"); EXT = "/tmp/fac-debye-scf-test/_diag_pfac"
FEIX = ["1s2 2s2 2p6 3s2 3p6"]
cases = {"D_FeIX_mu0": dict(configs=FEIX, debye_mode=1, ee_screen=1),
         **{f"D_FeIX_mu{m:g}": dict(configs=FEIX, mu=m, nps=1e-6, zps=8.0, debye_mode=1, ee_screen=1)
            for m in (1e-8, 1e-4, 0.1, 0.25, 1.0)},
         "D_He_mu0.1": dict(configs=["1s2"], mu=0.1, nps=1e-6, zps=24.0, debye_mode=1, ee_screen=1),
         **{f"F_FeIX_mu{m:g}_nps1e-9": dict(configs=FEIX, mu=m, nps=1e-9, zps=8.0, debye_mode=1, ee_screen=1) for m in (1e-4, 0.1, 1.0, 3.0)}}
import sys as _s
if len(_s.argv) > 2: cases = {k: v for k, v in cases.items() if k.startswith(_s.argv[2])}
run = sys.argv[1]
for name, spec in cases.items():
    wd = T / "out" / run / "diag" / name
    if wd.exists(): sys.exit(f"refusing to overwrite {wd}")
    wd.mkdir(parents=True)
    env = {"PYTHONPATH": EXT, "PATH": "/usr/bin:/bin", "FAC_D2SZ_DUMP": str(wd / "d2sz.dump")}
    rc = subprocess.run([sys.executable, str(T / "run_case.py"), json.dumps(dict(spec, expect_ext=EXT))],
                        cwd=wd, env=env, stdout=open(wd / "out.log", "w"), stderr=subprocess.STDOUT).returncode
    log = (wd / "out.log").read_text().splitlines()
    print(name, "rc", rc, [l for l in log if l.startswith(("CHECK", "DONE", "IMPORTED"))][-2:])
