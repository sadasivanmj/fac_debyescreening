#!/usr/bin/env python3
"""Run one PFAC case in the current directory (called by run_regression.py).

argv[1] = JSON spec: configs, mu, nps, zps, debye_mode, ee_screen,
config_energy, expect_ext.  Prints CHECK markers so the last completed
stage is visible even if FAC aborts.
"""
import json
import math
import sys

spec = json.loads(sys.argv[1])
from pfac import fac  # noqa: E402

if not fac.__file__.startswith(spec["expect_ext"]):
    sys.exit(f"WRONG EXTENSION IMPORTED: {fac.__file__}")
print("IMPORTED", fac.__file__, flush=True)

RB, H = 0.52917721067, 27.211386018
fac.SetAtom("Fe")
fac.SetOptimizePrint(1)
if spec.get("debye_mode") is not None:
    fac.SetOption("orbital:debye_mode", spec["debye_mode"])
if spec.get("ee_screen") is not None:
    fac.SetOption("radial:ee_screen", spec["ee_screen"])
fac.Config(*spec["configs"], group="g1")
mu = spec.get("mu", 0.0)
if mu > 0:
    nps = spec["nps"]
    tps = (1.0 / mu) ** 2 * 4.0 * math.pi * nps * RB**3 * H
    fac.PlasmaScreen(spec["zps"], nps, tps, 1, 0.0, 1)
    print(f"PLASMA zps={spec['zps']} nps={nps} tps={tps!r} lambda={1/mu}", flush=True)
if spec.get("config_energy"):
    print("CHECK 1: ConfigEnergy(0)", flush=True)
    fac.ConfigEnergy(0)
print("CHECK 2: OptimizeRadial", flush=True)
fac.OptimizeRadial(["g1"])
fac.GetPotential("final.pot")
if spec.get("config_energy"):
    print("CHECK 3: ConfigEnergy(1)", flush=True)
    fac.ConfigEnergy(1)
print("CHECK 4: Structure", flush=True)
fac.Structure("e.lev.b", ["g1"])
print("CHECK 5: tables", flush=True)
fac.MemENTable("e.lev.b")
fac.PrintTable("e.lev.b", "e.lev", 1)
print("DONE", flush=True)
