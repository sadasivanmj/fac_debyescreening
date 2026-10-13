# Minimal Fe XXVI (1s1) reproducer with Debye nuclear screening.
#
# Usage: python diag_fe_xxvi.py --mu 0.1 --outdir OUT [--ee 0|1] [--hxs-iter0 8]
#
# The same FAC sequence as the Fe XVII-XXV runs: free-ion OptimizeRadial, then
# PlasmaScreen in steps of 0.05 with OptimizeRadial after each step.
# It prints what FAC actually used: mps, dps (so mu = 1/dps), the Debye nuclear
# term ZPS(r) read back from the potential file, and the ground-state energy.
# For one electron the expected ZPS(r) is Z*(1 - exp(-mu r)).

import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "fe-xvii-xxv-ground-ip"))
from fe_ground_common import (DENSITY_A3, MU_STEP, Z_IRON, git_commit,
                              read_level_table, read_potential_header,
                              temperature_for_mu)

parser = argparse.ArgumentParser()
parser.add_argument("--mu", type=float, required=True)
parser.add_argument("--ee", type=int, default=0, help="radial:ee_screen (irrelevant for one electron)")
parser.add_argument("--hxs-iter0", type=float, default=8.0)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

os.makedirs(args.outdir)          # fails if it exists
os.chdir(args.outdir)

from pfac import fac
if not os.path.realpath(fac.__file__).startswith(os.path.realpath(REPO) + os.sep):
    sys.exit(f"wrong FAC extension imported: {fac.__file__}")
print("FAC extension:", fac.__file__, flush=True)

fac.SetAtom("Fe")
fac.Config("1s1", group="g0")
fac.OptimizeRadial(["g0"])        # free ion; fixes the grid

steps = []
if args.mu > 0:
    fac.SetOption("orbital:debye_mode", 1)
    fac.SetOption("radial:ee_screen", args.ee)
    fac.SetOption("radial:hxs_iter0", args.hxs_iter0)
    steps = [round(MU_STEP * k, 10) for k in range(1, int(round(args.mu / MU_STEP)))] + [args.mu]
    for mu_step in steps:
        fac.PlasmaScreen(Z_IRON - 1, DENSITY_A3, temperature_for_mu(mu_step), 1, 0.0, 1)
        print("screening step mu =", mu_step, flush=True)
        fac.OptimizeRadial(["g0"])

fac.GetPotential("final.pot")
fac.Structure("levels.b", ["g0"])
fac.MemENTable("levels.b")
fac.PrintTable("levels.b", "levels.txt", 1)

# ---- read back ----
header, asymptotic = read_potential_header("final.pot")
rows = [l.split() for l in open("final.pot") if l.strip() and not l.startswith("#")]
rad = [float(r[1]) for r in rows]
zps = [float(r[17]) for r in rows]
dps = float(header["dps"])
mu_used = 1.0 / dps if int(header["mps"]) == 1 and dps > 0 else 0.0
# compare ZPS(r) with Z(1-exp(-mu r)) at a few radii (nuclear charge Z = 26)
checks = []
for k in (100, 500, 900, 1300):
    expected = Z_IRON * (1.0 - math.exp(-mu_used * rad[k]))
    checks.append({"r": rad[k], "ZPS": zps[k], "expected_Z(1-exp(-mu r))": expected})
nele, e0, levels = read_level_table("levels.txt")
out = {
    "ion": "Fe XXVI", "config": "1s1", "electrons": nele,
    "mu_requested": args.mu, "mu_used": mu_used, "mps": int(header["mps"]),
    "ee_screen_option": args.ee if args.mu > 0 else 0,
    "hxs_iter0": args.hxs_iter0 if args.mu > 0 else 0,
    "asymptotic_charge": asymptotic, "ZPS_checks": checks,
    "E0_eV": e0, "n_levels": len(levels),
    "commit": git_commit(), "fac_extension": fac.__file__,
}
json.dump(out, open("diag.json", "w"), indent=2)
print(json.dumps(out, indent=2))
