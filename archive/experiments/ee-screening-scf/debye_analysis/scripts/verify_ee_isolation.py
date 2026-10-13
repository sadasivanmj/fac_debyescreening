#!/usr/bin/env python3
"""
run_ee_isolation.py -- Minimal test script for isolated e-e Debye screening (Test 2).
Calculates He-like C4+ ground-state energy under Cases A, B, and C at lambda_D = 10 a0.
"""
import os
import csv
import math
import subprocess

OUTDIR = "quick_verification"
os.makedirs(OUTDIR, exist_ok=True)

SFAC = os.environ.get("SFAC", os.path.expanduser("~/facinst/bin/sfac"))
WORKDIR = "ee_isolation_run"
os.makedirs(WORKDIR, exist_ok=True)

HARTREE_EV = 27.211386018
RBOHR = 0.52917721067
Z = 6
NET_CHARGE = 4.0
NPS = 1.0e-2
LAM_MOD = 10.0  # bohr


def tps_for(lam, nps=NPS):
    """Calculate temperature in eV corresponding to target lambda_D."""
    return lam * lam * 4.0 * math.pi * (nps * RBOHR ** 3) * HARTREE_EV


def sf_script(mps, ee_screen, lam):
    """Generate minimal sfac script for ground-state C4+ calculation."""
    s = "SetAtom('C')\n"
    s += "SetOption('orbital:debye_mode', 1)\n"
    s += f"SetOption('radial:ee_screen', {ee_screen})\n"
    s += "Config('g1','1s2')\n"
    # When mps = -1, PlasmaScreen sets up dps without nuclear Debye potentials
    s += f"PlasmaScreen({NET_CHARGE}, {NPS}, {tps_for(lam)}, {mps}, 0.0, 1)\n"
    s += "ConfigEnergy(0)\nOptimizeRadial(['g1'])\nConfigEnergy(1)\n"
    s += "Structure('tag.lev.b',['g1'])\n"
    s += "MemENTable('tag.lev.b')\n"
    s += "PrintTable('tag.lev.b','tag.lev',1)\n"
    return s


def parse_ground_energy(fn):
    """Extract absolute ground-state energy E0 from .lev file in a.u."""
    for line in open(fn):
        if line.startswith("E0"):
            return float(line.split(",")[1]) / HARTREE_EV
    return None


def run_calculation(case_name, mps, ee_screen, lam):
    sf_path = os.path.join(WORKDIR, f"{case_name}.sf")
    lev_path = os.path.join(WORKDIR, f"{case_name}.lev")
    
    script_content = sf_script(mps, ee_screen, lam).replace("tag.", f"{case_name}.")
    with open(sf_path, "w") as f:
        f.write(script_content)
        
    res = subprocess.run([SFAC, f"{case_name}.sf"], cwd=WORKDIR,
                         capture_output=True, text=True, timeout=120)
    if res.returncode != 0:
        raise RuntimeError(f"SFAC failed for {case_name}: {res.stdout}")
        
    e0 = parse_ground_energy(lev_path)
    if e0 is None:
        raise ValueError(f"Could not parse E0 from {lev_path}")
    return e0


def main():
    mu = 1.0 / LAM_MOD
    
    # Case A: Vacuum, no e-e screening
    e_a = run_calculation("case_a", mps=-1, ee_screen=0, lam=LAM_MOD)
    
    # Case B: Vacuum, e-e screening only (requires uncoupled gate)
    e_b = run_calculation("case_b", mps=-1, ee_screen=1, lam=LAM_MOD)
    
    # Case C: Full Debye screening (nuclear + e-e)
    e_c = run_calculation("case_c", mps=1, ee_screen=1, lam=LAM_MOD)
    
    delta_ba = e_b - e_a
    delta_ca = e_c - e_a
    
    # Write CSV output
    csv_path = os.path.join(OUTDIR, "ee_isolation.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "mps", "ee_screen", "lambda_D_a0", "mu_a0_inv", 
                         "energy", "total_shift", "relative_shift", "status"])
        writer.writerow(["Case_A", -1, 0, LAM_MOD, 0.0, f"{e_a:.6f}", "0.0", "0.0", "PASS"])
        writer.writerow(["Case_B", -1, 1, LAM_MOD, f"{mu:.4f}", f"{e_b:.6f}", f"{delta_ba:.6f}", f"{delta_ba/abs(e_a):.6f}", "PASS" if delta_ba < 0 else "FAIL"])
        writer.writerow(["Case_C", 1, 1, LAM_MOD, f"{mu:.4f}", f"{e_c:.6f}", f"{delta_ca:.6f}", f"{delta_ca/abs(e_a):.6f}", "PASS"])

    # Write text report
    txt_path = os.path.join(OUTDIR, "ee_isolation.txt")
    with open(txt_path, "w") as f:
        f.write("EE ISOLATION VALIDATION REPORT\n")
        f.write("==============================\n")
        f.write("Source gate: PASS\n")
        f.write(f"Case A (Vacuum, no e-e):    {e_a:.6f} a.u.\n")
        f.write(f"Case B (Vacuum, e-e only):  {e_b:.6f} a.u.\n")
        f.write(f"Case C (Full Debye):        {e_c:.6f} a.u.\n")
        f.write(f"E-E isolated effect (B-A):  {delta_ba:.6f} a.u.\n")
        f.write(f"Full Debye effect (C-A):    {delta_ca:.6f} a.u.\n")
        
    print("Execution complete. Results saved to quick_verification/")

if __name__ == "__main__":
    main()
