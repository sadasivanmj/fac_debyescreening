#!/usr/bin/env python3
"""
chen2018_kr34_reproduce.py
==========================
Reproduces Table II of Chen et al., Phys. Plasmas 25, 072120 (2018) for Kr34+.
Specifically targets the 1s2p 1P1 excitation energy in a.u.

Models:
    Model A (FAC 'b'): e-n screened, e-e unscreened (ee=0)
    Model B (FAC 'g'): e-n screened, e-e screened (ee=1)
"""

import os
import re
import math
import subprocess
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# --------------------------------------------------------------------------
# FAC CONVENTIONS & CONSTANTS
# --------------------------------------------------------------------------
SFAC       = os.environ.get("SFAC", "sfac")
WORKDIR    = "chen2018_kr34_run"
NPS        = 1.0e-2
RBOHR      = 0.52917721067
HARTREE_EV = 27.211386018

def tps_for(lam):
    """eV temperature that makes dps == lam (bohr), given nps in A^-3 and ups=0."""
    return lam * lam * 4.0 * math.pi * (NPS * RBOHR**3) * HARTREE_EV

# --------------------------------------------------------------------------
# CHEN ET AL (2018) TABLE II DATA (He-like Kr34+)
# Quantities are tabulated as direct Excitation Energies in a.u.
# --------------------------------------------------------------------------
# We omit lambda_D = 2.0 because Chen et al. do not provide FAC ('b' or 'g') 
# values for it, only MCDF ('a') values.
LAMBDAS = [1.0e6, 20.0, 15.0, 10.0, 7.5, 5.0]

# 1s2p 1P1 FAC values
# Model A corresponds to 'b' in the paper
REF_EXC_A = [482.1981, 482.1428, 482.1426, 482.1421, 482.1412, 482.0457]

# Model B corresponds to 'g' in the paper
# Note: At lambda_D = infinity, screening is 0, so Model A == Model B physically.
# We use the free ion 'b' value as the reference for 'g' at infinity.
REF_EXC_B = [482.1981, 482.1418, 482.1414, 482.1403, 482.1390, 482.0424]

# --------------------------------------------------------------------------
# FAC DRIVING FUNCTIONS
# --------------------------------------------------------------------------
def sf_input(tag, lam, ee):
    """Generates the FAC script tailored for Kr34+."""
    tps = tps_for(lam)
    
    # Matching the simpler CI basis typically used in Chen's FAC calculations
    lines = [
        "SetAtom('Kr')",
        "SetOption('orbital:debye_mode', 1)",
        f"SetOption('radial:ee_screen', {ee})",
        "Config('g1', '1s2')",
        "Config('g2', '1s1 2*1')",
        "Config('g3', '1s1 3*1')",
    ]
    groups = "'g1', 'g2', 'g3'"
    
    lines.extend([
        # Nuclear charge for Kr is 36.0
        f"PlasmaScreen(36.0, {NPS}, {tps}, 1, 0.0, 1)",
        "ConfigEnergy(0)",
        "OptimizeRadial(['g1'])",
        "ConfigEnergy(1)",
        f"Structure('{tag}.lev.b', [{groups}])",
        f"PrintTable('{tag}.lev.b', '{tag}.lev', 1)"
    ])
    
    return "\n".join(lines) + "\n"

def parse_lev(fn):
    """Extracts ground state E0 (eV) and level dictionary."""
    E0 = None
    lev = []
    if not os.path.exists(fn):
        return None, None
        
    for L in open(fn):
        if L.startswith("E0"):
            E0 = float(L.split(",")[1]) / HARTREE_EV
        f = L.split()
        if len(f) > 6 and re.match(r"^\d+$", f[0]):
            lev.append(dict(i=int(f[0]), e=float(f[2]) / HARTREE_EV,
                            vnl=f[4], tj=int(f[5])))
    return E0, lev

def get_1P1_excitation_energy(lev):
    """Finds the 1s2p 1P1 excitation energy (in a.u.) relative to ground."""
    if not lev: return None
    # 2p orbital is vnl='201', 1P1 has 2J=2 (tj=2). 
    # It is the upper state of the two J=1 levels in 1s2p.
    c = sorted([l for l in lev if l["vnl"] == "201" and l["tj"] == 2],
               key=lambda l: l["e"])
    return c[1]["e"] if len(c) > 1 else None

def run_fac(tag, lam, ee):
    sf_path = os.path.join(WORKDIR, tag + ".sf")
    with open(sf_path, "w") as f:
        f.write(sf_input(tag, lam, ee))
        
    r = subprocess.run([SFAC, tag + ".sf"], cwd=WORKDIR, capture_output=True, text=True)
    if r.returncode != 0:
        return None, None
        
    E0, lev = parse_lev(os.path.join(WORKDIR, tag + ".lev"))
    exc = get_1P1_excitation_energy(lev)
    return E0, exc

# --------------------------------------------------------------------------
# MAIN EXECUTION & ERROR CALCULATION
# --------------------------------------------------------------------------
def main():
    os.makedirs(WORKDIR, exist_ok=True)
    print(f"Running FAC calculations for Kr34+ in {WORKDIR}/ ...\n")
    
    results = []
    errors_A = []
    errors_B = []
    
    print("=" * 115)
    print(f"{'lam_D':<7} | {'Ref A (b)':<10} | {'FAC A':<10} | {'Err_A (%)':<10} | {'Ref B (g)':<10} | {'FAC B':<10} | {'Err_B (%)':<10}")
    print("-" * 115)
    
    for i, lam in enumerate(LAMBDAS):
        lam_str = f"{lam:.1f}" if lam != 1.0e6 else "inf"
        
        ref_A = REF_EXC_A[i]
        ref_B = REF_EXC_B[i]
        
        # Model A run (ee=0)
        tag_A = f"Kr34_A_L{lam_str}".replace(".", "p")
        _, exc_A = run_fac(tag_A, lam, ee=0)
        
        # Model B run (ee=1)
        tag_B = f"Kr34_B_L{lam_str}".replace(".", "p")
        _, exc_B = run_fac(tag_B, lam, ee=1)
        
        errA_pct = None
        if exc_A and ref_A:
            errA = exc_A - ref_A
            errA_pct = 100 * errA / abs(ref_A)
            errors_A.append(errA_pct)
            
        errB_pct = None
        if exc_B and ref_B:
            errB = exc_B - ref_B
            errB_pct = 100 * errB / abs(ref_B)
            errors_B.append(errB_pct)
            
        def fmt(val): return f"{val:.4f}" if val is not None else "N/A"
        def fmt_e(val): return f"{val:+.5f}" if val is not None else "N/A"
        
        print(f"{lam_str:<7} | {fmt(ref_A):<10} | {fmt(exc_A):<10} | {fmt_e(errA_pct):<10} | {fmt(ref_B):<10} | {fmt(exc_B):<10} | {fmt_e(errB_pct):<10}")
        
        results.append({
            'lam': lam,
            'ref_A': ref_A, 'fac_A': exc_A, 
            'ref_B': ref_B, 'fac_B': exc_B
        })
        
    print("=" * 115)
    
    if errors_A:
        max_err_A = max(abs(e) for e in errors_A)
        rms_A = math.sqrt(sum(e**2 for e in errors_A)/len(errors_A))
        print(f"Model A (ee=0) -> Max |Error|: {max_err_A:.5f} %  |  RMS Error: {rms_A:.5f} %")
    if errors_B:
        max_err_B = max(abs(e) for e in errors_B)
        rms_B = math.sqrt(sum(e**2 for e in errors_B)/len(errors_B))
        print(f"Model B (ee=1) -> Max |Error|: {max_err_B:.5f} %  |  RMS Error: {rms_B:.5f} %")
        
    if HAS_MATPLOTLIB:
        generate_plots(results)

def generate_plots(results):
    # Filter out the infinity point for plotting against a standard linear axis
    finite_res = [r for r in results if r['lam'] != 1.0e6]
    
    lam_vals = [r['lam'] for r in finite_res]
    
    ref_A = [r['ref_A'] for r in finite_res]
    fac_A = [r['fac_A'] for r in finite_res]
    
    ref_B = [r['ref_B'] for r in finite_res]
    fac_B = [r['fac_B'] for r in finite_res]

    # Plot 1: Excitation Energy vs Debye Length
    plt.figure(figsize=(9, 6))
    plt.plot(lam_vals, ref_A, 'k--', label="Ref 'b' (Model A)")
    plt.plot(lam_vals, fac_A, 'b^', alpha=0.6, label='FAC Model A')
    plt.plot(lam_vals, ref_B, 'k-', label="Ref 'g' (Model B)")
    plt.plot(lam_vals, fac_B, 'ro', alpha=0.6, label='FAC Model B')
    
    plt.xlabel(r'Debye Length $\lambda_D$ ($a_0$)')
    plt.ylabel(r'Excitation Energy $\Delta E$ (a.u.)')
    plt.title(r'Kr$^{34+}$ $1s^2 \ ^1S_0 \rightarrow 1s2p \ ^1P_1$ Transition Energy')
    # Invert X axis so stronger screening (smaller lambda) is on the right
    plt.gca().invert_xaxis()
    plt.legend()
    plt.grid(True, ls=':', alpha=0.7)
    plt.savefig('chen2018_kr34_reproduction.png', dpi=300)
    plt.close()
    
    # Plot 2: E-E Screening Increment
    plt.figure(figsize=(9, 6))
    diff_ref = [r['ref_B'] - r['ref_A'] for r in finite_res]
    diff_fac = [r['fac_B'] - r['fac_A'] for r in finite_res]
    
    plt.plot(lam_vals, diff_ref, 'k-', label="Ref 'g' - 'b'")
    plt.plot(lam_vals, diff_fac, 'gD', label='FAC Model B - A')
    
    plt.xlabel(r'Debye Length $\lambda_D$ ($a_0$)')
    plt.ylabel('Energy Difference (a.u.)')
    plt.title('Kr$^{34+}$ Electron-Electron Screening Increment')
    plt.gca().invert_xaxis()
    plt.legend()
    plt.grid(True, ls=':', alpha=0.7)
    plt.savefig('chen2018_kr34_ee_increment.png', dpi=300)
    plt.close()
    print("\nSaved plots to: chen2018_kr34_reproduction.png, chen2018_kr34_ee_increment.png")

if __name__ == "__main__":
    main()
