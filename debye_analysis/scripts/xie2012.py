#!/usr/bin/env python3
"""
xie2012_table2_reproduce.py
===========================
Reproduces Table 2 of Xie et al., Eur. Phys. J. D (2012) 66: 125.
Calculates ground state (1s^2 1S0) and excited state (1s2p 1P1) energies 
for He-like C4+ in Debye-Huckel plasmas.
"""

import os
import re
import sys
import math
import subprocess
try:
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# --------------------------------------------------------------------------
# FAC CONVENTIONS & CONSTANTS
# --------------------------------------------------------------------------
SFAC       = os.environ.get("SFAC", "sfac")
WORKDIR    = "xie2012_run"
NPS        = 1.0e-2
RBOHR      = 0.52917721067
HARTREE_EV = 27.211386018

def tps_for(lam):
    """eV temperature that makes dps == lam (bohr), given nps in A^-3 and ups=0."""
    return lam * lam * 4.0 * math.pi * (NPS * RBOHR**3) * HARTREE_EV

# --------------------------------------------------------------------------
# XIE ET AL (2012) TABLE 2 DATA (He-like C4+)
# Quantities are tabulated as positive B = -E (binding energy in a.u.)
# --------------------------------------------------------------------------
MU_VALUES = [
    0, 0.01, 0.02, 0.05, 0.1, 0.125, 0.167, 0.2, 0.25, 0.3, 0.333, 
    0.4, 0.5, 0.6, 0.667, 0.7, 0.75, 0.8, 0.87, 0.9, 0.952, 1.0, 1.073
]

REF_GROUND_A = [
    32.4176, 32.2978, 32.1783, 31.8215, 31.2324, 30.9412, 30.4563, 
    30.0790, 29.5137, 28.9557, 28.5877, 27.8614, 26.7963, 25.7576, 
    25.0767, 24.7457, 24.2492, 23.7594, 23.0841, 22.7986, 22.3089, 
    21.8629, None
]
REF_GROUND_B = [
    32.4176, 32.3079, 32.1982, 31.8711, 31.3306, 31.0633, 30.6181, 
    30.2717, 29.7523, 29.2394, 28.9010, 28.2328, 27.2522, 26.2948, 
    25.6666, 25.3611, 24.9029, 24.4505, 23.8264, 23.5621, 23.1089, 
    22.6959, 22.0774
]

REF_EXCITED_A = [
    21.1023, 20.9827, 20.8638, 20.5116, 19.9388, 19.6590, 19.1984, 
    18.8446, 18.3219, 17.8149, 17.4853, 16.8470, 15.9389, 15.0883, 
    14.5511, 14.2963, 13.9231, 13.5661, None, None, None, 
    None, None
]
REF_EXCITED_B = [
    21.1023, 20.9927, 20.8836, 20.5604, 20.0339, 19.7763, 19.3518, 
    19.0253, 18.5421, 18.0724, 17.7665, 17.1717, 16.3227, 15.5208, 
    15.0098, 14.7657, 14.4054, 14.0566, 13.5876, 13.3935, 13.0673, 
    12.7781, 12.3489
]


# --------------------------------------------------------------------------
# FAC DRIVING FUNCTIONS
# --------------------------------------------------------------------------
def sf_input(tag, mu, ee, max_n):
    """Generates the FAC script dynamically based on max_n."""
    lam = 1.0e6 if mu == 0 else 1.0/mu
    tps = tps_for(lam)
    
    lines = [
        "SetAtom('C')",
        "SetOption('orbital:debye_mode', 1)",
        f"SetOption('radial:ee_screen', {ee})",
        "Config('g1', '1s2')",
        "Config('g2', '1s1 2*1')"
    ]
    groups = ["'g1'", "'g2'"]
    
    if max_n >= 3:
        lines.append("Config('g3', '1s1 3*1')")
        groups.append("'g3'")
    if max_n >= 4:
        lines.append("Config('g4', '1s1 4*1')")
        groups.append("'g4'")
    if max_n >= 5:
        lines.append("Config('g5', '1s1 5*1')")
        groups.append("'g5'")
        
    # Scale correlation configurations based on available orbitals
    if max_n >= 3:
        lines.append("Config('g_corr', '2*2', '2*1 3*1', '3*2')")
    else:
        lines.append("Config('g_corr', '2*2')")
    groups.append("'g_corr'")
    
    groups_str = ", ".join(groups)
    
    lines.extend([
        f"PlasmaScreen(4.0, {NPS}, {tps}, 1, 0.0, 1)",
        "ConfigEnergy(0)",
        "OptimizeRadial(['g1'])",
        "ConfigEnergy(1)",
        f"Structure('{tag}.lev.b', [{groups_str}])",
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
    """Finds the 1s2p 1P1 excitation energy relative to ground."""
    if not lev: return None
    c = sorted([l for l in lev if l["vnl"] == "201" and l["tj"] == 2],
               key=lambda l: l["e"])
    return c[1]["e"] if len(c) > 1 else None

def run_fac(tag, mu, ee):
    # Step down max_n if highly excited states leave the bound spectrum
    for max_n in [5, 4, 3, 2]:
        current_tag = f"{tag}_n{max_n}"
        sf_path = os.path.join(WORKDIR, current_tag + ".sf")
        with open(sf_path, "w") as f:
            f.write(sf_input(current_tag, mu, ee, max_n))
            
        r = subprocess.run([SFAC, current_tag + ".sf"], cwd=WORKDIR, capture_output=True, text=True)
        if r.returncode == 0:
            E0, lev = parse_lev(os.path.join(WORKDIR, current_tag + ".lev"))
            exc = get_1P1_excitation_energy(lev)
            if exc is not None:
                return E0, exc, max_n
    return None, None, None

# --------------------------------------------------------------------------
# MAIN EXECUTION & ERROR CALCULATION
# --------------------------------------------------------------------------
def main():
    os.makedirs(WORKDIR, exist_ok=True)
    print(f"Running FAC calculations with dynamic n_max fallback in {WORKDIR}/ ...\n")
    
    results = []
    errors_A = []
    errors_B = []
    
    print("=" * 125)
    print(f"{'mu':<6} | {'lam_D':<7} | {'Ref dE_A':<10} | {'FAC dE_A':<10} | {'Err_A (%)':<10} | {'Ref dE_B':<10} | {'FAC dE_B':<10} | {'Err_B (%)':<10} | {'Max n'}")
    print("-" * 125)
    
    for i, mu in enumerate(MU_VALUES):
        lam = 1.0e6 if mu == 0 else 1.0/mu
        lam_str = f"{lam:.2f}" if mu != 0 else "inf"
        
        ref_dE_A = REF_GROUND_A[i] - REF_EXCITED_A[i] if (REF_GROUND_A[i] and REF_EXCITED_A[i]) else None
        ref_dE_B = REF_GROUND_B[i] - REF_EXCITED_B[i] if (REF_GROUND_B[i] and REF_EXCITED_B[i]) else None
        
        tag_A = f"C4_A_mu{mu}".replace(".", "p")
        E0_A, exc_A, n_A = run_fac(tag_A, mu, ee=0)
        
        tag_B = f"C4_B_mu{mu}".replace(".", "p")
        E0_B, exc_B, n_B = run_fac(tag_B, mu, ee=1)
        
        errA_pct = None
        if exc_A and ref_dE_A:
            errA = exc_A - ref_dE_A
            errA_pct = 100 * errA / abs(ref_dE_A)
            errors_A.append(errA_pct)
            
        errB_pct = None
        if exc_B and ref_dE_B:
            errB = exc_B - ref_dE_B
            errB_pct = 100 * errB / abs(ref_dE_B)
            errors_B.append(errB_pct)
            
        def fmt(val): return f"{val:.5f}" if val is not None else "N/A"
        def fmt_e(val): return f"{val:+.3f}" if val is not None else "N/A"
        
        # Display the lowest n_max successfully used across both models
        used_n = min([n for n in (n_A, n_B) if n is not None], default="-")
        
        print(f"{mu:<6} | {lam_str:<7} | {fmt(ref_dE_A):<10} | {fmt(exc_A):<10} | {fmt_e(errA_pct):<10} | {fmt(ref_dE_B):<10} | {fmt(exc_B):<10} | {fmt_e(errB_pct):<10} | {used_n}")
        
        results.append({
            'mu': mu,
            'lam': lam,
            'ref_A': ref_dE_A, 'fac_A': exc_A, 
            'ref_B': ref_dE_B, 'fac_B': exc_B
        })
        
    print("=" * 125)
    
    if errors_A:
        max_err_A = max(abs(e) for e in errors_A)
        rms_A = math.sqrt(sum(e**2 for e in errors_A)/len(errors_A))
        print(f"Model A (ee=0) -> Max |Error|: {max_err_A:.3f} %  |  RMS Error: {rms_A:.3f} %")
    if errors_B:
        max_err_B = max(abs(e) for e in errors_B)
        rms_B = math.sqrt(sum(e**2 for e in errors_B)/len(errors_B))
        print(f"Model B (ee=1) -> Max |Error|: {max_err_B:.3f} %  |  RMS Error: {rms_B:.3f} %")
        
    if HAS_MATPLOTLIB:
        generate_plots(results)

def generate_plots(results):
    # Safe extraction filtering out N/A entries
    mu_A = [r['mu'] for r in results if r['ref_A'] and r['fac_A']]
    ref_A = [r['ref_A'] for r in results if r['ref_A'] and r['fac_A']]
    fac_A = [r['fac_A'] for r in results if r['ref_A'] and r['fac_A']]
    
    mu_B = [r['mu'] for r in results if r['ref_B'] and r['fac_B']]
    ref_B = [r['ref_B'] for r in results if r['ref_B'] and r['fac_B']]
    fac_B = [r['fac_B'] for r in results if r['ref_B'] and r['fac_B']]

    # Plot 1: Excitation Energy vs Mu
    plt.figure(figsize=(9, 6))
    if mu_A:
        plt.plot(mu_A, ref_A, 'k--', label='Ref Model A')
        plt.plot(mu_A, fac_A, 'b^', alpha=0.6, label='FAC Model A')
    if mu_B:
        plt.plot(mu_B, ref_B, 'k-', label='Ref Model B')
        plt.plot(mu_B, fac_B, 'ro', alpha=0.6, label='FAC Model B')
    plt.xlabel(r'Screening parameter $\mu$ ($a_0^{-1}$)')
    plt.ylabel(r'Excitation Energy $\Delta E$ (a.u.)')
    plt.title(r'C$^{4+}$ $1s^2 \ ^1S_0 \rightarrow 1s2p \ ^1P_1$ Transition Energy')
    plt.legend()
    plt.grid(True, ls=':', alpha=0.7)
    plt.savefig('xie2012_fig3_reproduction.png', dpi=300)
    plt.close()
    
    # Plot 2: E-E Screening Increment (Model B - Model A)
    plt.figure(figsize=(9, 6))
    
    # Strictly matched x and y sets to avoid Dimension errors
    mu_ref_diff = [r['mu'] for r in results if r['ref_B'] is not None and r['ref_A'] is not None]
    diff_ref = [r['ref_B'] - r['ref_A'] for r in results if r['ref_B'] is not None and r['ref_A'] is not None]
    
    mu_fac_diff = [r['mu'] for r in results if r['fac_B'] is not None and r['fac_A'] is not None]
    diff_fac = [r['fac_B'] - r['fac_A'] for r in results if r['fac_B'] is not None and r['fac_A'] is not None]
    
    if mu_ref_diff:
        plt.plot(mu_ref_diff, diff_ref, 'k-', label='Reference $\Delta E_B - \Delta E_A$')
    if mu_fac_diff:
        plt.plot(mu_fac_diff, diff_fac, 'gD', label='FAC $\Delta E_B - \Delta E_A$')
        
    plt.xlabel(r'Screening parameter $\mu$ ($a_0^{-1}$)')
    plt.ylabel('Energy Difference (a.u.)')
    plt.title('Effect of Electron-Electron Screening Increment')
    plt.legend()
    plt.grid(True, ls=':', alpha=0.7)
    plt.savefig('xie2012_ee_increment.png', dpi=300)
    plt.close()
    
    # Plot 3: Percentage Error vs Mu
    plt.figure(figsize=(9, 6))
    if mu_A:
        err_A = [100 * (f - r)/r for f, r in zip(fac_A, ref_A)]
        plt.plot(mu_A, err_A, 'b-^', label='Model A % Error')
    if mu_B:
        err_B = [100 * (f - r)/r for f, r in zip(fac_B, ref_B)]
        plt.plot(mu_B, err_B, 'r-o', label='Model B % Error')
        
    plt.xlabel(r'Screening parameter $\mu$ ($a_0^{-1}$)')
    plt.ylabel('Percentage Error (%)')
    plt.title('FAC Numerical Truncation Error relative to Xie 2012')
    plt.axhline(0, color='k', linewidth=1)
    plt.legend()
    plt.grid(True, ls=':', alpha=0.7)
    plt.savefig('xie2012_error_analysis.png', dpi=300)
    plt.close()
    
    print("\nSaved plots to: xie2012_fig3_reproduction.png, xie2012_ee_increment.png, xie2012_error_analysis.png")

if __name__ == "__main__":
    main()
