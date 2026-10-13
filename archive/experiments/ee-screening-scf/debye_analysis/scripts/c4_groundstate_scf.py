#!/usr/bin/env python3
"""
chen_groundstate.py
Calculates the ground-state energy of He-like C4+ in a Debye-Huckel plasma 
using a modified Flexible Atomic Code (FAC)[cite: 1], with full CI up to n=6[cite: 2].
"""

import os
import math
import subprocess
import csv
import pathlib
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# --- CONSTANTS & CONVENTIONS[cite: 1] ---
NPS = 1.0e-2
RBOHR = 0.52917721067
HARTREE_EV = 27.211386018

def tps_for(lam):
    """eV temperature that makes dps == lam (bohr), given nps in A^-3 and ups=0[cite: 1]."""
    return lam * lam * 4.0 * math.pi * (NPS * RBOHR**3) * HARTREE_EV

# --- CHEN ET AL. (2018) TABLE I FAC REFERENCE DATA[cite: 2] ---
CHEN_TABLE_I_FAC = {
    "inf": -32.4771,
    100.0: -32.3768,
    20.0:  -31.9081,
    12.5:  -31.5601,
    10.0:  -31.3298,
    5.0:   -30.1973,
    3.333: -29.0950,
    2.5:   -28.0242,
    2.0:   -26.9818
}

def generate_sf_input(tag, lam_a0, max_n):
    """Generates the .sf script with dynamic CI up to max_n[cite: 2]."""
    tps = tps_for(lam_a0)
    
    lines = [
        "SetAtom('C')",
        "SetOption('orbital:debye_mode', 1)",
        "SetOption('radial:ee_screen', 1)",
        "Config('g1', '1s2')"
    ]
    
    groups = ["'g1'"]
    
    # Dynamically build the configuration basis up to max_n
    if max_n >= 2:
        for n1 in range(1, max_n + 1):
            cfgs = []
            for n2 in range(n1, max_n + 1):
                if n1 == 1 and n2 == 1:
                    continue # already handled by g1
                cfgs.append(f"'{n1}*1 {n2}*1'" if n1 != n2 else f"'{n1}*2'")
            if cfgs:
                group_name = f"g_n{n1}"
                lines.append(f"Config('{group_name}', {', '.join(cfgs)})")
                groups.append(f"'{group_name}'")

    lines.extend([
        f"PlasmaScreen(4.0, {NPS}, {tps}, 1, 0.0, 1)",
        "ConfigEnergy(0)",
        "OptimizeRadial(['g1'])",
        "ConfigEnergy(1)",
        f"Structure('{tag}.lev.b', [{', '.join(groups)}])",
        f"PrintTable('{tag}.lev.b', '{tag}.lev', 1)"
    ])
    
    return "\n".join(lines) + "\n"

def parse_e0(lev_file):
    """Parses the total ground state energy E0 from the .lev file[cite: 1]."""
    if not os.path.exists(lev_file):
        return None
    with open(lev_file, 'r') as f:
        for line in f:
            if line.startswith("E0"):
                e0_ev = float(line.split(',')[1])
                return e0_ev / HARTREE_EV
    return None

def main():
    sfac_bin = os.environ.get("SFAC", "sfac")
    work_dir = pathlib.Path("c4_ground_results")
    work_dir.mkdir(exist_ok=True)
    
    results = []
    
    print(f"Using FAC binary: {sfac_bin}")
    print(f"Output directory: {work_dir}/\n")
    print("-" * 85)
    print(f"{'lambda_D (a0)':<13} | {'Calculated E0':<15} | {'Chen E0 (Lit)':<15} | {'Abs Diff':<10} | {'Rel Err %':<10} | {'Max n'}")
    print("-" * 85)

    for lam_label, chen_e0 in CHEN_TABLE_I_FAC.items():
        lam_num = 1.0e6 if lam_label == "inf" else float(lam_label)
        tag_base = f"C4_ground_L{lam_label}".replace(".", "p")
        
        e0_calc = None
        best_n = None
        
        # Fallback loop: start at n=6 (full CI), step down if FAC fails[cite: 1, 2]
        for current_n in range(6, 0, -1):
            tag = f"{tag_base}_n{current_n}"
            sf_path = work_dir / f"{tag}.sf"
            sf_path.write_text(generate_sf_input(tag, lam_num, current_n))
            
            r = subprocess.run([sfac_bin, sf_path.name], cwd=work_dir, capture_output=True, text=True, timeout=120)
            
            if r.returncode == 0:
                lev_path = work_dir / f"{tag}.lev"
                val = parse_e0(lev_path)
                if val is not None:
                    e0_calc = val
                    best_n = current_n
                    break
        
        if e0_calc is None:
            print(f"{str(lam_label):<13} | {'FAILED':<15} | {chen_e0:<15.4f} | {'-':<10} | {'-':<10} | -")
            continue
            
        abs_err = e0_calc - chen_e0
        rel_err = 100.0 * abs_err / abs(chen_e0)
        
        print(f"{str(lam_label):<13} | {e0_calc:<15.4f} | {chen_e0:<15.4f} | {abs_err:<10.4f} | {rel_err:<10.4f} | n={best_n}")
            
        results.append({
            "lambda_D_label": lam_label,
            "lambda_D_a0": lam_num,
            "E0_fac_au": e0_calc,
            "chen_fac_au": chen_e0,
            "absolute_error_au": abs_err,
            "relative_error_percent": rel_err,
            "max_n_used": best_n
        })

    print("-" * 85)

    if results:
        csv_path = work_dir / "c4_ground_comparison.csv"
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        
        max_abs = max(abs(r["absolute_error_au"]) for r in results)
        max_rel = max(abs(r["relative_error_percent"]) for r in results)
        print(f"\nMaximum Absolute Error: {max_abs:.6f} a.u.")
        print(f"Maximum Relative Error: {max_rel:.6f} %")

        if HAS_MATPLOTLIB:
            plot_results(results, work_dir)

def plot_results(results, out_dir):
    finite_res = [r for r in results if r["lambda_D_label"] != "inf"]
    inf_res = next((r for r in results if r["lambda_D_label"] == "inf"), None)
    
    lams = [r["lambda_D_a0"] for r in finite_res]
    e0_calc = [r["E0_fac_au"] for r in finite_res]
    e0_chen = [r["chen_fac_au"] for r in finite_res]
    
    plt.figure(figsize=(8, 6))
    plt.plot(lams, e0_calc, 'b-', label="Calculated (Modified FAC with CI)")
    plt.plot(lams, e0_chen, 'ro', label="Chen et al. (Table I)")
    
    if inf_res:
        plt.axhline(y=inf_res["E0_fac_au"], color='b', linestyle='--', label=r"Calc $\lambda_D = \infty$ limit")
        plt.axhline(y=inf_res["chen_fac_au"], color='r', linestyle=':', label=r"Chen $\lambda_D = \infty$ limit")

    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.xlabel(r"Debye Length $\lambda_D$ ($a_0$) [Log Scale]")
    plt.ylabel(r"Ground-State Energy $E_0$ (a.u.)")
    plt.title(r"He-like C$^{4+}$ Ground State vs Debye Length")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_dir / "c4_ground_energy_ci.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
