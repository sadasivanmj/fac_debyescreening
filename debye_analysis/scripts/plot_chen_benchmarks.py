#!/usr/bin/env python3
"""
plot_chen_benchmarks.py
Reads existing FAC data from chen2018_results.json and plots reproduction curves 
and errors against Chen et al. (2018) Table II.
"""

import os
import json
import csv
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ENVIRONMENT & CONSTANTS
# =============================================================================
RESULTS_FILE = "chen2018_results.json"
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

TARGET_LEVELS = [
    "1s2s 3S1",
    "1s2p 3P1",
    "1s2p 1P1",
    "1s3p 3P1",
    "1s3p 1P1"
]

LAMBDAS = [float('inf'), 20.0, 15.0, 10.0, 7.5, 5.0]

# =============================================================================
# REFERENCE DATA (CHEN ET AL. 2018, TABLE II, C4+)
# Transcribed directly from your chen2018_reproduce.py file
# =============================================================================
CHEN_REF = {
    "en_only": {
        "1s2s 3S1": {np.inf: 11.1297, 20.0: 11.0882, 15.0: 11.0877, 10.0: 11.0862, 7.5: 11.0842, 5.0: 11.0056},
        "1s2p 3P1": {np.inf: 11.3987, 20.0: 11.3472, 15.0: 11.3466, 10.0: 11.3450, 7.5: 11.3426, 5.0: 11.2426},
        "1s2p 1P1": {np.inf: 11.5051, 20.0: 11.4844, 15.0: 11.4838, 10.0: 11.4820, 7.5: 11.4794, 5.0: 11.3777},
        "1s3p 3P1": {np.inf: 13.0612, 20.0: 13.0139, 15.0: 13.0125, 10.0: 13.0085, 7.5: 13.0033, 5.0: 12.9629},
        "1s3p 1P1": {np.inf: 13.1004, 20.0: 13.0529, 15.0: 13.0512, 10.0: 13.0466, 7.5: 13.0402, 5.0: 12.9936},
    },
    "en_ee": {
        # Chen only reported the 'g' model for the 1s2p 1P1 level for C4+
        "1s2s 3S1": {np.inf: np.nan, 20.0: np.nan, 15.0: np.nan, 10.0: np.nan, 7.5: np.nan, 5.0: np.nan},
        "1s2p 3P1": {np.inf: np.nan, 20.0: np.nan, 15.0: np.nan, 10.0: np.nan, 7.5: np.nan, 5.0: np.nan},
        "1s2p 1P1": {np.inf: np.nan, 20.0: 11.4771, 15.0: 11.4762, 10.0: 11.4688, 7.5: 11.4687, 5.0: 11.3685},
        "1s3p 3P1": {np.inf: np.nan, 20.0: np.nan, 15.0: np.nan, 10.0: np.nan, 7.5: np.nan, 5.0: np.nan},
        "1s3p 1P1": {np.inf: np.nan, 20.0: np.nan, 15.0: np.nan, 10.0: np.nan, 7.5: np.nan, 5.0: np.nan},
    }
}

# =============================================================================
# DATA PARSING
# =============================================================================
def load_fac_data():
    with open(RESULTS_FILE, 'r') as f:
        data = json.load(f)["results"]
    
    results = []
    
    for model_key, ee_flag in [("en_only", 0), ("en_ee", 1)]:
        for lam in LAMBDAS:
            # Construct the JSON key (e.g., C_ee0_L20)
            lam_str = "1e06" if lam == float('inf') else str(lam).replace(".", "p")
            # Remove trailing 'p0' for integers like 20.0 -> 20
            if lam_str.endswith("p0"): lam_str = lam_str[:-2]
                
            tag = f"C_ee{ee_flag}_L{lam_str}"
            
            if tag not in data:
                print(f"Warning: {tag} not found in JSON.")
                continue
                
            run_data = data[tag]["lev"]
            
            for level in TARGET_LEVELS:
                fac_e = run_data.get(level, np.nan)
                chen_e = CHEN_REF[model_key][level].get(lam, np.nan)
                
                if not np.isnan(fac_e) and not np.isnan(chen_e):
                    err = fac_e - chen_e
                    abs_err = abs(err)
                    pct_err = 100 * abs_err / abs(chen_e)
                    signed_pct = 100 * err / chen_e
                else:
                    err = abs_err = pct_err = signed_pct = np.nan
                    
                results.append({
                    "model": model_key, "level": level, "lambda_D_a0": lam,
                    "fac_energy": fac_e, "chen_energy": chen_e,
                    "signed_error": err, "absolute_error": abs_err,
                    "signed_percent_error": signed_pct, "absolute_percent_error": pct_err,
                    "status": "PASS" if not np.isnan(fac_e) else "FAIL"
                })
                
    return results

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================
def plot_reproduction(model, data):
    plt.figure(figsize=(8, 6))
    finite_lams = [l for l in LAMBDAS if l != float('inf')]
    
    plotted_any = False
    for level in TARGET_LEVELS:
        x_fin, y_fac, y_chen = [], [], []
        inf_fac, inf_chen = None, None
        
        for row in data:
            if row["model"] == model and row["level"] == level and row["status"] == "PASS":
                if row["lambda_D_a0"] == float('inf'):
                    inf_fac = row["fac_energy"]
                    inf_chen = row["chen_energy"]
                else:
                    x_fin.append(row["lambda_D_a0"])
                    y_fac.append(row["fac_energy"])
                    # Only plot Chen if it's not NaN
                    if not np.isnan(row["chen_energy"]):
                        y_chen.append(row["chen_energy"])
                    
        if x_fin:
            p = plt.plot(x_fin, y_fac, marker='o', label=f"FAC {level}", linestyle='-')
            color = p[0].get_color()
            if len(y_chen) == len(x_fin):
                plt.plot(x_fin, y_chen, marker='s', label=f"Chen {level}", linestyle='--', color=color, alpha=0.7)
            if inf_fac is not None:
                plt.axhline(inf_fac, color=color, linestyle=':', alpha=0.5)
            plotted_any = True
            
    if not plotted_any:
        plt.close()
        return

    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.xlabel(r"Debye Length $\lambda_D$ ($a_0$)")
    plt.ylabel("Excitation Energy (a.u.)")
    plt.title(f"Reproduction Curves: {model}\n(Dotted lines = Unscreened limits)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, f"chen_reproduction_{model}.pdf"), dpi=300)
    plt.savefig(os.path.join(FIGDIR, f"chen_reproduction_{model}.png"), dpi=300)
    plt.close()

def plot_error(model, data):
    plt.figure(figsize=(8, 6))
    
    plotted_any = False
    for level in TARGET_LEVELS:
        x_fin, y_err = [], []
        for row in data:
            if row["model"] == model and row["level"] == level and row["status"] == "PASS":
                if row["lambda_D_a0"] != float('inf') and not np.isnan(row["absolute_percent_error"]):
                    x_fin.append(row["lambda_D_a0"])
                    y_err.append(row["absolute_percent_error"])
                    
        if x_fin:
            plt.plot(x_fin, y_err, marker='o', label=level)
            plotted_any = True
            
    if not plotted_any:
        plt.close()
        return

    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.xlabel(r"Debye Length $\lambda_D$ ($a_0$)")
    plt.ylabel("Absolute Percentage Error (%) [Lower is Better]")
    plt.title(f"Numerical Error vs Chen et al.: {model}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, f"chen_error_{model}.pdf"), dpi=300)
    plt.savefig(os.path.join(FIGDIR, f"chen_error_{model}.png"), dpi=300)
    plt.close()

# =============================================================================
# MAIN EXECUTION
# =============================================================================
def main():
    print(f"Reading data from {RESULTS_FILE}...")
    results = load_fac_data()
    
    # Save CSV
    csv_file = "chen_reproduction.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        for r in results:
            writer.writerow(r)
            
    # Generate Plots
    print("Generating Plots...")
    plot_reproduction("en_only", results)
    plot_error("en_only", results)
    plot_reproduction("en_ee", results)
    plot_error("en_ee", results)
    
    print("Done! Check the 'figures/' folder and 'chen_reproduction.csv'.")

if __name__ == "__main__":
    main()
