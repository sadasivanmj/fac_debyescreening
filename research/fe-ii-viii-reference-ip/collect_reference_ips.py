# Collect the Fe II - Fe IX ground-state energies, form the ionization potentials
#   I_q(mu) = [E0(Fe q+1) - E0(Fe q)] / 27.211386 Ha * 27.211386 eV/Ha
# (FAC writes E0 in eV, so the IP is the difference of the E0 values) and compare
# them with the paper's Table 2 and with NIST.
#
# Usage:  python collect_reference_ips.py results/RUN_ID
# Writes energies.csv and ionization_potentials.csv into that folder.  An IP is
# only computed from two runs that are "ok" (exit code 0, FAC checks passed, no
# "Maximum iteration reached" in the log) with identical settings.

import csv
import json
import math
import os
import sys

HARTREE_EV = 27.211386
IONS = ["Fe II", "Fe III", "Fe IV", "Fe V", "Fe VI", "Fe VII", "Fe VIII", "Fe IX"]
MODELS = ["gref", "ref"]
MUS = ["0.0", "0.1", "0.25"]

# Comparison values only (not used in the calculation).
# Deprince et al. A&A 643, A57 (2020), Table 2: IP of the named ion at mu = 0, 0.1, 0.25.
PAPER = {"Fe II": (18.89, 14.20, 7.78), "Fe III": (32.93, 25.60, 16.42),
         "Fe IV": (55.37, 44.65, 29.64), "Fe V": (74.92, 61.56, 42.83),
         "Fe VI": (98.20, 82.20, 59.79), "Fe VII": (123.60, 104.96, 78.86),
         "Fe VIII": (151.14, 129.85, 99.98)}
NIST = {"Fe II": 16.19920, "Fe III": 30.651, "Fe IV": 54.91, "Fe V": 75.0,
        "Fe VI": 98.985, "Fe VII": 124.98, "Fe VIII": 151.060}

run_dir = sys.argv[1]


def read_run(ion, model, mu):
    key = ion.lower().replace(" ", "_")
    folder = os.path.join(run_dir, "runs", key, f"{model}_mu{mu}")
    log = os.path.join(run_dir, "logs", f"{key}_{model}_mu{mu}.log")
    row = {"ion": ion, "model": model, "mu": mu, "status": "missing", "E0_eV": ""}
    if not os.path.exists(os.path.join(folder, "manifest.json")):
        return row, None
    m = json.load(open(os.path.join(folder, "manifest.json")))
    text = open(log).read() if os.path.exists(log) else ""
    status = m["status"]
    if "Maximum iteration reached in OptimizeRadial" in text:
        status = "FAILED: SCF not converged"
    e = m["E0_eV"]
    if status == "ok" and not (e is not None and math.isfinite(e)):
        status = "FAILED: no finite energy"
    if status == "ok" and abs(m["mu_used"] - float(mu)) > 1e-6:
        status = "FAILED: mu mismatch"
    row.update({"status": status, "E0_eV": e, "E0_Hartree": e / HARTREE_EV if e is not None else "",
                "electrons": m["electrons"], "n_configurations": len(m["configurations"]),
                "n_levels": m["n_levels"], "mu_used": m["mu_used"], "hxs_iter0": m["hxs_iter0"],
                "debye_mode": m["debye_mode"], "ee_screen": m["ee_screen"], "commit": m["commit"]})
    return row, m


runs, manifests, rows = {}, {}, []
for ion in IONS:
    for model in MODELS:
        for mu in MUS:
            row, m = read_run(ion, model, mu)
            runs[(ion, model, mu)] = row
            manifests[(ion, model, mu)] = m
            rows.append(row)
fields = ["ion", "model", "mu", "status", "electrons", "n_configurations", "n_levels",
          "mu_used", "hxs_iter0", "debye_mode", "ee_screen", "E0_eV", "E0_Hartree", "commit"]
with open(os.path.join(run_dir, "energies.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

out = []
for k, ion in enumerate(IONS[:-1]):
    upper = IONS[k + 1]
    for model in MODELS:
        ip0 = None
        for i, mu in enumerate(MUS):
            a, b = runs[(ion, model, mu)], runs[(upper, model, mu)]
            row = {"ionization": f"{ion} -> {upper}", "ion": ion, "model": model, "mu": mu,
                   "IP_calc_eV": "", "IP_paper_eV": PAPER[ion][i],
                   "IP_NIST_eV": NIST[ion] if mu == "0.0" else "",
                   "err_vs_paper_eV": "", "err_vs_NIST_eV": "",
                   "shift_calc_eV": "", "shift_paper_eV": round(PAPER[ion][i] - PAPER[ion][0], 4),
                   "status": "ok"}
            if a["status"] != "ok" or b["status"] != "ok":
                row["status"] = "not computed (missing or failed energy)"
            else:
                ip = (b["E0_eV"] / HARTREE_EV - a["E0_eV"] / HARTREE_EV) * HARTREE_EV
                if mu == "0.0":
                    ip0 = ip
                row["IP_calc_eV"] = round(ip, 4)
                row["err_vs_paper_eV"] = round(ip - PAPER[ion][i], 4)
                if mu == "0.0":
                    row["err_vs_NIST_eV"] = round(ip - NIST[ion], 4)
                if ip0 is not None:
                    row["shift_calc_eV"] = round(ip - ip0, 4)
            out.append(row)
fields = ["ionization", "model", "mu", "IP_calc_eV", "IP_paper_eV", "IP_NIST_eV",
          "err_vs_paper_eV", "err_vs_NIST_eV", "shift_calc_eV", "shift_paper_eV", "status"]
with open(os.path.join(run_dir, "ionization_potentials.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(out)

n_ok = sum(1 for r in rows if r["status"] == "ok")
print(f"{n_ok}/{len(rows)} energies usable, {sum(1 for r in out if r['IP_calc_eV'] != '')}/{len(out)} IPs computed")
