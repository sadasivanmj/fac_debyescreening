# Collect the 18 reference-configuration runs into two CSV tables and compute
# the ionization potentials IP(Fe X) = E0(Fe X+1) - E0(Fe X).
# Usage:  python3 collect_reference_results.py   (after run_reference_suite.sh)

import csv
import os

HARTREE_EV = 27.211386   # FAC writes E0 in eV; this converts to Hartree
IONS = ["fe_ii", "fe_iii", "fe_iv"]
MODELS = ["ref", "gref"]
MU_VALUES = ["0.0", "0.1", "0.25"]


def read_run(ion, model, mu):
    folder = f"reference_runs/{ion}/{model}_mu{mu}"
    run = {"ion_key": ion, "model": model, "mu": float(mu), "status": "missing"}
    if os.path.exists(folder + "/summary.txt"):
        for line in open(folder + "/summary.txt"):
            key, value = line.split(" = ", 1)
            run[key] = value.strip()
    log = open(folder + ".log").read() if os.path.exists(folder + ".log") else ""
    # FAC prints this when an orbital optimization did not converge
    if "Maximum iteration reached in OptimizeRadial" in log:
        run["status"] = "FAILED: SCF not converged"
    if "exit code 0" not in log and run["status"] == "ok":
        run["status"] = "FAILED: non-zero exit code"
    return run


runs = {}
with open("reference_energies.csv", "w", newline="") as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(["ion", "electrons", "model", "mu_requested", "mu_used",
                     "E0_eV", "E0_Hartree", "status"])
    for ion in IONS:
        for model in MODELS:
            for mu in MU_VALUES:
                run = read_run(ion, model, mu)
                runs[(ion, model, mu)] = run
                energy_ev = run.get("E0_eV")
                energy_ha = float(energy_ev) / HARTREE_EV if run["status"] == "ok" else ""
                writer.writerow([run.get("ion", ion), run.get("electrons", ""), model, mu,
                                 run.get("mu_used", ""), energy_ev if run["status"] == "ok" else "",
                                 energy_ha, run["status"]])

steps = [("fe_ii", "fe_iii", "Fe II -> Fe III"), ("fe_iii", "fe_iv", "Fe III -> Fe IV")]
with open("reference_ionization_potentials.csv", "w", newline="") as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(["ionization", "model", "mu", "IP_eV", "IP_Hartree", "dIP_eV", "status"])
    for lower, upper, label in steps:
        for model in MODELS:
            ip_at_zero = None
            for mu in MU_VALUES:
                a, b = runs[(lower, model, mu)], runs[(upper, model, mu)]
                if a["status"] != "ok" or b["status"] != "ok":
                    writer.writerow([label, model, mu, "", "", "", "not computed (missing energy)"])
                    continue
                ip = float(b["E0_eV"]) - float(a["E0_eV"])
                if mu == "0.0":
                    ip_at_zero = ip
                shift = round(ip - ip_at_zero, 4) if ip_at_zero is not None else ""
                writer.writerow([label, model, mu, round(ip, 4), ip / HARTREE_EV, shift, "ok"])
                print(f"{label:16s} {model:4s} mu={mu:4s} IP = {ip:9.4f} eV")
