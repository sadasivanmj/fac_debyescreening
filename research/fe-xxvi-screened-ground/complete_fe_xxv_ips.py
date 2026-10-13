# Complete the Fe XXV ionization potentials of Deprince et al. (2019), Table 4.
#
# Usage: python complete_fe_xxv_ips.py OLD_RUN NEW_RUN OUT_DIR
#   OLD_RUN : the validated run with the 36 good energies (results/run_20261013_minimal)
#   NEW_RUN : the rerun with the patched build (results/run_20261013_patched_rerun)
#
# 1. Consistency: every energy that was valid in OLD_RUN must be reproduced by NEW_RUN.
# 2. Fe XXVI: its energy must not depend on mu_ee (one electron).
# 3. IP(Fe XXV) = E0(Fe XXVI) - E0(Fe XXV), with the validated Fe XXV energies
#    from OLD_RUN and the Fe XXVI energy of the same mu_ne from NEW_RUN.

import csv
import json
import os
import sys

HARTREE_EV = 27.211386
CASES = {"A": (0.10, 0.00), "B": (0.10, 0.10), "C": (0.25, 0.00), "D": (0.25, 0.25)}
TABLE4 = {"A": 8766.28, "B": 8768.90, "C": 8661.14, "D": 8667.86}
IONS = ["xvii", "xviii", "xix", "xx", "xxi", "xxii", "xxiii", "xxiv", "xxv", "xxvi"]

old_run, new_run, out_dir = sys.argv[1:4]
os.makedirs(out_dir)


def manifest(run, ion, case):
    return json.load(open(os.path.join(run, "runs", "fe_" + ion, "case_" + case, "manifest.json")))


# ---- 1. old energies reproduced by the patched build ----
max_diff = 0.0
n_compared = 0
for case in CASES:
    for ion in IONS:
        old = manifest(old_run, ion, case)
        if old["status"] != "ok":
            continue
        new = manifest(new_run, ion, case)
        assert new["status"] == "ok", (ion, case, new["status"])
        max_diff = max(max_diff, abs(new["E0_eV"] - old["E0_eV"]))
        n_compared += 1
print(f"{n_compared} old energies reproduced by the patched build; max |dE| = {max_diff:.3e} eV")

# ---- 2. Fe XXVI independent of mu_ee ----
e26 = {}
for mu_ne in (0.10, 0.25):
    cases = [c for c, (m, _) in CASES.items() if m == mu_ne]
    vals = [manifest(new_run, "xxvi", c) for c in cases]
    for v in vals:
        assert v["status"] == "ok" and abs(v["mu_ne_used"] - mu_ne) < 1e-6 * mu_ne
    spread = max(v["E0_eV"] for v in vals) - min(v["E0_eV"] for v in vals)
    e26[mu_ne] = vals[0]["E0_eV"]
    print(f"Fe XXVI mu_ne={mu_ne}: E0 = {e26[mu_ne]:.5f} eV, spread over mu_ee = {spread:.3e} eV")

# ---- 3. the four IPs ----
rows = []
for case, (mu_ne, mu_ee) in CASES.items():
    m25 = manifest(old_run, "xxv", case)
    assert m25["status"] == "ok" and m25["mu_ne_requested"] == mu_ne and m25["mu_ee_requested"] == mu_ee
    e25 = m25["E0_eV"]
    d_hartree = (e26[mu_ne] - e25) / HARTREE_EV
    ip = d_hartree * HARTREE_EV
    rows.append({"case": case, "mu_ne": mu_ne, "mu_ee": mu_ee,
                 "E0_FeXXV_eV": e25, "E0_FeXXVI_eV": e26[mu_ne],
                 "E0_FeXXV_Ha": e25 / HARTREE_EV, "E0_FeXXVI_Ha": e26[mu_ne] / HARTREE_EV,
                 "IP_calc_eV": round(ip, 4), "IP_table4_eV": TABLE4[case],
                 "diff_eV": round(ip - TABLE4[case], 4),
                 "diff_percent": round(100 * (ip - TABLE4[case]) / TABLE4[case], 4)})
by = {r["case"]: r for r in rows}
for r in rows:
    partner = {"A": "B", "B": "A", "C": "D", "D": "C"}[r["case"]]
    lo, hi = sorted([r["case"], partner])
    r["ee_effect_calc_eV"] = round(by[hi]["IP_calc_eV"] - by[lo]["IP_calc_eV"], 4)
    r["ee_effect_table4_eV"] = round(TABLE4[hi] - TABLE4[lo], 4)
with open(os.path.join(out_dir, "fe_xxv_ip_table4.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r["case"], r["IP_calc_eV"], r["IP_table4_eV"], r["diff_eV"], r["diff_percent"],
          "ee effect", r["ee_effect_calc_eV"], "vs", r["ee_effect_table4_eV"])
