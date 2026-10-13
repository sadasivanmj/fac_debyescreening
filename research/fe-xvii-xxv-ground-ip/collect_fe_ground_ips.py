# Collect the ground-state energies of Fe XVII-XXVI, form the ionization
# potentials IP(q) = E0(q+1) - E0(q) and compare them with Table 4 of
# Deprince et al., A&A 626, A83 (2019).
#
# Usage:  python collect_fe_ground_ips.py results/RUN_ID
# Writes ground_energies.csv, ip_differences.csv and table4_ip_comparison.csv
# into that folder.  An IP is only computed when both energies come from runs
# with status "ok" and identical settings.

import csv
import json
import math
import os
import sys

HARTREE_EV = 27.211386  # as specified for the energy differences

# Ion stages in order; each IP uses the ion itself and the next ion up.
IONS = ["XVII", "XVIII", "XIX", "XX", "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI"]
CASES = {"A": (0.10, 0.00), "B": (0.10, 0.10), "C": (0.25, 0.00), "D": (0.25, 0.25)}

# Published Table 4 values (eV), used only for the comparison.
TABLE4 = {
    "XVII":   {"A": 1190.72, "B": 1214.75, "C": 1089.22, "D": 1147.80},
    "XVIII":  {"A": 1287.21, "B": 1308.56, "C": 1185.55, "D": 1237.68},
    "XIX":    {"A": 1389.10, "B": 1407.92, "C": 1287.41, "D": 1333.09},
    "XX":     {"A": 1503.54, "B": 1519.58, "C": 1401.59, "D": 1440.80},
    "XXI":    {"A": 1619.15, "B": 1632.53, "C": 1517.02, "D": 1549.75},
    "XXII":   {"A": 1727.82, "B": 1738.54, "C": 1625.57, "D": 1651.80},
    "XXIII":  {"A": 1880.58, "B": 1888.61, "C": 1778.76, "D": 1798.43},
    "XXIV":   {"A": 1974.41, "B": 1979.79, "C": 1872.58, "D": 1885.77},
    "XXV":    {"A": 8766.28, "B": 8768.90, "C": 8661.14, "D": 8667.86},
}

run_dir = sys.argv[1]

# ---- read the manifests ----
energies = {}   # (ion, case) -> manifest dict, only for usable runs
rows = []
for case, (mu_ne, mu_ee) in CASES.items():
    for ion in IONS:
        path = os.path.join(run_dir, "runs", "fe_" + ion.lower(), "case_" + case, "manifest.json")
        row = {"case": case, "mu_ne": mu_ne, "mu_ee": mu_ee, "ion": "Fe " + ion,
               "E0_eV": "", "E0_Hartree": "", "status": "missing", "ground_2J": "",
               "ground_parity": "", "n_levels": "", "mu_ne_used": "", "hxs_iter0": "",
               "ee_screen_option": "", "asymptotic_charge": "", "commit": ""}
        if os.path.exists(path):
            m = json.load(open(path))
            e = m.get("E0_eV")
            ok = m.get("status") == "ok" and e is not None and math.isfinite(e)
            # requested settings must match the case and be what FAC used
            ok = ok and m["mu_ne_requested"] == mu_ne and m["mu_ee_requested"] == mu_ee
            ok = ok and abs(m["mu_ne_used"] - mu_ne) < 1e-6 * mu_ne
            row.update({"E0_eV": e, "E0_Hartree": e / HARTREE_EV if e is not None else "",
                        "status": "ok" if ok else "bad: " + str(m.get("status")),
                        "ground_2J": m["ground_level"]["two_j"] if m.get("ground_level") else "",
                        "ground_parity": m["ground_level"]["parity"] if m.get("ground_level") else "",
                        "n_levels": m.get("n_levels"), "mu_ne_used": m["mu_ne_used"],
                        "hxs_iter0": m["hxs_iter0"], "ee_screen_option": m["ee_screen_option"],
                        "asymptotic_charge": m["asymptotic_charge"], "commit": m["commit"]})
            if ok:
                energies[(ion, case)] = m
        rows.append(row)

with open(os.path.join(run_dir, "ground_energies.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# ---- ionization potentials ----
diffs = []
ip_value = {}   # (ion, case) -> IP in eV
for case in CASES:
    for k, ion in enumerate(IONS[:-1]):
        upper = IONS[k + 1]
        low, high = energies.get((ion, case)), energies.get((upper, case))
        if low is None or high is None:
            continue
        # both energies must come from the same screening settings
        same = all(low[key] == high[key] for key in
                   ("mu_ne_requested", "mu_ee_requested", "hxs_iter0", "ee_screen_option",
                    "debye_mode_option", "commit", "density_A3"))
        if not same:
            continue
        ip = (high["E0_eV"] / HARTREE_EV - low["E0_eV"] / HARTREE_EV) * HARTREE_EV
        ip_value[(ion, case)] = ip
        diffs.append({"case": case, "ion": "Fe " + ion, "E0_lower_charge_eV": low["E0_eV"],
                      "E0_next_ion_eV": high["E0_eV"],
                      "delta_E_Hartree": (high["E0_eV"] - low["E0_eV"]) / HARTREE_EV,
                      "IP_eV": ip})

with open(os.path.join(run_dir, "ip_differences.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["case", "ion", "E0_lower_charge_eV", "E0_next_ion_eV",
                                      "delta_E_Hartree", "IP_eV"])
    w.writeheader()
    w.writerows(diffs)

# ---- comparison with Table 4 ----
def ee_effect(ion, mu_case_no_ee, mu_case_ee):
    a, b = ip_value.get((ion, mu_case_no_ee)), ip_value.get((ion, mu_case_ee))
    return None if a is None or b is None else b - a

out = []
for ion in IONS[:-1]:
    for case in CASES:
        calc = ip_value.get((ion, case))
        pub = TABLE4[ion][case]
        # e-e effect at the same mu_ne: B-A for mu_ne=0.10, D-C for mu_ne=0.25
        partner = {"A": ("A", "B"), "B": ("A", "B"), "C": ("C", "D"), "D": ("C", "D")}[case]
        ee_calc = ee_effect(ion, *partner)
        ee_pub = TABLE4[ion][partner[1]] - TABLE4[ion][partner[0]]
        # expected trends: e-e screening raises the IP; a larger mu_ne lowers it
        trend_ee = None if ee_calc is None else (ee_calc > 0) == (ee_pub > 0)
        other = {"A": "C", "C": "A", "B": "D", "D": "B"}[case]
        if (ion, case) in ip_value and (ion, other) in ip_value:
            lo, hi = sorted([case, other], key=lambda c: CASES[c][0])  # lo = smaller mu_ne
            trend_mu = ip_value[(ion, lo)] > ip_value[(ion, hi)]
        else:
            trend_mu = None
        out.append({
            "ion": "Fe " + ion, "case": case, "mu_ne": CASES[case][0], "mu_ee": CASES[case][1],
            "IP_calc_eV": "" if calc is None else round(calc, 4),
            "IP_table4_eV": pub,
            "diff_eV": "" if calc is None else round(calc - pub, 4),
            "diff_percent": "" if calc is None else round(100 * (calc - pub) / pub, 4),
            "ee_effect_calc_eV": "" if ee_calc is None else round(ee_calc, 4),
            "ee_effect_table4_eV": round(ee_pub, 4),
            "ee_trend_reproduced": "" if trend_ee is None else trend_ee,
            "mu_ne_trend_reproduced": "" if trend_mu is None else trend_mu,
        })

with open(os.path.join(run_dir, "table4_ip_comparison.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    w.writerows(out)

n_ok = sum(1 for r in rows if r["status"] == "ok")
print(f"{n_ok}/40 ground-state energies usable, {len(ip_value)}/36 IPs computed")
done = [r for r in out if r["diff_eV"] != ""]
if done:
    d = sorted(abs(r["diff_eV"]) for r in done)
    print(f"|IP - Table 4|: median {d[len(d)//2]:.2f} eV, max {d[-1]:.2f} eV")
