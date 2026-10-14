#!/usr/bin/env python3
"""Analysis of the DensityToSZ A/B experiment from SAVED energies only (no FAC import, no calculations).

Variant A: Coulomb direct potential inside DensityToSZ (screened Slater integrals and nuclear screening kept).
Variant B: screened Yukawa direct potential (the production implementation).
Both with the upstream jmaxk fix, identical configurations and settings (Fe II, Fe III, Fe IV; gref; mu = 0, 0.1, 0.25).

Reads   results/variant_A/energies.csv, results/variant_B/energies.csv (+ manifests),
        results/variant_A_diagnosis/variant_a_diagnosis.csv (small-mu ladder, supplementary),
        results/physical_checks/variant_{A,B}/energies.csv (Fe XVII-XXVI, supplementary).
Writes  energies.csv, ip_comparison.csv, analysis_numbers.json, experiment_tables.md and the SVG plots in plots/.

Exit code 0 = every expected energy present and valid; 2 = problems (listed); --allow-incomplete returns 0 but still reports.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE.parent / "deprince_2020_fe_ii_viii"))

import analysis_support as an          # noqa: E402
import reference_data as ref           # noqa: E402
import svg_plot                        # noqa: E402

ION_ORDER = ["Fe II", "Fe III", "Fe IV"]
CASES = {"mu0p00": 0.0, "mu0p10": 0.1, "mu0p25": 0.25}
IONIZATIONS = [("Fe II", "Fe III"), ("Fe III", "Fe IV")]
HA = an.HARTREE_EV


def load_variant(res: Path, variant: str):
    d = res / f"variant_{variant}"
    records = an.load_energy_records(d)
    expected = {(i, "reference_configurations", "gref", c) for i in ION_ORDER for c in CASES}
    usable, problems = an.check_integrity(records, expected, d)
    problems = [f"variant {variant}: {p}" for p in problems]
    return records, usable, problems


def energy(usable, ion, case):
    r = usable.get((ion, "reference_configurations", "gref", case))
    return None if r is None else r["total_energy"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=HERE / "results")
    ap.add_argument("--allow-incomplete", action="store_true")
    args = ap.parse_args()
    res = args.results_dir.resolve()
    (res / "plots").mkdir(exist_ok=True)

    recs, usable, problems = {}, {}, []
    for v in "AB":
        recs[v], usable[v], p = load_variant(res, v)
        problems += p

    # ------------------------------------------------------------------ energies.csv (all records, failures included)
    fields = ["variant", "ion", "electron_count", "configuration_model", "orbital_optimization_model", "case", "mu_ne", "mu_ee",
              "direct_potential", "total_energy", "energy_unit", "total_energy_hartree_converted", "ground_level_identifier",
              "convergence_status", "status", "exit_code", "elapsed_time", "fac_executable_or_module", "fac_source_commit",
              "build_provenance", "run_directory", "log_path"]
    rows = []
    for v in "AB":
        for r in recs[v]:
            m = json.loads((res / f"variant_{v}" / r["run_directory"] / "manifest.json").read_text())
            rows.append({"variant": v, **{k: r.get(k, "") for k in fields[1:]}, "direct_potential": m.get("direct_potential", ""),
                         "run_directory": f"variant_{v}/{r['run_directory']}", "log_path": f"variant_{v}/{r['log_path']}"})
    rows.sort(key=lambda r: (r["variant"], ION_ORDER.index(r["ion"]), r["mu_ne"]))
    an.write_csv(res / "energies.csv", fields, rows)

    E = {(v, ion, c): energy(usable[v], ion, c) for v in "AB" for ion in ION_ORDER for c in CASES}

    # --------------------------------------------------------------------------- delta E = E_B - E_A
    dE = {}
    for ion in ION_ORDER:
        for c in CASES:
            a, b = E[("A", ion, c)], E[("B", ion, c)]
            dE[(ion, c)] = None if a is None or b is None else b - a
            if c == "mu0p00" and dE[(ion, c)] is not None and abs(dE[(ion, c)]) > 1e-6:
                problems.append(f"unscreened control differs between variants for {ion}: dE = {dE[(ion, c)]} eV")

    # --------------------------------------------------------------------------- IPs per variant
    def ip(v, lo, hi, c):
        a, b = E[(v, lo, c)], E[(v, hi, c)]
        if a is None or b is None:
            return None
        lo_rec, hi_rec = usable[v][(lo, "reference_configurations", "gref", c)], usable[v][(hi, "reference_configurations", "gref", c)]
        value, why = an.ionization_potential(lo_rec, hi_rec)     # refuses incompatible records
        if value is None:
            problems.append(f"variant {v} {lo}->{hi} {c}: {why}")
        return value

    IP = {(v, lo, c): ip(v, lo, hi, c) for v in "AB" for lo, hi in IONIZATIONS for c in CASES}
    hi_of = dict(IONIZATIONS)
    out_rows = []
    for lo, hi in IONIZATIONS:
        for c, mu in CASES.items():
            paper, nist = ref.PAPER_TABLE2[lo][c], (ref.NIST_IP[lo] if mu == 0 else None)
            row = {"ionization": f"{lo} -> {hi}", "case": c, "mu": mu, "ip_paper_eV": paper, "ip_nist_eV": nist if nist is not None else ""}
            for v in "AB":
                x, x0 = IP[(v, lo, c)], IP[(v, lo, "mu0p00")]
                row[f"IP_{v}_eV"] = "" if x is None else round(x, 6)
                row[f"err_vs_paper_{v}_eV"] = "" if x is None else round(x - paper, 6)
                row[f"rel_err_vs_paper_{v}_percent"] = "" if x is None else round(100 * (x - paper) / paper, 4)
                row[f"shift_{v}_eV"] = "" if (x is None or x0 is None) else round(x - x0, 6)
                row[f"err_vs_nist_{v}_eV"] = "" if (x is None or nist is None) else round(x - nist, 6)
            row["shift_paper_eV"] = round(paper - ref.PAPER_TABLE2[lo]["mu0p00"], 6)
            for v in "AB":
                s = row[f"shift_{v}_eV"]
                row[f"shift_err_{v}_eV"] = "" if s == "" else round(s - row["shift_paper_eV"], 6)
            a, b = IP[("A", lo, c)], IP[("B", lo, c)]
            if a is not None and b is not None:
                dI = b - a
                dI_from_dE = dE[(hi, c)] - dE[(lo, c)]
                row.update({"delta_I_BminusA_eV": round(dI, 6), "delta_I_from_dE_eV": round(dI_from_dE, 6),
                            "identity_residual_eV": abs(dI - dI_from_dE),
                            "abs_err_change_A_to_B_eV": round(abs(b - paper) - abs(a - paper), 6)})
                if abs(dI - dI_from_dE) > 1e-8:
                    problems.append(f"identity dI = dE(hi) - dE(lo) violated for {lo} {c}: {dI} vs {dI_from_dE}")
            else:
                row["delta_I_BminusA_eV"] = ""
            row["status"] = "A and B available" if (a is not None and b is not None) else \
                ("B only (A did not converge)" if b is not None else "missing")
            out_rows.append(row)
    ip_fields = ["ionization", "case", "mu", "IP_A_eV", "IP_B_eV", "delta_I_BminusA_eV", "delta_I_from_dE_eV", "identity_residual_eV",
                 "ip_paper_eV", "err_vs_paper_A_eV", "err_vs_paper_B_eV", "rel_err_vs_paper_A_percent", "rel_err_vs_paper_B_percent",
                 "abs_err_change_A_to_B_eV", "ip_nist_eV", "err_vs_nist_A_eV", "err_vs_nist_B_eV",
                 "shift_A_eV", "shift_B_eV", "shift_paper_eV", "shift_err_A_eV", "shift_err_B_eV", "status"]
    an.write_csv(res / "ip_comparison.csv", ip_fields, out_rows)

    # ---------------------------------------------------------- supplementary: small-mu ladders (A and B converge)
    ip0 = {lo: IP[("B", lo, "mu0p00")] for lo, _ in IONIZATIONS}

    def build_ladder(subdir: str, with_shifts: bool):
        path = res / subdir / "variant_a_diagnosis.csv"
        scan = list(csv.DictReader(open(path))) if path.exists() else []
        S = {(r["variant"], r["ion"], float(r["mu"])): (float(r["E0_eV"]) if r["status"] == "ok" else None) for r in scan}
        rows_ = []
        for mu in sorted({k[2] for k in S}):
            row = {"mu": mu}
            for ion in ION_ORDER:
                a, b = S.get(("A", ion, mu)), S.get(("B", ion, mu))
                row[f"dE_{ion}"] = None if a is None or b is None else b - a
                row[f"A_converged_{ion}"] = a is not None
            for lo, hi in IONIZATIONS:
                a2, a3 = S.get(("A", lo, mu)), S.get(("A", hi, mu))
                b2, b3 = S.get(("B", lo, mu)), S.get(("B", hi, mu))
                if None not in (a2, a3, b2, b3):
                    row[f"IP_A_{lo}"], row[f"IP_B_{lo}"] = a3 - a2, b3 - b2
                    row[f"dI_{lo}"] = (b3 - b2) - (a3 - a2)
                    row[f"dI_{lo}_from_dE"] = (b3 - a3) - (b2 - a2)
                    if abs(row[f"dI_{lo}"] - row[f"dI_{lo}_from_dE"]) > 1e-8:
                        problems.append(f"ladder identity violated {lo} mu={mu} ({subdir})")
                else:
                    row[f"dI_{lo}"] = None
                for v in "AB":
                    x = row.get(f"IP_{v}_{lo}")
                    row[f"shift{v}_{lo}"] = (x - ip0[lo]) if (with_shifts and x is not None and ip0[lo] is not None) else None
            rows_.append(row)
        return rows_

    sup_rows = build_ladder("variant_A_diagnosis", True)             # FAC default radial:csi = 0.55
    sup_rows0 = build_ladder("variant_A_diagnosis_csi0", False)      # radial:csi = 0 in both builds

    # ---------------------------------------------------- supplementary: Fe XVII-XXVI direct-potential contribution
    phys = {}
    for v in "AB":
        pth = res / "physical_checks" / f"variant_{v}" / "energies.csv"
        if pth.exists():
            phys[v] = {(r["ion"], r["case"]): float(r["total_energy"]) for r in csv.DictReader(open(pth)) if r["status"] == "ok"}
    phys_rows = []
    if len(phys) == 2:
        for ion in sorted({k[0] for k in phys["A"]}, key=lambda s: len(s)):
            for mu_lbl, (nuc, full) in (("0.1", ("A", "B")), ("0.25", ("C", "D"))):
                try:
                    d_slater = phys["A"][(ion, full)] - phys["A"][(ion, nuc)]
                    d_total = phys["B"][(ion, full)] - phys["B"][(ion, nuc)]
                    phys_rows.append({"ion": ion, "mu": mu_lbl, "dE_ee_screening_slater_only_A_build": d_slater,
                                      "dE_ee_screening_total_B_build": d_total, "dE_direct_potential": d_total - d_slater,
                                      "nuclear_only_identical_in_both_builds": phys["A"][(ion, nuc)] == phys["B"][(ion, nuc)]})
                except KeyError:
                    pass

    numbers = {"problems": problems, "dE_BminusA_eV": {f"{i}@{c}": dE[(i, c)] for i in ION_ORDER for c in CASES},
               "IP": {f"{v}:{lo}@{c}": IP[(v, lo, c)] for v in "AB" for lo, _ in IONIZATIONS for c in CASES},
               "ladder_default_csi": sup_rows, "ladder_csi0": sup_rows0, "physical_checks": phys_rows}
    (res / "analysis_numbers.json").write_text(json.dumps(numbers, indent=2, default=str))

    # ----------------------------------------------------------------------------------------------- tables
    f2 = lambda x, d=3: an.fmt(x, d)
    t = ["## Total energies and A/B differences (eV)", "",
         an.md_table(["Ion", "mu", "E(A)", "E(B)", "dE = E(B) - E(A)"],
                     [[ion, CASES[c], f2(E[("A", ion, c)], 4), f2(E[("B", ion, c)], 4), f2(dE[(ion, c)], 4)] for ion in ION_ORDER for c in CASES]), "",
         "(A is missing wherever its SCF did not converge; nothing is substituted.)", "",
         "## Ionization potentials (eV)", "",
         an.md_table(["Transition", "mu", "IP(A)", "IP(B)", "paper", "err A", "err B", "|err| change A->B", "NIST", "shift A", "shift B", "shift paper"],
                     [[r["ionization"], r["mu"], f2(r["IP_A_eV"] or None, 4) if r["IP_A_eV"] != "" else "—",
                       f2(r["IP_B_eV"], 4) if r["IP_B_eV"] != "" else "—", r["ip_paper_eV"],
                       f2(r["err_vs_paper_A_eV"], 4) if r["err_vs_paper_A_eV"] != "" else "—",
                       f2(r["err_vs_paper_B_eV"], 4) if r["err_vs_paper_B_eV"] != "" else "—",
                       f2(r.get("abs_err_change_A_to_B_eV", ""), 4) if r.get("abs_err_change_A_to_B_eV", "") != "" else "—",
                       r["ip_nist_eV"] if r["ip_nist_eV"] != "" else "—",
                       f2(r["shift_A_eV"], 4) if r["shift_A_eV"] != "" else "—", f2(r["shift_B_eV"], 4) if r["shift_B_eV"] != "" else "—",
                       r["shift_paper_eV"]] for r in out_rows]), "",
         "## Supplementary: small-mu ladder where variant A converges, FAC default radial:csi (single SCF step from the free ion)", "",
         an.md_table(["mu", "dE Fe II", "dE Fe III", "dE Fe IV", "dI Fe II", "dI Fe III", "A converged (II/III/IV)", "shift A Fe II", "shift B Fe II"],
                     [[r["mu"], f2(r["dE_Fe II"], 3), f2(r["dE_Fe III"], 3), f2(r["dE_Fe IV"], 3), f2(r["dI_Fe II"], 3), f2(r["dI_Fe III"], 3),
                       "/".join("y" if r[f"A_converged_{i}"] else "n" for i in ION_ORDER), f2(r["shiftA_Fe II"], 3), f2(r["shiftB_Fe II"], 3)]
                      for r in sup_rows]), "",
         "## Supplementary: the same ladder with radial:csi = 0 in both builds (removes FAC's Coulomb-designed asymptotic-charge adjustment)", "",
         an.md_table(["mu", "dE Fe II", "dE Fe III", "dE Fe IV", "dI Fe II", "dI Fe III", "A converged (II/III/IV)"],
                     [[r["mu"], f2(r["dE_Fe II"], 3), f2(r["dE_Fe III"], 3), f2(r["dE_Fe IV"], 3), f2(r["dI_Fe II"], 3), f2(r["dI_Fe III"], 3),
                       "/".join("y" if r[f"A_converged_{i}"] else "n" for i in ION_ORDER)] for r in sup_rows0]), "",
         "## Supplementary: Fe XVII - Fe XXVI, electron-electron screening effect split into Slater and direct-potential parts (eV)", "",
         an.md_table(["Ion", "mu", "dE(e-e screening), A build (Slater only)", "dE(e-e screening), B build (Slater + direct)", "direct-potential part"],
                     [[r["ion"], r["mu"], f2(r["dE_ee_screening_slater_only_A_build"], 4), f2(r["dE_ee_screening_total_B_build"], 4), f2(r["dE_direct_potential"], 4)]
                      for r in phys_rows]), ""]
    (res / "experiment_tables.md").write_text("\n".join(t))

    # ------------------------------------------------------------------------------------------------ plots
    P = res / "plots"
    mus_main = [0.0, 0.1, 0.25]
    # 1 total-energy change vs mu (ladders where both converge; main grid points where A is absent are marked)
    ser = []
    for k, ion in enumerate(ION_ORDER):
        col = svg_plot.PALETTE[k]
        ser.append({"name": f"{ion}, FAC default csi", "x": [0.0] + [r["mu"] for r in sup_rows], "y": [0.0] + [r[f"dE_{ion}"] for r in sup_rows], "color": col})
        ser.append({"name": f"{ion}, csi = 0", "x": [0.0] + [r["mu"] for r in sup_rows0], "y": [0.0] + [r[f"dE_{ion}"] for r in sup_rows0], "color": col,
                    "dashed": True, "marker": "s", "hollow": True})
    ser.append({"name": "A did not converge (main grid)", "x": [0.1, 0.25], "y": [0.0, 0.0], "marker": "x", "line": False, "color": "#000000"})
    svg_plot.line_plot(P / "01_total_energy_change_A_to_B.svg", "Total-energy change from the screened direct potential", "mu (a.u.)",
                       "E(B) - E(A) (eV)", ser, note="Lines are broken where variant A did not converge; x marks mu = 0.1, 0.25 (A failed, no value). Dashed: radial:csi = 0.")
    # 2 IPs vs published
    for lo, hi in IONIZATIONS:
        ser = [{"name": "B (screened direct)", "x": mus_main, "y": [IP[("B", lo, c)] for c in CASES]},
               {"name": "A (Coulomb direct)", "x": mus_main, "y": [IP[("A", lo, c)] for c in CASES], "color": "#D55E00", "marker": "s"},
               {"name": "Deprince et al. 2020", "x": mus_main, "y": [ref.PAPER_TABLE2[lo][c] for c in CASES], "color": "#009E73", "marker": "^", "dashed": True},
               {"name": "NIST (mu=0)", "x": [0.0], "y": [ref.NIST_IP[lo]], "color": "#000000", "marker": "x", "line": False}]
        svg_plot.line_plot(P / f"02_ip_{lo.replace(' ', '_').lower()}_vs_published.svg", f"{lo} ionization potential", "mu (a.u.)", "IP (eV)", ser,
                           note="Variant A converged only at mu = 0 (single point); no line is drawn through failed runs.")
    # 3 delta I vs mu
    ser = []
    for k, (lo, _) in enumerate(IONIZATIONS):
        col = svg_plot.PALETTE[k]
        ser.append({"name": f"dI {lo}, FAC default csi", "x": [0.0] + [r["mu"] for r in sup_rows], "y": [0.0] + [r[f"dI_{lo}"] for r in sup_rows], "color": col})
        ser.append({"name": f"dI {lo}, csi = 0", "x": [0.0] + [r["mu"] for r in sup_rows0], "y": [0.0] + [r[f"dI_{lo}"] for r in sup_rows0], "color": col,
                    "dashed": True, "marker": "s", "hollow": True})
    svg_plot.line_plot(P / "03_delta_I_vs_mu.svg", "Direct-potential contribution to the IPs, dI = I(B) - I(A)", "mu (a.u.)", "dI (eV)", ser,
                       note="Defined only where variant A converges (csi default: mu <= 0.01 Fe II, <= 0.02 Fe III).")
    # 4 screening shifts
    for lo, _ in IONIZATIONS:
        ser = [{"name": "B", "x": [0.0] + mus_main[1:] + [r["mu"] for r in sup_rows], "y": [0.0] + [IP[("B", lo, c)] - IP[("B", lo, "mu0p00")] for c in list(CASES)[1:]] + [r[f"shiftB_{lo}"] for r in sup_rows], "line": False},
               {"name": "A", "x": [0.0] + [r["mu"] for r in sup_rows], "y": [0.0] + [r[f"shiftA_{lo}"] for r in sup_rows], "color": "#D55E00", "marker": "s", "line": False},
               {"name": "Deprince et al. 2020", "x": mus_main, "y": [ref.PAPER_TABLE2[lo][c] - ref.PAPER_TABLE2[lo]["mu0p00"] for c in CASES], "color": "#009E73", "marker": "^", "dashed": True}]
        svg_plot.line_plot(P / f"04_shift_{lo.replace(' ', '_').lower()}.svg", f"{lo}: screening-induced IP shift", "mu (a.u.)", "IP(mu) - IP(0) (eV)", ser,
                           note="B and A points are separate calculations; the published shifts are connected only for reference.")
    # 5/6 residuals vs the paper
    for key, title, ylabel, fname in (("err", "Absolute IP residual against Deprince et al. 2020", "calc - paper (eV)", "05_abs_residual_vs_paper.svg"),
                                      ("rel", "Relative IP residual against Deprince et al. 2020", "(calc - paper)/paper (%)", "06_rel_residual_vs_paper.svg")):
        ser = []
        for lo, _ in IONIZATIONS:
            for v, mk in (("B", "o"), ("A", "s")):
                vals = []
                for c in CASES:
                    x = IP[(v, lo, c)]
                    p = ref.PAPER_TABLE2[lo][c]
                    vals.append(None if x is None else ((x - p) if key == "err" else 100 * (x - p) / p))
                ser.append({"name": f"{lo} {v}", "x": mus_main, "y": vals, "marker": mk, "hollow": v == "A",
                            "line": v == "B"})
        svg_plot.line_plot(P / fname, title, "mu (a.u.)", ylabel, ser, hlines=[(0.0, "0")],
                           note="Variant A has a value only at mu = 0 (identical to B there).")

    # 7 independent radial-potential comparison
    pv = res / "potential_validation.csv"
    if pv.exists():
        data = [r for r in csv.DictReader(open(pv)) if r["mu"] in ("0.25",) and r["density"] in ("valence", "multishell")]
        for dens in ("valence", "multishell"):
            sel = [r for r in data if r["density"] == dens]
            r_, fac = [float(r["r"]) for r in sel], [float(r["rV_FAC"]) for r in sel]
            cf, dr = [float(r["rV_closed_form"]) for r in sel], [abs(float(r["residual_abs"])) for r in sel]
            svg_plot.line_plot(P / f"07a_potential_{dens}_mu0.25.svg", f"r*V_H(r), {dens} density, mu = 0.25: FAC vs independent reference",
                               "r (a0)", "r V_H (a.u.)", [{"name": "closed-form reference", "x": r_, "y": cf, "marker": None},
                                                         {"name": "FAC DensityToSZ", "x": r_, "y": fac, "marker": "o", "line": False, "color": "#D55E00"}],
                               logx=True, logy=True, ylim=(1e-30, 1e3), note="Both curves overlap; values below 1e-30 not drawn; see 07b for the residual.")
            svg_plot.line_plot(P / f"07b_residual_{dens}_mu0.25.svg", f"|FAC - reference| as a function of radius, {dens}, mu = 0.25",
                               "r (a0)", "|residual| of r V_H (a.u.)", [{"name": "absolute residual", "x": r_, "y": dr, "marker": None}],
                               logx=True, logy=True, ylim=(1e-14, 1.0),
                               note="Points below 1e-14 are not drawn (far tail, where r V_H itself is < 1e-14). Floor ~1e-8 of the peak = FAC's own quadrature floor.")

    print(f"{sum(1 for v in 'AB' for k in usable[v])}/18 usable energies (A: {len(usable['A'])}, B: {len(usable['B'])})")
    for p in problems:
        print("PROBLEM:", p)
    print("wrote", res / "energies.csv", res / "ip_comparison.csv", res / "experiment_tables.md", "and plots/")
    return 0 if (args.allow_incomplete or not problems) else 2


if __name__ == "__main__":
    sys.exit(main())
