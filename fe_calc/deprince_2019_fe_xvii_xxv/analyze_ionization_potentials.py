#!/usr/bin/env python3
"""Ionization potentials of Fe XVII - Fe XXV from saved total energies (Paper B).

Reads results/energies.csv and the run manifests written by calculate_energies.py.
It never imports FAC and never launches a calculation; a missing energy is reported,
not recomputed, replaced or guessed.

    I_q(case) = [E(Fe q+1, case) - E(Fe q, case)]       (same configuration model, same case, same build)

IP(Fe XXV) uses the matching Fe XXVI energy.  The two configuration models (ground_only,
reference_only) are analysed separately and never mixed.  Only models present in the energy
file are analysed (use --models to restrict); the baseline is ground_only.

Outputs (in --output-dir, default = the results directory):
    ionization_potentials.csv   one row per ionization x model x case
    comparison_summary.md       literature comparison (Table 2, Table 4, NIST), e-e screening
                                effects, shifts and error statistics

Exit code 0 = complete and valid, 2 = missing / invalid / incompatible energies
(the tables are still written, with the problems listed; use --allow-incomplete to get 0).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE))

import analysis_support as an        # noqa: E402
import configurations as cfg         # noqa: E402
import reference_data as ref         # noqa: E402

CSV_FIELDS = ["ionization", "ion", "configuration_model", "case", "mu_ne", "mu_ee", "ip_calc_eV",
              "ip_published_eV", "published_source", "diff_vs_published_eV", "rel_diff_vs_published_percent",
              "ip_nist_eV", "diff_vs_nist_eV", "rel_diff_vs_nist_percent",
              "shift_vs_unscreened_calc_eV", "shift_vs_unscreened_published_eV", "shift_diff_eV",
              "ee_effect_calc_eV", "ee_effect_published_eV", "ee_effect_diff_eV",
              "energy_lower_eV", "energy_upper_eV", "status"]
IP_IONS = cfg.ION_ORDER[:-1]        # Fe XVII .. Fe XXV (Fe XXVI is only the daughter ion)
EE_PARTNER = {"B": "A", "D": "C"}   # e-e effect = fully screened minus nuclear-only


def expected_keys(models: list[str]) -> set[tuple]:
    return {(ion, model, cfg.ORBITAL_OPTIMIZATION, case)
            for ion in cfg.ION_ORDER for model in models for case in cfg.CASES}


def compute_rows(usable: dict, models: list[str]) -> list[dict]:
    ip: dict[tuple, float] = {}          # (ion, model, case) -> IP in eV
    rows = []
    for model in models:
        for k, ion in enumerate(IP_IONS):
            upper_ion = cfg.ION_ORDER[k + 1]
            for case, (mu_ne, mu_ee) in cfg.CASES.items():
                lo = usable.get((ion, model, cfg.ORBITAL_OPTIMIZATION, case))
                hi = usable.get((upper_ion, model, cfg.ORBITAL_OPTIMIZATION, case))
                pub, source = ref.published(ion, case)
                row = {"ionization": f"{ion} -> {upper_ion}", "ion": ion, "configuration_model": model, "case": case,
                       "mu_ne": mu_ne, "mu_ee": mu_ee, "ip_published_eV": pub, "published_source": source,
                       "ip_nist_eV": ref.NIST_IP[ion] if case == "unscreened" else "", "status": "ok"}
                if lo is None or hi is None:
                    row["status"] = "NOT COMPUTED: missing or invalid energy"
                    rows.append(row)
                    continue
                value, why = an.ionization_potential(lo, hi)
                if value is None:
                    row["status"] = f"NOT COMPUTED: {why}"
                    rows.append(row)
                    continue
                ip[(ion, model, case)] = value
                row.update({"ip_calc_eV": round(value, 4), "energy_lower_eV": lo["total_energy"],
                            "energy_upper_eV": hi["total_energy"],
                            "diff_vs_published_eV": round(value - pub, 4),
                            "rel_diff_vs_published_percent": round(100 * (value - pub) / pub, 3)})
                if case == "unscreened":
                    nist = ref.NIST_IP[ion]
                    row.update({"diff_vs_nist_eV": round(value - nist, 4),
                                "rel_diff_vs_nist_percent": round(100 * (value - nist) / nist, 3)})
                rows.append(row)
    # derived quantities that need two IPs of the same model
    for r in rows:
        if r["status"] != "ok":
            continue
        ion, model, case = r["ion"], r["configuration_model"], r["case"]
        if case != "unscreened" and (ion, model, "unscreened") in ip:
            r["shift_vs_unscreened_calc_eV"] = round(ip[(ion, model, case)] - ip[(ion, model, "unscreened")], 4)
            r["shift_vs_unscreened_published_eV"] = round(r["ip_published_eV"] - ref.TABLE2_MU0[ion], 4)
            r["shift_diff_eV"] = round(r["shift_vs_unscreened_calc_eV"] - r["shift_vs_unscreened_published_eV"], 4)
        if case in EE_PARTNER and (ion, model, EE_PARTNER[case]) in ip:
            calc = ip[(ion, model, case)] - ip[(ion, model, EE_PARTNER[case])]
            pub = ref.TABLE4[ion][case] - ref.TABLE4[ion][EE_PARTNER[case]]
            r.update({"ee_effect_calc_eV": round(calc, 4), "ee_effect_published_eV": round(pub, 4),
                      "ee_effect_diff_eV": round(calc - pub, 4)})
    return rows


def report(rows: list[dict], problems: list[str], results_dir: Path, commits: set[str], models: list[str]) -> str:
    out = ["# Fe XVII - Fe XXV ionization potentials: reduced configuration models",
           "",
           "Generated by `analyze_ionization_potentials.py` from the saved energies; no FAC calculation was run.",
           "",
           "**These are reduced-model benchmarks, not a reproduction of the published MCDF calculation.** `ground_only` is the",
           "ground configuration alone; `reference_only` is the listed reference configurations (orbitals optimised on the",
           "ground configuration). Neither includes the single/double excitations of the paper's active space, extended-average-level",
           "optimisation, Breit or QED. The models are never mixed. Published values: Deprince et al., A&A 626, A83 (2019),",
           "Table 2 (isolated ion) and Table 4 (cases A-D, convention (mu_ne, mu_ee)); NIST applies to the isolated ion only.",
           "", "- Energy source: `energies.csv` in the results directory (with the run manifests)",
           f"- FAC source commit(s) in the energy records: {', '.join(sorted(c[:10] for c in commits))}", ""]
    out += ["## Input integrity", ""]
    out += [f"- PROBLEM: {p}" for p in problems] or ["- All expected energies present, unique, converged and consistent with their manifests."]
    out += [f"- PROBLEM: {r['ionization']} {r['configuration_model']} {r['case']}: {r['status']}" for r in rows if r["status"] != "ok"]
    out += [""]
    for model in models:
        mrows = [r for r in rows if r["configuration_model"] == model]
        out += [f"## Model `{model}`", "", "### IPs (eV)", "",
                an.md_table(["Ionization", "case (mu_ne, mu_ee)", "IP calc", "published", "source", "diff", "diff (%)",
                             "NIST", "diff vs NIST", "shift vs unscreened (calc / publ.)"],
                            [[r["ionization"], f"{r['case']} ({r['mu_ne']}, {r['mu_ee']})", an.fmt(r.get("ip_calc_eV")),
                              an.fmt(r["ip_published_eV"]), r["published_source"], an.fmt(r.get("diff_vs_published_eV")),
                              an.fmt(r.get("rel_diff_vs_published_percent"), 3),
                              an.fmt(r["ip_nist_eV"] if r["ip_nist_eV"] != "" else None, 2), an.fmt(r.get("diff_vs_nist_eV")),
                              (f"{an.fmt(r.get('shift_vs_unscreened_calc_eV'))} / {an.fmt(r.get('shift_vs_unscreened_published_eV'))}"
                               if r["case"] != "unscreened" else "")] for r in mrows]), ""]
        ee = [r for r in mrows if "ee_effect_calc_eV" in r]
        out += ["### Electron-electron screening effect (fully screened minus nuclear-only: B - A, D - C)", "",
                an.md_table(["Ionization", "mu", "calc (eV)", "published (eV)", "calc - published (eV)"],
                            [[r["ionization"], r["mu_ne"], an.fmt(r["ee_effect_calc_eV"]), an.fmt(r["ee_effect_published_eV"]),
                              an.fmt(r["ee_effect_diff_eV"])] for r in ee]), "",
                "The nuclear-only minus fully-screened difference (A - B, C - D) is the same quantity with the opposite sign.", ""]
        stats = []
        for label, sel, key in (("IP vs published, all cases", mrows, "diff_vs_published_eV"),
                                ("IP vs Table 2 (mu = 0)", [r for r in mrows if r["case"] == "unscreened"], "diff_vs_published_eV"),
                                ("IP vs Table 4 (A-D)", [r for r in mrows if r["case"] != "unscreened"], "diff_vs_published_eV"),
                                ("IP vs NIST (mu = 0)", [r for r in mrows if r["case"] == "unscreened"], "diff_vs_nist_eV"),
                                ("shift vs unscreened, calc - published", [r for r in mrows if r["case"] != "unscreened"], "shift_diff_eV"),
                                ("e-e effect, calc - published", ee, "ee_effect_diff_eV")):
            s = an.summary_stats([(f"{r['ionization']} {r['case']}", r.get(key)) for r in sel])
            stats.append([label, s["n"], an.fmt(s["mean_abs"]), an.fmt(s["max_abs"]), s["max_label"], an.fmt(s["mean_signed"])])
        out += ["### Error statistics", "",
                an.md_table(["Comparison", "n", "mean abs (eV)", "max abs (eV)", "where", "mean signed (eV)"], stats), ""]
    out += ["## Notes", "",
            "- Table 4 columns B and D equal Table 2 at mu = 0.1 and 0.25; the paper's isolated-ion value (Table 2, mu = 0) is the reference for the unscreened case.",
            "- IP(Fe XXV) = E(Fe XXVI) - E(Fe XXV) with matching screening; the Fe XXVI energy does not depend on mu_ee.",
            "- A constant offset between calculated and published values across the screening cases points to the unscreened energies, not to the screening.",
            ""]
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=HERE / "results")
    ap.add_argument("--output-dir", type=Path, default=None)
    ap.add_argument("--models", nargs="+", choices=cfg.CONFIGURATION_MODELS, default=None,
                    help="models to analyse (default: every model that has records, at least ground_only)")
    ap.add_argument("--allow-incomplete", action="store_true", help="exit 0 even if energies are missing or invalid")
    args = ap.parse_args()
    results_dir = args.results_dir.resolve()
    out_dir = (args.output_dir or results_dir).resolve()

    try:
        records = an.load_energy_records(results_dir)
    except an.AnalysisError as e:
        print("ERROR:", e)
        return 2
    models = args.models or [m for m in cfg.CONFIGURATION_MODELS if any(r["configuration_model"] == m for r in records)] or ["ground_only"]
    usable, problems = an.check_integrity(records, expected_keys(models), results_dir)
    rows = compute_rows(usable, models)
    commits = {r["fac_source_commit"] for r in usable.values()}

    out_dir.mkdir(parents=True, exist_ok=True)
    an.write_csv(out_dir / "ionization_potentials.csv", CSV_FIELDS, rows)
    (out_dir / "comparison_summary.md").write_text(report(rows, problems, results_dir, commits, models))

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"models {models}: {len(usable)}/{len(expected_keys(models))} usable energies, {n_ok}/{len(rows)} IPs computed")
    for p in problems:
        print("PROBLEM:", p)
    for r in rows:
        if r["status"] != "ok":
            print("PROBLEM:", r["ionization"], r["configuration_model"], r["case"], r["status"])
    print("wrote", out_dir / "ionization_potentials.csv", "and", out_dir / "comparison_summary.md")
    incomplete = bool(problems) or n_ok != len(rows)
    return 0 if (args.allow_incomplete or not incomplete) else 2


if __name__ == "__main__":
    sys.exit(main())
