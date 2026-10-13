#!/usr/bin/env python3
"""Cross-model analysis of the deprince_t2_sfac.py results.

Reads results/ip_comparison.csv (written by the driver) and writes
    results/ip_lowering.csv    IP(mu) - IP(0): FAC vs Deprince et al.
    results/model_overview.md  all models side by side, IP lowering table,
                               error statistics
"""

import csv
import statistics
from collections import defaultdict
from pathlib import Path

import deprince_t2_sfac as drv

RES = drv.HERE / "results"
MU = drv.MU_VALUES
ORDER = [drv.ROMAN[i] for i in drv.TARGET_IONS]


def main():
    rows = list(csv.DictReader(open(RES / "ip_comparison.csv")))
    ip = {}
    for r in rows:
        if r["status"] == "ok":
            ip[(r["model"], r["ion"], float(r["mu"]))] = float(r["FAC_IP_eV"])
    status = {(r["model"], r["ion"], float(r["mu"])): r["status"] for r in rows}
    pub = {drv.ROMAN[i]: [float(x) for x in v[2:]] for i, v in drv.REFERENCE.items()}
    nist = {drv.ROMAN[i]: float(v[0]) for i, v in drv.REFERENCE.items()}
    models = [m for m in drv.MODELS if any(k[0] == m for k in status)]

    low = []
    for m in models:
        for ion in ORDER:
            for mu in MU[1:]:
                a, b = ip.get((m, ion, 0.0)), ip.get((m, ion, mu))
                p = pub[ion][MU.index(mu)] - pub[ion][0]
                rec = {"model": m, "ion": ion, "mu": mu, "pub_dIP_eV": round(p, 4)}
                if a is None or b is None:
                    rec["status"] = "FAILED/NOT RUN"
                else:
                    d = b - a
                    rec.update(status="ok", FAC_dIP_eV=round(d, 4),
                               diff_eV=round(d - p, 4),
                               rel_pct=round(100 * (d - p) / p, 2))
                low.append(rec)
    with open(RES / "ip_lowering.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["model", "ion", "mu", "status", "pub_dIP_eV",
                                "FAC_dIP_eV", "diff_eV", "rel_pct"])
        w.writeheader()
        w.writerows(low)

    def cell(m, ion, mu):
        v = ip.get((m, ion, mu))
        return f"{v:.3f}" if v is not None else status.get((m, ion, mu), "NOT RUN")

    md = ["# Model overview (computed FAC values vs quoted references)", ""]
    for mu in MU:
        md += [f"## IP (eV) at μ = {mu:g}", "",
               "| Ion | NIST | Deprince | " + " | ".join(models) + " |",
               "|---|---:|---:|" + "---:|" * len(models)]
        for ion in ORDER:
            n = f"{nist[ion]:g}" if mu == 0 else "—"
            md.append(f"| {ion} | {n} | {pub[ion][MU.index(mu)]:.2f} | "
                      + " | ".join(cell(m, ion, mu) for m in models) + " |")
        md.append("")
    md += ["## Plasma-induced IP lowering ΔIP = IP(μ) − IP(0) (eV)", "",
           "| Ion | μ | Deprince | " + " | ".join(models) + " |",
           "|---|---:|---:|" + "---:|" * len(models)]
    L = {(r["model"], r["ion"], r["mu"]): r for r in low}
    for ion in ORDER:
        for mu in MU[1:]:
            p = pub[ion][MU.index(mu)] - pub[ion][0]
            cells = []
            for m in models:
                r = L[(m, ion, mu)]
                cells.append(f"{r['FAC_dIP_eV']:.3f}" if r["status"] == "ok"
                             else "FAILED")
            md.append(f"| {ion} | {mu:g} | {p:.2f} | " + " | ".join(cells) + " |")
    md += ["", "## Error statistics over Fe II–Fe VIII", "",
           "| model | quantity | n | mean signed (eV) | mean abs (eV) | "
           "mean abs rel (%) | max abs rel (%) |",
           "|---|---|---:|---:|---:|---:|---:|"]
    for m in models:
        for label, ref, mus in (
                ("IP vs NIST, μ=0", lambda i, mu: nist[i], [0.0]),
                ("IP vs Deprince, μ=0", lambda i, mu: pub[i][0], [0.0]),
                ("IP vs Deprince, μ=0.1", lambda i, mu: pub[i][1], [0.1]),
                ("IP vs Deprince, μ=0.25", lambda i, mu: pub[i][2], [0.25]),
                ("ΔIP vs Deprince, μ=0.1,0.25",
                 lambda i, mu: pub[i][MU.index(mu)] - pub[i][0], [0.1, 0.25])):
            d, rel = [], []
            for ion in ORDER:
                for mu in mus:
                    if label.startswith("ΔIP"):
                        r = L[(m, ion, mu)]
                        if r["status"] != "ok":
                            continue
                        v = r["FAC_dIP_eV"]
                    else:
                        v = ip.get((m, ion, mu))
                        if v is None:
                            continue
                    t = ref(ion, mu)
                    d.append(v - t)
                    rel.append(100 * abs(v - t) / abs(t))
            if d:
                md.append(f"| {m} | {label} | {len(d)} | {statistics.mean(d):+.3f} | "
                          f"{statistics.mean(map(abs, d)):.3f} | "
                          f"{statistics.mean(rel):.2f} | {max(rel):.2f} |")
    (RES / "model_overview.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
