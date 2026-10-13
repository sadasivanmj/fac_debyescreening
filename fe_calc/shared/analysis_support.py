"""Support code for the IP-analysis programs (standard library only).

Imports nothing from FAC and runs no calculations: it only reads the saved
energies.csv and run manifests written by the energy-calculation programs.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

HARTREE_EV = 27.211386          # eV per Hartree, as specified for the IP definition

# Fields that must be identical for the two energies of one IP.  Besides the
# screening condition and the model, the same FAC build and source must have produced both.
COMPATIBILITY_FIELDS = ("configuration_model", "orbital_optimization_model", "case",
                        "mu_ne", "mu_ee", "energy_unit", "fac_source_commit", "build_provenance")


class AnalysisError(Exception):
    """Raised for input problems that make a result untrustworthy."""


def load_energy_records(results_dir: Path) -> list[dict]:
    path = Path(results_dir) / "energies.csv"
    if not path.exists():
        raise AnalysisError(f"{path} not found; run calculate_energies.py first")
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for key in ("mu_ne", "mu_ee", "total_energy", "total_energy_hartree_converted", "elapsed_time"):
                r[key] = float(r[key]) if r.get(key) not in (None, "") else math.nan
            r["electron_count"] = int(r["electron_count"])
            rows.append(r)
    return rows


def record_key(r: dict) -> tuple:
    return (r["ion"], r["configuration_model"], r["orbital_optimization_model"], r["case"])


def check_integrity(records: list[dict], expected: set[tuple], results_dir: Path) -> tuple[dict, list[str]]:
    """Index the usable records and list every problem found.

    Returns (usable, problems).  usable maps record_key -> record for records that are
    status 'ok', finite, unique and consistent with their manifest.  Anything else is
    reported in problems; nothing is replaced or guessed.
    """
    problems: list[str] = []
    by_key: dict[tuple, list[dict]] = {}
    for r in records:
        by_key.setdefault(record_key(r), []).append(r)
    usable: dict[tuple, dict] = {}
    for key, rs in by_key.items():
        if key not in expected:
            problems.append(f"unexpected record {key}")
            continue
        if len(rs) > 1:
            problems.append(f"duplicated record {key}: {len(rs)} rows")
            continue
        r = rs[0]
        if r["status"] != "ok":
            problems.append(f"invalid run {key}: status={r['status']} convergence={r['convergence_status']}")
            continue
        if not math.isfinite(r["total_energy"]):
            problems.append(f"non-finite energy {key}")
            continue
        manifest_path = _resolve_manifest(results_dir, r)
        if manifest_path is None:
            problems.append(f"manifest missing for {key}")
            continue
        m = json.loads(manifest_path.read_text())
        if m.get("status") != "ok" or abs(m.get("total_energy", math.nan) - r["total_energy"]) > 1e-9:
            problems.append(f"energies.csv and manifest disagree for {key}")
            continue
        usable[key] = r
    for key in sorted(expected - set(by_key)):
        problems.append(f"missing energy {key}")
    return usable, problems


def _resolve_manifest(results_dir: Path, r: dict) -> Path | None:
    """run_directory is stored relative to the results directory (portable); accept absolute too."""
    p = Path(r["run_directory"])
    for cand in ((p / "manifest.json"), (Path(results_dir) / p / "manifest.json")):
        if cand.exists():
            return cand
    return None


def ionization_potential(lower: dict, upper: dict) -> tuple[float | None, str]:
    """IP in eV from two compatible energy records; (None, reason) if they are not compatible."""
    for field in COMPATIBILITY_FIELDS:
        if lower[field] != upper[field]:
            return None, f"incompatible {field}: {lower[field]!r} vs {upper[field]!r}"
    if upper["electron_count"] != lower["electron_count"] - 1:
        return None, "adjacent charge states expected (electron counts differ by != 1)"
    d_hartree = upper["total_energy"] / HARTREE_EV - lower["total_energy"] / HARTREE_EV
    return d_hartree * HARTREE_EV, ""


def summary_stats(pairs: list[tuple[str, float]]) -> dict:
    """pairs = (label, signed difference).  Returns n, mean/max absolute difference and where."""
    vals = [(l, d) for l, d in pairs if d is not None and math.isfinite(d)]
    if not vals:
        return {"n": 0, "mean_abs": math.nan, "max_abs": math.nan, "max_label": "", "mean_signed": math.nan}
    worst = max(vals, key=lambda p: abs(p[1]))
    return {"n": len(vals), "mean_abs": sum(abs(d) for _, d in vals) / len(vals),
            "max_abs": abs(worst[1]), "max_label": worst[0],
            "mean_signed": sum(d for _, d in vals) / len(vals)}


def fmt(x, digits: int = 2) -> str:
    return "—" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{digits}f}"


def md_table(headers: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
