#!/usr/bin/env python3
"""Total energies of Fe II - Fe IX for Deprince et al. (2020), reduced configuration model.

This program calculates and saves TOTAL ENERGIES ONLY.  Ionization potentials are
computed by analyze_ionization_potentials.py from the saved energies.

One run = one ion x one orbital model (gref / ref) x one screening value mu, executed
in a fresh Python process.  Each run directory results/runs/<ion>/<model>/<case>/ holds
input.json (what was run), run.log, the raw FAC outputs (levels.txt, final.pot) and
manifest.json.  results/energies.csv is rebuilt from all manifests.

Examples
  python calculate_energies.py --preflight
  python calculate_energies.py --check-only
  python calculate_energies.py --ions fe_ii --models gref ref --mu 0 0.1     # a few runs
  python calculate_energies.py                                              # all 48 runs
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE))

import configurations as cfg                      # noqa: E402
import fac_run_support as sup                     # noqa: E402
import validate_inputs                            # noqa: E402

# Plasma settings (validated Fe II - Fe IV sequence).  Debye length lambda_D = sqrt(T/(4 pi n))
# in atomic units with ups = 0, so T = 4 pi n / mu^2.  The density is small enough that the
# ion-sphere radius lies far outside the atom.
DENSITY_A3 = 1e-6               # Angstrom^-3
BOHR_RADIUS_A = 0.52917721067   # FAC's value (consts.h)
HARTREE_EV_FAC = 27.211386018   # FAC's value, only used to turn T into eV for PlasmaScreen
MU_STEP = 0.05                  # screening is switched on in steps of 0.05
HXS_ITER0 = 8                   # validated exchange-ramp offset for warm restarts
DEFAULT_TIMEOUT_S = 1800.0


def temperature_ev(mu: float) -> float:
    return 4.0 * math.pi * DENSITY_A3 * BOHR_RADIUS_A**3 / mu**2 * HARTREE_EV_FAC


# ------------------------------------------------------------------------- the worker

def worker(run_dir: Path) -> int:
    """Runs inside the fresh process: the FAC sequence for one run.  Writes
    levels.txt, final.pot and worker_result.json into run_dir (the current directory)."""
    spec = json.loads((run_dir / "input.json").read_text())
    ion, model, mu = spec["ion"], spec["orbital_optimization_model"], spec["mu"]
    configs, electrons = spec["configurations"], spec["electron_count"]

    # Import FAC from the repository build only (never a globally installed pfac).
    build = sup.find_build_dir()
    if build is None:
        print("no build/lib.* found; run fe_calc/shared/build_pfac.sh")
        return 3
    sys.path.insert(0, str(build))
    from pfac import fac                          # the only place FAC is imported
    if not os.path.realpath(fac.__file__).startswith(str(sup.REPO_ROOT.resolve()) + os.sep):
        print("wrong FAC extension imported:", fac.__file__)
        return 3
    print("FAC extension:", fac.__file__, flush=True)

    fac.SetAtom("Fe")
    for name, value in spec.get("extra_options", []):    # diagnostic hook (--extra-option); empty in production runs
        fac.SetOption(name, value)
    if model == "ref":
        fac.Config(*configs, group="ref")
        optimize_groups, structure_groups = ["ref"], ["ref"]
    else:                                         # gref
        fac.Config(configs[0], group="g0")
        fac.Config(*configs[1:], group="ref")
        optimize_groups, structure_groups = ["g0"], ["g0", "ref"]

    # Free-ion SCF first; it also fixes the radial grid, so no ion-sphere cut-off is applied
    # once the plasma is switched on.  For mu = 0 this is the whole calculation: PlasmaScreen
    # is never called and no screening option is set.
    fac.OptimizeRadial(optimize_groups)
    steps: list[float] = []
    if mu > 0:
        fac.SetOption("orbital:debye_mode", 1)    # screen the electron-nucleus interaction
        fac.SetOption("radial:ee_screen", 1)      # screen the e-e Slater integrals AND the SCF direct potential
        fac.SetOption("radial:hxs_iter0", HXS_ITER0)
        zion = cfg.Z_IRON - electrons
        steps = [round(MU_STEP * k, 10) for k in range(1, int(round(mu / MU_STEP)))] + [mu]
        for mu_step in steps:
            fac.PlasmaScreen(zion, DENSITY_A3, temperature_ev(mu_step), 1, 0.0, 1)
            print(f"screening step mu = {mu_step:g}", flush=True)
            fac.OptimizeRadial(optimize_groups)

    fac.GetPotential("final.pot")
    fac.Structure("levels.b", structure_groups)   # all J and both parities
    fac.MemENTable("levels.b")
    fac.PrintTable("levels.b", "levels.txt", 1)

    # ---- read back what FAC actually did ----
    pot = sup.parse_potential(Path("final.pot"))
    table = sup.parse_level_table(Path("levels.txt"))
    ground = sup.ground_level(table)
    problems: list[str] = []
    checks: dict[str, bool] = {}

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks[name] = bool(ok)
        if not ok:
            problems.append(f"{name}: {detail}")

    check("electron_count_in_level_table", table["nele"] == electrons, f"NELE={table['nele']}")
    check("e0_finite", sup.is_finite(table["e0_ev"]), str(table["e0_ev"]))
    check("level_table_complete", table["nlev_header"] == len(table["levels"]),
          f"header {table['nlev_header']} vs parsed {len(table['levels'])}")
    check("level_table_referenced_to_lowest_level", ground is not None and abs(ground["energy_rel_ev"]) < 1e-9)
    if mu > 0:
        check("plasma_model_is_debye", pot["mps"] == 1, f"mps={pot['mps']}")
        check("nuclear_mu_applied", pot["dps"] > 0 and abs(pot["dps"] * mu - 1) < 1e-6, f"dps={pot['dps']}")
        check("no_ion_sphere_cut", pot["ips"] == 0, f"ips={pot['ips']}")
        # r*(Vc+U)+ZPS at large r is 0 only if BOTH the nuclear term and the e-e direct potential are screened
        if spec.get("direct_potential", "screened") == "screened":
            check("ee_direct_potential_screened", abs(pot["asymptotic_charge"]) < 1e-3, f"asymptotic charge {pot['asymptotic_charge']:.5f}")
        else:   # DensityToSZ experiment, variant A: Coulomb direct potential, so only the nucleus is screened
            check("ee_direct_potential_coulomb", abs(pot["asymptotic_charge"] - (electrons - 1.0)) < 1e-3,
                  f"asymptotic charge {pot['asymptotic_charge']:.5f}, expected {electrons - 1.0}")
    else:
        check("unscreened_path", pot["mps"] == -1, f"mps={pot['mps']} (expected -1: PlasmaScreen never called)")
        check("free_ion_asymptotic_charge", abs(pot["asymptotic_charge"] + (cfg.Z_IRON - electrons + 1)) < 1e-3,
              f"{pot['asymptotic_charge']:.5f}")

    result = {
        "total_energy_ev": (table["e0_ev"] + ground["energy_rel_ev"]) if (ground and sup.is_finite(table["e0_ev"])) else None,
        "ground_level": ground, "n_levels": len(table["levels"]), "nlev_header": table["nlev_header"],
        "blocks_seen": sorted({(lv["parity"], lv["two_j"]) for lv in table["levels"]}),
        "mu_used": (1.0 / pot["dps"]) if (pot["mps"] == 1 and pot["dps"] > 0) else 0.0,
        "mu_steps": steps, "mps": pot["mps"], "asymptotic_charge": pot["asymptotic_charge"],
        "fac_file": str(Path(os.path.realpath(fac.__file__)).relative_to(sup.REPO_ROOT.resolve())), "checks": checks, "problems": problems,
    }
    Path("worker_result.json").write_text(json.dumps(result, indent=2))
    return 0 if not problems else 4


# --------------------------------------------------------------------- the orchestrator

def plan_runs(ions: list[str], models: list[str], mus: list[float], direct_potential: str = "screened",
              extra_options: list | None = None) -> list[dict]:
    runs = []
    for ion in ions:
        for model in models:
            for mu in mus:
                electrons, configs, source = cfg.IONS[ion]
                runs.append({"ion": ion, "orbital_optimization_model": model, "mu": mu,
                             "case": cfg.case_name(mu), "electron_count": electrons,
                             "configurations": configs, "configuration_source": source,
                             "direct_potential": direct_potential, "extra_options": extra_options or []})
    return runs


def run_dir_of(results_dir: Path, run: dict) -> Path:
    return results_dir / "runs" / cfg.slug(run["ion"]) / run["orbital_optimization_model"] / run["case"]


def build_manifest(run: dict, rel_dir: str, exit_code: int, elapsed: float, timed_out: bool,
                   run_dir: Path) -> dict:
    """Assemble the manifest from the worker result, the log and the process outcome."""
    log_path = run_dir / "run.log"      # renamed to run.FAILED.log below if the run did not succeed
    log_text = log_path.read_text(errors="replace") if log_path.exists() else ""
    scan = sup.scan_log(log_text)
    result_file = run_dir / "worker_result.json"
    result = json.loads(result_file.read_text()) if result_file.exists() else None

    if timed_out:
        status, conv = "timeout", "unknown"
    elif scan["max_iteration_reached"]:
        status, conv = "not_converged", "not_converged"
    elif exit_code != 0 and result is None:
        status, conv = "failed", "unknown"
    elif result is None or result["problems"] or result["total_energy_ev"] is None:
        status, conv = "check_failed", "converged"
    elif exit_code != 0:
        status, conv = "failed", "converged"
    else:
        status, conv = "ok", "converged"

    energy = result["total_energy_ev"] if result else None
    g = sup.git_info()
    gl = result["ground_level"] if result else None
    ground_id = (f"ILEV={gl['ilev']};2J={gl['two_j']};parity={gl['parity']};VNL={gl['vnl']};{gl['labels']}"
                 if gl else "")
    ok = status == "ok"
    final_log = sup.finalize_log_name(run_dir, ok)
    manifest = {
        "paper": cfg.PAPER, "ion": run["ion"], "ion_order": cfg.ION_ORDER.index(run["ion"]),
        "charge_state": cfg.charge_state(run["ion"]), "electron_count": run["electron_count"],
        "configuration_model": cfg.CONFIGURATION_MODEL,
        "orbital_optimization_model": run["orbital_optimization_model"],
        "case": run["case"], "case_order": list(cfg.CASES).index(run["case"]),
        "mu_ne": run["mu"], "mu_ee": run["mu"], "direct_potential": run["direct_potential"], "extra_options": run["extra_options"],
        "total_energy": energy, "energy_unit": "eV",
        "total_energy_hartree_converted": energy / sup.HARTREE_EV if energy is not None else None,
        "ground_level_identifier": ground_id,
        "convergence_status": conv, "status": status, "exit_code": exit_code,
        "elapsed_time": round(elapsed, 3),
        "fac_executable_or_module": result["fac_file"] if result else "",
        "fac_source_commit": g["fac_source_commit"], "repo_head_commit": g["commit"], "build_provenance": sup.build_provenance_summary(),
        "run_directory": rel_dir, "log_path": f"{rel_dir}/{final_log.name}",
        # details
        "git_branch": g["branch"], "git_modified_tracked_files": g["modified_tracked_files"],
        "configurations": run["configurations"], "configuration_source": run["configuration_source"],
        "options": {"debye_mode": 1 if run["mu"] > 0 else 0, "ee_screen": 1 if run["mu"] > 0 else 0,
                    "hxs_iter0": HXS_ITER0 if run["mu"] > 0 else 0, "density_A3": DENSITY_A3,
                    "temperature_eV": temperature_ev(run["mu"]) if run["mu"] > 0 else None},
        "mu_used": result["mu_used"] if result else None, "mu_steps": result["mu_steps"] if result else [],
        "n_levels": result["n_levels"] if result else None,
        "blocks_seen_parity_2J": result["blocks_seen"] if result else None,
        "asymptotic_charge": result["asymptotic_charge"] if result else None,
        "checks": result["checks"] if result else {}, "problems": result["problems"] if result else [],
        "suspicious_log_lines": scan["suspicious_lines"],
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def run_one(run: dict, results_dir: Path, timeout_s: float, overwrite: bool) -> tuple[str, dict | None]:
    """Execute one run; returns (action, manifest).  Completed runs are never overwritten."""
    run_dir = run_dir_of(results_dir, run)
    rel = str(run_dir.relative_to(results_dir))
    existing = sup.read_manifest(run_dir)
    if existing and existing.get("status") == "ok" and not overwrite:
        return "skipped (already ok)", existing
    if run_dir.exists():
        sup.archive_previous(run_dir, "superseded" if existing and existing.get("status") == "ok" else "failed-attempt")
    run_dir.mkdir(parents=True)
    spec = {k: run[k] for k in ("ion", "orbital_optimization_model", "mu", "electron_count", "configurations", "direct_potential", "extra_options")}
    (run_dir / "input.json").write_text(json.dumps(spec, indent=2))
    cmd = [sup.python_executable(), str(Path(__file__).resolve()), "--worker", str(run_dir)]
    exit_code, elapsed, timed_out = sup.execute_worker(cmd, run_dir, timeout_s)
    manifest = build_manifest(run, rel, exit_code, elapsed, timed_out, run_dir)
    return manifest["status"], manifest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", type=Path, default=HERE / "results")
    ap.add_argument("--ions", nargs="+", default=cfg.ION_ORDER, help="e.g. fe_ii fe_iii (default: Fe II - Fe IX)")
    ap.add_argument("--models", nargs="+", choices=cfg.ORBITAL_MODELS, default=list(cfg.ORBITAL_MODELS))
    ap.add_argument("--mu", nargs="+", type=float, default=sorted(cfg.CASES.values()), help="screening values (0, 0.1, 0.25)")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_S, help="per-run timeout in seconds")
    ap.add_argument("--jobs", type=int, default=1, help="concurrent runs (each is a separate single-thread process)")
    ap.add_argument("--direct-potential", choices=("screened", "coulomb"), default="screened",
                    help="what the build's DensityToSZ is expected to do (default: screened). 'coulomb' is only for the "
                         "DensityToSZ A/B experiment, variant A, and only changes the expected asymptotic charge")
    ap.add_argument("--extra-option", nargs=2, action="append", default=[], metavar=("NAME", "VALUE"),
                    help="DIAGNOSTIC ONLY: FAC option set before the SCF (value parsed as a number), recorded in the manifest")
    ap.add_argument("--overwrite", action="store_true", help="re-run completed runs (the old directory is moved aside, not deleted)")
    ap.add_argument("--check-only", action="store_true", help="validate the configuration data and stop")
    ap.add_argument("--preflight", action="store_true", help="check build, options, provenance, output directory and stop")
    ap.add_argument("--skip-preflight", action="store_true", help="do not run the preflight before calculating")
    ap.add_argument("--dry-run", action="store_true", help="print the planned runs and stop")
    ap.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker:
        return worker(args.worker.resolve())

    problems = validate_inputs.validate()
    if args.check_only:
        for p in problems:
            print("PROBLEM:", p)
        print("configuration check:", "FAILED" if problems else "passed")
        return 1 if problems else 0
    if problems:
        print("configuration validation failed:", *problems, sep="\n  ")
        return 1

    results_dir = args.results_dir.resolve()
    rows = sup.preflight_checks(results_dir=results_dir)
    if args.preflight or not args.skip_preflight:
        for status, name, detail in rows:
            print(f"[{status}] {name}: {detail}")
        if any(s == "FAIL" for s, _, _ in rows):
            print("preflight FAILED")
            return 2
        if args.preflight:
            print("preflight passed")
            return 0

    ions = [cfg.normalize_ion(i) for i in args.ions]
    plan = plan_runs(ions, args.models, args.mu, args.direct_potential, [[n, float(v)] for n, v in args.extra_option])
    unknown = [m for m in args.mu if m not in cfg.CASES.values()]
    if unknown:
        print("unsupported mu values:", unknown, "(use 0, 0.1, 0.25)")
        return 1
    print(f"{len(plan)} planned run(s)")
    if args.dry_run:
        for r in plan:
            print("  ", r["ion"], r["orbital_optimization_model"], r["case"])
        return 0

    results_dir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    counts: dict[str, int] = {}

    def task(run: dict) -> str:
        action, manifest = run_one(run, results_dir, args.timeout, args.overwrite)
        detail = f" E={manifest['total_energy']:.4f} eV" if manifest and manifest.get("total_energy") is not None else ""
        print(f"{run['ion']:8s} {run['orbital_optimization_model']:5s} {run['case']}: {action}{detail}", flush=True)
        return action

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for action in pool.map(task, plan):
            counts[action] = counts.get(action, 0) + 1
    path = sup.write_energies_csv(results_dir)
    print(f"finished in {time.monotonic() - start:.0f} s: {counts}; wrote {path}")
    bad = sum(v for k, v in counts.items() if k not in ("ok", "skipped (already ok)"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
