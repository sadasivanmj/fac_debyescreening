#!/usr/bin/env python3
"""FAC calculation of Deprince et al. (2019), Tables 2 and 4.

Calculates E0 for Fe XVII--Fe XXV and the Fe XXVI daughter ion.  The RAS-like
configuration list is built as all unique one- and two-electron promotions from
the published reference configurations into 3s, 3p, and 3d.

Important: for Table 4's ee-screened cases, use a FAC/PFAC build where the
screened electron-electron direct potential is included in DensityToSZ (the
`ee-screening-scf` implementation), not only in the Slater integrals.
Run --preflight first to check that screening is applied to one-electron Fe XXVI.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

Z = 26.0
NPS_A3 = 1.0e-2
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018
UPS = 0.0
VXF = 1

IONS = ["FeXVII", "FeXVIII", "FeXIX", "FeXX", "FeXXI",
        "FeXXII", "FeXXIII", "FeXXIV", "FeXXV", "FeXXVI"]
ELECTRONS = {
    "FeXVII": 10, "FeXVIII": 9, "FeXIX": 8, "FeXX": 7,
    "FeXXI": 6, "FeXXII": 5, "FeXXIII": 4, "FeXXIV": 3,
    "FeXXV": 2, "FeXXVI": 1,
}

# Reference configurations from Deprince et al. (2019), A&A 626, A83, Sec. 2.1.
# [1s], [2s], or [2p] in the paper means a hole in that subshell.
REFS = {
    "FeXVII": [
        "1s2 2s2 2p6",
        "1s2 2s2 2p5 3s1",
        "1s2 2s2 2p5 3p1",
        "1s1 2s2 2p6 3s1",
        "1s1 2s2 2p6 3p1",
    ],
    "FeXVIII": [
        "1s2 2s2 2p5",
        "1s1 2s2 2p6",
    ],
    "FeXIX": [
        "1s2 2s2 2p4", "1s2 2s1 2p5", "1s2 2p6",
        "1s1 2s2 2p5", "1s1 2s1 2p6",
    ],
    "FeXX": [
        "1s2 2s2 2p3", "1s2 2s1 2p4", "1s2 2p5",
        "1s1 2s2 2p4", "1s1 2s1 2p5", "1s1 2p6",
    ],
    "FeXXI": [
        "1s2 2s2 2p2", "1s2 2s1 2p3", "1s2 2p4",
        "1s1 2s2 2p3", "1s1 2s1 2p4", "1s1 2p5",
    ],
    "FeXXII": [
        "1s2 2s2 2p1", "1s2 2s1 2p2", "1s2 2p3",
        "1s1 2s2 2p2", "1s1 2s1 2p3", "1s1 2p4",
    ],
    "FeXXIII": [
        "1s2 2s2", "1s2 2s1 2p1", "1s2 2p2",
        "1s1 2s2 2p1", "1s1 2s1 2p2", "1s1 2p3",
    ],
    # (1s + 2s + 2p)^3: all eight allowed shell-occupation configurations.
    "FeXXIV": [
        "1s2 2s1", "1s2 2p1", "1s1 2s2", "1s1 2s1 2p1",
        "1s1 2p2", "2s2 2p1", "2s1 2p2", "2p3",
    ],
    "FeXXV": ["1s2", "1s1 2s1", "1s1 2p1"],
    # The paper stops at Fe XXV; Fe XXVI is needed as the IP daughter ion.
    "FeXXVI": ["1s1"],
}

# Deprince et al. (2019), Table 2.
T2_PAPER = {
    "FeXVII":  [1260.58, 1214.75, 1147.80],
    "FeXVIII": [1357.09, 1308.56, 1237.68],
    "FeXIX":   [1459.12, 1407.92, 1333.09],
    "FeXX":    [1573.48, 1519.58, 1440.80],
    "FeXXI":   [1689.13, 1632.53, 1549.75],
    "FeXXII":  [1797.82, 1738.54, 1651.80],
    "FeXXIII": [1950.49, 1888.61, 1798.43],
    "FeXXIV":  [2044.34, 1979.79, 1885.77],
    "FeXXV":   [8836.74, 8768.90, 8667.86],
}
NIST = {
    "FeXVII": 1262.7, "FeXVIII": 1357.8, "FeXIX": 1460.0,
    "FeXX": 1575.6, "FeXXI": 1687.0, "FeXXII": 1798.4,
    "FeXXIII": 1950.4, "FeXXIV": 2045.759, "FeXXV": 8828.1875,
}

# Deprince et al. (2019), Table 4. Order: (mu_ne, mu_ee) =
# (0.1,0), (0.1,0.1), (0.25,0), (0.25,0.25).
T4_PAPER = {
    "FeXVII":  [1190.72, 1214.75, 1089.22, 1147.80],
    "FeXVIII": [1287.21, 1308.56, 1185.55, 1237.68],
    "FeXIX":   [1389.10, 1407.92, 1287.41, 1333.09],
    "FeXX":    [1503.54, 1519.58, 1401.59, 1440.80],
    "FeXXI":   [1619.15, 1632.53, 1517.02, 1549.75],
    "FeXXII":  [1727.82, 1738.54, 1625.57, 1651.80],
    "FeXXIII": [1880.58, 1888.61, 1778.76, 1798.43],
    "FeXXIV":  [1974.41, 1979.79, 1872.58, 1885.77],
    "FeXXV":   [8766.28, 8768.90, 8661.14, 8667.86],
}

# FAC needs only five unique conditions to produce both tables.  Table 2 uses
# mu_ne=mu_ee; Table 4 additionally uses mu_ee=0 at each nonzero mu_ne.
CASES = [
    {"key": "mu0",          "mu_ne": 0.0,  "mu_ee": 0.0},
    {"key": "ne01_ee0",     "mu_ne": 0.1,  "mu_ee": 0.0},
    {"key": "ne01_ee01",    "mu_ne": 0.1,  "mu_ee": 0.1},
    {"key": "ne025_ee0",    "mu_ne": 0.25, "mu_ee": 0.0},
    {"key": "ne025_ee025",  "mu_ne": 0.25, "mu_ee": 0.25},
]
CASE_BY_KEY = {c["key"]: c for c in CASES}
T2_CASE = {0.0: "mu0", 0.1: "ne01_ee01", 0.25: "ne025_ee025"}
T4_CASES = ["ne01_ee0", "ne01_ee01", "ne025_ee0", "ne025_ee025"]

SHELLS = ("1s", "2s", "2p", "3s", "3p", "3d")
CAPACITY = {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10}
TARGETS = ("3s", "3p", "3d")
TOKEN = re.compile(r"^(\d+[spdf])(\d*)$")


def parse_config(config: str) -> tuple[int, ...]:
    occ = {s: 0 for s in SHELLS}
    for token in config.split():
        m = TOKEN.fullmatch(token)
        if not m:
            raise ValueError(f"Unrecognized FAC shell token {token!r} in {config!r}")
        shell, number = m.groups()
        if shell not in occ:
            raise ValueError(f"Shell {shell!r} is outside the n<=3 active set")
        occ[shell] = int(number) if number else 1
        if occ[shell] > CAPACITY[shell]:
            raise ValueError(f"Overfilled shell {shell!r} in {config!r}")
    return tuple(occ[s] for s in SHELLS)


def format_config(occ: tuple[int, ...]) -> str:
    return " ".join(f"{shell}{n}" for shell, n in zip(SHELLS, occ) if n > 0)


def expand_sd(refs: list[str]) -> list[str]:
    """Unique reference configs plus all 1-/2-electron promotions into n=3.

    This is an explicit occupation-level RAS construction for FAC.  It should
    be treated as a reproducible implementation of the paper's stated recipe;
    it is not guaranteed to reproduce GRASP's exact CSF selection/coupling.
    """
    all_occ = set()
    target_indices = [SHELLS.index(s) for s in TARGETS]

    for ref in refs:
        base = parse_config(ref)
        all_occ.add(base)

        # Singles.
        for src, count in enumerate(base):
            if count == 0:
                continue
            for dst in target_indices:
                if src == dst:
                    continue
                final = list(base)
                final[src] -= 1
                final[dst] += 1
                if all(0 <= final[i] <= CAPACITY[SHELLS[i]] for i in range(len(SHELLS))):
                    all_occ.add(tuple(final))

        # Doubles: choose two electrons from occupied shells, and promote them
        # to 3s/3p/3d. The two electrons may come from the same shell.
        src_shells = [i for i, count in enumerate(base) if count > 0]
        source_pairs = [p for p in itertools.combinations_with_replacement(src_shells, 2)
                        if all(p.count(i) <= base[i] for i in set(p))]
        for srcs in source_pairs:
            for dsts in itertools.combinations_with_replacement(target_indices, 2):
                for dests in set(itertools.permutations(dsts)):
                    if any(src == dst for src, dst in zip(srcs, dests)):
                        continue
                    final = list(base)
                    for src in srcs:
                        final[src] -= 1
                    for dst in dests:
                        final[dst] += 1
                    if all(0 <= final[i] <= CAPACITY[SHELLS[i]] for i in range(len(SHELLS))):
                        all_occ.add(tuple(final))

    return [format_config(occ) for occ in sorted(all_occ)]


def config_groups(ion: str) -> tuple[list[str], list[str]]:
    refs = REFS[ion]
    if ion == "FeXXVI":
        return refs, []  # ground-state daughter; no excited levels needed for E0
    all_configs = expand_sd(refs)
    ref_set = {format_config(parse_config(c)) for c in refs}
    sd = [c for c in all_configs if c not in ref_set]
    return refs, sd


def validate_reference_configs() -> None:
    for ion in IONS:
        n = ELECTRONS[ion]
        for config in REFS[ion]:
            occ = parse_config(config)
            if sum(occ) != n:
                raise ValueError(f"{ion} config has {sum(occ)} electrons, expected {n}: {config}")
        refs, sd = config_groups(ion)
        for config in refs + sd:
            if sum(parse_config(config)) != n:
                raise ValueError(f"Generated config has wrong electron count for {ion}: {config}")
        print(f"{ion}: {n} electrons; {len(refs)} references; {len(sd)} generated SD configs")


def temperature_eV(mu_ne: float) -> float:
    """Set T so lambda_D=1/mu_ne, at fixed density and ups=0."""
    density_au = NPS_A3 * RBOHR_ANG**3
    t_hartree = (1.0 / mu_ne)**2 * 4.0 * math.pi * density_au * (1.0 + UPS)
    return t_hartree * HARTREE_EV


def read_e0(path: Path) -> float:
    with path.open() as f:
        for line in f:
            if line.strip().startswith("E0"):
                try:
                    return float(line.split(",", 1)[1].strip())
                except (IndexError, ValueError) as exc:
                    raise RuntimeError(f"Cannot parse E0 line in {path}: {line.rstrip()}") from exc
    raise RuntimeError(f"No E0 header in {path}")


def worker(ion: str, case_key: str, outdir_arg: str) -> None:
    from pfac import fac

    case = CASE_BY_KEY[case_key]
    mu_ne, mu_ee = case["mu_ne"], case["mu_ee"]
    run_dir = Path(outdir_arg) / "runs" / case_key / ion
    run_dir.mkdir(parents=True, exist_ok=True)
    level_bin = run_dir / "levels.lev.b"
    level_txt = run_dir / "levels.lev"

    fac.SetAtom("Fe")
    fac.SetOption("orbital:debye_mode", 0)
    fac.SetOption("radial:ee_screen", 1 if mu_ee > 0 else 0)
    if mu_ne > 0:
        fac.SetOption("orbital:debye_mode", 1)
        charge = Z - ELECTRONS[ion]
        tps = temperature_eV(mu_ne)
        # PFAC binding: PlasmaScreen(zps, nps[A^-3], tps[eV], m, ups, vxf).
        # m=1 is the Debye model. In this table, mu_ee is either zero or mu_ne.
        fac.PlasmaScreen(charge, NPS_A3, tps, 1, UPS, VXF)
        print(f"{ion}: mu_ne={mu_ne:g}, mu_ee={mu_ee:g}, lambda_D={1/mu_ne:g} a0, T={tps:.8g} eV")
    else:
        print(f"{ion}: unscreened, mu_ne=0, mu_ee=0")

    refs, sd = config_groups(ion)
    for cfg in refs:
        fac.Config(cfg, group="ref")
    for cfg in sd:
        fac.Config(cfg, group="sd")

    fac.ConfigEnergy(0)
    fac.OptimizeRadial(["ref"])
    fac.ConfigEnergy(1)

    groups = ["ref", "sd"] if sd else ["ref"]
    fac.Structure(str(level_bin), groups)
    fac.MemENTable(str(level_bin))
    fac.PrintTable(str(level_bin), str(level_txt), 1)  # energies in eV

    e0 = read_e0(level_txt)
    metadata = {
        "ion": ion, "electrons": ELECTRONS[ion], "case": case_key,
        "mu_ne": mu_ne, "mu_ee": mu_ee,
        "temperature_eV": temperature_eV(mu_ne) if mu_ne > 0 else None,
        "reference_configurations": refs, "n_sd_configurations": len(sd),
        "E0_eV": e0, "binary_table": str(level_txt),
    }
    (run_dir / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"RESULT {ion} {case_key} E0_eV={e0:.10f}")


def run_one(script: Path, outdir: Path, ion: str, case_key: str) -> float:
    run_dir = outdir / "runs" / case_key / ion
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "worker.log"
    result = subprocess.run(
        [sys.executable, str(script), "--worker", ion, case_key, str(outdir)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    log_path.write_text(result.stdout or "")
    if result.returncode != 0:
        tail = "\n".join((result.stdout or "").splitlines()[-30:])
        raise RuntimeError(f"FAC failed for {ion}/{case_key}; see {log_path}\n{tail}")
    energy_file = run_dir / "levels.lev"
    return read_e0(energy_file)


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def error_fields(calc: float, ref: float) -> tuple[float, float]:
    delta = calc - ref
    return delta, abs(delta) / abs(ref) * 100.0


def write_results(outdir: Path, energies: dict[tuple[str, str], float]) -> None:
    energy_rows = []
    for case in CASES:
        for ion in IONS:
            key = (ion, case["key"])
            if key in energies:
                energy_rows.append([ion, ELECTRONS[ion], case["key"], case["mu_ne"],
                                    case["mu_ee"], f"{energies[key]:.10f}"])
    write_csv(outdir / "energies.csv",
              ["ion", "electrons", "case", "mu_ne_a0^-1", "mu_ee_a0^-1", "E0_eV"],
              energy_rows)

    # Table 2: use (mu_ne, mu_ee)=(mu,mu), plus unscreened baseline.
    t2_rows = []
    for ion in IONS[:-1]:
        idx = IONS.index(ion)
        daughter = IONS[idx + 1]
        for mu, paper_value in zip((0.0, 0.1, 0.25), T2_PAPER[ion]):
            case_key = T2_CASE[mu]
            if (ion, case_key) not in energies or (daughter, case_key) not in energies:
                continue
            ip = energies[(daughter, case_key)] - energies[(ion, case_key)]
            delta_paper, rel_paper = error_fields(ip, paper_value)
            nist = NIST.get(ion) if mu == 0.0 else None
            delta_nist = ip - nist if nist is not None else None
            t2_rows.append([ion, daughter, mu, case_key,
                            f"{energies[(ion, case_key)]:.10f}",
                            f"{energies[(daughter, case_key)]:.10f}",
                            f"{ip:.8f}", f"{paper_value:.8f}", f"{delta_paper:.8f}",
                            f"{rel_paper:.5f}", "" if nist is None else f"{nist:.8f}",
                            "" if delta_nist is None else f"{delta_nist:.8f}",
                            "CHECK_HLIKE_SCREENING" if ion == "FeXXV" and mu > 0 else ""])
    write_csv(outdir / "table2_errors.csv",
              ["ion", "daughter_ion", "mu", "case", "E0_ion_eV", "E0_daughter_eV",
               "IP_FAC_eV", "IP_Deprince_eV", "delta_vs_Deprince_eV",
               "abs_relative_error_vs_Deprince_pct", "NIST_unscreened_eV", "delta_vs_NIST_eV",
               "validation_note"], t2_rows)

    # Table 4: independently screened nuclear/electron-electron terms.
    t4_rows = []
    for ion in IONS[:-1]:
        idx = IONS.index(ion)
        daughter = IONS[idx + 1]
        for condition_index, case_key in enumerate(T4_CASES):
            if (ion, case_key) not in energies or (daughter, case_key) not in energies:
                continue
            ip = energies[(daughter, case_key)] - energies[(ion, case_key)]
            ref = T4_PAPER[ion][condition_index]
            delta, rel = error_fields(ip, ref)
            case = CASE_BY_KEY[case_key]
            note = "CHECK_HLIKE_SCREENING" if ion == "FeXXV" else ""
            t4_rows.append([ion, daughter, case["mu_ne"], case["mu_ee"], case_key,
                            f"{energies[(ion, case_key)]:.10f}",
                            f"{energies[(daughter, case_key)]:.10f}",
                            f"{ip:.8f}", f"{ref:.8f}", f"{delta:.8f}", f"{rel:.5f}", note])
    write_csv(outdir / "table4_errors.csv",
              ["ion", "daughter_ion", "mu_ne", "mu_ee", "case", "E0_ion_eV",
               "E0_daughter_eV", "IP_FAC_eV", "IP_Deprince_eV",
               "delta_eV", "abs_relative_error_pct", "validation_note"], t4_rows)

    # Screened-H-like control: Fe XXVI has one electron, so mu_ee must not affect it.
    checks = []
    e0 = energies.get(("FeXXVI", "mu0"))
    if e0 is not None:
        for mu, case_off, case_on in [(0.1, "ne01_ee0", "ne01_ee01"),
                                      (0.25, "ne025_ee0", "ne025_ee025")]:
            off = energies.get(("FeXXVI", case_off))
            on = energies.get(("FeXXVI", case_on))
            if off is not None:
                shift = off - e0
                first_order = Z * mu * HARTREE_EV
                ok = shift > 0.5 * first_order
                checks.append({"mu_ne": mu, "E0_shift_FeXXVI_eV": shift,
                               "first_order_Z_mu_Ha_eV": first_order,
                               "screening_shift_check": "PASS" if ok else "FAIL_CHECK_BUILD_OR_POTENTIAL"})
            if off is not None and on is not None:
                checks.append({"mu_ne": mu, "E0_ee_on_minus_off_eV": on - off,
                               "ee_screening_one_electron_check":
                                   "PASS" if abs(on - off) < 1.0e-3 else "FAIL_EE_SCREEN_AFFECTS_ONE_ELECTRON"})
    (outdir / "screening_checks.json").write_text(json.dumps(checks, indent=2) + "\n")

    # Human-readable run summary and caveats.
    notes = [
        "FAC Deprince Tables 2 and 4 calculation summary",
        "",
        f"Energies calculated: {len(energies)} ion/case pairs.",
        "Table 2 results: table2_errors.csv",
        "Table 4 results: table4_errors.csv",
        "Ground-state total energies: energies.csv",
        "One-electron screening diagnostics: screening_checks.json",
        "",
        "Model: FAC ConfigEnergy(0) -> OptimizeRadial(['ref']) -> ConfigEnergy(1) -> Structure(['ref','sd']).",
        "The SD configurations are generated at occupation level as unique single/double promotions into 3s, 3p, 3d.",
        "This follows the paper's written RAS recipe but may not reproduce GRASP/RATIP's exact CSF selection or EAL optimization.",
        "Table 4 requires a FAC build that screens the electron-electron direct potential in DensityToSZ as well as Slater integrals.",
        "The current fix-debye-ee-screening branch alone is insufficient for that; use the ee-screening-scf implementation.",
        "The screened Fe XXVI control must shift positively with mu_ne. If screening_checks.json says FAIL, the Fe XXV IP rows are not physically validated.",
        "For each ion, the orbitals are optimized on the reference configurations; the Structure calculation includes reference + generated SD configurations.",
        "No Breit or QED corrections are added by this script; it is an FAC comparison, not an exact GRASP92/RATIP reproduction.",
    ]
    (outdir / "README.txt").write_text("\n".join(notes) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description="FAC reproduction of Deprince et al. Tables 2 and 4")
    parser.add_argument("--check-only", action="store_true", help="validate reference configurations and print RAS sizes only")
    parser.add_argument("--preflight", action="store_true", help="run only Fe XXVI at mu=0 and mu_ne=0.1, mu_ee=0")
    parser.add_argument("--ions", help="comma-separated ion labels, e.g. FeXXV,FeXXVI; default: all")
    parser.add_argument("--cases", help="comma-separated case keys; default: all five cases")
    parser.add_argument("--worker", nargs=3, metavar=("ION", "CASE", "OUTDIR"), help=argparse.SUPPRESS)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.worker:
        ion, case_key, outdir = args.worker
        worker(ion, case_key, outdir)
        return

    print("Validating Deprince Table-1 reference configurations and generated SD lists...")
    validate_reference_configs()
    if args.check_only:
        return

    print("IMPORTANT: full Table 4 requires the ee-screening-scf implementation "
          "(screened direct potential in DensityToSZ).")
    print("Run --preflight first; do not launch all cases until Fe XXVI's "
          "one-electron screening shift passes.")

    script = Path(__file__).resolve()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = script.parent / f"fac_deprince_tables2_4_{stamp}"
    outdir.mkdir()

    if args.preflight:
        requested_ions = ["FeXXVI"]
        requested_cases = ["mu0", "ne01_ee0"]
    else:
        requested_ions = args.ions.split(",") if args.ions else IONS
        requested_cases = args.cases.split(",") if args.cases else [c["key"] for c in CASES]

    bad_ions = set(requested_ions) - set(IONS)
    bad_cases = set(requested_cases) - set(CASE_BY_KEY)
    if bad_ions or bad_cases:
        raise ValueError(f"Unknown ions={sorted(bad_ions)} or cases={sorted(bad_cases)}")

    (outdir / "run_plan.json").write_text(json.dumps({
        "ions": requested_ions, "cases": requested_cases,
        "binary": "PFAC currently imported by this Python environment",
        "screening_conditions": [CASE_BY_KEY[k] for k in requested_cases],
    }, indent=2) + "\n")

    energies = {}
    for case_key in requested_cases:
        for ion in requested_ions:
            run_dir = outdir / "runs" / case_key / ion
            log = run_dir / "worker.log"
            print(f"\nRunning {ion}, {case_key}; log: {log}", flush=True)
            energies[(ion, case_key)] = run_one(script, outdir, ion, case_key)
            print(f"  E0 = {energies[(ion, case_key)]:.8f} eV", flush=True)

    if args.preflight:
        e0 = energies[("FeXXVI", "mu0")]
        e1 = energies[("FeXXVI", "ne01_ee0")]
        shift = e1 - e0
        expected = Z * 0.1 * HARTREE_EV
        print("\nFe XXVI one-electron screening preflight")
        print(f"E0(mu=0)       = {e0:.8f} eV")
        print(f"E0(mu_ne=0.1)  = {e1:.8f} eV")
        print(f"Energy shift   = {shift:.8f} eV")
        print(f"First-order Z*mu*Eh = {expected:.8f} eV")
        if shift < 0.5 * expected:
            print("FAIL: screened nuclear potential may not be reaching the one-electron ion.")
            print("Do not run/interpret the full tables until the FAC source/build is audited.")
        else:
            print("PASS: one-electron screening shift is nonzero and has the expected sign.")
        print(f"Outputs: {outdir}")
        return

    write_results(outdir, energies)
    print(f"\nResults saved in: {outdir}")
    print("Read screening_checks.json before interpreting Fe XXV screened IPs.")
    if set(requested_ions) != set(IONS) or set(requested_cases) != set(CASE_BY_KEY):
        print("NOTE: partial run; only IP pairs for which both adjacent ions were calculated are written.")


if __name__ == "__main__":
    main()

