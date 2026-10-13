#!/usr/bin/env python3
"""Calculate Fe II ionization potential with FAC's sfac executable, not pfac.

Creates and runs one native SFAC input file per ion and screening parameter.
Reference configurations are those listed in Deprince et al. Table 1. They are
reference configurations only; this does not generate the paper's full SD active
space / CSF expansion.

Run with:
    SFAC=/path/to/sfac python fe_ip_sfac.py
If SFAC is unset, ~/facinst/bin/sfac is used, matching the archive examples.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN_ROOT = HERE / "fe_ip_sfac_runs"
SFAC = os.path.expanduser(os.environ.get("SFAC", "~/facinst/bin/sfac"))

# Plasma conditions convention used by debye_analysis/scripts in the ZIP.
# PlasmaScreen arguments: (zps, nps[A^-3], tps[eV], m, ups, vxf).
NPS_A3 = 1.0e-2
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018
UPS = 0.0
VXF = 1

# The paper's Table 1 reference configurations.
# [n l] denotes a hole in the corresponding subshell.
CONFIGS = {
    "FeII": [
        "1s2 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d7 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s2",
    ],
    "FeIII": [
        "1s2 2s2 2p6 3s2 3p6 3d6",
        "1s2 2s2 2p6 3s2 3p5 3d7",
        "1s2 2s2 2p6 3s2 3p5 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p5 3d5 4s2",
        "1s2 2s2 2p5 3s2 3p6 3d7",
        "1s2 2s2 2p5 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p5 3s2 3p6 3d5 4s2",
        "1s1 2s2 2p6 3s2 3p6 3d7",
        "1s1 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s1 2s2 2p6 3s2 3p6 3d5 4s2",
    ],
}

# Ion net charges passed as zps in the ZIP's SFAC examples.
ZNET = {"FeII": 1.0, "FeIII": 2.0}
# None means ordinary free-ion path. Positive values are mu [a0^-1].
MU_VALUES = [None, 0.1, 0.25]


def label_mu(mu: float | None) -> str:
    if mu is None:
        return "free"
    return "mu" + format(mu, ".6g").replace(".", "p").replace("-", "m")


def tps_for_lambda(lambda_d_a0: float) -> float:
    """Return T [eV] matching lambda_D [a0] at fixed NPS and UPS=0."""
    density_au = NPS_A3 * RBOHR_ANG**3
    temperature_hartree = (
        lambda_d_a0**2 * 4.0 * math.pi * density_au * (1.0 + UPS)
    )
    return temperature_hartree * HARTREE_EV


def sf_input(ion: str, mu: float | None) -> str:
    tag = f"{ion}_{label_mu(mu)}"
    # A short group name avoids FAC's 16-character group-name limit.
    cfg_line = "Config('g1', " + ", ".join(repr(c) for c in CONFIGS[ion]) + ")"

    lines = [
        "SetAtom('Fe')",
    ]

    if mu is not None:
        lam = 1.0 / mu
        tps = tps_for_lambda(lam)
        lines += [
            "SetOption('orbital:debye_mode', 1)",
            "SetOption('radial:ee_screen', 1)",
        ]

    lines.append(cfg_line)

    if mu is not None:
        lines.append(
            f"PlasmaScreen({ZNET[ion]:.1f}, {NPS_A3:.10g}, "
            f"{tps:.14g}, 1, {UPS:.1f}, {VXF})"
        )

    lines += [
        "ConfigEnergy(0)",
        "OptimizeRadial(['g1'])",
        "ConfigEnergy(1)",
        f"Structure('{tag}.lev.b', ['g1'])",
        f"MemENTable('{tag}.lev.b')",
        f"PrintTable('{tag}.lev.b', '{tag}.lev', 1)",
        "",
    ]
    return "\n".join(lines)


def read_e0_eV(level_file: Path) -> float:
    """Read E0 from SFAC PrintTable verbose output (format flag 1; energy in eV)."""
    with level_file.open() as f:
        for line in f:
            if line.strip().startswith("E0"):
                try:
                    return float(line.split(",", 1)[1].strip())
                except (IndexError, ValueError) as exc:
                    raise RuntimeError(f"Could not parse E0 line in {level_file}: {line}") from exc
    raise RuntimeError(f"No E0 header found in {level_file}")


def run_one(ion: str, mu: float | None) -> float:
    tag = f"{ion}_{label_mu(mu)}"
    workdir = RUN_ROOT / tag
    workdir.mkdir(parents=True, exist_ok=True)

    sf_file = workdir / f"{tag}.sf"
    sf_file.write_text(sf_input(ion, mu))

    if os.path.sep not in SFAC and shutil.which(SFAC) is None:
        raise FileNotFoundError(
            f"Cannot find SFAC executable '{SFAC}'. Set it explicitly, for example:\n"
            "  SFAC=/absolute/path/to/sfac python fe_ip_sfac.py"
        )
    if os.path.sep in SFAC and not (Path(SFAC).is_file() and os.access(SFAC, os.X_OK)):
        raise FileNotFoundError(
            f"SFAC executable is not found or not executable: {SFAC}\n"
            "Set SFAC=/absolute/path/to/sfac and rerun."
        )

    print(f"\n[{ion}, mu={'free' if mu is None else mu}] running sfac in {workdir}", flush=True)
    proc = subprocess.run(
        [SFAC, sf_file.name],
        cwd=workdir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.stdout:
        print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n")

    if proc.returncode != 0:
        log = workdir / f"{tag}.log"
        log.write_text(proc.stdout or "")
        raise RuntimeError(
            f"sfac failed for {tag} (return code {proc.returncode}). Log: {log}"
        )

    level_file = workdir / f"{tag}.lev"
    if not level_file.exists():
        raise RuntimeError(
            f"sfac returned success but did not create {level_file}. "
            f"Check output above and the executable build."
        )

    e0 = read_e0_eV(level_file)
    print(f"{ion} E0 = {e0:.8f} eV", flush=True)
    return e0


def main() -> None:
    results = []
    for mu in MU_VALUES:
        e_feii = run_one("FeII", mu)
        e_feiii = run_one("FeIII", mu)
        ip = e_feiii - e_feii
        results.append((mu, e_feii, e_feiii, ip))

    print("\nFe II ionization potential: Fe II -> Fe III + e- (eV)")
    print(f"{'mu (a0^-1)':>12} {'E0 Fe II (eV)':>18} "
          f"{'E0 Fe III (eV)':>18} {'IP Fe II (eV)':>16}")
    for mu, e2, e3, ip in results:
        mu_label = "free" if mu is None else f"{mu:g}"
        print(f"{mu_label:>12} {e2:18.6f} {e3:18.6f} {ip:16.6f}")

    print(f"\nSFAC input/output files are in: {RUN_ROOT}")
    print("Note: this uses Table 1 reference configurations, not the paper's full SD active space.")


if __name__ == "__main__":
    main()

