# Shared helper for the Fe XVII-XXVI ground-configuration calculations.
#
# One call to run_ground() does one ion, one screening case, one Python
# process (FAC keeps its orbitals and potential in global memory).
#
# Screening model (see README.md):
#   mu_ne : screening of the electron-nucleus interaction, -Z exp(-mu r)/r
#   mu_ee : screening of the electron-electron interaction, exp(-mu r12)/r12
# FAC derives both from ONE Debye length (dps = 1/mu).  The nuclear term is
# switched on by orbital:debye_mode = 1 and the e-e term by radial:ee_screen = 1,
# so only mu_ee = 0 (e-e Coulomb, unscreened) and mu_ee = mu_ne are possible.
# These are exactly the two choices of Deprince et al. (2019), Table 4.

import argparse
import json
import math
import os
import subprocess
import sys

Z_IRON = 26
DENSITY_A3 = 1e-6                 # plasma density in Angstrom^-3 (as in the Fe II-IV cases)
BOHR_RADIUS_A = 0.52917721067     # FAC's value (consts.h)
HARTREE_EV = 27.211386018         # FAC's value (consts.h), used for T in eV only
MU_STEP = 0.05                    # screening is switched on in steps of 0.05
HXS_ITER0 = 8                     # validated exchange-ramp offset for warm restarts

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def count_electrons(config):
    # "1s2 2s2 2p6" -> 2 + 2 + 6
    return sum(int(shell[2:]) for shell in config.split())


def temperature_for_mu(mu):
    # Debye length lambda_D = sqrt(T / (4 pi n)) (ups = 0), so T = 4 pi n / mu^2.
    density_au = DENSITY_A3 * BOHR_RADIUS_A**3
    return 4.0 * math.pi * density_au / mu**2 * HARTREE_EV


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", REPO, "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:
        return "unknown"


def read_level_table(filename):
    """Return (nele, E0 in eV, list of (ilev, energy_rel_eV, parity, vnl, two_j))."""
    nele = None
    e0 = None
    levels = []
    for line in open(filename):
        if line.startswith("NELE"):
            nele = int(line.split("=")[1])
        elif line.startswith("E0"):
            e0 = float(line.split(",")[1])     # "E0 = block, energy"
        else:
            p = line.split()
            # level row: ILEV IBASE ENERGY P VNL 2J <labels>
            if len(p) >= 6 and p[0].isdigit() and "E" in p[2]:
                levels.append((int(p[0]), float(p[2]), int(p[3]), int(p[4]), int(p[5])))
    return nele, e0, levels


def read_potential_header(filename):
    """Return (header dict, r*(Vc+U)+ZPS at the last grid point)."""
    header = {}
    last_row = None
    for line in open(filename):
        if line.startswith("#"):
            parts = line[1:].split("=")
            if len(parts) == 2:
                header[parts[0].strip()] = parts[1].strip()
        elif line.strip():
            last_row = line.split()
    asymptotic = float(last_row[4]) + float(last_row[5]) + float(last_row[17])
    return header, asymptotic


def run_ground(ion_name, config):
    parser = argparse.ArgumentParser(description=f"{ion_name} ground configuration {config}")
    parser.add_argument("--mu-ne", type=float, required=True, help="nuclear-electron mu (a.u.)")
    parser.add_argument("--mu-ee", type=float, required=True, help="electron-electron mu (a.u.)")
    parser.add_argument("--outdir", required=True, help="new output folder (must not exist)")
    args = parser.parse_args()
    mu_ne, mu_ee = args.mu_ne, args.mu_ee

    nele = count_electrons(config)
    zion = Z_IRON - nele            # charge of the ion, e.g. Fe XVII has zion = 16
    if mu_ne < 0 or mu_ee < 0 or (mu_ee != 0 and mu_ee != mu_ne):
        sys.exit("FAC supports mu_ee = 0 or mu_ee = mu_ne only (one Debye length)")
    if mu_ee > 0 and mu_ne == 0:
        sys.exit("mu_ee > 0 needs mu_ne > 0")

    os.makedirs(args.outdir)        # fails if the folder already exists
    os.chdir(args.outdir)

    from pfac import fac
    # Make sure this is the locally built extension and not a global install.
    if not os.path.realpath(fac.__file__).startswith(os.path.realpath(REPO) + os.sep):
        sys.exit(f"wrong FAC extension imported: {fac.__file__}")
    print("FAC extension:", fac.__file__, flush=True)

    fac.SetAtom("Fe")
    fac.Config(config, group="g0")

    # Free-ion SCF first.  It also fixes the radial grid, so no ion-sphere
    # cut-off is applied when the plasma is switched on afterwards.
    fac.OptimizeRadial(["g0"])

    steps = []
    if mu_ne > 0:
        fac.SetOption("orbital:debye_mode", 1)               # screen the electron-nucleus term
        fac.SetOption("radial:ee_screen", 1 if mu_ee > 0 else 0)  # e-e: Yukawa or plain Coulomb
        fac.SetOption("radial:hxs_iter0", HXS_ITER0)
        # Switch the screening on in steps of MU_STEP (warm SCF restarts).
        steps = [round(MU_STEP * k, 10) for k in range(1, int(round(mu_ne / MU_STEP)))] + [mu_ne]
        for mu_step in steps:
            fac.PlasmaScreen(zion, DENSITY_A3, temperature_for_mu(mu_step), 1, 0.0, 1)
            print(f"screening step mu = {mu_step:g}", flush=True)
            fac.OptimizeRadial(["g0"])

    fac.GetPotential("final.pot")
    fac.Structure("levels.b", ["g0"])
    fac.MemENTable("levels.b")
    fac.PrintTable("levels.b", "levels.txt", 1)

    # ---- read back what FAC actually did ----
    header, asymptotic = read_potential_header("final.pot")
    mps = int(header["mps"])
    dps = float(header["dps"])
    ips = int(header["ips"])
    mu_ne_used = 1.0 / dps if (mps == 1 and dps > 0) else 0.0
    nele_file, e0, levels = read_level_table("levels.txt")

    # The ground level is the one with the lowest energy; levels.txt energies
    # are relative to E0, so the ground level must have energy 0 and E0 is
    # its total energy in eV (negative).
    ground = min(levels, key=lambda lv: lv[1]) if levels else None
    energy_ev = e0 + ground[1] if (e0 is not None and ground) else None

    problems = []
    if energy_ev is None or not math.isfinite(energy_ev):
        problems.append("no finite ground-state energy in levels.txt")
    if nele_file != nele:
        problems.append(f"levels.txt has NELE = {nele_file}, expected {nele}")
    if ground and abs(ground[1]) > 1e-9:
        problems.append("level table is not referenced to its lowest level")
    if mu_ne > 0:
        if mps != 1 or abs(mu_ne_used - mu_ne) > 1e-6 * mu_ne:
            problems.append(f"nuclear screening not applied: mps = {mps}, dps = {dps}")
        if ips != 0:
            problems.append("ion-sphere cut-off active")
        # r*(Vc+U)+ZPS at large r: 0 when nucleus and e-e are screened,
        # N-1 when only the nucleus is screened, -(zion+1) for the free ion.
        expected = 0.0 if mu_ee > 0 else nele - 1.0
        if abs(asymptotic - expected) > 1e-3:
            problems.append(f"asymptotic charge {asymptotic:.5f}, expected {expected:.5f}")
    else:
        if abs(asymptotic + (zion + 1)) > 1e-3:
            problems.append(f"free-ion asymptotic charge {asymptotic:.5f}")
    status = "ok" if not problems else "FAILED: " + "; ".join(problems)

    manifest = {
        "ion": ion_name, "Z": Z_IRON, "electrons": nele, "config": config,
        "mu_ne_requested": mu_ne, "mu_ee_requested": mu_ee,
        "mu_ne_used": mu_ne_used,
        "ee_screen_option": (1 if mu_ee > 0 else 0) if mu_ne > 0 else 0,
        "debye_mode_option": 1 if mu_ne > 0 else 0,
        "hxs_iter0": HXS_ITER0 if mu_ne > 0 else 0,
        "plasma_model_mps": mps, "debye_length_dps": dps, "ion_sphere_ips": ips,
        "asymptotic_charge": asymptotic,
        "density_A3": DENSITY_A3, "temperature_eV": temperature_for_mu(mu_ne) if mu_ne > 0 else None,
        "mu_steps": steps,
        "ground_level": {"ilev": ground[0], "parity": ground[2], "vnl": ground[3], "two_j": ground[4]} if ground else None,
        "n_levels": len(levels),
        "E0_eV": energy_ev, "energy_unit": "eV (total energy, negative)",
        "hartree_eV": HARTREE_EV,
        "commit": git_commit(), "fac_extension": fac.__file__,
        "outdir": os.getcwd(), "status": status,
    }
    with open("manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2), flush=True)
    if problems:
        sys.exit("calculation failed: " + status)
