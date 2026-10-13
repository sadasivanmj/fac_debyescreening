
#!/usr/bin/env python3
from pathlib import Path
import math
import subprocess
import sys

from pfac import fac

OUT = Path(__file__).resolve().parent

# ---------------------------------------------------------
# Deprince et al. reference configurations
# ---------------------------------------------------------

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

# Net ionic charges
NET_CHARGE = {"FeII": 1.0, "FeIII": 2.0}

# Debye parameters following the project examples.
# nps is in Angstrom^-3, tps in eV.
NPS = 1.0e-2
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018
UPS = 0.0
VXF = 1

# None: ordinary free-ion calculation, without PlasmaScreen.
# 1e-6: long-Debye-length control using the screened implementation.
# 0.1 and 0.25: parameters used in the Deprince comparison.
MU_VALUES = [None, 0.1, 0.25]


def mu_tag(mu):
    if mu is None:
        return "free"
    return "mu_" + format(mu, ".6g").replace(".", "p").replace("-", "m")


def temperature_for_lambda(lambda_d):
    """Temperature in eV giving the target Debye length in Bohr radii."""
    density_au = NPS * RBOHR_ANG**3

    t_hartree = (
        lambda_d**2
        * 4.0 * math.pi
        * density_au
        * (1.0 + UPS)
    )

    return t_hartree * HARTREE_EV


def calculate_ion(ion, mu):
    """Calculate one ion in a fresh FAC process."""
    tag = mu_tag(mu)
    group = f"{ion}_{tag}_ref"

    lev_bin = OUT / f"{ion}_{tag}.lev.b"
    lev_txt = OUT / f"{ion}_{tag}.lev"

    fac.SetAtom("Fe")

    # Only activate screening for finite, positive mu.
    if mu is not None:
        if mu <= 0.0:
            raise ValueError("Screened calculations require mu > 0.")

        # Electron-nucleus and electron-electron Debye screening.
        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 0)

    for cfg in CONFIGS[ion]:
        fac.Config(cfg, group=group)

    if mu is not None:
        lambda_d = 1.0 / mu
        tps = temperature_for_lambda(lambda_d)

        # PFAC Python argument order:
        # PlasmaScreen(zps, nps, tps, m, ups, vxf)
        # m=1 selects Debye screening.
        fac.PlasmaScreen(
            NET_CHARGE[ion],
            NPS,
            tps,
            1,
            UPS,
            VXF,
        )

        print(
            f"{ion}: mu={mu:g}, lambda_D={lambda_d:g} a0, "
            f"T={tps:.6g} eV"
        )
    else:
        print(f"{ion}: ordinary free-ion calculation")

    fac.ConfigEnergy(0)
    fac.OptimizeRadial([group])
    fac.ConfigEnergy(1)

    fac.Structure(str(lev_bin), [group])
    fac.MemENTable(str(lev_bin))

    # Verbose format: energy values in eV.
    fac.PrintTable(str(lev_bin), str(lev_txt), 1)

    print(f"Saved {lev_txt}")


def read_e0_eV(filename):
    """Read the E0 reference energy from a verbose FAC level table."""
    with open(filename) as f:
        for line in f:
            if line.strip().startswith("E0"):
                return float(line.split(",")[1])

    raise RuntimeError(f"E0 not found in {filename}")


def main():
    script = str(Path(__file__).resolve())
    results = []

    for mu in MU_VALUES:
        tag = mu_tag(mu)
        energies = {}

        for ion in ("FeII", "FeIII"):
            # Isolate FAC's global state between separate ion calculations.
            mu_arg = "free" if mu is None else str(mu)

            subprocess.run(
                [
                    sys.executable, script,
                    "--worker", ion, mu_arg,
                ],
                check=True,
            )

            table = OUT / f"{ion}_{tag}.lev"
            energies[ion] = read_e0_eV(table)

        ip = energies["FeIII"] - energies["FeII"]

        results.append((mu, energies["FeII"], energies["FeIII"], ip))

    print("\n" + "=" * 76)
    print("Fe II ionization potential: Fe II -> Fe III + e-")
    print("Energy values in eV")
    print("=" * 76)
    print(f"{'mu (a0^-1)':>12} {'E0 Fe II':>18} "
          f"{'E0 Fe III':>18} {'IP Fe II':>16}")

    for mu, e2, e3, ip in results:
        label = "free" if mu is None else f"{mu:g}"
        print(f"{label:>12} {e2:18.6f} {e3:18.6f} {ip:16.6f}")

    print("=" * 76)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--worker":
        ion = sys.argv[2]
        mu = None if sys.argv[3] == "free" else float(sys.argv[3])

        if ion not in CONFIGS:
            raise ValueError(f"Unknown ion: {ion}")

        calculate_ion(ion, mu)
    else:
        main()

