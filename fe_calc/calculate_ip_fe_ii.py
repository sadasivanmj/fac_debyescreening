
from pathlib import Path
import subprocess
import sys
import csv
import math


# --------------------------------------------------
# Project and output directory
# --------------------------------------------------

project_dir = Path(__file__).resolve().parent.parent
output_dir = project_dir / "fe_calc" / "fe_energy"
output_dir.mkdir(parents=True, exist_ok=True)

script_path = Path(__file__).resolve()


# --------------------------------------------------
# Configuration model
# --------------------------------------------------

configurations = {
    "fe_ii": [
        "1s2 2s2 2p6 3s2 3p6 3d6 4s1",
        "1s2 2s2 2p6 3s2 3p6 3d7",
        "1s2 2s2 2p6 3s2 3p6 3d5 4s2",
    ],

    "fe_iii": [
        "1s2 2s2 2p6 3s2 3p6 3d6",
        "1s2 2s2 2p6 3s2 3p6 3d5 4s1",
        "1s2 2s2 2p6 3s2 3p6 3d4 4s2",
    ],
}

electron_count = {
    "fe_ii": 25,
    "fe_iii": 24,
}

MU_VALUES = [0.0, 0.1, 0.25]
MU_STEP = 0.01

# Plasma parameters
DENSITY_A3 = 1e-6
BOHR_RADIUS_A = 0.52917721067
HARTREE_EV = 27.211386018

HXS_ITER0 = 8


# --------------------------------------------------
# Convert screening parameter to temperature
# --------------------------------------------------

def temperature_ev(mu):
    return (
        4.0 * math.pi * DENSITY_A3 * BOHR_RADIUS_A**3
        / mu**2 * HARTREE_EV
    )


# --------------------------------------------------
# Output filename helper
# --------------------------------------------------

def mu_tag(mu):
    return f"{mu:.2f}".replace(".", "p")


# --------------------------------------------------
# Calculate energy for one ion and one screening value
# --------------------------------------------------

def calculate_energy(ion, mu):

    from pfac import fac

    # Refuse a stale pfac from site-packages: the screened SCF needs the Debye RadialBound bracket and
    # radial:hxs_iter0, which only this checkout's build has.  Build with fe_calc/shared/build_pfac.sh and
    # `source fe_calc/shared/use_local_pfac.sh`.
    if not Path(fac.__file__).resolve().is_relative_to(project_dir / "build"):
        sys.exit(f"wrong FAC extension imported: {fac.__file__}\n"
                 "run: source fe_calc/shared/use_local_pfac.sh")

    fac.SetAtom("Fe")

    # Register the expanded configuration basis
    for config in configurations[ion]:
        fac.Config(config, group="g")

    # Start with the unscreened free-ion orbitals.
    # No ConfigEnergy(0)/(1): they store the free-ion configuration energies and Structure() then shifts the
    # screened configurations back onto them, which cancels the screening in E0 (the IP stayed at the
    # unscreened value for every mu).  Same sequence as the validated Fe II-IX study.
    fac.OptimizeRadial(["g"])

    # Gradually introduce Debye screening
    if mu > 0.0:

        fac.SetOption("orbital:debye_mode", 1)
        fac.SetOption("radial:ee_screen", 1)
        fac.SetOption("radial:hxs_iter0", HXS_ITER0)

        # Ionic charge: Fe II = +1, Fe III = +2
        zion = 26 - electron_count[ion]

        steps = []
        step = MU_STEP

        while step < mu - 1e-10:
            steps.append(round(step, 10))
            step += MU_STEP

        steps.append(mu)

        for mu_step in steps:
            fac.PlasmaScreen(
                zion,
                DENSITY_A3,
                temperature_ev(mu_step),
                1,
                0.0,
                1,
            )

            fac.OptimizeRadial(["g"])

            print(
                f"{ion}: screening mu = {mu_step:.2f}",
                flush=True,
            )

    # Calculate the structure with the screened orbitals
    tag = mu_tag(mu)

    en_file = output_dir / f"{ion}_{tag}.en"
    txt_file = output_dir / f"{ion}_{tag}.txt"

    fac.Structure(str(en_file), ["g"])
    fac.MemENTable(str(en_file))
    fac.PrintTable(str(en_file), str(txt_file), 1)

    print(f"Saved: {txt_file}", flush=True)


# --------------------------------------------------
# Read total energy from FAC's E0 line
# --------------------------------------------------

def read_energy(filename):

    with open(filename, "r") as file:
        for line in file:
            if line.strip().startswith("E0"):
                value = line.split("=", 1)[1].split(",", 1)[1]
                return float(value.strip())

    raise ValueError(f"E0 not found in {filename}")


# --------------------------------------------------
# Worker mode: one ion and one screening value
# --------------------------------------------------

if len(sys.argv) == 4 and sys.argv[1] == "--run":

    ion = sys.argv[2]
    mu = float(sys.argv[3])

    if ion not in configurations:
        raise ValueError(f"Unknown ion: {ion}")

    calculate_energy(ion, mu)
    sys.exit(0)


# --------------------------------------------------
# Run all calculations in independent FAC processes
# --------------------------------------------------

for mu in MU_VALUES:
    for ion in configurations:
        print(f"\nCalculating {ion}, mu = {mu:.2f}")

        subprocess.run(
            [
                sys.executable,
                str(script_path),
                "--run",
                ion,
                str(mu),
            ],
            check=True,
            cwd=project_dir,
        )


# --------------------------------------------------
# Calculate the ionization potentials
# --------------------------------------------------

results = []

for mu in MU_VALUES:

    tag = mu_tag(mu)

    E_FeII = read_energy(output_dir / f"fe_ii_{tag}.txt")
    E_FeIII = read_energy(output_dir / f"fe_iii_{tag}.txt")

    IP = E_FeIII - E_FeII

    results.append({
        "mu": mu,
        "E_FeII_eV": E_FeII,
        "E_FeIII_eV": E_FeIII,
        "IP_eV": IP,
    })


# --------------------------------------------------
# Print results
# --------------------------------------------------

IP_free = results[0]["IP_eV"]

print("\n============== RESULTS ==============")
print(f"{'mu':>8} {'E(Fe II) [eV]':>18} "
      f"{'E(Fe III) [eV]':>18} {'IP [eV]':>12} {'Delta IP [eV]':>15}")

for result in results:
    print(
        f"{result['mu']:8.2f} "
        f"{result['E_FeII_eV']:18.6f} "
        f"{result['E_FeIII_eV']:18.6f} "
        f"{result['IP_eV']:12.6f} "
        f"{result['IP_eV'] - IP_free:15.6f}"
    )


# --------------------------------------------------
# Save summary CSV
# --------------------------------------------------

summary_file = output_dir / "ip_screening_summary.csv"

with open(summary_file, "w", newline="") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "mu",
            "E_FeII_eV",
            "E_FeIII_eV",
            "IP_eV",
        ],
    )

    writer.writeheader()
    writer.writerows(results)

print(f"\nSummary saved to: {summary_file}")
print(f"All results saved in: {output_dir}")

