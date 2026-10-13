# Fe IV ground-state energy from the Deprince et al. (2020) Table 1 reference
# configurations, with Debye-Hueckel plasma screening.
#
# Usage:  python3 fe_iv_reference.py MODEL MU OUTPUT_FOLDER
#   MODEL = ref  : orbitals optimized on the average of all reference configurations
#   MODEL = gref : orbitals optimized on the ground configuration only
#   MU    = screening parameter mu = 1/lambda_D in atomic units (0, 0.1, 0.25)
#
# Run one calculation per Python process: FAC keeps its orbitals and
# potential in global memory, so separate runs must not share a process.

import math
import os
import sys

from pfac import fac

ION_NAME = "Fe IV"
NUM_ELECTRONS = 23

# Table 1 of Deprince et al. (2020), hole notation written out in full.
# The first configuration is the ground configuration.
REFERENCE_CONFIGS = [
    "1s2 2s2 2p6 3s2 3p6 3d5",
    "1s2 2s2 2p6 3s2 3p5 3d6",
    "1s2 2s2 2p5 3s2 3p6 3d6",
    "1s1 2s2 2p6 3s2 3p6 3d6",
]

# Plasma settings.  In Debye mode (PlasmaScreen m = 1, ups = 0) FAC uses the
# Debye length lambda_D = sqrt(T / (4 pi n)) in atomic units, with the
# density given in A^-3 and T in eV.  We keep the density fixed and choose T
# so that lambda_D = 1/mu.  The density is small so that the ion-sphere
# radius is far outside the atom.
DENSITY_A3 = 1e-6
BOHR_RADIUS_A = 0.52917721067     # FAC's value (consts.h)
HARTREE_EV = 27.211386018         # FAC's value, used only to convert T to eV
MU_STEP = 0.05                    # screening is switched on in steps of 0.05
HXS_ITER0 = 8                     # validated exchange-ramp offset for warm restarts


def count_electrons(config):
    # "1s2 2s2 2p6 ..." -> 2 + 2 + 6 + ...
    return sum(int(shell[2:]) for shell in config.split())


def temperature_for_mu(mu):
    density_au = DENSITY_A3 * BOHR_RADIUS_A**3
    return 4.0 * math.pi * density_au / mu**2 * HARTREE_EV


model = sys.argv[1]
mu = float(sys.argv[2])
output_folder = sys.argv[3]

if model not in ("ref", "gref"):
    sys.exit(f"unknown model {model!r}; use ref or gref")
for config in REFERENCE_CONFIGS:
    if count_electrons(config) != NUM_ELECTRONS:
        sys.exit(f"{config} has {count_electrons(config)} electrons, expected {NUM_ELECTRONS}")

os.makedirs(output_folder)        # fails if the folder already exists
os.chdir(output_folder)
print("FAC extension:", fac.__file__, flush=True)

fac.SetAtom("Fe")
if model == "ref":
    fac.Config(*REFERENCE_CONFIGS, group="ref")
    optimize_groups = ["ref"]
    structure_groups = ["ref"]
else:
    fac.Config(REFERENCE_CONFIGS[0], group="g0")
    fac.Config(*REFERENCE_CONFIGS[1:], group="ref")
    optimize_groups = ["g0"]
    structure_groups = ["g0", "ref"]

# Free-ion SCF first.  It also fixes the radial grid, so no ion-sphere
# cut-off is applied when the plasma is switched on afterwards.
fac.OptimizeRadial(optimize_groups)

if mu > 0:
    # Screen both the electron-nucleus (debye_mode 1) and the electron-electron
    # interaction (ee_screen 1: Slater integrals and SCF direct potential).
    fac.SetOption("orbital:debye_mode", 1)
    fac.SetOption("radial:ee_screen", 1)
    # Start each warm SCF with (almost) full exchange; without this the first
    # iteration has no exchange and the 3d orbital can become unbound.
    fac.SetOption("radial:hxs_iter0", HXS_ITER0)
    zion = 26 - NUM_ELECTRONS
    # rounded so that 3 * 0.05 is exactly 0.15
    steps = [round(MU_STEP * k, 10) for k in range(1, int(round(mu / MU_STEP)))] + [mu]
    for mu_step in steps:
        fac.PlasmaScreen(zion, DENSITY_A3, temperature_for_mu(mu_step), 1, 0.0, 1)
        print(f"screening step mu = {mu_step:g}", flush=True)
        fac.OptimizeRadial(optimize_groups)

fac.GetPotential("final.pot")
fac.Structure("levels.b", structure_groups)
fac.MemENTable("levels.b")
fac.PrintTable("levels.b", "levels.txt", 1)

# Read back what FAC actually used.  The potential header has mps (plasma
# model, 1 = Debye) and dps (the Debye length in a0).  The last grid point
# gives the asymptotic charge r*(Vc + U) + ZPS: 0 when both screenings are on,
# -(zion + 1) for the free ion.
header = {}
last_row = None
for line in open("final.pot"):
    if line.startswith("#"):
        parts = line[1:].split("=")
        if len(parts) == 2:
            header[parts[0].strip()] = parts[1].strip()
    elif line.strip():
        last_row = line.split()
plasma_model = int(header["mps"])
debye_length = float(header["dps"])
asymptotic_charge = float(last_row[4]) + float(last_row[5]) + float(last_row[17])

# FAC writes level energies in eV; E0 is the lowest level.
energy_ev = None
electrons_in_file = None
for line in open("levels.txt"):
    if line.startswith("E0"):
        energy_ev = float(line.split(",")[1])
    if line.startswith("NELE"):
        electrons_in_file = int(line.split("=")[1])

problems = []
if energy_ev is None or not math.isfinite(energy_ev):
    problems.append("no finite E0 in levels.txt")
if electrons_in_file != NUM_ELECTRONS:
    problems.append(f"levels.txt has NELE = {electrons_in_file}")
if mu > 0:
    if plasma_model != 1 or abs(debye_length * mu - 1) > 1e-6:
        problems.append(f"screening not applied: mps = {plasma_model}, dps = {debye_length}")
    if abs(asymptotic_charge) > 1e-3:
        problems.append(f"r*V at r_max = {asymptotic_charge}, e-e screening missing")
    if int(header["ips"]) != 0:
        problems.append("ion-sphere cut-off active")
elif abs(asymptotic_charge + (26 - NUM_ELECTRONS + 1)) > 1e-3:
    problems.append(f"free-ion r*V at r_max = {asymptotic_charge}")

status = "ok" if not problems else "FAILED: " + "; ".join(problems)
with open("summary.txt", "w") as summary:
    summary.write(f"ion = {ION_NAME}\n")
    summary.write(f"electrons = {NUM_ELECTRONS}\n")
    summary.write(f"model = {model}\n")
    summary.write(f"mu_requested = {mu}\n")
    summary.write(f"mu_used = {1 / debye_length if debye_length > 0 else 0.0}\n")
    summary.write(f"plasma_model_mps = {plasma_model}\n")
    summary.write(f"hxs_iter0 = {HXS_ITER0 if mu > 0 else 0}\n")
    summary.write(f"E0_eV = {energy_ev}\n")
    summary.write(f"status = {status}\n")
print(open("summary.txt").read(), flush=True)
if problems:
    sys.exit("calculation failed: " + status)
