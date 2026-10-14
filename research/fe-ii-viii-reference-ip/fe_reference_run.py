# Ground-state energy of one Fe ion (Fe II - Fe IX) from the reference
# configurations of Deprince et al. (reduced model: no single/double excitations),
# with Debye-Hueckel plasma screening.  Same FAC sequence as the validated
# Fe II - Fe IV cases in research/debye-screening/cases/fe/.
#
# Usage:  python fe_reference_run.py ION MODEL MU OUTPUT_FOLDER
#   ION   = "Fe II" ... "Fe IX"   (quote it) or fe_ii ... fe_ix
#   MODEL = ref  : orbitals optimized on the average of all reference configurations
#           gref : orbitals optimized on the ground configuration only
#   MU    = screening parameter mu = 1/lambda_D in atomic units (0, 0.1, 0.25)
#
# Run one calculation per Python process: FAC keeps its orbitals and potential
# in global memory.  The screened runs use debye_mode 1 (electron-nucleus) and
# ee_screen 1 (electron-electron) with the same mu.

import json
import math
import os
import subprocess
import sys

from reference_configs import IONS, count_electrons

DENSITY_A3 = 1e-6                 # small density: the ion sphere lies far outside the atom
BOHR_RADIUS_A = 0.52917721067     # FAC's value (consts.h)
HARTREE_EV = 27.211386018         # FAC's value, used only to convert T to eV
MU_STEP = 0.05                    # screening is switched on in steps of 0.05
HXS_ITER0 = 8                     # validated exchange-ramp offset for warm restarts
Z_IRON = 26

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def temperature_for_mu(mu):
    # lambda_D = sqrt(T / (4 pi n)) = 1/mu  =>  T = 4 pi n / mu^2   (atomic units)
    density_au = DENSITY_A3 * BOHR_RADIUS_A**3
    return 4.0 * math.pi * density_au / mu**2 * HARTREE_EV


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", REPO, "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


ion_arg, model, mu, output_folder = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
ion_name = "Fe " + ion_arg.replace("_", " ").split()[-1].upper()      # "fe_ii" or "Fe II" -> "Fe II"
if ion_name not in IONS:
    sys.exit(f"unknown ion {ion_arg!r}")
if model not in ("ref", "gref"):
    sys.exit(f"unknown model {model!r}; use ref or gref")
num_electrons, configs = IONS[ion_name]
for config in configs:
    if count_electrons(config) != num_electrons:
        sys.exit(f"{config} has {count_electrons(config)} electrons, expected {num_electrons}")

os.makedirs(output_folder)        # fails if the folder already exists
os.chdir(output_folder)

from pfac import fac
if not os.path.realpath(fac.__file__).startswith(os.path.realpath(REPO) + os.sep):
    sys.exit(f"wrong FAC extension imported: {fac.__file__}")
print("FAC extension:", fac.__file__, flush=True)

fac.SetAtom("Fe")
if model == "ref":
    fac.Config(*configs, group="ref")
    optimize_groups = ["ref"]
    structure_groups = ["ref"]
else:
    fac.Config(configs[0], group="g0")
    fac.Config(*configs[1:], group="ref")
    optimize_groups = ["g0"]
    structure_groups = ["g0", "ref"]

# Free-ion SCF first.  It also fixes the radial grid, so no ion-sphere
# cut-off is applied when the plasma is switched on afterwards.
fac.OptimizeRadial(optimize_groups)

steps = []
if mu > 0:
    fac.SetOption("orbital:debye_mode", 1)       # screen the electron-nucleus term
    fac.SetOption("radial:ee_screen", 1)         # screen the electron-electron term
    fac.SetOption("radial:hxs_iter0", HXS_ITER0)
    zion = Z_IRON - num_electrons
    steps = [round(MU_STEP * k, 10) for k in range(1, int(round(mu / MU_STEP)))] + [mu]
    for mu_step in steps:
        fac.PlasmaScreen(zion, DENSITY_A3, temperature_for_mu(mu_step), 1, 0.0, 1)
        print(f"screening step mu = {mu_step:g}", flush=True)
        fac.OptimizeRadial(optimize_groups)

fac.GetPotential("final.pot")
fac.Structure("levels.b", structure_groups)
fac.MemENTable("levels.b")
fac.PrintTable("levels.b", "levels.txt", 1)

# ---- read back what FAC actually used ----
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
# r*(Vc+U)+ZPS at the last grid point: 0 with both screenings on, -(zion+1) for the free ion
asymptotic_charge = float(last_row[4]) + float(last_row[5]) + float(last_row[17])

# FAC writes level energies in eV relative to E0, the energy of the lowest level.
energy_ev = None
electrons_in_file = None
n_levels = 0
lowest_relative = None
for line in open("levels.txt"):
    if line.startswith("E0"):
        energy_ev = float(line.split(",")[1])
    elif line.startswith("NELE"):
        electrons_in_file = int(line.split("=")[1])
    else:
        p = line.split()
        if len(p) >= 6 and p[0].isdigit() and "E" in p[2]:
            n_levels += 1
            e = float(p[2])
            lowest_relative = e if lowest_relative is None else min(lowest_relative, e)

problems = []
if energy_ev is None or not math.isfinite(energy_ev):
    problems.append("no finite E0 in levels.txt")
if electrons_in_file != num_electrons:
    problems.append(f"levels.txt has NELE = {electrons_in_file}")
if lowest_relative is None or abs(lowest_relative) > 1e-9:
    problems.append("level table is not referenced to its lowest level")
if mu > 0:
    if plasma_model != 1 or abs(debye_length * mu - 1) > 1e-6:
        problems.append(f"screening not applied: mps = {plasma_model}, dps = {debye_length}")
    if abs(asymptotic_charge) > 1e-3:
        problems.append(f"r*V at r_max = {asymptotic_charge}, e-e screening missing")
    if int(header["ips"]) != 0:
        problems.append("ion-sphere cut-off active")
elif abs(asymptotic_charge + (Z_IRON - num_electrons + 1)) > 1e-3:
    problems.append(f"free-ion r*V at r_max = {asymptotic_charge}")

status = "ok" if not problems else "FAILED: " + "; ".join(problems)
manifest = {
    "ion": ion_name, "electrons": num_electrons, "model": model,
    "configurations": configs, "n_levels": n_levels,
    "mu_requested": mu, "mu_used": 1 / debye_length if debye_length > 0 and plasma_model == 1 else 0.0,
    "mu_steps": steps, "plasma_model_mps": plasma_model,
    "debye_mode": 1 if mu > 0 else 0, "ee_screen": 1 if mu > 0 else 0,
    "hxs_iter0": HXS_ITER0 if mu > 0 else 0, "density_A3": DENSITY_A3,
    "temperature_eV": temperature_for_mu(mu) if mu > 0 else None,
    "asymptotic_charge": asymptotic_charge,
    "E0_eV": energy_ev, "energy_unit": "eV (total energy of the lowest level)",
    "commit": git_commit(), "fac_extension": fac.__file__, "status": status,
}
with open("manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)
print(json.dumps({k: v for k, v in manifest.items() if k != "configurations"}, indent=2), flush=True)
if problems:
    sys.exit("calculation failed: " + status)
