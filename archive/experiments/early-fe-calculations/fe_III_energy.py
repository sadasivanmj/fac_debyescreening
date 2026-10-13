
from pathlib import Path
from pfac import fac

# --------------------------------------------------
# Fe III: reference configurations from Deprince et al.
# Table 1. Unscreened calculation: mu = 0.
# --------------------------------------------------

OUT = Path("fe_calc")
OUT.mkdir(exist_ok=True)

GROUP = "FeIII_ref"
LEVEL_BIN = OUT / "FeIII_mu0.lev.b"
LEVEL_TXT = OUT / "FeIII_mu0.lev"

fac.SetAtom("Fe")

# Fe III has 24 electrons.
# [n l] in the paper denotes a hole in that subshell.
configs = [
    # Ground reference
    "1s2 2s2 2p6 3s2 3p6 3d6",

    # [3p] references
    "1s2 2s2 2p6 3s2 3p5 3d7",
    "1s2 2s2 2p6 3s2 3p5 3d6 4s1",
    "1s2 2s2 2p6 3s2 3p5 3d5 4s2",

    # [2p] references
    "1s2 2s2 2p5 3s2 3p6 3d7",
    "1s2 2s2 2p5 3s2 3p6 3d6 4s1",
    "1s2 2s2 2p5 3s2 3p6 3d5 4s2",

    # [1s] references
    "1s1 2s2 2p6 3s2 3p6 3d7",
    "1s1 2s2 2p6 3s2 3p6 3d6 4s1",
    "1s1 2s2 2p6 3s2 3p6 3d5 4s2",
]

print("Configuring Fe III references...")
for cfg in configs:
    fac.Config(cfg, group=GROUP)

print("Optimizing radial orbitals...")
fac.ConfigEnergy(0)
fac.OptimizeRadial([GROUP])

print("Calculating energy levels...")
fac.ConfigEnergy(1)
fac.Structure(str(LEVEL_BIN), [GROUP])

fac.MemENTable(str(LEVEL_BIN))

# Verbose output format: energies in eV.
fac.PrintTable(str(LEVEL_BIN), str(LEVEL_TXT), 1)

# Read the total ground-state energy reference from the table header.
e0 = None
nele = None
nlev = None

with LEVEL_TXT.open() as f:
    for line in f:
        line = line.strip()

        if line.startswith("E0"):
            e0 = float(line.split(",")[1])

        elif line.startswith("NELE"):
            nele = int(line.split("=")[1])

        elif line.startswith("NLEV"):
            nlev = int(line.split("=")[1])

if e0 is None:
    raise RuntimeError(f"Could not read E0 from {LEVEL_TXT}")

print("\n========== Fe III RESULTS ==========")
print(f"Electron count          : {nele}")
print(f"Number of levels         : {nlev}")
print(f"Ground-state E0          : {e0:.8f} eV")
print(f"Binary output            : {LEVEL_BIN}")
print(f"Text output              : {LEVEL_TXT}")
print("====================================")

