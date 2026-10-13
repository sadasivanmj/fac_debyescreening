import math
from pfac import fac

# Screening parameters matching the project examples
NPS = 1.0e-2          # density in Angstrom^-3
LAMBDA_D = 10.0       # Debye length in Bohr radii
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018

density_au = NPS * RBOHR_ANG**3
tps = (
    LAMBDA_D**2
    * 4.0 * math.pi
    * density_au
    * HARTREE_EV
)

fac.SetAtom("Fe")

# Electron-nucleus screening ON; electron-electron screening ON (only change vs fe_calc/fe_iii_diag.py)
fac.SetOption("orbital:debye_mode", 1)
fac.SetOption("radial:ee_screen", 1)

# Fe III has 24 electrons
fac.Config("1s2 2s2 2p6 3s2 3p6 3d6", group="g1")

# Python binding signature:
# PlasmaScreen(zps, nps, tps, m, ups, vxf)
fac.PlasmaScreen(2.0, NPS, tps, 1, 0.0, 1)

print("CHECK 1: Before ConfigEnergy(0)", flush=True)
fac.ConfigEnergy(0)

print("CHECK 2: Before OptimizeRadial", flush=True)
fac.OptimizeRadial(["g1"])

print("CHECK 3: Before ConfigEnergy(1)", flush=True)
fac.ConfigEnergy(1)

print("CHECK 4: Before Structure", flush=True)
fac.Structure("fe_calc/debug_FeIII.lev.b", ["g1"])

print("CHECK 5: Before MemENTable", flush=True)
fac.MemENTable("fe_calc/debug_FeIII.lev.b")
fac.PrintTable(
    "fe_calc/debug_FeIII.lev.b",
    "fe_calc/debug_FeIII.lev",
    1,
)

print("Fe III diagnostic calculation completed.", flush=True)
