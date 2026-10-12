import math
from pfac import fac

# Minimal nuclear-only Debye screening test: He-like C4+
LAMBDA_D = 10.0       # Debye length in Bohr radii
NPS = 1.0e-2          # density in Angstrom^-3
RBOHR_ANG = 0.52917721067
HARTREE_EV = 27.211386018

density_au = NPS * RBOHR_ANG**3
tps = LAMBDA_D**2 * 4.0 * math.pi * density_au * HARTREE_EV

fac.SetAtom("C")

fac.SetOption("orbital:debye_mode", 1)
fac.SetOption("radial:ee_screen",0)

fac.Config("g1", "1s2")

# Python wrapper order: zps, nps, tps, m, ups, vxf
fac.PlasmaScreen(4.0, NPS, tps, 1, 0.0, 1)

print("CHECK 1: Before ConfigEnergy(0)", flush=True)
fac.ConfigEnergy(0)

print("CHECK 2: Before OptimizeRadial", flush=True)
fac.OptimizeRadial(["g1"])

print("CHECK 3: Before Structure", flush=True)
fac.ConfigEnergy(1)
fac.Structure("fe_calc/C4_test.lev.b", ["g1"])

fac.MemENTable("fe_calc/C4_test.lev.b")
fac.PrintTable("fe_calc/C4_test.lev.b", "fe_calc/C4_test.lev", 1)

print("Minimal screened C4+ calculation completed.", flush=True)
