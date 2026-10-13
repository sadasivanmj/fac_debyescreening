#!/usr/bin/env python3
"""
Nuclear-only Debye-Huckel screening in FAC (pfac), Chen et al. PoP 25, 072120 (2018) model "b".

Requires the one-hunk sign patch to faclib/orbital.c (fac-debye-nuclear-screening.patch).

FAC has no direct lambda_D entry point.  The Debye length actually used by the code is
    dps = sqrt(tps / (4*pi*nps*(1+ups)))                      [orbital.c:3380]
with  nps  given in A^-3  (converted internally by RBOHR**3)
      tps  given in eV    (converted internally by HARTREE_EV)
      ups  = 0 selects the electron-only Debye length.
This module inverts that relation.

IMPORTANT: use lambda_D = 1e6 a0 as the free-ion control, NOT a run without
PlasmaScreen().  Switching mps>=0 on changes how FAC builds the direct potential U
(the mps<0 path uses the parametrised N1 = N - |ihx| form, VT*r -> -(Z-N1); the
mps>=0 path builds the full Hartree-Slater potential, VT*r -> -(Z-N)).  That is a
constant, lambda-independent offset (~4e-3 a.u. in E0, ~1e-2 a.u. in excitation
energies for He-like C) which cancels exactly only if the reference run is taken in
the same mode.
"""
import math

RBOHR = 0.52917721067      # consts.h
HARTREE_EV = 27.211386018  # consts.h

def debye_params(lambda_d_a0, nps_per_A3=1.0e-2):
    """Return (nps [A^-3], tps [eV]) giving the requested lambda_D in bohr, with ups=0."""
    n_au = nps_per_A3 * RBOHR**3
    tps_au = lambda_d_a0**2 * 4.0 * math.pi * n_au
    return nps_per_A3, tps_au * HARTREE_EV

def plasma_screen_args(lambda_d_a0, z_net, nps_per_A3=1.0e-2):
    """Arguments for PlasmaScreen(zps, nps, tps, m, ups, vxf) with m=1 (Debye)."""
    nps, tps = debye_params(lambda_d_a0, nps_per_A3)
    return dict(zps=float(z_net), nps=nps, tps=tps, m=1, ups=0.0, vxf=1)

FREE_LAMBDA = 1.0e6   # a0; use this for the lambda -> infinity control

if __name__ == "__main__":
    import sys
    lam = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    z_net = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    p = plasma_screen_args(lam, z_net)
    print("PlasmaScreen(%.10g, %.10g, %.14g, %d, %.1f, %d)"
          % (p['zps'], p['nps'], p['tps'], p['m'], p['ups'], p['vxf']))
