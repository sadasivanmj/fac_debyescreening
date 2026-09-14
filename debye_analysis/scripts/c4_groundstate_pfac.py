# Ground state (1s2 1S0) energy of He-like C4+ (Z=6, N=2) vs Debye length.
# Uses the patched FAC options confirmed in this project:
#   orbital:debye_mode  -> electron-nucleus Debye screening of the mean potential
#   radial:ee_screen    -> Yukawa screening of the two-electron Slater integrals
# Benchmark: Chen, Ma, Hu & Wang, Phys. Plasmas 25, 072120 (2018), Table I,
# FAC column, for C4+.
#
# Debye length definition (matches this project's convention and the
# InitializePS() formula in the patched orbital.c):
#   lambda_D [a0] = sqrt( T[Hartree] / (4 pi n_e[a0^-3] (1 + Z*)) )
# Z* defaults to the ion net charge (Z - N_bound) if not set explicitly.
# This is the same default the patched code falls back to when ups < 0.
#
# n_e is fixed at a convenient value. Only lambda_D affects the physics here,
# so the exact n_e is not important. T is solved for each target lambda_D.
#
# Open items, not verified by code inspection alone:
#  - The exact text layout of the "E0 = ..." line in the ASCII .lev file is
#    assumed to match the FAC 1.1.5 format seen in the manual demo output.
#    Check this against your installed FAC version if parsing fails.
#  - This script was not executed here (no FAC build available in this
#    environment). Run it in your own FAC install and check the printed
#    values before trusting them.

from pfac.fac import *
import math

HARTREE_EV = 27.211386245988
RBOHR_CM = 0.529177210903e-8   # Bohr radius, cm

Z = 6                 # nuclear charge, carbon
N_BOUND = 2           # He-like: 2 bound electrons
ZSTAR = Z - N_BOUND   # net ion charge, used as (1 + Z*) in the Debye formula
NE_CM3 = 1.0e21       # fixed electron density, cm^-3 (arbitrary choice)

# Debye lengths from Chen Table I for C4+, in a0. None = free ion (unscreened).
DEBYE_LENGTHS = [None, 100.0, 20.0, 12.5, 10.0, 5.0, 10.0/3.0, 2.5, 2.0]

# Chen Table I, FAC column, ground-state energy (a.u.) for C4+, same order.
CHEN_FAC_C4 = [-32.4771, -32.3768, -31.9081, -31.5601,
               -31.3298, -30.1973, -29.0950, -28.0242, -26.9818]


def temperature_ev(lam_d_a0, ne_cm3, zstar):
    """Solve for T [eV] that gives the target Debye length at fixed density."""
    ne_a0 = ne_cm3 * RBOHR_CM**3
    t_hartree = (lam_d_a0**2) * 4.0 * math.pi * ne_a0 * (1.0 + zstar)
    return t_hartree * HARTREE_EV


def read_e0(fname):
    """Read the absolute reference energy E0 (Hartree) from an ASCII .lev file."""
    with open(fname) as f:
        for line in f:
            if line.strip().startswith('E0'):
                # expected form: "E0 = <ilev0>, <energy>"
                return float(line.split(',')[1])
    raise RuntimeError('E0 line not found in ' + fname)


def ground_state_energy(lam_d_a0, tag):
    group = 'g_' + tag
    lev_b = 'c4_' + tag + '.lev.b'
    lev_a = 'c4_' + tag + '.lev'

    SetAtom('C')
    Config(group, '1s2')

    if lam_d_a0 is None:
        # explicit off: clears any screening left over from a prior call
        SetOption('orbital:debye_mode', '', 0, 0.0)
        SetOption('radial:ee_screen', '', 0, 0.0)
        PlasmaScreen(-1, 0, 0.0, 0.0, 0.0, 0.0)
    else:
        t_ev = temperature_ev(lam_d_a0, NE_CM3, ZSTAR)
        SetOption('orbital:debye_mode', '', 1, 0.0)
        SetOption('radial:ee_screen', '', 1, 0.0)
        # PlasmaScreen(m, vxf, zps, nps, tps, ups)
        # m=1: Debye model. zps=0 lets FAC derive the net charge from the
        # configuration. ups=-1 lets FAC default Z* to that same net charge.
        PlasmaScreen(1, 0, 0.0, NE_CM3, t_ev, -1.0)

    ConfigEnergy(0)
    OptimizeRadial([group])
    ConfigEnergy(1)
    Structure(lev_b, [group])
    PrintTable(lev_b, lev_a, 0)
    return read_e0(lev_a)


def main():
    print('GROUND STATE 1s2 1S0, He-like C4+')
    print('Debye length in a0. Energy in a.u. Diff = FAC - Chen.')
    header = '{:>8} {:>14} {:>14} {:>10} {:>8}'.format(
        'lam_D', 'FAC', 'Chen_FAC', 'Diff', '%Diff')
    print(header)
    for lam_d, chen in zip(DEBYE_LENGTHS, CHEN_FAC_C4):
        tag = 'inf' if lam_d is None else str(lam_d).replace('.', 'p')
        label = 'inf' if lam_d is None else '{:g}'.format(lam_d)
        e = ground_state_energy(lam_d, tag)
        diff = e - chen
        pct = 100.0 * abs(diff) / abs(chen)
        row = '{:>8} {:14.4f} {:14.4f} {:10.4f} {:8.3f}'.format(
            label, e, chen, diff, pct)
        print(row)


if __name__ == '__main__':
    main()
