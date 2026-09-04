#!/usr/bin/env python3
"""
Minimal He-like C4+ nuclear-only Debye test.   Usage:  python3 test_debye_c4.py 20
Prints ground-state total energy, level energies, and excitation energies (a.u.).
Set element/Z_net at the top for Ne8+, Ar16+, Kr34+.
"""
import sys, re
from pfac import fac
from debye_fac import plasma_screen_args, FREE_LAMBDA, HARTREE_EV

ELEM, ZNET = 'C', 4.0
lam = sys.argv[1] if len(sys.argv) > 1 else 'free'
lam_d = FREE_LAMBDA if lam in ('free', 'inf') else float(lam)
tag = 'c_%s' % lam

fac.SetAtom(ELEM)
fac.SetOption('orbital:debye_mode', 1)      # 1 = screen Z(r) only  -> V_en = -Z e^{-r/lam}/r
fac.Config('g1', '1s2')
fac.Config('g2', '1s1 2s1', '1s1 2p1')
fac.Config('g3', '1s1 3s1', '1s1 3p1', '1s1 3d1')

p = plasma_screen_args(lam_d, ZNET)
fac.PlasmaScreen(p['zps'], p['nps'], p['tps'], p['m'], p['ups'], p['vxf'])

fac.ConfigEnergy(0)
fac.OptimizeRadial(['g1'])
fac.ConfigEnergy(1)
fac.GetPotential('%s.pot' % tag)            # cols 2,3,18 = r, Z(r), ZPS(r)
fac.Structure('%s.lev.b' % tag, ['g1', 'g2', 'g3'])
fac.MemENTable('%s.lev.b' % tag)
fac.PrintTable('%s.lev.b' % tag, '%s.lev' % tag, 1)

for line in open('%s.lev' % tag):
    if line.startswith('E0'):
        print('lambda_D = %-10s  E(ground) = %.8f a.u.' % (lam_d, float(line.split(',')[1]) / HARTREE_EV))
    f = line.split()
    if len(f) > 6 and re.match(r'^\d+$', f[0]):
        print('  %3s  2J=%s  vnl=%s  E_exc = %12.6f a.u.  %s'
              % (f[0], f[5], f[4], float(f[2]) / HARTREE_EV, f[-1]))
