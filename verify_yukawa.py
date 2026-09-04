#!/usr/bin/env python3
"""Potential-level check, independent of any atomic energy.
Reads a GetPotential() dump and prints V_screened/V_Coulomb, which must equal exp(-r/lambda_D).
Usage:  python3 verify_yukawa.py c_20.pot 20"""
import sys, math
fn, lam = sys.argv[1], float(sys.argv[2])
rows = [[float(x) for x in L.split()] for L in open(fn)
        if not L.startswith('#') and len(L.split()) == 18]
print('%9s %13s %13s %11s %11s %10s' % ('r(a0)', 'V_Coulomb', 'V_screened', 'ratio', 'exp(-r/lam)', 'rel.err'))
for target in (1.0, 5.0, 10.0, 20.0):
    i = min(range(len(rows)), key=lambda k: abs(rows[k][1] - target))
    r, Z, ZPS = rows[i][1], rows[i][2], rows[i][17]
    vc, vs = -Z / r, (-Z + ZPS) / r          # VT_nuclear = (-Z + ZPS)/r
    ref = math.exp(-r / lam)
    print('%9.4f %13.7f %13.7f %11.7f %11.7f %10.1e' % (r, vc, vs, vs / vc, ref, vs / vc / ref - 1))
