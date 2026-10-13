# SUPPLEMENT, not a FAC result.
# FAC cannot apply the Debye nuclear screening to a one-electron ion (see
# README.md), so this independent solver gives the 1s Dirac energy of a single
# electron in V(r) = -Z exp(-mu r)/r (point nucleus) as a cross-check on the
# Fe XXVI helper energy.  Atomic units; E is measured from the rest energy.
# Plain RK4 on a logarithmic grid, outward/inward integration, bisection.
#
# Usage: python supplement_fe_xxvi_dirac_yukawa.py

import math

C = 137.035999084          # speed of light (a.u.)
Z = 26.0
HARTREE_EV = 27.211386
KAPPA = -1                 # 1s1/2
N_STEPS = 40000


def deriv(r, y, e, mu):
    # P = r*g, Q = r*f; dy/dr
    p, q = y
    v = -Z * math.exp(-mu * r) / r
    return (-KAPPA * p / r + (2 * C + (e - v) / C) * q,
            KAPPA * q / r - (e - v) / C * p)


def integrate(r_from, r_to, y, e, mu):
    # RK4 in t = ln r, so dr = r dt; keeps the step proportional to r
    t0, t1 = math.log(r_from), math.log(r_to)
    h = (t1 - t0) / N_STEPS
    p, q = y
    t = t0
    for _ in range(N_STEPS):
        def f(t, p, q):
            r = math.exp(t)
            dp, dq = deriv(r, (p, q), e, mu)
            return r * dp, r * dq
        k1 = f(t, p, q)
        k2 = f(t + h / 2, p + h / 2 * k1[0], q + h / 2 * k1[1])
        k3 = f(t + h / 2, p + h / 2 * k2[0], q + h / 2 * k2[1])
        k4 = f(t + h, p + h * k3[0], q + h * k3[1])
        p += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        q += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t += h
        n = math.hypot(p, q)
        if n > 1e100:
            p, q = p / n, q / n      # keep the numbers finite (only the ratio matters)
    return p, q


def mismatch(e, mu, r_match=0.05, r_max=30.0):
    gamma = math.sqrt(KAPPA**2 - (Z / C)**2)
    r0 = 1e-7
    # small-r behaviour: P ~ r^gamma, Q/P = (kappa + gamma) C / Z
    p1, q1 = integrate(r0, r_match, (1.0, (KAPPA + gamma) * C / Z), e, mu)
    kap = math.sqrt(-e * (e + 2 * C**2)) / C        # decay constant outside the potential
    p2, q2 = integrate(r_max, r_match, (1.0, -kap * C / (e + 2 * C**2)), e, mu)
    # both solutions must be proportional at r_match; use the normalised cross product
    return (p1 * q2 - p2 * q1) / (math.hypot(p1, q1) * math.hypot(p2, q2))


def energy_ev(mu):
    e_c = C**2 * (math.sqrt(1 - (Z / C)**2) - 1)    # point-Coulomb 1s
    # the Yukawa screening only raises the energy, by at most ~ Z*mu
    lo, hi = e_c - 1.0, e_c + 1.5 * Z * mu + 1.0
    f_lo = mismatch(lo, mu)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f_mid = mismatch(mid, mu)
        if f_lo * f_mid <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi) * HARTREE_EV


if __name__ == "__main__":
    exact = C**2 * (math.sqrt(1 - (Z / C)**2) - 1) * HARTREE_EV
    print(f"mu=0 analytic point-nucleus 1s: {exact:.4f} eV")
    for mu in (0.0, 0.1, 0.25):
        print(f"mu={mu:<5} Dirac-Yukawa 1s: {energy_ev(mu):.4f} eV")
