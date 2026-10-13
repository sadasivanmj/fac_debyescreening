#!/usr/bin/env python3
"""
quick_verify_stage2_v3.py -- TEST 1 (large-Debye limit, self-consistent)
and TEST 2 (e-e screening isolation), rebuilt on the sfac-driver pattern
from chen2018_reproduce.py (sf_input / parse_lev / run_one).

CHANGES FROM v2:
1. Energy extraction matches parse_lev() exactly (absolute ground energy).
2. Units: .lev energies converted from eV to a.u.
3. TEST 1 compares two large, finite lambda_D values on the SAME mps=1 path.
4. TEST 2's case B (e-e only) uses mps=0, not mps=1, properly isolating e-e.
5. Pre-check confirms mps=0 and mps=1 match closely when screening is off.

FIXED IN THIS VERSION:
Test 2's failure condition now correctly evaluates (E_b - E_ref0) rather 
than (E_b - E_a) to ensure the e-e screening actually lowered the energy 
relative to the unscreened vacuum state.
"""
import os
import csv
import math
import subprocess

OUTDIR = "/data/project/fac/quick_verification"
os.makedirs(OUTDIR, exist_ok=True)

SFAC = os.environ.get("SFAC", "sfac")
WORKDIR = "quick_stage2_run"
TIMEOUT_S = 300          # generous but hard; one run per case, no retry

RBOHR = 0.52917721067
HARTREE_EV = 27.211386018

Z = 6
NET_CHARGE = 4.0          # He-like C4+ : znet = Z - N = 6 - 2 = 4
NPS = 1.0e-2              # A^-3, arbitrary -- tps compensates

# TEST 1: two large, finite lambda_D values, same mps=1 path
LAM1 = 1.0e5
LAM2 = 1.0e6
TOL_REL_TEST1 = 0.05      # loose: this is a leading-order slope check

# TEST 2: one moderate lambda_D
LAM_MOD = 10.0
PRECHECK_TOL_AU = 1e-6    # mps=0 vs mps=1 with all screening off must match this tightly


def tps_for(lam, nps=NPS):
    """eV temperature giving dps == lam (bohr), ups=0."""
    return lam * lam * 4.0 * math.pi * (nps * RBOHR ** 3) * HARTREE_EV


def sf_input(mps, ee_screen, lam, zps=NET_CHARGE, nps=NPS):
    """Minimal ground-state-only sfac script. mps: 1 = nuclear Debye path,
    0 = plasma path active but nuclear branch off. ee_screen: 0/1."""
    s = "SetAtom('C')\n"
    s += "SetOption('orbital:debye_mode', 1)\n"
    s += "SetOption('radial:ee_screen', %d)\n" % ee_screen
    s += "Config('g1','1s2')\n"
    s += "PlasmaScreen(%.10g, %.10g, %.14g, %d, 0.0, 1)\n" % (zps, nps, tps_for(lam, nps), mps)
    s += "ConfigEnergy(0)\nOptimizeRadial(['g1'])\nConfigEnergy(1)\n"
    s += "GetPotential('tag.pot')\n"
    s += "Structure('tag.lev.b',['g1'])\n"
    s += "MemENTable('tag.lev.b')\n"
    s += "PrintTable('tag.lev.b','tag.lev',1)\n"
    return s


def parse_lev(fn):
    """Identical logic to chen2018_reproduce.py:parse_lev."""
    E0 = None
    for L in open(fn):
        if L.startswith("E0"):
            E0 = float(L.split(",")[1]) / HARTREE_EV
            break
    return E0


def run_case(tag, mps, ee_screen, lam, zps=NET_CHARGE, nps=NPS):
    os.makedirs(WORKDIR, exist_ok=True)
    sf_path = os.path.join(WORKDIR, tag + ".sf")
    script = sf_input(mps, ee_screen, lam, zps, nps).replace("tag.", tag + ".")
    open(sf_path, "w").write(script)
    try:
        r = subprocess.run([SFAC, tag + ".sf"], cwd=WORKDIR,
                           capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT", f"{tag}: exceeded {TIMEOUT_S}s, not retried"
    if r.returncode != 0:
        msg = [l for l in r.stdout.split("\n") if "Error" in l or "Max iteration" in l]
        return None, "FAIL", (msg[0].strip() if msg else f"{tag}: sfac returned {r.returncode}")
    lev_path = os.path.join(WORKDIR, tag + ".lev")
    if not os.path.exists(lev_path):
        return None, "FAIL", f"{tag}: no .lev file produced"
    E0 = parse_lev(lev_path)
    if E0 is None:
        with open(lev_path) as f:
            head = "".join(f.readlines()[:15])
        return None, "FAIL", (f"{tag}: no 'E0' header line found in {tag}.lev.\n"
                              f"First 15 lines:\n{head}")
    return E0, "PASS", ""


def analytic_leading_nuclear_shift(lam):
    """dE_leading = N*Z*kappa = N*Z/lam."""
    N = 2
    return N * Z / lam


def test1_large_lambda():
    e1, s1, m1 = run_case("t1_lam1", mps=1, ee_screen=0, lam=LAM1)
    if s1 != "PASS":
        return s1, None, m1
    e2, s2, m2 = run_case("t1_lam2", mps=1, ee_screen=0, lam=LAM2)
    if s2 != "PASS":
        return s2, None, m2

    delta_e = e2 - e1
    analytic = analytic_leading_nuclear_shift(LAM2) - analytic_leading_nuclear_shift(LAM1)
    residual = delta_e - analytic
    rel_err = abs(residual) / abs(analytic) if analytic != 0 else float("inf")
    status = "PASS" if rel_err < TOL_REL_TEST1 else "FAIL"
    row = [e1, e2, delta_e, analytic, residual, rel_err, status]
    return status, row, ""


def test2_ee_isolation():
    # Pre-check: does mps=0 reproduce mps=1 to tight tolerance when screening is off?
    e_ref1, s_r1, m_r1 = run_case("t2_precheck_mps1", mps=1, ee_screen=0, lam=LAM_MOD)
    if s_r1 != "PASS":
        return s_r1, None, m_r1
    e_ref0, s_r0, m_r0 = run_case("t2_precheck_mps0", mps=0, ee_screen=0, lam=LAM_MOD)
    if s_r0 != "PASS":
        return s_r0, None, m_r0
    precheck_diff = abs(e_ref0 - e_ref1)
    if precheck_diff > PRECHECK_TOL_AU:
        return "FAIL", None, (
            "OBSERVATION: mps=0 and mps=1 give different energies "
            f"({e_ref0:.10f} vs {e_ref1:.10f} a.u., both nominally unscreened)\n"
            "LIKELY CAUSE: the mps=0 and mps=1 code paths use a different\n"
            "  zeroth-order Hamiltonian convention (same class of issue as the\n"
            "  documented off-vs-mps=1 offset), so case B below would not be a\n"
            "  clean isolation even with EEScreenMu patched.\n"
            "NEXT CHECK: compare pot->zps and the asymptotic VT*r values FAC\n"
            "  actually uses for mps=0 vs mps=1 at the same nominal conditions."
        )

    e_a, s_a, m_a = run_case("t2_a_nuclear_only", mps=1, ee_screen=0, lam=LAM_MOD)
    if s_a != "PASS":
        return s_a, None, m_a
    e_b, s_b, m_b = run_case("t2_b_ee_only", mps=0, ee_screen=1, lam=LAM_MOD)
    if s_b != "PASS":
        return s_b, None, m_b
    e_c, s_c, m_c = run_case("t2_c_both", mps=1, ee_screen=1, lam=LAM_MOD)
    if s_c != "PASS":
        return s_c, None, m_c

    b_minus_ref = e_b - e_ref0
    c_minus_a = e_c - e_a
    c_minus_b = e_c - e_b
    row = [e_a, e_b, e_c, b_minus_ref, c_minus_a, c_minus_b]

    # Check if e-e screening properly lowered the energy vs unscreened state
    if b_minus_ref >= 0:
        return "FAIL", row, (
            "OBSERVATION: (B - Ref) >= 0, i.e. isolated e-e screening did not lower\n"
            "  the total energy below the unscreened baseline.\n"
            "LIKELY CAUSE: either mu=0 in case B still (EEScreenMu patch not\n"
            "  actually built/installed into the binary sfac is calling), or a\n"
            "  sign error in GetYkDebye/DensityToSZ.\n"
            "NEXT CHECK: add `fprintf(stderr,\"mu=%g\\n\",_mu);` at the top of\n"
            "  GetYk right after `double _mu = EEScreenMu();`, rebuild, and\n"
            "  confirm case B actually prints mu=0.1 (1/lambda_D), not mu=0."
        )
    return "PASS", row, ""


def main():
    t1_status, t1_row, t1_msg = test1_large_lambda()
    t2_status, t2_row, t2_msg = test2_ee_isolation()

    with open(os.path.join(OUTDIR, "large_lambda.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_lam1", "E_lam2", "Delta_E", "analytic_leading_term",
                    "residual", "relative_error", "status"])
        if t1_row:
            w.writerow(t1_row)
        else:
            w.writerow([f"# {t1_status}: {t1_msg}"])

    with open(os.path.join(OUTDIR, "ee_isolation.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_A_nuclear_only", "E_B_ee_only", "E_C_both",
                    "B_minus_Ref", "C_minus_A", "C_minus_B"])
        if t2_row:
            w.writerow(t2_row)
        else:
            w.writerow([f"# {t2_status}: {t2_msg}"])

    overall = "PRELIMINARY PASS" if t1_status == "PASS" and t2_status == "PASS" else "PRELIMINARY FAIL"
    lines = [
        "FAST VERIFICATION",
        f"Large-Debye limit: {t1_status}",
        f"E-E isolation:      {t2_status}",
        "",
        f"Overall: {overall}",
        "",
        "Preliminary verification only.",
    ]
    if t1_msg:
        lines += ["", f"TEST 1 note: {t1_msg}"]
    if t2_msg:
        lines += ["", f"TEST 2 note: {t2_msg}"]

    with open(os.path.join(OUTDIR, "summary_stage2.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
