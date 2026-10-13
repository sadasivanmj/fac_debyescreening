#!/usr/bin/env python3
"""
quick_verify_stage2.py -- TEST 1 (large-Debye free-ion limit) and
TEST 2 (e-e screening isolation), for the patched FAC in your build tree.

HONESTY NOTE (read before trusting any PASS in the output):
This script actually calls `pfac`, which requires a *compiled FAC binary*
with the full faclib source tree (structure.c, config.c, coulomb.c, the
angular-momentum/Racah modules, etc.) and the Structure()/CI machinery that
solves for real bound-state energies. Only two source files (orbital.c,
radial.c) were made available for review, there is no build system, and
this sandbox has no gfortran -- so this script has NOT been executed
end-to-end against a working FAC here, and no energy numbers below were
invented to fill in the expected CSV columns. Run it in the environment
where you actually build FAC. If pfac cannot be imported, the script
reports SKIP rather than guessing.
"""
import csv
import math
import os
import subprocess
import sys
import textwrap

OUTDIR = "/home/claude/quick_verification"
os.makedirs(OUTDIR, exist_ok=True)

TIMEOUT_S = 120          # hard per-run timeout, one FAC run per case, no retry
Z = 6
N = 2                    # He-like C4+
LAMBDA_LARGE = 1.0e6      # a0, "asymptotic" test point for TEST 1
LAMBDA_MODERATE = 10.0    # a0, for TEST 2
TOL_REL_TEST1 = 1e-2      # residual/leading_term must be < this to PASS
                          # (loose on purpose: this is a leading-order check,
                          #  not an exact one -- see analytic derivation below)


# ---------------------------------------------------------------------------
# Analytic expectations, derived (not assumed) from the sign convention
# confirmed against your actual orbital.c in the prior verification pass:
#   V_en_screened(r) = -Z exp(-r/lambda_D) / r
#   V_ee_screened(r12) = exp(-r12/lambda_D) / r12   (mu = 1/lambda_D)
# ---------------------------------------------------------------------------

def analytic_leading_nuclear_shift(Z, N, lambda_D):
    """TEST 1. Perturbation added to H by turning on nuclear screening:
       dV(r) = V_screened(r) - V_unscreened(r)
             = -Z*exp(-r/lambda_D)/r - (-Z/r)
             =  Z*(1 - exp(-r/lambda_D))/r
    Small-kappa expansion (kappa = 1/lambda_D):
       (1 - exp(-kappa r))/r = kappa - kappa^2 r/2 + O(kappa^3 r^2)
    First-order PT, summed over N electrons:
       dE = N*Z*kappa - (Z*kappa^2/2)*sum<r_i> + O(kappa^3)
    Leading term as lambda_D -> infinity (kappa -> 0):
       dE_leading = N*Z*kappa = N*Z/lambda_D   (POSITIVE)
    Sign check: dV(r) > 0 pointwise (screening always weakens the bare
    attraction), so its expectation value over any state is positive ->
    the screened ion is LESS bound than the free ion, i.e. E_screened is
    algebraically LARGER (less negative) than E_unscreened. This is the
    opposite sign from +/-NZ/lambda_D asserted without derivation -- it is
    derived here as strictly +N*Z/lambda_D, never negative, regardless of Z.
    """
    kappa = 1.0 / lambda_D
    return N * Z * kappa


def analytic_ee_sign():
    """TEST 2, sign only (the magnitude needs orbitals FAC would compute).
       dV_ee(r12) = exp(-mu r12)/r12 - 1/r12 = (exp(-mu r12) - 1)/r12 < 0
       for all r12 > 0, mu > 0. First-order PT over electron pairs gives a
       NEGATIVE expectation value: turning on e-e screening alone lowers
       (more negative) the total energy relative to case A.
       => expect (B - A) < 0.
    Nuclear-only screening does the opposite (see analytic_leading_nuclear_shift):
       => expect the nuclear contribution alone to be > 0.
    Combined case C mixes a positive (nuclear) and negative (e-e) contribution;
    the NET sign of (C - A) is not determined by this argument alone and must
    come from the actual calculation. (C - B) is expected to carry the same
    sign as the nuclear-only contribution (positive), since it isolates adding
    nuclear screening on top of an already e-e-screened Hamiltonian.
    """
    return {
        "B_minus_A": "expected < 0 (screening reduces e-e repulsion -> more bound)",
        "nuclear_only_contribution": "expected > 0 (see TEST 1 derivation)",
        "C_minus_A": "sign not determined analytically -- competing terms, needs the run",
        "C_minus_B": "expected > 0 (same sign as the isolated nuclear contribution)",
    }


# ---------------------------------------------------------------------------
# FAC drivers. These generate a small pfac script and run it as a subprocess
# with a hard timeout. They are written against the API surface documented
# in your own forensic report (Section 13.2: SetOption/PlasmaScreen/OptimizeRadial)
# plus the standard pfac Structure()/MemENTable() pattern for reading a level
# energy. The exact config-building calls (Config(), etc.) depend on how your
# build reads/writes level tables -- adjust GET_ENERGY_PY below to match your
# own working scripts if the placeholder does not match your pfac version.
# ---------------------------------------------------------------------------

RUN_TEMPLATE = textwrap.dedent("""
    import sys
    try:
        from pfac import fac
    except ImportError:
        print("SKIP: pfac not importable in this environment")
        sys.exit(3)

    fac.SetAtom('C')
    fac.Config('He-like', '1s2')

    {screening_lines}

    fac.OptimizeRadial(['He-like'])
    fac.Structure('stage2.lev.b', ['He-like'])
    fac.MemENTable('stage2.lev.b')
    fac.PrintTable('stage2.lev.b', 'stage2.lev', 1)

    # NOTE: adjust this if your level-table column layout differs.
    with open('stage2.lev') as f:
        lines = [l for l in f if l.strip() and not l.startswith('#')]
    e0 = float(lines[0].split()[2])
    print("ENERGY=%.12f" % e0)
""")

NUCLEAR_OFF = ""  # mps defaults to 0 / unscreened
NUCLEAR_LARGE_LAMBDA = f"""
    fac.SetOption('orbital:debye_mode', 1)
    fac.PlasmaScreen({Z - N}, 1e-30, {LAMBDA_LARGE**2 * 4*math.pi*1e-30*27.211386018:.6e}, 1, 0.0, 1)
"""
NUCLEAR_MODERATE = f"""
    fac.SetOption('orbital:debye_mode', 1)
    fac.PlasmaScreen({Z - N}, 1e-4, {LAMBDA_MODERATE**2 * 4*math.pi*1e-4*27.211386018:.6e}, 1, 0.0, 1)
"""
EE_ONLY_MODERATE = f"""
    fac.SetOption('radial:ee_screen', 1)
    fac.PlasmaScreen(0, 1e-4, {LAMBDA_MODERATE**2 * 4*math.pi*1e-4*27.211386018:.6e}, 0, 0.0, 1)
    # Requires the EEScreenMu decoupling in ee_screen_decouple.patch --
    # without it, mps must be 1 for e-e screening to activate at all,
    # which would silently also screen the nucleus.
"""
BOTH_MODERATE = f"""
    fac.SetOption('orbital:debye_mode', 1)
    fac.SetOption('radial:ee_screen', 1)
    fac.PlasmaScreen({Z - N}, 1e-4, {LAMBDA_MODERATE**2 * 4*math.pi*1e-4*27.211386018:.6e}, 1, 0.0, 1)
"""


def run_fac_case(screening_lines, label):
    script = RUN_TEMPLATE.format(screening_lines=screening_lines)
    try:
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True, text=True, timeout=TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", None, f"{label}: exceeded {TIMEOUT_S}s, not retried")

    if "SKIP: pfac not importable" in result.stdout:
        return ("SKIP", None, f"{label}: pfac not available in this environment")
    if result.returncode != 0:
        return ("FAIL", None, f"{label}: nonzero exit -- stderr: {result.stderr[-500:]}")

    for line in result.stdout.splitlines():
        if line.startswith("ENERGY="):
            return ("PASS", float(line.split("=", 1)[1]), "")
    return ("FAIL", None, f"{label}: no ENERGY= line in output")


def test1_large_lambda():
    rows = []
    status_a, e_unscreened, msg_a = run_fac_case(NUCLEAR_OFF, "A (screening off)")
    if status_a != "PASS":
        return status_a, rows, msg_a
    status_b, e_large, msg_b = run_fac_case(NUCLEAR_LARGE_LAMBDA, f"B (lambda_D={LAMBDA_LARGE:g})")
    if status_b != "PASS":
        return status_b, rows, msg_b

    delta_e = e_large - e_unscreened
    leading = analytic_leading_nuclear_shift(Z, N, LAMBDA_LARGE)
    residual = delta_e - leading
    rel_err = abs(residual) / abs(leading) if leading != 0 else float("inf")
    status = "PASS" if rel_err < TOL_REL_TEST1 else "FAIL"
    rows.append([e_unscreened, e_large, delta_e, leading, residual, rel_err, status])
    return status, rows, ""


def test2_ee_isolation():
    rows = []
    status_a, e_a, msg_a = run_fac_case(NUCLEAR_OFF, "A (no screening)")
    if status_a != "PASS":
        return status_a, rows, msg_a
    status_b, e_b, msg_b = run_fac_case(EE_ONLY_MODERATE, "B (e-e only)")
    if status_b != "PASS":
        return status_b, rows, msg_b
    status_c, e_c, msg_c = run_fac_case(BOTH_MODERATE, "C (both)")
    if status_c != "PASS":
        return status_c, rows, msg_c

    rows.append([e_a, e_b, e_c, e_b - e_a, e_c - e_a, e_c - e_b])
    sign_check = (e_b - e_a) < 0
    status = "PASS" if sign_check else "FAIL"
    return status, rows, ("" if sign_check else
        "OBSERVATION: (B-A) >= 0\n"
        "LIKELY CAUSE: e-e screening not actually active in case B (check the\n"
        "  EEScreenMu decoupling was applied and mu>0), or a sign error remains\n"
        "  in GetYkDebye/DensityToSZ.\n"
        "NEXT CHECK: dump V_ee at fixed r1,r2 for case B and compare against\n"
        "  exp(-r12/lambda_D)/r12 the way the nuclear potential was checked.")


def main():
    t1_status, t1_rows, t1_msg = test1_large_lambda()
    t2_status, t2_rows, t2_msg = test2_ee_isolation()

    with open(os.path.join(OUTDIR, "large_lambda.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_unscreened", "E_large_lambda", "Delta_E",
                    "analytic_leading_term", "residual", "relative_error", "status"])
        for r in t1_rows:
            w.writerow(r)
        if not t1_rows:
            w.writerow([f"# {t1_status}: {t1_msg}"])

    with open(os.path.join(OUTDIR, "ee_isolation.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_A_no_screen", "E_B_ee_only", "E_C_both",
                    "B_minus_A", "C_minus_A", "C_minus_B"])
        for r in t2_rows:
            w.writerow(r)
        if not t2_rows:
            w.writerow([f"# {t2_status}: {t2_msg}"])

    overall = "PRELIMINARY PASS" if t1_status == "PASS" and t2_status == "PASS" else \
              ("PRELIMINARY SKIP" if "SKIP" in (t1_status, t2_status) else "PRELIMINARY FAIL")

    summary = []
    summary.append("FAST VERIFICATION")
    summary.append(f"Large-Debye limit: {t1_status}")
    summary.append(f"E-E isolation:      {t2_status}")
    summary.append("")
    summary.append(f"Overall: {overall}")
    summary.append("")
    summary.append("Preliminary verification only.")
    if t1_msg:
        summary.append("")
        summary.append(f"TEST 1 note: {t1_msg}")
    if t2_msg:
        summary.append("")
        summary.append(f"TEST 2 note: {t2_msg}")

    with open(os.path.join(OUTDIR, "summary_stage2.txt"), "w") as f:
        f.write("\n".join(summary) + "\n")

    print("\n".join(summary))


if __name__ == "__main__":
    main()
