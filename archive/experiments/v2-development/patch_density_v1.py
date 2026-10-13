#!/usr/bin/env python3
"""Patch faclib/orbital.c:DensityToSZ with a Debye (Yukawa) screened direct potential.

Isolated-worktree tool: refuses to run anywhere except /tmp/fac-debye-scf-test.

Math (FAC conventions, see REPORT): d = radial density (int d dr = N), the grid
has unit steps in rho with dr/drho = dr_drho, and z stores r*V_direct.
Coulomb:  z = int_0^r d dr' + r int_r^inf d/r' dr'
Yukawa:   z = e^{-mu r} int_0^r d sinhc(mu r') dr'
              + r sinhc(mu r) int_r^inf (d/r') e^{-mu r'} dr'
with mu = EEScreenMu() (1/lambda_D when radial:ee_screen=1 and Debye mode,
0 otherwise).  For mu == 0 the original statements are kept verbatim.

Steps: verify location and the exact original function text (sha256), refuse
a second application, back up orbital.c, apply exact-string replacements
inside the function boundary, re-read and verify, write patch_density.diff.
"""

import difflib
import hashlib
import os
import sys
from pathlib import Path

ROOT = Path("/tmp/fac-debye-scf-test")
TARGET = ROOT / "faclib" / "orbital.c"
BACKUP = ROOT / "faclib" / "orbital.c.pre_debye_scf"
DIFF = ROOT / "patch_density.diff"
MARKER = "DEBYE_SCF_DIRECT"
SIG = ("int DensityToSZ(POTENTIAL *pot, double *d, double *z,\n"
       "\t\tdouble *zx, double *jps, int md) {\n")
# sha256 of the original function text (signature .. closing brace) at c4a016a
ORIG_SHA = "0c0980116d7c5e0fe95600e3fbcd1570b740dfc3f0fa488eb2a2b67cfa03fe24"

HELPERS = """/* DEBYE_SCF_DIRECT: Yukawa-screened direct potential helpers.
** sinhc(x) = sinh(x)/x.  _ee_wsinhc(w, x) = w*sinhc(x) evaluated without
** overflow: for x >= 700, sinh(x)/x = e^x/(2x) to double precision, so the
** product is formed as exp(log|w| + x - log(2x)); w == 0 gives exactly 0
** (never 0*inf). */
extern double EEScreenMu(void);

static double _ee_shx(double x) {
  if (x < 1e-3) {
    double x2 = x*x;
    return 1.0 + x2/6.0*(1.0 + x2/20.0);
  }
  return sinh(x)/x;
}

static double _ee_wsinhc(double w, double x) {
  if (w == 0) return 0.0;
  if (x < 700.0) return w*_ee_shx(x);
  return copysign(exp(log(fabs(w)) + x - log(2.0*x)), w);
}

"""

REPL = [
    ("  int i, im, maxrp, k;\n  double nx, rx, a, vxs;\n",
     "  int i, im, maxrp, k;\n  double nx, rx, a, vxs;\n"
     "  double mu = EEScreenMu(); /* DEBYE_SCF_DIRECT */\n"),
    ("""  im = maxrp-1;
  for (i = 0; i <= im; i++) {
    _dwork[i] = d[i]*pot->dr_drho[i];
  }
  _dwork1[0] = d[0]*pot->rad[0]/3.0;
  NewtonCotesIP(_dwork1, _dwork, 0, im, -1, 0);
  for (i = 0; i <= im; i++) {
    _dwork[i] = (d[i]/pot->rad[i])*pot->dr_drho[i];
  }
  _dwork2[im] = 0.0;
  NewtonCotesIP(_dwork2, _dwork, 0, im, -1, -1);
  for (i = 0; i <= im; i++) {
    z[i] = _dwork1[i] + pot->rad[i]*_dwork2[i];
  }\x20\x20
""",
     """  im = maxrp-1;
  if (mu > 0) {
    /* DEBYE_SCF_DIRECT: r*V for the kernel exp(-mu r12)/r12,
    ** z = e^{-mu r} int_0^r d sinhc(mu r') dr'
    **   + r sinhc(mu r) int_r^inf (d/r') e^{-mu r'} dr' */
    for (i = 0; i <= im; i++) {
      _dwork[i] = _ee_wsinhc(d[i]*pot->dr_drho[i], mu*pot->rad[i]);
    }
    _dwork1[0] = _ee_wsinhc(d[0]*pot->rad[0]/3.0, mu*pot->rad[0]);
    NewtonCotesIP(_dwork1, _dwork, 0, im, -1, 0);
    for (i = 0; i <= im; i++) {
      _dwork[i] = (d[i]/pot->rad[i])*pot->dr_drho[i]*exp(-mu*pot->rad[i]);
    }
    _dwork2[im] = 0.0;
    NewtonCotesIP(_dwork2, _dwork, 0, im, -1, -1);
    for (i = 0; i <= im; i++) {
      rx = pot->rad[i];
      z[i] = exp(-mu*rx)*_dwork1[i] + _ee_wsinhc(rx*_dwork2[i], mu*rx);
    }
  } else {
  for (i = 0; i <= im; i++) {
    _dwork[i] = d[i]*pot->dr_drho[i];
  }
  _dwork1[0] = d[0]*pot->rad[0]/3.0;
  NewtonCotesIP(_dwork1, _dwork, 0, im, -1, 0);
  for (i = 0; i <= im; i++) {
    _dwork[i] = (d[i]/pot->rad[i])*pot->dr_drho[i];
  }
  _dwork2[im] = 0.0;
  NewtonCotesIP(_dwork2, _dwork, 0, im, -1, -1);
  for (i = 0; i <= im; i++) {
    z[i] = _dwork1[i] + pot->rad[i]*_dwork2[i];
  }\x20\x20
  }
"""),
]


def die(msg):
    sys.exit("patch_density.py: ABORT: " + msg)


def function_span(src):
    if src.count(SIG) != 1:
        die(f"DensityToSZ signature found {src.count(SIG)} times (expected 1)")
    a = src.index(SIG)
    b = src.index("\n}\n", a) + 3
    return a, b


def main():
    if Path.cwd().resolve() != ROOT:
        die(f"run from {ROOT} (cwd is {Path.cwd().resolve()})")
    if ROOT not in TARGET.resolve().parents:
        die("target outside the isolated worktree")
    src = TARGET.read_text()
    if MARKER in src:
        die("patch already applied (marker present); refusing to apply twice")
    a, b = function_span(src)
    func = src[a:b]
    sha = hashlib.sha256(func.encode()).hexdigest()
    if len(sys.argv) > 1 and sys.argv[1] == "--print-hash":
        print(sha)
        return
    if sha != ORIG_SHA:
        die(f"DensityToSZ differs from the expected c4a016a version (sha {sha})")
    new = func
    for old, rep in REPL:
        if new.count(old) != 1:
            die(f"expected block found {new.count(old)} times:\n{old}")
        new = new.replace(old, rep)
    out = src[:a] + HELPERS + new + src[b:]
    if BACKUP.exists():
        die(f"backup {BACKUP} already exists; inspect it before re-running")
    BACKUP.write_text(src)
    TARGET.write_text(out)
    # re-read and verify
    chk = TARGET.read_text()
    a2, b2 = function_span(chk)
    if chk[:a] != src[:a] or chk[b2:] != src[b:]:
        die("text outside the function and helper block changed (BUG); "
            f"original preserved in {BACKUP}")
    if chk.count(MARKER) != 3 or "EEScreenMu()" not in chk[a2:b2]:
        die("verification of patched function failed")
    diff = difflib.unified_diff(src.splitlines(True), chk.splitlines(True),
                                "a/faclib/orbital.c", "b/faclib/orbital.c")
    DIFF.write_text("".join(diff))
    print(f"patched {TARGET}; backup {BACKUP}; diff {DIFF}")


if __name__ == "__main__":
    main()
