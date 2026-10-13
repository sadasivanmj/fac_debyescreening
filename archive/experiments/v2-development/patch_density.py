#!/usr/bin/env python3
"""v2: numerically stable Debye-Hueckel screened direct potential in DensityToSZ.

Isolated-worktree tool (/tmp/fac-debye-scf-test only).  v1 of this script is
preserved as patch_density_v1.py (its diff as patch_density_v1.diff and its
patched source as faclib/orbital.c.debye_scf_v1).

v2 is built from the ORIGINAL c4a016a source (faclib/orbital.c.pre_debye_scf).
Accepted current states of faclib/orbital.c: the original, or exactly v1.

FAC conventions (verified in the audit): d = radial density (int d dr = N),
grid with unit steps in rho, dr_drho = dr/drho, z = r V_direct, jps = the
outer integral at r_0.  With E(x) = (1 - e^{-2x})/(2x) = e^{-x} sinh(x)/x:
    z(r)  = At(r) + r E(mu r) Bt(r)
    At(r) = int_0^r d(r') E(mu r') e^{-mu(r-r')} dr'      (= e^{-mu r} A(r))
    Bt(r) = int_r^inf (d(r')/r') e^{-mu(r'-r)} dr'         (= e^{mu r} B(r))
All factors are bounded; cumulative integrals use NewtonCotesIP on blocks of
exponent span <= 100 (a single block when mu (r_im - r_0) <= 100).
"""

import difflib
import hashlib
import sys
from pathlib import Path

ROOT = Path("/tmp/fac-debye-scf-test")
TARGET = ROOT / "faclib" / "orbital.c"
ORIGINAL = ROOT / "faclib" / "orbital.c.pre_debye_scf"
V1_SAVED = ROOT / "faclib" / "orbital.c.debye_scf_v1"
DIFF = ROOT / "patch_density.diff"
DIFF_V1V2 = ROOT / "patch_density_v1_to_v2.diff"
MARKER = "DEBYE_SCF_V2"
SHA_ORIGINAL_FILE = "59d695abcde1ab6f13a57cf9c00fcf0689ab8af1d3d895637292d28f0391c40c"
SHA_V1_FILE = "bbef452c2463102190083f81b591136189e6bfb3cbed32dcc6872811d2ebe9b7"
SHA_ORIGINAL_FUNC = "0c0980116d7c5e0fe95600e3fbcd1570b740dfc3f0fa488eb2a2b67cfa03fe24"
SIG = ("int DensityToSZ(POTENTIAL *pot, double *d, double *z,\n"
       "\t\tdouble *zx, double *jps, int md) {\n")

HELPERS = r"""/* DEBYE_SCF_V2: Debye-Hueckel (Yukawa) screened direct potential.
** With E(x) = (1 - exp(-2x))/(2x) = exp(-x) sinh(x)/x  (E(0) = 1, 0 < E <= 1),
**   z(r)  = r V(r) = At(r) + r E(mu r) Bt(r),
**   At(r) = exp(-mu r) int_0^r d sinhc(mu r') dr'
**         = int_0^r d(r') E(mu r') exp(-mu (r - r')) dr',
**   Bt(r) = exp(mu r) int_r^inf (d/r') exp(-mu r') dr'
**         = int_r^inf (d/r') exp(-mu (r' - r)) dr'.
** Every factor is bounded, so no intermediate value can overflow.  The
** cumulative integrals use NewtonCotesIP on blocks whose exponent span
** mu*(r_end - r_start) is <= _EE_SPAN: inside a block the integrand is
** scaled by exp(-mu |R - r'|) relative to the block reference R and the
** cumulative result is rescaled by exp(mu |R - r|) <= exp(_EE_SPAN).  When
** mu*(r_im - r_0) <= _EE_SPAN this is one block over the whole range.
** _EE_SPAN = 100 keeps scaled integrands far above the absolute 1e-100
** below which uvip3s (interpolation.c) silently falls back to the
** trapezoid rule; a span of 300 measurably degraded the quadrature.
** Returns B(r_0) = int_{r_0}^inf (d/r') exp(-mu r') dr', the screened
** analogue of the Coulomb jps. */
#define _EE_SPAN 100.0
extern double EEScreenMu(void);

static double _ee_E(double x) {
  if (x == 0) return 1.0;
  return -expm1(-2.0*x)/(2.0*x);
}

static double _ee_screened_sz(int im, double mu, double *rad, double *drdrho,
			      double *d, double *z, double *w,
			      double *at, double *bt) {
  int i, s, e;
  double R, L;

  at[0] = d[0]*rad[0]/3.0*_ee_E(mu*rad[0]);
  s = 0;
  while (s < im) {
    for (e = s+1; e < im && mu*(rad[e+1]-rad[s]) <= _EE_SPAN; e++);
    R = rad[e];
    for (i = s; i <= e; i++) {
      w[i] = d[i]*_ee_E(mu*rad[i])*drdrho[i]*exp(-mu*(R-rad[i]));
    }
    at[s] *= exp(-mu*(R-rad[s]));
    NewtonCotesIP(at, w, s, e, -1, 0);
    for (i = s; i <= e; i++) {
      at[i] *= exp(mu*(R-rad[i]));
    }
    s = e;
  }
  bt[im] = 0.0;
  e = im;
  while (e > 0) {
    for (s = e-1; s > 0 && mu*(rad[e]-rad[s-1]) <= _EE_SPAN; s--);
    L = rad[s];
    for (i = s; i <= e; i++) {
      w[i] = (d[i]/rad[i])*drdrho[i]*exp(-mu*(rad[i]-L));
    }
    bt[e] *= exp(-mu*(rad[e]-L));
    NewtonCotesIP(bt, w, s, e, -1, -1);
    for (i = s; i <= e; i++) {
      bt[i] *= exp(mu*(rad[i]-L));
    }
    e = s;
  }
  for (i = 0; i <= im; i++) {
    z[i] = at[i] + rad[i]*_ee_E(mu*rad[i])*bt[i];
  }
  return exp(-mu*rad[0])*bt[0];
}
/* DEBYE_SCF_V2_END */

"""

OLD_DECL = "  int i, im, maxrp, k;\n  double nx, rx, a, vxs;\n"
NEW_DECL = OLD_DECL + """  double mu = EEScreenMu(); /* DEBYE_SCF_V2 */

  if (!(mu >= 0.0) || isinf(mu)) { /* DEBYE_SCF_V2: NaN, negative or inf */
    printf("invalid Debye screening parameter in DensityToSZ: %g %g\\n",
\t   mu, pot->dps);
    Abort(1);
  }
"""
OLD_BLOCK = """  im = maxrp-1;
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
"""
COULOMB = OLD_BLOCK.split("  im = maxrp-1;\n", 1)[1]
NEW_BLOCK = ("  im = maxrp-1;\n"
             "  if (mu > 0) { /* DEBYE_SCF_V2: screened; _dwork2[0] = jps */\n"
             "    _dwork2[0] = _ee_screened_sz(im, mu, pot->rad, pot->dr_drho, d, z,\n"
             "\t\t\t\t _dwork, _dwork1, _dwork2);\n"
             "  } else {\n" + COULOMB + "  }\n")


def die(msg):
    sys.exit("patch_density.py: ABORT: " + msg)


def sha(t):
    return hashlib.sha256(t.encode()).hexdigest()


def span(src):
    if src.count(SIG) != 1:
        die(f"DensityToSZ signature found {src.count(SIG)} times")
    a = src.index(SIG)
    return a, src.index("\n}\n", a) + 3


def main():
    if Path.cwd().resolve() != ROOT:
        die(f"run from {ROOT}")
    for p in (TARGET, ORIGINAL):
        if ROOT not in p.resolve().parents:
            die(f"{p} outside the isolated worktree")
    orig = ORIGINAL.read_text()
    if sha(orig) != SHA_ORIGINAL_FILE:
        die("original backup does not match the recorded c4a016a orbital.c")
    cur = TARGET.read_text()
    if MARKER in cur:
        die("v2 already applied; refusing to apply twice")
    state = {SHA_ORIGINAL_FILE: "original", SHA_V1_FILE: "v1"}.get(sha(cur))
    if state is None:
        die("faclib/orbital.c is neither the original nor v1; inspect it first")
    a, b = span(orig)
    func = orig[a:b]
    if sha(func) != SHA_ORIGINAL_FUNC:
        die("original DensityToSZ differs from the expected version")
    for old in (OLD_DECL, OLD_BLOCK):
        if func.count(old) != 1:
            die(f"expected block found {func.count(old)} times:\n{old}")
    new_func = func.replace(OLD_DECL, NEW_DECL).replace(OLD_BLOCK, NEW_BLOCK)
    out = orig[:a] + HELPERS + new_func + orig[b:]
    if state == "v1" and sha(V1_SAVED.read_text()) != SHA_V1_FILE:
        die("v1 state not preserved in orbital.c.debye_scf_v1")
    TARGET.write_text(out)
    chk = TARGET.read_text()
    a2, b2 = span(chk)
    f2 = chk[a2:b2]
    ok = (chk[:a] == orig[:a] and chk[b2:] == orig[b:]
          and chk[a:a2] == HELPERS
          and "double mu = EEScreenMu();" in f2 and "Abort(1);" in f2
          and "_ee_screened_sz(im, mu" in f2
          and ("} else {\n" + COULOMB) in f2
          and f2.replace(NEW_DECL, OLD_DECL).replace(NEW_BLOCK, OLD_BLOCK) == func)
    if not ok:
        TARGET.write_text(cur)
        die("verification failed; previous orbital.c restored")
    DIFF.write_text("".join(difflib.unified_diff(
        orig.splitlines(True), chk.splitlines(True),
        "a/faclib/orbital.c", "b/faclib/orbital.c")))
    v1 = V1_SAVED.read_text()
    DIFF_V1V2.write_text("".join(difflib.unified_diff(
        v1.splitlines(True), chk.splitlines(True),
        "v1/faclib/orbital.c", "v2/faclib/orbital.c")))
    print(f"v2 applied (previous state: {state}); diff vs original: {DIFF}; "
          f"diff v1->v2: {DIFF_V1V2}")


if __name__ == "__main__":
    main()
