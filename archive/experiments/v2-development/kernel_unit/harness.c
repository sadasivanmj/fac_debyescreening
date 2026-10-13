/* Kernel unit harness: original Coulomb loops, v1 arithmetic, v2 kernel.
** Uses the real NewtonCotesIP from the worktree libfac.a.  Input: grid file
** (n, then r drdrho d per line); output columns: r z_coul z_v1 z_v2 jps info. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int NewtonCotesIP(double *r, double *x, int i0, int i1, int m, int id);
static void Abort(int r) { exit(r); }
#include "v1_helpers.c"
#include "v2_helpers.c"
int main(int argc, char **argv) {
  double mu = atof(argv[2]); int n, i, im;
  FILE *f = fopen(argv[1], "r"); if (fscanf(f, "%d", &n) != 1) return 2;
  double *r = malloc(8*n), *dr = malloc(8*n), *d = malloc(8*n), *w = malloc(8*n),
         *A = malloc(8*n), *B = malloc(8*n), *zc = malloc(8*n), *z1 = malloc(8*n), *z2 = malloc(8*n);
  for (i = 0; i < n; i++) if (fscanf(f, "%lf %lf %lf", r+i, dr+i, d+i) != 3) return 3;
  im = n-1;
  /* original Coulomb statements (verbatim structure) */
  for (i = 0; i <= im; i++) w[i] = d[i]*dr[i];
  A[0] = d[0]*r[0]/3.0; NewtonCotesIP(A, w, 0, im, -1, 0);
  for (i = 0; i <= im; i++) w[i] = (d[i]/r[i])*dr[i];
  B[im] = 0.0; NewtonCotesIP(B, w, 0, im, -1, -1);
  for (i = 0; i <= im; i++) zc[i] = A[i] + r[i]*B[i];
  double jc = B[0];
  /* v1 arithmetic (mu > 0 branch of patch v1) */
  for (i = 0; i <= im; i++) w[i] = _ee_wsinhc(d[i]*dr[i], mu*r[i]);
  A[0] = _ee_wsinhc(d[0]*r[0]/3.0, mu*r[0]); NewtonCotesIP(A, w, 0, im, -1, 0);
  for (i = 0; i <= im; i++) w[i] = (d[i]/r[i])*dr[i]*exp(-mu*r[i]);
  B[im] = 0.0; NewtonCotesIP(B, w, 0, im, -1, -1);
  for (i = 0; i <= im; i++) z1[i] = exp(-mu*r[i])*A[i] + _ee_wsinhc(r[i]*B[i], mu*r[i]);
  double j1 = B[0];
  /* v2 kernel */
  double j2 = _ee_screened_sz(im, mu, r, dr, d, z2, w, A, B);
  printf("# mu %.17g jps_coul %.17g jps_v1 %.17g jps_v2 %.17g\n", mu, jc, j1, j2);
  for (i = 0; i <= im; i++) printf("%.17g %.17g %.17g %.17g\n", r[i], zc[i], z1[i], z2[i]);
  return 0;
}
