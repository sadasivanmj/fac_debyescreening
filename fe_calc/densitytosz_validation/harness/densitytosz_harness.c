/* Calls the REAL compiled DensityToSZ() from libfac.a on FAC's own radial grid.
**
** usage:  densitytosz_harness grid    MU            > grid.txt
**         densitytosz_harness eval    MU DENSITYFILE [MAXRP] > out.txt
**         (grid mode: densitytosz_harness grid MU - [MAXRP])
**
** MU = 0 : plain FAC (no PlasmaScreen, no ee_screen)   -> original Coulomb path
** MU > 0 : radial:ee_screen = 1 and PlasmaScreen with a Debye length 1/MU (same T, n mapping as the
**          energy calculations), so EEScreenMu() = MU exactly as in a real screened SCF.
** MU < 0 : "ee_screen = 0 with the plasma on" (EEScreenMu() = 0, nuclear screening still configured),
**          used to check that disabling the direct-potential screening is a pure Coulomb DensityToSZ call.
**
** The grid is that of a free-ion Fe III-like calculation (N = 24, N1 = 23), built by SetOrbitalRGrid before
** the plasma is switched on (so no ion-sphere cut-off), exactly as in the energy calculations.
** DENSITYFILE: n, then n values of the radial density D(r) (integral of D dr = N, FAC's convention).
** Output: header lines starting with '#', then  index r dr_drho D rV_FAC  per grid point.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "init.h"
#include "nucleus.h"
#include "radial.h"
#include "orbital.h"

#ifndef NO_EESCREENMU      /* the pristine upstream library has no screening: build with -DNO_EESCREENMU */
extern double EEScreenMu(void);
#else
static double EEScreenMu(void) { return 0.0; }
#endif

#define DENSITY_A3 1e-6
#define BOHR_A 0.52917721067
#define HARTREE_EV 27.211386018

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: grid|eval MU [DENSITYFILE]\n"); return 2; }
  int do_eval = (strcmp(argv[1], "eval") == 0);
  double mu = atof(argv[2]);
  int i, n;

  InitFac0();                               /* what 'import pfac' does at module initialisation */
  InitFac1();
  SetAtom("Fe", -1, -1, -1, -1, -1);
  ReinitRadial(0);                          /* initialises the potential for the atom, as the first FAC calculation does */
  int maxrp_arg = (argc > 4) ? atoi(argv[4]) : -1;      /* optional finer/coarser grid for the convergence study */
  SetRadialGrid(maxrp_arg, -1, -1, -1, -1);
  POTENTIAL *pot = RadialPotential();
  pot->N = 24.0;
  pot->N1 = 23.0;
  SetOrbitalRGrid(pot);                     /* free-ion grid, mps = -1 */

  if (mu != 0.0) {
    if (mu > 0) SetOption("radial:ee_screen", "1", 1, 1.0);
    else SetOption("radial:ee_screen", "0", 0, 0.0);
    double m = fabs(mu);
    double tps = 4.0*M_PI*DENSITY_A3*pow(BOHR_A, 3)/(m*m)*HARTREE_EV;   /* lambda_D = 1/mu */
    PlasmaScreen(1, 1, 23.0, DENSITY_A3, tps, 0.0);
    InitializePS(pot);                      /* sets dps = sqrt(T/(4 pi n)) */
  }
  n = pot->maxrp;
  printf("# maxrp %d ips %d mps %d dps %.17g ahx %.17g EEScreenMu %.17g ups %.17g\n",
         n, pot->ips, pot->mps, pot->dps, pot->ahx, EEScreenMu(), pot->ups);

  double *d = calloc(n, sizeof(double)), *z = calloc(n, sizeof(double)), *zx = calloc(n, sizeof(double));
  double jps = 0.0;
  if (do_eval) {
    FILE *f = fopen(argv[3], "r");
    int nf;
    if (!f || fscanf(f, "%d", &nf) != 1 || nf != n) { fprintf(stderr, "density file does not match grid (%d)\n", n); return 3; }
    for (i = 0; i < n; i++) if (fscanf(f, "%lf", d + i) != 1) return 4;
    fclose(f);
    int im = DensityToSZ(pot, d, z, zx, &jps, 0);
    printf("# DensityToSZ returned %d jps %.17g\n", im, jps);
  }
  for (i = 0; i < n; i++)
    printf("%d %.17g %.17g %.17g %.17g\n", i, pot->rad[i], pot->dr_drho[i], d[i], z[i]);
  return 0;
}
