/* DEBYE_SCF_V2: Debye-Hueckel (Yukawa) screened direct potential.
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
