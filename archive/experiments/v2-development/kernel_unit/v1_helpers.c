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

