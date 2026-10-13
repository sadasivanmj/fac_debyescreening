#!/usr/bin/env python3
"""Synthetic kernel tests: original Coulomb, v1, v2 vs scipy-quad reference."""
import numpy as np, subprocess, math, json
from scipy.integrate import quad

def grid(n=2000, rmin=1e-6, rmax=1500.0, ag=25.0):
    bg = (n-1 - ag*(math.sqrt(rmax)-math.sqrt(rmin)))/math.log(rmax/rmin)
    rho0 = ag*math.sqrt(rmin) + bg*math.log(rmin)
    r = np.empty(n); r[0] = rmin; x = rmin
    for i in range(1, n):
        target = rho0 + i
        for _ in range(100):
            f = ag*math.sqrt(x) + bg*math.log(x) - target
            x_new = x - f/(ag/(2*math.sqrt(x)) + bg/x)
            if x_new <= 0: x_new = x/2
            if abs(x_new-x) < 1e-15*x: x = x_new; break
            x = x_new
        r[i] = x
    r[-1] = rmax if abs(r[-1]-rmax) < 1e-6*rmax else r[-1]
    drdrho = r/(ag*np.sqrt(r)/2 + bg)
    return r, drdrho

def dens(r, a, N):  # radial density of a 1s-type cloud, int d dr = N
    return N*4*a**3*r**2*np.exp(-2*a*r)

def E(x): return 1.0 if x == 0 else -math.expm1(-2*x)/(2*x)

def zref(rv, a, N, mu):
    f = lambda t: N*4*a**3*t*t*math.exp(-2*a*t)
    i1 = quad(lambda t: f(t)*E(mu*t)*math.exp(-mu*(rv-t)), 0, rv, limit=400,
              points=[min(rv, 1/a)], epsabs=0, epsrel=1e-13)[0]
    i2 = quad(lambda t: f(t)/t*math.exp(-mu*(t-rv)), rv, np.inf, limit=400,
              epsabs=0, epsrel=1e-13)[0]
    return i1 + rv*E(mu*rv)*i2

r, dr = grid()
out = {}
for cname, a, N in [("compact", 13.0, 2.0), ("valence", 1.0, 18.0), ("diffuse", 0.05, 1.0)]:
    d = dens(r, a, N)
    fn = f"grid_{cname}.txt"
    with open(fn, "w") as fh:
        fh.write(f"{len(r)}\n")
        for i in range(len(r)): fh.write(f"{r[i]:.17g} {dr[i]:.17g} {d[i]:.17g}\n")
    idx = np.unique(np.geomspace(1, len(r)-1, 70).astype(int))
    for mu in [0.0, 1e-8, 1e-4, 0.1, 0.25, 1.0, 5.0]:
        txt = subprocess.run(["./harness", fn, repr(mu)], capture_output=True, text=True).stdout
        lines = txt.splitlines(); hdr = lines[0].split()
        A = np.array([[float(v) for v in l.split()] for l in lines[1:]])
        zc, z1, z2 = A[:,1], A[:,2], A[:,3]
        ref = np.array([zref(r[i], a, N, mu) for i in idx])
        sc = np.max(np.abs(ref))
        def err(z):
            if not np.all(np.isfinite(z[idx])): return float('nan')
            return float(np.max(np.abs(z[idx]-ref))/sc)
        rec = dict(case=cname, mu=mu, max_mu_r=mu*r[-1],
                   nonfinite_v1=int(np.sum(~np.isfinite(z1))), nonfinite_v2=int(np.sum(~np.isfinite(z2))),
                   err_coul=err(zc) if mu == 0 else None, err_v1=err(z1), err_v2=err(z2),
                   v2_minus_v1=float(np.nanmax(np.abs(z2-z1))/sc) if np.all(np.isfinite(z1)) else None,
                   v2_minus_coul=float(np.max(np.abs(z2-zc))/sc),
                   jps=[float(hdr[4]), float(hdr[6]), float(hdr[8])])
        out[f"{cname}_{mu}"] = rec
        fe = lambda v, f=".2e": "n/a" if v is None else format(v, f)
        print(f"{cname:8s} mu={mu:<7g} mu*rmax={mu*r[-1]:8.1f} nonfinite v1/v2={rec['nonfinite_v1']}/{rec['nonfinite_v2']} "
              f"err vs ref: coul={fe(rec['err_coul'])} v1={fe(rec['err_v1'])} v2={fe(rec['err_v2'])} "
              f"|v2-v1|={fe(rec['v2_minus_v1'], '.1e')} |v2-coul|={fe(rec['v2_minus_coul'], '.1e')}", flush=True)
json.dump(out, open("unit_results.json", "w"), indent=1)
