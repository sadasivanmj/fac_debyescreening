#!/usr/bin/env python3
"""v1 vs v2(span S) vs analytic Yukawa potential of a 1s cloud, all grid points."""
import numpy as np, subprocess, json
def zan(r,a,N,mu):
    b=2*a
    if mu==0: return N*(1-np.exp(-b*r)*(1+b*r/2))
    D=b*b-mu*mu
    return N*(b**4/D**2*(np.exp(-mu*r)-np.exp(-b*r)) - b**3*r*np.exp(-b*r)/(2*D))
res={}
print(f"{'case':8s} {'mu':>7s} {'coul':>9s} {'v1':>9s} " + " ".join(f"{'S='+s:>9s}" for s in ['300','150','100','60','30']) + "   (max|z-z_an|/N; nan = non-finite)")
for cname,a,N in [("diffuse006",0.06,1.0)]:
    for mu in [0.0,1e-8,1e-4,0.1,0.25,1.0,5.0]:
        cells=[]; row={}
        for S in ['300','150','100','60','30']:
            A=np.loadtxt(subprocess.run([f"./harness_S{S}",f"grid_{cname}.txt",repr(mu)],capture_output=True,text=True).stdout.splitlines()[1:])
            r,zc,z1,z2=A.T; ref=zan(r,a,N,mu)
            e=lambda z: float(np.max(np.abs(z-ref))/N) if np.all(np.isfinite(z)) else float('nan')
            row[S]=e(z2)
        row['coul']=e(zc) if mu==0 else None; row['v1']=e(z1)
        res[f"{cname}_{mu}"]=row
        f=lambda v: "     n/a" if v is None else f"{v:9.1e}"
        print(f"{cname:8s} {mu:7g} {f(row['coul'])} {f(row['v1'])} "+" ".join(f(row[S]) for S in ['300','150','100','60','30']),flush=True)
json.dump(res,open("span_scan.json","w"),indent=1)
