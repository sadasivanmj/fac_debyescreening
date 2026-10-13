#!/usr/bin/env python3
"""
chen2018_reproduce.py
=====================

Reproduce Chen, Ma, Hu & Wang, Phys. Plasmas 25, 072120 (2018), Table II,
for He-like C4+, Ne8+, Ar16+ and Kr34+ in a Debye-Hueckel plasma.

Two physical models are run:

    model "b"  (ee=0):  V_en = -Z exp(-r/lambda_D)/r     V_ee = 1/r12
    model "g"  (ee=1):  V_en = -Z exp(-r/lambda_D)/r     V_ee = exp(-r12/lambda_D)/r12

REQUIREMENTS
------------
Both FAC patches must be applied and built:

    patch -p1 < fac-debye-nuclear-screening.patch     # e-n sign fix   (orbital.c)
    patch -p1 < fac-debye-ee-screening.patch          # e-e Yukawa     (radial.c, orbital.c)
    ./configure --prefix=$HOME/facinst && make -j8

Set SFAC below to the resulting sfac binary.

CONVENTIONS THAT MATTER
-----------------------
1. FAC has no direct lambda_D input.  It derives
       dps = sqrt(tps / (4*pi*nps*(1+ups)))
   with nps in A^-3 and tps in eV.  This script inverts that, with ups = 0.

2. The lambda_D -> infinity control MUST be a run at lambda_D = 1e6 a0 with the
   plasma path active (mps=1).  A run without PlasmaScreen() uses a different
   zeroth-order Hamiltonian (VT*r -> -(Z-N1) instead of -(Z-N)) and carries a
   constant offset of ~1e-2 a.u. in the excitation energies.  Using it as the
   reference contaminates every shift.

3. Screening is off by default.  Both SetOption() calls are required.
"""

import os, re, sys, math, json, subprocess

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
SFAC     = os.environ.get("SFAC", "sfac")   # path to the patched sfac binary
WORKDIR  = "chen2018_run"
NPS      = 1.0e-2                            # A^-3; any value works, tps compensates
FREE_LAM = 1.0e6                             # a0; the lambda -> infinity control
LAMBDAS  = [FREE_LAM, 20.0, 15.0, 10.0, 7.5, 5.0]
MODELS   = [0, 1]                            # 0 -> Chen "b", 1 -> Chen "g"

RBOHR      = 0.52917721067
HARTREE_EV = 27.211386018

IONS = {                # element symbol, net charge of the He-like ion
    "C":  4.0,
    "Ne": 8.0,
    "Ar": 16.0,
    "Kr": 34.0,
}
ION_ORDER = ["C", "Ne", "Ar", "Kr"]
ION_LABEL = {"C": "C4+", "Ne": "Ne8+", "Ar": "Ar16+", "Kr": "Kr34+"}

# --------------------------------------------------------------------------
# Chen et al. (2018) Table II, transcribed from the paper and checked line by line
#   a = MCDF, electron-nucleus screening only
#   b = FAC,  electron-nucleus screening only
#   f = MCDF, electron-nucleus AND e-e screening
#   g = FAC,  electron-nucleus AND e-e screening
# Column order: [free ion, 20, 15, 10, 7.5, 5] a0.  None = not reported.
# f and g are reported by the paper for two levels only.
#
# NOTE ON THREE NIST VALUES: the paper's NIST column for Ar16+ 1s2p 3P2,
# Ar16+ 1s3p 3P1 and Kr34+ 1s3p 1P1 is 114.8890, 135.2165 and 567.2068
# respectively.  Transcriptions that give 115.1025, 135.2872 and 567.5555 for
# those three have picked up the adjacent "b" or "a" free-ion entry instead.
# --------------------------------------------------------------------------
CHEN = {
 "C": {
  "1s2s 3S1": {"NIST":10.9866,
     "a":[10.9036,10.8975,10.8926,10.8800,10.8622,10.8137],
     "b":[11.1297,11.0882,11.0877,11.0862,11.0842,11.0056]},
  "1s2p 3P1": {"NIST":11.1865,
     "a":[11.1016,11.0966,11.0925,11.0819,11.0671,11.0261],
     "b":[11.3987,11.3472,11.3466,11.3450,11.3426,11.2426]},
  "1s2p 1P1": {"NIST":11.3151,
     "a":[11.2372,11.2318,11.2274,11.2160,11.2001,11.1563],
     "b":[11.5051,11.4844,11.4838,11.4820,11.4794,11.3777],
     "f":[None,   11.2309,11.2242,11.2112,11.1916,11.1443],
     "g":[None,   11.4771,11.4762,11.4688,11.4687,11.3685]},
  "1s2p 3P2": {"NIST":11.1871,
     "a":[11.1028,11.0977,11.0936,11.0830,11.0682,11.0272],
     "b":[11.3987,11.3475,11.3469,11.3453,11.3429,11.2430]},
  "1s3p 3P1": {"NIST":12.9920,
     "a":[12.9042,12.8887,12.8768,12.8448,12.8012,12.6829],
     "b":[13.0612,13.0139,13.0125,13.0085,13.0033,12.9629]},
  "1s3p 1P1": {"NIST":13.0282,
     "a":[12.9432,12.9272,12.9149,12.8818,12.8363,12.7149],
     "b":[13.1004,13.0529,13.0512,13.0466,13.0402,12.9936]},
 },
 "Ne": {
  "1s2s 3S1": {"NIST":33.2610,
     "a":[33.1828,33.1770,33.1728,33.1598,33.1422,33.0928],
     "b":[33.4011,33.3584,33.3580,33.3569,33.3553,33.2778]},
  "1s2p 3P1": {"NIST":33.6189,
     "a":[33.6250,33.6203,33.6105,33.6061,33.5918,33.4660],
     "b":[33.8290,33.7762,33.7757,33.7745,33.7727,33.6756]},
  "1s3p 3P1": {"NIST":39.3872,
     "a":[39.3042,39.2895,39.2787,39.2470,39.2040,39.0853],
     "b":[39.4539,39.4073,39.4063,39.4037,39.4000,39.3594]},
  "1s3p 1P1": {"NIST":39.4603,
     "a":[39.3793,39.3643,39.3538,39.3209,39.2772,39.1574],
     "b":[39.5300,39.4834,39.4823,39.4793,39.4751,39.4309]},
 },
 "Ar": {
  "1s2s 3S1": {"NIST":114.0754,
     "a":[114.0262,114.0196,114.0154,114.0038,113.9864,113.9357],
     "b":[114.2159,114.1722,114.1720,114.1712,114.1701,114.0938]},
  "1s2p 3P1": {"NIST":114.7878,
     "a":[114.7419,114.7366,114.7334,114.7242,114.7104,114.6698],
     "b":[115.0051,114.9509,114.9506,114.9490,114.9481,114.8528]},
  "1s2p 1P1": {"NIST":115.3775,
     "a":[115.3388,115.3357,115.3323,115.3230,115.3088,115.2671],
     "b":[115.6022,115.5488,115.5484,115.5471,115.5460,115.4508]},
  "1s2p 3P2": {"NIST":114.8890,
     "a":[114.8427,114.8378,114.8345,114.8253,114.8114,114.7706],
     "b":[115.1025,115.0487,115.0484,115.0472,115.0461,114.9516]},
  "1s3p 3P1": {"NIST":135.2165,
     "a":[135.1670,135.1520,135.1413,135.1113,135.0687,134.9479],
     "b":[135.2872,135.2410,135.2404,135.2381,135.2361,135.1959]},
  "1s3p 1P1": {"NIST":135.3789,
     "a":[135.3323,135.3177,135.3068,135.2765,135.2334,135.1112],
     "b":[135.4521,135.4061,135.4055,135.4031,135.4008,135.3592]},
 },
 "Kr": {
  "1s2s 3S1": {"NIST":476.9807,
     "a":[477.2644,477.2589,477.2504,477.2423,477.2253,477.1769],
     "b":[477.1289,477.0831,477.0830,477.0825,477.0819,477.0035]},
  "1s2p 3P1": {"NIST":478.7028,
     "a":[479.0517,479.0474,479.0427,479.0344,479.0211,478.9831],
     "b":[478.9663,478.9086,478.9084,478.9079,478.9070,478.8067]},
  "1s2p 1P1": {"NIST":481.9492,
     "a":[482.2991,482.2946,482.2901,482.2812,482.2673,482.2279],
     "b":[482.1981,482.1428,482.1426,482.1421,482.1412,482.0457],
     "f":[None,    482.2936,482.2888,482.2792,482.2648,482.2243],
     "g":[None,    482.1418,482.1414,482.1403,482.1390,482.0424]},
  "1s2p 3P2": {"NIST":481.0825,
     "a":[481.4235,481.4190,481.4105,481.4057,481.3920,481.3529],
     "b":[481.3254,481.2698,481.2696,481.2691,481.2683,481.1730]},
  "1s3p 3P1": {"NIST":566.2595,
     "a":[566.6075,566.5937,566.5775,566.5529,566.5109,566.3923],
     "b":[566.3631,566.3153,566.3149,566.3140,566.3124,566.2702]},
  "1s3p 1P1": {"NIST":567.2068,
     "a":[567.5555,567.5415,567.5355,567.5000,567.4573,567.3366],
     "b":[567.3063,567.2599,567.2596,567.2586,567.2568,567.2167]},
 },
}

LEVELS = ["1s2s 3S1", "1s2p 3P1", "1s2p 1P1", "1s2p 3P2", "1s3p 3P1", "1s3p 1P1"]

# --------------------------------------------------------------------------
# FAC driving
# --------------------------------------------------------------------------
def tps_for(lam, nps=NPS):
    """eV temperature that makes dps == lam (bohr), given nps in A^-3 and ups=0."""
    return lam * lam * 4.0 * math.pi * (nps * RBOHR**3) * HARTREE_EV


def sf_input(tag, elem, znet, lam, ee, with_n3=True):
    cfgs = ["Config('g1','1s2')", "Config('g2','1s1 2s1','1s1 2p1')"]
    grps = "'g1','g2'"
    if with_n3:
        cfgs.append("Config('g3','1s1 3s1','1s1 3p1','1s1 3d1')")
        grps += ",'g3'"
    s  = "SetAtom('%s')\n" % elem
    s += "SetOption('orbital:debye_mode', 1)\n"       # V_en = -Z exp(-r/lam)/r
    s += "SetOption('radial:ee_screen', %d)\n" % ee   # V_ee screened when ee=1
    s += "\n".join(cfgs) + "\n"
    s += "PlasmaScreen(%.10g, %.10g, %.14g, 1, 0.0, 1)\n" % (znet, NPS, tps_for(lam))
    s += "ConfigEnergy(0)\nOptimizeRadial(['g1'])\nConfigEnergy(1)\n"
    s += "GetPotential('%s.pot')\n" % tag
    s += "Structure('%s.lev.b',[%s])\n" % (tag, grps)
    s += "MemENTable('%s.lev.b')\n" % tag
    s += "PrintTable('%s.lev.b','%s.lev',1)\n" % (tag, tag)
    return s


def run_one(tag, elem, znet, lam, ee, with_n3=True):
    """Returns (E_ground_au, {level: E_exc_au}, note) or (None, None, error_text)."""
    path = os.path.join(WORKDIR, tag + ".sf")
    open(path, "w").write(sf_input(tag, elem, znet, lam, ee, with_n3))
    r = subprocess.run([SFAC, tag + ".sf"], cwd=WORKDIR,
                       capture_output=True, text=True, timeout=7200)
    if r.returncode != 0:
        msg = [l for l in r.stdout.split("\n") if "Error" in l or "Max iteration" in l]
        return None, None, (msg[0].strip() if msg else "sfac returned %d" % r.returncode)
    E0, lev = parse_lev(os.path.join(WORKDIR, tag + ".lev"))
    return E0, identify(lev), ("n<=2 only" if not with_n3 else "")


def parse_lev(fn):
    E0 = None
    lev = []
    for L in open(fn):
        if L.startswith("E0"):
            E0 = float(L.split(",")[1]) / HARTREE_EV
        f = L.split()
        if len(f) > 6 and re.match(r"^\d+$", f[0]):
            lev.append(dict(i=int(f[0]), e=float(f[2]) / HARTREE_EV,
                            vnl=f[4], tj=int(f[5])))
    return E0, lev


def identify(lev):
    """Map FAC levels to the Chen level labels.

    vnl encodes the outer orbital: 200 = 2s, 201 = 2p, 301 = 3p.
    Within 1s2p there are two 2J=2 levels; the lower is 3P1, the upper 1P1
    (same for 1s3p).  2J=4 in 1s2p is 3P2.
    """
    def pick(vnl, tj, which=0):
        c = sorted([l for l in lev if l["vnl"] == vnl and l["tj"] == tj],
                   key=lambda l: l["e"])
        return c[which]["e"] if len(c) > which else None
    return {
        "1s2s 3S1": pick("200", 2, 0),
        "1s2p 3P1": pick("201", 2, 0),
        "1s2p 1P1": pick("201", 2, 1),
        "1s2p 3P2": pick("201", 4, 0),
        "1s3p 3P1": pick("301", 2, 0),
        "1s3p 1P1": pick("301", 2, 1),
    }


# --------------------------------------------------------------------------
# main sweep
# --------------------------------------------------------------------------
def sweep(cache="chen2018_results.json", force=False):
    os.makedirs(WORKDIR, exist_ok=True)
    if os.path.exists(cache) and not force:
        return json.load(open(cache))
    res, fails = {}, []
    for elem in ION_ORDER:
        znet = IONS[elem]
        for ee in MODELS:
            for lam in LAMBDAS:
                tag = "%s_ee%d_L%s" % (elem, ee, ("%g" % lam).replace(".", "p").replace("+", ""))
                E0, lv, note = run_one(tag, elem, znet, lam, ee, with_n3=True)
                if E0 is None:
                    # Under strong screening the n=3 states can leave the bound
                    # spectrum: with the nucleus screened and e-e unscreened the
                    # asymptotic one-electron potential is +1/r (repulsive).
                    # Retry with n<=2 rather than losing the whole point.
                    fails.append(dict(tag=tag, elem=elem, ee=ee, lam=lam,
                                      stage="n<=3", error=note))
                    E0, lv, note = run_one(tag + "_n2", elem, znet, lam, ee, with_n3=False)
                    if E0 is None:
                        fails.append(dict(tag=tag, elem=elem, ee=ee, lam=lam,
                                          stage="n<=2", error=note))
                        continue
                res[tag] = dict(elem=elem, ee=ee, lam=lam, E0=E0, lev=lv, note=note)
                print("  %-22s E0=%14.6f  %s" % (tag, E0, note), flush=True)
    out = dict(results=res, failures=fails)
    json.dump(out, open(cache, "w"), indent=1)
    return out


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
LAMLAB = ["inf", "20", "15", "10", "7.5", "5"]


def get(res, elem, ee, j, level):
    tag = "%s_ee%d_L%s" % (elem, ee, ("%g" % LAMBDAS[j]).replace(".", "p").replace("+", ""))
    r = res.get(tag)
    return None if r is None else r["lev"].get(level)


def report(data, out="chen2018_report.txt"):
    res, fails = data["results"], data["failures"]
    L = []
    A = L.append

    A("=" * 100)
    A("REPRODUCTION OF CHEN et al., Phys. Plasmas 25, 072120 (2018), TABLE II")
    A("All energies in atomic units.  ee=0 -> Chen model 'b'.  ee=1 -> Chen model 'g'.")
    A("Reference for all shifts: lambda_D = %g a0, run on the same mps=1 path." % FREE_LAM)
    A("=" * 100)

    # ---- 1. absolute excitation energies -------------------------------
    A("\n\n1. ABSOLUTE EXCITATION ENERGIES\n")
    for elem in ION_ORDER:
        A("--- He-like %s" % ION_LABEL[elem])
        A("%-11s %5s %12s %12s %12s %12s %12s" %
          ("level", "lam", "this(ee=0)", "Chen b", "this(ee=1)", "Chen g", "NIST"))
        for lv in LEVELS:
            if lv not in CHEN[elem]:
                continue
            c = CHEN[elem][lv]
            for j in range(6):
                v0, v1 = get(res, elem, 0, j, lv), get(res, elem, 1, j, lv)
                b = c["b"][j]
                g = c.get("g", [None] * 6)[j]
                A("%-11s %5s %12s %12s %12s %12s %12s" % (
                    lv if j == 0 else "", LAMLAB[j],
                    "%12.4f" % v0 if v0 else "n/a",
                    "%12.4f" % b if b else "-",
                    "%12.4f" % v1 if v1 else "n/a",
                    "%12.4f" % g if g else "-",
                    "%12.4f" % c["NIST"] if j == 0 else ""))
        A("")

    # ---- 2. screening shifts -------------------------------------------
    A("\n2. SCREENING SHIFTS  dE(lam) = E_exc(lam) - E_exc(inf)\n")
    A("   The e-n-only shift is compared with Chen 'a' (MCDF, same physical model).")
    dev = []
    for elem in ION_ORDER:
        A("--- He-like %s" % ION_LABEL[elem])
        A("%-11s %5s %11s %11s %11s   %8s" %
          ("level", "lam", "this(ee=0)", "Chen a", "Chen b", "this/a"))
        for lv in LEVELS:
            if lv not in CHEN[elem]:
                continue
            c = CHEN[elem][lv]
            base = get(res, elem, 0, 0, lv)
            for j in range(1, 6):
                v = get(res, elem, 0, j, lv)
                da = c["a"][j] - c["a"][0]
                db = c["b"][j] - c["b"][0]
                if v is None or base is None:
                    A("%-11s %5s %11s %11.4f %11.4f   %8s" %
                      (lv if j == 1 else "", LAMLAB[j], "n/a", da, db, "-"))
                    continue
                do = v - base
                dev.append((do - da) / abs(da))
                A("%-11s %5s %11.4f %11.4f %11.4f   %8.2f" %
                  (lv if j == 1 else "", LAMLAB[j], do, da, db, do / da))
        A("")
    if dev:
        dev.sort()
        med = dev[len(dev) // 2]
        within = sum(1 for t in dev if abs(t) < 0.05)
        A("   Summary vs Chen 'a':  n=%d   median deviation %+.1f%%   within 5%%: %d/%d"
          % (len(dev), 100 * med, within, len(dev)))

    # ---- 3. kappa^2 scaling --------------------------------------------
    A("\n\n3. SCALING TEST   shift(lam=10) / shift(lam=20)   [theory: 4.00]\n")
    A("%-6s %-11s %10s %10s %10s" % ("ion", "level", "this work", "Chen a", "Chen b"))
    for elem in ION_ORDER:
        for lv in LEVELS:
            if lv not in CHEN[elem]:
                continue
            c = CHEN[elem][lv]
            base = get(res, elem, 0, 0, lv)
            v20, v10 = get(res, elem, 0, 1, lv), get(res, elem, 0, 3, lv)
            mine = "%10.2f" % ((v10 - base) / (v20 - base)) if None not in (base, v20, v10) else "       n/a"
            A("%-6s %-11s %s %10.2f %10.2f" % (
                ION_LABEL[elem], lv, mine,
                (c["a"][3] - c["a"][0]) / (c["a"][1] - c["a"][0]),
                (c["b"][3] - c["b"][0]) / (c["b"][1] - c["b"][0])))

    # ---- 4. e-e screening increment ------------------------------------
    A("\n\n4. ELECTRON-ELECTRON SCREENING INCREMENT  E_exc(ee=1) - E_exc(ee=0)\n")
    A("   Compared with Chen (g - b) and (f - a) where the paper reports them.")
    A("%-6s %-11s %5s %12s %12s %12s" %
      ("ion", "level", "lam", "this work", "Chen g-b", "Chen f-a"))
    for elem in ION_ORDER:
        for lv in LEVELS:
            if lv not in CHEN[elem]:
                continue
            c = CHEN[elem][lv]
            has = "g" in c
            for j in range(1, 6):
                v0, v1 = get(res, elem, 0, j, lv), get(res, elem, 1, j, lv)
                mine = "%12.5f" % (v1 - v0) if None not in (v0, v1) else "         n/a"
                gb = "%12.5f" % (c["g"][j] - c["b"][j]) if has and c["g"][j] else "           -"
                fa = "%12.5f" % (c["f"][j] - c["a"][j]) if has and c["f"][j] else "           -"
                A("%-6s %-11s %5s %s %s %s" %
                  (ION_LABEL[elem] if j == 1 else "", lv if j == 1 else "",
                   LAMLAB[j], mine, gb, fa))

    # ---- 5. total energies ---------------------------------------------
    A("\n\n5. GROUND-STATE TOTAL ENERGIES (a.u.)\n")
    A("%-6s %5s %16s %16s" % ("ion", "lam", "E0 (ee=0)", "E0 (ee=1)"))
    for elem in ION_ORDER:
        for j in range(6):
            t = lambda ee: res.get("%s_ee%d_L%s" % (
                elem, ee, ("%g" % LAMBDAS[j]).replace(".", "p").replace("+", "")))
            r0, r1 = t(0), t(1)
            A("%-6s %5s %16s %16s" % (
                ION_LABEL[elem] if j == 0 else "", LAMLAB[j],
                "%16.6f" % r0["E0"] if r0 else "n/a",
                "%16.6f" % r1["E0"] if r1 else "n/a"))

    # ---- 6. failures ----------------------------------------------------
    A("\n\n6. RUNS THAT DID NOT CONVERGE\n")
    if not fails:
        A("   none")
    for f in fails:
        A("   %-26s ee=%d lam=%-8g stage=%-6s  %s"
          % (f["tag"], f["ee"], f["lam"], f["stage"], f["error"]))
    A("")
    A("   A failure at 'n<=3' followed by no 'n<=2' entry means the n=2 levels were")
    A("   recovered with a reduced configuration set and the n=3 levels are absent.")
    A("   With the nucleus screened and e-e unscreened the asymptotic one-electron")
    A("   potential is +1/r, so only a finite number of states stay bound.")

    text = "\n".join(L)
    open(out, "w").write(text)
    print(text)
    return text


if __name__ == "__main__":
    force = "--force" in sys.argv
    print("Running %d FAC calculations ..." % (len(ION_ORDER) * len(MODELS) * len(LAMBDAS)))
    data = sweep(force=force)
    report(data)
