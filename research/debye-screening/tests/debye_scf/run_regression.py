#!/usr/bin/env python3
"""Unscreened (mu = 0) regression of the screened DensityToSZ against a baseline build.

usage: python3 run_regression.py <baseline_pfac_dir> <test_pfac_dir> <out_dir>

Each *_pfac_dir is a directory containing the built `pfac` package (e.g.
build/lib.linux-x86_64-cpython-310 of a `make pfac` build); nothing is installed.
Runs the same cases with both builds (run_case.py checks which extension was
imported) and compares E0, every level energy and the final.pot arrays.
Exits non-zero if any compared value differs.
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = {
    "He_free": dict(configs=["1s2"]),
    "FeIX_free": dict(configs=["1s2 2s2 2p6 3s2 3p6"]),
    "FeIII_free_CE": dict(configs=["1s2 2s2 2p6 3s2 3p6 3d6"], config_energy=True),
    # ee_screen on but no plasma: EEScreenMu() = 0, so the Coulomb path must run
    "FeIX_free_ee1": dict(configs=["1s2 2s2 2p6 3s2 3p6"], debye_mode=1, ee_screen=1),
}


def levels(p):
    t = Path(p).read_text().splitlines()
    i = next(k for k, l in enumerate(t) if l.startswith("E0"))
    rows = [l.split() for l in t[i + 1:] if l.strip() and l.split()[0].isdigit()]
    return float(t[i].split(",")[1]), [float(r[2]) for r in rows]


def pot_arrays(p):
    rows = [l.split() for l in open(p) if l.strip() and not l.startswith("#")]
    return [[float(x) for x in r] for r in rows]


def run(name, spec, lib, wd):
    wd.mkdir(parents=True)
    spec = dict(spec, expect_ext=str(lib))
    env = {"PYTHONPATH": str(lib), "PATH": "/usr/bin:/bin", "HOME": str(Path.home())}
    with open(wd / "out.log", "w") as log:
        rc = subprocess.run([sys.executable, str(HERE / "run_case.py"), json.dumps(spec)],
                            cwd=wd, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    imp = next((l.split()[1] for l in open(wd / "out.log") if l.startswith("IMPORTED")), None)
    return rc, imp


def main():
    base, test, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve(), Path(sys.argv[3])
    if out.exists():
        sys.exit(f"refusing to overwrite {out}")
    res, fail = [], False
    for name, spec in CASES.items():
        r = {"case": name}
        for tag, lib in (("base", base), ("test", test)):
            r[tag + "_rc"], r[tag + "_imported"] = run(name, spec, lib, out / tag / name)
        if r["base_rc"] or r["test_rc"]:
            fail = True
            res.append(r)
            continue
        e0b, lb = levels(out / "base" / name / "e.lev")
        e0t, lt = levels(out / "test" / name / "e.lev")
        pb = pot_arrays(out / "base" / name / "final.pot")
        pt = pot_arrays(out / "test" / name / "final.pot")
        r.update(E0_base=e0b, E0_test=e0t, nlev=(len(lb), len(lt)),
                 max_abs_E0=abs(e0b - e0t),
                 max_abs_lev=max((abs(a - b) for a, b in zip(lb, lt)), default=0.0),
                 pot_shape=(len(pb), len(pb[0]) if pb else 0),
                 max_abs_pot=max(abs(a - b) for ra, rb in zip(pb, pt) for a, b in zip(ra, rb)))
        mx = max((max(abs(a), abs(b)) for ra, rb in zip(pb, pt) for a, b in zip(ra, rb)), default=1)
        r["max_rel_pot"] = r["max_abs_pot"] / mx if mx else 0.0
        if (r["max_abs_E0"] or r["max_abs_lev"] or r["max_abs_pot"] or len(lb) != len(lt)
                or len(pb) != len(pt)):
            fail = True
        res.append(r)
        print(json.dumps(r), flush=True)
    (out / "regression.json").write_text(json.dumps(res, indent=1))
    print("REGRESSION", "FAIL" if fail else "PASS (bitwise identical)")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
