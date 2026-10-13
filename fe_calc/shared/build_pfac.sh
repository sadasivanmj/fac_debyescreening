#!/bin/bash
# Build the local PFAC extension (nothing is installed) and record how it was built.
#
# Usage:  fe_calc/shared/build_pfac.sh            (from anywhere inside the repository)
# Env:    CC                 C compiler for ./configure (default: cc)
#         FAC_CONFIGURE_ARGS extra arguments for ./configure
#                            (macOS arm64: --build=aarch64-apple-darwin20.0.0 --host=aarch64-apple-darwin20.0.0)
#         FAC_CFLAGS         flags written to the Makefiles
#                            (default: -O0 -g -fPIC -ftrivial-auto-var-init=zero)
#
# Why -ftrivial-auto-var-init=zero: PotentialHX1() in faclib/radial.c returns
# Max(jmax, jmaxk), but jmaxk is only assigned when mps != 1, so in Debye mode it is
# uninitialised stack garbage.  On macOS arm64 this crashes some screened runs
# (SIGBUS in PotentialHX).  The flag starts every local at 0, which gives
# Max(jmax, 0) = jmax as intended.  It is a build WORKAROUND; the source defect is
# NOT fixed.  See fe_calc/README.md.
#
# Writes build/BUILD_INFO.json (flags, compiler, commit, sha256 of the extension).
# Never runs make install or pip install.

set -e
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
FAC_CFLAGS="${FAC_CFLAGS:--O0 -g -fPIC -ftrivial-auto-var-init=zero}"

CFLAGS='-O2 -fPIC' ./configure $FAC_CONFIGURE_ARGS > build_configure.log 2>&1
# ./configure writes "CFLAGS = -O2 -fPIC" into these Makefiles; replace it.
for mk in faclib/Makefile python/Makefile; do
    sed -i.bak "s|^CFLAGS = -O2 -fPIC.*|CFLAGS = $FAC_CFLAGS|" "$mk" && rm -f "$mk.bak"
    grep -q "^CFLAGS = $FAC_CFLAGS" "$mk" || { echo "could not set CFLAGS in $mk"; exit 1; }
done
make clean > /dev/null 2>&1 || true
rm -rf build
make pfac > build_make.log 2>&1

SO="$(ls build/lib.*/pfac/fac*.so | head -1)"
python3 - "$REPO" "$SO" "$FAC_CFLAGS" <<'PY'
import hashlib, json, platform, subprocess, sys, datetime
repo, so, flags = sys.argv[1:4]
def run(*cmd):
    try:
        return subprocess.check_output(cmd, cwd=repo, text=True, stderr=subprocess.STDOUT).strip()
    except Exception as e:
        return f"unavailable ({e})"
cc = run("sh", "-c", "grep -m1 '^CC = ' faclib/Makefile | cut -d= -f2-").strip()
info = {
    "built_utc": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z",
    "git_commit": run("git", "rev-parse", "HEAD"),
    "git_dirty": bool(run("git", "status", "--porcelain", "--untracked-files=no")),
    "cflags": flags,
    "cc": cc,
    "cc_version": run("sh", "-c", f"{cc} --version | head -1"),
    "platform": platform.platform(), "python": sys.version.split()[0],
    "extension": so,
    "extension_sha256": hashlib.sha256(open(f"{repo}/{so}", "rb").read()).hexdigest(),
    "workaround_jmaxk_zero_init": "-ftrivial-auto-var-init=zero" in flags,
}
json.dump(info, open(f"{repo}/build/BUILD_INFO.json", "w"), indent=2)
print(json.dumps(info, indent=2))
PY
