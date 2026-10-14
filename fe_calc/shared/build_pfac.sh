#!/bin/bash
# Build the local PFAC extension from THIS checkout (nothing is installed) and record how it was built.
#
# Usage:  fe_calc/shared/build_pfac.sh            (from anywhere inside the repository, in the Python
#                                                  environment you will run FAC from, e.g. `conda activate fac_mac`)
# Env:    FAC_CC             C compiler (default: Apple clang on macOS, else the inherited CC).  The inherited
#                            CC is ignored on macOS because the conda cross-compiler wrapper that
#                            conda sets there (arm64-apple-darwin20.0.0-clang) makes ./configure fail.
#         FAC_CONFIGURE_ARGS extra ./configure arguments (default on macOS arm64:
#                            --build=aarch64-apple-darwin20.0.0 --host=aarch64-apple-darwin20.0.0)
#         FAC_CFLAGS         flags written to the Makefiles (default: -O2 -fPIC)
#
# The extension is built for the first `python` on PATH, so it matches the active conda environment.
#
# jmaxk: PotentialHX1() in faclib/radial.c used to return Max(jmax, jmaxk) with jmaxk uninitialised when
# mps == 1 (Debye mode).  Upstream fixed it (M. F. Gu, FAC #851, `jmaxk = 0;`; patch-id identical to the
# fix in this repository).  The build therefore uses ordinary flags.  The script REFUSES to build with
# -ftrivial-auto-var-init=zero (it would mask any other uninitialised variable) and refuses to build if the
# source fix is missing.
#
# Writes build/BUILD_INFO.json, then imports the fresh extension with PYTHONPATH set and checks that
# `pfac` resolves to this checkout (and warns if a plain `import pfac` resolves elsewhere).
# Never runs make install or pip install.

set -e
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
FAC_CFLAGS="${FAC_CFLAGS:--O2 -fPIC}"

case "$FAC_CFLAGS" in
  *-ftrivial-auto-var-init*) echo "refusing: FAC_CFLAGS contains -ftrivial-auto-var-init (it masks uninitialised variables); the jmaxk defect is fixed in the source"; exit 1;;
esac
grep -Eq 'jmaxk[[:space:]]*=[[:space:]]*0;' faclib/radial.c || { echo "refusing: faclib/radial.c lacks the upstream 'jmaxk = 0;' fix (FAC #851)"; exit 1; }

if [ -n "${FAC_CC:-}" ]; then
    export CC="$FAC_CC"
elif [ "$(uname -s)" = "Darwin" ] && [ -x /Library/Developer/CommandLineTools/usr/bin/clang ]; then
    export CC=/Library/Developer/CommandLineTools/usr/bin/clang
fi
if [ -z "${FAC_CONFIGURE_ARGS+x}" ] && [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ]; then
    FAC_CONFIGURE_ARGS="--build=aarch64-apple-darwin20.0.0 --host=aarch64-apple-darwin20.0.0"
fi

fail() { echo "build failed; last lines of $1:"; tail -15 "$1"; exit 1; }
CFLAGS="$FAC_CFLAGS" ./configure ${FAC_CONFIGURE_ARGS:-} > build_configure.log 2>&1 || fail build_configure.log
# ./configure writes "CFLAGS = -O2 -fPIC" into these Makefiles; replace it when other flags were requested.
for mk in faclib/Makefile python/Makefile; do
    sed -i.bak "s|^CFLAGS = .*|CFLAGS = $FAC_CFLAGS|" "$mk" && rm -f "$mk.bak"
    grep -q "^CFLAGS = $FAC_CFLAGS\$" "$mk" || { echo "could not set CFLAGS in $mk"; exit 1; }
done
make clean > /dev/null 2>&1 || true
rm -rf build
make pfac > build_make.log 2>&1 || fail build_make.log

SO="$(ls build/lib.*/pfac/fac*.so | head -1)"
LIBDIR="$(dirname "$(dirname "$SO")")"
python - "$REPO" "$SO" "$LIBDIR" "$FAC_CFLAGS" <<'PY'
import hashlib, json, os, platform, subprocess, sys, datetime
repo, so, libdir, flags = sys.argv[1:5]
def run(*cmd, **kw):
    try:
        return subprocess.check_output(cmd, cwd=repo, text=True, stderr=subprocess.STDOUT, **kw).strip()
    except Exception as e:
        return f"unavailable ({e})"
cc = run("sh", "-c", "grep -m1 '^CC = ' faclib/Makefile | cut -d= -f2-").strip()
blob = open(f"{repo}/{so}", "rb").read()
radial = open(f"{repo}/faclib/radial.c", errors="replace").read()
info = {
    "built_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "git_commit": run("git", "rev-parse", "HEAD"),
    "git_branch": run("git", "rev-parse", "--abbrev-ref", "HEAD"),
    "git_dirty": bool(run("git", "status", "--porcelain", "--untracked-files=no")),
    "cflags": flags,
    "cc": cc,
    "cc_version": run("sh", "-c", f"{cc} --version | head -1"),
    "platform": platform.platform(), "python": sys.version.split()[0],
    "python_executable": sys.executable, "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
    "extension": so,
    "extension_sha256": hashlib.sha256(blob).hexdigest(),
    "workaround_jmaxk_zero_init": "-ftrivial-auto-var-init=zero" in flags,
    "jmaxk_source_fix": "jmaxk = 0;" in radial,
    "features_in_extension": {k: k.encode() in blob for k in
        ("radial:ee_screen", "radial:hxs_iter0", "orbital:debye_mode", "No Debye bracket",
         "invalid Debye screening parameter in DensityToSZ")},
}
# import the fresh build in a clean subprocess and confirm where `pfac` comes from
env = dict(os.environ, PYTHONPATH=f"{repo}/{libdir}")
code = "import pfac, pfac.fac as f; print(pfac.__file__); print(f.__file__)"
out = subprocess.check_output([sys.executable, "-c", code], env=env, cwd="/", text=True).split()
info["import_with_pythonpath"] = out
info["import_resolves_to_this_checkout"] = all(os.path.realpath(p).startswith(os.path.realpath(repo) + os.sep) for p in out)
plain = subprocess.run([sys.executable, "-c", "import pfac; print(pfac.__file__)"], cwd="/", text=True,
                       capture_output=True, env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"})
info["plain_import_without_pythonpath"] = plain.stdout.strip() or plain.stderr.strip().splitlines()[-1:]
json.dump(info, open(f"{repo}/build/BUILD_INFO.json", "w"), indent=2)
print(json.dumps(info, indent=2))
if not info["import_resolves_to_this_checkout"]:
    sys.exit("ERROR: the fresh build does not import from this checkout")
if not all(info["features_in_extension"].values()):
    sys.exit("ERROR: the build lacks required Debye features: " + str([k for k, v in info["features_in_extension"].items() if not v]))
if isinstance(info["plain_import_without_pythonpath"], str) and not os.path.realpath(info["plain_import_without_pythonpath"]).startswith(os.path.realpath(repo)):
    print("\nNOTE: a plain `import pfac` resolves to another installation:\n  " + info["plain_import_without_pythonpath"] +
          "\n  Use the local build:  source fe_calc/shared/use_local_pfac.sh   (never pip/make install over it)")
PY
echo
echo "To use this build in the current shell:  source $REPO/fe_calc/shared/use_local_pfac.sh"
echo "tracked-file changes after the build (should be empty):"; git status --porcelain --untracked-files=no
