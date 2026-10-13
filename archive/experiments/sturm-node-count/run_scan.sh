#!/bin/sh
# Fe IX (1s2 2s2 2p6 3s2 3p6) Debye scan with the STURM_DIAG build of this branch.
# usage: tests/sturm_diag/run_scan.sh <pfac_build_lib_dir> <out_dir>   (nothing is installed)
# Inputs are those of the original failure (run_diag.py case F_FeIX_mu3_nps1e-9):
# PlasmaScreen(zps=8, nps=1e-9, tps from mu, m=1, ups=0), debye_mode 1, ee_screen 1.
set -e
L=$(cd "$1" && pwd); OUT=$2; HERE=$(cd "$(dirname "$0")/.." && pwd)
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT"; exit 1; }
for m in 1.5 2.0 2.25 2.5 2.75 3.0; do
  mkdir -p "$OUT/mu$m"
  (cd "$OUT/mu$m" && PYTHONPATH=$L python3 "$HERE/debye_scf/run_case.py" \
    "{\"configs\":[\"1s2 2s2 2p6 3s2 3p6\"],\"mu\":$m,\"nps\":1e-9,\"zps\":8.0,\"debye_mode\":1,\"ee_screen\":1,\"expect_ext\":\"$L\"}" \
    > out.log 2>&1 || true)
done
