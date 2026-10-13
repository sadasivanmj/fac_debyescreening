#!/bin/bash
# Run the Fe XVII-XXVI ground-configuration calculations for the four Table 4
# screening cases, one fresh Python process per run, sequentially.
#
# Usage:  ./run_fe_ground_configs.sh [RUN_ID]
# Everything goes to results/RUN_ID/ (runs/, logs/); RUN_ID defaults to a time
# stamp and the script refuses to reuse an existing one.
# Set PYTHON to choose the interpreter; the locally built pfac in the
# repository's build/ folder is used (nothing is installed).

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PYTHON="${PYTHON:-python3}"
RUN_ID="${1:-run_$(date +%Y%m%d_%H%M%S)}"
OUT="$HERE/results/$RUN_ID"

BUILD="$(ls -d "$REPO"/build/lib.* 2>/dev/null | head -1)"
[ -n "$BUILD" ] || { echo "no build/lib.* in $REPO; run 'make pfac' first"; exit 1; }
export PYTHONPATH="$BUILD"
"$PYTHON" -c "from pfac import fac; print('using', fac.__file__)" || exit 1

[ -e "$OUT" ] && { echo "$OUT already exists; choose another RUN_ID"; exit 1; }
mkdir -p "$OUT/runs" "$OUT/logs"
echo "ion,case,mu_ne,mu_ee,exit_code,seconds,outdir" > "$OUT/run_summary.csv"

# case name : mu_ne : mu_ee   (Deprince et al. 2019, Table 4)
CASES="A:0.10:0.00 B:0.10:0.10 C:0.25:0.00 D:0.25:0.25"
IONS="xvii xviii xix xx xxi xxii xxiii xxiv xxv xxvi"

failed=0
for case in $CASES; do
    IFS=: read name mu_ne mu_ee <<< "$case"
    for ion in $IONS; do
        folder="$OUT/runs/fe_${ion}/case_${name}"
        log="$OUT/logs/fe_${ion}_case_${name}.log"
        mkdir -p "$(dirname "$folder")"
        echo "running Fe ${ion} case ${name} (mu_ne=${mu_ne}, mu_ee=${mu_ee})"
        start=$SECONDS
        "$PYTHON" -u "$HERE/fe_${ion}_ground.py" --mu-ne "$mu_ne" --mu-ee "$mu_ee" \
            --outdir "$folder" > "$log" 2>&1
        rc=$?
        echo "fe_${ion},${name},${mu_ne},${mu_ee},${rc},$((SECONDS-start)),${folder}" >> "$OUT/run_summary.csv"
        if [ $rc -ne 0 ]; then
            echo "  FAILED (exit code $rc), see $log"
            failed=$((failed+1))
        fi
    done
done
echo "finished: $failed failed run(s); results in $OUT"
exit $((failed > 0))
