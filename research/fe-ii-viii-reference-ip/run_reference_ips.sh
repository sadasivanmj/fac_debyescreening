#!/bin/bash
# Run Fe II - Fe IX (reference configurations only) at mu = 0, 0.1, 0.25 for both
# orbital models, one fresh Python process per run, sequentially.
# Usage:  ./run_reference_ips.sh [RUN_ID]      (refuses to reuse an existing RUN_ID)
# Output: results/RUN_ID/runs/<ion>/<model>_mu<mu>/ and results/RUN_ID/logs/.
# The locally built pfac in the repository's build/ folder is used (nothing is installed).

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
echo "ion,model,mu,exit_code,seconds" > "$OUT/run_summary.csv"

failed=0
cd "$HERE"
for ion in fe_ii fe_iii fe_iv fe_v fe_vi fe_vii fe_viii fe_ix; do
    for model in gref ref; do
        for mu in 0.0 0.1 0.25; do
            echo "running $ion $model mu=$mu"
            start=$SECONDS
            "$PYTHON" -u fe_reference_run.py "$ion" "$model" "$mu" "$OUT/runs/$ion/${model}_mu${mu}" \
                > "$OUT/logs/${ion}_${model}_mu${mu}.log" 2>&1
            rc=$?
            echo "$ion,$model,$mu,$rc,$((SECONDS-start))" >> "$OUT/run_summary.csv"
            if [ $rc -ne 0 ]; then
                echo "  FAILED (exit code $rc), see $OUT/logs/${ion}_${model}_mu${mu}.log"
                failed=$((failed+1))
            fi
        done
    done
done
echo "finished: $failed failed run(s); results in $OUT"
exit $((failed > 0))
