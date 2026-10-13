#!/bin/bash
# Run the 18 reference-configuration calculations (3 ions x 3 mu x 2 models)
# with the FAC build of this repository, then collect the summary tables.
# Usage:  research/debye-screening/cases/fe/run_reference_suite.sh
# The scripts stay in cases/fe; the runs and the two CSV tables are written to
# research/debye-screening/results/reference-configurations.
# Nothing is installed: Python imports pfac from the repository's build folder
# (build it first with ./configure and make pfac in the repository root).

CASES="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$CASES/../../../.." && pwd)"
OUT="$REPO/research/debye-screening/results/reference-configurations"
export PYTHONPATH="$REPO/build/lib.linux-x86_64-cpython-310"
python3 -c "from pfac import fac; print('using', fac.__file__)" || exit 1

cd "$OUT"
if [ -e reference_runs ]; then
    echo "$OUT/reference_runs already exists; move it away before rerunning"
    exit 1
fi
mkdir reference_runs

for ion in fe_ii fe_iii fe_iv; do
    for model in ref gref; do
        for mu in 0.0 0.1 0.25; do
            folder="reference_runs/${ion}/${model}_mu${mu}"
            mkdir -p "reference_runs/${ion}"
            echo "running ${ion} ${model} mu=${mu}"
            python3 -u "$CASES/${ion}_reference.py" "$model" "$mu" "$folder" > "${folder}.log" 2>&1
            echo "exit code $?" >> "${folder}.log"
        done
    done
done

python3 "$CASES/collect_reference_results.py"
