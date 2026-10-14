#!/bin/bash
# Source this file (do not execute it) to make `import pfac` resolve to the build of THIS checkout:
#     cd /Users/sadasivanmj/Documents/projects/fac_debyescreening
#     source fe_calc/shared/use_local_pfac.sh
# It puts build/lib.*/ first on PYTHONPATH and prints what Python will import.  Nothing is installed.
_fac_repo="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"
_fac_lib="$(ls -d "$_fac_repo"/build/lib.*/ 2>/dev/null | head -1)"
if [ -z "$_fac_lib" ]; then
    echo "no build found in $_fac_repo/build; run fe_calc/shared/build_pfac.sh first"
else
    export PYTHONPATH="${_fac_lib%/}${PYTHONPATH:+:$PYTHONPATH}"
    python - <<'PY'
import os, pfac
print("pfac ->", pfac.__file__)
PY
fi
unset _fac_repo _fac_lib
