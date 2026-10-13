# Fe XXIII (4 electrons): ground configuration 1s2 2s2.
# Usage: python fe_xxiii_ground.py --mu-ne 0.1 --mu-ee 0.0 --outdir results/runs/<ion>/ne010_ee000
from fe_ground_common import run_ground

run_ground("Fe XXIII", "1s2 2s2")
