# Fe XX (7 electrons): ground configuration 1s2 2s2 2p3.
# Usage: python fe_xx_ground.py --mu-ne 0.1 --mu-ee 0.0 --outdir results/runs/<ion>/ne010_ee000
from fe_ground_common import run_ground

run_ground("Fe XX", "1s2 2s2 2p3")
