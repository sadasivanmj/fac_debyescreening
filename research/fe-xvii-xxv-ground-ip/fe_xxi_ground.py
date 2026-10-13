# Fe XXI (6 electrons): ground configuration 1s2 2s2 2p2.
# Usage: python fe_xxi_ground.py --mu-ne 0.1 --mu-ee 0.0 --outdir results/runs/<ion>/ne010_ee000
from fe_ground_common import run_ground

run_ground("Fe XXI", "1s2 2s2 2p2")
