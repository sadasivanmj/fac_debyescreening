# Fe XVIII (9 electrons): ground configuration 1s2 2s2 2p5.
# Usage: python fe_xviii_ground.py --mu-ne 0.1 --mu-ee 0.0 --outdir results/runs/<ion>/ne010_ee000
from fe_ground_common import run_ground

run_ground("Fe XVIII", "1s2 2s2 2p5")
