#!/usr/bin/env python3
"""Run all Paper B tests and write results/validation_report.md (every pass / fail / skip)."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "shared"))
from testing_support import run_and_report            # noqa: E402

sys.exit(run_and_report(HERE, HERE.parent / "results" / "validation_report.md",
                        "Validation report: Deprince 2019, Fe XVII - Fe XXVI (reduced ground-only / reference-only models)"))
