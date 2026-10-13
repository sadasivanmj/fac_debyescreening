"""Helpers for the test suites (standard library only; never imports FAC)."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import time
import traceback
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fac_run_support as sup   # noqa: E402


def write_synthetic_results(results_dir: Path, records: list[dict]) -> None:
    """Write energies.csv and one manifest per record, in the layout the calculators use.

    Each record needs: ion, ion_order, electron_count, configuration_model,
    orbital_optimization_model, case, case_order, mu_ne, mu_ee, total_energy; optional
    overrides for status, fac_source_commit and build_provenance.
    """
    results_dir = Path(results_dir)
    rows = []
    for r in records:
        rec = {"paper": "synthetic", "charge_state": 26 - r["electron_count"], "energy_unit": "eV",
               "total_energy_hartree_converted": r["total_energy"] / sup.HARTREE_EV,
               "ground_level_identifier": "ILEV=0", "convergence_status": "converged", "status": "ok",
               "exit_code": 0, "elapsed_time": 1.0, "fac_executable_or_module": "build/x.so",
               "fac_source_commit": "a" * 40, "build_provenance": "synthetic-build"}
        rec.update(r)
        slug = rec["ion"].lower().replace(" ", "_")
        rel = f"runs/{slug}/{rec['configuration_model']}/{rec['orbital_optimization_model']}/{rec['case']}"
        rec["run_directory"] = rel
        rec["log_path"] = rel + "/run.log"
        d = results_dir / rel
        d.mkdir(parents=True, exist_ok=True)
        (d / "manifest.json").write_text(json.dumps(rec))
        rows.append(rec)
    with open(results_dir / "energies.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sup.ENERGY_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run_script(script: Path, *args: str, env_extra: dict | None = None, timeout: float = 600) -> subprocess.CompletedProcess:
    import os
    env = dict(os.environ)
    env.update(env_extra or {})
    return subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True, timeout=timeout, env=env)


def poisoned_pfac_dir() -> tempfile.TemporaryDirectory:
    """A directory with a fake 'pfac' package that raises on import, to prove that a
    program never imports FAC."""
    tmp = tempfile.TemporaryDirectory()
    pkg = Path(tmp.name) / "pfac"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("raise ImportError('pfac must not be imported by the analysis program')\n")
    return tmp


class RecordingResult(unittest.TextTestResult):
    """Collects pass / fail / skip per test for the validation report."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.outcomes: list[tuple[str, str, str, float]] = []
        self._t0 = 0.0

    def startTest(self, test):
        self._t0 = time.monotonic()
        super().startTest(test)

    def _add(self, test, status, detail=""):
        self.outcomes.append((test.id().split(".", 1)[-1], status, detail, time.monotonic() - self._t0))

    def addSuccess(self, test):
        super().addSuccess(test)
        self._add(test, "PASS")

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self._add(test, "FAIL", str(err[1]).splitlines()[0][:300] if str(err[1]) else "")

    def addError(self, test, err):
        super().addError(test, err)
        self._add(test, "FAIL", "".join(traceback.format_exception_only(err[0], err[1])).strip()[:300])

    def addSkip(self, test, reason):
        super().addSkip(test, reason)
        self._add(test, "SKIP", reason)


def run_and_report(test_dir: Path, report_path: Path, title: str, extra_dirs: list[Path] | None = None) -> int:
    """Run every test_*.py in test_dir with unittest and write a compact Markdown report."""
    for d in [test_dir, *(extra_dirs or [])]:
        sys.path.insert(0, str(d))
    suite = unittest.defaultTestLoader.discover(str(test_dir), pattern="test_*.py", top_level_dir=str(test_dir))
    runner = unittest.TextTestRunner(resultclass=RecordingResult, verbosity=1, stream=sys.stderr)
    result = runner.run(suite)
    g = sup.git_info()
    lines = [f"# {title}", "", f"Run on commit `{g['commit'][:10]}` (FAC source `{g['fac_source_commit'][:10]}`), branch `{g['branch']}`.",
             f"Build: {sup.build_provenance_summary()}", "",
             "| Test | Result | Seconds | Note |", "|---|---|---|---|"]
    for name, status, detail, secs in result.outcomes:
        lines.append(f"| `{name}` | {status} | {secs:.1f} | {detail.replace('|', '/')} |")
    n_pass = sum(1 for o in result.outcomes if o[1] == "PASS")
    n_fail = sum(1 for o in result.outcomes if o[1] == "FAIL")
    n_skip = sum(1 for o in result.outcomes if o[1] == "SKIP")
    lines += ["", f"**{n_pass} passed, {n_fail} failed, {n_skip} skipped.**", ""]
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)
    Path(report_path).write_text("\n".join(lines))
    print(f"validation report: {report_path} ({n_pass} passed, {n_fail} failed, {n_skip} skipped)")
    return 0 if n_fail == 0 else 1
