"""Support code for the energy-calculation programs (standard library only).

This module never imports FAC (pfac).  The FAC calls live in the per-paper
``calculate_energies.py`` worker, so importing this module is always safe.

Contents: provenance of the build and source, parsing of FAC output files, log
scanning, safe run directories (nothing completed is overwritten), subprocess
execution with a timeout, and the aggregated energies.csv.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HARTREE_EV = 27.211386          # eV per Hartree, used for every conversion
REPO_ROOT = Path(__file__).resolve().parents[2]

# Tracked files that every build rewrites; they do not make the tree "dirty".
BUILD_ARTIFACTS = {"faclib/f2c.mod"}

# FAC prints this when an orbital optimisation hits its iteration limit.
NOT_CONVERGED_MARKER = "Maximum iteration reached in OptimizeRadial"
WARNING_PATTERN = re.compile(r"warning|error|abort|invalid|fail|nan\b|inf\b", re.IGNORECASE)

ENERGY_FIELDS = [
    "paper", "ion", "charge_state", "electron_count", "configuration_model",
    "orbital_optimization_model", "case", "mu_ne", "mu_ee",
    "total_energy", "energy_unit", "total_energy_hartree_converted",
    "ground_level_identifier", "convergence_status", "status", "exit_code",
    "elapsed_time", "fac_executable_or_module", "fac_source_commit",
    "build_provenance", "run_directory", "log_path",
]


# --------------------------------------------------------------------- provenance

def sha256_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(*args: str, strip: bool = True) -> str:
    try:
        out = subprocess.check_output(["git", "-C", str(REPO_ROOT), *args], text=True,
                                      stderr=subprocess.DEVNULL)
        return out.strip() if strip else out
    except Exception:
        return ""


def git_info() -> dict:
    """Commit, branch and the list of modified tracked files (build artifacts excluded)."""
    changed = [line[3:] for line in _git("status", "--porcelain", "--untracked-files=no", strip=False).splitlines()]
    changed = [p for p in changed if p not in BUILD_ARTIFACTS]
    # fac_source_commit: last commit that touched the FAC C/Python sources, so that committing
    # analysis scripts later does not change the identity of the FAC source that made an energy.
    source_changed = [line for line in _git("status", "--porcelain", "--untracked-files=no", "--", "faclib", "python", strip=False).splitlines()
                      if line[3:] not in BUILD_ARTIFACTS]
    return {"fac_source_commit": (_git("log", "-1", "--format=%H", "--", "faclib", "python", "setup.py") or "unknown")
            + ("+modified" if source_changed else ""),
            "commit": _git("rev-parse", "HEAD") or "unknown",
            "branch": _git("rev-parse", "--abbrev-ref", "HEAD") or "unknown",
            "modified_tracked_files": changed}


def find_build_dir(repo: Path = REPO_ROOT) -> Path | None:
    dirs = sorted((repo / "build").glob("lib.*"))
    return dirs[0] if dirs else None


def find_extension(repo: Path = REPO_ROOT) -> Path | None:
    build = find_build_dir(repo)
    if build is None:
        return None
    found = sorted((build / "pfac").glob("fac*.so"))
    return found[0] if found else None


def build_info(repo: Path = REPO_ROOT) -> dict:
    """What build/BUILD_INFO.json (written by shared/build_pfac.sh) says, or {}."""
    p = repo / "build" / "BUILD_INFO.json"
    return json.loads(p.read_text()) if p.exists() else {}


def build_provenance_summary(repo: Path = REPO_ROOT) -> str:
    """One compact string stored in every energy record."""
    so = find_extension(repo)
    info = build_info(repo)
    if so is None:
        return "no-build"
    sha = sha256_of(so)
    flags = info.get("cflags", "unknown-flags")
    verified = "BUILD_INFO-matches" if info.get("extension_sha256") == sha else "BUILD_INFO-MISSING-or-MISMATCH"
    return f"cflags=[{flags}]; cc={info.get('cc_version', 'unknown')}; ext_sha256={sha[:16]}; {verified}"


def preflight_checks(repo: Path = REPO_ROOT, results_dir: Path | None = None,
                     require_workaround: bool = True) -> list[tuple[str, str, str]]:
    """Return (status, check, detail) rows; status is PASS, WARN or FAIL."""
    rows: list[tuple[str, str, str]] = []

    def add(ok: bool, name: str, detail: str, warn_only: bool = False) -> None:
        rows.append(("PASS" if ok else ("WARN" if warn_only else "FAIL"), name, detail))

    so = find_extension(repo)
    add(so is not None, "extension built", str(so) if so else "no build/lib.*/pfac/fac*.so; run fe_calc/shared/build_pfac.sh")
    blob = so.read_bytes() if so else b""
    for opt in ("radial:ee_screen", "radial:hxs_iter0", "orbital:debye_mode"):
        add(opt.encode() in blob, f"option known to build: {opt}", "string found in extension")
    add(b"invalid Debye screening parameter in DensityToSZ" in blob,
        "stable Debye DensityToSZ (v2) compiled in", "v2 error string found in extension")

    def source(name: str) -> str:
        p = repo / "faclib" / name
        return p.read_text(errors="replace") if p.exists() else ""
    orbital, radial = source("orbital.c"), source("radial.c")
    add("_ee_screened_sz" in orbital, "source has screened direct potential (v2)", "faclib/orbital.c: _ee_screened_sz")
    add("debye_emax_eps" in orbital, "source has Debye RadialBound bracket", "faclib/orbital.c: debye_emax_eps")
    add("hxs_iter0" in radial, "source has radial:hxs_iter0", "faclib/radial.c")
    add("if (potential->mps != 1) potential->ZPS[m] = 0;" in radial,
        "source has the Fe XXVI one-electron nuclear-screening fix", "faclib/radial.c PotentialHX")

    if so is not None:
        newest = max((repo / "faclib" / n).stat().st_mtime for n in ("orbital.c", "radial.c") if (repo / "faclib" / n).exists())
        add(so.stat().st_mtime >= newest, "build not older than sources", "extension newer than orbital.c/radial.c")
    info = build_info(repo)
    add(bool(info), "build/BUILD_INFO.json present", "written by shared/build_pfac.sh", warn_only=True)
    if so is not None and info:
        add(info.get("extension_sha256") == sha256_of(so), "BUILD_INFO matches extension", "sha256 equal")

    # The jmaxk defect: in PotentialHX1, jmaxk is declared without initialiser and is only
    # assigned when mps != 1, but Max(jmax, jmaxk) is returned.
    defect = re.search(r"int\s+jmax,\s*jmaxk,\s*m;", radial) is not None
    flag = "-ftrivial-auto-var-init=zero" in info.get("cflags", "")
    if defect:
        add(flag or not require_workaround, "jmaxk uninitialised-variable defect",
            "source defect PRESENT (not fixed); workaround flag -ftrivial-auto-var-init=zero "
            + ("is in the build" if flag else "is MISSING from the build: screened runs may crash"))
    else:
        add(True, "jmaxk uninitialised-variable defect", "no longer present in source")

    if results_dir is not None:
        try:
            results_dir.mkdir(parents=True, exist_ok=True)
            probe = results_dir / ".write_test"
            probe.write_text("x")
            probe.unlink()
            add(True, "results directory writable", str(results_dir))
        except OSError as e:
            add(False, "results directory writable", str(e))
    g = git_info()
    add(not g["modified_tracked_files"], "working tree clean (tracked files)",
        f"commit {g['commit'][:10]} on {g['branch']}" + (f"; modified: {g['modified_tracked_files']}" if g["modified_tracked_files"] else ""),
        warn_only=True)
    return rows


# ------------------------------------------------------------------ FAC output files

def parse_level_table(path: Path) -> dict:
    """Parse the text table written by FAC's PrintTable.

    Level energies are in eV relative to E0, the energy of the lowest level
    (block, E0).  Returns nele, nlev (from the header), e0_ev and the level rows.
    """
    out = {"nele": None, "nlev_header": None, "e0_ev": None, "levels": []}
    for line in Path(path).read_text().splitlines():
        if line.startswith("NELE"):
            out["nele"] = int(line.split("=")[1])
        elif line.startswith("NLEV"):
            out["nlev_header"] = int(line.split("=")[1])
        elif line.startswith("E0"):
            out["e0_ev"] = float(line.split(",")[1])
        else:
            p = line.split()
            # level row: ILEV IBASE ENERGY P VNL 2J <labels...>
            if len(p) >= 6 and p[0].isdigit() and "E" in p[2]:
                out["levels"].append({"ilev": int(p[0]), "energy_rel_ev": float(p[2]),
                                      "parity": int(p[3]), "vnl": int(p[4]), "two_j": int(p[5]),
                                      "labels": " ".join(p[6:])})
    return out


def ground_level(table: dict) -> dict | None:
    """The lowest level of the table (all J and both parities are in one table)."""
    return min(table["levels"], key=lambda lv: lv["energy_rel_ev"]) if table["levels"] else None


def parse_potential(path: Path) -> dict:
    """Header values of final.pot and the asymptotic charge r(Vc+U)+ZPS at the last point."""
    header, last = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            parts = line[1:].split("=")
            if len(parts) == 2:
                header[parts[0].strip()] = parts[1].strip()
        elif line.strip():
            last = line.split()
    return {"mps": int(header["mps"]), "dps": float(header["dps"]), "ips": int(header["ips"]),
            "N": float(header["N"]),
            "asymptotic_charge": float(last[4]) + float(last[5]) + float(last[17]),
            "zps_last": float(last[17])}


def scan_log(text: str) -> dict:
    """Convergence marker and suspicious lines in a run log."""
    suspicious = [l.strip() for l in text.splitlines()
                  if WARNING_PATTERN.search(l) and "FAC extension" not in l][:20]
    return {"max_iteration_reached": NOT_CONVERGED_MARKER in text, "suspicious_lines": suspicious}


# -------------------------------------------------------------- run directory handling

def timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def read_manifest(run_dir: Path) -> dict | None:
    p = Path(run_dir) / "manifest.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except json.JSONDecodeError:
        return None


def archive_previous(run_dir: Path, reason: str) -> Path | None:
    """Move an existing run directory aside (never delete) and return its new path."""
    run_dir = Path(run_dir)
    if not run_dir.exists():
        return None
    target = run_dir.with_name(f"{run_dir.name}.{reason}-{timestamp()}")
    shutil.move(str(run_dir), str(target))
    return target


def execute_worker(cmd: list[str], run_dir: Path, timeout_s: float) -> tuple[int, float, bool]:
    """Run one worker process in its own fresh interpreter.

    stdout and stderr go to run_dir/run.log.  Returns (exit_code, seconds, timed_out).
    A fresh process per run matters: FAC keeps orbitals and potentials in global memory.
    """
    log = Path(run_dir) / "run.log"
    start = time.monotonic()
    with open(log, "w") as f:
        try:
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=timeout_s, cwd=run_dir)
            return proc.returncode, time.monotonic() - start, False
        except subprocess.TimeoutExpired:
            f.write(f"\nTIMEOUT after {timeout_s} s\n")
            return -9, time.monotonic() - start, True


def finalize_log_name(run_dir: Path, ok: bool) -> Path:
    """Successful logs stay run.log (git-ignored, verbose); failed ones become run.FAILED.log
    so the diagnostics can be tracked."""
    log = Path(run_dir) / "run.log"
    if not ok and log.exists():
        failed = Path(run_dir) / "run.FAILED.log"
        log.rename(failed)
        return failed
    return log


# ------------------------------------------------------------------------ aggregation

def collect_manifests(results_dir: Path) -> list[dict]:
    """Every run manifest below results_dir/runs (current ones only, archived copies skipped)."""
    out = []
    for m in sorted(Path(results_dir, "runs").glob("*/*/*/manifest.json")):
        data = json.loads(m.read_text())
        out.append(data)
    return out


def write_energies_csv(results_dir: Path) -> Path:
    rows = collect_manifests(results_dir)
    rows.sort(key=lambda r: (r["ion_order"], r["configuration_model"], r["orbital_optimization_model"], r["case_order"]))
    path = Path(results_dir) / "energies.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ENERGY_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def is_finite(x) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def python_executable() -> str:
    return sys.executable
