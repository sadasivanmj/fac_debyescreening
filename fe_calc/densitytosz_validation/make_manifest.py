#!/usr/bin/env python3
"""Write results/experiment_manifest.json: exact commits, the upstream jmaxk fix, the A/B source difference, build and
compiler provenance, library hashes and input definitions of the DensityToSZ experiment.  Fails if the two variants differ
anywhere in the FAC sources other than the intended line.

Usage: python make_manifest.py [--repo-a PATH] [--repo-ref PATH] [--repo-neg PATH]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "shared"))
sys.path.insert(0, str(HERE.parent / "deprince_2020_fe_ii_viii"))
import fac_run_support as sup      # noqa: E402
import configurations as cfg       # noqa: E402

UPSTREAM_FIX = "fb32c16a304a089522ce6ba67ca199527b27519c"


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _changed(diff: str) -> list[str]:
    """Only the added/removed lines (hunk offsets and blob ids differ when a patch is cherry-picked onto another base)."""
    return [l for l in diff.splitlines() if l[:1] in "+-" and not l.startswith(("+++", "---"))]


def build_of(repo: Path) -> dict:
    info = json.loads((repo / "build" / "BUILD_INFO.json").read_text()) if (repo / "build" / "BUILD_INFO.json").exists() else {}
    so = next(iter(sorted((repo / "build").glob("lib.*/pfac/fac*.so"))), None)
    return {"repository": str(repo), "head": git(repo, "rev-parse", "HEAD"), "branch": git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
            "fac_source_commit": git(repo, "log", "-1", "--format=%H", "--", "faclib", "python", "setup.py"),
            "build_info": info, "extension": str(so) if so else None, "extension_sha256": sha(so) if so else None,
            "libfac_a_sha256": sha(repo / "libfac.a") if (repo / "libfac.a").exists() else None,
            "workaround_flag_in_build": "-ftrivial-auto-var-init=zero" in info.get("cflags", ""),
            "radial_c_has_jmaxk_fix": "jmaxk = 0;" in (repo / "faclib" / "radial.c").read_text()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=HERE / "results")
    ap.add_argument("--repo-a", type=Path, default=sup.REPO_ROOT.parent / "fac_densitytosz_variant_a")
    ap.add_argument("--repo-ref", type=Path, default=sup.REPO_ROOT.parent / "fac_upstream_reference")
    ap.add_argument("--repo-neg", type=Path, default=sup.REPO_ROOT.parent / "fac_negative_control")
    args = ap.parse_args()
    B, A = sup.REPO_ROOT, args.repo_a.resolve()

    # upstream fix: hash, author, parent, diff (from the fetched official history, if present locally)
    up = {"hash": UPSTREAM_FIX}
    try:
        up.update({"author": git(B, "show", "-s", "--format=%an <%ae>", UPSTREAM_FIX), "parent": git(B, "show", "-s", "--format=%P", UPSTREAM_FIX),
                   "subject": git(B, "show", "-s", "--format=%s", UPSTREAM_FIX), "date": git(B, "show", "-s", "--format=%ad", UPSTREAM_FIX),
                   "diff": git(B, "show", "--format=", UPSTREAM_FIX), "official_repository": "https://github.com/flexible-atomic-code/fac"})
        up["is_ancestor_of_B_head"] = subprocess.run(["git", "-C", str(B), "merge-base", "--is-ancestor", UPSTREAM_FIX, "HEAD"]).returncode == 0
    except subprocess.CalledProcessError:
        up["note"] = "official history not present in this clone (fetch https://github.com/flexible-atomic-code/fac)"
    cherry = git(B, "log", "-1", "--format=%H", "--author=Ming Feng Gu", "--fixed-strings", "--grep=jmaxk fix (#851)")
    up["local_cherry_pick"] = {"hash": cherry, "author": git(B, "show", "-s", "--format=%an", cherry) if cherry else None,
                               "changed_lines_identical_to_upstream": bool(cherry) and _changed(git(B, "show", "--format=", cherry)) == _changed(up.get("diff", "")),
                               "note": "cherry-picked (new hash, upstream author and message preserved); the upstream hash is not an ancestor"}

    # source difference between the variants (FAC sources only)
    changed = git(B, "diff", "--name-only", "HEAD", git(A, "rev-parse", "HEAD"), "--", "faclib", "python", "setup.py").splitlines()
    diff = git(B, "diff", "HEAD", git(A, "rev-parse", "HEAD"), "--", "faclib", "python", "setup.py")
    ok = changed == ["faclib/orbital.c"] and diff.count("\n-  double mu") == 1 and diff.count("\n+  double mu = 0.0;") == 1
    if not ok:
        print("A and B differ in more than the intended line:", changed)
        return 1

    res = args.results_dir
    summary = json.loads((res / "potential_validation_summary.json").read_text())
    manifest = {
        "experiment": "DensityToSZ A/B: Coulomb (A) versus screened Yukawa (B) SCF direct potential, Fe II - Fe IV, gref",
        "upstream_jmaxk_fix": up,
        "variant_B": {**build_of(B), "description": "production screened DensityToSZ (v2); jmaxk fix cherry-picked; no workaround flag"},
        "variant_A": {**build_of(A), "description": "B + one line: double mu = 0.0 in DensityToSZ (Coulomb direct potential); screened Slater integrals and nuclear screening unchanged"},
        "pristine_upstream_reference": build_of(args.repo_ref.resolve()) if args.repo_ref.exists() else None,
        "negative_control_prefix_source": build_of(args.repo_neg.resolve()) if args.repo_neg.exists() else None,
        "source_difference_A_vs_B": {"files": changed, "diff": diff},
        "inputs": {"ions": cfg.ION_ORDER[:3], "configurations": {i: cfg.IONS[i][1] for i in cfg.ION_ORDER[:3]},
                   "orbital_model": "gref (OptimizeRadial on the ground configuration only; references added in Structure)",
                   "mu_values": [0.0, 0.1, 0.25], "plasma": {"density_A3": 1e-6, "ups": 0.0, "mu_step": 0.05, "hxs_iter0": 8,
                                                              "debye_mode": 1, "ee_screen": 1},
                   "calculator": "fe_calc/deprince_2020_fe_ii_viii/calculate_energies.py (unchanged physics; --direct-potential coulomb only changes the expected asymptotic charge)"},
        "potential_validation": {"all_passed": summary["all_passed"], "failures": summary["failures"], "tolerances": summary["tolerances"],
                                 "grid": summary["grid"]},
        "harness_executables_sha256": {p.name: sha(p) for p in sorted((res / "_bin").glob("harness_*"))},
        "host": {"platform": platform.platform(), "python": sys.version.split()[0], "numpy": __import__("numpy").__version__},
        "result_files_sha256": {p.name: sha(p) for p in sorted(res.glob("*.csv")) + sorted(res.glob("*.json")) if p.name != "experiment_manifest.json"},
    }
    (res / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2))
    print("wrote", res / "experiment_manifest.json")
    print("A vs B FAC-source difference:", changed, "| workaround flag in B/A builds:",
          manifest["variant_B"]["workaround_flag_in_build"], manifest["variant_A"]["workaround_flag_in_build"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
