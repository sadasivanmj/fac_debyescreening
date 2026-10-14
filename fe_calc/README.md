# fe_calc: reproducible FAC energies and IPs for the Deprince et al. iron-ion studies

Two independent studies, each with a **separate energy calculator** and a **separate IP-analysis
program**:

| directory | paper | ions |
|---|---|---|
| [`deprince_2020_fe_ii_viii/`](deprince_2020_fe_ii_viii/README.md) | Deprince et al., A&A 643, A57 (2020), DOI 10.1051/0004-6361/202038659 | Fe II - Fe VIII (+ Fe IX, needed for IP(Fe VIII)) |
| [`deprince_2019_fe_xvii_xxv/`](deprince_2019_fe_xvii_xxv/README.md) | Deprince et al., A&A 626, A83 (2019), DOI 10.1051/0004-6361/201935444 | Fe XVII - Fe XXV (+ Fe XXVI, needed for IP(Fe XXV)) |

**Scope: reduced configuration models only.** The papers' MCDF restricted active spaces (all single and
double excitations) are deliberately NOT implemented. Nothing here is an exact reproduction of the
published MCDF results.

## Order of use

```bash
fe_calc/shared/build_pfac.sh; source fe_calc/shared/use_local_pfac.sh   # build (once per source change) and select the local build
cd fe_calc/deprince_2020_fe_ii_viii                    # or deprince_2019_fe_xvii_xxv
python calculate_energies.py --check-only              # configuration validation (no FAC needed)
python calculate_energies.py --preflight               # build / option / provenance checks
python calculate_energies.py                           # energies only  -> results/energies.csv + runs/
python analyze_ionization_potentials.py                # IPs and literature comparison (no FAC, no calculations)
python tests/run_validation.py                         # tests -> results/validation_report.md
```

* `calculate_energies.py` computes and saves **total energies only**. `analyze_ionization_potentials.py`
  reads the saved `energies.csv` and run manifests; it never imports FAC, never launches a run and never
  recomputes or substitutes a missing energy (it exits with code 2 and lists every problem).
* `shared/` holds only generic helpers (provenance, FAC output parsing, run directories, CSV/statistics,
  configuration syntax checks, build script). Model definitions, reference data, run conditions, the FAC
  sequences and the comparison logic are separate in each paper directory.

## Build (one workflow, this checkout only)

```bash
cd /Users/sadasivanmj/Documents/projects/fac_debyescreening
conda activate fac_mac                       # the environment you run FAC from
fe_calc/shared/build_pfac.sh                 # builds build/lib.*/pfac for THIS python; installs nothing
source fe_calc/shared/use_local_pfac.sh      # every new shell: puts this build first on PYTHONPATH
python -c "import pfac; print(pfac.__file__)"   # must print .../fac_debyescreening/build/...
```

`build_pfac.sh` runs `./configure` and `make pfac` only (never `make install` or `pip install`), with
ordinary flags (`-O2 -fPIC`; macOS: Apple clang and the arm64 `--build/--host` arguments are chosen
automatically, because conda's cross-compiler wrapper in `CC` breaks `./configure`). It writes
`build/BUILD_INFO.json` (commit, branch, dirty flag, flags, compiler, Python and conda env, sha256 of the
extension, which Debye features are compiled in, and where `import pfac` resolves). Every energy record
stores the build provenance and `--preflight` checks that the extension matches `BUILD_INFO.json`.

**The `jmaxk` defect is fixed in the source.** `PotentialHX1()` (`faclib/radial.c`) returned
`Max(jmax, jmaxk)` with `jmaxk` unassigned when `mps == 1` (Debye mode): stack garbage, which crashed
screened runs (SIGBUS) on macOS arm64. Upstream fixed it (M. F. Gu, FAC #851: `jmaxk = 0;` before the
`mps` block; this repository has it as `603c0e53`, patch-id identical to upstream `fb32c16a`). The earlier
workaround flag `-ftrivial-auto-var-init=zero` is **no longer used and the build script refuses it**
(it would hide other uninitialised variables); the script also refuses to build if the fix is missing.
Evidence: `fe_calc/densitytosz_validation/` (negative control with the unfixed source and identical flags
versus the fixed source).

**A plain `import pfac` may still find an older installation in `site-packages`** (the one in `fac_mac`
predates the Debye SCF fixes). Always `source fe_calc/shared/use_local_pfac.sh` first, and never
`pip install` / `make install` over it.

## Validated screening code required

`--preflight` verifies that the build contains: the stable Debye direct potential `DensityToSZ` v2
(commit `00eb229b`, screening the SCF direct potential, not only the Slater integrals), `radial:hxs_iter0`
(`0650edb1`), the Debye `RadialBound` bracket (`8c8b1bfb`) and the Fe XXVI one-electron nuclear-screening fix
(`7f05cb9d`, same patch as `fbc69734`). At run time every screened run checks that
r(Vc+U)+ZPS tends to 0 when both interactions are screened and to N-1 when only the nucleus is, which proves
that the electron-electron direct potential, not just the Slater integrals, is screened.

## Output layout, git and provenance

Each paper directory has `results/` with `runs/<ion>/<model>/<case>/` (`input.json`, `run.log`, `levels.txt`,
`final.pot`, `worker_result.json`, `manifest.json`), `energies.csv`, `ionization_potentials.csv`,
`comparison_summary.md` and `validation_report.md`. Completed runs are never overwritten (`--overwrite` moves the
old directory aside); failed runs are kept as `*.failed-attempt-*` directories. Tracked in Git: manifests,
inputs, `energies.csv`, the IP tables and reports. Not tracked (regenerable, bulky): `run.log` of successful
runs, `levels.txt`, `final.pot`, `levels.b`, the build. Failure logs are named `run.FAILED.log` and are tracked.
`fac_source_commit` in every record is the last commit that touched `faclib/`, `python/` or `setup.py`.
