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
fe_calc/shared/build_pfac.sh                           # once: build the local PFAC extension (nothing is installed)
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

## Build and the `jmaxk` defect (important)

Use `shared/build_pfac.sh`. It runs `./configure` and `make pfac` only (never `make install` or `pip install`)
and writes `build/BUILD_INFO.json` (flags, compiler, commit, sha256 of the extension). Every energy record
stores the build provenance, and `--preflight` checks that the extension matches `BUILD_INFO.json`.

**Outstanding source defect (not fixed):** in `PotentialHX1()` (`faclib/radial.c`) the local `jmaxk` is only
assigned when `mps != 1`, but `Max(jmax, jmaxk)` is returned. In Debye mode (`mps = 1`) the result is
uninitialised stack garbage; with a garbage value of 0x55555555 `PotentialHX` reads far outside its arrays
and crashes (SIGBUS) on macOS arm64 (seen for Fe III `gref` at mu = 0.1 and Fe XXVI with `-O2`). On Linux the
garbage was harmless. The framework **does not change the C source**; it builds with the compiler flag
`-ftrivial-auto-var-init=zero` (plus `-O0 -g`), which zero-initialises locals so that `Max(jmax, 0) = jmax` as
intended. This is a workaround, not a fix: the defect is still in the source, the preflight reports it, and
a build without the flag fails the preflight. A source fix (initialise `jmaxk` at its declaration) should be a
separate, reviewed commit with its own regression test. Checks that support the workaround: all 18 Fe II-IV
energies agree with Linux-validated values, and 40 earlier Fe XVII-XXVI energies from builds without the flag
are reproduced exactly (max |dE| = 0 eV).

## Validated screening code required

`--preflight` verifies that the build contains: the stable Debye direct potential `DensityToSZ` v2
(commit `00eb229b`, screening the SCF direct potential, not only the Slater integrals), `radial:hxs_iter0`
(`0650edb1`), the Debye `RadialBound` bracket (`8c8b1bfb`) and the Fe XXVI one-electron nuclear-screening fix
(cherry-picked from `fbc69734` onto this branch). At run time every screened run checks that
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
