# Archive

These are old scripts, patches and experiments that are no longer part of the
normal workflow. They are kept because they explain how the current code and
results came about. Nothing here is used by the active scripts.

**Path warning:** some archived scripts still use their old relative paths,
for example `fe_calc/` or the repository root. Run them only from a copy of
the original layout.

| Item | Original purpose | Why archived | Replacement | Status |
|---|---|---|---|---|
| `experiments/early-fe-calculations/fe_II_energy.py`, `fe_III_energy.py` | First Fe II / Fe III μ = 0 reference-configuration runs (output `FeII_mu0.lev`, `FeIII_mu0.lev`) | Superseded | `research/debye-screening/cases/fe/fe_ii_reference.py`, `fe_iii_reference.py` (identical μ = 0 energies, e.g. Fe III −34542.4125 eV) | Superseded |
| `experiments/early-fe-calculations/fe_iii_diag.py` + `debug_FeIII.lev*` | Fe III screened diagnostic (λ = 10 a₀, nps = 0.01) | It sets `ee_screen 0`, so it never tested the screened SCF. The stored `debug_FeIII.lev` came from an older build and is about 1650 eV off | v2 tests in `research/debye-screening/tests/debye_scf/` | Historical |
| `experiments/early-fe-calculations/ip_potential_feii.py`, `screened_ip_potential_feii.py` + `FeII_*`, `FeIII_*` level files | Early free and screened Fe II IP scripts | Superseded by the validated screening chain and the reference suite | reference suite | Superseded |
| `experiments/early-fe-calculations/fe_ip_sfac.py` + `fe_ip_sfac_runs/` | Fe II IP with the `sfac` executable; contains the original μ = 0.1 `RadialBound` failure (`error.pot`) | It used an `sfac` built from the old `ee-screening-scf` branch | reference suite | Historical (first evidence of the bracket failure) |
| `experiments/deprince-t2-sfac/` | Earlier Table 2 study of Fe II–VIII with `sfac` (`ee-screening-scf` + local patches): driver, CSF counter, report, result CSVs, logs | Different code base; superseded by the validated v2 chain. Its CSF counter is correct (the shared-seen bug was in a later script) | reference suite; `research/full-active-space/` | Historical. Raw run folders (≈ 380 MB) stay untracked in the main working tree under `fe_calc/deprince_t2/` |
| `experiments/ee-screening-scf/` | Parallel line of work on branch `ee-screening-scf`: an earlier screened `DensityToSZ`, with Xie 2012, Fritzsche 2006, Chen 2018 Kr³⁴⁺ and C4+ studies | Its `DensityToSZ` is superseded by v2 (agreement ≤ 0.2 meV on the Fe reference cases); its results were not recomputed with v2 | `faclib/orbital.c` (v2); `research/debye-screening/` | Historical; merged for history. Its C code is not active (see its README) |
| `experiments/v2-development/` | Development of the v2 screened direct potential: reports, patch generator, kernel unit tests against analytic Yukawa potentials | Development finished; the result is in `faclib/orbital.c` | `research/debye-screening/patches/01_*.patch` | Historical (documents the 10⁻⁹ accuracy test and the v1 NaN cases) |
| `experiments/sturm-node-count/` + `patches/sturm_node_count_diagnostic.patch` | Diagnostic instrumentation of `RadialBound` (full-grid Sturm count) and the Fe IX μ = 1.5–3 scan | Diagnostic only, never meant for production | `research/debye-screening/patches/03_debye_radialbound_bracket.patch` | Historical. Branch `diagnostic/sturm-node-count` |
| `experiments/warm-start-continuation/` | Warm-start and cold-start experiments that found the exchange-free first SCF iteration | The cause was found and fixed by `hxs_iter0` | `research/debye-screening/patches/02_*.patch`, `results/hxs-iter0-validation/` | Historical. Branch `experiment/warm-start-continuation` |
| `experiments/audits-2026-10-12/` | Audit notes on the Fe II–VIII IPs and the Fe III anomaly, plus the change history of the investigation | Analysis notes, not code | `PROJECT_STATUS.md`, `docs/numerical-validation.md` | Historical |
| `scripts/verify_ip_deprince.py` | Reference-space IP driver (all J, both parities), with `--hxs-iter0`; produced the regression references | Replaced by simpler per-ion scripts that reproduce its energies exactly | `research/debye-screening/cases/fe/` | Superseded. Its first CSF counter had a shared-seen-set bug: it undercounted the active space. The bug affected only the CSF counts, not the energies, and was fixed on 2026-10-13 |
| `patches/fac-debye-nuclear-screening.patch` | One-hunk sign fix for nuclear-only Debye screening (Chen model) | No longer applies; the code in `faclib/orbital.c` has moved on | current `faclib/orbital.c` | Superseded |
| `patches/accidental_git_diff_2026-10-12.diff` | A `git diff` saved by mistake as a file called `pbcopy` (the macOS clipboard command) | Not a program | — | Historical |
| `source/orbital.c.orig` | An older working copy of `faclib/orbital.c` (earlier `ZPS` formula) | It was never used by the build | git history of `faclib/orbital.c` | Historical |
