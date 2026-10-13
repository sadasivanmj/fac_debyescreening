# DensityToSZ() validation and A/B experiment (Fe II, Fe III, Fe IV)

Tests whether FAC's screened SCF direct (Hartree) potential in `DensityToSZ()` is correct, and measures its effect on
Fe II -> III and Fe III -> IV ionization potentials. **Experimental; nothing here is merged into the production branch.**
Read `results/experiment_summary.md` first, then `results/densitytosz_validation.md`.

## Worktrees and builds (all `-O2 -fPIC`, no `-ftrivial-auto-var-init=zero`)

| worktree | branch / commit | role |
|---|---|---|
| `fac_densitytosz_experiment` | `experiment/densitytosz-validation` | variant **B** (production screened kernel) + the upstream `jmaxk` fix `fb32c16a` (cherry-picked, M. F. Gu) + this directory |
| `fac_densitytosz_variant_a` | `experiment/densitytosz-variant-a` | variant **A** = B + one line (`double mu = 0.0;` in `DensityToSZ`) |
| `fac_upstream_reference` | detached `fb32c16a` | pristine upstream, only for the bit-identity test of the Coulomb path |
| `fac_negative_control` | detached `bbbed5da` | pre-fix source, only for the `jmaxk` negative control |

Build each with `fe_calc/shared/build_pfac.sh` (env `FAC_CFLAGS='-O2 -fPIC'`, `CC`, `FAC_CONFIGURE_ARGS`); it writes `build/BUILD_INFO.json`.
`FE_CALC_REPO=<worktree>` makes the framework use another worktree's build and record its commit.

## Order of use

```bash
python potential_validation.py          # real DensityToSZ via the C harness vs independent references -> potential_validation.csv/.json
# A/B energies (18 runs): each variant from its own build
python ../deprince_2020_fe_ii_viii/calculate_energies.py --ions fe_ii fe_iii fe_iv --models gref --direct-potential screened --results-dir results/variant_B
FE_CALC_REPO=../../../fac_densitytosz_variant_a python ../deprince_2020_fe_ii_viii/calculate_energies.py --ions fe_ii fe_iii fe_iv --models gref --direct-potential coulomb --results-dir results/variant_A
python variant_a_diagnosis.py [--csi 0 --results-dir results/variant_A_diagnosis_csi0]   # why A fails; small-mu ladder
python first_order_check.py             # exact first-order Debye-Hueckel check of both builds
python analyze_densitytosz_experiment.py --allow-incomplete   # reads saved energies only (no FAC) -> energies.csv, ip_comparison.csv, plots/
python make_manifest.py                 # -> experiment_manifest.json (fails if A and B differ by more than the intended line)
python write_reports.py                 # -> densitytosz_validation.md, experiment_summary.md
python -m unittest discover -s tests    # unit tests of the evaluators and of the analysis (no FAC)
```

`harness/densitytosz_harness.c` + `build_harness.sh` compile a small C program against a build's `libfac.a`; it initialises FAC like
`import pfac` (`InitFac0/1`, `SetAtom`, `ReinitRadial`, `SetRadialGrid`, `SetOrbitalRGrid`, plasma via `PlasmaScreen` + `InitializePS`) and calls the real `DensityToSZ()`.

Results not tracked in Git (regenerable): `_bin/`, `_work/`, level tables, potentials, successful-run logs. Failure logs (`run.FAILED.log`) are tracked.
The first-run potential validation is kept as `potential_validation_summary_first_run.json` (see Section 9 of `densitytosz_validation.md`).
