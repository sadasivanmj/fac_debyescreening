# Debye-screening research

Scripts, tests and results for FAC with Debye-Hückel screening.

- **Screened interactions:** the electron-nucleus and the electron-electron
  interaction are both screened as e^(−μr)/r.
- **Screening parameter:** μ = 1/λ_D, in atomic units.
- **Where the C changes live:** [`../../docs/architecture.md`](../../docs/architecture.md).
- **What has been validated:** [`../../docs/numerical-validation.md`](../../docs/numerical-validation.md).

## Folders

| folder | contents |
|---|---|
| `patches/` | The C changes as patch files, in the order they were made (see below). |
| `scripts/` | `debye_fac.py`: helper that converts λ_D into `PlasmaScreen` arguments. `chen2018_reproduce.py`: Chen et al. (2018) He-like ions with `sfac`. |
| `tests/` | Small checks: `verify_yukawa.py`, `test_debye_c4.py`, `test_debye_screening_minimal.py`, `debye_scf/` (μ = 0 regression of the screened direct potential against an unmodified build) and `warm_start/` (`hxs_iter0` continuation test). |
| `cases/fe/` | Fe II–IV reference-configuration calculations (one script per ion, run script, result collector). |
| `results/reference-configurations/` | Their CSV tables, the write-up `REFERENCE_CONFIGURATION_VALIDATION.md` and the run folders. |
| `results/hxs-iter0-validation/` | Logs showing the Fe II / Fe III ladders converge with `hxs_iter0 = 8` and fail without it. |
| `results/chen2018/` | Chen et al. reproduction: report, JSON results, `sfac` input files. |
| `results/c4/` | C4+ test output (`c_20.*`, λ_D = 20 a₀) and the earlier C4+ ground-state study (`c4_ground_results/`). |
| `references/` | Reference paper(s); see `references/README.md`. |

The patches in `patches/`, in order:

1. `fac-debye-ee-screening.patch` — screened Slater integrals (`GetYkDebye`).
2. `01_v2_screened_direct_potential_DensityToSZ.patch`
3. `02_hxs_iter0_exchange_ramp.patch`
4. `03_debye_radialbound_bracket.patch`

All four are already in `faclib/`. Do **not** apply them again.

## Build first (nothing is installed)

From the repository root:

```bash
./configure --prefix=$PWD/_noinstall
make pfac
export PYTHONPATH=$PWD/build/lib.linux-x86_64-cpython-310   # use this build, not a global pfac
python3 -c "from pfac import fac; print(fac.__file__)"     # check which FAC is imported
```

## Commands

Run all of these from the repository root, with `PYTHONPATH` set as above.

```bash
# Fe II–IV reference configurations: 18 runs, about 1 minute.  It refuses to
# overwrite results/reference-configurations/reference_runs, so move that
# folder away first if you want a fresh run.
research/debye-screening/cases/fe/run_reference_suite.sh

# mu = 0 regression of the screened direct potential: compares this build with
# an unmodified FAC build (build c4a016a in another worktree first)
python3 research/debye-screening/tests/debye_scf/run_regression.py \
    <baseline>/build/lib.linux-x86_64-cpython-310 build/lib.linux-x86_64-cpython-310 /tmp/regression_out

# potential-level Yukawa check on the stored C4+ potential (lambda_D = 20 a0)
python3 research/debye-screening/tests/verify_yukawa.py research/debye-screening/results/c4/c_20.pot 20

# C4+ nuclear-only test; writes c_<lambda>.* into the current folder
python3 research/debye-screening/tests/test_debye_c4.py 20

# Chen et al. (2018): needs an sfac binary built from this tree (SFAC=...).
# Outputs go into the folder it is started from.
cd research/debye-screening/results/chen2018 && SFAC=/path/to/sfac python3 ../../scripts/chen2018_reproduce.py
```

## Notes

- **Nuclear-only scripts are older.** `debye_fac.py`, `test_debye_c4.py` and
  `chen2018_reproduce.py` were written for the nuclear-only screening model
  of Chen et al.
  - The nuclear sign patch they mention, `fac-debye-nuclear-screening.patch`,
    has since been superseded by the code in `faclib/`.
  - That patch is kept in `archive/patches/`.
- **Earlier Fe IP scripts are archived.** These are the `fe_calc/` scripts
  and the Deprince `sfac` study; see `../../archive/README.md`.
