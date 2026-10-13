# Screened Fe XXVI ground state and the Fe XXV ionization potentials

Completes Table 4 of Deprince et al., A&A 626, A83 (2019) for Fe XXV. Follows
`../fe-xvii-xxv-ground-ip/` (minimal ground-configuration benchmark).
Branch `research/fe-xxvi-screened-ground`, based on `research/fe-xvii-xxv-ground-ip` (`31fa9738`).

## Why FAC failed to give the screened Fe XXVI energy

The calculation did not crash: it silently ran **unscreened**.

* The runtime values were right: `mps = 1`, `dps = 10.0` (μ_ne = 0.1), 1 electron, `1s1`.
* `SetPotentialPS` computed the Debye nuclear term `ZPS(r) = Z(r)(1 − e^{−r/λ})` correctly.
* Then `PotentialHX` (`faclib/radial.c`) took its early return for fewer than
  two electrons (`jmax <= 0 || N+nqf < 1+EPS5`) and set `ZPS = 0` for all r. For
  N = 1 the potential is just `Vc + ZPS/r` (`SetPotentialVT`), so the nuclear
  screening was lost. For N ≥ 2 that branch is not taken.
* Evidence (`results/run_20261013_diag/unpatched_baseline/`): the potential file has
  ZPS = 0 at every radius where Z(1−e^{−μr}) is 12.5–26, and r·V → −26 (free ion)
  instead of 0. The small energy changes seen before (+0.12 eV at μ = 0.1) are
  therefore not the screened nuclear potential.
* Not the cause: the exchange ramp (`hxs_iter0`), the `RadialBound` bracket, the
  screening kernel, or the unit conversion. The radial solver never failed.

## Fix (C level, one line)

`faclib/radial.c`, in that early return: `if (potential->mps != 1) potential->ZPS[m] = 0;`
In Debye mode ZPS does not depend on the electrons, so it must be kept. Other
plasma models and N ≥ 2 are untouched. It is committed separately from the results.
The three validated patches (`DensityToSZ`, `hxs_iter0`, `RadialBound`) are unchanged.

## Regression

* Unscreened Fe XXVI (μ = 0): `final.pot` and `levels.txt` identical to the unpatched
  build (`cmp`/`diff`).
* Full rerun of all 40 ground-state runs on the patched build
  (`../fe-xvii-xxv-ground-ip/results/run_20261013_patched_rerun/`): all 40 valid, and
  the 36 energies that were valid before are reproduced with max |ΔE| = 0 eV.

## Results

Fe XXVI `1s1`, built as for the other ions (free-ion SCF, steps of 0.05 in μ,
`hxs_iter0 = 8`, `-O0` build, see the earlier README):

| μ_ne | E₀ (eV) | E₀ (Ha) | μ_ee = 0 vs μ_ee = μ_ne |
|---|---|---|---|
| 0 | −9277.69186 | | (control, unscreened) |
| 0.10 | −9207.14345 | −338.3563 | identical (spread 0 eV) |
| 0.25 | −9102.07055 | −334.4949 | identical (spread 0 eV) |

Cross-check with an independent solver for the same potential
(`../fe-xvii-xxv-ground-ip/supplement_fe_xxvi_dirac_yukawa.py`: point nucleus,
Dirac equation, V = −Z e^{−μr}/r; not FAC, `fac_vs_dirac_yukawa.txt`):
FAC − Dirac is 4.0530, 4.0530 and 4.0526 eV at μ = 0, 0.1, 0.25 (finite nucleus + QED
in FAC), and the screening shifts agree: 70.5484 vs 70.5485 eV (0.1) and
175.6213 vs 175.6218 eV (0.25).

### Fe XXV ionization potentials (`results/run_20261013_fe_xxv_ips/fe_xxv_ip_table4.csv`)

IP = E₀(Fe XXVI) − E₀(Fe XXV), Fe XXV energies from the validated `run_20261013_minimal`.

| Case | μ_ne | μ_ee | IP calc (eV) | Table 4 (eV) | diff (eV) | diff (%) |
|---|---|---|---|---|---|---|
| A | 0.10 | 0.00 | 8756.1057 | 8766.28 | −10.17 | −0.116 |
| B | 0.10 | 0.10 | 8758.8152 | 8768.90 | −10.08 | −0.115 |
| C | 0.25 | 0.00 | 8651.0649 | 8661.14 | −10.08 | −0.116 |
| D | 0.25 | 0.25 | 8657.7977 | 8667.86 | −10.06 | −0.116 |

Screening-induced changes: electron–electron screening raises the IP by 2.71 eV
(μ = 0.1; Table 4: 2.62) and 6.73 eV (μ = 0.25; Table 4: 6.72). Raising μ_ne from
0.1 to 0.25 lowers it by 105.04 eV (A→C; Table 4: 105.14) and 101.02 eV (B→D; 101.04).
The trends are reproduced.

### The constant offset

The 10 eV difference is the same in all four cases, so it does not come from the
screening. Extrapolating the two μ_ne values linearly to μ = 0 gives about 8826.1 eV
for these results (FAC's own unscreened IP is 8826.65 eV) and about 8836.4 eV for
Table 4. The paper's unscreened Fe XXV value is not available here, so I cannot say
where its extra ~10 eV comes from (single configuration here vs MCDF with correlation,
Breit and QED there, or the paper's own Fe XXV/XXVI energies). Nothing was tuned.
The other ions differ from Table 4 by at most 6.9 eV, so Fe XXV is the largest
discrepancy in the set.

Same limits as before: a minimal ground-configuration benchmark, not a RAS/MCDF reproduction.

## Files

| file | purpose |
|---|---|
| `diag_fe_xxvi.py` | minimal reproducer (Tests A/B/C); prints μ used, ZPS(r) read from the potential, E₀ |
| `complete_fe_xxv_ips.py` | consistency check, μ_ee independence, the four IPs |
| `results/run_20261013_diag/{unpatched_baseline,patched}/mu*/` | `diag.json`, `levels.txt`, log for μ = 0, 0.1, 0.25 |
| `results/run_20261013_fe_xxv_ips/` | IP table, independent solver output, FAC-vs-Dirac comparison |

Extension used: `build/lib.macosx-11.0-arm64-cpython-310/pfac/fac.cpython-310-darwin.so`
inside this worktree (built with `make pfac` at `-O0`; nothing installed); the scripts
refuse any other `pfac`. Commit, imports and options are in each `diag.json`/`manifest.json`.
