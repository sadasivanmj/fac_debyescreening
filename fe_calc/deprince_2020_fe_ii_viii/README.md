# Deprince et al. (2020): Fe II - Fe VIII ionization potentials, reduced reference-configuration model

* Paper: J. Deprince et al., *Plasma environment effects on K lines of astrophysical interest. IV. IPs, K thresholds,
  radiative rates, and Auger widths in Fe II-Fe VIII*, A&A 643, A57 (2020). Compared: **Table 2** (computed IPs at
  mu = 0, 0.1, 0.25 au) and the NIST isolated-ion IPs quoted in it. The local PDF copy is under
  `research/debye-screening/references/` (not committed; copyright).
* Fe IX (needed for IP(Fe VIII)): reference set from Paper III, A&A 635, A70 (2020), Table 1.

## Model

Only the listed reference configurations (`configurations.py`; `[nl]` = hole, closed core
`1s2 2s2 2p6 3s2 3p6`): Fe II 7, Fe III 10, Fe IV - Fe VIII 4 each, Fe IX 4 (`3p6, [3p]3d, [2p]3d, [1s]3d`).
**No single or double excitations are generated**, the paper's 12 268 - 85 798 CSFs are not reproduced.
Orbital models (same CI space, kept apart everywhere): `gref` = orbitals optimised on the ground configuration only
(primary), `ref` = orbitals optimised on the average of all reference configurations.
Screening: mu = 0 (plain unscreened FAC, `PlasmaScreen` never called), 0.1 and 0.25 with electron-nucleus and
electron-electron screening equal (`orbital:debye_mode 1`, `radial:ee_screen 1`, `radial:hxs_iter0 8`, density
1e-6 Angstrom^-3, T = 4 pi n / mu^2, mu switched on in steps of 0.05 after a free-ion SCF). 8 ions x 3 mu x 2 models = 48 runs.

## Commands

```bash
python calculate_energies.py --check-only | --preflight | --dry-run
python calculate_energies.py --ions fe_ii --models gref ref --mu 0 0.1     # a selection
python calculate_energies.py [--jobs 2] [--timeout 1800] [--overwrite]       # all 48 runs; resumable
python analyze_ionization_potentials.py                                      # -> ionization_potentials.csv, comparison_summary.md
python tests/run_validation.py                                               # -> results/validation_report.md
```

## Conventions

* FAC writes level energies in eV relative to `E0` (the energy of the lowest level); the total ground-state energy is
  `E0` + the lowest relative energy (0), in eV, negative. `energies.csv` stores it as `total_energy` with
  `energy_unit = eV` and a labelled converted `total_energy_hartree_converted` (1 Ha = 27.211386 eV).
* IP(q) = [E(q+1) - E(q)] / 27.211386 Ha x 27.211386 eV/Ha for the same model, mu and build (anything else is refused).
* Ground level: the lowest level of the single FAC level table (all J, both parities; `Structure` over all groups). The run
  checks NELE, finite E0, that the table is referenced to its lowest level and that the header NLEV equals the number of parsed
  levels, and records `ILEV/2J/parity` of the ground level and the (parity, 2J) blocks seen. It cannot prove independently that
  every symmetry allowed by the configurations was generated; FAC's default `Structure` covers all of them.

## Validation (`results/validation_report.md`: 27 passed, 0 failed, 0 skipped)

Configuration syntax/counts/duplicates; Fe II `ref` and `gref` at mu = 0 (distinct energies, identical configuration space and
986 levels); `mu = 0` path unscreened (`mps = -1`, free-ion asymptotic charge); Fe II at mu = 0.1 screened on both interactions
(asymptotic charge 0); Fe III with the same settings and the Fe II IP by independent subtraction (15.5203 / 10.2258 eV for `gref`
at mu = 0 / 0.1); the 18 Fe II-IV energies against the Linux-validated values (tolerance 1e-3 eV; platform and flags differ, so
bit-identity is not claimed); analysis on synthetic data (independent recomputation, missing/duplicated/failed/incompatible records,
no FAC import); resume, no-overwrite, timeout and retry; preflight.

## Results (all 48 runs converged; 42/42 IPs; see `comparison_summary.md`, `ionization_potentials.csv`)

| `gref` (primary) | n | mean abs (eV) | max abs (eV) |
|---|---|---|---|
| IP vs paper Table 2 | 21 | 1.93 | 7.15 (Fe III -> IV, mu = 0.25) |
| IP vs NIST (mu = 0) | 7 | 1.28 | 2.56 (Fe III) |
| screening shift vs paper | 14 | 0.58 | 2.31 |

`ref`: 2.85 / 16.39 eV vs paper, 2.99 / 11.00 eV vs NIST, shift 0.67 / 3.12 eV; the Fe III -> Fe IV `ref` value is an orbital-averaging
artefact (equal-weight average dominated by core-hole configurations). IPs fall with mu for all ions and both models; FAC is
about 1 eV below NIST for Fe V - VIII. The paper's own values also differ from NIST, so agreement with the paper and with NIST are
different tests. Likely sources of the remaining differences (not isolated here): the reduced CI space, single-configuration orbitals
instead of extended-average-level optimisation, and no Breit/QED.

## Status

Complete: all 48 energies, 42 IPs, analysis, validation. Not done (by design): the full single/double-excitation active space.
Known source defect and build workaround: see `../README.md`.
