# Deprince et al. (2019): Fe XVII - Fe XXV ionization potentials, reduced configuration models

* Paper: J. Deprince et al., *Plasma environment effects on K lines of astrophysical interest. II. Ionization potentials,
  K thresholds, radiative rates, and Auger widths in Ne- through He-like iron ions (Fe XVII-Fe XXV)*, A&A 626, A83 (2019).
  Compared: **Table 2** (fully screened, plus NIST for the isolated ion) and **Table 4** (independent screening, convention
  (mu_ne, mu_ee)). The reference configurations are in Sec. 2.1 of the paper (its Table 1 lists the mu values).
* Fe XXVI (hydrogen-like) is not in the paper; it is the daughter ion needed for IP(Fe XXV).

## Models (always kept apart)

`ground_only`: the ground configuration alone (the existing validated benchmark, primary). `reference_only`: the listed reference
configurations (`configurations.py`: Fe XVII 5, XVIII 2, XIX 5, XX-XXIII 6, XXIV 8 = all 3-electron occupations of 1s, 2s, 2p, XXV 3,
XXVI 1), with orbitals optimised on the ground configuration and the references added in `Structure`. **No single or double excitations
are generated**; neither model is the published MCDF calculation (no Breit, QED or extended-average-level optimisation).

## Screening conditions (nuclear and electron-electron controlled separately)

| case | mu_ne | mu_ee |
|---|---|---|
| unscreened | 0 | 0 (plain FAC; `PlasmaScreen` never called) |
| A | 0.1 | 0 |
| B | 0.1 | 0.1 |
| C | 0.25 | 0 |
| D | 0.25 | 0.25 |

mu_ne is applied with `orbital:debye_mode 1`; `radial:ee_screen` is 1 or 0. `mu_ee = 0` keeps the electron-electron repulsion as plain Coulomb.
FAC uses one Debye length, so only mu_ee = 0 or mu_ee = mu_ne is possible (the Table 4 cases). Other settings as in the Paper A study.
10 ions x 5 cases x 2 models = 100 runs (the baseline `ground_only` is 50).

## Commands

```bash
python calculate_energies.py --check-only | --preflight | --dry-run
python calculate_energies.py --ions fe_xvii --models ground_only --cases unscreened A B      # a selection
python calculate_energies.py --models ground_only reference_only [--jobs 2] [--overwrite]   # default: ground_only, all ions/cases
python analyze_ionization_potentials.py [--models ground_only]      # -> ionization_potentials.csv, comparison_summary.md
python tests/run_validation.py                                      # -> results/validation_report.md
```

Conventions (energy units, IP definition, ground-level selection and its limits) are as in `../deprince_2020_fe_ii_viii/README.md`.

## Validation (`results/validation_report.md`: 29 passed, 0 failed, 0 skipped)

Reference configurations and electron counts (Fe XXIV = 8 occupations); Fe XVII unscreened, A, B (raw logs/outputs inspected; unscreened
path `mps = -1`; nuclear-only asymptotic charge 9 = N-1, fully screened 0, so the e-e **direct potential** is screened; E(unscreened) < E(B) < E(A));
Fe XXVI at mu = 0 and with nuclear screening: energy shifts +70.5485 / +175.6218 eV agree with an independent point-nucleus Dirac-Yukawa solver
to 5 meV, and A = B, C = D exactly (no mu_ee dependence); Fe XXV IP from the matching Fe XXVI energy; 40 earlier validated `ground_only`
A-D energies reproduced (tolerance 1e-3 eV); physical ordering of all energies and positivity/screening-dependence of all IPs; analysis on
synthetic data; resume/overwrite/timeout; preflight.

## Results (100/100 runs converged; 90/90 IPs; see `comparison_summary.md`)

| | `ground_only` (n, mean / max abs, eV) | `reference_only` |
|---|---|---|
| IP vs published, all cases | 45: 2.84 / 10.17 (Fe XXV, A) | 45: 2.24 / 10.17 |
| IP vs Table 2 (mu = 0) | 9: 2.85 / 10.09 (Fe XXV) | 9: 2.27 / 10.09 |
| IP vs Table 4 (A-D) | 36: 2.83 / 10.17 | 36: 2.23 / 10.17 |
| IP vs NIST (mu = 0) | 9: 2.60 / 6.80 (Fe XXIII) | 9: 1.85 / 4.07 (Fe XX) |
| screening shift vs unscreened, calc - published | 36: 0.07 / 0.24 | 36: 0.07 / 0.32 |
| e-e effect (B-A, D-C), calc - published | 18: 0.06 / 0.15 | 18: 0.07 / 0.22 |

The screening effects, including the electron-electron effect, agree with the paper to 0.35 eV or better in every case; the remaining IP differences are
ion-dependent offsets that are nearly the same in all screening cases (so they come from the unscreened energies, not the screening). The Fe XXV
offset (about -10 eV; paper's unscreened value 8836.74 eV vs 8826.65 eV here and 8828.19 eV NIST) is the largest and is unexplained. Adding the reference
configurations lowers the energies of Fe XVII and Fe XX - XXV by 0.006 - 7.5 eV (the largest, 7.5 eV, is 2s2-2p2 mixing in Fe XXIII) and leaves Fe XVIII, XIX and XXVI unchanged (no mixing
with the ground level for J/parity reasons). Differences may come from the reduced configuration space, orbital optimisation, relativistic corrections
(Breit/QED) or the paper's own energies; these are not separated here.

## Status

Complete: all 100 energies (50 per model), 90 IPs, analysis, validation. Not done (by design): the single/double-excitation active space.
Known source defect and build workaround: see `../README.md`.
