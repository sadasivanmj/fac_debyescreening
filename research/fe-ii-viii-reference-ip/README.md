# Fe II – Fe VIII ionization potentials from a reduced configuration space

FAC ionization potentials (IPs) of Fe II → Fe III … Fe VIII → Fe IX at μ = 0, 0.1 and
0.25 a.u., using **only the reference configurations** of Deprince et al.
(no single or double excitations), compared with Table 2 of
A&A 643, A57 (2020) and with NIST.

**This is a reduced-model calculation, not a reproduction of the paper.** The paper
uses an MCDF active space of 12 268 – 85 798 configuration state functions (CSFs) per
ion plus extended-average-level optimisation, Breit and QED. Here FAC builds the CI
only among the listed reference configurations (29 to 1640 levels, see
`configuration_manifest.csv`), with FAC's own orbitals.

Branch `research/fe-ii-viii-reference-ip`, based on `research/fe-xxvi-screened-ground`
(`07062f73`). **No C code was changed.** Run: `results/run_20261013_reference/`.

## Model

* Configurations: `reference_configs.py` (Table 1 of the paper, holes written out; closed
  core `1s2 2s2 2p6 3s2 3p6`). The generated Fe II–IV lists are identical to the validated
  ones in `../debye-screening/cases/fe/`. Every electron count is checked before FAC is called.
* **Fe IX** (needed for IP(Fe VIII)) is not in Paper IV. Its reference configurations are
  in Paper III (Deprince et al., A&A 635, A70, Table 1): `3p6, [3p]3d, [2p]3d, [1s]3d`,
  i.e. the same pattern as Fe VIII–IV, used with the same reduced method. The Fe VIII IP is
  therefore computed, not pending.
* Two orbital models (same CI space; only the orbitals differ):
  `gref` = `OptimizeRadial` on the ground configuration only (**primary**), `ref` = on the
  equal-weight average of all reference configurations. `ref` is kept for completeness; its
  Fe III → Fe IV value is an averaging artefact (see the earlier
  `REFERENCE_CONFIGURATION_VALIDATION.md`).
* Screening: electron–nucleus (`orbital:debye_mode 1`) and electron–electron
  (`radial:ee_screen 1`) with the same μ; `PlasmaScreen(zion, 1e-6 Å⁻³, T, 1, 0, 1)`
  with T = 4πn/μ²; free-ion SCF first, then μ in steps of 0.05, `radial:hxs_iter0 = 8`.
  Same sequence as the validated Fe II–IV cases. μ = 0: no plasma options.
* IP = [E0(next ion) − E0(ion)] × 27.211386 eV/Ha / 27.211386, i.e. the difference of
  the E0 values. FAC writes E0 (energy of the lowest level) in eV, total energy, negative.
* Runtime checks per run (all passed): `mps = 1` and `dps·μ = 1`, no ion-sphere cut,
  asymptotic r·V = 0 (electron–electron screening active; −(zion+1) for the free ion),
  `NELE` correct, finite E0, level table referenced to its lowest level, and no
  "Maximum iteration reached in OptimizeRadial" in the log (checked by the collector).

## Results

All 48 runs (8 ions × 3 μ × 2 models) converged and passed every check; 42 of 42 IPs
computed (165 s in total). Regression: the 18 Fe II–IV energies are identical (max
|ΔE| = 0 eV) to the Linux-validated `reference_energies.csv`.

IP in eV, `gref` model (paper Table 2 in brackets); full data in `ionization_potentials.csv`:

| IP | μ = 0 | μ = 0.1 | μ = 0.25 | NIST | err vs NIST (μ=0) |
|---|---|---|---|---|---|
| Fe II | 15.52 (18.89) | 10.23 (14.20) | 4.69 (7.78) | 16.199 | −0.68 |
| Fe III | 28.09 (32.93) | 20.13 (25.60) | 9.27 (16.42) | 30.651 | −2.56 |
| Fe IV | 53.33 (55.37) | 42.92 (44.65) | 28.17 (29.64) | 54.91 | −1.58 |
| Fe V | 73.96 (74.92) | 60.84 (61.56) | 42.48 (42.83) | 75.0 | −1.04 |
| Fe VI | 97.99 (98.20) | 82.23 (82.20) | 60.20 (59.79) | 98.985 | −0.99 |
| Fe VII | 123.98 (123.60) | 105.57 (104.96) | 79.84 (78.86) | 124.98 | −1.00 |
| Fe VIII | 149.94 (151.14) | 128.83 (129.85) | 99.43 (99.98) | 151.060 | −1.13 |

Screening shifts ΔIP(μ) = IP(μ) − IP(0), calc / paper (eV):

| IP | μ = 0.1 | μ = 0.25 |
|---|---|---|
| Fe II | −5.29 / −4.69 | −10.83 / −11.11 |
| Fe III | −7.96 / −7.33 | −18.82 / −16.51 |
| Fe IV | −10.41 / −10.72 | −25.17 / −25.73 |
| Fe V | −13.11 / −13.36 | −31.47 / −32.09 |
| Fe VI | −15.76 / −16.00 | −37.80 / −38.41 |
| Fe VII | −18.42 / −18.64 | −44.14 / −44.74 |
| Fe VIII | −21.10 / −21.29 | −50.50 / −51.16 |

* **Against NIST (μ = 0):** `gref` is 0.7–2.6 eV low (mean 1.28 eV); FAC is below NIST for
  every ion, by about 1 eV for Fe V–VIII. The paper's own values are 0.1–2.7 eV from NIST
  (e.g. Fe II +2.7, Fe III +2.3, Fe VII −1.4, Fe VIII +0.1), so the paper is not the same as NIST either.
* **Against the paper:** |IP − paper| for `gref`: mean 1.93 eV, median 1.02 eV, maximum 7.15 eV
  (Fe III → Fe IV at μ = 0.25). It is within 1.3 eV for Fe V–IX at μ = 0 and 0.1 and
  within about 1 eV for Fe VI–IX at 0.25. `ref`: mean 2.85, median 0.93, max 16.4 eV
  (the Fe III artefact); for Fe IV–IX it is close to `gref` (0.2–0.5 eV lower).
* **Screening trends:** reproduced. In both models the IP decreases monotonically with μ
  for all seven ionizations, the shift grows nearly linearly with the effective charge
  (≈ −2.65 eV per unit Z_eff at μ = 0.1, close to the Debye–Hückel limit −27.21 μ Z_eff
  eV), and the shifts agree with the paper to 0.2–0.7 eV for Fe IV–VIII. They are least
  accurate for Fe II–III (up to 2.3 eV, Fe III at μ = 0.25).
* **Likely sources of discrepancy (not tested here):** the largest differences are for the
  lowest charges (Fe II–IV), where 3d–4s correlation and the 3d⁶ term structure matter most
  and the reduced CI space (1640 levels for Fe III against 31 501 CSFs in the paper) cannot
  provide it; orbital optimisation on one configuration instead of EAL; missing Breit and
  QED terms. The ~1 eV deficit that is nearly constant for Fe V–IX at μ = 0 and the
  smaller shifts at μ = 0.25 may also reflect a missing or different treatment of
  correlation, but this run cannot separate the causes. The full active space is the next test.

## Platform problem found and how it was handled (important)

FAC crashed (SIGBUS) on macOS arm64 for some screened runs, for example Fe III `gref` at
μ = 0.1 with the plain `-O0` build, and earlier for Fe XXVI with `hxs_iter0` at `-O2`.
Cause (from an `-O0 -g` build under lldb): in `PotentialHX1` (`faclib/radial.c`,
function returning `Max(jmax, jmaxk)`) the local `jmaxk` is only assigned inside a block that is
skipped when `mps == 1` (Debye mode), so the return value contains stack garbage (here `jmax =
1431655765 = 0x55555555`), and `PotentialHX` then reads `ut[m]` far outside the array
(`radial.c:2036`). On Linux the garbage happened to be harmless.

Following the instruction not to change unrelated source code, **the C code was not modified.**
The extension was built with `make pfac` and the compiler flag
`-ftrivial-auto-var-init=zero` (plus `-O0 -g`), set in `faclib/Makefile` and
`python/Makefile` (`CFLAGS = -O0 -g -fPIC -ftrivial-auto-var-init=zero`) after
`./configure` with Apple clang. This starts `jmaxk` at 0, so `Max(jmax, jmaxk) = jmax` as
intended. Checks: the five smoke tests and all 18 Fe II–IV energies reproduce the Linux
numbers exactly, and the 40 Fe XVII–XXVI energies of the previous task are reproduced
with max |ΔE| = 0 eV (a rerun kept untracked in
`../fe-xvii-xxv-ground-ip/results/run_20261013_zeroinit_check/`), so earlier results are
not affected. A source fix would be to initialise `jmaxk = jmax` (or 0) at its
declaration in `PotentialHX1`; it is **not applied** and should be decided separately.

## Files and rerun

| file | purpose |
|---|---|
| `reference_configs.py` | the configurations (run `python reference_configs.py` to print and check them) |
| `fe_reference_run.py` | one ion, model and μ: `python fe_reference_run.py fe_iii gref 0.1 OUTDIR` |
| `run_reference_ips.sh` | all 48 runs, one process each: `./run_reference_ips.sh [RUN_ID]` |
| `collect_reference_ips.py` | `python collect_reference_ips.py results/RUN_ID` → `energies.csv`, `ionization_potentials.csv` |
| `results/run_20261013_reference/` | `configuration_manifest.csv`, `energies.csv`, `ionization_potentials.csv`, `run_summary.csv`, `logs/`, `runs/*/*/{manifest.json,levels.txt}` |

Build (nothing installed): `./configure` with `CC=/Library/Developer/CommandLineTools/usr/bin/clang`
and `--build=aarch64-apple-darwin20.0.0 --host=aarch64-apple-darwin20.0.0`, edit the two
`CFLAGS` lines as above, then `make pfac`. The scripts refuse any `pfac` that is not inside
this repository. `final.pot` and `levels.b` are not committed (regenerable).
