# Fe XVII–XXV ionization potentials: minimal ground-configuration benchmark

Compares FAC ionization potentials (IPs), computed from **one ground
configuration per ion**, with Table 4 of Deprince et al., A&A 626, A83 (2019)
(μ_ne = 0.10 / 0.25 a.u., with and without electron–electron screening).

**This is a minimal ground-configuration benchmark, not a reproduction of the
paper.** The paper uses an MCDF restricted active space (single and double
excitations into n = 3 from several reference configurations), extended-average-level
optimisation, the Breit interaction and QED corrections. None of that is done here.

## Files

| file | purpose |
|---|---|
| `fe_xvii_ground.py` … `fe_xxvi_ground.py` | one ion each (a 3-line wrapper); same command line for all |
| `fe_ground_common.py` | the shared FAC sequence, runtime checks and manifest |
| `supplement_fe_xxvi_dirac_yukawa.py` | independent Dirac–Yukawa 1s solver (**not FAC**, see "Fe XXVI") |
| `run_fe_ground_configs.sh` | runs 4 cases × 10 ions sequentially, one fresh Python process each |
| `collect_fe_ground_ips.py` | energies → IPs → comparison with Table 4 (CSV) |
| `results/run_20261013_minimal/` | this run: `runs/`, `logs/`, `run_summary.csv`, `ground_energies.csv`, `ip_differences.csv`, `table4_ip_comparison.csv` |

```bash
python fe_xviii_ground.py --mu-ne 0.1 --mu-ee 0.0 --outdir results/x/fe_xviii/case_A
./run_fe_ground_configs.sh [RUN_ID]          # refuses to reuse an existing RUN_ID
python collect_fe_ground_ips.py results/RUN_ID
```

Results are in `results/<RUN_ID>/` instead of directly in `results/`, so that a
new run never overwrites an old one. `final.pot` and `levels.b` are not committed.

## Method

* Base: `integration/debye-screening-complete` (`ff08ede4`), branch
  `research/fe-xvii-xxv-ground-ip`. Fixes present in the C code:
  `00eb229b` stable `DensityToSZ` (`_ee_E`, `_ee_screened_sz`),
  `0650edb1` `radial:hxs_iter0`, `8c8b1bfb` Debye `RadialBound` bracket.
  The C code is unchanged.
* Sequence (same as the validated Fe II–IV cases): free-ion `OptimizeRadial`
  (fixes the grid, no ion-sphere cut-off), then `PlasmaScreen` in steps of
  μ = 0.05 up to the target with `OptimizeRadial` after each step,
  `radial:hxs_iter0 = 8`, `Structure` on the single ground configuration.
* Plasma: n = 10⁻⁶ Å⁻³, `ups = 0`, T = 4πn/μ² (in Hartree) so that λ_D = √(T/4πn) = 1/μ.
* **Independent control of the two screenings.** FAC derives both from one
  Debye length. `orbital:debye_mode = 1` screens the electron–nucleus term
  (−Z e^{−μr}/r) and `radial:ee_screen` (1/0) switches the Yukawa kernel for
  the electron–electron term (Slater integrals and SCF direct potential).
  `ee_screen = 0` leaves the electron–electron repulsion in place, as plain
  Coulomb. Only μ_ee = 0 or μ_ee = μ_ne is possible, which are exactly the Table 4 cases.
* Runtime checks in every run (written to `manifest.json`): plasma model
  `mps = 1`; `dps` from `final.pot` gives μ_ne,used = 1/dps (equal to the request
  within 10⁻⁶ relative); no ion-sphere cut (`ips = 0`); `NELE` equals the configuration's
  electron count; and the asymptotic charge r(V_c+U)+ZPS at the last grid point is
  0 when both are screened, N−1 when only the nucleus is, −(zion+1) for the free ion.
  There is no getter for `ee_screen`, so this asymptotic charge is the runtime evidence for it.
* Ground energy: `E0` + the lowest relative level energy (always level 0), in eV,
  negative. IP = E0(next ion) − E0(ion), converted with 1 Ha = 27.211386 eV.

## Results (run `run_20261013_minimal`)

* 36 of 40 ground-state energies are valid; 32 of 36 IPs computed. All 4 Fe XXVI runs fail
  the screening check (below), so the four Fe XXV IPs are **not** computed.
* |IP − Table 4| over the 32 IPs: median 1.39 eV, mean 1.92 eV, maximum 6.87 eV
  (Fe XXIII, 0.38 %). Fe XVII: ≤ 0.31 eV; Fe XX, XXI: ≤ 0.5 eV; Fe XXIV: ≤ 1.2 eV;
  Fe XVIII, XIX, XXII: 1.6–2.7 eV; Fe XXIII: 6.7–6.9 eV.
* The discrepancy is nearly the same in all four cases for each ion (e.g. Fe XXIII
  −6.87, −6.83, −6.73, −6.68 eV), so it is an ion-dependent offset, not a screening error.
* Screening trends: in all 32 cases the sign of the electron–electron effect (B−A, D−C)
  and the decrease of the IP with μ_ne agree with Table 4. The size of the e–e effect
  agrees with the paper to 0.01–0.2 eV (e.g. Fe XVII 24.09 vs 24.03 eV at μ = 0.1,
  58.73 vs 58.58 eV at 0.25).

Ground levels found (lowest level of the single configuration): Fe XVII ¹S₀,
XVIII ²P₃/₂, XIX ³P₂, XX ⁴S₃/₂, XXI ³P₀, XXII ²P₁/₂, XXIII ¹S₀, XXIV ²S₁/₂, XXV ¹S₀
(`two_j` and parity are in `ground_energies.csv`).

### Where the discrepancies probably come from (not tested here)

* Fe XXII–XXIV: the sign pattern (IP too large for Fe XXII, too small for Fe XXIII,
  smaller for XXIV) is what a missing 2s²–2p² mixing would give, because it lowers
  E(Fe XXIII) by several eV. This is consistent with the data but was not tested.
  Because the offset is the same in all four screening cases, it is a
  configuration-space effect and not a screening one.
* Remaining differences could come from the restricted configuration, the orbital
  optimisation (here: the single ground configuration, not EAL), Breit/QED, or
  numerical error. The data here cannot separate these; a run with the active space would.

## Two problems found

### 1. Fe XXVI cannot be screened by FAC (blocks the Fe XXV IP)

For a one-electron ion `PotentialHX` (`faclib/radial.c`, the branch
`N+nqf < 1+EPS5`) sets `ZPS = 0` and returns, so the Debye nuclear term is never
applied: `dps` is set to 1/μ, but the asymptotic charge stays at the free-ion
value −26. The check in the script catches it (`status = FAILED`), the logs are in
`results/run_20261013_minimal/logs/fe_xxvi_case_*.log`. The earlier
`fe_calc/fe_xxv.py` outputs (Fe XXVI E0 = −9277.574 eV at μ = 0.1, only 0.12 eV above
the free ion) are affected by the same problem.

The minimal missing capability is to keep the nuclear Debye term for N = 1.
No energy was substituted. For orientation only, `supplement_fe_xxvi_dirac_yukawa.py`
(independent, point nucleus, reproduces the analytic Dirac 1s energy at μ = 0) gives
1s energies of −9211.20 eV (μ = 0.1) and −9106.12 eV (μ = 0.25). Combined with the
FAC Fe XXV energies this would give Fe XXV IPs of about 8752 / 8755 / 8647 / 8654 eV
for cases A / B / C / D, which is 10–14 eV below Table 4 (the smaller number if the FAC Fe XXVI zero-point offset of about 4 eV is applied) (8766.28 / 8768.90 /
8661.14 / 8667.86) in all four cases, with the B−A and A−C differences matching the
paper to about 0.1 eV. This mixes two codes, is not part of the 36 IPs, and the
offset is unexplained.

### 2. `radial:hxs_iter0 = 8` crashes the optimised (`-O2`) macOS build

Platform: macOS arm64, Apple clang. With the default `-O2` build, the first
screened `OptimizeRadial` with `hxs_iter0 = 8` dies with SIGBUS in
`PotentialHX` (`faclib/radial.c:2034`, the loop `ut[m-1] < ut[m]`; located
with an AddressSanitizer build). The same input runs without `hxs_iter0`, and
runs at `-O0` and in one `-O2 -g` build, so this looks like undefined behaviour
that depends on code generation; the cause was not found. The C code was not changed.
The extension used here was built with `-O0 -g` (compiler flags only; edit
`CFLAGS = -O2 -fPIC` → `-O0 -g -fPIC` in `faclib/Makefile` and `python/Makefile`
after `./configure`). As a check, this build reproduces the Linux-validated Fe III
10-configuration energies in
`research/debye-screening/results/reference-configurations/reference_energies.csv`
exactly: −34542.4125 eV (μ = 0) and −33611.0422 eV (μ = 0.1).

Build (nothing installed): in the worktree, with `CC` set to Apple clang because the
conda cross-compiler wrapper failed in `configure`:

```bash
CC=/Library/Developer/CommandLineTools/usr/bin/clang CFLAGS='-O2 -fPIC' ./configure \
    --build=aarch64-apple-darwin20.0.0 --host=aarch64-apple-darwin20.0.0
# set CFLAGS = -O0 -g -fPIC in faclib/Makefile and python/Makefile
make pfac
```

The scripts refuse to run unless the imported `pfac.fac` lies inside this
repository (the conda environment has its own, different, global `pfac`).
