# Project status — Debye-screened FAC for Fe II–Fe IV

Last updated 2026-10-13. Detailed numbers are in
[`docs/numerical-validation.md`](docs/numerical-validation.md); this page is
the short version.

## Validated

These are all in `faclib/`; the patch files are in
`research/debye-screening/patches/`.

- **Screened SCF direct potential (v2)** — `DensityToSZ` in `faclib/orbital.c`.
  - With μ = 0, results are bitwise identical to unmodified FAC (c4a016a).
  - Against analytic Yukawa potentials it agrees to 10⁻¹⁰–10⁻⁹ for
    μ ≤ 5 a₀⁻¹.
  - It agrees to ≤ 0.2 meV with an independently written screened potential
    (`ee-screening-scf` branch).
- **Screened Slater integrals** — `GetYkDebye` / `EEScreenMu` in
  `faclib/radial.c`, from `fac-debye-ee-screening.patch`, which is already
  applied.
- **`radial:hxs_iter0`** — an option that starts FAC's exchange ramp later on
  warm SCF restarts.
  - The default (0) is bitwise unchanged.
  - The value 8 removes the Fe II / Fe III `RadialBound` failures, which came
    from the exchange-free first iteration.
- **Debye upper bracket in `RadialBound`** — a full-grid node count at
  E = −10⁻⁶ Ha.
  - It fixes the bracket failure for nodeless orbitals near threshold (3d in
    Fe II, μ = 0.25).
  - Non-Debye runs are untouched.

## Experimental

- **Full single/double active-space driver**
  (`research/full-active-space/scripts/verify_ip_deprince_full_as.py`).
  It works, but the production run is unfinished (see below).
- **OpenMP build** (`./configure --with-mpi=omp`).
  - It gives the same energies as the serial build at μ = 0.
  - It **segfaults in `PotentialHX`** at the first screened `OptimizeRadial`.
  - Not used for production.
- **`RadialBound` Debye bracket for orbitals with nodes** (nr > 0, e.g. Fe IX
  3s at μ = 1.5–2.25). Here the bracket is accepted, but a later search stage
  runs out of iterations. This is not fixed.

## Completed calculations

- **Fe II, Fe III, Fe IV reference configurations** (Deprince et al. 2020,
  Table 1).
  - μ = 0, 0.1, 0.25 and models `ref` / `gref`: 18 energies and 12 IPs, all
    converged.
  - The runs reproduce the earlier regression results exactly.
  - Location: `research/debye-screening/results/reference-configurations/`.
- **Chen et al. (2018) He-like ions and C4+ tests.** These results date from
  2026-09-04 and were not rerun in this cleanup.
  - Location: `research/debye-screening/results/chen2018/` and
    `results/c4/`.

## Incomplete calculations

**Full single/double active space**, all even-parity J blocks
(`research/full-active-space/`).

- **Converged with the serial build, 5 of 18 jobs:**
  - Fe II μ = 0 (`ref`, `gref`)
  - Fe II μ = 0.1 (`ref`)
  - Fe III μ = 0 (`ref`, `gref`)
- **Interrupted:** Fe II μ = 0.1 `gref`, stopped by hand.
- **Not started:** everything else, including every Fe IV job with the serial
  build.
- **Cross-checks with the OpenMP build:** Fe IV μ = 0 for `ref` and `gref`.
  These are not production results.
- **No IP from the full active space has been computed yet.**

## Known limitations

- **`ref` orbital-averaging artefact.** The `ref` orbitals come from an
  equal-weight average dominated by core-hole configurations.
  - In the reference space: Fe III → Fe IV IP is 19.65 eV at μ = 0 (paper
    32.93) and 0.025 eV at μ = 0.25.
  - In the full space, the `ref` ground state of Fe II becomes 3d⁵4s²
    (J = 5/2) instead of 3d⁶4s.
- **Published active-space sizes can't be reproduced.** The paper's stated
  recipe gives 15147 / 40189 / 75309 CSFs, against the published 12268 /
  31501 / 57189.
  - None of the restrictions I tested (frozen core shells, hole limits,
    J windows) reproduces them.
- **Orbitals are not the paper's.** FAC uses one local central potential plus
  CI, while the paper used MCDHF-EAL. So agreement with Table 2 is not
  expected to be exact.
- **Ion-sphere cut.** If `PlasmaScreen` is called before the first
  `OptimizeRadial`, FAC builds the grid with an ion-sphere cut-off.
  - The scripts therefore always do a free-ion SCF first.
- **Exchange is still local and unscreened.**
- **`EEScreenMu` returns 0 when T = 0**, i.e. silently unscreened.
- **PFAC crash with an empty list.** `fac.Structure(p, [])` crashes, because
  of an uninitialised pointer in `IntFromList`.
  - Use `fac.Structure(p, [-1])` to select all J.

## Next steps

1. Finish the full active space with the serial build:
   - The Fe IV jobs need up to 9 GB each (measured with OpenMP), so run them
     one at a time.
   - Serial runtime is estimated at about 1–1.5 h per Fe IV job; it has not
     been measured.
   - Restart instructions are in `research/full-active-space/README.md`.
2. Compute the 12 full-active-space IPs and compare them with the reference
   space and with Table 2.
3. Optionally:
   - investigate the OpenMP `PotentialHX` crash (it would cut Fe IV time
     about 3×);
   - extend the `RadialBound` Debye bracket to orbitals with nodes.
