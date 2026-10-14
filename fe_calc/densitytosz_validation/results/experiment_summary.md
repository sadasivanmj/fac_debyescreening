# DensityToSZ A/B experiment: effect of the screened SCF direct potential on Fe II / Fe III / Fe IV

*Reduced reference-configuration model (gref), not the paper's MCDF active space.* Generated from the saved results by `write_reports.py`.

## Setup

* Upstream `jmaxk` fix `fb32c16a30` (Ming Feng Gu, FAC PR #851, parent `5f8c8d23e8`) cherry-picked as `603c0e5358` (changed lines identical: True). Both builds use `-O2 -fPIC`, **no** `-ftrivial-auto-var-init=zero` (the workaround is gone). Negative control (identical flags, pre-fix source): 18 of 20 runs of a reproducer crash; post-fix 0 of 20 (`jmaxk_negative_control.json`). All 148 framework energies of the earlier sweeps are reproduced exactly by the fixed build (max |dE| = 0 eV).
* Variant B = production screened `DensityToSZ` (FAC source `603c0e5358`); variant A = B plus one line, `double mu = 0.0;` in `DensityToSZ` (`54f043349c`). FAC-source diff A vs B: ['faclib/orbital.c'] only. `radial:ee_screen` is on in both, so the screened Slater integrals and nuclear screening are identical.
* Fe II, Fe III, Fe IV, existing reference-configuration lists, `gref` orbitals, mu = 0, 0.1, 0.25, same sequence and settings; 18 runs.

## Outcome of the 18 energy runs

* **B: 9/9 converged and validated; A: 3/9 converged (all mu = 0 controls, identical to B to 0.0 eV); the six screened A runs failed** with `No Debye bracket in RadialBound (unbound): n=3 kappa=2` (Fe II and Fe III at the first step mu = 0.05; Fe IV at mu = 0.1). Logs preserved as `run.FAILED.log`.
* Therefore dE = E(B) - E(A), and dI, **cannot be measured at the published mu = 0.1 and 0.25**; it is measured below, at smaller mu, where A converges.

### Why A fails (evidence, not proof)

With a Coulomb direct potential the total electron potential energy `r V_tot` tends to **N-1 = 22-24** at large r (B: 0), because the nuclear term is screened but the electrons' is not. It turns repulsive at `r_c = ln(Z/(N-1))/mu`; measured 8.14 a0 (Fe II, mu = 0.01) and 6.15 a0 (Fe III, mu = 0.02), analytic 8.0 and 6.1. FAC reports no bound 3d state once that attractive region shrinks to ~4 a0 (Fe II fails from mu = 0.02, Fe III from 0.03; Fe IV survives to 0.05). So A is not a usable self-consistent model for these ions; it is a controlled artefact, not a physical alternative. (`variant_A_diagnosis/`.) Orbital-by-orbital energies were not extracted.

## Results (B at all mu; A only at mu = 0)

| Transition | mu | IP A | IP B | paper | B - paper | B rel % | NIST | B - NIST | shift B | shift paper | shift err B |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Fe II -> Fe III | 0.0 | 15.5203 | 15.5203 | 18.89 | -3.3697 | -17.83 | 16.1992 | -0.678 | 0.0 | 0.0 | 0.0 |
| Fe II -> Fe III | 0.1 | — | 10.2258 | 14.2 | -3.9742 | -27.98 | — | — | -5.2945 | -4.69 | -0.604 |
| Fe II -> Fe III | 0.25 | — | 4.6872 | 7.78 | -3.0928 | -39.75 | — | — | -10.833 | -11.11 | 0.2769 |
| Fe III -> Fe IV | 0.0 | 28.0942 | 28.0942 | 32.93 | -4.8358 | -14.68 | 30.651 | -2.556 | 0.0 | 0.0 | 0.0 |
| Fe III -> Fe IV | 0.1 | — | 20.1338 | 25.6 | -5.4662 | -21.35 | — | — | -7.9604 | -7.33 | -0.630 |
| Fe III -> Fe IV | 0.25 | — | 9.2698 | 16.42 | -7.1502 | -43.54 | — | — | -18.824 | -16.51 | -2.314 |

Unscreened baseline errors (identical in A and B, so unrelated to DensityToSZ): Fe II -> III 15.5203 eV (paper 18.89, NIST 16.1992); Fe III -> IV 28.0942 eV (paper 32.93, NIST 30.651). Screening-shift errors of B against the paper: Fe II -0.60 / +0.28 eV and Fe III -0.63 / -2.31 eV at mu = 0.1 / 0.25.

## Direct-potential effect where A converges

Single SCF steps from the free ion, A and B, same settings (supplementary: not the published mu). With FAC's default `radial:csi`, B-A looks large (e.g. Fe II dE = -0.173 eV and dI(Fe II) = -0.404 eV at mu = 0.01), but that is the `csi` artefact of Section 8 of `densitytosz_validation.md`. With `radial:csi = 0` in both builds the pure direct-potential effect is:

| mu | dE Fe II | dE Fe III | dE Fe IV | dI Fe II | dI Fe III | A converged II/III/IV |
|---|---|---|---|---|---|---|
| 0.001 | -0.000 | 0.000 | -0.000 | 0.000 | -0.000 | y/y/y |
| 0.003 | -0.001 | -0.000 | -0.000 | 0.001 | 0.000 | y/y/y |
| 0.005 | 0.060 | 0.002 | -0.001 | -0.058 | -0.003 | y/y/y |
| 0.01 | 0.089 | 0.077 | 0.051 | -0.011 | -0.026 | y/y/y |
| 0.02 | — | 0.076 | 0.071 | — | -0.005 | n/y/y |
| 0.03 | — | — | 0.065 | — | — | n/n/y |
| 0.05 | — | — | — | — | — | n/n/n |

The identities `dI(Fe II) = dE(Fe III) - dE(Fe II)` and `dI(Fe III) = dE(Fe IV) - dE(Fe III)` hold to < 1e-8 eV in every row. For the highly charged ions (Fe XVII - Fe XXVI, mu = 0.1 and 0.25, both variants converge) the direct-potential part of the energy is at most 0.082 eV, while the electron-electron screening effect itself is -2.7 to -120.6 eV at mu = 0.1 (larger at 0.25): it comes almost entirely from the screened Slater integrals (`results/physical_checks/`). This is the variational property: the direct potential only changes the orbitals.

## Answers

1. **Is `DensityToSZ()` mathematically correct?** Yes. It implements `r V = e^{-mu r} int D sinh(mu r')/(mu r') + r (sinh(mu r)/(mu r)) int D e^{-mu r'}/r'` (via the stable `E(x)` form) with the Coulomb limit bit-identical to the pristine library. One declared criterion is marginally missed (deep tail of an extreme density, below).
2. **Does it agree with an independent Yukawa Hartree calculation?** Yes: 1e-7..4e-7 of the peak for compact densities, 1e-8..3e-7 for valence/multishell, 1e-10..1e-8 for diffuse at mu <= 0.25, on the real grid; the floor is FAC's own quadrature and the same in the Coulomb path.
3. **Is the long-range behaviour correct?** Yes: slope `-mu` to 1e-13, deviation 1e-7..1e-8 down to 1e-110, no plateau/truncation/overflow. Deep-tail relative accuracy degrades only for the extreme `2a < mu` density (1.0047e-5 against the declared 1e-5 at values > 1e-17 of the peak; the one failed pre-declared criterion).
4. **Change of the Fe II, III, IV total energies?** At mu = 0 none (0.0 eV). At 0.1 and 0.25 **not measurable** (A did not converge). Where measurable (mu <= 0.03, `csi = 0`): -0.001 to 0.089 eV (Fe II), -0.000 to 0.077 eV (Fe III), -0.001 to 0.071 eV (Fe IV); for Fe XVII - XXVI at the published mu at most 0.08 eV.
5. **Change of both IPs?** At the published mu not measurable. Where measurable: |dI| <= 0.06 eV with `csi = 0` (<= 0.03 eV for mu >= 0.01) for both ionizations; the identity holds. With the default `csi` up to ~0.4 eV, which is the `csi` artefact.
6. **Does it explain a substantial part of the Fe III screening-shift discrepancy?** No evidence of that. The discrepancy of B is -0.63 eV (mu = 0.1) and -2.31 eV (0.25); the measurable direct-potential contribution is <= 0.03 eV, and in the highly charged analogues <= 0.08 eV in energy. The conclusion for mu = 0.25 itself rests on extrapolation, because A cannot be run there.
7. **Better, worse, or no difference vs Deprince?** At mu = 0: identical (both 3.4 and 4.8 eV below the paper). At the published mu: A cannot be compared (no converged energies); B reproduces the screening shifts to 0.3-0.6 eV for Fe II -> III and 0.6 / 2.3 eV for Fe III -> IV. Where both exist the difference is < 0.06 eV, far below any discrepancy: **little difference**. No statement that either is "more correct" follows from agreement with the literature; correctness is established by the potential-level tests.
8. **What remains unexplained?** (i) The unscreened IP errors (-3.37 / -4.84 eV vs paper; -0.68 / -2.56 vs NIST) and the Fe III -> IV shift errors: reduced CI space, single-configuration orbitals, no Breit/QED are candidates, none tested here. (ii) The small-mu first-order deficiency of B (up to -1.2%) tied to FAC's `csi`/`ihx` block: `csi = 0` removes it, the mechanism is not understood and has no effect at mu >= 0.1 for these ions. (iii) The exact reason A loses its 3d bound state (consistent with the repulsive tail, not proven orbital by orbital). (iv) The deep-tail accuracy for `2a < mu` densities. (v) FAC's grid-independent quadrature floor of 1e-8..1e-7.

## Limitations

* A cannot reach the published mu, so the A/B comparison at 0.1 and 0.25 is replaced by measurements at smaller mu plus the highly charged analogues.
* The ladders are single SCF steps from the free ion; SCF-convergence jitter in E(mu) is ~0.03 eV, comparable to the small dI values.
* `jmaxk` upstream fix is in both experimental builds only; the production branch has not been changed. Experiment branches are local; nothing was pushed or merged.
