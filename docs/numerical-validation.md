# Numerical validation

This page lists the checks the screened FAC has passed and the regression
numbers to compare against. Energies are in eV, which is what FAC prints.
Ionization potentials (IPs) are differences of E0.

## 1. Unscreened regression (μ = 0)

- **Test:** `research/debye-screening/tests/debye_scf/run_regression.py
  <baseline_build> <test_build> <out>`.
- **Cases:** He-like Fe, Fe IX, Fe III with configuration energies, and Fe IX
  with `ee_screen 1` but no plasma.
- **Passes if:** E0, all levels and the full potential arrays are
  **bitwise identical** to unmodified FAC (c4a016a).
- **Status:** passed for the v2 patch alone, with `hxs_iter0` and with the
  `RadialBound` fix.
- **Extra check:** the ten-configuration Fe III run gives E0 = −34542.4125 eV
  and 1640 levels, identical to unmodified FAC.

## 2. Screened direct-potential kernel

`archive/experiments/v2-development/kernel_unit/` compares the potential
with closed-form Yukawa potentials of model densities.

- **Range tested:** μ = 0–5 a₀⁻¹, r up to 1500 a₀.
- **v2 accuracy:** 10⁻¹⁰–10⁻⁹, never non-finite.
- **Exception:** a diffuse density at μ = 5 is only accurate to 3·10⁻⁶. There
  the screening length is shorter than the grid spacing, which limits any
  method.

## 3. Independent implementation

The `ee-screening-scf` branch has a separately written screened potential
(an `sfac` study, archived in `archive/experiments/deprince-t2-sfac/`). All
converged reference-space energies agree to ≤ 0.2 meV.

## 4. First-order Debye-Hückel limit (archived `sfac` study)

For small μ, ΔE ≈ μ[ZN − N(N−1)/2].

- He-, Li- and Ne-like Fe agree to 0.02–0.06 %.
- Fe II–IX agree to 0.2–1.2 %, which is first-order orbital relaxation.

## 5. Solver fixes

- **`hxs_iter0`.**
  - Default 0 is bitwise unchanged.
  - With 8, the Fe II ladder to μ = 0.1 and the Fe III ladder to μ = 0.25
    converge, and give the same orbital energies (to 1.5·10⁻⁶ Ha) wherever
    the default also converges.
  - Logs: `research/debye-screening/results/hxs-iter0-validation/`.
- **`RadialBound` Debye bracket.**
  - The 17 reference cases that already converged are bitwise unchanged.
  - `gref` Fe II μ = 0.25 now converges: E0 = −32293.7705 eV, identical to
    the `sfac` study.
  - The Sturm diagnostic showed 3d₃/₂ bound at −0.288 Ha in the failing
    potential (`archive/experiments/warm-start-continuation/E6_*.log`).

## 6. Regression references: reference-configuration IPs (eV)

The source is `research/debye-screening/results/reference-configurations/`.
All 18 E0 are reproduced exactly by the per-ion scripts.

| IP | model | μ = 0 | μ = 0.1 | μ = 0.25 |
|---|---|---:|---:|---:|
| Fe II → Fe III | ref | 18.8717 | 14.1881 | 8.7175 |
| Fe II → Fe III | gref | 15.5203 | 10.2258 | 4.6872 |
| Fe III → Fe IV | ref | 19.6518 | 11.1572 | 0.0253 |
| Fe III → Fe IV | gref | 28.0942 | 20.1338 | 9.2698 |

For comparison, Deprince et al. Table 2 gives:

- Fe II: 18.89 / 14.20 / 7.78 eV;
- Fe III: 32.93 / 25.60 / 16.42 eV;
- NIST (isolated ion): Fe II 16.19920 eV, Fe III 30.651 eV.

The `ref` Fe III values are an orbital-averaging artefact (see
`PROJECT_STATUS.md`). They are regression numbers, not physics.

## 7. Runtime checks in every production script

- correct electron count (`NELE`);
- `mps = 1` and `dps · μ = 1` to 10⁻⁶;
- asymptotic r·V = 0 when screened, −(zion + 1) when free;
- `ips = 0`: no ion-sphere cut;
- no "Maximum iteration reached in OptimizeRadial";
- finite E0.

The full active-space driver also checks that FAC's level count per J block
equals the counted CSFs, and that the final `HXS` = 1.
