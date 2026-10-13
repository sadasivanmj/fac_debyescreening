# Numerically stable Debye–Hückel SCF direct potential (`DensityToSZ`, v2)

Worktree `/tmp/fac-debye-scf-test` (Linux host `luna`, as authorised by the user), detached at `c4a016a`. The v1 report (`REPORT_debye_scf.md`) and v1 artifacts are preserved: `patch_density_v1.py`, `patch_density_v1.diff`, `faclib/orbital.c.debye_scf_v1`, and the `_patched_pfac/` build.

## 1. Initial findings
- **Routine semantics.** Source: `faclib/orbital.c:DensityToSZ`, `interpolation.c:NewtonCotesIP/uvip3s`, `orbital.c` grid setup.
  - `d` is the radial density, with ∫d dr = N (Fe IX 18.000, He-like 2.000).
  - `z` is r·V_direct; `jps` = ∫_{r0}^∞ d/r dr; `zx` is r·V_xc out, exchange density in.
  - The grid is ρ = a·r^q + b·ln r with exactly unit index steps (`del = 1`), and `dr_drho = dr/dρ`.
  - `pot->rad[0]` = r_min > 0 (FAC's `rmin`, here 10⁻⁸-scale), never 0. The origin segment is approximated by `d0·r0/3`.
  - `NewtonCotesIP(r, x, i0, i1, m, id)` is a cumulative cubic-spline quadrature over unit index spacing. `id ≥ 0` goes forward from a preset `r[i0]`, `id < 0` backward from `r[i1]`; `m = −1` means no swap.
- **`EEScreenMu()`** (`radial.c`) returns 1/`dps` [a₀⁻¹] only if `ee_screen = 1`, `mps = 1` and `dps > 0`, otherwise 0. `dps` = [T/(4πn(1+ups))]^{1/2} a₀ comes from `InitializePS`. It is not declared in any header, so `orbital.c` needs an `extern`.
- **Root causes addressed in v2:**
  - Writing A(r) = ∫d·sinhc can overflow when μ > 2κ for a diffuse density, even though e^{−μr}A stays finite. Demonstrated: v1 is NaN/inf for the valence (μ = 5) and diffuse (μ = 1, 5) synthetic densities.
  - A NaN μ (e.g. negative T gives NaN `dps`) silently took the Coulomb path.
- **Audit facts for the requested tests:**
  - `fe_calc/fe_iii_unscreened.py` **does not exist**; the existing unscreened script is `fe_calc/fe_III_energy.py`.
  - `fe_calc/fe_iii_diag.py` sets **`ee_screen 0`**, so it never runs the screened direct potential.
  - `PYTHONPATH=…/python` imports the **global** `pfac` (site-packages), not the build. Verified, so not used.
  - The committed `fe_calc/debug_FeIII.lev` (E0 −36202.93 eV, about 1650 eV below the free ion) was produced by a build with a different (old-sign) e–n screening and is not a valid reference.

## 2. Formulation
With E(x) = (1 − e^{−2x})/(2x) = e^{−x}·sinhc(x), E(0) = 1, 0 < E ≤ 1:

- z(r) = Ã(r) + r·E(μr)·B̃(r)
- Ã(r) = ∫₀ʳ d·E(μr′)·e^{−μ(r−r′)} dr′ = e^{−μr}A(r)
- B̃(r) = ∫ᵣ^∞ (d/r′)·e^{−μ(r′−r)} dr′ = e^{μr}B(r)

This is algebraically identical to the requested form and contains no sinh. E is evaluated as −expm1(−2x)/(2x), which avoids cancellation for small x and overflow for large x.

**Cumulative integration.** `NewtonCotesIP` cannot apply an r-dependent factor inside its running sum, so the integrals are accumulated in **blocks** of exponent span μΔr ≤ 100:
- Inside a block the integrand is scaled by e^{−μ|R − r′|} relative to the block reference R, and the cumulative result is rescaled by e^{μ|R − r|} ≤ e¹⁰⁰.
- Blocks chain exactly.
- One block covers everything when μ(r_im − r₀) ≤ 100.

**Why span 100 (measured, not assumed).** `uvip3s` falls back to the trapezoid rule when a 7-point window spans less than an **absolute** 10⁻¹⁰⁰. A span of 300 pushed the scaled integrands below that and degraded accuracy to 10⁻⁶–10⁻⁵. Spans of 100 and 150 restore the 10⁻⁹ level and agree with each other; spans of 60 and 30 lose slightly through extra block starts.

**Other details.**
- `jps` keeps its meaning: B(r₀) = e^{−μr₀}B̃(r₀).
- The origin term is d₀r₀/3 × E(μr₀).
- μ = 0 runs the **original statements verbatim**.
- NaN, negative or infinite μ prints `invalid Debye screening parameter in DensityToSZ` and calls `Abort(1)`, FAC's existing convention.

## 3. Source changes
- Only `faclib/orbital.c`: **+81 / −0 lines** relative to c4a016a.
  - Helpers `_ee_E` and `_ee_screened_sz`, and `extern EEScreenMu`, inserted before `DensityToSZ` (markers `DEBYE_SCF_V2`…`_END`).
  - Inside the function: the μ declaration with its validity check, and `if (mu > 0) {screened} else {original}`.
- Exchange-correlation, `ips`/`ups`/`vxs` and extension branches are unchanged.
- `patch_density.py` (v2) builds from the original backup. It checks the location, the full-file sha256 of the original/v1 states, the function sha256 and the block match counts, refuses double application, and verifies that everything outside the helper block and function is unchanged.
- Diffs: `patch_density.diff` (original → v2) and `patch_density_v1_to_v2.diff`. Source sha256 `8fe975cb…`.

## 4. Build
- `./configure --prefix=/tmp/fac-debye-scf-test/_noinstall` (earlier), then `make pfac`: rc 0, `orbital.c` recompiled, no `orbital.c` warnings.
- Extension: `build/lib.linux-x86_64-cpython-310/pfac/fac.cpython-310-x86_64-linux-gnu.so`, frozen as `_patched_v2_pfac/` (sha256 `6a5a7583…`).
- No install target was run; `_noinstall` does not exist.

## 5. Import verification
- Python `/home/sadasivanmj/miniconda3/envs/fac_env/bin/python3`.
- `pfac` → `/tmp/fac-debye-scf-test/_patched_v2_pfac/pfac/__init__.py`; extension → `…/_patched_v2_pfac/pfac/fac.cpython-310-x86_64-linux-gnu.so`. Asserted in every run.
- The global site-packages `pfac` was never imported (mtime 2026-08-29, unchanged).

## 6. Unscreened regression (pass condition 1: PASS)
- `fe_calc/fe_III_energy.py` (unmodified; 10 Table 1 references; run from fresh directories, so the committed outputs are untouched): baseline, v1 and v2 all give E0 = **−34542.41250000 eV**, NLEV 1640, and **identical level tables** (excluding TSess).
- Minimal cases (He-like, Fe IX, Fe III; Fe IX with `ee_screen 1` and no plasma): **all 18 `final.pot` columns × 1500 points identical to the baseline** (max abs diff 0), with identical headers and `e.lev`.
- Bit-identity is meaningful here because the μ = 0 path executes the original statements, built with the same compiler and flags.

## 7. Screened diagnostic (pass condition 2: PASS)
- **Fe III, 3d⁶, λ_D = 10 a₀ (μ = 0.1 a₀⁻¹), nps = 0.01 Å⁻³, zps = 2, cold `ConfigEnergy(0)` → `OptimizeRadial` → `ConfigEnergy(1)` → `Structure`.** The copy `tests/fe_iii_diag_ee1.py` differs from `fe_calc/fe_iii_diag.py` only in `ee_screen 0 → 1`.
  - Baseline **fails** at `ConfigEnergy(0)`, SCF iteration 2, `RadialBound a`, n = 3, κ = −1 (3s), grid point 230, trial energy −8.7×10⁻¹⁵⁴.
  - v1 and v2 **complete all stages** in 48 SCF iterations: **E0 = −33625.6792 eV**, with identical 34-level tables.
  - The reference value's provenance is verified: it came from the v1 run `tests/out/r6/patch/T4_FeIII_diag_ee1`, with identical inputs.
- The unmodified `fe_iii_diag.py` (`ee_screen 0`) fails identically on all builds. That is expected: e–n-only screening is not the code path under test.
- Other screened controls (`tests/out/r10` vs `r2`/`r6`): v2 equals v1 to 10⁻⁴ eV for He-like μ = 10⁻⁴…0.25 and Fe IX μ = 10⁻⁸…1.0. Fe IX at μ = 1 is multi-block. The baseline fails at Fe IX μ = 1.0 and Fe III ee1.

**Kernel validation against the closed-form Yukawa potential of a 1s cloud** (`tests/kernel_unit/final_unit_table.txt`; FAC-like grid with q = ½, r = 10⁻⁶…1500, 2000 points; compact a = 13, valence a = 1, diffuse a = 0.06; μ = 0…5):
- v2 errors match the Coulomb quadrature level: compact 2.4–2.9×10⁻⁹, valence 0.8–5.7×10⁻⁹, diffuse 0.9–75×10⁻¹⁰.
- Exception: diffuse μ = 5 at 3.0×10⁻⁶. The kernel decay length of 0.2 a₀ is below the local grid spacing there; every method shares this limit.
- **0 non-finite values.** v1 is non-finite in 3 of the 21 cases.
- `jps` matches v1 to ≤ 2×10⁻¹¹ relative.

**Tested ranges:** μ ∈ [0, 5] a₀⁻¹, μ·r_max up to 7500, up to 75 blocks, on synthetic densities; real SCF runs for μ ≤ 1 (Fe IX) and μ ≤ 0.25 (He-like). Stability is not claimed beyond these.

**Invalid μ** (negative T, so `dps` = NaN): baseline and v1 continue silently and fail later in `RadialBound` (1s); v2 stops immediately with `invalid Debye screening parameter in DensityToSZ: -nan -nan`.

## 8. Remaining issues and next discriminating tests
1. **`RadialBound` bracket failures** (pre-existing; not modified, per the brief). Seen previously at Fe IX μ = 3 (λ = 0.33 a₀), and in this diagnostic whenever the SCF potential is inconsistent. The potential at those failures is finite (v1 dumps).
   - Next test: a full-grid Sturm node count at E = −ε in that failing potential, to separate a genuinely unbound orbital from a bracket artefact.
2. **Ion-sphere truncation** (pre-existing, preserved). With `ups = 0`, `DensityToSZ` integrates only to `ips`; the Fe III diagnostic's nps = 0.01 puts the cut at 6.86 a₀. Fe IX E0 changes by ~0.05 eV between nps = 10⁻² and 10⁻⁶.
3. **Local exchange** in `zx` remains unscreened.
4. **`EEScreenMu` returns 0 for `dps = 0`** (T = 0), i.e. silently unscreened. That convention is shared with the CI Slater integrals and was not changed (outside `DensityToSZ`).
5. **Unexplained 1.7 meV** base–patch difference at Fe IX μ = 10⁻⁴ (5×10⁻⁸ relative); v1 and v2 agree there.

## 9. Isolation
- All work is in `/tmp/fac-debye-scf-test`. `~/Documents/projects/fac_debyescreening` does not exist on this host and was never accessed.
- No install command; global `pfac` unchanged; no remote jobs or processes touched (no `sfac` running).
- No `git reset`, `clean`, checkout or branch change. HEAD is still `c4a016a`. No file outside the worktree changed since it was created.
