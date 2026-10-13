# Debye–Hückel screened SCF direct potential in `DensityToSZ`: implementation and validation

Worktree `/tmp/fac-debye-scf-test`, detached at `c4a016a` (fix-debye-ee-screening), created 2026-10-12 on the Linux host `luna`. The brief's macOS paths do not exist here; Linux use was authorised by the user (option 2). Starting revision chosen because it has the Yukawa Slater integrals (`EEScreenMu`, `GetYkDebye`) but an unscreened `DensityToSZ`.

## 1. Initial source findings (c4a016a)
- `DensityToSZ(pot,d,z,zx,jps,md)` (`faclib/orbital.c:4753`) computes `z = ∫₀ʳ d dr′ + r∫ᵣ^∞ (d/r′)dr′`, i.e. **z = r·V_direct(r)**. `jps = ∫ d/r dr = V(0)` and the return value is `im`.
- `d` is the **radial** density: the dumps give ∫d dr = 18.0000 for Fe IX and 2.0000 for He-like. `nx = zx/(4πr²)` converts to a volume density. The origin term `d[0]r[0]/3` assumes d ∝ r².
- Grid: ρ = a·r^q + b·ln r with a step of exactly 1 per index (`orbital.c` ~3455: `del = 1` by construction of `b`). `dr_drho = dr/dρ`. `NewtonCotesIP` → `uvip3s` integrates with unit index spacing, so Σ f·dr_drho is ∫f dr. `id ≥ 0` accumulates forward from `r[i0]`, `id < 0` backward from `r[i1]`, and `m = −1` means no swap.
- `zx`: on input the exchange density (when `ahx ≠ 0`); on output r·V_xc (local exchange). Downstream (`PotentialHX1`) uses `u − ue`.
- Callers: `PotentialHX1` with `md = 0` (bound density `NPS`), `PotentialHX1` with `md = 1` (plasma `VPS`, skipped when `mps = 1`), and `AverageAtom` (`mps ≥ 3`). Only the `md = 0` Debye-mode call is affected.
- `EEScreenMu()` (`radial.c:10653`) returns 1/`dps` only if `radial:ee_screen = 1`, `mps = 1` and `dps > 0`, otherwise 0. `dps` = [T/(4πn(1+ups))]^{1/2} a₀ is set in `InitializePS`. It reads the global `potential`, which is the object every caller passes. Using it ties the SCF screening to the same switch as the CI Slater integrals. It is not declared in any header, so an `extern` is needed.
- Historical reference: branch `ee-screening-scf` commit 4dd4b40 has a screened `DensityToSZ` with a `_shx` helper (series for x < 10⁻⁶, argument clamped at 400). It was used only for comparison, not merged.

## 2. Derivation
The angular average of e^{−μ|r−r′|}/|r−r′| is sinh(μr<)e^{−μr>}/(μ r r′). With d = 4πr²ρ:

r·V_DH(r) = ∫ d(r′) sinh(μr<)e^{−μr>}/(μ r′) dr′ = **e^{−μr} ∫₀ʳ d·sinhc(μr′) dr′ + r·sinhc(μr) ∫ᵣ^∞ (d/r′) e^{−μr′} dr′**

At μ → 0 both kernels reduce to the original Coulomb expression; all factors are dimensionless or lengths, consistent with z. The candidate's formula is therefore correct. Its arithmetic is not: sinh(x) overflows for x ≳ 710, and where d = 0 the inner integrand becomes 0·∞ = NaN.

Implementation: `_ee_wsinhc(w,x)` = w·sinhc(x). It uses a series for x < 10⁻³, uses exp(ln|w| + x − ln 2x) for x ≥ 700 (exact to double precision), and returns 0 exactly when w = 0. The outer term is written as `_ee_wsinhc(r·B, μr)`. The μ = 0 branch keeps the original statements **verbatim**. Nothing else changed: the exchange-correlation, `ips`/`ups`/`vxs` and extension branches are untouched.

## 3. Changes
- Only file changed: `faclib/orbital.c`, 42 inserted lines and 0 deleted: the helpers plus `extern` before the function, `double mu = EEScreenMu();`, and an `if (mu > 0) {…} else {original}` block.
- Applied by `patch_density.py`. It checks location, signature uniqueness and the sha256 of the original function (`0c098011…`), refuses double application, and re-verifies the result.
- Backup: `faclib/orbital.c.pre_debye_scf`. Reversible diff: `patch_density.diff`.
- A temporary, env-guarded dump (marker `TEMP_DIAG_D2SZ`) was used only in the separate build `_diag_pfac`. It was removed; the clean rebuild is byte-identical to the pre-diagnostic patched build (sha256 `a18f80ec…`). A stale-object trap was caught: `cp -p` kept the old mtime, so `make` initially relinked the diagnostic object until the source was `touch`ed.

## 4. Build
- `./configure --prefix=/tmp/fac-debye-scf-test/_noinstall` then `make pfac` (rc 0 for baseline, patched, diag and final; no `orbital.c` warnings). No `make install*`, no `pip install`.
- The conda `fac_env/lib` contains no `libfac.a`, so the extension links the worktree's `libfac.a`.

## 5. Runtime verification
- Each test asserts `pfac.fac.__file__` is in `_baseline_pfac/` or `_patched_pfac/` (`PYTHONPATH` set per subprocess). A global `pfac` exists in `fac_env` site-packages and was never imported or modified (mtime 2026-08-29).
- The compiled `DensityToSZ` calls `EEScreenMu` and `_ee_wsinhc` (patched build); the baseline makes no such calls.
- Screening was active in every screened run: `final.pot` has `mps = 1` and `dps = 10.0 / 4.0 a₀` for μ = 0.1 / 0.25.

## 6. Results (`tests/out/r1…r6`, `tests/test_helpers.out`)
| Test | Result | Evidence | Status |
|---|---|---|---|
| T1 unscreened regression (He-like, Fe IX, Fe III free; Fe IX with ee_screen=1 and no plasma) | base vs patch bitwise identical | all 18 `final.pot` columns × 1500 points: max abs diff 0; `e.lev` (excl. TSess) and headers identical | **PASS** |
| T2 minimal screened controls (He-like, Fe IX; μ = 0.1, 0.25) | all stages complete on both builds | patched Fe IX 3p₃/₂ shifts −5.70 → −7.32 Ha (≈(N−1)μ, unscreened tail removed); patched E0 lower by 0.07 / 0.71 eV; He-like ΔE matches first + second order | **PASS** |
| T3 kernel vs independent Yukawa double integral (same density) | agrees to quadrature accuracy | max\|z_FAC − z_ref\|/max\|z\| = 4.8e-5 (μ = 0 calibration) … 9.5e-5 (μ = 1), 1.7e-4 (μ = 3); μ = 10⁻⁸ screened path differs from reference by 8.56e-4 abs, identical to the Coulomb calibration | **PASS** |
| T5a helper unit test | max rel err 6.1e-14 vs log-domain reference, x ∈ [0, 10⁴] | 0 non-finite where reference finite; candidate: 7 (e.g. 0·sinhc(720) = NaN) | **PASS** (candidate FAIL) |
| T5b large μr in a real SCF (Fe IX, μ = 3, full grid) | `DensityToSZ` finite, 224 points in the x ≥ 700 branch | the run then fails in **RadialBound a** (3s, λ_D = 0.33 a₀) | kernel **PASS**; run fails in the solver |
| T5c small μ energies | He-like 10⁻⁴, 10⁻³: both builds equal first order to 1e-4 eV | Fe IX 10⁻⁸: both 3–5e-4 eV from first order (grid offset from rps); Fe IX 10⁻⁴: patch − base = −1.7 meV | PASS; 1.7 meV unexplained |
| T5d Fe IX μ = 1.0 | base fails (RadialBound a, 3s); patch completes (E0 −26689.18 eV) | | **PASS** |
| T4 existing Fe III diagnostic settings (3d⁶, μ = 0.1, nps = 0.01, zps = 2, cold, ConfigEnergy) with ee_screen = 1 | **base fails** at ConfigEnergy(0), SCF iteration 2: `RadialBound a` 3s, grid point 230, E = −8.7e-154. **Patch completes** all stages, E0 = −33625.6792 eV, all orbitals bound (3d −0.77 Ha) | with ee_screen = 0 (original script) both fail identically, as expected | **PASS** for the screened direct potential |

## 7. Remaining problems
1. **`RadialBound` bracket failure (pre-existing, separate).** It still occurs when an orbital is near or above threshold, or in a cold-start transient (seen at Fe IX μ = 3). Its cause is the upper bracket requiring the next nodal state, which a Yukawa potential may not have. This was not changed here (no unrelated edits).
   - Next discriminating experiment: a full-grid Sturm node count at E = −ε for that failing case, to separate "genuinely unbound" from "bracket artefact".
2. **Ion-sphere truncation (pre-existing, preserved).** `DensityToSZ` still integrates only to `ips` when `ups = 0`; with nps = 0.01 the cut is at 6.86 a₀ (Fe III). The Fe III success above carries that artefact: E0 differs by ~0.05 eV between nps = 10⁻² and 10⁻⁶ for Fe IX.
3. **Local exchange is unscreened.** `zx` is a μ-independent local XC potential; only the direct term is Yukawa-screened.
4. **A 1.7 meV difference** at Fe IX μ = 10⁻⁴ between builds (5×10⁻⁸ relative) is not explained.

This is a code-path validation. It does not reproduce published Fe IPs.

## 8. Isolation
- All work is in `/tmp/fac-debye-scf-test`. `~/Documents/projects/fac_debyescreening` does not exist on this host and was never touched.
- No installation; global `pfac` unchanged.
- No remote access, no process signals; the Linux production results and the main repository working tree are unchanged.
- The only change outside the worktree is git's worktree metadata (`/data/project/fac/.git/worktrees/fac-debye-scf-test`). No `reset`, `clean` or checkout in any tree.
