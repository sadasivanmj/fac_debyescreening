# Deprince et al. (2020) Table 2 with native SFAC — report

Target: Table 2 of J. Deprince et al., *Plasma environment effects on K lines of
astrophysical interest IV*, A&A 643, A57 (2020): ionization potentials (IPs) of
Fe II–Fe VIII for Debye screening parameters μ = 0, 0.1 and 0.25 a₀⁻¹.

**Status: preliminary cross-code comparison, not an exact reproduction.** All
42 required IPs (7 ions × 3 μ) were computed with native `sfac` for six FAC
models; all runs succeeded. Values labelled *FAC* below were computed here.
NIST and *Deprince* values are quoted from Table 2 of the paper.

Everything is in `fe_calc/deprince_t2/`:

| path | content |
|---|---|
| `deprince_t2_sfac.py` | production driver (native SFAC only, no PFAC) |
| `diagnose_radial_failure.py` | focused diagnostic suite for the screened radial-solver failure |
| `count_csfs.py` | jj-CSF counting of the paper's active space |
| `analyze.py` | model overview, IP-lowering table, error statistics |
| `runs/<model>/<ion>/<mu>/` | `in.sf`, `sfac.log`, `e.lev`, `e.lev.b`, `final.pot`, `run.json` per case |
| `results/ip_comparison.csv` | FAC IPs, reference values, absolute and relative errors (all models) |
| `results/energies.csv` | E0, NELE, NLEV, ground J/parity, λ_D check, binary hash, status per case |
| `results/ip_lowering.csv`, `results/model_overview.md`, `results/summary.md` | derived tables |
| `results/csf_counts.csv` | CSF counts of the active spaces |
| `diagnostics/` | diagnostic runs, `diagnostics.md/json`, grid, warm-restart and small-μ checks |
| `patches/` | source patches used for the final binary |
| `provenance.json` | commits and sha256 of every binary used |

The original `fe_calc/fe_ip_sfac.py` and `fe_calc/fe_ip_sfac_runs/` were left
untouched.

---

## 1. Provenance: which code was actually run

* Repository `/data/project/fac`, branch `fix-debye-ee-screening`, HEAD
  `c4a016a`; the working tree was not modified (only new files under
  `fe_calc/deprince_t2/`). `ee-screening-scf` is at `ff9a5c9`.
* **The installed `~/facinst/bin/sfac` used by `fe_ip_sfac.py` is not a
  `fix-debye-ee-screening` build.** It was built on 2026‑09‑14 04:42 UTC,
  before the upstream merge `49c7b97` (Oct 9), and its `DensityToSZ` calls
  `EEScreenMu` and `sinh`, which only exist on `ee-screening-scf`. A clean
  rebuild of `ee-screening-scf` `ff9a5c9` is **bit-identical** (sha256
  `2789ca41…`) to the installed binary. The preliminary 18.8716 eV result and
  the failing μ = 0.1 run therefore both came from `ee-screening-scf`.
* Builds were made in separate git worktrees, `/data/project/fac_worktrees/ee`
  (c4a016a) and `/data/project/fac_worktrees/scf` (ff9a5c9), each with its own
  install prefix. `~/facinst` and the repository working tree were not
  touched. Unpatched copies of both builds are kept (`inst_unpatched/`).
* **Final binary** for all production runs:
  `/data/project/fac_worktrees/scf/inst/bin/sfac`, sha256 `3a77793e…`. It is
  `ee-screening-scf` ff9a5c9 plus:
  1. `patches/debye_sfac_fixes.patch` (this work, §4): the RadialBound bracket
     for Debye potentials and the opt-in `radial:hxs_iter0` option;
  2. upstream commit `d3a4ec3` (int→size_t overflow fix, #842). It is already
     on `fix-debye-ee-screening` but not on `ee-screening-scf`. Without it,
     `AllocHamMem` segfaults for blocks above ~8192 CSFs (Fe IV, VI and VII in
     the paper's active space). The `structure.c` hunk was applied by hand
     because the surrounding lines differ (`patches/scf_build_working_tree.diff`
     is the exact diff that was built).

  With default options the patches are a no-op for free ions: the patched
  build reproduces the unpatched free-ion E0 to every printed digit
  (diagnostics E1/E2).
* The archive **`fac_debyescreening.zip` does not exist** anywhere on this
  machine (full filesystem search for `*.zip` and `*debye*`), so it could not
  be inspected. The calling conventions were therefore verified against the
  patched source itself (`sfac.c:PPlasmaScreen`, `radial.c:PlasmaScreen`,
  `orbital.c:InitializePS`). They were cross-checked against the project's
  own example inputs (`c4_ground_results/*.sf`), which use the same
  `PlasmaScreen(zps, 0.01, T, 1, 0.0, 1)` form as `fe_ip_sfac.py`.

## 2. What the paper did (from the PDF)

* Hamiltonian (Eq. 1): H = Σᵢ[cα·p + βc² − Z e^{−μrᵢ}/rᵢ] + Σ_{i>j} e^{−μr_ij}/r_ij,
  with μ = 1/λ_D = (4π nₑ/kTₑ)^{1/2} in atomic units (Eq. 2).
* Table 1 reference configurations (`[nl]` = hole):

  | Ion | references | N_CSF |
  |---|---|---:|
  | Fe II | 3d⁶4s, [3p]3d⁷4s, [3p]3d⁶4s², [2p]3d⁷4s, [2p]3d⁶4s², [1s]3d⁷4s, [1s]3d⁶4s² | 12268 |
  | Fe III | 3d⁶, [3p]3d⁷, [3p]3d⁶4s, [3p]3d⁵4s², [2p]3d⁷, [2p]3d⁶4s, [2p]3d⁵4s², [1s]3d⁷, [1s]3d⁶4s, [1s]3d⁵4s² | 31501 |
  | Fe IV–VIII | 3dⁿ, [3p]3dⁿ⁺¹, [2p]3dⁿ⁺¹, [1s]3dⁿ⁺¹ (n = 5…1) | 57189, 81237, 85798, 71135, 43271 |

* Active space: all single and double (SD) excitations from the references
  into the n = 3 and 4s orbitals. MCDHF with the extended average level (EAL),
  weights ∝ (2J+1), including QED; GRASP2K + RATIP; 0 ≤ μ ≤ 0.25.
* NIST values with the uncertainties printed in Table 2: Fe II 16.19920(5),
  Fe III 30.651(12), Fe IV 54.91(4), Fe V 75.0(2), Fe VI 98.985(15),
  Fe VII 124.98(1), Fe VIII 151.060(12) eV. Central values are used for the
  errors; the uncertainties are in the `NIST_unc_eV` column.
* IP(Fe VIII) needs Fe IX. Fe IX is not in Table 1, so this work uses
  [Ar], [3p]3d, [2p]3d, [1s]3d. This is the same pattern as Fe IV–VIII (ground
  plus the 3p, 2p and 1s holes with the electron in 3d).

## 3. How FAC implements the screened Hamiltonian

| piece | where | `fix-debye-ee-screening` (c4a016a) | `ee-screening-scf` (ff9a5c9) |
|---|---|---|---|
| e–n: `orbital:debye_mode 1` | `orbital.c:SetPotentialIPS` | ZPS = Z(r)(1−e^{−r/λ}), added as +ZPS/r in `SetPotentialVT`, so the orbital potential has −Z e^{−r/λ}/r | same |
| e–e in Slater integrals: `radial:ee_screen 1` | `radial.c:GetYk → GetYkDebye` | Yukawa multipole kernel μ(2k+1) i_k(μr<) k_k(μr>) for every Slater integral (all `GetYk` call sites; the cache bypass is intentional) | same |
| e–e in the **SCF direct potential** | `orbital.c:DensityToSZ` | **unscreened Coulomb** | Yukawa-screened (k = 0 kernel, `_shx`, `exp(−μr)`) |
| Breit, QED | — | unscreened | unscreened |

The one-body CI matrix element is εₐδ_ab + ⟨a|−Z/r − Vc − U|b⟩
(`ResidualPotential`). The orbitals satisfy (T + Vc + U + ZPS/r)φ = εφ, so the
ZPS/r terms cancel and the Hamiltonian contains exactly the screened −Z e^{−μr}/r.
With both switches on, both screened interactions of Eq. (1) enter the final
energies on both branches. The branches differ only in the potential used to
generate the orbitals.

**Units and λ_D.** `PlasmaScreen(zps, nps[Å⁻³], tps[eV], m, ups, vxf)`
converts n → n·a₀³ and T → T/27.2114. For m = 1, `InitializePS` sets
λ_D = `dps` = [T/(4πn(1+ups))]^{1/2} a₀. With ups = 0 and
T = λ²·4πn·a₀³·27.2114 eV, the dumped potentials give `dps` = 10.000000 and
4.000000 a₀ for μ = 0.1 and 0.25. The driver checks `dps`·μ = 1 to 10⁻⁶ in
every run. In Debye mode 1, `zps` and `nps` have no other effect on the
energies: zps = 3 instead of 1, or nps = 10⁻⁴ instead of 10⁻⁶, changes E0 by
< 0.1 meV (diagnostics E4, E5). The exception is the ion-sphere cut-off
described next.

## 4. Diagnosis of the screened radial-solver failure

Original symptom (Fe II, μ = 0.1, `fe_ip_sfac.py`):
`Max iteration before finding correct nodes in RadialBound a: 0 0 279 3 2 -8.12838e-154`,
then `Error occured in RadialSolver: -2 3 2`. The failing orbital is n = 3,
κ = 2 (3d₃/₂) in `OptimizeRadial` → `SolveDirac` → `RadialBound`, phase a
(upper energy bracket). The diagnostics table is in `diagnostics/diagnostics.md`.
Each case is a separate native SFAC run with its input, log and FAC
`error.pot` kept.

**Not the reference configurations.** The failure is identical with the ground
configuration alone and with every incremental subset of the Table 1
references (A1–A3). It is not caused by core-hole configurations, so none
were removed.

There are three independent causes:

1. **Ion-sphere truncation of the screened Hartree potential (input level).**
   For m = 1, `InitializePS` still computes the ion-sphere radius
   r_ps = (3 zps/4πn)^{1/3} (5.44 a₀ for nps = 0.01 Å⁻³). When ups = 0,
   `DensityToSZ` integrates the density only out to that grid point
   (`ips = 278`) and copies the value outward. The screened Hartree charge u(r),
   which should decay like e^{−μr}, is frozen at 14.51 beyond 5.44 a₀. The
   electron therefore sees r·V → +14.5 at large r, a repulsive Coulomb tail
   with V > 0 from r ≈ 5.9 a₀ (C1). *Fix:* use nps = 10⁻⁶ Å⁻³ with T scaled to
   keep λ_D (r_ps = 117 a₀; r·V → 0). In the final protocol the grid is built
   in the free-ion step, so `ips = 0` and no truncation occurs at all.
2. **RadialBound assumes infinitely many bound states (code level).** Phase a
   raises E until the turning-point-limited solution has nr+1 nodes, i.e. until
   it reaches the next state of the same κ. The step is capped at ½|E|, so E
   creeps geometrically to 0⁻ (hence E = −8×10⁻¹⁵⁴). A Yukawa potential has
   finitely many bound states, so 4d₃/₂ may not exist even when 3d₃/₂ is bound.
   A first fix that used the turning-point node count at E = −ε gave **false
   negatives**: as E→0 the turning point stays finite (~3.5 a₀ for l = 2), and
   the outermost node of a weakly bound state can lie beyond it.
   *Fix (`patches/debye_sfac_fixes.patch`, `orbital.c`):* only when `mps == 1`
   and E > −ε (ε = `orbital:debye_emax_eps`, default 10⁻⁶ Ha), count the
   nodes of the regular solution over the **whole grid**. By Sturm's theorem
   this counts the bound states below −ε. If the count is exactly nr+1, use
   −ε as the upper bracket and skip the phase-5 upward search (which would
   step past 0). Otherwise return error −6 with
   "No Debye bracket in RadialBound (unbound | higher state bound)". The final
   log-derivative matching bisection is unchanged. The patch never touches
   non-Debye runs.
3. **SCF start-up transients (protocol level).** With the corrected solver, a
   cold start still fails, now with a *genuine* unbound 3d₃/₂ at SCF
   iteration 2 (D1). In the screened cold start, 3d₃/₂ jumps from −0.23 Ha
   (iteration 0) to −6.9 Ha (iteration 1), and the next potential then
   over-screens. The free-ion SCF goes through the same overshoot (maximum
   orbital-energy changes of 33, 20, 8 Ha in iterations 1–3), but its Coulomb
   tail keeps the valence orbitals bound, whereas a Yukawa tail does not. A warm start from converged orbitals still fails (D6) because
   `OptimizeILoop` switches exchange off at iteration 0
   (`hxs = hxs0·(1−e^{−0.75·iter})`). The exchange-free potential with no
   self-interaction removal does not bind 3d for Fe II at μ ≳ 0.09 (verified
   with the full-grid count, so this is not a solver artefact). Converged
   3d₃/₂ energies show the state is firmly bound: −0.594, −0.562, −0.529,
   −0.506, −0.436, −0.330 Ha at μ = 0.02, 0.05, 0.08, 0.10, 0.16, 0.25.
   Stronger damping (`SetOptimizeControl` stabilizer 0.5–0.9) does not help.
   *Fix:* (a) start every run with a free-ion `OptimizeRadial` and reach the
   target μ through warm-started steps of 0.05; (b) the opt-in option
   `radial:hxs_iter0` (patch, default 0) starts the exchange ramp as if that
   many iterations had passed. Value 8 is used, giving 99.75% exchange at the
   first warm iteration. FAC does not accept convergence while the ramp is
   active, so the fixed point keeps full exchange.

Not used as a fix, and why:
* Lowering μ, raising iteration limits, suppressing the error, dropping
  references: not needed, and they would not address any of the causes.
* e–n screening alone: the SCF Hartree tail is then unscreened, so an
  electron sees r·V → +(N−1) = +24 (B1). Even 3s fails. It is a diagnostic,
  not a usable model.
* `hxs_iter0` for **cold** starts changes the free-ion E0 by up to 0.30 eV
  (value 20), showing that FAC's SCF has some history dependence when started
  from its crude initial potential. It is therefore only used for warm steps.
  There the history dependence was measured to be ≤ 0.011 eV (§6).

`fix-debye-ee-screening` (unscreened `DensityToSZ`) cannot run the screened
SCF for Fe II–Fe V even with these fixes. The orbital potential then has an
asymptotic charge of +(N−1) (F1). It does converge for Fe VI–Fe IX at
μ = 0.1, and there the IPs differ from `ee-screening-scf` by only 0.02–0.19 eV
(e.g. `ref` Fe VI 81.739 vs 81.860 eV). The CI Hamiltonian is the same on both
branches; only the orbitals differ (valence ε shifted by ≈ +1.8 Ha). Free-ion
energies are identical on the two branches (all ions), so the upstream master
changes merged into `fix-debye-ee-screening` do not affect these energies.
`diagnostics/branch_ee/` contains the full comparison.

Side finding (not relevant for Fe II–IX): for a **one-electron ion** the e–n
Debye shift vanishes (Fe XXVI 1s: ΔE = 0.00009 eV instead of Zμ = 0.0708 eV
at μ = 10⁻⁴). `PotentialHX` zeroes `ZPS` whenever N < 1+ε
(`diagnostics/first_order/`).

## 5. Computational protocol and models

Every (model, ion, μ) is one `sfac` process in `runs/<model>/<ion>/<mu>/`:

1. `Config` groups (names `ref`, `g0`, `sd`, all ≤ 16 characters); occupations
   are written out in full, e.g. `[2p]3d⁷4s` → `1s2 2s2 2p5 3s2 3p6 3d7 4s1`.
2. Cold `OptimizeRadial` for the free ion, which also fixes the radial grid.
   All μ for an ion share the same grid.
3. μ > 0 only: `orbital:debye_mode 1`, `radial:ee_screen 1`,
   `radial:hxs_iter0 8`, then `PlasmaScreen(z_net, 1e-6, T(μₖ), 1, 0.0, 1)`
   with `OptimizeRadial` for μₖ = 0.05, 0.10, …, μ.
4. `Structure` (restricted to the ground parity/2J block for the SD models),
   `PrintTable(..., 1)`. E0 = energy of the lowest level, which FAC prints
   **in eV** (`dbase.c`: energy × HARTREE_EV). No unit conversion is applied.

There is no `ConfigEnergy` step. The paper has no equivalent, and for the
ground level it changes E0 by 0.0001 eV (diagnostic E3 reproduces the original
baseline with `ConfigEnergy`: −34561.2841 eV).

| model | orbitals optimized on | CI space |
|---|---|---|
| `ref` | statistical-weight average of all Table 1 references (FAC's analogue of an EAL over the references; FAC weights configurations by their degeneracy = Σ(2J+1)) | references, all J, both parities |
| `sdv` | as `ref` | + SD of the ground reference within {3s,3p,3d,4s}; ground (P, 2J) block |
| `sdp` | as `ref` | + SD of **every** reference into {1s,…,3d,4s} (the paper's active space); ground (P, 2J) block |
| `gref`, `gsdv`, `gsdp` | ground reference configuration only | as `ref`, `sdv`, `sdp` |

An IP is reported only if both runs succeeded and share the binary hash,
Debye settings, nps, ups, λ_D and `hxs_iter0`. Otherwise it is FAILED or
NOT RUN. Restarts reuse a case only if `run.json` says ok with the same input
hash and binary hash, and `e.lev` is newer than the run start. Otherwise the
outputs are deleted before the rerun.

### Active space and CSF counts

`count_csfs.py` enumerates the relativistic subshell occupations of every
nonrelativistic configuration and counts jj-coupled CSFs per J, as
N(J) = c(M = J) − c(M = J+1). FAC builds the same jj-coupled basis, so FAC's
level count equals this CSF count over the same blocks. This was checked: the
counted CSFs of the reference configurations equal FAC's NLEV for every ion
(986, 1640, 460, 536, 442, 267, 108, 29). A FAC level is one eigenvector of a
(J, parity) block. The number of levels equals the number of CSFs in the
blocks computed, not the number of nonrelativistic configurations, and FAC
counts each J once (not each M).

| Ion | published N_CSF | paper AS: configs | CSFs all J, both P | CSFs all J, ground P | CSFs ground (P,2J) = `sdp` block | SD of ground ref. within {3s,3p,3d,4s}, ground block¹ | references, all J (= FAC NLEV of `ref`) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Fe II | 12268 | 120 | 15147 (1.23×) | 7603 | 927 | 107 | 986 |
| Fe III | 31501 | 153 | 40189 (1.28×) | 20305 | 3363 | 267 | 1640 |
| Fe IV | 57189 | 153 | 75309 (1.32×) | 38157 | 8239 | 470 | 460 |
| Fe V | 81237 | 153 | 109657 (1.35×) | 55665 | 3008 | 151 | 536 |
| Fe VI | 85798 | 153 | 118162 (1.38×) | 60014 | 10591 | 365 | 442 |
| Fe VII | 71135 | 153 | 99555 (1.40×) | 50575 | 10083 | 280 | 267 |
| Fe VIII | 43271 | 152 | 61391 (1.42×) | 31173 | 5964 | 125 | 108 |
| Fe IX | — | 146 | 28697 | 14565 | 963 | 23 | 29 |

¹ FAC's `sdv` block additionally contains the even-parity [1s]-hole references, e.g. NLEV = 126 instead of 107 for Fe II. FAC's NLEV for the completed `sdp` runs equals the counted ground-block CSFs exactly (927, 3363, 8239, 3008, … ).

The paper's recipe as stated (SD from all references into n = 3 and 4s, all
J, both parities) gives 1.23–1.42 times the published counts. The published
numbers must include restrictions that the paper does not state (for example
on J, or on which shells may be excited). They could not be reproduced. For an
IP only the ground (J, parity) block matters, and in that block the paper's
active space has 927–10 591 CSFs. That is tractable, so `sdp`/`gsdp` include it
completely. What remains different from the paper is the orbital treatment:
FAC uses one local central potential optimized on an average configuration and
then CI. The paper used MCDHF-EAL over the full CSF expansion, with orbitals
variationally optimized on it.

## 6. Validation

* **Baseline.** Fresh `ref` runs give E0(Fe II) = −34561.2842 eV and
  E0(Fe III) = −34542.4125 eV, so IP = 18.8717 eV. This matches the
  preliminary −34561.2841 / −34542.4125 / 18.8716 eV (the 0.1 meV difference
  is the `ConfigEnergy` step).
* **Electron counts.** NELE from every `e.lev` equals 25…18 for Fe II…Fe IX.
  This is enforced in the driver, and every configuration is checked when it
  is generated.
* **Adjacent charge states.** IP(Fe X) = E0(Fe X+1) − E0(Fe X), with the same
  model, μ and binary for both ions; IP(Fe VIII) uses Fe IX.
* **Both screenings active.** `mps = 1` and `dps = 1/μ` are read from every
  final potential. e–e screening in the SCF shows up as r·V → 0 (versus +24
  with e–n alone). In the CI, the first-order tests below are only reproduced
  when both terms are present.
* **First-order Debye–Hückel limit.** With e^{−μr}/r ≈ 1/r − μ, the energy
  shift is ΔE = μ[ZN − N(N−1)/2]. FAC gives this to 0.02% for He-like and
  Li-like Fe and 0.06% for Ne-like Fe at μ = 10⁻⁴ and 10⁻³, which verifies the
  e–n and e–e implementation. For Fe II–IX the agreement is 0.2–1.2%. The
  remaining deviation is linear in μ, so it is first order. FAC orbitals are
  not variational for the ground level, so the orbital response contributes
  at first order; this is largest for the EAL-like `ref` orbitals. The
  resulting slope dIP/dμ is −(z+1) within 0–20% for Fe IV–VIII, but not for
  `ref` Fe II/III (`diagnostics/small_mu_check.txt`).
* **Large finite λ_D vs the free path.** μ = 10⁻⁴ (λ = 10⁴ a₀) differs from the
  free-ion run by the expected first-order shift (Fe II ground: +0.941 eV vs
  +0.952 eV analytic). There is no constant offset between the Debye and
  free code paths.
* **Radial grid.** 1500 → 3000 points changes IPs by ≤ 0.023 eV
  (`diagnostics/grid/grid_sensitivity.txt`).
* **SCF history.** One or three extra warm `OptimizeRadial` calls at μ = 0
  change E0 by ≤ 0.011 eV and IPs by ≤ 0.011 eV
  (`diagnostics/warm_restart_check.txt`).
* **Monotonicity.** IP(0) > IP(0.1) > IP(0.25) holds for every ion and model
  (`results/summary.md`).

RESULTS_PLACEHOLDER
