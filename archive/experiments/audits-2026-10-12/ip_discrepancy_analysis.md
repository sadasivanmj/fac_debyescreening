# Fe II–Fe VIII IP discrepancies vs Deprince et al. (2020): read-only analysis

Prepared 2026-10-12, ~13:05–13:20 UTC. No SFAC run, rebuild, edit, or process signal was made for this analysis (`pgrep -x sfac`: none running at 13:03 UTC).
Paths: `D = /data/project/fac/fe_calc/deprince_t2`, `A = /tmp/claude-1000/-data-project-fac/a8e208a0-910c-4f2e-b650-3503c715a8c5/scratchpad/audit`, `SRC = /data/project/fac_worktrees/scf` (the source tree of the production binary).

Labels: **[Verified]**, **[Supported interpretation]**, **[Untested hypothesis]**, **[Unresolved]**.

## 1. Executive Finding

1. **The Fe III → Fe IV anomaly in `ref`, `sdv` and `sdp` (≈19.7 eV vs 32.93 eV) comes from the orbital-optimization average, not from the CI space or the bookkeeping. [Verified]** The CI configuration sets are identical between `ref`/`gref`, `sdv`/`gsdv` and `sdp`/`gsdp`. Changing only the optimization group from the reference average to the ground configuration raises the Fe III IP by **+8.44 eV** (`ref`→`gref`) or **+8.05 eV** (`sdp`→`gsdp`) at μ = 0. That is 61–64% of the −13.2 eV gap.
2. **Correction to earlier statements:** FAC's reference average is **not** a statistical-weight average. Within a group, FAC weights every **jj configuration equally** (SRC `faclib/config.c`, `GetAverageConfig`: weight `1/n_cfgs` per configuration; `sweight` weights only groups, and only for `acfg_wmode=2`). This model reproduces the average occupancies printed in every `D/runs/ref/<ion>/mu0/final.pot` to ≤6×10⁻⁶. For Fe III the actual average has **ground weight 6.7% (5 of 75 jj configurations) and ⟨4s⟩ = 1.000**; the ground configuration 3d⁶ has no 4s. My earlier "1.9% / ⟨4s⟩ = 1.14" were wrong. The wrong description also appears in `D/REPORT.md`, the driver docstring, `A/claude_change_history.txt` and my previous audit answer; none of those were edited. **[Verified]**
3. **The averaged orbitals over-stabilize 4s-bearing configurations. [Verified observation; causal role = Supported interpretation]**
   - In `sdp` (averaged orbitals), the lowest Fe III 3d⁵4s level of the 2J = 8 block is only 0.058 eV above 3d⁶ at μ = 0 and becomes the lowest level of that block at μ = 0.1 and 0.25.
   - In `gsdp` it stays 9.06 / 7.38 / 4.78 eV above 3d⁶.
   - Fe III is the only ion whose n = 3 eigenvalues are shifted strongly by the averaging (3d₅/₂ by −0.254 Ha, versus −0.06 to −0.11 Ha for Fe IV–VIII).
4. **New protocol defect in 6 of 144 runs: `sdv` and `sdp` for Fe II at μ = 0.25 and Fe III at μ = 0.1 and 0.25. [Verified]** The lowest computed level changes configuration character (3d⁵4s² for Fe II, 3d⁵4s for Fe III). Because the SD models compute only the expected (parity, 2J) block, E0 for these runs may not be the true ground level: for Fe III, 3d⁵(⁶S)4s ⁷S₃ has 2J = 6. The affected IPs are `sdv`/`sdp` Fe II (μ = 0.1, 0.25) and Fe III (μ = 0.1, 0.25). The previous audit's 2J check did not catch this because the 2J stayed the same. **[Verified]** Whether a lower level exists in another J block is **[Untested hypothesis]**.
5. **The ground-optimized gap (`gsdp` 27.77 vs 32.93 eV, −5.16 eV) splits arithmetically into FAC − NIST = −2.885 eV and NIST − paper = −2.279 eV. [Verified arithmetic]** FAC lies below experiment and the paper lies above it, by similar amounts. Fe II behaves the same way: −0.694 and −2.691 eV. Attributing FAC's −2.9 eV to missing differential 3d correlation is **[Supported interpretation]**, not demonstrated.
6. **Screening-induced change ΔIP(μ) for Fe IV–VIII: |ΔIP| is 0.8–2.9% smaller than the paper's, in every model, with a model spread ≤ 0.4 percentage points. The mismatch decreases with ionic charge. [Verified]** A two-point fit puts the mismatch mainly in the term linear in μ: FAC's coefficient is 0.978–0.997 × the exact first-order value −(z+1)μ; the paper's is 1.003–1.011 ×. **[Verified arithmetic]** The mechanism is **[Unresolved]**: non-variational orbital response, the Coulomb local exchange in the screened SCF potential, and the paper's own EAL non-stationarity cannot be told apart with the existing data.
7. **The baseline error and the screening shift are separable. [Verified]** For Fe IV–VIII, models whose μ = 0 IPs differ by up to 1.4 eV give the same ΔIP to within 0.08 eV. For Fe II/III, ΔIP is strongly model-dependent (Fe III μ = 0.25: −16.06 to −19.63 eV vs −16.51 published). This tracks the 3d^n / 3d^{n−1}4s near-degeneracy.

## 2. Evidence and Data Provenance

| Artifact | Status / use |
|---|---|
| `A/claude_change_history.txt` | read; contains the incorrect weight description noted in §1.2 |
| `A/ip_audit.csv`, `A/raw.json` (12:48 UTC) | read; 126 IP rows, all `ok`; source of the E0 values below **[SUB: independent subtraction from raw e.lev]** |
| `D/results/{energies.csv,ip_comparison.csv,summary.md}` (12:48 UTC), `D/results/csf_counts.csv`, `D/diagnostics/diagnostics.md` | present, read |
| `D/runs/<model>/<ion>/<mu>/{in.sf,e.lev,final.pot,run.json}` | read directly in this analysis (average occupancies, eigenvalues, ground-level labels, per-configuration lowest levels) |
| `D/deprince_t2_sfac.py`, `D/patches/debye_sfac_fixes.patch`, `D/patches/scf_build_working_tree.diff` | read |
| `D/diagnostics/small_mu/runs/*/*/*/run.json` | read (μ = 10⁻⁴, 10⁻³; binary 3a77793e) |
| SRC `faclib/config.c` (`GetAverageConfig`), `faclib/orbital.c` (`SetPotentialIPS`, `DensityToSZ`), `faclib/radial.c` (`EEScreenMu`, `GetYkDebye`, `PotentialHX`) | read |

Provenance of all 144 production runs: binary sha256 `3a77793e4b5e6f2c…`, commit `ff9a5c9cc6623e3f6adb8062ba23a52e9aace582`, source-diff sha256 `ceac162e…` (all 144 `run.json` files; re-verified 12:53 UTC). **[Verified]** The analysis numbers were produced by read-only Python in this session. No helper script was saved; a temporary stdout capture `A/.analysis_stdout.tmp` was created and deleted within a minute.

Model definitions, checked against `in.sf` **[Verified]**:
- `ref`/`sdv`/`sdp` call `OptimizeRadial(['ref'])` (all Table-1 refs in one group).
- `gref`/`gsdv`/`gsdp` call `OptimizeRadial(['g0'])` (ground configuration only).
- The configuration union and the symmetry restriction are identical within each orbital pair, for all 8 ions.
- SD models call `Structure(p,[2J])` before `Structure`.
- Fe III config counts: `ref` 10, `sdv` 21, `sdp` 78 (`sdv` ⊂ `sdp`).
- No `ConfigEnergy` call anywhere.

## 3. Fe III → Fe IV: Quantitative Decomposition

IPs at μ = 0 (eV), from raw E0 subtraction **[Verified]**:

| model | IP | − Deprince (32.93) | − NIST (30.651) |
|---|---:|---:|---:|
| ref | 19.6518 | −13.278 | −10.999 |
| sdv | 19.7269 | −13.203 | −10.924 |
| sdp | 19.7215 | −13.209 | −10.930 |
| gref | 28.0942 | −4.836 | −2.557 |
| gsdv | 27.7852 | −5.145 | −2.866 |
| gsdp | 27.7663 | −5.164 | −2.885 |

Raw E0 (eV, μ = 0):
- Fe III: ref −34542.4125, gref −34555.7544, sdp −34543.3083, gsdp −34556.1514.
- Fe IV: ref −34522.7607, gref −34527.6602, sdp −34523.5868, gsdp −34528.3851.

Decomposition (μ = 0) **[Verified arithmetic]**: ref 19.652 → (orbital average → ground-optimized, CI fixed) **+8.442** → gref 28.094 → (CI refs → refs+SD) **−0.328** → gsdp 27.766 → residual **+5.164** to the paper = +2.885 (NIST − gsdp) + 2.279 (paper − NIST). The alternative path ref → sdp (CI, +0.070) → gsdp (orbitals, +8.045) gives the same total. The two effects are not additive: CI captures −0.90 eV of correlation with averaged orbitals and −0.40 eV with ground orbitals (E(sdp)−E(ref) vs E(gsdp)−E(gref)), so CI partially compensates for the poor orbitals.

Orbital penalty P = E0(avg-optimized) − E0(ground-optimized), same CI, μ = 0 / 0.1 / 0.25 (eV) **[Verified]**:

| ion | P_ref | P_sdp | avg ground weight | avg ⟨4s⟩ (ground) | avg 3d excess | 3d₅/₂ ε shift (Ha) |
|---|---|---|---:|---:|---:|---:|
| Fe II | 9.99 / 10.59 / 10.95 | 9.83 / 10.42 / 8.87 | 10.0% | 1.50 (1) | +0.40 | −0.142 |
| **Fe III** | **13.34 / 14.56 / 14.98** | 12.84 / 13.15 / 10.89 | **6.7%** | **1.00 (0)** | −0.07 | **−0.254** |
| Fe IV | 4.90 / 5.58 / 5.74 | 4.80 / 5.43 / 5.59 | 16.7% | 0 (0) | +0.83 | −0.059 |
| Fe V | 4.45 / 5.11 / 5.24 | 4.35 / 4.95 / 5.08 | 16.7% | 0 | +0.83 | −0.072 |
| Fe VI | 4.21 / 4.90 / 5.01 | 4.12 / 4.73 / 4.84 | 13.8% | 0 | +0.86 | −0.087 |
| Fe VII | 3.85 / 4.53 / 4.62 | 3.77 / 4.36 / 4.44 | 13.0% | 0 | +0.87 | −0.100 |
| Fe VIII | 3.59 / 4.25 / 4.33 | 3.52 / 4.07 / 4.14 | 11.8% | 0 | +0.88 | −0.109 |
| Fe IX | 3.39 / 4.06 / 4.13 | 3.33 / 3.88 / 3.94 | 9.1% | 0 | +0.91 | — |

Weights and occupancies come from the equal-jj-weight model, validated against `final.pot` (number of jj configurations per Table-1 reference: Fe III [5,8,10,10,8,10,10,4,5,5] = 75; Fe IV [5,10,10,5] = 30). Eigenvalue shifts are ref − gref from the `final.pot` headers.

Answers to the six questions:
1. *Why is Fe III penalized more than Fe IV?*
   - **[Verified]** Only Fe III's reference average combines a very small ground weight (6.7%) with a full 4s electron that is absent from its ground configuration. Its n = 3 eigenvalues are 4× more displaced than Fe IV's.
   - **[Supported interpretation]** Moving ~1 electron from 3d into 4s (plus ~0.93 core holes) leaves the 3d less screened and more contracted than optimal for 3d⁶, giving a large variational penalty. Fe IV–IX keep the extra electron in 3d, and their penalty falls smoothly from 4.9 to 3.4 eV with charge.
2. *How do the different references act?* Through the average potential only. The core-hole references enter the CI but do not change E0: ground-optimized runs with and without the references in the CI give identical E0 (`D/diagnostics/scratch_feiii/*_gonly` vs `*_optg_allref`, earlier session; and gref = ground-config-only E0). **[Verified]**
3. *Are the 4s-bearing core-hole references responsible?* **[Supported interpretation / Untested hypothesis]**. They are the only compositional difference between Fe III and Fe IV–IX, and the sdp level structure shows 4s over-stabilization (§1.3). Their effect has not been isolated from the effect of the small ground weight; §8 does that.
4. *Weights, orbitals, CI, or interaction?* The penalty is an orbital effect whose size is set by the jj-equal weights **[Verified: CI set fixed in the comparison]**. CI interacts with it: it partly compensates at μ = 0 and changes the ground character at μ > 0 in `sdv`/`sdp` **[Verified]**.
5. *Why is gsdp still 5.16 eV low?* 2.885 eV is FAC's own deficit relative to experiment and 2.279 eV is the paper's excess **[Verified arithmetic]**. FAC's deficit is consistent with missing differential 3d correlation (no 4p/4d/4f in either active space) and single-configuration orbitals **[Supported interpretation]**.
6. *Paper's restrictions, energy convention, relativity, QED?* Both codes use the Dirac–Coulomb Hamiltonian; FAC includes Breit and QED by default, and QED/Breit change Fe III/IV IPs by much less than 1 eV **[Supported interpretation]**. The paper's CSF selection cannot be reconstructed (§4). These cannot explain a 13 eV effect. Whether they explain part of the 2.28 eV paper − NIST excess is **[Unresolved]**.

## 4. Orbital Optimization Versus Configuration-Space Effects

Orbital effect, ΔIP = IP(avg) − IP(ground-opt) with identical CI, μ = 0 / 0.1 / 0.25 (eV) **[Verified]**:

| ion | ref − gref | sdv − gsdv | sdp − gsdp |
|---|---|---|---|
| Fe II | +3.35 / +3.96 / +4.03 | +3.02 / +2.74 / +2.02* | +3.01 / +2.73 / +2.02* |
| Fe III | −8.44 / −8.98 / −9.24 | −8.06 / −7.73* / −5.31* | −8.05 / −7.72* / −5.30* |
| Fe IV | −0.45 / −0.47 / −0.49 | −0.46 / −0.48 / −0.52 | −0.45 / −0.48 / −0.51 |
| Fe V | −0.24 / −0.21 / −0.23 | −0.24 / −0.23 / −0.25 | −0.23 / −0.22 / −0.25 |
| Fe VI | −0.36 / −0.37 / −0.39 | −0.36 / −0.38 / −0.40 | −0.35 / −0.37 / −0.39 |
| Fe VII | −0.26 / −0.28 / −0.30 | −0.27 / −0.30 / −0.31 | −0.25 / −0.29 / −0.30 |
| Fe VIII | −0.20 / −0.19 / −0.20 | −0.20 / −0.21 / −0.22 | −0.18 / −0.19 / −0.21 |

CI effect, same orbital protocol (eV) **[Verified]**:

| ion | sdv − ref | sdp − ref | gsdv − gref | gsdp − gref |
|---|---|---|---|---|
| Fe II | −0.34 / −1.25* / −2.03* | −0.35 / −1.26* / −2.04* | −0.01 / −0.03 / −0.03 | −0.02 / −0.03 / −0.03 |
| Fe III | +0.08 / +0.93* / +3.64* | +0.07 / +0.92* / +3.63* | −0.31 / −0.31 / −0.29 | −0.33 / −0.33 / −0.31 |
| Fe IV | −0.72 / −0.74 / −0.71 | −0.75 / −0.77 / −0.74 | −0.71 / −0.72 / −0.69 | −0.75 / −0.76 / −0.73 |
| Fe V | −0.70 / −0.72 / −0.69 | −0.74 / −0.76 / −0.73 | −0.70 / −0.70 / −0.68 | −0.74 / −0.75 / −0.72 |
| Fe VI | −0.87 / −0.88 / −0.85 | −0.92 / −0.94 / −0.91 | −0.87 / −0.87 / −0.84 | −0.93 / −0.94 / −0.91 |
| Fe VII | −1.05 / −1.06 / −1.03 | −1.13 / −1.14 / −1.11 | −1.04 / −1.04 / −1.01 | −1.14 / −1.14 / −1.10 |
| Fe VIII | −0.99 / −0.99 / −0.97 | −1.06 / −1.08 / −1.05 | −0.98 / −0.98 / −0.95 | −1.08 / −1.07 / −1.04 |

\* involves one of the six ground-character-changed runs (§1.4).

Trends **[Verified]**:
- For Fe IV–VIII both effects are nearly μ-independent. The orbital effect is −0.2 to −0.5 eV; the CI effect is −0.7 to −1.1 eV and grows in magnitude with charge. SD CI lowers these IPs, moving them *away* from the paper and NIST.
- For Fe II/III, the averaged-orbital CI effect is strongly μ-dependent (up to +3.6 eV), while the ground-orbital CI effect is small and constant. The μ-dependence coincides with the ground-character changes.

CSF counts (`D/results/csf_counts.csv`) **[Verified]**: the reconstruction of the stated active space gives 1.23–1.42 × the published N_CSF (all J, both parities). FAC NLEV equals the counted ground-block CSFs exactly for every `sdp` run (e.g. Fe III 3363, Fe VI 10591). What the paper restricted is **[Unresolved]**; only the ground (P, 2J) block is used here. This is a CI-space difference. Given the size of the CI effects above (≤1.1 eV for Fe IV–VIII, ≤0.33 eV for ground-optimized Fe III), it is unlikely to explain the Fe III anomaly **[Supported interpretation]**.

## 5. Baseline IP Error Versus Screening-Induced Change

ΔIP(μ) = IP(μ) − IP(0) in eV; the published ΔIP and the first-order Debye–Hückel limit −(z+1)μ·27.2114 eV (z+1 = charge of the product ion) are shown for reference. Each entry is ΔIP (error in eV, relative error in |ΔIP| in %) **[Verified]**:

| ion | μ | IP(0) / IP(μ) gsdp | paper ΔIP | 1st order | ref | sdp | gref | gsdp |
|---|---|---|---:|---:|---|---|---|---|
| Fe II | 0.1 | 15.505 / 10.197 | −4.690 | −5.442 | −4.684 (+0.01, −0.1%) | −5.590* (−0.90, +19.2%) | −5.295 (−0.60, +12.9%) | −5.308 (−0.62, +13.2%) |
| Fe II | 0.25 | 15.505 / 4.658 | −11.110 | −13.606 | −10.154 (+0.96, −8.6%) | −11.840* (−0.73, +6.6%) | −10.833 (+0.28, −2.5%) | −10.847 (+0.26, −2.4%) |
| Fe III | 0.1 | 27.766 / 19.802 | −7.330 | −8.163 | −8.495 (−1.16, +15.9%) | −7.641* (−0.31, +4.2%) | −7.960 (−0.63, +8.6%) | −7.965 (−0.63, +8.7%) |
| Fe III | 0.25 | 27.766 / 8.960 | −16.510 | −20.409 | −19.626 (−3.12, +18.9%) | −16.064* (+0.45, −2.7%) | −18.824 (−2.31, +14.0%) | −18.806 (−2.30, +13.9%) |
| Fe IV | 0.1 | 52.587 / 42.161 | −10.720 | −10.885 | −10.429 (+0.29, −2.7%) | −10.454 (+0.27, −2.5%) | −10.413 (+0.31, −2.9%) | −10.426 (+0.29, −2.7%) |
| Fe IV | 0.25 | 52.587 / 27.439 | −25.730 | −27.211 | −25.213 (+0.52, −2.0%) | −25.209 (+0.52, −2.0%) | −25.169 (+0.56, −2.2%) | −25.148 (+0.58, −2.3%) |
| Fe V | 0.1 | 73.212 / 60.098 | −13.360 | −13.606 | −13.089 (+0.27, −2.0%) | −13.110 (+0.25, −1.9%) | −13.111 (+0.25, −1.9%) | −13.114 (+0.25, −1.8%) |
| Fe V | 0.25 | 73.212 / 41.760 | −32.090 | −34.014 | −31.471 (+0.62, −1.9%) | −31.468 (+0.62, −1.9%) | −31.474 (+0.62, −1.9%) | −31.452 (+0.64, −2.0%) |
| Fe VI | 0.1 | 97.061 / 81.295 | −16.000 | −16.327 | −15.771 (+0.23, −1.4%) | −15.788 (+0.21, −1.3%) | −15.764 (+0.24, −1.5%) | −15.766 (+0.23, −1.5%) |
| Fe VI | 0.25 | 97.061 / 59.290 | −38.410 | −40.817 | −37.824 (+0.59, −1.5%) | −37.812 (+0.60, −1.6%) | −37.799 (+0.61, −1.6%) | −37.770 (+0.64, −1.7%) |
| Fe VII | 0.1 | 122.848 / 104.430 | −18.640 | −19.048 | −18.435 (+0.21, −1.1%) | −18.450 (+0.19, −1.0%) | −18.416 (+0.22, −1.2%) | −18.418 (+0.22, −1.2%) |
| Fe VII | 0.25 | 122.848 / 78.739 | −44.740 | −47.620 | −44.174 (+0.57, −1.3%) | −44.156 (+0.58, −1.3%) | −44.142 (+0.60, −1.3%) | −44.109 (+0.63, −1.4%) |
| Fe VIII | 0.1 | 148.858 / 127.758 | −21.290 | −21.769 | −21.096 (+0.19, −0.9%) | −21.110 (+0.18, −0.8%) | −21.103 (+0.19, −0.9%) | −21.099 (+0.19, −0.9%) |
| Fe VIII | 0.25 | 148.858 / 98.391 | −51.160 | −54.423 | −50.505 (+0.65, −1.3%) | −50.489 (+0.67, −1.3%) | −50.501 (+0.66, −1.3%) | −50.467 (+0.69, −1.4%) |

(`sdv`/`gsdv` agree with `sdp`/`gsdp` to ≤0.03 eV for Fe IV–VIII; the full set is in the session output summarised here.)

Two-point fit ΔIP = aμ + bμ², ratio a / [−(z+1)·27.2114] **[Verified arithmetic]**:

| ion | paper | ref | gref | gsdp |
|---|---:|---:|---:|---:|
| Fe II | 0.892 | 0.937 | 1.091 | 1.094 |
| Fe III | 0.957 | 1.093 | 1.010 | 1.012 |
| Fe IV | 1.011 | 0.979 | 0.978 | 0.980 |
| Fe V | 1.008 | 0.987 | 0.989 | 0.990 |
| Fe VI | 1.006 | 0.992 | 0.992 | 0.992 |
| Fe VII | 1.005 | 0.995 | 0.993 | 0.994 |
| Fe VIII | 1.003 | 0.996 | 0.997 | 0.997 |

The quadratic coefficients b agree between FAC and the paper to within about 20% for Fe IV–VIII (e.g. Fe VI 42.8 vs 42.4).

Answers:
1. *Is the 1–3% mismatch systematic?* **[Verified]** Yes for Fe IV–VIII: same sign at both μ, magnitude falling with charge (Fe IV 2.0–2.9% → Fe VIII 0.8–1.4%).
2. *Does it depend on orbital optimization?* **[Verified]** No, not for Fe IV–VIII (≤0.4 percentage points).
3. *Does it depend on configuration space?* **[Verified]** No, not for Fe IV–VIII (≤0.1 percentage points).
4. *Can source evidence explain it?* **[Unresolved]** The deficit sits mainly in the linear term. For exactly stationary wave functions the linear term of E(μ) is μ[ZN − N(N−1)/2], so dIP/dμ|₀ = −(z+1) for any such method. FAC reproduces this to 0.02% for He/Li-like Fe and 0.06% for Ne-like Fe (`D/diagnostics/first_order`), so the Yukawa integrals themselves are implemented correctly. FAC orbitals are not stationary for the ground-level CI energy, so orbital response can enter at first order. Two features of the screened SCF potential that only affect the orbitals are verified in the source: the local exchange term in `DensityToSZ` (`zx = ahx*XCPotential(tps,nx,vxm)*rx`) carries no μ-screening, and the asymptotic self-interaction tail logic in `PotentialHX` (`_csi`, `ihx`) was designed for a Coulomb tail. Their quantitative effect is unknown.
5. *Method difference or open hypothesis?* **[Unresolved]** The paper's own fitted slopes deviate from 1 too (1.003–1.011, and 0.89/0.96 for Fe II/III), so neither code satisfies the stationary limit exactly. Calling the 1–3% a "method difference" is a hypothesis.
6. **Small-μ direct slopes** (`D/diagnostics/small_mu`, μ = 10⁻³) range from 0.69 to 0.98 for gref Fe III–VIII and are erratic for ref Fe II/III (−0.45, 2.69). Warm-restart history noise is ~0.01 eV against shifts of only 0.08–0.22 eV, so these slopes are too noisy to decide anything. **[Verified noise level]**

A model can reproduce the screening shift with a badly wrong baseline: `ref` Fe IV has a −2.49 eV baseline error but its ΔIP is within 0.29 / 0.52 eV of the paper. **[Verified]** For Fe II/III this does not hold, because the ground-state character itself depends on μ and on the model.

## 6. Screening Implementation and Numerical Protocol

Statements from the brief, checked against the evidence:

| statement | finding |
|---|---|
| μ > 0 inputs use `debye_mode 1`, `ee_screen 1`, `hxs_iter0 8` | **[Verified]** every μ > 0 `in.sf` (`A/raw.json`) |
| ladder of increasing μ | **[Verified]** 2 steps for μ = 0.1 (0.05, 0.1), 5 for μ = 0.25 (0.05…0.25), after a cold free-ion optimization |
| `nps = 1e-6` is a density, not μ | **[Verified]** last `PlasmaScreen(zps, 1e-06, T, 1, 0.0, 1)`; e.g. Fe III μ = 0.25: T = 8.107×10⁻⁴ eV, λ = 4.000 a₀ |
| smallest μ run = 10⁻⁴ | **[Verified]** `D/diagnostics/small_mu`, `D/diagnostics/first_order` |
| no run used `debye_mode 0` | **[Verified]** production inputs, and `D/diagnose_radial_failure.py` passes `en` = 1 in every case |
| "debye_mode 1 derives λ from μ" | **Correction [Verified, SRC orbital.c InitializePS]:** FAC computes λ = [T/(4πn(1+ups))]^{1/2} from T and n. The driver chooses T to give λ = 1/μ, and `final.pot` gives dps·μ = 1 in all μ > 0 runs. zps and nps insensitivity was tested for Fe II ground at μ = 0.1 only (diagnostics E4, E5). |
| "the screened SCF patch changes the DensityToSZ truncation" | **Correction [Verified]:** none of my patches touches `DensityToSZ`; its `ips` truncation is still in the source. Production avoids it by protocol: nps = 10⁻⁶, and the grid is built in the free step, so `ips = 0` in **every** μ > 0 `final.pot`. The Yukawa direct screening in `DensityToSZ` itself comes from branch commit 4dd4b40. |
| int-overflow fix in the production binary | **[Verified]** `mmalloc.c`/`structure.c` hunks are in diff `ceac162e…` |
| debye_mode 0 sign | **[Verified]** pre-#840 sign in this binary; not exercised |

Could a known implementation detail explain the observations?
- **Fe III baseline anomaly:** no. The same screening code is used in every model, and the anomaly exists at μ = 0, where `PlasmaScreen` is not called at all. **[Verified]**
- **Fe IV–VIII baseline errors (−1 to −3 eV at μ = 0):** these are free-ion errors, so screening code is irrelevant to them. **[Verified]**
- **Residual ΔIP mismatch:** possibly, via the orbital-only features listed in §5.4. **[Untested hypothesis]**

Solver and convergence defects (RadialBound bracket, cold-start transient) were fixed or bypassed before production. All 144 runs converged; IP sensitivity to grid is ≤0.023 eV and to warm restarts ≤0.011 eV (earlier diagnostics). These numerical effects are 2–3 orders of magnitude smaller than the Fe III discrepancy. **[Verified]**

## 7. Competing Explanations and Their Evidence

| explanation | for | against | status |
|---|---|---|---|
| E0 extraction, units or pairing error | — | independent audit; E0 = minimum level; NELE; same binary | rejected **[Verified]** |
| screening settings differ between adjacent ions | — | identical options, λ, zps rule | rejected **[Verified]** |
| CI space (SD truncation, CSF counts) causes the Fe III anomaly | — | CI-fixed orbital swap gives +8.0 to +8.4 eV; ground-orbital CI effect ≤0.33 eV | rejected as the main cause **[Verified]** |
| orbital average weights (jj-equal, small ground weight, ⟨4s⟩ = 1) | penalty 13.3 vs 4.9 eV; n = 3 eigenvalue shifts 4× larger; 4s-bearing levels over-stabilized | not isolated from the ground-weight effect | **[Supported interpretation]**; §8 tests it |
| missing differential correlation (FAC −2.9 eV vs NIST) | gsdp below NIST for all ions; no 4p/4d/4f correlation | no FAC calculation with a larger space | **[Supported interpretation]** |
| paper's own model error (+2.3 eV vs NIST) | NIST − paper = −2.28 eV (Fe III), −2.69 (Fe II) | unknown paper details | **[Verified arithmetic]**; cause **[Unresolved]** |
| block-restricted E0 is not the true ground (6 runs) | ground character changes; ⁷S₃ (2J = 6) not computed | not computed | defect **[Verified]**; size **[Untested]** |
| implementation error in the screened SCF affecting ΔIP | ΔIP linear term 1–3% low | He/Li/Ne first-order exact; model-independent | **[Unresolved]** |

## 8. Proposed Single Controlled Test (not executed)

**Question:** is the excess Fe III orbital penalty (13.34 eV vs 4.90 eV for Fe IV) caused by the six 4s-bearing references in the orbital average?

1. **Start from** `D/runs/ref/FeIII/mu0/in.sf` (production `ref` protocol, binary `3a77793e`).
2. **Change only the optimization group, not the CI set.** Split the existing `Config('ref', …)` into two groups:
   - `Config('o', '1s2 2s2 2p6 3s2 3p6 3d6', '1s2 2s2 2p6 3s2 3p5 3d7', '1s2 2s2 2p5 3s2 3p6 3d7', '1s1 2s2 2p6 3s2 3p6 3d7')`. These are the four references without 4s; they follow the same pattern as Fe IV–VIII.
   - `Config('x', '1s2 2s2 2p6 3s2 3p5 3d6 4s1', '1s2 2s2 2p6 3s2 3p5 3d5 4s2', '1s2 2s2 2p5 3s2 3p6 3d6 4s1', '1s2 2s2 2p5 3s2 3p6 3d5 4s2', '1s1 2s2 2p6 3s2 3p6 3d6 4s1', '1s1 2s2 2p6 3s2 3p6 3d5 4s2')`. These are the six 4s-bearing references, with strings taken from the existing `in.sf`.
   - `OptimizeRadial(['o'])`; `Structure('e.lev.b', ['o','x'])`. The CI set is then identical to `ref`.
3. **Repeat Fe III only.** Fe IV's `ref` references contain no 4s, so the existing `ref` Fe IV E0 (−34522.7607 eV at μ = 0) is the correct partner. Run μ = 0 first; μ = 0.1 and 0.25 with the identical ladder are optional extensions.
4. **Keep fixed:** binary, grid policy (cold free-ion optimization first), all-J CI, no `ConfigEnergy`, and for μ > 0 the ladder, `debye_mode 1`, `ee_screen 1`, `hxs_iter0 8`, nps = 10⁻⁶.
5. **Compare:**
   - P′ = E0(test) − E0(gref Fe III) = E0(test) + 34555.7544 eV, against P_ref = 13.342 eV and the Fe IV–V pattern (4.90, 4.45 eV).
   - IP′ = −34522.7607 − E0(test), against ref 19.652 and gref 28.094 eV.
   - The 3d₅/₂ eigenvalue shift (`final.pot`), against −0.254 Ha (ref) and ≈ −0.06 Ha (Fe IV).
   - The average occupancies in `final.pot`, which should be ⟨4s⟩ = 0 and ground weight 5/25 = 20% under the jj-equal model.
6. **If the 4s hypothesis is correct:** P′ ≈ 5–6 eV, IP′ ≈ 26–28 eV, and an eigenvalue shift ≲ 0.1 Ha.
7. **If the cause is elsewhere** (small ground weight or core-hole fraction in general): P′ stays ≳ 11 eV, IP′ ≲ 22 eV. If CI were the cause, P′ would change little, but CI is already excluded by the fixed-CI comparison.
8. **Confounders:**
   - Removing six references also changes the ground weight (6.7% → 20%) and the core-hole fractions. The test isolates "4s-bearing references" as a group, not 4s occupancy alone.
   - A second, optional run giving the `x` group a small nonzero weight (`OptimizeRadial(['o','x'],[w,1-w])`, `acfg_wmode` aware) would separate the effects but changes two things at once.
   - SCF history dependence is ≈0.01 eV, which is negligible.
9. **Decision criterion:**
   - Support: P′ ≤ 6.5 eV (≥ 80% of the 8.44 eV excess over Fe IV removed).
   - Weaken: P′ ≥ 11 eV (≤ 28% removed).
   - In between: mixed cause; follow up with the weight-scan variant.

**Is it sufficient?** It decides the most important open question, the cause of the Fe III anomaly. It does **not** address the ~5.2 eV gap between ground-optimized FAC and the paper, the block-restriction defect in six runs (that needs an all-J recomputation of those `sdv`/`sdp` cases), or the 1–3% ΔIP mismatch. A larger sweep is not justified until this test has run.

## 9. Conclusions and Unresolved Questions

- **[Verified]** The Fe III ≈19.7 eV IP in `ref`, `sdv` and `sdp` is an orbital-average artefact. Ground-optimized orbitals with the same CI give 27.8–28.1 eV.
- **[Verified]** FAC weights jj configurations equally within a group. For Fe III this means ground weight 6.7% and ⟨4s⟩ = 1.0 in the average. Earlier "statistical-weight" descriptions were wrong.
- **[Verified]** Six `sdv`/`sdp` runs (Fe II at μ = 0.25, Fe III at μ = 0.1 and 0.25) have a ground-character change inside a symmetry-restricted block; their E0 and IPs are not reliable ground-state values.
- **[Verified arithmetic / Supported interpretation]** The remaining gsdp − paper gap of −5.16 eV is half FAC below experiment (−2.89 eV) and half the paper above experiment (+2.28 eV).
- **[Verified]** For Fe IV–VIII, ΔIP is model-independent and 0.8–2.9% smaller in magnitude than published, mostly in the linear coefficient. **[Unresolved]** cause.
- **[Unresolved]** The paper's CSF restrictions, its Fe III ground-state composition, and whether its EAL also over-stabilizes 4s configurations.

**Most important unresolved scientific question:** is the excess Fe III orbital-averaging penalty (13.3 eV vs ~4.9 eV for Fe IV), which drives the Fe III → Fe IV anomaly, caused by the 4s-bearing reference configurations in FAC's equal-jj-weight average?
**Follow-up most likely to answer it:** the single Fe III run of §8. Optimize on the four non-4s references while keeping the full 10-reference CI set, and compare P′ with 13.34 eV (ref) and ~5 eV (the Fe IV–V pattern).
