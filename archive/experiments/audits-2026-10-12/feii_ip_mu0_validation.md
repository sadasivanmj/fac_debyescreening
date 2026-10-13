# Fe II ionization potential at μ = 0: validation against Deprince et al. Table 1

Date: 2026-10-12 (UTC). Read-only validation. No FAC/SFAC calculation was run, because the existing inputs match the specification exactly (§2). No existing file was modified.

## 1. Result

| quantity | value |
|---|---|
| E0(Fe III), raw FAC output | −34542.4125 eV |
| E0(Fe II), raw FAC output | −34561.2842 eV |
| IP = E0(Fe III) − E0(Fe II) | −34542.4125 − (−34561.2842) = **18.8717 eV** |
| Deprince et al. (calculated, μ = 0) | 18.89 eV |
| FAC − Deprince | −0.0183 eV (**0.097 %** absolute) |
| NIST experimental IP (as quoted in Deprince Table 2) | 16.19920 eV |
| FAC − NIST | +2.6725 eV (+16.50 %) |

The two E0 values are printed in `e.lev` to 9 significant figures (1×10⁻⁴ eV resolution), so the IP carries a ±0.0001 eV rounding uncertainty.

18.8717 eV is a **FAC-calculated** value: the difference between two calculated total energies, with orbitals optimized on the average of the Table 1 reference configurations. 18.89 eV is Deprince et al.'s **calculated** (MCDHF/EAL) value. The NIST 16.19920 eV is the **experimental/tabulated** ionization energy. The close FAC–Deprince agreement is between two models; both lie about 2.7 eV above experiment.

## 2. Configuration verification

Hole notation [nl] was expanded by removing one electron from the filled subshell (1s2 2s2 2p6 3s2 3p6 core), with 3d/4s occupations as given. Every configuration was compared string-for-string with `in.sf` (all MATCH, same order, no extra or missing entries):

**Fe II** (`/data/project/fac/fe_calc/deprince_t2/runs/ref/FeII/mu0/in.sf`, group `ref`)

| # | Table 1 | FAC string | e⁻ |
|---|---|---|---|
| 1 | 3d⁶4s | `1s2 2s2 2p6 3s2 3p6 3d6 4s1` | 25 |
| 2 | [3p]3d⁷4s | `1s2 2s2 2p6 3s2 3p5 3d7 4s1` | 25 |
| 3 | [3p]3d⁶4s² | `1s2 2s2 2p6 3s2 3p5 3d6 4s2` | 25 |
| 4 | [2p]3d⁷4s | `1s2 2s2 2p5 3s2 3p6 3d7 4s1` | 25 |
| 5 | [2p]3d⁶4s² | `1s2 2s2 2p5 3s2 3p6 3d6 4s2` | 25 |
| 6 | [1s]3d⁷4s | `1s1 2s2 2p6 3s2 3p6 3d7 4s1` | 25 |
| 7 | [1s]3d⁶4s² | `1s1 2s2 2p6 3s2 3p6 3d6 4s2` | 25 |

**Fe III** (`/data/project/fac/fe_calc/deprince_t2/runs/ref/FeIII/mu0/in.sf`, group `ref`)

| # | Table 1 | FAC string | e⁻ |
|---|---|---|---|
| 1 | 3d⁶ | `1s2 2s2 2p6 3s2 3p6 3d6` | 24 |
| 2 | [3p]3d⁷ | `1s2 2s2 2p6 3s2 3p5 3d7` | 24 |
| 3 | [3p]3d⁶4s | `1s2 2s2 2p6 3s2 3p5 3d6 4s1` | 24 |
| 4 | [3p]3d⁵4s² | `1s2 2s2 2p6 3s2 3p5 3d5 4s2` | 24 |
| 5 | [2p]3d⁷ | `1s2 2s2 2p5 3s2 3p6 3d7` | 24 |
| 6 | [2p]3d⁶4s | `1s2 2s2 2p5 3s2 3p6 3d6 4s1` | 24 |
| 7 | [2p]3d⁵4s² | `1s2 2s2 2p5 3s2 3p6 3d5 4s2` | 24 |
| 8 | [1s]3d⁷ | `1s1 2s2 2p6 3s2 3p6 3d7` | 24 |
| 9 | [1s]3d⁶4s | `1s1 2s2 2p6 3s2 3p6 3d6 4s1` | 24 |
| 10 | [1s]3d⁵4s² | `1s1 2s2 2p6 3s2 3p6 3d5 4s2` | 24 |

## 3. Protocol (identical for both ions)

Both `in.sf` files contain exactly:
```
SetAtom('Fe')
Config('ref', <7 or 10 configurations above>)
OptimizeRadial(['ref'])
GetPotential('final.pot')
Structure('e.lev.b', ['ref'])
MemENTable('e.lev.b')
PrintTable('e.lev.b', 'e.lev', 1)
```
- μ = 0: no `PlasmaScreen` and no `SetOption`. `final.pot` shows `mps = -1` and zps = nps = dps = 0 for both ions.
- No `ConfigEnergy` call. The original preliminary baseline (`fe_calc/fe_ip_sfac.py`, with `ConfigEnergy(0)`/`(1)`) gave −34561.2841 eV for Fe II and 18.8716 eV for the IP. That 0.0001 eV difference is the `ConfigEnergy` step (diagnostic `fe_calc/deprince_t2/diagnostics/E3_free_installed_CE`).
- Orbitals: optimized on the FAC average of the `ref` group, in which FAC weights each jj configuration equally. CI: all J and both parities of the reference configurations (NLEV 986 for Fe II, 1640 for Fe III).
- `sfac.log` is empty for both runs (0 bytes, no warnings). `run.json`: `status ok`, `returncode 0`.

## 4. Ground state and units

| | Fe II | Fe III |
|---|---|---|
| `e.lev` header | `E0 = 0, -3.45612842E+04`, NELE 25, NLEV 986 | `E0 = 0, -3.45424125E+04`, NELE 24, NLEV 1640 |
| ground level (index 0) | relative energy 0, parity 0 (even), 2J = 9, `3d6.4s1` | relative energy 0, parity 0 (even), 2J = 8, `3d6` |
| check | minimum over all 986 levels = index 0 | minimum over all 1640 levels = index 0 |

The ground levels are a⁶D₉/₂ (Fe II) and ⁵D₄ (Fe III), the expected terms. E0 is the absolute energy of the lowest level, in **eV**: FAC's `PrintTable` (verbose) writes `energy × HARTREE_EV`, and `final.pot` gives `EUNIT = 27.2113860`. No further unit conversion is applied.

## 5. Binary and provenance

Both `run.json` files record:
- `sfac`: `/data/project/fac_worktrees/scf/inst/bin/sfac`
- sha256 `3a77793e4b5e6f2c42efe2c358cfca09c911fd629895ccde48c15eac1894e229`, equal to the binary currently on disk (re-hashed during this validation)
- source commit `ff9a5c9cc6623e3f6adb8062ba23a52e9aace582` (ee-screening-scf) plus the working-tree diff sha256 `ceac162ead6ee5d5159e8ed02d1a4577791be1cd987e1ce2585fa49a0267d7a0`

The patches in that diff (Debye-only RadialBound bracket, opt-in `radial:hxs_iter0`, int-overflow fix) are inert for these μ = 0 runs. The patched and unpatched builds give identical free-ion Fe II E0 (diagnostics E1/E2).

Output files used: `.../runs/ref/{FeII,FeIII}/mu0/{in.sf,e.lev,final.pot,sfac.log,run.json}` under `/data/project/fac/fe_calc/deprince_t2/`.

## 6. Discrepancies

None between the existing runs and the Table 1 specification. The only differences found are:
- the 0.0001 eV `ConfigEnergy` effect relative to the original preliminary baseline (§3);
- the 0.0183 eV (0.097 %) difference from the published 18.89 eV;
- the +2.67 eV difference of both calculations from the experimental NIST value.

This agreement is specific to Fe II. The same `ref` protocol gives 19.65 eV for Fe III → Fe IV, against 32.93 eV published; see `ip_discrepancy_analysis.md`. The 18.87 eV agreement therefore does not validate the `ref` orbital-averaging protocol in general.
