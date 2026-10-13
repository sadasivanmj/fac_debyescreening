# Fe II → Fe VIII ionization potentials at μ = 0 (`ref` model): validation against Deprince et al. (2020)

Date: 2026-10-12 (UTC). Read-only audit of existing raw outputs. **No FAC/SFAC calculation was launched and no project file was modified**; no `sfac` process was running at the time of the audit. Run directory root: `R = /data/project/fac/fe_calc/deprince_t2/runs/ref/<ion>/mu0/`.

## 1. Summary

| Transition | E0 lower (eV) | E0 upper (eV) | IP FAC (eV) | Deprince | ΔIP (eV) | ε_rel | NIST | FAC − NIST (eV) | verification |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Fe II → Fe III | −34561.2842 | −34542.4125 | **18.8717** | 18.89 | −0.0183 | 0.10 % | 16.19920 | +2.6725 | all checks pass |
| Fe III → Fe IV | −34542.4125 | −34522.7607 | **19.6518** | 32.93 | **−13.2782** | **40.32 %** | 30.651 | −10.9992 | all checks pass |
| Fe IV → Fe V | −34522.7607 | −34469.8762 | **52.8845** | 55.37 | −2.4855 | 4.49 % | 54.91 | −2.0255 | all checks pass |
| Fe V → Fe VI | −34469.8762 | −34396.1570 | **73.7192** | 74.92 | −1.2008 | 1.60 % | 75.0 | −1.2808 | all checks pass |
| Fe VI → Fe VII | −34396.1570 | −34298.5258 | **97.6312** | 98.20 | −0.5688 | 0.58 % | 98.985 | −1.3538 | all checks pass |
| Fe VII → Fe VIII | −34298.5258 | −34174.8058 | **123.7200** | 123.60 | +0.1200 | 0.10 % | 124.98 | −1.2600 | all checks pass |

ΔIP = IP_FAC − IP_Deprince and ε_rel = 100·|ΔIP|/IP_Deprince. Arithmetic: IP = E0(upper) − E0(lower). For example, Fe III → Fe IV: −34522.7607 − (−34542.4125) = 19.6518 eV. E0 is printed to 1×10⁻⁴ eV, so each IP carries ±0.0001 eV rounding.

Deprince values are **calculated** (MCDHF/EAL). NIST values are **experimental/tabulated**. Agreement with one does not imply agreement with the other: Fe II agrees with Deprince to 0.10 % but lies 2.67 eV above NIST.

**Reproduced within 0.6 %:** Fe II, Fe VI, Fe VII. **Moderate deviation:** Fe V (1.6 %), Fe IV (4.5 %). **Large discrepancy:** Fe III → Fe IV (−13.28 eV, 40 %).

## 2. Configuration and electron-count verification [Verified]

Each `in.sf` contains a single group `ref`. Its strings were compared one by one, in order, with the Table 1 configurations expanded on a closed `1s2 2s2 2p6 3s2 3p6` core (a hole [nl] removes one electron from that subshell). All match, with no extra or missing entries.

| Ion | `in.sf` configurations | count | e⁻ (all) |
|---|---|---:|---:|
| Fe II | `…3p6 3d6 4s1`, `…3p5 3d7 4s1`, `…3p5 3d6 4s2`, `2p5 … 3d7 4s1`, `2p5 … 3d6 4s2`, `1s1 … 3d7 4s1`, `1s1 … 3d6 4s2` | 7 | 25 |
| Fe III | `3d6`, `3p5 3d7`, `3p5 3d6 4s1`, `3p5 3d5 4s2`, `2p5 3d7`, `2p5 3d6 4s1`, `2p5 3d5 4s2`, `1s1 3d7`, `1s1 3d6 4s1`, `1s1 3d5 4s2` | 10 | 24 |
| Fe IV | `3d5`, `3p5 3d6`, `2p5 3d6`, `1s1 3d6` | 4 | 23 |
| Fe V | `3d4`, `3p5 3d5`, `2p5 3d5`, `1s1 3d5` | 4 | 22 |
| Fe VI | `3d3`, `3p5 3d4`, `2p5 3d4`, `1s1 3d4` | 4 | 21 |
| Fe VII | `3d2`, `3p5 3d3`, `2p5 3d3`, `1s1 3d3` | 4 | 20 |
| Fe VIII | `3d1`, `3p5 3d2`, `2p5 3d2`, `1s1 3d2` | 4 | 19 |

The full strings are in `R/in.sf`; unlisted subshells are closed.

## 3. Protocol, executable and run status [Verified]

All seven `in.sf` files are identical apart from the configuration list:
```
SetAtom('Fe')
Config('ref', <Table 1 references>)
OptimizeRadial(['ref'])
GetPotential('final.pot')
Structure('e.lev.b', ['ref'])
MemENTable('e.lev.b')
PrintTable('e.lev.b', 'e.lev', 1)
```
- **No screening:** no `PlasmaScreen`, `SetOption` or `ConfigEnergy`. Every `final.pot` has `mps = -1` and zps = nps = dps = 0, so there is no residual screening.
- **Same executable for all seven:** `run.json` records `/data/project/fac_worktrees/scf/inst/bin/sfac`, sha256 `3a77793e4b5e6f2c…`. This equals the binary on disk, re-hashed now. Commit `ff9a5c9`, source diff `ceac162e…`. Every adjacent pair uses the same executable.
- **Clean completion:** `returncode 0`, `status ok`, `sfac.log` empty (0 bytes; no error or convergence message) for all seven ions.

## 4. Ground-state identification and units [Verified]

| Ion | `e.lev` E0 header | NELE | NLEV | lowest level | parity | 2J | label |
|---|---|---:|---:|---|---:|---:|---|
| Fe II | `0, -3.45612842E+04` | 25 | 986 | index 0, rel. 0 | even | 9 | 3d6.4s1 |
| Fe III | `0, -3.45424125E+04` | 24 | 1640 | index 0, rel. 0 | even | 8 | 3d6 |
| Fe IV | `0, -3.45227607E+04` | 23 | 460 | index 0, rel. 0 | even | 5 | 3d5 |
| Fe V | `0, -3.44698762E+04` | 22 | 536 | index 0, rel. 0 | even | 0 | 3d4 |
| Fe VI | `0, -3.43961570E+04` | 21 | 442 | index 0, rel. 0 | even | 3 | 3d3 |
| Fe VII | `0, -3.42985258E+04` | 20 | 267 | index 0, rel. 0 | even | 4 | 3d2 |
| Fe VIII | `0, -3.41748058E+04` | 19 | 108 | index 0, rel. 0 | even | 3 | 3d1 |

- The E0 header is the absolute energy of the minimum-energy level over **all computed levels**: every J and both parities of the reference-configuration CI. Its units are eV (`PrintTable` verbose writes energy × HARTREE_EV; `final.pot` EUNIT = 27.2114).
- The 2J values match the expected ground terms (a⁶D₉/₂, ⁵D₄, ⁶S₅/₂, ⁵D₀, ⁴F₃/₂, ³F₂, ²D₃/₂).
- **Limit of the claim:** this is the lowest level *within the reference configuration space*. Configurations outside Table 1 (e.g. Fe II 3d⁷ or Fe III 3d⁵4s) are not included, so this is not a claim about the true ground state of the ion.

## 5. Fe III → Fe IV discrepancy

**Present in the raw energies [Verified].** −34522.7607 − (−34542.4125) = 19.6518 eV, against Deprince's 32.93 eV. Priorities 1–3 are excluded by §2–§4: the configuration strings, electron counts, ground-state selection, units, protocol and executable are all correct and identical.

**Orbital optimization produces it [Verified from existing raw outputs].** The existing `gref` runs (`…/runs/gref/<ion>/mu0/`) use **the same CI configuration set** and the same binary, but optimize the orbitals on the ground configuration only (`OptimizeRadial(['g0'])`). The orbital penalty P = E0(ref) − E0(gref) is:

| Ion | Fe II | Fe III | Fe IV | Fe V | Fe VI | Fe VII | Fe VIII |
|---|---:|---:|---:|---:|---:|---:|---:|
| E0 gref (eV) | −34571.2747 | −34555.7544 | −34527.6602 | −34474.3260 | −34400.3706 | −34302.3757 | −34178.3917 |
| P (eV) | 9.9905 | **13.3419** | 4.8995 | 4.4498 | 4.2136 | 3.8499 | 3.5859 |

IP(ref) − IP(gref) = P(upper) − P(lower). For Fe III → Fe IV this is 4.8995 − 13.3419 = **−8.4424 eV**: with only the orbital-optimization group changed, the IP moves from 19.65 to 28.09 eV. The same arithmetic gives Fe II → Fe III **+3.3514 eV** (P(Fe III) > P(Fe II)). The ref model's close agreement for Fe II is therefore partly the result of this penalty difference; gref gives 15.52 eV. For Fe IV–Fe VIII the P differences are only −0.20 to −0.45 eV.

**Why Fe III is penalized most: supported interpretation, not proven.**
- FAC weights each jj configuration of a group equally (`fac_worktrees/scf/faclib/config.c`, `GetAverageConfig`). This rule reproduces the average occupancies in every `final.pot` (see `ip_discrepancy_analysis.md`).
- For Fe III this gives a ground weight of 5/75 = 6.7 % and an average 4s occupancy of 1.00, while the ground 3d⁶ has no 4s. Fe IV–VIII averages contain no 4s.
- That the six 4s-bearing Fe III references cause the excess penalty is an **untested hypothesis**. The controlled one-run test is specified in `ip_discrepancy_analysis.md` §8 and has not been run.

**What remains after orbital optimization [Verified arithmetic; cause not established].** The gref Fe III IP, 28.09 eV, is still 4.84 eV below Deprince and 2.56 eV below NIST. Deprince's own value is 2.28 eV above NIST. Attributing FAC's remaining deficit to missing correlation or other approximations is a possible explanation only.

Other transitions: the −1.2 to −2.5 eV deviations for Fe IV–VI are free-ion model differences from Deprince's MCDHF/EAL calculation with SD active space. They are not traced to any input or protocol defect here, and their cause is unresolved. Every FAC `ref` IP for Fe IV–VIII lies 1.3–2.0 eV below NIST, apart from Fe IV at −2.03 eV.

## 6. Verified facts vs possible explanations

- **Verified:**
  - All seven inputs match Table 1, with correct electron counts.
  - The protocol and executable are identical across all seven runs.
  - μ = 0 with no screening; all runs completed cleanly.
  - E0 values are the lowest computed level, in eV.
  - The six IPs above.
  - The Fe III → Fe IV discrepancy is in the raw energies and is reduced by 8.44 eV by changing only the orbital-optimization group.
  - FAC's average uses equal jj-configuration weights.
- **Supported interpretation:** the Fe III average, dominated by core-hole and 4s-bearing configurations, gives orbitals poorly suited to the 3d⁶ ground state.
- **Untested hypothesis:** the 4s-bearing references are the specific cause.
- **Unresolved:** the remaining gaps to Deprince and to NIST.

## 7. Integrity

No calculation was launched and no source, input, driver, binary or result file was read for writing. The only file created is this report, written with no-clobber. Fe VIII → Fe IX and the screened cases were not examined, as instructed.
