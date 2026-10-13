# Analysis from branch `ee-screening-scf`

This folder holds the research files of branch `ee-screening-scf`, brought in
by a git merge into `integration/debye-screening-complete`:

- **4dd4b40** "Implement self-consistent e-e Debye screening in DensityToSZ
  and restructure analysis scripts"
- **ff9a5c9** "analysis files", 2026-10-09

## What it is

That branch had its **own** screened SCF direct potential. It multiplies the
Coulomb integrands by `_shx(μr)` and `exp(−μr)` inside `DensityToSZ`. This
work was done with that implementation:

| folder | contents |
|---|---|
| `debye_analysis/scripts/` | analysis scripts: C4+ ground state (several variants), Chen et al. 2018 benchmark and Kr³⁴⁺ extension, Xie et al. 2012, Fritzsche (2006) runs and plots, Yukawa-potential and e–e-isolation checks |
| `debye_analysis/outputs/` | their outputs and raw runs |
| `xie2012_run/`, `xie2012_muc_run/` | `sfac` inputs and outputs of the Xie 2012 reproduction |
| `fritzsche_2006_results/` | Fritzsche (2006) comparison results |
| `chen2018_kr34_run/` | Chen et al. 2018 Kr³⁴⁺ runs |
| `c4_ground_results/` | C4+ ground-state results of that branch |
| `root-files/` | figures and FAC failure dumps that were at the repository root; `ver` renamed to `quick_verify_stage2.py` (its own docstring name) |

## Status

- **Superseded implementation.** That screened `DensityToSZ` has been
  replaced by the validated v2 implementation (`faclib/orbital.c`, see
  `research/debye-screening/patches/01_*.patch`).
  - In the merge, the conflict in `faclib/orbital.c` was resolved by keeping
    v2 completely. None of this branch's C code is active.
- **Close agreement.** For the Fe II–IV reference configurations, the two
  implementations agree to ≤ 0.2 meV (`docs/numerical-validation.md`).
  - The results here are therefore probably still representative, but they
    have **not** been recomputed with v2.
- **Duplicates dropped.** Files that this branch had only moved (the Chen
  2018 inputs, `debye_fac.py`, the two old patches, …) were byte-identical to
  the copies already in `research/debye-screening/`. Only those copies are
  kept.
- **Old paths.** The scripts still use the paths of that branch
  (`debye_analysis/...`), so run them from a copy of the old layout.
