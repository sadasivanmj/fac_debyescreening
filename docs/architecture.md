# Where the Debye screening lives in FAC

The physics: every Coulomb interaction is replaced by a Yukawa one,

    −Z/r → −Z e^(−μr)/r          (electron–nucleus)
    1/r₁₂ → e^(−μ r₁₂)/r₁₂       (electron–electron)

with μ = 1/λ_D, and λ_D = √(T / 4πn(1+ups)) in atomic units.

## Switching it on (Python)

```python
fac.SetOption("orbital:debye_mode", 1)   # electron-nucleus screening
fac.SetOption("radial:ee_screen", 1)     # electron-electron screening
fac.SetOption("radial:hxs_iter0", 8)     # warm SCF restarts (see below)
fac.PlasmaScreen(zion, n_A3, T_eV, 1, 0.0, 1)   # m = 1 Debye, ups = 0
```

`PlasmaScreen` converts n from Å⁻³ to a₀⁻³ and T from eV to Hartree.
`InitializePS` (`orbital.c`) then sets `dps` = λ_D. The scripts check
`dps · μ = 1` in the dumped potential (`GetPotential`).

## The pieces

| piece | file / function | what it does |
|---|---|---|
| nuclear screening | `orbital.c: SetPotentialIPS` | adds ZPS(r) = Z(1 − e^(−r/λ)), so the orbital potential becomes −Z e^(−μr)/r |
| μ used by the e–e code | `radial.c: EEScreenMu` | returns 1/dps when `ee_screen = 1`, `mps = 1` and `dps > 0`; otherwise 0 |
| screened Slater integrals | `radial.c: GetYk → GetYkDebye` | Yukawa multipole kernel μ(2k+1) i_k(μr<) k_k(μr>) for all Slater integrals (`fac-debye-ee-screening.patch`) |
| screened SCF direct potential (v2) | `orbital.c: DensityToSZ`, helpers `_ee_E`, `_ee_screened_sz` | r·V(r) = Ã(r) + r E(μr) B̃(r), with E(x) = (1 − e^(−2x))/2x, so every factor stays bounded; integrated in blocks with μΔr ≤ 100; aborts if μ is NaN, negative or infinite; for μ = 0 the original Coulomb code runs unchanged (`01_*.patch`) |
| exchange ramp on restarts | `radial.c: OptimizeILoop`, option `radial:hxs_iter0` | FAC ramps exchange in as 1 − e^(−0.75·iter) and restarts at iter = 0 in every `OptimizeRadial`; x > 0 starts the ramp at iter + x (`02_*.patch`) |
| Debye upper bracket | `orbital.c: RadialBound`, option `orbital:debye_emax_eps` | if `mps = 1` and the upward search reaches E > −ε, count nodes over the whole grid; nr+1 nodes → use −ε as the upper bracket, otherwise return −6 with a clear message (`03_*.patch`) |

## Things that are *not* screened or changed

- **Exchange** (local exchange potential), Breit and QED terms are unscreened.
- **Ion-sphere cut.** With m = 1, `InitializePS` still computes an ion-sphere
  radius. If the radial grid is built *after* `PlasmaScreen`, `DensityToSZ`
  is cut there.
  - The scripts always run a free-ion `OptimizeRadial` first, so the grid is
    fixed with no cut (`ips = 0`).

## Checking a run

The potential header written by `GetPotential` has `mps`, `dps`, `ips`, `HXS`.
The last grid point gives the asymptotic charge r(V_c + U) + ZPS:

- `0` when both screenings are active;
- −(zion + 1) for a free ion;
- +(N − 1) when only the nucleus is screened.
