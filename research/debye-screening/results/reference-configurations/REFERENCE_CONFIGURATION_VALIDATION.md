# Fe II–Fe IV reference-configuration calculations with Debye screening

## Purpose

Before the full single/double-excitation active space is run, these
calculations establish a clean, reproducible baseline. Ground-state energies
of Fe II, Fe III and Fe IV are computed from the reference configurations of
Deprince et al., A&A 643, A57 (2020), Table 1. They use μ = 0, 0.1 and
0.25 a.u. and two orbital models. Twelve ionization potentials (IPs) follow
from them. The run also checks that the validated FAC patches reproduce
earlier results.

**Scope.** Only the reference configurations enter the calculation. The
full single/double-excitation active space of the paper has *not* been run
here, so these numbers are not a reproduction of Deprince et al. Table 2.

## Code

- **Branch:** `feature/deprince-reference-validation`, based on
  `fix/debye-radialbound` (commit `8c8b1bfb`).
- **C code:** unchanged from that branch. Relative to `c4a016ad` it contains
  three patches:
  - the v2 Debye-screened SCF direct potential (`DensityToSZ`);
  - the `radial:hxs_iter0` exchange-ramp option;
  - the Debye-only upper energy bracket in `RadialBound`.
- **Build:** `./configure --prefix=$PWD/_noinstall && make pfac`. The serial
  build is used and nothing is installed; the scripts import `pfac` from
  `build/lib.linux-x86_64-cpython-310`.

## Reference configurations (Table 1, holes written out)

| Ion | electrons | configurations |
|---|---:|---|
| Fe II | 25 | 3d⁶4s, [3p]3d⁷4s, [3p]3d⁶4s², [2p]3d⁷4s, [2p]3d⁶4s², [1s]3d⁷4s, [1s]3d⁶4s² |
| Fe III | 24 | 3d⁶, [3p]3d⁷, [3p]3d⁶4s, [3p]3d⁵4s², [2p]3d⁷, [2p]3d⁶4s, [2p]3d⁵4s², [1s]3d⁷, [1s]3d⁶4s, [1s]3d⁵4s² |
| Fe IV | 23 | 3d⁵, [3p]3d⁶, [2p]3d⁶, [1s]3d⁶ |

The explicit shell occupations are in the scripts. Each script checks every
electron count before calling FAC.

## Orbital models

The CI space is the same in both models: all reference configurations, all
J, both parities. Only the orbitals differ.

- **`ref`:** `OptimizeRadial` on the average of *all* reference
  configurations. FAC weights each configuration equally in this average. The
  core-hole configurations therefore dominate it, and the valence orbitals are
  poor (see Limitations).
- **`gref`:** `OptimizeRadial` on the ground configuration only. The other
  references are added only in `Structure`.

## Screening settings

- **Plasma model:** `PlasmaScreen(zion, 1e-6 Å⁻³, T, 1, 0.0, 1)`. That is
  Debye mode m = 1 and ups = 0, with λ_D = √(T / 4πn) in atomic units.
- **Temperature:** T is chosen so that λ_D = 1/μ.
- **Density:** chosen small, so that the ion-sphere radius lies far outside
  the atom.
- **Electron–nucleus screening:** `orbital:debye_mode 1`.
- **Electron–electron screening:** `radial:ee_screen 1`. This screens the
  Slater integrals and the SCF direct potential.
- **Sequence:**
  1. A free-ion `OptimizeRadial` runs first. It also fixes the radial grid, so
     no ion-sphere cut-off is applied.
  2. Screening is switched on in steps of μ = 0.05. Each step starts from the
     previous orbitals, with `radial:hxs_iter0 = 8`.
- **μ = 0:** the default unscreened sequence; no screening options are set.
- **Runtime checks** (from `final.pot`), each of which stops a run if it fails:
  - `mps = 1` and `dps · μ = 1` to 10⁻⁶.
  - The asymptotic r·V is 0, which shows that e–e screening is active.
  - `ips = 0`: no ion-sphere cut.
  - For the free ion, r·V → −(zion+1).
  - `NELE` in the level file is correct, and E0 is finite.
  - The log contains no "Maximum iteration reached in OptimizeRadial".

## How to run

The scripts are in `research/debye-screening/cases/fe/`. Runs and tables are
written to this folder (`research/debye-screening/results/reference-configurations/`).
Commands from the repository root:

```bash
research/debye-screening/cases/fe/run_reference_suite.sh   # 18 runs, about 1 minute;
                                                          # refuses to overwrite reference_runs/
# or a single run, from this results folder:
cd research/debye-screening/results/reference-configurations
PYTHONPATH=../../../../build/lib.linux-x86_64-cpython-310 \
    python3 ../../cases/fe/fe_iii_reference.py gref 0.1 my_test_run
python3 ../../cases/fe/collect_reference_results.py   # rebuilds the two CSV tables
```

Each run writes the following to `reference_runs/<ion>/<model>_mu<mu>/`:

- the FAC level table `levels.txt` and binary `levels.b`;
- the potential `final.pot`;
- `summary.txt`, with the ion, electron count, model, requested and used μ,
  `hxs_iter0`, E0 and status.

Its log is `reference_runs/<ion>/<model>_mu<mu>.log`.

## Results

All 18 calculations converged and passed every check.

- **Units:** FAC writes level energies in eV, so the IPs are differences of
  E0 in eV. Hartree values use 1 Ha = 27.211386 eV.
- **Data files:** `reference_energies.csv` and
  `reference_ionization_potentials.csv`.

### Ground-state energies E0 (eV)

| ion | model | μ = 0 | μ = 0.1 | μ = 0.25 |
|---|---|---:|---:|---:|
| Fe II | ref | −34561.2842 | −33625.2303 | −32282.8180 |
| Fe II | gref | −34571.2747 | −33635.8239 | −32293.7705 |
| Fe III | ref | −34542.4125 | −33611.0422 | −32274.1005 |
| Fe III | gref | −34555.7544 | −33625.5981 | −32289.0833 |
| Fe IV | ref | −34522.7607 | −33599.8850 | −32274.0752 |
| Fe IV | gref | −34527.6602 | −33605.4643 | −32279.8135 |

### Ionization potentials (eV); ΔIP = IP(μ) − IP(0)

| IP | model | μ = 0 | μ = 0.1 | μ = 0.25 | ΔIP(0.1) | ΔIP(0.25) |
|---|---|---:|---:|---:|---:|---:|
| Fe II → Fe III | ref | 18.8717 | 14.1881 | 8.7175 | −4.6836 | −10.1542 |
| Fe II → Fe III | gref | 15.5203 | 10.2258 | 4.6872 | −5.2945 | −10.8331 |
| Fe III → Fe IV | ref | 19.6518 | 11.1572 | 0.0253 | −8.4946 | −19.6265 |
| Fe III → Fe IV | gref | 28.0942 | 20.1338 | 9.2698 | −7.9604 | −18.8244 |

The published values are given for context only; these runs are not their
equivalent.

| | μ = 0 | μ = 0.1 | μ = 0.25 |
|---|---:|---:|---:|
| Deprince Table 2, Fe II → Fe III | 18.89 | 14.20 | 7.78 |
| Deprince Table 2, Fe III → Fe IV | 32.93 | 25.60 | 16.42 |
| NIST (isolated ion), Fe II / Fe III | 16.19920 / 30.651 | | |

## Comparison with previous results

These runs reproduce the earlier validated reference-space results:

- **Energies:** all 18 E0 values are identical to every printed digit to
  `verify_ip_deprince_out_hxs8_rb/energies.csv` (branch
  `fix/debye-radialbound`, same C code).
- **IPs:** all 12 IPs agree with the regression table (18.87, 14.18, 8.718 /
  15.520, 10.226, 4.687 / 19.652, 11.157, 0.025 / 28.094, 20.134, 9.270 eV)
  to the quoted precision.
- **Independent cross-check:** the earlier `ee-screening-scf` runs, made with
  a separately implemented screened potential, agree to ≤ 0.0002 eV.

## Known limitations

- **Reference configurations only.** No single or double excitations are
  included. The paper's active space has 12 268 / 31 501 / 57 189 CSFs; the
  reference spaces here have 986 / 1640 / 460 (all J, both parities).
- **`ref` orbital-averaging artefact (kept, not corrected):**
  - The `ref` Fe III → Fe IV IP is 13.3 eV below the paper at μ = 0. It falls
    to 0.025 eV at μ = 0.25.
  - Cause: the equal-weight average over the core-hole configurations gives
    poor valence orbitals.
  - Its ΔIP should not be interpreted physically.
  - The `gref` model is the more meaningful reference-space result.
- **`gref` offsets:** `gref` is 3.1–4.0 eV (Fe II) and 4.8–7.2 eV (Fe III)
  below the paper. The remaining gap needs correlation (the active space) and
  better orbital optimization.
- **OpenMP:** only the serial FAC build is used. An isolated OpenMP build
  (`--with-mpi=omp`) segfaults in `PotentialHX` at the first screened
  `OptimizeRadial`. It has not been investigated further.
- **Full active space not run:** the full single/double-excitation active
  space has not been run as part of this stage.
