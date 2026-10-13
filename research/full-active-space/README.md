# Full single/double active space for Fe II–Fe IV (UNFINISHED)

**Goal:** see how much of the gap between the FAC reference-space IPs and
Deprince et al. (2020) Table 2 comes from the small configuration space,
compared with the orbitals or the screening.

**Status:** unfinished. 5 of the 18 serial jobs have converged, and no IP has
been computed yet. See [`../../PROJECT_STATUS.md`](../../PROJECT_STATUS.md).

## Model

- **Active space.** For each Table 1 reference configuration, all
  configurations within one or two electron moves among
  {1s, 2s, 2p, 3s, 3p, 3d, 4s}. Each reference is expanded separately, and
  the union is taken.
  - This is the paper's recipe as written, but it gives **more CSFs than
    published**, and the restriction that would explain the difference could
    not be found (see the table below).
- **CI.** FAC `Structure` over all **even-parity** J blocks, using
  `fac.Structure(0, [-1])`. E0 is the lowest level of all blocks.
  - Odd parity needs a 2p or 3p hole, so it lies tens of eV higher.
- **Orbitals.** Same `ref` / `gref` models and the same screening sequence as
  the reference-configuration suite.
  - The extra configurations enter the CI only, not `OptimizeRadial`.
  - For Fe IV, the 4s orbital is not optimized: it is a spectator in the SCF
    potential.

| ion | configurations | CSFs (both parities) | even CSFs | largest even J block | published CSFs |
|---|---:|---:|---:|---:|---:|
| Fe II | 120 | 15147 | 7603 | 1854 (2J = 5) | 12268 |
| Fe III | 153 | 40189 | 20305 | 4489 (2J = 4) | 31501 |
| Fe IV | 153 | 75309 | 38157 | 8239 (2J = 5) | 57189 |

FAC's level counts per J block equal these counted CSFs exactly; the script
checks this in every job.

## Script

`scripts/verify_ip_deprince_full_as.py` is copied unchanged from branch
`feature/full-active-space-deprince` (commit 7ba19fd).

- One FAC process per job.
- Each converged job is saved at once as an atomic `result.json` checkpoint.
- A job is reused only if its input hash matches. The hash covers the
  configurations, options, extension sha256, script sha256 and git HEAD.
- Failed attempts are renamed to `.attemptN` and kept.
- Jobs are scheduled under a memory budget.

Example, from the repository root with the serial build (`make pfac`):

```bash
export PYTHONPATH=$PWD/build/lib.linux-x86_64-cpython-310
python3 research/full-active-space/scripts/verify_ip_deprince_full_as.py --audit-only \
        --out /tmp/fas_audit                                  # active-space counts only
python3 research/full-active-space/scripts/verify_ip_deprince_full_as.py --space sd \
        --jobs 3 --mem-budget-gb 11 --out research/full-active-space/results/<new-folder>
```

Because the git HEAD is part of the hash, the old checkpoints below are **not**
reused automatically from this branch. A rerun into their folder renames
them to `.attemptN` rather than deleting them. A new output folder is
simpler.

## Results so far (`results/`, copied from `/tmp/fac-full-as-worktree`)

| folder | what | status |
|---|---|---|
| `serial-partial-2026-10-13/runs/sd/` | **production attempt, serial build** | 5 converged; Fe II μ = 0.1 `gref` interrupted |
| `serial-partial-2026-10-13/runs/ref/` | reference-space rerun with the same driver | 18/18 converged; identical to the reference suite |
| `openmp-pilots-2026-10-13/benchmarks/` | Fe II and Fe IV μ = 0 pilots, thread scaling, segfault reproduction | pilots converged |
| `openmp-pilots-2026-10-13/openmp-production-stopped/` | first production attempt with the OpenMP build | stopped after the screened jobs segfaulted |

Converged full-space energies (eV; serial build unless marked):

| ion | model | μ | E0 | ground 2J | wall time | peak memory |
|---|---|---:|---:|---:|---:|---:|
| Fe II | ref | 0 | −34564.7462 | 5 | 47 s | 1.2 GB |
| Fe II | gref | 0 | −34571.6567 | 9 | 48 s | 1.2 GB |
| Fe II | ref | 0.1 | −33629.6069 | 5 | 43 s | 1.2 GB |
| Fe III | ref | 0 | −34548.1596 | 6 | 535 s | 3.1 GB |
| Fe III | gref | 0 | −34556.1514 | 8 | 536 s | 3.1 GB |
| Fe IV | ref | 0 | −34523.5868 (OpenMP, 4 threads) | 5 | 1397 s | 9.0 GB |
| Fe IV | gref | 0 | −34528.3851 (OpenMP, 4 threads) | 5 | 1397 s | 9.0 GB |

- **Agreement with other runs.** Where both builds were run (Fe II μ = 0;
  Fe III μ = 0 `ref`), the OpenMP and serial energies are identical. For the `ref`
  model, the 2J = 9 block of Fe II and the 2J = 5 block of Fe IV reproduce the
  earlier ground-symmetry-block results of the archived `sfac` study
  (−34561.8279 vs −34561.8278 eV, and −34523.5868 eV).
- **Different ground state with `ref` orbitals.** In the full space the Fe II
  ground state becomes 3d⁵4s² (J = 5/2), and Fe III becomes 2J = 6. This is
  the `ref` orbital-averaging artefact, kept as computed.
- **Large files stay local.** Level tables (`e.lev`), binaries and
  potentials are kept on disk but git-ignored. Manifests, results, job
  definitions, active-space lists and logs are committed.

## Next step

Rerun the production with the serial build into a new folder. Fe IV runs one
job at a time, at up to 9 GB. Then compute the 12 IPs from `ip_sd.csv`.
