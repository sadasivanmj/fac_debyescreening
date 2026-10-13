# Development notes: branches and how the work was done

All work is in the repository at `/data/project/fac`. Each step was done on
its own branch in a separate git worktree. The main working tree
(`fix-debye-ee-screening`, c4a016a) was never modified.

| branch | tip | contents |
|---|---|---|
| `fix-debye-ee-screening` | c4a016a | starting point: upstream FAC + screened Slater integrals + nuclear screening (unscreened SCF direct potential) |
| `feature/screened-scf-v2` | ee326b0 | 00eb229: v2 screened `DensityToSZ` + μ = 0 regression test; ee326b0: `verify_ip_deprince.py` |
| `fix/hxs-iter0` | f12e3bf | 0650edb: `radial:hxs_iter0`; f12e3bf: `--hxs-iter0` flag for the validation script |
| `fix/debye-radialbound` | 8c8b1bf | Debye upper bracket in `RadialBound` — **validated solver state** |
| `feature/deprince-reference-validation` | 3ae5c5b | Fe II–IV reference-configuration suite (student-style scripts) |
| `chore/repository-organization` | this branch | research files moved into `research/`, `docs/`, `archive/` (no code or result changes) |
| `feature/full-active-space-deprince` | 7ba19fd | full single/double active-space driver (copied to `research/full-active-space/scripts/`) |
| `diagnostic/sturm-node-count` | 0bdc147 | Sturm node-count instrumentation (archived as a patch) |
| `experiment/warm-start-continuation` | c5f3d1d | warm-start experiments (archived) |
| `ee-screening-scf`, `fix-debye-nuclear-screening` | — | older screening attempts; not used |

Rules followed throughout:

- **Building.** Build with `./configure --prefix=$PWD/_noinstall && make pfac`
  only. Python always uses `PYTHONPATH=<worktree>/build/lib...`, and each
  script prints `fac.__file__` to prove which build it used.
- **No C changes without a regression.** Every C change was checked with the
  μ = 0 bitwise regression before it was committed.
- **Failed runs are kept.** They are renamed (`.attemptN`) or archived, never
  deleted.

Lessons worth remembering:

- **`PYTHONPATH=<repo>/python` is wrong.** It imports the globally installed
  `pfac`, not the build in the repository.
- **FAC's SCF restarts its exchange ramp on every `OptimizeRadial`.** With
  screening, that exchange-free first iteration can unbind 3d; this is why
  `hxs_iter0` exists.
- **`RestorePotential` clears the orbital table**, so it is not a warm start.
- **`fac.Structure(p, [])` crashes PFAC.** Use `[-1]` for "all J".
- **The OpenMP build** needs `fac.InitializeMPI(n)` to run in parallel at all
  (about 2× speed-up for Fe II). It also segfaults on screened runs.
- **Python's `csv.writer` writes `\r\n`** unless `lineterminator="\n"` is
  set.
