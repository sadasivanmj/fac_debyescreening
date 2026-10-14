CLAUDE CHANGE HISTORY - FAC DEBYE-SCREENING INVESTIGATION (Deprince et al. 2020 Table 2)
=====================================================================================
Audit written: 2026-10-12 ~12:53-13:00 UTC (date -u at start of documentation task: 2026-10-12 12:52:46 UTC).
Fe VI job status at audit time: NOT running. `pgrep -x sfac` returned no process at 12:52:46 UTC.
  All three Fe VI reruns (runs/sdp/FeVI/{mu0,mu0p1,mu0p25}) had already completed with status ok
  at 12:48 UTC (see logs_rerun_sdp_FeVI.txt and their run.json files).
  One helper process of mine was still alive: `python3 throttle_sfac.py 3` (pid 2662028). It only
  pauses/resumes sweep sfac processes and exits by itself after 10 min with no sweep sfac process.
  It was not touched during this documentation task.
Author of this file: Claude (Claude Code session a8e208a0-910c-4f2e-b650-3503c715a8c5, Linux host).
Scope: everything below refers to work done in THIS session on the Linux machine, path /data/project/fac.
  I have NO record of, and did not perform, any macOS work (see section 6).

Source-type labels used below:
  [RAW]   verified from raw FAC output files (e.lev, final.pot, in.sf, sfac.log, error.pot)
  [CSV]   taken from a CSV/JSON/markdown written by a script
  [SUB]   independently verified subtraction of ground-state energies (audit.py, raw files)
  [SESS]  recorded in this session's transcript (tool outputs seen by me). The transcript file is
          /home/sadasivanmj/.claude/projects/-data-project-fac/a8e208a0-910c-4f2e-b650-3503c715a8c5.jsonl
  [SRC]   read directly from source code / git objects
  [INF]   my inference / interpretation (not a measurement)
  Times are UTC. Exact times come from file mtimes or tool output; others are approximate ordering
  reconstructed from the session transcript and are marked "~".


1. PROJECT OBJECTIVE AND ENVIRONMENT
------------------------------------
Objective (from the user's first task): reproduce Table 2 of J. Deprince et al., A&A 643, A57 (2020)
("Plasma environment effects on K lines ... IV"), i.e. ionization potentials (IPs) of Fe II-Fe VIII
for Debye screening parameter mu = 0, 0.1, 0.25 a0^-1, with native SFAC (no PFAC). Required:
diagnose the screened radial-solver failure first, write a reproducible native-SFAC driver, a
focused diagnostic test, CSV and Markdown outputs, and compute relative errors vs the paper and vs
NIST. Follow-up tasks: (2) a read-only audit of the IP discrepancies; (3) this change history.

Reference values (Deprince et al. Table 2, quoted, not tuned to):
  Transition          mu=0     mu=0.1   mu=0.25
  Fe II  -> Fe III    18.89    14.20     7.78
  Fe III -> Fe IV     32.93    25.60    16.42
  Fe IV  -> Fe V      55.37    44.65    29.64
  Fe V   -> Fe VI     74.92    61.56    42.83
  Fe VI  -> Fe VII    98.20    82.20    59.79
  Fe VII -> Fe VIII  123.60   104.96    78.86
  Fe VIII-> Fe IX    151.14   129.85    99.98
  NIST (from the paper): 16.19920(5), 30.651(12), 54.91(4), 75.0(2), 98.985(15), 124.98(1), 151.060(12) eV.

Environment:
  Host: Linux 6.8.0-142-generic, 16 cores, 15 GB RAM, 4 GB swap. User sadasivanmj.
  Repository: /data/project/fac, branch fix-debye-ee-screening, HEAD c4a016ad8cf0458366714b013727d915e1fdec80
    ("fe_ip", author Sadasivan M J, 2026-10-12 13:00 +0530). Untracked at session start:
    c4_ground_results/, fe_calc/fe_ip_sfac_runs/ (and the Deprince PDF).
  Paper PDF: /data/project/fac/Fe2_Fe8 in plasma environment_Deprince et al.pdf (text extracted with
    pdftotext into scratchpad/deprince.txt) [SESS].
  The archive fac_debyescreening.zip requested by the user was NOT FOUND anywhere on the
    filesystem (`find / -xdev -iname '*.zip'` and `-iname '*debye*'`). It was never inspected [SESS].


2. CHRONOLOGICAL TIMELINE
-------------------------
(UTC. "mt" = file mtime evidence. Order otherwise from the session transcript.)
~07:40  Read-only inspection: git log/branches; fe_calc/fe_ip_sfac.py; fe_calc/fe_ip_sfac_runs/
        {FeII_free,FeIII_free,FeII_mu0p1}. Original failure log FeII_mu0p1.log:
        "Max iteration before finding correct nodes in RadialBound a: 0 0 279 3 2 -8.12838e-154 /
         Error occured in RadialSolver: -2 3 2 0.0000000000E+00" [RAW].
~07:41  Found installed ~/facinst/bin/sfac is bit-identical to /data/project/fac/sfac/sfac, mtime
        2026-09-14 04:43; libfac.a 2026-09-14 04:42; many faclib sources newer than libfac.a [SESS].
~07:42  objdump of DensityToSZ in the installed binary: it calls EEScreenMu and sinh -> binary contains the
        ee-screening-scf DensityToSZ change, not fix-debye-ee-screening code [SESS].
07:43   scratchpad/orbital_scf.c = `git show ee-screening-scf:faclib/orbital.c` (read copy, mt).
~07:44  `git worktree add --detach /data/project/fac_worktrees/ee c4a016a` and
        `git worktree add --detach /data/project/fac_worktrees/scf ee-screening-scf` (ff9a5c9).
07:46   Both built: ./configure --prefix=<worktree>/inst; make -j8; make install (build_status.txt:
        "ee exit=0", "scf exit=0"; mt 07:46). Result: scf clean build sha256 2789ca41... == installed
        binary. ee clean build sha256 279d393e... [SESS].
~07:47-07:55  Diagnostics with installed binary in fe_calc/deprince_t2/diagnostics/scratch_nps
        (nps scan, potential dumps, free vs e-n-only) -> cause 1 (ion-sphere truncation) identified.
~07:55  Patch v1 (RadialBound, orbital.c) written in ee worktree, copied to scf; unpatched binaries
        copied to inst_unpatched/; rebuild. Tests: still "unbound" at -1e-6, -1e-4, -1e-2, -5e-2.
~07:58  Option radial:hxs_iter0 v1 added (radial.c); rebuild. Tests: free-ion E0 changed by 0.30 eV
        with hxs_iter0=20; screened cold starts still failed -> radial.c change reverted
        (`git checkout faclib/radial.c` in both worktrees).
~08:00  Patch v2 of RadialBound (full-grid Sturm node count, skip phase 5); rebuild; tests.
08:02   hxs_iter0 re-added (v2, identical arithmetic, different comment); rebuild (mt 08:02:
        ee/inst sfac 23b5644a..., scf/inst sfac 3b88ec91..., later copied to scf/inst_fix_no842).
~08:03-08:06  path-dependence / ladder tests (diagnostics/scratch_nps/pathdep, ladder3).
~08:06  Driver fe_calc/deprince_t2/deprince_t2_sfac.py created (models ref, sdv).
08:08   ref sweep, 24 runs, all ok (logs_ref_sweep.txt mt 08:08), binary 3b88ec91.
        Fe III anomaly observed (IP 19.65 eV). Ground-only tests (diagnostics/scratch_feiii).
~08:09  Models gref, gsdv added to driver; sweep of gref,sdv,gsdv (logs_sweep2.txt mt 08:10), 72 ok.
08:10   diagnose_radial_failure.py written (mt 08:08) and run (bhxa3j0w3, mt 08:10).
08:10-08:11  count_csfs.py written and run (results/csf_counts.csv mt 08:11).
~08:12  Models sdp,gsdp added; mu=0 test sweep (logs_sdp_mu0.txt mt 08:20): FeIV, FeVI, FeVII
        segfaulted (sfac exit -11). gdb backtrace: mmalloc <- AllocHamMem <- ConstructHamilton.
08:21   Backport of upstream d3a4ec3 into scf worktree (mmalloc.c via git apply --include;
        structure.c AllocHamMem by hand); rebuild -> FINAL binary scf/inst/bin/sfac sha256 3a77793e...
        (mt 08:21:26). scf_final_working_tree.diff written (mt 08:21).
~08:2x  Standalone FeIV paper-AS retest started (task b1cpcwifo); later killed by me (rc 143, ~08:39).
08:32   Fast models (ref,sdv,gref,gsdv) rerun with final binary: 96 ok (logs_final_fast.txt mt 08:32).
~08:3x  sdp,gsdp sweep started with final binary, --jobs 6 (logs_final_sdp.txt, task bxwqvqr5z).
08:35   Diagnostics rerun with final binary (logs_diagnostics.txt); patches/ and provenance.json
        written (mt 08:35).
~08:3x-08:5x  small_mu, branch_ee, grid, warm_restart, first_order diagnostics (see section 5).
08:39   deprince_t2_sfac.py last edited (monotonicity check + docstring; mt 08:39); verified that
        generated in.sf for ref/FeII/mu0.1 is byte-identical to the existing run input.
~09:00-09:12  Memory pressure (13-14 GB used). 09:12 I sent SIGTERM to the three sdp FeVI processes
        (pids 2654276/2654277/2654323, ~34 min of compute each; recorded "sfac exit -15",
        t=2042-2049 s). 09:12:48 started throttle_sfac.py 3 (pauses extras with SIGSTOP).
09:13   REPORT.md last edited (mt 09:13). It is INCOMPLETE: still contains literal "RESULTS_PLACEHOLDER".
10:03-11:20  gsdp FeVI, FeVII, FeVIII completed; sdp FeVII resumed at 11:20 (monitor logs).
~11:21  Rerun of sdp FeVI launched: `python3 deprince_t2_sfac.py --models sdp --ions FeVI --jobs 3`
        (task b1qqa6pnn). Its 3 sfac processes were paused by throttle until 12:00:43-12:01:16.
12:01   Main sdp/gsdp sweep finished: 45/48 printed ok, 3 failed = the stopped FeVI
        ("runs: 141 ok, 3 failed" across all models; logs_final_sdp.txt mt 12:01).
12:13   Read-only audit: scratchpad/audit/audit.py written and run -> raw.json, ip_audit.csv
        (FeV/FeVI sdp pending at that time).
12:48   FeVI reruns finished (runtime_s 5199.7-5249.5 s wall, which includes ~40 min paused);
        driver regenerated results/{energies.csv,ip_comparison.csv,summary.md} (mt 12:48).
12:48   audit.py rerun (same script) -> raw.json and ip_audit.csv OVERWRITTEN with the complete set
        (126/126 IP rows ok). This overwrite happened in the previous task, before this one.
12:52   This documentation task started.


3. FILES CHANGED AND EXACT DIFFS
--------------------------------
3.0 Authorship summary
  * By me (Claude): everything under /data/project/fac/fe_calc/deprince_t2/ (new directory),
    /data/project/fac_worktrees/ (new worktrees, builds, patch files), the git worktree metadata
    entries in /data/project/fac/.git/worktrees/{ee,scf}, the scratchpad files.
  * NOT changed by me: the main repository working tree (no tracked file modified, no commit, no
    branch switch, no stash), ~/facinst (installed binary untouched), fe_calc/fe_ip_sfac.py,
    fe_calc/fe_ip_sfac_runs/, c4_ground_results/, any Python environment, PFAC.
  * Pre-existing user work (authorship: Sadasivan M J per git): commits c768c24 "Fixing Screening",
    4dd4b40 (on ee-screening-scf), 49c7b97 merge, c4a016a "fe_ip" (adds fe_ip_sfac.py, PFAC scripts,
    .DS_Store files). A stash "stash@{0}: WIP on ee-screening-scf: ff9a5c9" existed before I started;
    I did not touch it.

3.1 /data/project/fac_worktrees/ee/faclib/orbital.c  (worktree of c4a016a; uncommitted)
3.2 /data/project/fac_worktrees/scf/faclib/orbital.c (worktree of ff9a5c9; uncommitted)
  Original: git blob at the respective commit (recoverable: `git -C <wt> diff`, or `git show
  c4a016a:faclib/orbital.c`, `git show ff9a5c9:faclib/orbital.c`).
  History of my edits (both trees received the same edits; ee edited first, then patch applied to scf):
   v1 (~07:55, SUPERSEDED): in RadialBound phase a, inside the upward energy scan, when
      pot->mps==1 and e > -_debye_emax_eps: set e=-eps, compute bqp=DpDr(...), count nodes with
      IntegrateRadial(p,e,pot,0,bqp,i2,1.0,2) up to the turning point i2; if nodes>nr break, else print
      "Orbital unbound in Debye potential in RadialBound" and return -6. Added static
      _debye_emax_eps=1e-6 and SetOption "orbital:debye_emax_eps".
      Tests: Fe II ground mu=0.1 (nps 1e-2, 1e-6), eps 1e-4,1e-2,5e-2: all returned -6 [SESS].
      Later shown to give FALSE NEGATIVES near E=0 (turning-point-limited count misses a node
      beyond the classical turning point). Superseded by v2. The exact v1 text is recoverable only
      from the session transcript (not stored as a file).
   v2 (~08:00, FINAL, present now): same gate, but count nodes of the regular solution over the
      whole grid: nodes = IntegrateRadial(p, e, pot, 0, 0.0, pot->maxrp-1, 1.0, 0); accept
      e=-eps as upper bracket only if nodes == nr+1 (sets debye_emax=1), else print
      "No Debye bracket in RadialBound (unbound|higher state bound): ..." and return -6; when
      debye_emax is set, phase 5 (upward search for nodes==nr+1) is skipped (wrapped in
      `if (!debye_emax) { ... }`). Plus option orbital:debye_emax_eps (default 1e-6 Ha).
      Exact diff: /data/project/fac/fe_calc/deprince_t2/patches/debye_sfac_fixes.patch (orbital.c
      part, +32 lines). Non-Debye runs (mps != 1) are untouched by construction.
  Tests after v2: cold start Fe II mu=0.1 still -6 (genuine transient, see 5.D1); warm ladders OK;
      free-ion E0 identical to unpatched (5.E1/E2). Present in both worktrees now; compiled into
      ee/inst (23b5644a) and scf/inst (3a77793e, the FINAL binary used by all 144 production runs).

3.3 /data/project/fac_worktrees/{ee,scf}/faclib/radial.c (uncommitted)
   hxs_iter0 v1 (~07:58, REVERTED): replaced `ahx = exp(-iter*0.75);` in OptimizeILoop by
      `ahx = exp(-(iter+_hxs_iter0)*0.75);`, added `static double _hxs_iter0 = 0.0;` and SetOption
      "radial:hxs_iter0". Tests: free Fe II ground cold start E0 = -34571.2747 (default),
      -34571.2550 (2), -34570.9771 (20) eV, i.e. path dependence up to 0.30 eV; screened cold
      starts with 2/4/20 still failed at SCF iteration 2. Reverted with `git checkout faclib/radial.c`.
   hxs_iter0 v2 (08:02, FINAL, present now): same arithmetic, comment "starts the exchange ramp as
      if that many iterations had passed, for warm-started optimizations"; used ONLY in warm
      steps (value 8) after a cold free-ion OptimizeRadial. Diff in patches/debye_sfac_fixes.patch
      (radial.c part, +9/-1). Default 0 = stock behaviour.
   NOTE: I did NOT modify PotentialHX, PotentialHX1, SetPotential or DensityToSZ. The ramp lives in
      OptimizeILoop. Those functions were only READ (section 5). No change to jmaxk anywhere.

3.4 /data/project/fac_worktrees/scf/faclib/mmalloc.c and structure.c (08:21, uncommitted, FINAL)
   Purpose: upstream fix d3a4ec3 ("fix int overflow in large-size calculation (int->size_t) #842")
      is on fix-debye-ee-screening but not on ee-screening-scf; without it AllocHamMem/mmalloc
      segfaulted for Hamiltonian blocks > ~8192 CSFs.
   mmalloc.c: `git apply --include='faclib/mmalloc.c' ../upstream_d3a4ec3_int_overflow.patch`
      (22 lines changed, identical to upstream).
   structure.c (AllocHamMem): the upstream hunk did not apply (context differs); replaced by hand:
      -  t = hdim*(hdim+1)/2;              +  t = (size_t)hdim * (size_t)(hdim+1)/2;
      -  h->hsize = t + hdim*jp + jp;      +  h->hsize = t + (size_t)hdim*(size_t)jp + (size_t)jp;
   Full exact diff of the scf tree that produced the final binary:
      patches/scf_build_working_tree.diff (sha256 ceac162e...; also /data/project/fac_worktrees/
      scf_final_working_tree.diff). Present in the scf worktree now. The ee tree already has it.
   Tests: FeIV paper-AS no longer segfaults; small-model IPs identical before/after (ref Fe II 18.872,
      Fe III 19.652 etc. [CSV]).

3.5 New files created by me (all new; nothing pre-existing overwritten except as noted)
   /data/project/fac/fe_calc/deprince_t2/
     deprince_t2_sfac.py  driver (edited several times via python string replacement and Edit:
        +gref/gsdv models (~08:09), +sdp/gsdp + Config chunking (~08:12), +monotonicity check and
        docstring fix (08:39), removed unused CORE constant). Edits never changed generated in.sf
        for existing cases (checked once, ref/FeII/mu0.1) [SESS].
     diagnose_radial_failure.py, count_csfs.py, analyze.py (never run), throttle_sfac.py,
     REPORT.md (incomplete; contains RESULTS_PLACEHOLDER), provenance.json, patches/*,
     logs_*.txt, throttle.log, results/*, runs/<model>/<ion>/<mu>/*, diagnostics/*,
     __pycache__/ (created by imports).
     Overwritten by my own scripts: results/* (by each driver invocation), diagnostics/diagnostics.
     {json,md} (rerun 08:35), runs/* for fast models (rerun with final binary at 08:32; the
     reuse logic deletes outputs before rerunning).
   /data/project/fac_worktrees/: ee/, scf/ worktrees with build products, inst/, inst_unpatched/,
     scf/inst_fix_no842/, build logs, build_status.txt, debye_sfac_fixes.patch (stale copy of v1+...
     from 08:02 - see 3.6), upstream_d3a4ec3_int_overflow.patch, scf_final_working_tree.diff.
   Scratchpad (/tmp/claude-1000/-data-project-fac/a8e208a0-910c-4f2e-b650-3503c715a8c5/scratchpad):
     deprince.txt, orbital_scf.c, audit/audit.py, audit/raw.json, audit/ip_audit.csv, and this file.
3.6 Caveat: /data/project/fac_worktrees/debye_sfac_fixes.patch (mt 08:02) is a copy written at 08:02;
   the authoritative final diff is fe_calc/deprince_t2/patches/debye_sfac_fixes.patch (mt 08:35),
   whose sha256 ef42820f... equals the current ee worktree diff. I did not re-verify the 08:02 copy
   against it in this task.


4. COMMANDS EXECUTED AND OUTCOMES
---------------------------------
All commands ran on the Linux host via Claude Code's Bash tool (no PFAC, no macOS). Full verbatim
commands and outputs are in the session transcript (.jsonl, see header). Key ones:
 Git inspection: git log/branch/status; git show c768c24, 4dd4b40, c3be2b4, d3a4ec3; git diff between
   branches; git merge-base --is-ancestor checks. Outcome: see section 9.
 Worktrees/builds: git worktree add --detach (ee: c4a016a; scf: ff9a5c9);
   ./configure --prefix=<wt>/inst && make -j8 && make install (logs <wt>/build.log,
   build_patched.log, scf/build_842.log). All rebuilds returned rc=0.
 Binary inspection: sha256sum; strings | grep radial:/orbital:; objdump -d (DensityToSZ calls).
 Diagnostics: direct `sfac in.sf` runs in diagnostics/scratch_*; GetPotential dumps; gdb -batch
   backtrace for the segfault.
 Production: python3 deprince_t2_sfac.py --models ... --jobs N (logs_ref_sweep.txt, logs_sweep2.txt,
   logs_sdp_mu0.txt, logs_final_fast.txt, logs_final_sdp.txt, logs_rerun_sdp_FeVI.txt).
 Process control (affected running calculations):
   kill 2649285 2649286 (standalone FeIV retest, ~08:39);
   kill 2654276 2654277 2654323 (sdp FeVI, 09:12, SIGTERM);
   throttle_sfac.py 3 (SIGSTOP/SIGCONT, log throttle.log: paused 2662023/4/7 at 09:12:48, later
   paused rerun pids 2701590/3/9 and resumed them 12:00:43-12:01:16).
 Audit: python3 scratchpad/audit/audit.py scratchpad/audit (12:13 and 12:48).
 PFAC: never built, installed, imported or run by me. Python environments not modified.


5. DIAGNOSTIC EXPERIMENTS AND NEGATIVE RESULTS
----------------------------------------------
Binaries: INST = ~/facinst/bin/sfac (2789ca41, == clean ee-screening-scf ff9a5c9);
 SCF0 = scf/inst_unpatched (2789ca41); FIX = scf/inst (3a77793e final; earlier tests used the
 pre-backport 3b88ec91 which differs only by the int-overflow fix); EE1 = ee/inst (23b5644a, c4a016a+patch).
 Diagnostic suite outputs: fe_calc/deprince_t2/diagnostics/<case>/{in.sf,sfac.log,d.lev,final.pot,
 error.pot}; summary diagnostics/diagnostics.{md,json} (rerun 08:35 with FINAL binary for FIX cases).
 All are reproducible from existing artifacts (inputs kept) [RAW].

A1-A3  Q: do core-hole refs trigger the failure? INST, Fe II mu=0.1, nps=0.01, cold, Table-1 refs,
       ground only, and incremental subsets 1-2..1-7. Expected: maybe a specific ref fails.
       Actual: ALL fail identically (RadialSolver -2, n=3 kappa=2). Establishes: failure independent of
       refs. Does not establish the cause.
B1/B2  Q: e-n only vs e-n+e-e. INST mu=0.1. Both fail; B1 asymptotic r*V = +24 (unscreened Hartree
       tail), fails on 3s; B2 r*V=+14.5. Establishes e-n-only is not a viable SCF model in this code.
C1/C2  Q: role of nps. INST, ground, mu=0.1, nps 1e-2 vs 1e-6 (T scaled). C1: u(r) frozen beyond
       ips=278, r=5.4413 a0 (= ion-sphere radius), r*V -> +14.51 (repulsive tail). C2: r*V -> 0
       but STILL fails. Establishes cause 1 (DensityToSZ truncates at ips when ups==0); shows it is not
       the only cause.
(scratch) Patch v1 eps scan: all "unbound"; later recognised as false negatives of v1.
(scratch) Converged 3d3/2 energies (ground config, FIX/scratch builds): -0.594, -0.562, -0.529 Ha at
       mu=0.02,0.05,0.08 (cold), -0.506 (0.10), -0.436 (0.16), -0.330 (0.25) via ladder [RAW in
       diagnostics/scratch_nps/{muscan,ladder*}].
D1  FIX cold mu=0.1: fails with "No Debye bracket ... (unbound)" at SCF iteration 2; potential at
    that iteration has V_eff(l=2) > 0 everywhere [SESS]. Establishes a genuine SCF transient.
D2/D3/D4  FIX free-start + warm ladder (hxs_iter0=8): ok. Fe II ground mu=0.1 E0=-33635.8239,
    mu=0.25 -32293.7705; all refs mu=0.25 -32282.8180 eV [RAW].
D5  SCF0 (unpatched) ladder: fails (RadialBound a).  D6 FIX ladder with hxs_iter0=0: fails.
    Together: both the solver patch and the ramp offset are needed for the ladder.
E1/E2  free Fe II refs, SCF0 vs FIX: identical E0 -34561.2842 eV. Patch inert for free ions.
E3  INST with ConfigEnergy(0)/(1) (original protocol): -34561.2841 eV, reproduces original baseline.
E4  zps=3 instead of 1 (mu=0.1 ladder): E0 unchanged (-33635.8239).  E5 nps=1e-4: unchanged.
E6  mu=1e-4 ladder (lambda=1e4 a0), Fe II ground: E0=-34570.3335 vs free -34571.2747; shift +0.941 eV
    vs first-order analytic +0.952 eV.
F1  EE1 ladder mu=0.1: fails (r*V -> +24; unscreened DensityToSZ).  F2 EE1 free: -34561.2842 (identical).
Other studies (not in the A-F table):
 * scratch_feiii: Fe II-V optimized on ground config only vs refs CI: IP 15.520/28.094/53.334 eV;
   including refs in CI with ground-optimized orbitals gives identical E0 (refs do not mix).
 * scratch_segv: sdp FeIV mu0 segfault reproduction; gdb backtrace mmalloc<-AllocHamMem; ulimit -s
   unlimited did not help; after backport the run proceeded (I killed that standalone test at ~17 min,
   rc 143, because the identical case ran inside the sweep and completed there).
 * branch_ee (driver with EE1, ref/gref, mu 0,0.1,0.25): free E0 identical to FIX for all ions;
   Fe II-V and all mu=0.25 cases FAILED; Fe VI-IX mu=0.1 converged; IP differences vs FIX
   0.02-0.19 eV (e.g. ref FeVI 81.739 vs 81.860) [CSV diagnostics/branch_ee/results].
 * grid: SetRadialGrid(3000,1.05,100) vs default 1500 points: IP changes <= 0.023 eV
   [diagnostics/grid/grid_sensitivity.txt].
 * warm_restart: extra warm OptimizeRadial at mu=0: E0 changes <= 0.011 eV, IP <= 0.011 eV
   [diagnostics/warm_restart_check.txt].
 * small_mu (mu=1e-4, 1e-3; ref,gref): per-ion E0 shifts within 0.2-1.2% of mu*[Z N - N(N-1)/2];
   IP slopes for ref Fe II/III deviate strongly (wrong sign at 1e-4 for ref Fe II)
   [diagnostics/small_mu_check.txt]. NEGATIVE/CONTRADICTORY observation: unexplained beyond my
   interpretation (orbital non-stationarity) [INF].
 * first_order: H-, He-, Li-, Ne-like Fe, mu=1e-4/1e-3, with/without QED+Breit: He/Li within 0.02%,
   Ne 0.06% of first order; H-like shift ~0 (0.00009 vs 0.0708 eV) -> PotentialHX zeroes ZPS for N<1
   (read in source) [RAW diagnostics/first_order].
 * debye_mode=0: NO calculation was run with debye_mode=0 in this session. Only source inspection.
 * "mu = 1e-6 control": I did not run mu=1e-6. The value 1e-6 in this project is the plasma density
   nps = 1e-6 A^-3 (C2 and production), not mu. Smallest mu run: 1e-4.


6. macOS PFAC CRASH INVESTIGATION
---------------------------------
I did not perform any macOS work in this session, and found no evidence of it on this host:
 * `grep -rIl -E 'SIGBUS|EXC_BAD_ACCESS|jmaxk.*uninit'` over /data/project/fac and the session
   temp dir: no matches. The only macOS traces are .DS_Store files committed in c4a016a by the user.
 * I never ran PFAC, never imported pfac, never inspected a macOS Python environment or PFAC
   extension path, never saw a backtrace with SIGBUS/EXC_BAD_ACCESS, and did not inspect or claim
   anything about an uninitialized `jmaxk` in PotentialHX1. (I read PotentialHX1 while studying the
   Linux failure; I made no finding about jmaxk.)
 * PFAC scripts present in the repo (user's, not run by me): test_debye_screening_minimal.py,
   fe_calc/fe_iii_diag.py, fe_calc/screened_ip_potential_feii.py, fe_calc/fe_II_energy.py,
   fe_calc/fe_III_energy.py, fe_calc/ip_potential_feii.py. The three I grepped use
   debye_mode=1 with ee_screen=0 and ConfigEnergy(0)/(1).
 * "macOS ref results" and the "previously recorded macOS branch fix-debye-ee-screening" cannot be
   verified from anything available to me. MISSING: crash logs, backtraces, macOS build provenance.


7. LINUX SFAC SWEEP AND RUNNING JOBS
------------------------------------
Driver: /data/project/fac/fe_calc/deprince_t2/deprince_t2_sfac.py. Standard grid: 6 models x 8 ions
 (Fe II..Fe IX) x 3 mu = 144 runs (the user's "48 cases" = sdp+gsdp x 8 ions x 3 mu).
Protocol per run (one sfac process, own directory runs/<model>/<ion>/<mutag>/): Config groups
 'ref' (and 'g0','sd'), cold OptimizeRadial on ['ref'] (ref,sdv,sdp) or ['g0'] (gref,gsdv,gsdp) for
 the free ion; for mu>0: SetOption orbital:debye_mode 1, radial:ee_screen 1, radial:hxs_iter0 8,
 then PlasmaScreen(zps=26-N, 1e-6, T(mu_k), 1, 0.0, 1) + OptimizeRadial for mu_k = 0.05,0.10,...,mu;
 GetPotential('final.pot'); for SD models Structure(parity,[2J_ground]); Structure('e.lev.b', all
 groups); MemENTable; PrintTable(...,1). No ConfigEnergy in any production input.
Models: ref = Table-1 refs, orbitals on statistical-weight average; sdv = + SD of ground ref in
 {3s,3p,3d,4s}; sdp = + SD of every ref into {1s..3d,4s} (paper's active space), ground (P,2J)
 block; gref/gsdv/gsdp = same CI but orbitals optimized on the ground configuration only.
 Fe IX refs (not in paper): [Ar], [3p]3d, [2p]3d, [1s]3d (my choice).
Status at 12:52 UTC (from run.json, all 144 files): 144 status ok, 144 binary 3a77793e,
 144 git_commit ff9a5c9, 144 source_diff_sha256 ceac162e... [RAW/CSV].
 History: sdp FeVI (3) stopped by me at 09:12 (recorded failed, exit -15), rerun ok at 12:48
 (E0 -34398.4667, -33499.2228, -32206.4715 eV; NLEV 10591; ground 2J 3; dps 10.0/4.0 a0) [RAW].
 Earlier non-final runs (binary 3b88ec91) were superseded by reruns with 3a77793e in the same dirs.
Validation performed on all completed runs by audit.py [SUB]: NELE correct; E0 header equals the
 minimum-energy level (relative energy 0, same index); for mu>0: debye_mode 1, ee_screen 1,
 lambda from last PlasmaScreen = 1/mu to 1e-9, final.pot mps=1, dps*mu=1 to 1e-6, zps=net charge;
 mu=0: no PlasmaScreen/screening options; no ConfigEnergy; adjacent-ion pairs share binary hash.
 Result: "consistency issues: none". Metadata (run.json/CSV) vs raw: E0 values agree.
NOT checked: per-run convergence messages in sfac.log beyond exit code; full level tables beyond
 the ground row; results of analyze.py (never run); REPORT.md results section (not written).


8. IP AUDIT AND NUMERICAL DISCREPANCIES
---------------------------------------
Provenance of the table below: scratchpad/audit/ip_audit.csv, written 12:48 UTC by audit.py, which
 parses e.lev (E0 header and level rows), in.sf and final.pot directly and computes
 IP = E0(next ion) - E0(ion) with same model, same mu, same binary hash [SUB]. These are CURRENT
 values. They agree with the driver's results/ip_comparison.csv [CSV] (spot checks [SESS]).
 Historical values: ref/sdv/gref/gsdv first computed with binary 3b88ec91 (08:08-08:10) were identical
 to 4 decimals to the 3a77793e reruns where compared [SESS]. The "19.72 eV sdp Fe III" figure mentioned
 by the user corresponds to current sdp Fe III mu=0 = 19.7215 eV (sdv: 19.7269).

ref, sdp, gsdp (all six models are in the CSV):
model  ion     mu       E0_ion(eV)   E0_next(eV)     IP_FAC   IP_pub   err(eV)   rel(%)
ref    FeII    0.0     -34561.2842   -34542.4125    18.8717    18.89    -0.018    -0.10
ref    FeII    0.1     -33625.2303   -33611.0422    14.1881    14.20    -0.012    -0.08
ref    FeII    0.25    -32282.8180   -32274.1005     8.7175     7.78    +0.938   +12.05
ref    FeIII   0.0     -34542.4125   -34522.7607    19.6518    32.93   -13.278   -40.32
ref    FeIII   0.1     -33611.0422   -33599.8850    11.1572    25.60   -14.443   -56.42
ref    FeIII   0.25    -32274.1005   -32274.0752     0.0253    16.42   -16.395   -99.85
ref    FeIV    0.0     -34522.7607   -34469.8762    52.8845    55.37    -2.486    -4.49
ref    FeIV    0.1     -33599.8850   -33557.4292    42.4558    44.65    -2.194    -4.91
ref    FeIV    0.25    -32274.0752   -32246.4036    27.6716    29.64    -1.968    -6.64
ref    FeV     0.0     -34469.8762   -34396.1570    73.7192    74.92    -1.201    -1.60
ref    FeV     0.1     -33557.4292   -33496.7988    60.6304    61.56    -0.930    -1.51
ref    FeV     0.25    -32246.4036   -32204.1556    42.2480    42.83    -0.582    -1.36
ref    FeVI    0.0     -34396.1570   -34298.5258    97.6312    98.20    -0.569    -0.58
ref    FeVI    0.1     -33496.7988   -33414.9389    81.8599    82.20    -0.340    -0.41
ref    FeVI    0.25    -32204.1556   -32144.3481    59.8075    59.79    +0.018    +0.03
ref    FeVII   0.0     -34298.5258   -34174.8058   123.7200   123.60    +0.120    +0.10
ref    FeVII   0.1     -33414.9389   -33309.6538   105.2851   104.96    +0.325    +0.31
ref    FeVII   0.25    -32144.3481   -32064.8020    79.5461    78.86    +0.686    +0.87
ref    FeVIII  0.0     -34174.8058   -34025.0667   149.7391   151.14    -1.401    -0.93
ref    FeVIII  0.1     -33309.6538   -33181.0112   128.6426   129.85    -1.207    -0.93
ref    FeVIII  0.25    -32064.8020   -31965.5683    99.2337    99.98    -0.746    -0.75
sdp    FeII    0.0     -34561.8278   -34543.3083    18.5195    18.89    -0.370    -1.96
sdp    FeII    0.1     -33625.7892   -33612.8595    12.9297    14.20    -1.270    -8.95
sdp    FeII    0.25    -32285.2509   -32278.5717     6.6792     7.78    -1.101   -14.15
sdp    FeIII   0.0     -34543.3083   -34523.5868    19.7215    32.93   -13.209   -40.11
sdp    FeIII   0.1     -33612.8595   -33600.7790    12.0805    25.60   -13.520   -52.81
sdp    FeIII   0.25    -32278.5717   -32274.9141     3.6576    16.42   -12.762   -77.72
sdp    FeIV    0.0     -34523.5868   -34471.4499    52.1369    55.37    -3.233    -5.84
sdp    FeIV    0.1     -33600.7790   -33559.0963    41.6827    44.65    -2.967    -6.65
sdp    FeIV    0.25    -32274.9141   -32247.9863    26.9278    29.64    -2.712    -9.15
sdp    FeV     0.0     -34471.4499   -34398.4667    72.9832    74.92    -1.937    -2.59
sdp    FeV     0.1     -33559.0963   -33499.2228    59.8735    61.56    -1.687    -2.74
sdp    FeV     0.25    -32247.9863   -32206.4715    41.5148    42.83    -1.315    -3.07
sdp    FeVI    0.0     -34398.4667   -34301.7568    96.7099    98.20    -1.490    -1.52
sdp    FeVI    0.1     -33499.2228   -33418.3007    80.9221    82.20    -1.278    -1.55
sdp    FeVI    0.25    -32206.4715   -32147.5734    58.8981    59.79    -0.892    -1.49
sdp    FeVII   0.0     -34301.7568   -34179.1619   122.5949   123.60    -1.005    -0.81
sdp    FeVII   0.1     -33418.3007   -33314.1557   104.1450   104.96    -0.815    -0.78
sdp    FeVII   0.25    -32147.5734   -32069.1342    78.4392    78.86    -0.421    -0.53
sdp    FeVIII  0.0     -34179.1619   -34030.4868   148.6751   151.14    -2.465    -1.63
sdp    FeVIII  0.1     -33314.1557   -33186.5904   127.5653   129.85    -2.285    -1.76
sdp    FeVIII  0.25    -32069.1342   -31970.9483    98.1859    99.98    -1.794    -1.79
gsdp   FeII    0.0     -34571.6567   -34556.1514    15.5053    18.89    -3.385   -17.92
gsdp   FeII    0.1     -33636.2088   -33626.0115    10.1973    14.20    -4.003   -28.19
gsdp   FeII    0.25    -32294.1224   -32289.4642     4.6582     7.78    -3.122   -40.13
gsdp   FeIII   0.0     -34556.1514   -34528.3851    27.7663    32.93    -5.164   -15.68
gsdp   FeIII   0.1     -33626.0115   -33606.2100    19.8015    25.60    -5.798   -22.65
gsdp   FeIII   0.25    -32289.4642   -32280.5042     8.9600    16.42    -7.460   -45.43
gsdp   FeIV    0.0     -34528.3851   -34475.7982    52.5869    55.37    -2.783    -5.03
gsdp   FeIV    0.1     -33606.2100   -33564.0491    42.1609    44.65    -2.489    -5.57
gsdp   FeIV    0.25    -32280.5042   -32253.0657    27.4385    29.64    -2.201    -7.43
gsdp   FeV     0.0     -34475.7982   -34402.5863    73.2119    74.92    -1.708    -2.28
gsdp   FeV     0.1     -33564.0491   -33503.9514    60.0977    61.56    -1.462    -2.38
gsdp   FeV     0.25    -32253.0657   -32211.3060    41.7597    42.83    -1.070    -2.50
gsdp   FeVI    0.0     -34402.5863   -34305.5258    97.0605    98.20    -1.139    -1.16
gsdp   FeVI    0.1     -33503.9514   -33422.6565    81.2949    82.20    -0.905    -1.10
gsdp   FeVI    0.25    -32211.3060   -32152.0160    59.2900    59.79    -0.500    -0.84
gsdp   FeVII   0.0     -34305.5258   -34182.6777   122.8481   123.60    -0.752    -0.61
gsdp   FeVII   0.1     -33422.6565   -33318.2265   104.4300   104.96    -0.530    -0.50
gsdp   FeVII   0.25    -32152.0160   -32073.2771    78.7389    78.86    -0.121    -0.15
gsdp   FeVIII  0.0     -34182.6777   -34033.8199   148.8578   151.14    -2.282    -1.51
gsdp   FeVIII  0.1     -33318.2265   -33190.4681   127.7584   129.85    -2.092    -1.61
gsdp   FeVIII  0.25    -32073.2771   -31974.8859    98.3912    99.98    -1.589    -1.59

Derived observations (from the same raw energies) [SUB]:
 * Orbital-averaging penalty P = E0(ref) - E0(gref) at mu=0: FeII 9.990, FeIII 13.342, FeIV 4.899,
   FeV 4.450, FeVI 4.214, FeVII 3.850, FeVIII 3.586, FeIX 3.390 eV. IP(ref)-IP(gref) for Fe III =
   P(FeIV)-P(FeIII) = -8.442 eV (mu=0). Statistical-weight composition of the ref average: Fe III
   ground weight 1.9%, <4s>=1.14 while ground 3d6 has no 4s (only ion with this mismatch).
 * sdp-gsdp for FeVI: 4.120, 4.729, 4.835 eV (mu 0, 0.1, 0.25).
 * Screening-only part dIP = IP(mu)-IP(0) for Fe IV-VIII is model-independent and 1-3% smaller in
   magnitude than Deprince (e.g. FeVIII mu=0.25: -50.47..-50.51 vs -51.16). For Fe II/III dIP is
   model dependent (Fe III mu=0.25: -16.06 sdp, -18.81 gsdp, -19.63 ref, vs -16.51).
 * Monotonic decrease IP(0)>IP(0.1)>IP(0.25) holds for every ion and model (results/summary.md).


9. GIT, BUILD, AND EXECUTABLE PROVENANCE
----------------------------------------
Main repo /data/project/fac: branch fix-debye-ee-screening @ c4a016ad8cf0458366714b013727d915e1fdec80;
 tracked files unmodified by me; untracked additions include fe_calc/deprince_t2/.
 Local branch master = d40ba31 (STALE; origin/master = 5f8c8d2). ee-screening-scf = ff9a5c9cc6623e3f6
 adb8062ba23a52e9aace582 (verified). Stash stash@{0} pre-existing.
Worktrees (git metadata in /data/project/fac/.git/worktrees):
 /data/project/fac_worktrees/scf: HEAD ff9a5c9 (detached), modified: faclib/mmalloc.c, orbital.c,
   radial.c, structure.c (56+/11-). `git diff -- faclib sfac | sha256sum` = ceac162ead6ee5d5159e8ed02d
   1a4577791be1cd987e1ce2585fa49a0267d7a0 (VERIFIED at 12:53; equals patches/scf_build_working_tree.diff
   and the value recorded in all 144 run.json).
 /data/project/fac_worktrees/ee: HEAD c4a016a (detached), modified orbital.c, radial.c (40+/1-);
   diff sha256 ef42820f... = patches/debye_sfac_fixes.patch.
Executables (sha256 prefix, mtime UTC):
 ~/facinst/bin/sfac                       2789ca41138defd1  2026-09-14 04:43  (not built by me)
 scf/inst_unpatched/bin/sfac              2789ca41138defd1  2026-10-12 07:46  clean ff9a5c9 build
 ee/inst_unpatched/bin/sfac               279d393eb7150854  2026-10-12 07:46  clean c4a016a build
 ee/inst/bin/sfac                         23b5644a22fbb653  2026-10-12 08:02  c4a016a + patch
 scf/inst_fix_no842/bin/sfac              3b88ec91c3757787  2026-10-12 08:02  ff9a5c9 + patch (no #842)
 scf/inst/bin/sfac  (FINAL, all 144 runs) 3a77793e4b5e6f2c  2026-10-12 08:21  ff9a5c9 + ceac162e diff
 Exact source for 3a77793e: ff9a5c9 + diff with sha256 ceac162e (sources' mtimes 08:01-08:21 precede
 binary mtime 08:21:26). For 2789ca41 (installed): equals a clean build of ff9a5c9 bit-for-bit, so
 its effective source is known even though it was built (2026-09-14) from a then-uncommitted tree.
 Build compiler: environment default at build time (original installed build used a conda
 x86_64-conda-linux-gnu-cc per the repo's config.log); not OpenMP (no libgomp linked).
Active workers: none at audit time. During the sweep, /proc/<pid>/exe of the FeVI reruns was
 /data/project/fac_worktrees/scf/inst/bin/sfac (checked at ~12:13 and ~12:16, from recorded elapsed times 52:00 and 54:42 [SESS]).
Mode-0 sign provenance [SRC]: upstream c3be2b4 (#840, 2026-09-04) changed debye_mode 0 to
 ZPS = -zps*(x-1). fix-debye-ee-screening and origin/master contain it; ee-screening-scf (and hence
 INST and FINAL binaries) and local master do not (ZPS = zps*(x-1)). debye_mode 1 code is identical
 in both trees.
macOS provenance: NOT AVAILABLE to me (no files, no logs).


10. VERIFIED FINDINGS
---------------------
V1 Installed sfac = ee-screening-scf build, not fix-debye-ee-screening (sha256 equality).
V2 Original failure cause 1: ion-sphere truncation (ips) in DensityToSZ with ups=0 -> repulsive tail
   (C1 vs C2 potentials).
V3 Cause 2: RadialBound phase-a upper bracket cannot exist when the next same-kappa state is unbound;
   turning-point node counts give false negatives near E=0 (v1 vs v2 behaviour).
V4 Cause 3: cold-start SCF transient and exchange-free iteration 0 unbind 3d at mu>=~0.09 for Fe II;
   free-start warm ladder with hxs_iter0=8 converges (D1-D6).
V5 Patched final binary reproduces unpatched free-ion energies exactly (E1/E2) and the original
   baseline (E3, 18.8716/18.8717 eV).
V6 lambda_D = 1/mu internally (dps in final.pot), zps and nps irrelevant in debye_mode 1 (E4/E5).
V7 Int-overflow segfault in AllocHamMem for >~8192-CSF blocks on ee-screening-scf; fixed by d3a4ec3
   backport.
V8 All 144 production runs ok and internally consistent (audit, section 7).
V9 Fe III ref/sdv/sdp IP anomaly (~19.7 eV) quantitatively equals the orbital-averaging penalty
   difference (section 8).
V10 debye_mode 0 sign differs between ee-screening-scf-derived binaries and upstream (source only;
   not exercised in any run).
V11 One-electron ions lose e-n Debye screening (first_order H-like test + source).
V12 Grid (<=0.023 eV) and warm-restart (<=0.011 eV) sensitivities of IPs are small.


11. HYPOTHESES AND UNRESOLVED QUESTIONS
---------------------------------------
H1 [INF] Remaining 1-3 eV deviations from Deprince for Fe IV-VIII and the ~1-3% smaller screening
   lowering are due to FAC's single local-potential orbitals vs MCDHF-EAL and missing correlation;
   not tested.
H2 [INF] Fe III penalty is driven by the 4s-bearing core-hole references (<4s>=1.14 vs 0). The
   one-factor test (optimize Fe III on refs without 4s-bearing configs) was PROPOSED, NOT RUN.
H3 [INF] Small-mu IP-slope deviations (ref Fe II/III) reflect first-order orbital response of
   non-variational orbitals; not proven.
H4 Why the published CSF counts are 1.23-1.42x smaller than my reconstruction of the stated active
   space (results/csf_counts.csv): unknown.
H5 hxs_iter0 history dependence: 0.30 eV free-ion E0 change for cold starts with value 20 is not
   understood mechanistically.
H6 The macOS crash (if any) is entirely outside my evidence.


12. CURRENT PROJECT STATE
-------------------------
 * 144/144 runs ok; results/{energies.csv,ip_comparison.csv,summary.md} regenerated 12:48.
 * analyze.py exists but has NEVER been run; results/ip_lowering.csv and model_overview.md do NOT exist.
 * REPORT.md incomplete (RESULTS_PLACEHOLDER at line ~317; written before the paper-AS results and
   before the read-only audit findings, e.g. it does not mention the debye_mode-0 sign issue).
 * Worktrees contain uncommitted patches (sections 3.1-3.4); no commits or branches were created.
 * throttle_sfac.py (pid 2662028) still alive at 12:53; will exit by itself.
 * No sfac running.


13. KNOWN RISKS AND PRESERVATION CONSTRAINTS
--------------------------------------------
 * Rerunning the driver with a different binary will rerun cases (reuse requires matching sha256) and
   DELETE the existing outputs of those cases before running them.
 * Rebuilding in fac_worktrees/scf changes the binary hash -> all 144 runs become "incompatible" pairs.
 * debye_mode=0 with INST or FINAL binaries uses the pre-#840 sign.
 * sdp/gsdp FeIV/VI/VII need ~2-3.5 GB RAM each and 45-90 min; >3 in parallel risks OOM on this host.
 * git worktree removal would delete the only copies of the build trees; patch files are in
   fe_calc/deprince_t2/patches/.
 * The session scratchpad (/tmp/...) is session-specific and may not persist.


14. EVIDENCE INDEX
------------------
 Project:        /data/project/fac/fe_calc/deprince_t2/  (driver, REPORT.md, provenance.json, patches/,
                 logs_*.txt, throttle.log, results/, runs/, diagnostics/)
 Runs:           runs/<model>/<ion>/<mutag>/{in.sf,sfac.log,e.lev,e.lev.b,final.pot,run.json}
 Diagnostics:    diagnostics/diagnostics.{md,json}, diagnostics/<case>/, small_mu_check.txt,
                 warm_restart_check.txt, grid/grid_sensitivity.txt, branch_ee/, first_order/,
                 scratch_nps/ (including iter/, ladder*/, pathdep/, muscan/, potcmp/, v3/, hxs*/),
                 scratch_feiii/, scratch_segv/
 CSF counts:     results/csf_counts.csv
 Builds:         /data/project/fac_worktrees/{ee,scf}/ (build.log, build_patched.log, build_842.log,
                 inst*/), build_status.txt, scf_final_working_tree.diff, upstream_d3a4ec3_int_overflow.patch
 Original work:  /data/project/fac/fe_calc/fe_ip_sfac.py, fe_calc/fe_ip_sfac_runs/ (untouched)
 Paper text:     scratchpad/deprince.txt
 Audit:          scratchpad/audit/{audit.py,raw.json,ip_audit.csv}
 Task logs:      /tmp/claude-1000/-data-project-fac/a8e208a0-910c-4f2e-b650-3503c715a8c5/tasks/*.output
 Transcript:     /home/sadasivanmj/.claude/projects/-data-project-fac/a8e208a0-910c-4f2e-b650-3503c715a8c5.jsonl

Information requested but UNAVAILABLE: fac_debyescreening.zip contents; any macOS crash data,
 backtraces, PFAC extension path or macOS environment; exact text of RadialBound patch v1 and
 hxs_iter0 v1 as files (only in the transcript); exact start times of some steps (approximate);
 a mu=1e-6 calculation (never run).


TOP 5 UNRESOLVED QUESTIONS (ranked)
-----------------------------------
1. Which orbital-optimization protocol should define the FAC result for comparison with Deprince
   (EAL-like ref average vs ground-optimized), given the ion-dependent penalty that drives the
   Fe III error? (H2 test not yet run.)
2. Is the remaining 1-3% systematic difference in the screening-induced IP lowering (Fe IV-VIII) a
   method difference (local potential vs MCDHF) or an implementation issue?
3. Does the macOS crash (outside my evidence) share a cause with any Linux finding?
4. What active-space restrictions produce the paper's published CSF counts?
5. What is the origin of the SCF history dependence (0.30 eV with hxs_iter0=20 on cold starts;
   0.01 eV on warm restarts), and does it bias screened vs free comparisons?
