# Audit handoff: cloud session -> Claude Code CLI on Aitken

Written 2026-07-29 when the Enzo performance-audit effort moved from a
Claude Code cloud session to a CLI session running directly on an Aitken
front end. This file is the memory transfer: everything a fresh session
needs that is not obvious from the code itself. Read `audit/status.yml`
(the ledger) alongside it; the two together are the project state.

## What this effort is

A performance and parallel-scaling audit of Enzo in the FOGGIE regime
(deep AMR zooms, 10+ levels, ~1.4M concentrated star particles), followed
by staged remediation. Key documents, all in `audit/`:

- `PERFORMANCE_AUDIT.md` - the audit itself (tiers T0-T3 of findings).
- `enzotestplan.html` / `enzoaudit.html` - HTML renders.
- `status.yml` - the remediation ledger. **Update it with every fix.**
- `build_dashboard.py` - renders the ledger to `dashboard.html`.
  `dashboard.html` is committed only by the GitHub Actions bot
  (`audit-dashboard.yml`) - do not commit it by hand, it causes rebase
  conflicts.

## Standing policies (from the user, still in force)

1. **No PRs to main until further notice** (2026-07-28). All work lands
   directly on the `enzo-performance` branch.
2. Every audit fix updates `audit/status.yml` in the same push.
3. Commits authored as Claude use `noreply@anthropic.com`.
4. The bench_0 production-reference comparison is **iced** - the user is
   re-running bench0 with a different Enzo branch. Do not resume it until
   they ask. When revived: archive under `results/production-ref-*` on
   `bench-results` and add a yardstick entry to the ledger.

## Branches

- `enzo-performance` - the working branch: audit docs, CI workflows,
  bench harness, and all code fixes.
- `bench-results` - one-way results archive written from Aitken
  (`results/<id>/state.json`, comparison verdicts, log tails). The
  measured production noise floor lives at
  `results/t19-manual/noise_floor_A1_vs_A2.json`. Never merge these
  branches into each other.
- CI (`.github/workflows/enzo-performance-ci.yml`) is scoped to pushes /
  PRs on `enzo-performance` only, with `paths-ignore` for `audit/**`,
  `run/foggie_bench/**`, and `**.md`. It builds Enzo, runs smoke decks
  (SedovBlastAMR np4, GravityTest, TestOrbit, TabularFeedbackTest with an
  SN-actually-fired assertion), a gravity force-accuracy gate pinned to
  the measured baseline (RMS/max/tan = 0.15/0.70/0.05), serial bitwise
  determinism pairs (hydro and star-feedback), and an ASan/UBSan job.
  It keeps running on every push - nothing to migrate.

## Ledger state in brief (authoritative copy: status.yml)

Done: T0.8, T1.9 (RebuildHierarchy + FastSiblingLocator O(S^2) removal -
validated bitwise in CI and within 10x noise floor at production scale),
T1.8 (instrumentation - production-stamped 2026-07-29 by
t18-instrumentation-r7 at the corrected 1x128 mil_ait config; first
production numbers in the ledger note and the r7 state.json), CI.1-CI.6
(CI.6 closed as superseded by the CLI-session migration).

In flight:
- **T1.6** - RNG seeding in star_feedback6.F is done; allocatable
  arrays + early return still to do.

Next queue (user-approved order): T1.1 (StarParticleFindAll guard) ->
T1.2/T1.3 -> T1.13 (sanitizer-gated) -> finish T1.6 -> T1.5.

Open question for the user: the T1.9 production stamp was measured at
the old 512-rank rom_ait config; its serial-bitwise CI evidence stands,
but re-stamping at 1x128 mil_ait is cheap now that the envelope floor
exists, if wanted.

## Key technical facts (hard-won, do not re-learn)

- **T1.14 / noise floor**: parallel Enzo is NOT run-to-run bitwise
  reproducible (message-arrival order changes FP accumulation;
  CONFIG_BITWISE_IDENTICALITY covers only gravity/photon paths). Measured
  at the anchor at 512 ranks / rom_ait (a superseded config - production
  is 1 node x 128 ranks mil_ait; the r7 A1/A2 pair re-measures the floor
  there): relative mass diffs grow ~1 decade per root step (1e-11 -> 1e-7
  by step 5), refinement structure diverges by step 4, +-1 stochastic
  star event, ~2.2% wall-clock noise. Therefore all production A/B gating
  MUST use `compare_runs.py --noise-floor <archived floor json>` (gates
  become max(class tol, 10x floor)), with a floor measured at the same
  rank count / node model. Absolute tolerances only suit serial runs.
- The anchor snapshot and its properties: `run/foggie_bench/anchor.md`.
- Bench-run preparation gotchas are all encoded in
  `run/foggie_bench/make_bench_run.py`: restart syntax is
  `-r RD0016/RD0016` (snapshot is a directory), shadow-dir with symlinks
  and an edited parameter file, `NumberOfOutputsBeforeExit = 0` forced
  (else Enzo exits at the first output), `InitialCycleNumber` (not
  CycleNumber), devel walltime clamps to 2:00:00, group_list s3128.
- Toolchain on Aitken (also in `Make.mach.nasa-aitken-milan-mpich` and
  `pbs_pleiades.template.sh`): comp-intel/2020.4.304,
  hdf5/1.8.18_serial, MPICH 4.0.3 at
  `/u/jtumlins/installs/mpich-4.0.3/usr/local` (launcher MUST be its
  mpiexec - MPT's mpiexec fails with sgicheckppversion), grackle at
  `/u/jtumlins/grackle/grackle-3.3.1-dev/build` (lib64), anaconda at
  `/home1/jtumlins/anaconda3`. System python is 3.6: no `text=` /
  `capture_output=` subprocess kwargs.
- The devel queue has a 2 h walltime cap and a small per-user job limit -
  submitting baseline and candidate simultaneously can bounce with
  "would exceed queue generic's per-user limit". Submit sequentially or
  absorb the rejection and retry.
- PBS commands are not in cron/bare-shell PATH; interactive shells get
  them from profile scripts. (`/PBS/bin`.)
- Builds: `make machine-nasa-aitken-milan-mpich && make grackle-yes &&
  make opt-high && make -j8`. Cached per-sha builds live in
  `~/foggie_bench_root/builds/<sha12>/`. When building a historical ref,
  overlay the tip's `Make.mach.nasa-aitken-milan-mpich` into the old tree
  first (old refs predate the machine file / have stale paths).
  Cosmetic known issue: that machine file's MACH_FILE line says
  nasa-aitken-rome.

## The t18-instrumentation bench (T1.8 production stamp) - current state

Baseline = `7f22798^` (= 42dcbbd, tip before T1.8); candidate =
`8aa475f` (T1.8, CI green). 5 root steps from the anchor restart,
devel queue, class C0. Attempt history (full records under
`results/t18-instrumentation-r*/` on bench-results): r1 stale grackle
path in the committed machine file; r2 runner self-update lag; r3
placeholder restart path; r4 qsub not in cron PATH; r5 bench_A (job
24914598) died in 24 s writing its first dump - "Disk quota exceeded",
because foggie_bench_root then lived on /home1 proper (since moved to
/home1/jtumlins/nobackup). r6 was queued for the cron runner but never
ran (cron decommissioned).

r5 also exposed a config error carried since the first benches: runs
were prepared at 512 ranks on 4 x rom_ait (the make_bench_run.py
defaults), but FOGGIE production runs on **1 Milan node, 128 ranks
(mil_ait)**. The defaults now match production, and because a noise
floor is config-specific, the 512-rank t19-manual floor cannot gate
128-rank runs.

**Resolved by t18-instrumentation-r7** (2026-07-29, CLI session):
three baseline runs (A1/A2/A3) plus candidate B at 1 x 128 mil_ait.
A single-pair floor proved statistically fragile - identical-code
step-5 baryon diffs scattered 5.8e-9 to 2.9e-8 between pairs, and the
first (unluckily tight) pair tripped a false FAIL on the candidate -
so the floor is now the all-pairs envelope built with
`merge_noise_floors.py` (new harness script). Against it the T1.8
candidate PASSED class C0 with 3x margin at the worst gate.
Archived under `results/t18-instrumentation-r7/` on bench-results
(state.json there has the full story + the first production
RHperf/DtLimiter/MGSolver numbers). T1.8 is done in the ledger.
Lesson recorded in anchor.md: always measure floors from 3+ baseline
runs and gate against the envelope; bench jobs may use the normal
queue (user approval 2026-07-29), which allows concurrent runs.

## The cron runner is decommissioned by this migration

Remove the crontab line on the front end (`crontab -e`, delete the
`hpc_runner_cron.sh` entry) so the runner and the interactive session
never race. `hpc_runner.py` remains usable as a manual one-tick script
(`sh run/foggie_bench/hpc_runner_cron.sh` = one tick) if batch mode is
ever wanted again; `queue.yml` entry `t18-instrumentation-r6` is moot
once r5 is salvaged manually - remove it or leave it, nothing consumes
it with cron off. Mark CI.6 accordingly (already noted in the ledger).

## Aitken-local layout

Root moved 2026-07-29 from /home1/jtumlins/foggie_bench_root to
`/home1/jtumlins/nobackup/foggie_bench_root` (= /nobackupnfs1/jtumlins/
foggie_bench_root) after a bench dump blew the home quota (r5 failure).

    .../foggie_bench_root/repo      clone of enzo-foggie (enzo-performance),
                                    pushes via SSH deploy key
    .../foggie_bench_root/results   clone tracking bench-results
    .../foggie_bench_root/builds/   per-sha cached builds
    .../foggie_bench_root/runs/     bench run dirs (r7 is the live one)
    ~/.ssh/config                   Host github.com -> ssh.github.com:443
                                    with the foggie_bench_deploy key
    ~/foggie_bench_cron.log         runner log (historical; cron is off)
