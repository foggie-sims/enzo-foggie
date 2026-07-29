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
CI.1-CI.5.

In flight:
- **T1.8 (instrumentation)** - code landed (`RH_PERF` in
  RebuildHierarchy, DtLimiter histogram in SetLevelTimeStep /
  Grid_ComputeTimeStep, MGSolver stats in MultigridSolver, per-root-cycle
  reduction+print in EvolveHierarchy). CI green (run 14). **Missing only
  the production A/B stamp** - that is the t18-instrumentation bench
  (see below). When it passes: flip T1.8 to done and report the first
  production RHperf/DtLimiter/MGSolver numbers from the candidate's
  `enzo_bench.log`.
- **T1.6** - RNG seeding in star_feedback6.F is done; allocatable
  arrays + early return still to do.
- **CI.6** - the pull-model cron runner; superseded by this migration
  (see below).

Next queue after the T1.8 stamp (user-approved order): T1.1
(StarParticleFindAll guard) -> T1.2/T1.3 -> T1.13 (sanitizer-gated) ->
finish T1.6 -> T1.5.

## Key technical facts (hard-won, do not re-learn)

- **T1.14 / noise floor**: parallel Enzo is NOT run-to-run bitwise
  reproducible (message-arrival order changes FP accumulation;
  CONFIG_BITWISE_IDENTICALITY covers only gravity/photon paths). At the
  production anchor (512 ranks, 5 root steps): relative mass diffs grow
  ~1 decade per root step (1e-11 -> 1e-7 by step 5), refinement structure
  diverges by step 4, +-1 stochastic star event, ~2.2% wall-clock noise.
  Therefore all production A/B gating MUST use
  `compare_runs.py --noise-floor <archived floor json>` (gates become
  max(class tol, 10x floor)). Absolute tolerances only suit serial runs.
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
`8aa475f` (T1.8, CI green). 5 root steps, 512 ranks, rom_ait/128,
devel, class C0 with the archived noise floor. Attempt history (full
records under `results/t18-instrumentation-r*/` on bench-results):
r1 stale grackle path in the committed machine file; r2 runner
self-update lag; r3 placeholder restart path; r4 qsub not in cron PATH;
r5 bench_A submitted (job 24914598) but bench_B bounced off the
per-user queue limit and the entry failed. r6 was queued for the cron
runner just before this migration.

**Salvage path for the CLI session** (cheapest route to the stamp):
`~/foggie_bench_root/runs/t18-instrumentation-r5/` has both run dirs
fully prepared. If bench_A (baseline) completed (check
`bench_A/bench_exit_status` == exit=0), just `qsub bench_B/launch.sh`
when a devel slot is free, and when B finishes run:

    python3 run/foggie_bench/compare_runs.py \
        ~/foggie_bench_root/runs/t18-instrumentation-r5/bench_A \
        ~/foggie_bench_root/runs/t18-instrumentation-r5/bench_B \
        --class C0 \
        --noise-floor ~/foggie_bench_root/results/results/t19-manual/noise_floor_A1_vs_A2.json

(Adjust the results-clone path as found on disk.) Then: archive the
comparison + candidate `enzo_bench.log` instrumentation output to
`bench-results` under `results/t18-instrumentation-r5/`, flip T1.8 to
done in the ledger, and report the first RHperf/DtLimiter/MGSolver
numbers.

## The cron runner is decommissioned by this migration

Remove the crontab line on the front end (`crontab -e`, delete the
`hpc_runner_cron.sh` entry) so the runner and the interactive session
never race. `hpc_runner.py` remains usable as a manual one-tick script
(`sh run/foggie_bench/hpc_runner_cron.sh` = one tick) if batch mode is
ever wanted again; `queue.yml` entry `t18-instrumentation-r6` is moot
once r5 is salvaged manually - remove it or leave it, nothing consumes
it with cron off. Mark CI.6 accordingly (already noted in the ledger).

## Aitken-local layout

    ~/foggie_bench_root/repo      clone of enzo-foggie (enzo-performance),
                                  pushes via SSH deploy key
    ~/foggie_bench_root/results   clone tracking bench-results
    ~/foggie_bench_root/builds/   per-sha cached builds
    ~/foggie_bench_root/runs/     bench run dirs (r5 is the live one)
    ~/.ssh/config                 Host github.com -> ssh.github.com:443
                                  with the foggie_bench_deploy key
    ~/foggie_bench_cron.log       runner log (historical once cron is off)
