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

Done: T0.1, T0.2, T0.6, T0.8, T1.1, T1.2, T1.3, T1.5, T1.6, T1.8, T1.9,
T1.13, CI.1-CI.6. T1.11 is wontfix (superseded - see below).

The user-approved Tier-1 queue and three Tier-0 items completed
2026-07-30 (this CLI session). Every code item carries three-level
evidence: local serial bitwise A/B (h5diff), branch-CI green (incl.
sanitizers + serial determinism), and a production bench stamp at the
1x128 mil_ait anchor gated by the envelope noise floor (archives under
`results/<id>/` on bench-results, each with a narrative state.json).
Wall-clock at the anchor: T1.5 is the standout (~7%, the eliminated
per-subcycle Grackle ComputeCoolingTime); the others are individually
at or below the ~1.6% node-to-node scatter.

Two Tier-0 items closed as measured negatives - both worth NOT
repeating:
- **T0.1**: production already runs RebuildHierarchyCycleSkip=2 on all
  levels (the audit's recommendation predates the audit). Going to
  skip=3 changes the answer far beyond the noise envelope for -0.6%
  wall. Keep skip=2.
- **T0.2**: LoadBalancing=4 now runs correctly at production scale
  (needed the T1.13 partitioner fixes + the key rework), passes C0,
  and is wall-neutral: non-rebuild phases ~2% faster, rebuild ~30 s
  slower (balancer + migration). Keep LoadBalancing=1 until T2.1 makes
  the work weights particle-aware.

In flight: instrumentation pass 3 (sha ed40e892) - CI + anchor bench
running as of this writing. See "Where the time actually goes" below.

## Where the time actually goes (measured, not estimated)

The audit's priorities were estimates. Since T1.8 the branch carries
real per-section timers, and `run/foggie_bench/time_budget.py <rundir>`
prints the attribution from any bench run's `performance.out`.

Instrumentation pass 2 (`results/instr2-r1` on bench-results, 5 root
steps at the anchor, 2273 s):

    ChemistryCooling      192.6 s   8.5%   <- Grackle; NOT the elephant
    RebuildHierarchy      178.3 s   7.8%
    SolveHydroEquations   111.1 s   4.9%
    SetBoundaryConditions 107.1 s   4.7%
    PrepareDensityField    79.9 s   3.5%
    SetLevelTimeStep       51.3 s   2.3%   (incl. the dt allreduce)
    SolveForPotential      51.2 s   2.3%
    Group_WriteAllData     50.7 s   2.2%
    StarParticleHandler    12.5 s   0.5%   <- after the Tier-1 batch
    ---------------------------------------
    attributed             37.2%
    unattributed           62.8%

    cross-cutting (inside the above, so they attribute but do not
    partition):
    CommReceiveHandler    100.1 s   4.4%
      of which pure wait   87.8 s   3.9%   (MPI_Waitsome)

Three findings that redirect the audit:
1. Chemistry is 8.5%, not the missing majority - the "Grackle
   elephant" hypothesis is falsified. Do not start cooling work on
   the strength of the audit text alone.
2. The receive pump is 88% pure MPI wait. Comm-flavoured sections are
   ~15% of wall as an upper bound.
3. The star subsystem, the audit's headline Tier-1 target, is now
   0.5% of wall. That work is done.

Passes 3 and 4 closed the accounting to **99.8%**, and the answer was a
surprise: **62.7% of wall is a single `MPI_Allreduce`** - the one inside
`CommunicationUpdateStarParticleCount`, called from
`StarParticleFinalize` once per subcycle of every level (4327 times in a
5-root-step run). It is ~100% barrier wait, not communication: the
payload is 1-8 kB among 128 ranks on ONE node against 320 ms measured
per call. Chemistry (8.5%), rebuild (8.0%) and hydro (5.0%) are the
next largest and are an order of magnitude smaller.

**Full analysis, proposed fix, and its three failure modes are in
`audit/SPFinalize_Edits.md`.** That work is deliberately DEFERRED behind
T2.1 (see that file's "Relationship to T2.1").

Two traps for anyone reading performance.out directly:
- `Level_NN` lines are NOT nested - the level timer brackets stop before
  the recursive EvolveLevel call, so they are exclusive per-level times
  summing to ~87% of Total. They are a different decomposition axis and
  must never be added to the section budget.
- Several sections are cross-cutting (`SolveForPotential` is started
  inside PrepareDensityField.C; `SPCommUpdateCount` inside SPFinalize;
  the comm pump inside several). `time_budget.py` knows which; coverage
  meaningfully over 100% means a new overlap needs classifying.

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
- **Header changes do NOT trigger recompiles on Aitken.** `makedepend`
  is not installed, so the Makefile's `dep` target fails silently
  ("Error 127 (ignored)") and `DEPEND` stays a 0-byte file - there are
  no header->object dependencies at all. Editing a `.C` file is safe
  (timestamp rule), but after touching ANY header you must
  `make clean` and rebuild, or you get stale objects. The failure mode
  is usually a link error (undefined reference to a new global), which
  is benign; the dangerous case is a changed struct layout, which links
  fine and is undefined behaviour at runtime. Discovered 2026-07-30
  when a new `EXTERN` in global_data.h failed to link because enzo.C
  (which defines the storage) was never rebuilt.
- **The compiler is not on PATH in a non-interactive shell.** A tool-run
  `make` dies with `icpc: command not found` on every object, which looks
  alarming but is only a missing module environment. Prefix any scripted
  build with:

      source /usr/share/modules/init/bash
      module use -a /nasa/modulefiles/testing
      module load comp-intel/2020.4.304 hdf5/1.8.18_serial

- **`enzo.exe` does not run on the login node**, and this is pre-existing,
  not a symptom of your change: MPICH's `libmpifort` wants
  `libgfortran.so.3`, which is present on compute nodes but not on the
  front end. `ldd enzo.exe | grep "not found"` reports it for any build,
  including known-good ones - check a previous executable before
  investigating. Smoke tests therefore have to go through PBS; the login
  node can only build and link.
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
