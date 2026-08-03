# StarParticleFinalize: the per-subcycle Allreduce barrier

> **CLOSED 2026-08-02 — DO NOT IMPLEMENT. Measured ceiling is 4.6%.**
>
> Before building the (intricate, corruption-prone) deferred-sync design
> below, we measured the upper bound directly: an env-gated throwaway build
> that skips the collective entirely, identical binary otherwise
> (`results/spf-measure-r1` on bench-results).
>
> | | wall | SPCommUpdateCount |
> |---|---|---|
> | keep | 1704.6 s | 832.1 s |
> | skip | 1626.2 s | 0.0 s |
>
> **832 s of barrier removed; 78 s of wall saved.** ~90% relocated:
> `MPIWaitReceive` +389 s, `SetBoundaryConditions` +360 s, the dt reduction
> +164 s. Real work sections moved by less than 3 s.
>
> **Why the analysis below was wrong.** It argued from `max(Σ) ≤ Σ max` that
> removing the per-subcycle barrier would let fluctuating imbalance cancel
> across subcycles. That requires this to be the *only* synchronisation
> point. It is not — `SetBoundaryConditions` exchanges ghost zones every
> subcycle and the dt reduction is global — so ranks re-synchronise anyway
> and the wait simply moves to the next barrier. The 62% was never "the cost
> of a collective"; it was *where real work imbalance happened to be paid*.
>
> **The transferable lesson:** only changes that reduce the imbalance (or the
> work) help. Moving where it is paid does nothing. This is why T0.3 worked
> (it removed indivisible oversized grids) and why particle-weighted
> balancing failed (it worsened the imbalance).
>
> The rest of this document is retained as the record of the analysis and of
> the code facts it established, which remain accurate.


**Status: CLOSED as a measured negative (see banner).** Originally deferred 2026-07-30 in
favour of experimenting with particle-aware load balancing (T2.1) first,
because the two interact — see "Relationship to T2.1" below.

This document is the pickup point for that work. It records what was
measured, what the code actually does, the proposed change, and the three
ways it can go wrong. Anyone (human or agent) resuming this should be able
to work from this file plus the archived bench results.

## The finding

At the production anchor (1 Milan node, 128 ranks, 5 root steps from
RD0016 cycle 2110), with the full Tier-1 remediation batch already
landed:

| Section | seconds | % of 2214.5 s wall |
|---|---|---|
| `SPFinalize` (StarParticleFinalize) | 1387.7 | **62.7%** |
| ↳ `SPCommUpdateCount` (the Allreduce) | 1387.7 | 62.7% |
| `SPInitialize` (StarParticleInitialize) | 5.6 | 0.3% |
| ↳ `SPFindTotalParticles` | 3.5 | 0.2% |
| ↳ `SPRecordStarCount` | 2.1 | 0.1% |
| ChemistryCooling | 187.7 | 8.5% |
| RebuildHierarchy | 176.1 | 8.0% |
| SolveHydroEquations | 109.9 | 5.0% |
| SetBoundaryConditions | 105.2 | 4.7% |
| everything else | ~230 | ~10% |

`SPFinalize` and `SPCommUpdateCount` agree to 0.1 s. **The entire cost of
StarParticleFinalize is the single `MPI_Allreduce` inside
`CommunicationUpdateStarParticleCount`.** Every other branch of Finalize
early-returns because `AllStars` is NULL (a consequence of the T1.1 guard,
which is correct for FOGGIE: no Star-framework maker is enabled).

Archive: `results/instr4-r1/` on the `bench-results` branch (state.json,
performance.out, time_budget.txt). Preceding passes: `results/instr2-r1`,
`results/instr3-r1`.

## It is barrier wait, not communication

- 1387.7 s over 4327 subcycles = **320 ms per call**.
- The payload is `2*NumberOfGrids` ints (~1–8 kB) among 128 ranks **on a
  single node**, i.e. shared-memory MPI. True cost: tens of microseconds.
- So ~100% of the measured time is ranks waiting for each other, not
  moving data. Four orders of magnitude apart.

Subcycle counts for the 5-root-step run (from `Level[N]: dt =` lines in
`enzo_bench.log`): L6 132, L7 291, L8 559, L9 1124, L10 2140; 4327 total
across all levels. Every one of them pays this barrier.

## The imbalance fluctuates (this is the favourable case)

Per-root-cycle spread of `SPCommUpdateCount` across the 128 ranks
(mean / stddev / min / max, seconds):

    293.1  22.6  188.8  338.3
    256.5  16.6  205.1  296.1
    252.5  20.3  178.1  284.3
    271.2  24.1  130.3  311.0
    314.3  21.9  237.2  354.3

**No rank is ever near zero.** If one rank were persistently the slowest —
doing real work while the other 127 waited — that rank would show ~0 s at
the barrier. Instead the minimum is 130–237 s. So the slowest rank varies
from subcycle to subcycle: the imbalance is fluctuating, not fixed.

Why that matters: wall time is currently `Σ_subcycles max_rank(work)`.
Remove the per-subcycle barrier and it tends toward
`max_rank(Σ_subcycles work)`. Since `max(Σ) ≤ Σ max`, and the gap widens
exactly when the argmax varies, this is the regime where barrier removal
pays. If instead one rank had been permanently slow, the wait would simply
relocate downstream and only real load balancing (T2.1) would help.

Rough scale: summed non-barrier work is ~810 s of the 2214 s. Reducing
4327 collectives to ~55 (one per level per root cycle) could approach a
**2x** wall-clock improvement — against ~7% for the entire Tier-1 batch.

## Why it is unnecessary on almost every subcycle

FOGGIE production sets `StarFormationOncePerRootGridTimeStep = 1`.

- `EvolveLevel.C:415-427`: `SetMakeStars()` is called **only** when
  `level == 0` and that parameter is set, flagging every grid on every
  level, once per root grid timestep.
- `Grid_StarParticleHandler.C` (H2REG branch): runs the maker only if
  `this->MakeStars || !StarFormationOncePerRootGridTimeStep`, then sets
  `this->MakeStars = 0`.

So **new star particles can only appear in the first `StarParticleHandler`
call on each grid per root cycle.** The real job of
`CommunicationUpdateStarParticleCount` is `SetNewParticleIndex` — assigning
global indices to newly created particles. On the other ~99.9% of the 4327
subcycles it Allreduces counts that provably have not changed, and then
assigns no indices.

Relevant code:
- `StarParticleFinalize.C:83` — the call site (now wrapped in
  `TIMER_START("SPCommUpdateCount")`).
- `CommunicationUpdateStarParticleCount.C:94` — the `MPI_Allreduce` over
  `2*GridCount` ints.
- Same file, ~:123-150 — the post-reduce loop: `SetNewParticleIndex` for
  owned grids; for non-owned grids it updates `NumberOfStarParticles`,
  `NumberOfOtherParticles` and calls `SetNumberOfParticles`.

## Proposed change (REVISED 2026-08-02 - supersedes the earlier sketch)

### Why the first design was wrong

The original sketch keyed the skip on `StarFormationOncePerRootGridTimeStep`,
reasoning that new stars can only appear in the first subcycle of each level
per root cycle. **That is FOGGIE-specific, not a fix.** Not all FOGGIE
production runs set that parameter, and with it unset star creation can occur
on any subcycle, so no predicate built from replicated state can decide the
question. Predicting from parameters is the wrong axis.

### The right axis: sync before consumers, not on a cadence

The Allreduce's real job is `SetNewParticleIndex` (Grid_SetNewParticleIndex.C):
walk each grid and assign `ParticleNumber` to particles still flagged
`INT_UNDEFINED`, using running global counters. So it is needed only

  (a) when new particles exist anywhere, and
  (b) before something consumes the resulting indices or counts.

Condition (b) is the useful one, because it is parameter-independent.
Investigation of the consumers:

| consumer | needs the sync? | evidence |
|---|---|---|
| particle migration between grids/ranks | **No** | `Grid_CommunicationSendParticles.C` copies `ParticleNumber` into an opaque `buffer[].id` and restores it; `INT_UNDEFINED` travels intact and is never interpreted |
| `RebuildHierarchy` | **No** | it performs its own `CommunicationSyncNumberOfParticles` |
| data output | **YES** | `Grid_Group_WriteGrid.C` writes `ParticleNumber` to disk; unassigned IDs would corrupt the dump |

Ordering within an EvolveLevel subcycle (line numbers as of `86616544`):

    398  SetLevelTimeStep      <- global dt reduction, ranks in lockstep here
    447  StarParticleInitialize
    716  StarParticleHandler   <- particles created here
    796  StarParticleFinalize  <- the 62% Allreduce
    851  recurse to level+1
    867  OutputFromEvolveLevel <- the real consumer
    1002 RebuildHierarchy      <- does its own count sync

### Design

Move the count sync from *every* `StarParticleFinalize` call to *before data
output*, gated by a parameter defaulting to current behaviour.

- Output is rare in production (at most once per root cycle), versus 4327
  subcycles per 5 root steps, so this is a ~1000x reduction in occurrences
  rather than the ~2x that syncing at rebuild cadence would give.
- The call site moves to an event that is **already global and already
  synchronised**, which is why it should be cheap: our own data shows position
  dominates cost far more than the collective itself. `SetLevelTimeStep`'s
  `CommunicationMinValue` runs at the *same per-subcycle cadence* and costs
  **50 s against 1388 s** - a 28x difference - purely because it sits at the
  start of a subcycle where ranks are still in lockstep, while the star-count
  Allreduce sits at the end and absorbs all the imbalanced work.
- **No rank-divergence risk.** There is no predicate to disagree about: every
  rank syncs at the same global event. This removes the deadlock hazard that
  killed the parameter-based design and that T2.1 demonstrated for real.

### The bookkeeping hazard this creates

`SetNewParticleIndex` assigns IDs from running counters
(`NumberOfStarParticles`, `NumberOfOtherParticles`). `StarParticleInitialize`
recomputes `NumberOfOtherParticles` from `NumberOfStarParticles` **every
subcycle**, so deferring the sync leaves those globals stale in between.

That is tolerable **only if nothing consumes the stale values**. The one
consumer identified is `SetNewParticleIndex` itself, which runs inside the
deferred sync - so counters and assignment move together and stay mutually
consistent. **This is the single most important thing to verify before
trusting the change**: if any other code path reads `NumberOfStarParticles` or
`NumberOfOtherParticles` between syncs and acts on it, IDs can collide, and a
duplicate particle ID is a silent data-corruption bug rather than a crash.

Verification plan for that specific hazard: assert ID uniqueness across the
whole hierarchy at each output (cheap, and catches collisions directly rather
than inferring their absence).

## Risks (in order of severity)

1. **DEADLOCK — the skip condition must not diverge across ranks.**
   A collective is all-or-nothing: if rank A enters the Allreduce and rank
   B skips it, the run hangs. The obvious local test ("did *my* grids
   create particles?") is exactly this bug.
   **Specific trap:** `MakeStars` itself is *not* consistent across ranks.
   `Grid_StarParticleHandler` returns early on the
   `MyProcessorNumber != ProcessorNumber` check **before** reaching
   `this->MakeStars = 0`, so non-owning ranks never clear the flag. Do not
   use `MakeStars` in the condition.
   If a genuinely data-dependent condition is ever needed, it must be
   agreed via its own collective — which defeats the purpose, since the
   cost here is the synchronization, not the payload. Prefer a
   conservative replicated-state predicate.

2. **STALE COUNTS.** The Allreduce also refreshes *non-owned* grids'
   particle counts as a side effect (`SetNumberOfParticles`, and the
   running `NumberOfStarParticles` / `NumberOfOtherParticles` globals).
   Particles migrate between grids at rebuild, which with
   `RebuildHierarchyCycleSkip = 2` happens every 2 subcycles — more often
   than once per root cycle. Before skipping, confirm nothing between
   rebuilds consumes those counts (candidates: output paths,
   `FindTotalNumberOfParticles` consumers, load-balance work estimates,
   anything reading `Grid::NumberOfParticles` for a non-owned grid).
   If something does, the fix may need to keep a cheap count sync at
   rebuild cadence while dropping the per-subcycle one.

3. **WAIT MIGRATION.** Removing a barrier does not remove imbalance. Some
   of the recovered time will reappear at the next synchronization point
   (`MPIWaitReceive` in the boundary exchange, the dt Allreduce in
   `SetLevelTimeStep`, RebuildHierarchy's collectives). The fluctuating-
   imbalance evidence above says a large fraction should genuinely
   disappear, but only a bench proves it. **Expect the measured win to be
   smaller than 62.7%.**

## Validation plan

Same three-level chain used for every item this session:

1. Local serial A/B on `run/StarParticle/TabularFeedbackTest` — must be
   bitwise identical (h5diff on every dump, byte-identical feedback logs).
   Serial runs cannot exercise the deadlock, so this proves only that the
   physics is untouched.
2. Branch CI (build, smoke decks incl. np=4, gravity gate, both serial
   determinism pairs, ASan/UBSan). The np=4 decks are the first real
   multi-rank exercise of the skip condition.
3. Production bench at the anchor gated by the envelope noise floor
   (`results/t18-instrumentation-r7/noise_floor_envelope_A123.json`),
   using `compare_runs.py --class C0 --noise-floor ...`, plus
   `time_budget.py` on the result. **Watch specifically whether the
   recovered `SPCommUpdateCount` time reappears in `MPIWaitReceive` /
   `SetLevelTimeStep` — that is risk 3 made visible.**

A run that hangs at a collective is the deadlock signature; it will
consume its full walltime rather than failing fast.

## Relationship to T2.1 (why this was deferred)

The two changes attack the same underlying problem from opposite ends:

- **This change** removes the *sites where imbalance is paid*.
- **T2.1** (particle-aware load-balance work estimate) reduces *the
  imbalance itself*. The balancer currently weights grids by cell count
  only (`ComputeTime[i] = float(NumberOfCells)`), ignoring particles
  entirely — and FOGGIE concentrates ~1.4M star particles in a small
  region, so particle work is precisely what the estimate is blind to.

Doing T2.1 first is the better experiment: if the imbalance shrinks, the
barrier cost shrinks with it, and re-running this bench afterwards
measures how much of the 62.7% was ever attributable to imbalance versus
to something else. It also avoids attributing a T2.1 win to this change or
vice versa. Note T0.2 (Hilbert balancing) measured as wall-neutral
*because* it balanced cells-only weights — the same blindness.

## Reproducing the measurement

    python3 run/foggie_bench/time_budget.py <rundir>       # section budget
    grep "^SPCommUpdateCount " <rundir>/performance.out     # per-rank spread
    grep -c '^Level\[[0-9]*\]: dt = ' <rundir>/enzo_bench.log  # subcycle count

The timers themselves landed in instrumentation passes 2–4 (shas
`dc98ab91`, `ed40e892`, `f62e1945`); they are permanent, cost 0.05% or
less, and are validated bitwise.

## Addendum 2026-07-30: the deadlock risk is not theoretical

While implementing T2.1 (particle-aware load balancing) the exact hazard
described under risk 1 was hit for real, which makes it worth recording
concretely here.

`CommunicationLoadBalanceGrids` runs from `RebuildHierarchy.C:606`,
**before** the `CommunicationCollectParticles` sync at line 628. At that
point a newly created subgrid's particle count is known only to its
owning rank; non-owners hold zero or a stale value. Feeding that count
into the work estimate made each rank compute a different `ComputeTime`,
so ranks disagreed about which grids to move and the transfers
mismatched — "MPI communication error" a few seconds into the run, at
every non-zero weight, while the weight-0 control ran clean.

Two transferable lessons for the SPFinalize work:

1. **Cell counts are replicated metadata; particle counts are not.** Any
   collective decision built on particle state needs an explicit
   agreement step first. In T2.1 that is one `MPI_Allreduce` (owner
   reports, others report zero, MPI_SUM) before the estimate is formed.
2. **Serial tests cannot catch this class of bug**, and neither can a
   single-rank bitwise check. It first appears at np>1. Validate any
   skip-the-collective change at np=4 minimum before spending a
   128-rank bench slot.

Also note the affordability argument that made the extra reduction
acceptable in T2.1, since it applies in reverse here: `RebuildHierarchy`
measures max/mean = 1.01, i.e. ranks are already in lockstep inside it,
so an added collective there costs little. The per-subcycle collective
in `StarParticleFinalize` is expensive for precisely the opposite
reason — it sits where ranks arrive at wildly different times.
