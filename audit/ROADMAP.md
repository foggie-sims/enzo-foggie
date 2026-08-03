# FOGGIE/Enzo performance roadmap

**Status as of 2026-08-03.** This is the *measured* plan. The original
audit (`PERFORMANCE_AUDIT.md`) ranked its recommendations by inspection;
several of those rankings have since been overturned by benchmarking at
the production anchor. Where the two disagree, this file wins — and says
why. Per-item status lives in `status.yml` (rendered to `dashboard.html`).

Anchor for every number below: halo 008508 RD0016, 5 root steps,
**1 Aitken Milan node × 128 ranks**, node-to-node scatter ~1.6%.

---

## 1. Where the time actually goes

Instrumentation passes 2–4 closed the wall-clock accounting from 37% to
**99.8%** (`run/foggie_bench/time_budget.py` on any bench run).

At the control configuration (~2270 s):

| section | s | % |
|---|---|---|
| `SPCommUpdateCount` — one `MPI_Allreduce` | ~1410 | **62%** |
| ChemistryCooling (Grackle) | 190 | 8.4% |
| RebuildHierarchy | 180 | 8.0% |
| SolveHydroEquations | 110 | 4.9% |
| SetBoundaryConditions | 107 | 4.7% |
| everything else | ~270 | ~12% |

**62% of wall is one collective**, in `CommunicationUpdateStarParticleCount`,
called once per subcycle of every level (4327 times per 5-step run). Its
payload is 1–8 kB among 128 ranks on a single node — tens of microseconds
of actual communication against 320 ms measured per call. It is **not
communication; it is ranks waiting for each other.**

## 2. Why they wait: granularity, not work estimation

The imbalance is structural and indivisible:

- L10 carries 195 grids over 128 ranks = **1.52 grids/rank**.
- Grid sizes span 24 → 15,552 cells (648× spread).
- A fair share is 846,880/128 = 6,616 cells, so **the largest single grid
  is 2.35× a fair share** — and a grid cannot be split across ranks.
- Measured per-subcycle max/mean work: **2.71×**. The two agree.

No work-estimate weighting can fix an oversized atom. This is the single
most important finding of the effort, and it explains three separate
negative results below.

## 3. Completed, with measured effect

| item | outcome |
|---|---|
| **T0.3** subgrid sizing | **−25.8% wall** — the largest win. See below. |
| **T1.5** cooling/prologue gating | −7% wall (removed a per-subcycle Grackle solve) |
| T1.1, T1.2, T1.3, T1.6, T1.13 | correctness/waste fixes, individually at or below scatter |
| T1.8 | instrumentation (0.05% overhead) — enabled everything else |
| T1.9 | O(S²) sibling removal (validated earlier) |
| T0.6 | removed a 5-int `MPI_Bcast` from every subcycle of every level |
| T0.8, CI.1–CI.6 | build guards and CI harness |

### T0.3 — the audit's recommendation was backwards
*Upstreamed as **PR #73** (branch `subgrid-autoadjust-floor`, commit
`5f507d80`), an isolated branch off `main` carrying only this change. PR
text: `foggie_bench_root/PR_subgrid_autoadjust_floor.md`.*

The audit said `MinimumSubgridEdge=16, MaximumSubgridSize~2.6e5`, i.e.
*bigger* grids, to cut per-grid overhead and O(S²) sibling searches. That
rationale predates T1.9, which removed the sibling cost. Measured:

- audit's version as written: **+70% wall** (3853 s)
- opposite direction (`MaximumSubgridSize=400`): −24%
- **proper fix: `SubgridSizeAutoAdjustMinimum=512`: −25.8%** (1685.7 s
  two-sample mean vs 2272.5 s), barrier −43%

Root cause was a single constant. `DetermineSubgridSizeExtrema` already
computes the correct cap — `NumberOfCells/(NumberOfProcessors ×
OptimalSubgridsPerProcessor)`, which is 364–413 cells at L9/L10 here —
then discards it via `max(., MINIMUM_SIZE)` with a hard-coded
`#define MINIMUM_SIZE 2000` dating from lower core counts. That floor is
now the parameter `SubgridSizeAutoAdjustMinimum` (default 2000 =
bit-identical to prior behaviour).

**512 is an optimum, not a direction.** Floor 256 balances *better*
(barrier 762 s, chem skew 1.87) yet runs slower: per-grid overhead
overtakes the balance gain (SetBC 133→153 s, PrepareDensityField
119→143 s, CommReceiveHandler 152→182 s, for only 27 s of barrier).

C2 correctness: baryon 4.5e-10, dark matter 4.0e-16 (machine precision),
star count −0.075%.

## 4. Measured negatives — do not repeat these

| item | result | why |
|---|---|---|
| **T2.1** particle-weighted balancing | **+65% to +196% wall** | Total chemistry work is *constant* across weights (192/189/192/193 s) while its skew explodes 2.85→14.88. Weighting particles doesn't change how much work exists, it relocates cells onto the wrong ranks. Directly measured particle work is only 0.6% of wall. |
| **T0.1** deeper rebuild skip | answer changes far beyond noise for −0.6% wall | Production *already* runs `RebuildHierarchyCycleSkip=2` everywhere. At skip=3 particle drift (2.4 cells) exceeds the 2-cell flagging buffer, as the drift budget predicts. |
| **T0.2** `LoadBalancing=4` (Hilbert) | wall-neutral | Locality gain (~2% off non-rebuild phases) offset by balancer cost. It was balancing *cells-only* weights — the wrong quantity, and bounded by granularity anyway. Path is now correct and crash-free; revisit only after §5.1. |
| **T1.11** Hilbert bbox normalisation | wontfix | Superseded by the 126-bit integer key rework: anisotropic distortion plus key churn as the bounding box moves. |

## 5. Prioritised remaining work

### 5.1 T2.10 — signature reuse *(RETRACTED as a granularity lever, 2026-08-02)*

**An earlier version of this roadmap claimed T2.10 was the highest-value
remaining item, on the grounds that a cap of 400 "should yield ~2117
grids at L10" but only produced 436, supposedly because Berger-Rigoutsos
cannot split a solidly flagged region. Both halves of that were wrong.**

- A geometric bisection fallback **already exists**:
  `IdentifyNewSubgridsBySignature.C` calls
  `LargeAxisRatioCheck(StrongestDim, GridEnds, 0.0)` when no inflection
  point yields sufficiently thick grids, which splits along the long
  axis unconditionally.
- The 2117 figure was a unit error. `MaximumSubgridSize` is compared
  against the ProtoSubgrid size in **parent** cells, while the child
  grid created is 8x larger (RefineBy=2 in 3D) — see the telltale
  commented-out `//size *= 8;` in `ProtoSubgrid_AcceptableGrid.C`.

Measured directly from the 20-step dumps, the cap is enforced exactly:

| config | L10 grids | largest grid | 8 x cap |
|---|---|---|---|
| floor 2000 | 208 | 15,120 cells | 16,000 |
| floor 512 | 378 | 4,032 cells | 4,096 |

Largest grids sit at 95-98% of the limit. The splitter is working
correctly and there is **no hidden granularity headroom** behind it.

Consequence: the granularity avenue is essentially exhausted. It is
controlled entirely by the cap, and the sweep already located the
optimum (512; 256 balances better but loses to per-grid overhead).
T2.10's actual content — signature reuse and the multi-dimensional
inflection search — remains a legitimate but ordinary optimisation of
the regridding cost, not a route to better load balance. RebuildHierarchy
is ~8% of wall, so it is bounded accordingly.

### 5.2 `SPFinalize` — CLOSED, measured ceiling 4.6% *(2026-08-02)*

Measured before building: skipping the per-subcycle star-count Allreduce
entirely removed **832 s of barrier and saved 78 s of wall (4.6%)**. About
90% relocated to `MPIWaitReceive` (+389 s), `SetBoundaryConditions`
(+360 s) and the dt reduction (+164 s); real work was unchanged.

The 62% figure was never the cost of a collective — it was where real
work imbalance happened to be paid. Remove that barrier and the next one
absorbs it. See `SPFinalize_Edits.md` for the full record and for why the
`max(Σ) ≤ Σ max` argument did not hold (it assumed a single
synchronisation point; ghost-zone exchange and the dt reduction are
others).

**The transferable lesson for everything remaining:** only changes that
reduce the *imbalance itself* — or the work — pay off. Changes that move
where imbalance is paid do not. T0.3 succeeded because it eliminated
indivisible oversized grids; particle-weighted balancing failed because
it worsened the imbalance.

### 5.3 Measured-work load balancing (T2.1's remaining half) *(measured 2026-08-03)*

**Why chemistry is the only target left.** After T0.3 every section whose
cost is proportional to cell count is already well balanced — rebuild
1.01, boundaries 1.18, density 1.20, hydro 1.23 max/mean. `ChemistryCooling`
alone sits at **2.13**, because Grackle's cost per cell tracks gas state,
which is precisely what a cell-count weight cannot see.

**Design finding: the obvious implementation does not work.** The
`InputGridWork` overload is not sufficient by itself, because grid objects
do **not** survive `RebuildHierarchy` — it deletes the old subgrids
(`:370`), creates new ones from the flagging (`:449`), and only then
balances them (`:606`). At the moment the balancer needs an estimate, the
grids it is balancing did not exist during the cycle that was measured.
Per-grid history has no referent; work is a property of a *region of
space*, not of an object, so an estimate can only be carried across a
rebuild **spatially**.

The implementation that follows from this is a **per-cell work field**,
handled like density or temperature: filled during the subcycle, then
carried onto the new grids by machinery that already exists —
`InterpolateFieldValues` (`:546`) from the parent, `CopyZonesFromGrid`
(`:562`) from the old same-level grids, both of which run *before* the
balancer at `:606`, with the old grids still alive until `:634`. Freshly
refined regions inherit the parent's measured density rather than a
level average. Two constraints: it must be skipped by the hydro solver
(it is a diagnostic density, not a fluid quantity, and advecting it would
smear the signal), and it is filled uniformly per grid at `work/cells`,
since Grackle is called per grid — which is exactly the resolution the
balancer needs, because grids are its atoms.

**Measured ceiling** (`GridWorkMapOutput=1`, offline replay of a greedy
balancer scored against *actual* work, `harness/analyze_workmap.py`):

| balancer input | chemistry imbalance |
|---|---|
| cells (today) | 3.89× |
| measured work | 3.36× |
| perfect foreknowledge | 3.35× |

The slow-evolution premise holds decisively: interval-to-interval
correlation is **0.99** at every level from L5 down, which is 98% of all
chemistry work. The design is then essentially optimal — it recovers
**98%** of the room between the cells-based estimate and perfect
foreknowledge.

**But the room itself is only 13.9%.** The chemistry critical path goes
680.3 s → 587.5 s over 5 root steps, against a bench wall of ~1670 s.
Because this model sums per-level maxima (different ranks can be the
slowest at different levels), it overstates what any single rank
experiences, so **~5% of wall is an optimistic bound and ~3–5% is the
honest range** — not the 15–25% first estimated from the raw skew.
Measured-work balancing is worth doing, but it is a single-digit item:
the estimate was never the main problem.

**What the oracle floor points at instead.** The limit is granularity
again, and it is worst exactly where the work is. Weighted by share of
total chemistry, the oracle floor per level is: L10 **2.97×** (45.9% of
chemistry), L9 **5.26×** (21.2%), L8 **4.85×** (12.0%), against L6
**1.01×** (7.8%, effectively perfect). Two thirds of all chemistry sits
at levels whose *best achievable* imbalance is 3–5×, set by the largest
single grid's share of that level's work. The
natural successor to T0.3 is therefore to make
`DetermineSubgridSizeExtrema` cap **predicted work** rather than cell
count — splitting dense regions finer and leaving diffuse ones coarse.
The same work field would drive both the sizing and the balancing. T0.3
shrank the atom uniformly for 25–41%; this is the targeted version of the
same move, and on these numbers it is the larger lever of the two.

*Measurement history, since two of these were wrong before they were
right.* The first run recorded only at root-step boundaries and captured
**5.9%** of chemistry — grids die at every rebuild and take their
accumulated time with them, biasing hardest against the deepest levels,
which rebuild most. Recording moved to the top of `RebuildHierarchy`,
before destruction: capture went to ~100% (23038 s against a section-timer
scale of 22812 s), across 876 rebuild intervals. The offline replay then
had two bugs of its own — `argmin` on an all-zero load array piles every
zero-weight grid onto rank 0, and the two dump sites interleave, so some
"intervals" were a rebuild dump followed immediately by a root-step dump
with near-zero elapsed work. Fixed by tie-breaking on grid count and
dropping degenerate pairs (50 of 1628). The tell was a self-contradiction:
82× imbalance reported alongside 0.99 correlation. The per-cell field
design makes the *sampling* failure structurally impossible, since work
then accumulates in space rather than on an object.

**Not a constraint:** subgrids are *not* pinned to their parent's rank.
`LoadBalancing=1` balances over the full rank range (`StartProc=0`,
`EndProc=NumberOfProcessors`), and measurement confirms it — L10's 350
grids occupy all 128 ranks, with only 9 of 350 on a rank that also owns
the root tile they sit in and 171 on ranks owning no root tile at all.
(`LoadBalancing=2` *would* restrict to a `CoresPerNode` window.) The cost
of that freedom is that nearly every parent-child exchange crosses ranks,
which is what T0.2's Hilbert ordering tried to recover — and plausibly why
it came out wall-neutral, since a space-filling curve hands the
concentrated expensive region to a single rank.

The rank-agreement `Allreduce` added for T2.1 is a prerequisite for any
data-dependent weighting and is already in place.

### 5.3b T0.10 — Grackle is built without vectorisation *(new, 2026-08-02)*

ChemistryCooling is the **largest single compute item** in the measured
budget — 8.5% of wall — and it is entirely outside Enzo's build. Grackle
is linked as a prebuilt external library, so no Enzo compiler flag
reaches it. Its build config:

    CMAKE_BUILD_TYPE            = Release
    CMAKE_C_FLAGS_RELEASE       = -O3 -DNDEBUG
    CMAKE_Fortran_FLAGS_RELEASE = -O3 -DNDEBUG -O3
    -march / -x / -xHost        = ZERO occurrences

So its Fortran kernels target the generic x86-64 baseline (SSE2), not the
AVX2 the Milan nodes support. This is the **same defect the audit flagged
for Pleiades** (5.1.6, "a plausible 1.5-3x slowdown") — except that one is
now moot (Pleiades decommissioned) while this one is live, on the biggest
compute item we have.

It is also built with `/usr/bin/gfortran` while Enzo uses Intel
2020.4.304, which is why the executable carries a `libgfortran.so.3`
dependency.

Experiment: rebuild Grackle with `-O3 -march=core-avx2`, relink, bench.
Ceiling is ~8.5% of wall if the kernels fully vectorise; realistically
less. Touches no Enzo source. **C1, not C0** — it changes chemistry
results at roundoff, so it needs the noise-floor comparison rather than a
bitwise gate.

### 5.4 OPEN ISSUE — T0.3's systematic star-particle deficit

> **Update 2026-08-03 — the padding mechanism is refuted, the deficit is
> still unexplained.** The buffer-zone test proposed below has been run
> (`t03-buffer-r1`, 20 steps, floor 512 with `NumberOfBufferZones=3`;
> note FOGGIE already runs 2, not the Enzo default of 1, so this was
> 2→3). It fails on both counts. Restoring refined padding did **not**
> close the deficit — it widened slightly, −0.373% vs −0.334% at step 20
> — and it erased the speedup entirely: on the same metric (exclusive
> `Level_NN` time summed over 20 root steps) the control is 8179.1 s,
> floor 512 with buffer 2 is 5821.4 s (−28.8%), and buffer 3 is
> **8665.5 s, +5.9% — slower than the control**, because the extra
> refined volume costs more than the granularity gain returns.
>
> So the hypothesis stated below — that tighter grids refine ~15% less
> unflagged padding and thereby change how many marginal cells cross the
> star-formation threshold — is **measured and wrong**. Restoring the
> padding does not restore the stars. PR #73's disclosure stands as
> written *except* for that mechanism claim, which should be corrected.
>
> Untested candidates: the SF density threshold interacting with
> cell-size-dependent density estimates, or changed subgrid boundaries
> altering where feedback is deposited.

The 20-step de-risk run (`results/t03-long-r1`) confirmed the speedup
holds (24.8% at 20 steps vs 25.8% at 5) but surfaced an unresolved
problem: **`SubgridSizeAutoAdjustMinimum = 512` forms systematically
fewer star particles.**

| root step | stars (floor 2000) | stars (floor 512) | deficit |
|---|---|---|---|
| 1 | 1,418,504 | 1,418,304 | 0.014% |
| 5 | 1,423,057 | 1,421,994 | 0.075% |
| 10 | 1,429,207 | 1,426,902 | 0.161% |
| 15 | 1,435,493 | 1,431,866 | 0.253% |
| 20 | 1,441,921 | 1,437,105 | 0.334% |

Why this is a real effect and not noise:

- **The sign never flips.** Chaotic divergence wanders both ways; a
  consistent deficit at every checkpoint means a systematic cause.
- **The rate is constant** — 0.015/0.017/0.018/0.016 % per root step —
  i.e. linear, with no sign of saturating over this window.
- It is **~1000x the identical-code noise floor**, which is +-1 star.
- Conservation is impeccable throughout (dark matter at machine
  precision, baryon 1e-9), so nothing is leaking. This is a change in
  how many stars *form*, not a bookkeeping error.

Likely mechanism: smaller subgrids fit the flagged region more tightly
and refine 15% less unflagged padding, which changes how many marginal
cells cross the star-formation density threshold.

**Status: not a blocker for short runs or testing; unquantified for long
production campaigns.** Naive linear extrapolation gives ~1.7% at 100
root steps and ~8% at 500, but 20 steps is far too short to know whether
it saturates, and that arithmetic should not be trusted. Do not commit a
long campaign to this setting on the strength of the speedup alone.

Testable next step: a floor=512 arm with `NumberOfBufferZones` raised by
one, to see whether restoring the refined volume closes the deficit and
how much of the 24.8% survives.

### 5.5 Smaller / unsequenced
- **T0.7** compiler flags — SCOPE REDUCED. Its headline benefit was
  Pleiades-only and is moot (decommissioned); Aitken is already
  vectorised and Grackle is unreachable (see 5.3b), so `-ip` /
  `-heap-arrays` reach only ~5% of wall. Variants built and queued;
  expect low single digits.
- **T0.4, T0.5, T0.9**, T1.7, T1.10, T1.12 — untouched.
- Re-stamp **T1.9** at 1×128 (its production stamp predates the config fix).
- Longer-baseline bench (10+ root steps) before committing a production
  campaign to the T0.3 setting.

## 6. Methodology notes that cost real time to learn

- **Noise floors are config-specific and must be multi-pair.** Single-pair
  floors scatter ~5× at step 5 and produced a false FAIL. Use
  `merge_noise_floors.py` over ≥3 baseline runs.
- **The production config is 1 node × 128 ranks `mil_ait`.** Earlier benches
  ran 512 ranks on 4× `rom_ait` — an accidental `make_bench_run.py`
  default, never a production configuration.
- **`compare_runs.py` mass sums were wrong until 2026-07-31** (commit
  `92e436d0`). The HDF5 grid groups carry no attributes, so cell volume
  silently defaulted to 1.0: no volume weighting, no particle `dx³`, and
  refined regions double-counted. It reported impossible 2.4% dark-matter
  differences. Same-structure (C0/C1) verdicts are unaffected — the errors
  cancel — but every **C2** verdict needs the corrected metric.
- **Header edits do not trigger recompiles here** (`makedepend` absent,
  `DEPEND` is 0 bytes). After touching any header, `make clean` first. The
  benign symptom is a link error; the dangerous one is a changed struct
  layout, which links fine and is UB at runtime.
- **Validate collective-affecting changes at np≥4.** Serial tests and
  single-rank bitwise checks cannot see rank divergence.
- **Bench runs are ~50 GB each** and filled the quota mid-session. Archive
  and prune dumps; conclusions live in `bench-results`.
