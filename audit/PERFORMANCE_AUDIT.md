# Enzo Performance & Scaling Audit for FOGGIE

**Scope:** `src/enzo` at commit `4ceb000` (foggie-sims/enzo-foggie).
**Focus:** the FOGGIE regime — deep AMR zooms (10+ levels, nested must-refine regions,
`MultiRefineRegion` tracking) with 10⁶–10⁷ star particles concentrated in a few galaxies
on the finest levels — where Enzo is known to scale poorly.
**Method:** six parallel code audits (hierarchy construction, time integration, MPI
communication, load balancing, gravity/particles, star formation/feedback), each reading
the relevant sources line-by-line, plus an architectural comparison against the public
RAMSES and Nyx repositories and the ART code as described in its method papers.
All file references are to `src/enzo/` unless otherwise noted; line numbers are as of
`4ceb000`.

---

## 1. Executive summary

Enzo's poor scaling in the FOGGIE regime is not one bottleneck but the product of three
architectural multipliers, each of which the FOGGIE configuration maximizes:

1. **The subcycle multiplier.** Level L executes ~2^L times per root-grid step
   (`EvolveLevel.C:391`), and every subcycle is a globally synchronous phase: a timestep
   allreduce, ~14 communication rounds, star-particle bookkeeping collectives, and (by
   default) a full rebuild of the entire sub-hierarchy. At L_max = 10–11 this is
   ~2000 level-executions and **~4×10⁴ global collectives per root step** before any
   feedback-related communication.

2. **The replicated-hierarchy multiplier.** Every MPI rank holds a full copy of the grid
   hierarchy — not just metadata, but an instantiated `grid` object (~4 kB) per remote
   grid (`CommunicationShareGrids.C:176-241`) — and every hot loop iterates over the
   *global* grid list, filtering ownership inside the call. Because the auto-tuner
   targets 16 grids per rank (`DetermineSubgridSizeExtrema.C:45-53`), grid count grows
   linearly with rank count, so per-rank metadata memory and per-rank loop overhead grow
   as O(N_ranks), and aggregate bookkeeping as O(N_ranks²). Adding processors makes
   this component of the code *slower*.

3. **The concentrated-particle multiplier.** The load balancer's work estimate is cell
   count only — the particle count is computed and then explicitly discarded
   (`CommunicationLoadBalanceGrids.C:84`). Every per-particle cost (deposits, feedback
   table lookups, particle Courant scans, rebuild churn) lands on the handful of ranks
   that own the galaxies, and the bulk-synchronous structure makes everyone else wait.
   Meanwhile particle deposits for gravity recur once per ancestor level per subcycle,
   i.e. **~2^(L+1) CIC deposits per particle per root step** (`DepositParticleMassField.C`),
   and the star feedback path performs 6 yield-table integrations per star per subcycle
   with no age cutoff before the calls (`star_feedback6.F`).

On top of these, the audit found a large number of specific defects: O(N²) loops on hot
paths, a serialized O(N_ranks)-step transpose inside the gravity FFT, per-star global
allreduce storms in feedback spheres, a chaining mesh that degenerates to a single
bucket for zoom geometries, Hilbert load-balancing keys that run out of resolution at
exactly level ~11, dense allreduces of sparse data, several latent correctness bugs, and
a June 2025 build-flag change that disabled vectorization on Pleiades.

The good news: a substantial fraction of the total cost can be recovered with
**parameter-file and build-flag changes alone** (§7.1), and another large fraction with
**localized, low-risk code fixes** (§7.2). The architectural gaps relative to
RAMSES/Nyx/ART (§6) — distributed metadata, composite gravity solves, threading — are
long-term projects, but the comparison also shows Enzo's patch-based design is not
inherently the problem; the implementation predates the scales FOGGIE runs at.

---

## 2. Why the FOGGIE regime is Enzo's worst case

A useful mental model for the cost of one root-grid step:

```
Cost ≈ Σ_L  2^L  ×  [ C_sync(N_ranks)  +  C_meta(N_grids_L)  +  C_local(work_L) ]
```

- `C_sync` — per-subcycle collectives and communication rounds. Latency-bound; grows
  with log(N_ranks) per collective but the *count* of collectives is fixed per subcycle,
  so it is multiplied by ~2^(L_max+1) ≈ 2000–4000 per root step.
- `C_meta` — per-rank traversal of the replicated grid list (~10 loops per subcycle over
  all grids on the level, plus the rebuild machinery). Independent of rank count;
  proportional to *total* grid count, which itself is proportional to rank count under
  the auto-tuner. This is the anti-scaling term.
- `C_local` — actual physics on locally owned grids. This is the only term that shrinks
  as ranks are added, and only if the load balancer distributes it — which, for
  particle-dominated work, it currently does not attempt.

FOGGIE's configuration pushes all three terms to their worst case simultaneously: deep
hierarchy (maximizes the 2^L multiplier), high rank counts (maximizes N_grids via the
auto-tuner and collective latency), and extreme spatial concentration of particles and
fine grids (maximizes load imbalance within C_local and defeats the spatial-locality
assumptions in the chaining mesh and Hilbert balancer).

---

## 3. Findings: architecture and AMR machinery

### 3.1 Replicated hierarchy metadata (the fundamental non-scalable structure)

Every rebuild, `CommunicationShareGrids.C` (called from `RebuildHierarchy.C:454` once
per level) does an `MPI_Allgather` + `MPI_Allgatherv` of a ~80-byte `PackedGrid` for
every new subgrid, and then **every rank instantiates a full `HierarchyEntry` + `grid`
object for every remote grid** (`CommunicationShareGrids.C:176-241`), including heap
allocation of `CellLeftEdge`/`CellWidth` arrays per dimension
(`Grid_PrepareGridDerivedQuantities.C:39-40`). At ~4 kB per grid per rank, a
160,000-grid hierarchy costs **~640 MB per rank of pure remote-grid metadata**,
replicated identically everywhere. This is the deep-zoom memory wall, and it is also
what makes every `for (grid1 = 0; grid1 < NumberOfGrids; ...)` loop O(total) rather
than O(local).

Downstream, ownership is filtered *inside* the per-grid calls
(e.g. `Grid_SolveHydroEquations.C:64-65` returns early for non-local grids), so a rank
owning zero grids on level 10 still executes ~10 loops of length N_10, with a virtual
call per grid per loop, at every one of the 1024 level-10 subcycles.

**Recommendations:** (i) short term, build a per-rank `LocalGrids[]` index array once
per `EvolveLevel` entry and drive all local-only loops from it — mechanical and low
risk; (ii) medium term, replace the remote-grid `grid` object with a lightweight proxy
(rank, dims, edges, processor, parent) and compute `CellLeftEdge`/`CellWidth` on demand
(uniform spacing makes them O(1)); (iii) long term, share grid metadata only with the
neighborhood that can overlap (the chaining mesh can determine this), moving from
Allgatherv to a sparse exchange. Item (iii) is the Enzo equivalent of RAMSES's "locally
essential tree" and is the structural fix.

### 3.2 RebuildHierarchy: frequency and internal costs

- **Frequency.** The hierarchy below level L is rebuilt after every subcycle except the
  last (`EvolveLevel.C:973-974`), and `RebuildHierarchyCycleSkip[level]` — the only
  throttle — defaults to 1 (never skip) for every level. At L_max = 10 this is ~2000
  rebuild-level iterations per root step; rebuilds are commonly 30–50% of wall clock in
  deep Enzo zooms. **Setting `RebuildHierarchyCycleSkip[L] = 2–4` for L ≥ 4 is the
  single highest-value zero-code change available.**
- **O(N²) parent grouping.** `RebuildHierarchy.C:203-242` groups grids by parent with a
  naive double loop — O(N²/2) pointer comparisons per level per rebuild, unfiltered by
  ownership, executed identically on every rank (~10⁸ ops/level/rebuild/rank at FOGGIE
  scale). Fix: sort or hash by parent pointer, and restrict to locally-owned parents
  (the downstream `MoveAllParticles` already no-ops for remote grids).
- **Dead work and stack hazards.** Five to seven `MAX_NUMBER_OF_SUBGRIDS` (=100,000)
  pointer arrays live on the stack (~4–5 MB of frame in a routine reached through the
  level recursion); `RebuildHierarchy.C:439-441` zeroes 100,000 entries per level per
  rebuild and the filled values are then overwritten before use (dead code);
  `RebuildHierarchy.C:372-375` reads uninitialized memory when
  `SubgridSizeAutoAdjust = FALSE`.
- **Debug collectives.** `RebuildHierarchy.C:425-426` issues two unconditional blocking
  `MPI_Reduce` calls per level per rebuild whose only consumer is an `if (debug)`
  printf. At high rank counts this is tens of seconds of pure synchronization per root
  step. Guard them.

### 3.3 The chaining mesh (FastSiblingLocator) collapses for zoom geometries

`FastSiblingLocatorInitialize.C:44` sizes the sibling-search chaining mesh as
`min(TopGridDims/4, 128)` cells **over the whole domain**, for every level
(`FastSiblingLocatorEntireDomain = TRUE` by default). For a FOGGIE box with a ~10⁻³
zoom region, all grids on levels ≳6 fall into one or two mesh cells, so sibling search
degenerates to O(N_level²) — the accelerator provides zero acceleration exactly where
it is needed. Compounding this, the duplicate check inside
`Grid_FastSiblingLocatorFindSiblings.C:117-141` is a linear rescan of the accumulated
sibling list per candidate (O(S²), before the processor filter and overlap check), and
the mesh construction performs up to millions of individual 16-byte `new`/`delete`
pairs plus a 16 MB pointer-array clear per call.

**Recommendations:** size the mesh to the bounding box of the grids actually being
inserted with O(1) mean occupancy (the bounding-box machinery already exists in
`LoadBalanceHilbertCurve.C:89-96`); replace the dedup rescan with a visit-stamp on the
grid object (O(1), ~10 lines); arena-allocate the chain links. As an immediate
parameter-level experiment, `FastSiblingLocatorEntireDomain = 0` (mesh over the
RefineRegion) should recover most of the deep-level locality, at the cost of pushing
root grids onto an unaccelerated overflow list — worth measuring.

Note also that every overlap loop carries an `#ifdef FAST_SIB` / `#else` pair where the
fallback is a literal all-pairs O(N²) `CheckForOverlap`
(e.g. `PrepareDensityField.C:328-343`, `SetBoundaryConditions.C:164-202`). **Verify
`CONFIG_FAST_SIB = yes` in the production build** — mis-building without it is a silent
~100× regression — and consider deleting the fallback branches.

### 3.4 Subgrid sizing: deep-level grids are mostly ghost zones

With `SubgridSizeAutoAdjust = TRUE` (default), `DetermineSubgridSizeExtrema.C:45-53`
computes `MaximumSubgridSize = N_cells/(N_ranks × 16)`, clamped below by hard-coded
floors `MINIMUM_SIZE = 2000` (parent cells) and `MINIMUM_EDGE = 4`. On deep FOGGIE
levels the floors bind: a 4-parent-cell-edge subgrid becomes a 14³ child grid after
refinement and ghost zones — **81% ghost cells**; even at the 2000-cell ceiling (~31³
child) it is 46% ghost. Consequences: 2–5× inflated memory and boundary-exchange
volume per unit of real work, an exploding message count (message count, not bandwidth,
limits `SetBoundaryConditions` at high rank counts), and grid counts that scale with
rank count (§3.1). The hydro sweeps compound this: `Grid_SolvePPM_DE.C:91` and
`Grid_xEulerSweep.C:174-175` sweep the full transverse ghost extent, a 2.1× flop
overhead at the auto-tuned grid size vs 1.4× at 32³.

**Recommendations:** for FOGGIE, set `SubgridSizeAutoAdjust = 0`,
`MinimumSubgridEdge = 16` (ghost fraction 81% → ~28%), `MaximumSubgridSize` in the
10⁵–2.6×10⁵ range, and benchmark; in code, promote the two floors to runtime
parameters and make the sizing ghost-aware (enforce edge ≥ 2×NumberOfGhostZones as a
floor rationale). Fewer, fatter grids simultaneously improve memory, message counts,
sibling-list sizes, CIC cache locality, and the behavior of the greedy load balancer.

### 3.5 Grid creation (Berger–Rigoutsos) inefficiencies

`IdentifyNewSubgridsBySignature.C` never reuses signatures across bisection splits
(each tree node pays ≥3 full projections of its volume plus 1–2 full copies;
O(V·K) for the skewed trees FOGGIE's nested geometry produces), and the
inflection-point split search has been reduced to the single longest dimension —
`for (i = 0; i < 1; i++)` at line 126, with the original `MAX_DIMENSION` loop commented
out — falling through to a blind bisection when that one dimension has no acceptable
crossing, which systematically produces worse-fitting grids and therefore more of them.
Additionally any proto-subgrid with an edge ≥ 516 cells is force-split regardless of
efficiency (an artifact of the Fortran work-array limit `MAX_ANY_SINGLE_DIRECTION`).

**Recommendations:** slice the parent's signature arrays into children instead of
recomputing (removes ~2/3 of projection cost); restore the multi-dimension inflection
search; carry `NumberFlagged` into children.

### 3.6 Flagging-field costs

Each `CellFlaggingMethod` makes its own independent full-grid passes (≥15 passes per
grid per rebuild for a typical FOGGIE method set), `Grid_FlagBufferZones.C` copies the
full temp buffer back per buffer-zone iteration, and
`Grid_SetFlaggingFieldStaticRegions.C:51` scans all 1000 `MAX_STATIC_REGIONS` slots for
every grid on every level even when only a handful are defined. The must-refine
particle path is the FOGGIE-critical one and is covered in §3.7.

### 3.7 Must-refine particle flagging ships dense fields (FOGGIE-critical)

For remote grids, any rank holding contributing particles computes the **entire dense
`float` flagging field for that grid (ghost zones included) and sends it**
(`Grid_SetParticleMassFlaggingField.C:166-171`); the owner then sums the incoming
buffers element-wise. This is invoked for every level ≤
`MustRefineParticlesRefineToLevel` — i.e. all the nested levels of a FOGGIE zoom — once
per level per rebuild. The source itself concedes the problem
(`RebuildHierarchy.C:407-411`: "there will be scaling issues if one tries to do
nonlocal refinement on large numbers of particles").

**Recommendations:** send a sparse representation — a bitfield (32× smaller), a list of
flagged cell indices, or simply the contributing particle positions for owner-side
deposit. Also exploit the existing `AllLocal` fast path (taken when
`level > MaximumStaticSubgridLevel`) by defining static refine regions down to the
deepest level where particle locality holds.

---

## 4. Findings: parallelism

### 4.1 Load balancing (the FOGGIE scaling wall)

1. **The work estimate ignores particles entirely.**
   `CommunicationLoadBalanceGrids.C:79-87` fetches the per-grid particle count and
   discards it; `ComputeTime[i] = float(NumberOfCells)`. The Hilbert path
   (`LoadBalanceHilbertCurve.C:126-132, 419-425`) does the same. A disk grid with 10⁵
   star particles and a CGM grid with zero are identical to the balancer. Every
   per-particle cost in this report therefore concentrates on a few ranks, unmitigated.
   The infrastructure for a fix already exists: the radiative-transfer balancer feeds a
   measured per-grid work array through the `float *InputGridWork` overload of
   `LoadBalanceHilbertCurve` (`CommunicationLoadBalancePhotonGrids.C:113-115`).
   **Recommendation:** add `ComputeWorkEstimate() = CellsTotal + w_p·N_particles +
   w_s·N_stars` (new parameters, default 0 → behavior unchanged), calibrate the weights
   against per-level timings; longer term, feed measured per-grid times through
   `InputGridWork`. This is the single largest structural payoff available.

2. **The default balancer (mode 1) gives up on first failure.** The greedy loop
   (`CommunicationLoadBalanceGrids.C:119-217`) aborts balancing for the *whole level*
   the first time the smallest grid on the most-loaded rank exceeds half the current
   imbalance (`:197-198`) — exactly the deep-zoom case where one rank holds one big
   grid. There is no grid-splitting fallback. The loop is also serial, replicated on
   every rank, moves one grid per O(N_ranks + N_grids) iteration, and picks first-fit
   rather than best-fit. **Recommendation:** switch production to `LoadBalancing = 4`
   (Hilbert) now; fix mode 1 to mark stuck ranks and continue.

3. **Hilbert keys run out of resolution at level ~11.** `HilbertCurve3D.C` extracts 19
   bits of *domain-absolute* coordinate; for a 256³ root grid the key resolution equals
   the cell size at exactly L = 11, so deeper grids get tied keys, the sort loses
   locality, and the fuzzy-boundary refinement silently no-ops. The fix is **already
   written but commented out**: bounding-box normalization at
   `LoadBalanceHilbertCurve.C:112-113` (and :405-406, :642-643;
   `BoundingBoxWidthInv` is computed and unused). Re-enabling it recovers ~500× key
   resolution for a 10⁻³ zoom. Also widen the key comparison (the two-word key is
   collapsed into a 53-bit double while carrying 57 bits).

4. **Per-level-only balancing; root grids frozen.** Levels are balanced in isolation;
   the only cross-level-aware balancer (`CommunicationLoadBalanceRootGrids`, which
   weights root tiles by descendant cells) **returns immediately for any nested-IC run**
   (`StaticRefineRegionLevel[0] != INT_UNDEFINED` guard at `:60-64`), so FOGGIE root
   tiles keep their initialization-time cyclic assignment
   (`CommunicationPartitionGrid.C:551`, `gridcounter % NumberOfProcessors`) forever —
   including the tiles under the zoom that parent the entire hierarchy.
   **Recommendation:** parameterize the zoom early-return and weight root placement by
   expected subgrid load; longer term, a multi-constraint SFC partition over all levels.

5. **Grids are born on the parent's rank, filled, then migrated.** New subgrids are
   created on the parent's processor, populated by interpolation, and only afterwards
   load-balanced — so full field payloads (plus particles on non-static levels) cross
   the network every rebuild, and the few ranks parenting the zoom are transient memory
   hot spots. **Recommendation:** assign the destination processor *before*
   `InterpolateFieldValues` (geometry and work estimates are known right after
   `CommunicationShareGrids`); the interpolation itself then performs the transfer.

6. **Modes 2/3 are a silent no-op on fresh zoom runs** (`CoresPerNode` stays 1 because
   `DetermineNumberOfNodes()` is unreachable before the early returns; each rank then
   "balances" against itself with no warning). Add a guard/failure and make
   `CoresPerNode` a readable parameter.

Latent bugs in the balancers: an out-of-bounds `HilbertData[-1]` read in
`LoadBalanceHilbertCurveRootGrids.C` (missing `BlockDivisions[i] == 0` guard); `int`
overflow of `TotalWork` above 2.1×10⁹ ghost-inclusive cells (garbage partitions);
unbounded partition `while` loops; a division by zero when a rank receives zero work.

### 4.2 Communication layer

1. **Per-star collective storms (very high impact for FOGGIE).**
   `Star_FindFeedbackSphere.C:113-170` performs **one 7-element `MPI_Allreduce` per
   radius shell per star** while growing feedback spheres
   (`StarParticleAddFeedback.C:131,164`); with hundreds of forming stars this is
   10³–10⁴ global allreduces per level-subcycle. The active-particle feedback-zone
   machinery similarly runs a full three-phase communication round over all grids *per
   particle* (`ConstructFeedbackZone.C:119-141`). **Recommendation:** batch across
   stars — compute all local shells, then one packed allreduce (or pipelined
   `MPI_Iallreduce`); one three-phase round covering all (particle, grid) pairs.

2. **`CommunicationReceiveHandler` is O(N²).** The drain loop
   (`CommunicationReceiveHandler.C:74-343`) rescans the entire receive-record table
   (two full scans, plus a third on active-particle completions) after every
   `MPI_Waitsome` return; with trickling completions this is O(N²) per round, ~14
   rounds per subcycle. With `GRIDS_PER_LOOP = 100000` (batching effectively disabled
   in `SetBoundaryConditions.C`, `PrepareDensityField.C`, `UpdateFromFinerGrids.C`,
   `CopyZonesFromOldGrids.C`, `DepositParticleMassFlaggingField.C`), N reaches
   N_grids × N_siblings. **Recommendation:** dispatch only the completed indices
   `MPI_Waitsome` already returns, keep an explicit deferred list for
   dependency-blocked records, and lower `GRIDS_PER_LOOP` to ~512–1024 (also bounds
   receive-buffer memory and removes the unchecked overflow of
   `MAX_RECEIVE_BUFFERS = 150000`, which currently has **no bounds check anywhere** and
   fails as silent corruption).

3. **The gravity FFT transpose is a serialized ring.** `CommunicationTranspose.C`
   variants 0/1 loop over all processor "jumps" with a blocking send/wait per step —
   O(N_ranks) serialized latency steps per transpose, ×6 transposes per root-grid
   potential solve (~13 `CommunicationTranspose` calls in
   `CommunicationParallelFFT.C`); variant 2 pipelines only 128 at a time. On top of
   this the FFT itself is a 1-D slab decomposition (`CommunicationParallelFFT.C:56-88`)
   — parallelism capped at N_root slabs (256 ranks for a 256³ root grid), with all
   other ranks idle but still exchanging zero-length messages.
   **Recommendations:** replace the ring with a single `MPI_Alltoallv` (the routing
   tables already contain the counts/displacements); ensure `UnigridTranspose = 2` (the
   default — verify it is not overridden); medium term, adopt a pencil decomposition or
   an external FFT (heFFTe/FFTW-MPI), or gather the root solve onto a
   min(P, N_root) sub-communicator. Also note the bundled `s90` serial FFT is used
   unless a vendor macro is set (`select_fft.F90:86`) — several × slower than
   FFTW/MKL, with two extra full-array transposes per call in `wrapper3d.F90`.

4. **Dense collectives of sparse data.** `CommunicationSyncNumberOfParticles` allreduces
   `45 × N_grids` ints (~1.8–18 MB at FOGGIE grid counts) with ~1/N_ranks occupancy, up
   to four times per `CommunicationCollectParticles`, several times per rebuild level.
   The stride grew from 33 to 45 in this fork (PR #54), which also added **12 full
   linear scans of every grid's particle array per sync**
   (`Grid_ReturnNumberOfParticlesOfThisType` per type) to feed an I/O diagnostic.
   `CommunicationShareActiveParticles` uses `MPI_Allgatherv` where `Alltoallv` is
   correct (every rank receives every moving particle — a factor-N_ranks waste; the
   per-destination counts are computed and thrown away).
   **Recommendations:** one-pass 12-bin histogram (or incrementally maintained
   `ParticleTypeCount[]`, which already exists in `Grid.h`); gate the extra 12 slots on
   output steps; sparse owner-keyed allgatherv for the sync; convert
   ShareActiveParticles to the Alltoallv pattern already used correctly by
   `CommunicationShareParticles`.

5. **Collective fusion and stragglers.** Numerous adjacent small reductions can be
   packed: the two timestep allreduces in `SetLevelTimeStep.C:69,88`; the two rebuild
   sums; the two per-level star-count allreduces
   (`RecordTotalStarParticleCount` + `CommunicationUpdateStarParticleCount`); the
   five-in-a-row `CommunicationAllSumValues` in
   `RecalibrateMBHFeedbackThermalRadius.C`, `RecalibrateAccretingMass.C`,
   `CalculateSubtractionParameters.C`; and `StarParticleFindAll`'s global star-list
   allgatherv (see §5.1). Two unconditional barriers bracketing `EvolveLevel`
   (`EvolveHierarchy.C:508,554`) exist purely to time the `Evtime` log and should be
   removed or guarded.

6. **Send-buffer management.** `CommunicationBufferedSend.C` scans the whole buffer
   table per send taking the *last* free slot (O(B) per message), reclaims only every
   30th call, and hard-exits at 40,000 outstanding buffers. Fix: free-list (O(1)),
   byte-budget-driven reclamation, soft in-flight cap with `MPI_Waitany`.

7. **Message-payload trims.** Halo exchanges send all ~20–30 baryon fields even when
   the consumer needs ~8 (`Grid_CommunicationSendRegion.C` with
   `SendField == ALL_FIELDS`) — a ~3× bandwidth saving available via per-call field
   masks. `int TransferSize` overflows at 2.1×10⁹ elements — within a factor of ~1.1
   for a whole-grid move of a 400³, 30-field grid, doubled by `NEW_AND_OLD`; promote to
   64-bit.

### 4.3 No threading

The only OpenMP in the tree is in `ffte4X.F90`. Everything else — hydro, gravity, CIC,
chemistry, feedback, and all the replicated-metadata loops — is flat MPI, so the
replicated hierarchy is paid once per *core* and every collective spans all cores.
**Recommendation:** thread the *grid loop* (not the kernels): a
`#pragma omp parallel for schedule(dynamic)` over the local grid list in
`EvolveLevel.C:505-587, 666-765` covers `SolveHydroEquations`, chemistry, particle
updates, and feedback with minimal code change (prerequisite: make the few `static`
scratch arrays, e.g. `Grid_FastSiblingLocatorFindSiblings.C:41`, thread-local). Running
1 rank per NUMA domain × 16–64 threads divides replicated-hierarchy memory and
collective width by the same factor — the cheapest big lever against §3.1 short of
distributing the metadata. This is precisely the path ART took in the early 2000s and
RAMSES-yOMP takes today (§6).

---

## 5. Findings: physics modules

### 5.1 Star particles and feedback (FOGGIE-specific path: `star_maker2.F` + `star_feedback6.F`)

The recurring theme: almost every star-particle code path runs once per level per
subcycle and scales with the *total* particle or grid count of the hierarchy, not with
local work. At 2^11 subcycles per root step with 10⁶–10⁷ concentrated particles, every
term is multiplied simultaneously.

1. **`StarParticleFindAll` scans the whole hierarchy every level-subcycle.**
   `StarParticleInitialize` → `StarParticleFindAll` walks **all levels
   0..MAX_DEPTH**, scanning every particle (dark matter included) to find new
   `Star`-framework objects, then allgathers the global star list (~1.3 kB per
   `StarBuffer`). For FOGGIE's configuration (no Star-framework particle types in use)
   this is ~10¹⁰ particle-type comparisons per root step **to find nothing**. The
   Active-Particle framework has exactly the right guard
   (`ActiveParticleInitialize.C:45` returns when no types are enabled);
   `StarParticleFindAll` has none. **Recommendation (highest-leverage single change in
   the star subsystem, ~hours):** skip `StarParticleFindAll` and its dependents
   entirely when no Star-framework maker (`POP3_STAR`, `STAR_CLUSTER`, `MBH_PARTICLE`,
   `SINGLE_STAR`, `StarParticleRadiativeFeedback`) is enabled. Also: never enable
   `StarParticleRadiativeFeedback` at FOGGIE star counts — it promotes every star
   particle to a `Star` object, making this path O(N²) with a
   1.3 kB × N_stars allgatherv per subcycle.

2. **Yield-table integration without an age gate.** `star_feedback6.F` performs six
   `integrate_yields` calls per star particle per subcycle (`sne_number`, `sne_mass`,
   `sne_II_metal`, `sne_Ia_metal`), each with a linear metallicity search; the
   sub-threshold skip (`nsn < nsn_timestep`) is applied only *after* all six, and there
   is no age cutoff before the calls — a 10 Gyr old particle pays the full cost to
   return zero. On level-11 grids each particle pays this ~2000× per root step.
   **Recommendations:** (i) test `age > tabAge(ntabAge)` immediately after the
   particle-type checks and skip; (ii) call `sne_number` first, test, then the rest;
   (iii) fuse the four routines onto one shared index triple; (iv) cache the
   metallicity index per particle. Also, **the code default
   `StarFeedbackSNePerTimestepLimit = 1.0e-6` contradicts the documented default of
   1e-3** — and the docs explicitly warn that values below 1e-3 significantly slow the
   code. Reconcile; production runs inheriting the code default are running 1000×
   below the recommended threshold.

3. **Per-cell waste in the feedback prologue.** In `Grid_StarParticleHandler.C`:
   - `IdentifySpeciesFields` — 12 linear field searches — is called **inside the
     innermost cell loop** of the `mu_field` construction (`:1783`, and again in the
     MOM_STAR branch at `:1719`): 300–400 wasted comparisons per cell, a 20–30×
     overhead on that loop, on every grid at every subcycle regardless of whether the
     grid has stars. Hoist it (minutes).
   - A per-cell `fopen`/`fprintf`/`fclose` diagnostic (`:1801-1829`) has all ranks
     appending to the **same file** whenever μ goes out of range — a latent
     filesystem-serializing job-killer, and redundant (the Fortran side already clamps
     μ). Remove or gate + rank-suffix.
   - `ComputeCoolingTime` — a full Grackle solve over the whole grid — runs on **every
     grid at every subcycle** whenever star formation is on (`:861-862`), but its
     output is consumed only when `StarMakerThermalCrit == 1` (`star_maker2.F:212`).
     This is a hidden near-doubling of the cooling budget. Gate it.
   - `AllocateNewParticles(0.25·N_cells)` allocates ~6 MB of particle buffers per grid
     per subcycle even when zero stars form; the full-grid temperature field, DM-field
     copy, and species→fraction round trip are likewise unconditional. Right-size and
     gate on `MakeStars`.

4. **`star_feedback6.F` internals.**
   - 14 full-grid automatic (stack) arrays per call — ~38 MB of stack for a 70³ grid at
     `-r8` — with no `-heap-arrays` in any machine file: a stack-overflow class of
     crash whose likely symptom is "memory issues" at high optimization (see item 6).
     Make them `allocatable` or add `-heap-arrays`; structurally, replace the dense
     per-cell accumulators with a compact list of explosion cells (typically ≪1% of
     cells), which also removes the unconditional full-grid zeroing and second-phase
     sweep (add an `nfbcells == 0` early return meanwhile — one line).
   - Loop-invariant `pow()` calls (real-literal exponents forcing libm `pow`) are
     recomputed per neighbor cell — e.g. the unresolved-momentum expression (3 `pow`s)
     identically 27 times per explosion cell. Hoist; use integer exponents and
     `sqrt()`.
   - The (2r+1)³ feedback stencil is traversed ~13 times per explosion cell across
     bookkeeping passes, one of which (`mom_inj`/`energy_inj`) ignores the distance
     mask. Fuse the fusable passes; the cost grows 4.6×/12.7× if
     `StarFeedbackDistRadius` is raised to 2/3.
   - `draw_stochastic` uses Fortran `RANDOM_NUMBER` with **no seed call** (unlike
     `star_maker_ssn.F`), so ranks likely draw identical sequences — correlated
     stochastic SNe across the volume. A statistics bug worth fixing alongside.
   - The two `momentum_convert` calls mutate the shared velocity arrays over the
     stencil, so overlapping feedback regions from adjacent explosion cells interact
     mid-conversion — a correctness concern in starbursts and a blocker for threading
     the explosion loop.

5. **Particle bookkeeping at scale.**
   - The particle Courant scan (`Grid_ComputeTimeStep.C:379-398`) does three passes with
     a division per particle; restructure to a single fused max-reduction and one
     division (~30× cheaper on million-particle grids).
   - Rebuild-time particle churn: every structural change reallocates and copies all
     ~11 SoA arrays; `TransferSubgridParticles` computes the cell-index assignment
     twice (CountOnly then real). Adopt capacity/size semantics with geometric growth,
     in-place compaction, and reuse of the count-pass assignment.
   - Two O(N_grids) star-count allreduces per level per subcycle
     (`RecordTotalStarParticleCount`, `CommunicationUpdateStarParticleCount`), each
     preceded by O(N_particles) type scans — merge into one, skip when nothing formed,
     maintain counts incrementally.

6. **Build-flag regression.** Commit `1bc879d` (June 2025) dropped
   `-xAVX -ip -ipo` from `Make.mach.pleiades-mpich` to work around memory issues —
   very plausibly the stack overflows of item 4. The result is scalar SSE2 code for
   all Fortran kernels on Pleiades, a plausible 1.5–3× slowdown, and no inlining of the
   small helpers called millions of times per step. **Recommendation:** fix the stack
   allocation root cause, then restore `-xAVX`/`-march` and `-ip`; reintroduce `-ipo`
   last and separately. Best effort-to-benefit ratio in the audit.

### 5.2 Gravity

1. **Exponential particle-deposit multiplicity.** `DepositParticleMassField` recurses
   over the entire subtree below each grid, so every particle is CIC-deposited into
   every ancestor level at every subcycle of that level: Σ_L 2^L ≈ 2^(L_max+1)
   deposits per particle per root step (~2×10¹⁰ CIC scatters for 10⁷ particles at
   L_max = 10), each deposit costing ~5 array passes (mass-rescale copy with per-call
   `new`/`delete`, position drift out and back, scatter).
   **Recommendations:** set `MaximumGravityRefinementLevel` explicitly (e.g.
   L_max − 2/3) — the code's one lever, and check whether FOGGIE production decks set
   it; drift once per recursion instead of per ancestor; hoist the rescale buffer;
   long term, restrict the child's *already-built* mass field to the parent
   (O(N_cells), particle-count independent — the RAMSES/Nyx approach).

2. **Per-patch Dirichlet solves — the known level-boundary force artifact.** Each
   subgrid's potential is solved independently with boundary values tri-linearly
   interpolated from the parent, iterated `PotentialIterations` (default 4) times with
   sibling-boundary averaging, plus one more solve in `EvolveLevel` — 5+ multigrid
   solves per grid per subcycle, with no coarse-fine flux matching and no composite
   solve. The buffer-region source mass is injected with **nearest-neighbor**
   interpolation (`INTERPOLATE_LINEAR` explicitly disabled in
   `Grid_CopyParentToGravitatingFieldBoundary.C:136`), the previous solution is
   discarded every subcycle (`delete`/`new` in `Grid_PreparePotentialField.C:80-83`),
   and the root-grid FFT uses the continuum −1/k² Green's function while the levels use
   a 2nd-order finite-difference Laplacian — mismatched operators at the 0/1 interface
   (`INFLUENCE2`, the discrete form, is compiled out). The multigrid itself coarsens
   non-2^j+1 dimensions with the nesting assertion commented out
   (`MultigridSolver.C:72-80`), degrading V-cycle convergence; the convergence
   tolerance *loosens* with grid size (`Grid_SolveForPotential.C:76`); and the
   unconverged fallback is up to 200 silent Gauss–Seidel sweeps per grid, then a fatal
   error.
   **Recommendations:** short term — enable `INTERPOLATE_LINEAR`, reuse the previous
   potential as the initial guess, test `INFLUENCE2`, instrument multigrid iteration
   counts (top suspect for unexplained level-time spikes), fix the tolerance to a fixed
   RMS criterion, and benchmark `PotentialIterations = 2` with capped V-cycles on later
   iterations. Long term — a composite multilevel solve (RAMSES octree multigrid /
   AMReX MLMG equivalent) is the correct fix for both cost and the well-known spurious
   heating of particles crossing refinement boundaries; at 10+ levels this accuracy
   issue is material, not cosmetic.

3. **Node-level particle kernels.** `cic_deposit.F`/`cic_interp.F` are scalar (scatter
   cannot vectorize; no directives), no OpenMP, and particles are never spatially
   sorted, so million-particle grids thrash L1 on every scatter. Morton-sorting
   particles at rebuild (~50 lines, insert after `CommunicationCollectParticles`) is
   typically worth 2–4× on deposit/interpolation. Note also
   `GravitatingMassField` sizing (+12 cells per dimension) makes the three gravity
   arrays 1.6–2.9× the baryon-field size per grid — another argument for larger
   subgrids (§3.4). Verify `CONFIG_PFLOAT_16` is *not* set: the CIC kernel bodies are
   compiled empty under 16-byte positions.

### 5.3 Timestep control

- `SetLevelTimeStep` takes a global minimum per level — one anomalous cell anywhere
  halves the level dt and doubles the subcycle count for that level *and everything
  below*. The `dtAcceleration` scan includes ghost zones
  (`Grid_ComputeTimeStep.C:412-422`), so interpolated boundary accelerations can set
  the global step; restrict it to the active region (one line). There is no cap or
  diagnostic on the subcycle ratio.
- **Recommendation (cheap, high diagnostic value):** log, once per level-step, which
  constraint on which grid set the timestep. Feedback-heated cells shrinking level dt
  is invisible today and is a plausible driver of anomalous subcycle counts.
- Chemistry (`solve_rate_cool.F`, itmax = 10⁴ per-cell substeps) creates order-of-
  magnitude within-grid cost variance that no cell-count-based balancer can see —
  supporting the measured-work balancing recommendation (§4.1.1).

### 5.4 High-frequency bookkeeping

`OutputFromEvolveLevel` polls the filesystem (4 `access()` calls + a broadcast) at
every finest-level subcycle (~1024×/root step) under `FileDirectedOutput`; `TIMER_WRITE`
issues one `MPI_Gather` per timer key per root step (`TimingCycleSkip = 1` default) and
`TIMER_START` performs two string-keyed map lookups per grid per subcycle. Throttle the
polling, set `TimingCycleSkip = 10`, pack the gathers, and use pre-registered timer
handles.

### 5.5 Instrumentation gap

The default build reports **one number** for the entire rebuild subsystem. A complete
16-slot sub-phase breakdown already exists (`RHperf` in `RebuildHierarchy.C`) but is
disabled by `#define NO_RH_PERF` at line 92. Enabling it (or converting the slots to
`TIMER_*` regions so they flow into `performance.out` and
`performance_tools.py`) is trivial and is a prerequisite for measuring everything else
in this report. Similarly, uncomment/promote the multigrid iteration-count diagnostics
(`MultigridSolver.C:217,234`) and add the per-level dt-limiter line (§5.3).

---

## 6. Comparison: RAMSES, Nyx, and ART

| Property | Enzo (this fork) | RAMSES | Nyx (AMReX) | ART/cART |
|---|---|---|---|---|
| Refinement unit | arbitrary rectangular patch | oct (2³ cells) in a fully-threaded tree | fixed-size box (`max_grid_size`, `blocking_factor`) | cell/oct in a refinement tree |
| Hierarchy metadata | fully replicated `grid` objects on every rank | distributed; each rank holds local octs + "locally essential" boundary | compact `BoxArray` + `DistributionMapping` (metadata cheap, field data local-only) | distributed local trees + buffer zone |
| Regridding | global `RebuildHierarchy` per subcycle per level | incremental, local oct refine/derefine; no global rebuild | parallel Berger–Rigoutsos tagging → new BoxArray; cheap remap | incremental, local |
| Load-balance unit & driver | whole patch; cell-count-only greedy (default) | oct segments along Peano–Hilbert curve, cost-weighted | boxes via SFC or knapsack over per-box costs; can weight particles; dual-grid option separates particle/fluid balance | SFC segments weighted by measured cost |
| Particle load balance | none (follows grid owner; particles invisible to balancer) | implicit — particles live in oct lists, follow the SFC domain | explicit — particle container redistribution; separate particle distribution maps | implicit via SFC |
| Gravity on fine levels | independent per-patch relaxation with interpolated Dirichlet BCs | multigrid solved across the octree (composite) | MLMG composite multilevel solve with coarse-fine sync/flux registers | multilevel relaxation on the tree |
| Threading / accelerators | none (flat MPI); CUDA for PPM only, off | MPI (+ RAMSES-yOMP hybrid branch; cuRAMSES GPU work) | MPI + OpenMP tiling + full GPU (CUDA/HIP) offload | hybrid MPI+OpenMP since early 2000s |
| Subcycling | factor-2 W-cycle, uncapped ratio, no alternatives | per-level dt, factor 2, incremental tree | runtime-selectable: none / standard / user pattern / **optimal** subcycling | factor-2 |

What the comparison implies for the FOGGIE pain points specifically:

- **Deep hierarchy:** RAMSES and ART never rebuild — refinement is an incremental local
  tree operation, so the 2^L × rebuild multiplier that dominates Enzo's deep levels
  simply does not exist. Nyx regrids like Enzo but the operation is a parallel
  clustering over distributed tags producing a compact box list — no per-rank
  instantiation of remote grid objects and no O(N²) loops. Enzo's equivalents are
  §3.1–3.3, and the mitigations listed there (skip-cycles, local grid lists, proxy
  metadata, neighborhood sharing) close most of the practical gap without abandoning
  patches.
- **Concentrated star particles:** in RAMSES/ART, particles are attached to the local
  tree and *automatically* follow the cost-weighted SFC decomposition — an oct-granular
  unit of balance means a single galaxy spreads across many ranks. Nyx makes particle
  cost an explicit input to the distribution map and can even balance particles on a
  different grid than the fluid. Enzo's balancer discards the particle count (§4.1.1);
  fixing that estimate, plus Hilbert mode with working keys (§4.1.3) and smaller
  effective balance granularity via larger-but-fewer grids per level, is the Enzo-shaped
  version of what the other codes get architecturally.
- **Gravity at level boundaries:** RAMSES (octree multigrid) and Nyx (MLMG with flux
  registers) both solve a composite problem across levels; Enzo's per-patch Dirichlet
  scheme (§5.2.2) is the original 1998 design and is the biggest *accuracy* gap for a
  10-level zoom, in addition to its cost.
- **Node-level:** Nyx demonstrates that block-structured AMR maps cleanly onto
  OpenMP tiling and GPUs; RAMSES-yOMP and cuRAMSES show the octree codes retrofitting
  the same hybrid parallelism Enzo lacks entirely (§4.3). The scaling walls FOGGIE hits
  are being actively attacked in every one of these communities; none of them ships a
  flat-MPI, fully-replicated-metadata configuration at FOGGIE scale.
- **Honest caveats:** octrees pay ~2× memory per cell in tree overhead and worse
  cache/vector behavior per cell than patch sweeps; AMReX's fixed box sizes waste some
  refinement volume relative to Enzo's fitted patches. Enzo's per-patch solvers are
  genuinely fast *per cell* — the deficit is in the machinery around them.

References: RAMSES — Teyssier (2002), github.com/ramses-organisation/ramses;
Nyx — Almgren et al. (2013), github.com/AMReX-Astro/Nyx, AMReX docs (MLMG, subcycling
modes, dual-grid); ART — Kravtsov, Klypin & Khokhlov (1997), Rudd, Zentner & Kravtsov
(2008) (no canonical public repo; hybrid MPI+OpenMP versions distributed by the group).

---

## 7. Recommendations, tiered

### 7.1 Tier 0 — parameter files and build flags (no code changes; try first)

| Change | Addresses | Expected effect |
|---|---|---|
| `RebuildHierarchyCycleSkip[L] = 2–4` for L ≥ 4 | §3.2 | Directly divides the dominant rebuild + collective count; highest-value single knob |
| `LoadBalancing = 4` | §4.1.2 | Escapes the mode-1 abort pathology; spatial locality preserved |
| `SubgridSizeAutoAdjust = 0`, `MinimumSubgridEdge = 16`, `MaximumSubgridSize ≈ 2.6×10⁵` | §3.4 | Ghost fraction 81%→~28%; fewer grids, fewer messages, smaller metadata |
| Set `MaximumGravityRefinementLevel = L_max − 2` (validate forces) | §5.2.1 | Truncates the exponential deposit multiplier |
| `LoadBalancingMinLevel = 1`; verify `UnigridTranspose = 2`, `ParallelRootGridIO = 1`, `ParallelParticleIO = 1` | §4.2.3, I/O | Avoids known slow paths (note the `MustRefineParticlesCreateParticles == 1` × `ParallelRootGridIO` constraint in `CosmologySimulationInitialize.C:693`) |
| `TimingCycleSkip = 10`; keep `FileDirectedOutput` polling throttled | §5.4 | Removes ~10³ filesystem ops + ~20 gathers per root step |
| Restore `-xAVX`/`-march` and `-ip` on Pleiades **after** the §5.1.4 stack fix; add `-heap-arrays 64` | §5.1.6 | 1.5–3× on all Fortran kernels |
| Verify `CONFIG_FAST_SIB = yes` and `CONFIG_PFLOAT_16` not set in production builds | §3.3, §5.2.3 | Guards against silent 100× regressions |
| Check in a FOGGIE production parameter deck under `run/` | §4.1 (L11) | Nothing in-tree records production settings today; makes tuning auditable |

### 7.2 Tier 1 — trivial-to-low-effort code fixes (minutes to hours each)

1. Guard `StarParticleFindAll` when no Star-framework type is enabled (§5.1.1) — the
   single biggest star-path win.
2. Hoist `IdentifySpeciesFields` out of the `mu_field` cell loop; skip `mu_field` on
   starless grids (§5.1.3).
3. Delete or gate the per-cell μ file write (§5.1.3).
4. Age cutoff + call reordering in `star_feedback6.F`; reconcile the
   `StarFeedbackSNePerTimestepLimit` default vs documentation (§5.1.2).
5. Gate `ComputeCoolingTime` on `StarMakerThermalCrit`; right-size
   `AllocateNewParticles` (§5.1.3).
6. Make the `star_feedback6.F` grid arrays allocatable + early return when no feedback
   (§5.1.4); seed the RNG.
7. Guard the two debug reductions in `RebuildHierarchy.C:425-426`; remove the `Evtime`
   barriers; fuse adjacent small allreduces (§4.2.5).
8. Enable `RH_PERF` (or convert to TIMER regions); add the per-level dt-limiter
   diagnostic; instrument multigrid iteration counts (§5.5).
9. Visit-stamp dedup in `FastSiblingLocatorFindSiblings`; hash/sort the
   parent-grouping loop (§3.2–3.3).
10. Lower `GRIDS_PER_LOOP` to ~512–1024; bounds-check `CommunicationReceiveIndex`
    (§4.2.2).
11. Re-enable the Hilbert bounding-box normalization (3 commented-out lines) (§4.1.3).
12. Single-pass particle Courant reduction (§5.1.5); single-pass particle-type
    histogram in `CommunicationSyncNumberOfParticles` (§4.2.4).
13. Fix the small latent bugs: the `case 21` fall-through in
    `CommunicationReceiveHandler.C:298-305` (live for any Active-Particle feedback
    user), the sink-path double free in `Grid_DepositParticlePositions.C`, the
    `CreateSUBlingList` leak, the Hilbert OOB/overflow items (§4.1), the
    `Grid_CommunicationReceiveRegion.C:422` wrong-array zeroing, and the uninitialized
    `APTotalNumber`.

### 7.3 Tier 2 — medium efforts (days to weeks)

1. Particle-aware (then measured-work) load-balance estimate via the existing
   `InputGridWork` overload (§4.1.1) — the structural fix for FOGGIE's imbalance.
2. Batch the per-star feedback-sphere allreduces and per-particle feedback-zone rounds
   (§4.2.1).
3. Level-adaptive chaining mesh (§3.3).
4. O(N) receive-handler dispatch + real batching (§4.2.2).
5. `MPI_Alltoallv` transpose; `Alltoallv` for ShareActiveParticles; sparse
   sync-of-particle-counts (§4.2.3–4.2.4).
6. Local-grid-list driving of all local-only loops in `EvolveLevel`,
   `UpdateFromFinerGrids`, `PrepareDensityField` (§3.1).
7. Assign subgrid processors before interpolation (§4.1.5).
8. Sparse must-refine flagging messages (§3.7).
9. Morton-sort particles at rebuild; drift-once deposit restructure; scratch-pool the
   PPM sweep temporaries and potential-solve `rhs` (§5.2.1, §5.2.3, §3.2).
10. Signature reuse + restored multi-dim inflection search in Berger–Rigoutsos (§3.5).
11. OpenMP over the local grid loop; 1 rank per NUMA domain (§4.3).
12. Gravity accuracy bundle: `INTERPOLATE_LINEAR`, potential reuse as initial guess,
    `INFLUENCE2` test, fixed-RMS tolerance, capped fallback (§5.2.2).

### 7.4 Tier 3 — architectural (months; sequence after Tiers 0–2 are measured)

1. **Distributed hierarchy metadata** — proxy objects first, then neighborhood-limited
   metadata exchange (§3.1). This is the change that removes the O(N_ranks) per-rank
   scaling term everywhere and is what RAMSES/ART/AMReX all have in some form.
2. **Composite multilevel gravity solve** (§5.2.2) — the correct long-term fix for both
   the solve cost and the refinement-boundary force artifacts.
3. **Pencil-decomposed or external FFT** for the root-grid solve (§4.2.3).
4. **Cross-level (multi-constraint) load balancing** including root-tile placement for
   zooms (§4.1.4).
5. Sub-communicators per level for the fine-level collectives (§2), and/or optimal
   subcycling à la Nyx.

### 7.5 Measurement plan

Before and after each tier: (1) enable the rebuild sub-phase timers and per-level dt
diagnostics (Tier 1.8); (2) record per-level, per-rank {max, mean, min} of cells,
particles, and wall time each root step (the data exists in
`OutputLevelInformation`/`CollectGridInformation`; only the reduction is missing);
(3) track the four headline numbers this audit predicts will move: wall-clock share of
`RebuildHierarchy`, time-in-collectives per root step, max/mean per-rank work on the
two finest levels, and total star-feedback time per root step. A representative
restart at production scale run for ~10 root steps per configuration is sufficient for
all of the above.

---

## Appendix A — correctness issues found in passing (fix regardless of performance)

| Issue | Location |
|---|---|
| `case 21` falls through into `case 22` — wrong-argument active-particle send on every feedback-zone receive (live for AccretingParticle/GalaxyParticle/SmartStar users) | `CommunicationReceiveHandler.C:298-305` |
| No bounds check on `CommunicationReceiveIndex` vs `MAX_RECEIVE_BUFFERS` — silent corruption past ~150k receives | 12+ call sites, e.g. `Grid_CommunicationSendRegion.C:264` |
| `int` overflow of `TransferSize` for large grid moves; `MPI_Arg` is 32-bit | `Grid_CommunicationSendRegion.C:77-78` |
| Use-after-free + double free of `ParticleMassPointerSink`; leak of `ParticleMassTemp` (sink-particle path) | `Grid_DepositParticlePositions.C:300-372, 255-257/493-495` |
| Wrong array zeroed (`BaryonField` instead of `OldBaryonField`) | `Grid_CommunicationReceiveRegion.C:422` |
| `Grids`/`ChildGrids` leak on `FluxCorrection == 0` early return | `CreateSUBlingList.C:78-83` |
| Uninitialized `APTotalNumber` accumulation | `CommunicationCollectParticles.C:221` |
| `cSndRcv` routing tables never freed; stale if decomposition changes | `CommunicationTranspose.C:339,361-362` |
| Hilbert balancer: OOB `HilbertData[-1]` read; `int` overflow of `TotalWork`; unbounded partition loops; `MinWork == 0` division | `LoadBalanceHilbertCurveRootGrids.C:146-179`, `LoadBalanceHilbertCurve.C:66,131,141-144,188,254` |
| Uninitialized read of `NumberOfCells[]` when `SubgridSizeAutoAdjust = FALSE` | `RebuildHierarchy.C:372-375` |
| Off-by-one sibling-list overflow guard (`>` vs `>=`) | `Grid_FastSiblingLocatorFindSiblings.C:120,136,170` |
| Unseeded `RANDOM_NUMBER` → correlated stochastic SNe across ranks | `star_feedback6.F:2043-2081` |
| Overlapping feedback stencils mutate shared velocity arrays mid-conversion | `star_feedback6.F:1146,1292` |
| Unbalanced `LCAPERF_START`/`STOP` on the Grackle path | `Grid_MultiSpeciesHandler.C:31,39` |
| `StarFeedbackSNePerTimestepLimit` code default (1e-6) contradicts documented default (1e-3) | `SetDefaultGlobalValues.C:628` vs `doc/manual/.../star_particles.rst:337` |
