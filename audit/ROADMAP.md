# FOGGIE/Enzo performance roadmap

**Status as of 2026-07-31.** This is the *measured* plan. The original
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

### 5.1 T2.10 — forced bisection in Berger-Rigoutsos *(highest expected value)*
A cap of 400 should yield ~2117 grids at L10; only 436 appear.
`ProtoSubgrid_AcceptableGrid` rejects oversized grids, but Berger-Rigoutsos
splits on signature zero-crossings and a **solidly flagged region** — the
galaxy centre — offers none, so it returns one oversized grid regardless
of the cap. Adding a geometric bisection fallback when a grid exceeds the
cap would unlock the granularity we are currently asking for and not
getting. Given 195→436 grids bought 25.8%, this is the clearest remaining
lever.

### 5.2 `SPFinalize` — skip the per-subcycle collective
Full analysis and design in **`SPFinalize_Edits.md`**. Still the largest
single line item (62%), but **re-measure first**: T0.3 already cut the
barrier 43% without touching this code, because the barrier is a
*symptom* of granularity. Its remaining value shrinks as §5.1 lands, and
it is the riskiest change on the board — a rank-divergent skip condition
**hangs** the run, a failure mode T2.1 demonstrated for real.

### 5.3 Measured-work load balancing (T2.1's remaining half)
Per-grid elapsed time from the previous cycle, via the existing
`LoadBalanceHilbertCurve(grid*[], int, float *InputGridWork, int*)`
overload. Captures chemistry variation without guessing a ratio. Sequence
**after** §5.1 — it is bounded by granularity until grids can actually be
subdivided. The rank-agreement `Allreduce` added for T2.1 is a
prerequisite for any data-dependent weighting and is already in place.

### 5.4 Smaller / unsequenced
- **T0.7** compiler flags (`-march`, `-ip`) — unblocked now that T1.6 moved
  the big arrays off the stack; could move the compute-bound remainder.
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
