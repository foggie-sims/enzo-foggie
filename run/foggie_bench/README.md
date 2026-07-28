# FOGGIE-bench: abbreviated production runs for validating enzo-performance changes

This directory implements the "B3" benchmark of the performance test plan
(`audit/enzotestplan.html`, sections 2.1-2.2): a **short restart window from a
real FOGGIE production snapshot**, run twice (baseline commit vs. candidate
commit) and compared. CI on the `enzo-performance` branch proves changes
correct on small idealized problems; this harness is the check that catches
what CI structurally cannot — deep nesting, must-refine particle machinery,
thousands of grids per level, load-balancer behavior, and particle churn at
production scale.

**Status: harness written in advance of machine access; every script prints
what it is doing and fails loudly. Expect to shake it down on the first real
restart — treat the first baseline/candidate pair as a trial of the harness
itself.**

## The protocol

1. **Pick the anchor snapshot once** and record it here (halo, snapshot name,
   filesystem path, rank count). Prefer a z ~ 2 output of a halo with an
   established galaxy (>= 1e6 star particles, L_max reached). This snapshot is
   the fixed reference for all A/B comparisons until deliberately re-anchored.

2. **Baseline run:** build the *baseline* commit (normally the tip of
   `enzo-performance` before your change), then:

       python3 make_bench_run.py --restart /path/to/DD1234/DD1234 \
           --enzo /path/to/baseline/enzo.exe --nsteps 5 --out bench_baseline

   This creates `bench_baseline/` with an edited copy of the restart parameter
   file (StopCycle = CycleNumber + nsteps, per-root-step data dumps on, dt-based
   dumps off so both runs dump at identical cycles) and a `launch.sh` from the
   PBS template. Submit it.

3. **Candidate run:** same, with the candidate build and `--out bench_candidate`.
   Same rank count, same node type, same PBS resources — vary nothing but the
   code.

4. **Compare:**

       python3 compare_runs.py bench_baseline bench_candidate --class C0

   `--class` selects the gate per the test plan's change classes:
   - `C0`: bitwise field/particle identity expected (h5diff over dumps)
     *plus* everything below. NOTE: per ledger item T1.14, parallel Enzo is
     not run-to-run reproducible, so C0 bitwise A/B is only meaningful when
     baseline and candidate produce deterministic outputs at these settings —
     in practice C0 at production scale means "compare the invariants and
     statistics, expect them at machine precision", not literal bitwise.
   - `C1`: conservation invariants at machine precision; global sums and
     per-level statistics within tight tolerances; timing report informational.
   - `C2`: report everything; gates only the conservation invariants; the
     physics comparison is a human review against the seed-ensemble bands
     (test plan section 3.3) once those exist.

   The comparator emits `comparison.json` plus a human-readable verdict.

5. **Record the result** in `audit/status.yml` on the relevant item (`tests:`
   entry with the comparison verdict/location), so the dashboard reflects
   HPC-validated status, not just CI status.

## What the comparator checks

From each run directory (parsed from the dumps Enzo wrote):

- **Conservation invariants** (gate for every class): total baryon mass,
  total metal mass, total stellar mass + star particle count, total DM mass,
  per dump. Baseline vs candidate at matching cycles.
- **Per-level structure**: grid count, cell count, mean grid size, ghost-zone
  fraction, particle count per level (from the `.hierarchy` files). C1
  tolerance: identical unless the change is supposed to alter gridding.
- **Global field sums**: sum(rho*dV), sum(E*rho*dV), momentum components,
  metal mass by field — grid-aware, level-masked (finest data wins where
  levels overlap).
- **Timing KPIs** (informational, from `performance.out` if EnzoTiming is on):
  wall per root step, RebuildHierarchy share, per-level times. This is where
  the audit's predicted speedups (e.g. T1.9's O(N^2) -> O(N) rebuild loops)
  should appear at production grid counts.
- **Run health**: both runs reached StopCycle; no ENZO_FAIL in logs.

## Files

| File | Purpose |
|---|---|
| `make_bench_run.py` | Prepare an abbreviated-run directory from a restart |
| `compare_runs.py` | A/B comparison and gate |
| `pbs_pleiades.template.sh` | PBS script template (modules per `src/qsub_compile_enzo.sh`) |
| `anchor.md` | Record of the chosen anchor snapshot(s) — fill in on first use |

## Caveats and known limitations

- Requires `h5py` for field sums (module load python3 with h5py on Pleiades,
  or a conda env); structure-only comparisons work without it.
- `make_bench_run.py` edits a *copy* of the restart parameter file and never
  touches the production output directory; the restart data itself is
  referenced in place (Enzo reads it read-only).
- Wall-clock comparisons need same-node-type, ideally same-reservation runs;
  treat < 10% differences as noise until repeats say otherwise (test plan 1.3).
- The dt/output cadence override assumes cycle-based dumping works on the
  production parameter set; if the production deck drives outputs by redshift
  the harness's edits may need extending — first-use shakedown item.
