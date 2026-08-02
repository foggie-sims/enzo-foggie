# FOGGIE-bench anchor snapshot(s)

The fixed reference snapshot used for all baseline/candidate comparisons.
Change only deliberately (a new anchor invalidates comparability with prior
results, so date any change).

## Primary anchor

| Field | Value |
|---|---|
| Halo | 008508 |
| Run | H2regulated / H2mech_tab_cont_ff |
| Snapshot | RD0016 |
| Restart parameter file | `/home1/jtumlins/nobackup/halo_008508/H2regulated/H2mech_tab_cont_ff_DONE/RD0016/RD0016` |
| Rank count | 128 |
| Node model | mil_ait (1 node, 128 cores) |
| Star particles | ~1.42 million (from the T1.9 bench comparisons) |
| Anchored on | 2026-07-29, established via the T1.9 manual benches (A1/A2/B); parallel config corrected to match production 2026-07-29 |

The 1 node x 128 ranks mil_ait configuration matches the FOGGIE production
jobs. Earlier benches (T1.9 validation, t19-manual floor) ran at 512 ranks
on 4 x rom_ait - an accidental default from `make_bench_run.py`, never a
production configuration.

## Measured properties at this anchor

- **Run-to-run noise floor at 512 ranks / rom_ait (superseded config)**
  (identical code, 5 root steps): relative mass differences grow ~1 decade
  per root step, 1e-11 at step 1 to 1e-7 by step 5; refinement structure
  diverges by step 4; +-1 stochastic star formation event; wall-clock noise
  ~2.2%. Archived at `results/t19-manual/noise_floor_A1_vs_A2.json` on the
  `bench-results` branch. Valid only for gating 512-rank runs.
- **Noise floor at 1 x 128 mil_ait (production config)**: measured
  2026-07-29 from three baseline runs (t18-instrumentation-r7 A1/A2/A3,
  all pairwise comparisons merged with `merge_noise_floors.py`). Envelope:
  relative mass diffs grow ~1 decade per root step to (baryon/metal/stellar)
  2.9e-8 / 1.4e-7 / 2.6e-8 at step 5; the MESH (per-level grid/cell
  counts) stayed bit-identical in every identical-code pair through all 5
  steps; the observed structure-channel noise is integer events - +-1
  stochastic star and single particles hopping a refinement level with
  the mesh unchanged (seen from step 4 on). compare_runs.py gates the
  mesh exactly and the migration/star channels as counts. Wall-clock
  noise ~0.5%. Archived at
  `results/t18-instrumentation-r7/noise_floor_envelope_A123.json` on
  `bench-results`. **SUPERSEDED 2026-07-31 - use
  `noise_floor_envelope_A123_fixed.json` instead.** The original was
  computed with the broken `compare_runs` mass metric (no cell volume, no
  particle dx^3, refined regions double-counted; fixed in commit
  92e436d0). Corrected envelope at DD0399: baryon 8.1e-15, dark matter
  exactly 0.0, metal 7.0e-08, stellar 1.8e-08, +-1 star, 1 particle
  migration. DM matching bit-for-bit between identical-code runs is the
  proof the correction is right. C0/C1 gating is now far sharper; C2
  items change the answer by design and will legitimately exceed the
  floor, so treat it there as a reference scale, not pass/fail.
  Single-pair floors
  scatter ~5x at step 5; always gate against a multi-pair envelope. A noise
  floor is config-specific: gate runs only with a floor measured at the
  same rank count and node model.


## Anchor fragility (2026-08-02)

The anchor snapshot lives inside a **user-managed production directory**,
not a protected location. It has already moved once: the parent run was
renamed `H2mech_tab_cont_ff` -> `H2mech_tab_cont_ff_DONE`, which broke
every bench submission until the path was updated here.

If that directory is ever deleted, **every prior bench result loses its
baseline** and cross-comparability with the whole audit is gone. The
snapshot is 8.6 GB. Consider preserving a copy under
`foggie_bench_root/` (or another location not subject to production
housekeeping) before relying on it for further work.

## Notes

- Runs restart with `-r RD0016/RD0016` from a shadow snapshot directory
  prepared by `make_bench_run.py`; the production dump is read-only.
- The queue entry field `restart:` must point at the parameter file INSIDE
  the snapshot directory (`.../RD0016/RD0016`), not the directory itself.
