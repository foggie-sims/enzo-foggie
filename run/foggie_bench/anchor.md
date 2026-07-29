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
| Restart parameter file | `/home1/jtumlins/nobackup/halo_008508/H2regulated/H2mech_tab_cont_ff/RD0016/RD0016` |
| Rank count | 512 |
| Node model | rom_ait (128 cores/node) |
| Star particles | ~1.42 million (from the T1.9 bench comparisons) |
| Anchored on | 2026-07-29, established via the T1.9 manual benches (A1/A2/B) |

## Measured properties at this anchor

- **Run-to-run noise floor** (identical code, 5 root steps, 512 ranks):
  relative mass differences grow ~1 decade per root step, 1e-11 at step 1 to
  1e-7 by step 5; refinement structure diverges by step 4; +-1 stochastic
  star formation event; wall-clock noise ~2.2%. Archived at
  `results/t19-manual/noise_floor_A1_vs_A2.json` on the `bench-results`
  branch; all production A/B gating uses it via `--noise-floor`.

## Notes

- Runs restart with `-r RD0016/RD0016` from a shadow snapshot directory
  prepared by `make_bench_run.py`; the production dump is read-only.
- The queue entry field `restart:` must point at the parameter file INSIDE
  the snapshot directory (`.../RD0016/RD0016`), not the directory itself.
