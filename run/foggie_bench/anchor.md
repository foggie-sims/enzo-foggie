# FOGGIE-bench anchor snapshot(s)

Record the fixed reference snapshot(s) used for all baseline/candidate
comparisons here. Fill in on first use; change only deliberately (a new
anchor invalidates comparability with prior results, so date it).

| Field | Value |
|---|---|
| Halo | *TBD (e.g. 008508 / Tempest)* |
| Snapshot | *TBD (e.g. DD1234, z ~ 2)* |
| Path on machine | *TBD (/nobackup/...)* |
| Rank count | *TBD (production count)* |
| Node model | *TBD (rom_ait / bro)* |
| Star particles | *TBD* |
| L_max reached | *TBD* |
| Anchored on | *date, by whom* |

Notes:
- Choose an output with an established galaxy (>= 1e6 star particles) so the
  star-particle machinery, must-refine flagging, and load balancer are all
  genuinely stressed.
- A second, smaller anchor (earlier snapshot or smaller halo, ~64-128 ranks)
  is useful as a cheap shakedown target before burning time at full scale.
