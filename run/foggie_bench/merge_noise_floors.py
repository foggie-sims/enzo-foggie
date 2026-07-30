#!/usr/bin/env python3
"""Merge several identical-code comparison.json files into one envelope floor.

    python3 merge_noise_floors.py floor_A1_vs_A2.json floor_A1_vs_A3.json \
        floor_A2_vs_A3.json --out noise_floor_envelope.json

A noise floor measured from a single identical-code pair is one sample of a
chaotically growing divergence and can be unluckily small at any given dump
(observed in t18-instrumentation-r7: step-5 baryon diffs of 5.8e-9 vs 2.3e-8
between pairs of the same three baseline runs). With N baseline runs, all
N(N-1)/2 pairwise comparisons are floor samples; this script takes, per dump
and per metric, the maximum relative difference across the pairs, the maximum
star-count delta, and the earliest dump at which any pair's refinement
structure diverged. The result is a valid --noise-floor input for
compare_runs.py, gating against the measured envelope rather than one draw.
"""

import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("floors", nargs="+",
                    help="comparison.json files from identical-code pairs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if len(args.floors) < 2:
        print("note: merging a single floor is an identity operation")

    pairs = [json.load(open(p)) for p in args.floors]
    merged = {
        "class": pairs[0].get("class"),
        "merged_from": args.floors,
        "dumps_compared": sorted(set().union(*[p.get("dumps_compared", [])
                                               for p in pairs])),
        "failures": [],
        "per_dump": {},
    }

    # Earliest dump at which any pair's per-level structure diverged; keep a
    # single representative failure line (compare_runs.py parses the dump
    # name off the front of that message).
    struct_break = None
    for p in pairs:
        for msg in p.get("failures", []):
            if "per-level structure differs" in msg:
                d = msg.split(":")[0].strip()
                if struct_break is None or d < struct_break:
                    struct_break = d
    if struct_break is not None:
        merged["failures"].append(
            f"{struct_break}: per-level structure differs")

    for dump in merged["dumps_compared"]:
        entries = [p["per_dump"][dump] for p in pairs
                   if dump in p.get("per_dump", {})]
        out = {}
        # Envelope of every rel_* metric present in any pair.
        for e in entries:
            for k, v in e.items():
                if k.startswith("rel_"):
                    out[k] = max(out.get(k, 0.0), v)
        # Star-count floor: compare_runs.py reconstructs the delta from the
        # sums_* blocks, so copy them from the pair with the largest delta.
        best = None
        for e in entries:
            sa, sb = e.get("sums_baseline"), e.get("sums_candidate")
            if sa and sb:
                d = abs(sa["star_count"] - sb["star_count"])
                if best is None or d > best[0]:
                    best = (d, sa, sb)
        if best:
            out["star_count_diff"] = best[0]
            out["sums_baseline"] = best[1]
            out["sums_candidate"] = best[2]
        merged["per_dump"][dump] = out

    with open(args.out, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"Merged {len(pairs)} pair floors -> {args.out}")
    for dump in merged["dumps_compared"]:
        e = merged["per_dump"][dump]
        rels = "  ".join(f"{k[4:]}={v:.3e}" for k, v in sorted(e.items())
                         if k.startswith("rel_"))
        print(f"  {dump}: {rels}  dstar={e.get('star_count_diff', 0)}")
    if struct_break:
        print(f"  structure diverges from {struct_break}")


if __name__ == "__main__":
    main()
