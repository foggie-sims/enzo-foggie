#!/usr/bin/env python3
"""Simulate a work-capped Berger-Rigoutsos, offline, from a T2.1 work map.

Berger-Rigoutsos currently stops splitting on CELL count
(ProtoSubgrid_AcceptableGrid.C: size > MaximumSubgridSize).  The measured
oracle floor of ~3.35x exists because that leaves a few grids holding a
large share of a level's WORK, and a grid cannot be split across ranks.

This asks what would change if the cap were on work instead, using grids
and measured work that actually occurred, so the answer costs a few
minutes rather than an Enzo implementation.

Two effects are modelled, because balance alone is not the objective -
the buffer-zone arm had the best balance of any arm measured and was the
slowest.

  gain  the chemistry+gravity critical path (sum over levels and
        intervals of the busiest rank), recomputed after splitting
  cost  per-grid overhead - ghost zones, boundary exchange, receive
        handling - which is what defeated floor 256

The overhead is calibrated from the T0.3 scan rather than assumed: at
5 root steps and 128 ranks, floor 512 -> 256 grew L10 from 378 to 504
grids (+33%) and wall from 1670.7 s to 1775.0 s (+6.2%), i.e. roughly
0.19% of wall per 1% more grids, even though that arm balanced BETTER.
That elasticity is crude - it assumes grid counts at other levels moved
similarly - so it is reported separately from the gain rather than
folded in silently.

MODELLING ASSUMPTION, and it is optimistic: a grid whose work exceeds the
cap is split into equal pieces of equal work.  Real B-R cuts on geometry,
so the pieces will be uneven and more of them will be needed to get under
a given cap.  Treat the resulting floor as a bound, not a forecast.  A
pessimistic variant (--skew) splits into deliberately uneven pieces to
show how fast the conclusion degrades.
"""

import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np

# From the T0.3 scan: fractional wall cost per fractional increase in
# grid count, measured on an arm that balanced better and still lost.
OVERHEAD_ELASTICITY = 0.19
BENCH_WALL = 1670.7          # s, floor 512, 5 root steps, 128 ranks


def load(directory):
    recs = defaultdict(list)
    files = sorted(glob.glob(os.path.join(directory, "gridwork_rank*.txt")))
    if not files:
        sys.exit("no gridwork_rank*.txt in %s" % directory)
    for path in files:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            chem = float(f[5])
            grav = float(f[6]) if len(f) >= 14 else 0.0
            recs[(int(f[0]), int(f[2]))].append(chem + grav)
    return {k: np.array(v) for k, v in recs.items()}, len(files)


def greedy_max(weights, nranks):
    """Busiest rank under largest-first assignment; ties -> fewest grids.

    Weights are the actual work here, so this is the oracle: the best any
    assignment policy could do with these grids.  That is the right thing
    to measure, because the question is whether the GRIDS allow balance,
    not whether the estimator is good.
    """
    load = np.zeros(nranks)
    count = np.zeros(nranks)
    for w in np.sort(weights)[::-1]:
        r = int(np.lexsort((count, load))[0])
        load[r] += w
        count[r] += 1
    return load.max()


def split_to_cap(work, cap, skew=0.0):
    """Split any grid above the cap; return the new work array."""
    if cap <= 0:
        return work
    out = []
    for w in work:
        if w <= cap:
            out.append(w)
            continue
        k = int(np.ceil(w / cap))
        if skew <= 0:
            out.extend([w / k] * k)
        else:
            # uneven pieces: linearly graded, mean preserved
            g = np.linspace(1.0 - skew, 1.0 + skew, k)
            out.extend(list(w * g / g.sum()))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="run directory holding gridwork_rank*.txt")
    ap.add_argument("--skew", type=float, default=0.0,
                    help="0 = equal pieces (optimistic); 0.5 = quite uneven")
    ap.add_argument("--wall", type=float, default=BENCH_WALL,
                    help="wall clock of the baseline run, for the net estimate")
    args = ap.parse_args()

    recs, nranks = load(args.path)
    levels = sorted(set(l for _, l in recs))

    base_path = 0.0
    base_grids = 0
    base_ideal = 0.0
    for key, w in recs.items():
        if w.sum() <= 0:
            continue
        base_path += greedy_max(w, nranks)
        base_grids += len(w)
        base_ideal += w.sum() / nranks

    print("=== work-capped Berger-Rigoutsos, simulated on measured work ===")
    print("    %d ranks, %d level-intervals, %d grid-instances"
          % (nranks, len(recs), base_grids))
    print("    pieces are %s"
          % ("equal (optimistic bound)" if args.skew <= 0
             else "uneven, skew %.2f" % args.skew))
    print()
    print("    baseline critical path %8.1f s   (even share %.1f s"
          "  ->  %.2fx)"
          % (base_path, base_ideal, base_path / base_ideal))
    print()
    print("  cap    grids   vs base   crit path   floor   path saved"
          "   est. overhead    net")
    print("  " + "-" * 78)

    # caps expressed as a multiple of the per-rank even share at each
    # level-interval, so one number applies across levels of very
    # different absolute cost
    for mult in (4.0, 3.0, 2.0, 1.5, 1.0, 0.75, 0.5, 0.25):
        path = 0.0
        grids = 0
        for key, w in recs.items():
            if w.sum() <= 0:
                continue
            cap = mult * w.sum() / nranks
            nw = split_to_cap(w, cap, args.skew)
            path += greedy_max(nw, nranks)
            grids += len(nw)

        saved = base_path - path
        grid_growth = (grids - base_grids) / base_grids
        overhead = OVERHEAD_ELASTICITY * grid_growth * args.wall
        net = saved - overhead
        print("  %4.2f %8d  %+7.1f%%  %9.1f s  %6.2fx  %9.1f s  %13.1f s %+7.1f s"
              % (mult, grids, 100 * grid_growth, path, path / base_ideal,
                 saved, overhead, net))

    print()
    print("  cap is a multiple of the even share: 1.00 means no grid may")
    print("  hold more work than one rank's fair portion of its level.")
    print()
    print("  'net' is the gain minus the modelled per-grid overhead.  It is")
    print("  the number that matters: floor 256 improved balance and lost")
    print("  6.2%% of wall to overhead, and the buffer-zone arm balanced best")
    print("  of all and was slowest.  A positive net here is necessary for")
    print("  this to be worth implementing, not sufficient - the overhead")
    print("  elasticity is calibrated from a single pair of runs.")


if __name__ == "__main__":
    main()
