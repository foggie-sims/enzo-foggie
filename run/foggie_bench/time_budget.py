#!/usr/bin/env python3
"""Attribute a FOGGIE-bench run's wall time across performance.out sections.

    python3 time_budget.py <run_dir> [<run_dir> ...]

Enzo's performance.out carries one block per timed cycle: a `Total` line,
one `Level_NN` line per level, and one line per named TIMER section (see
the instrumentation passes recorded in audit/status.yml T1.8 / T0.6).
Column 1 is the mean-over-ranks time for that cycle; this script sums it
over cycles and reports each section as a share of `Total`.

Two things make naive summing misleading, so they are handled here:

  - Level_NN lines are a SEPARATE DECOMPOSITION AXIS, not part of the
    section budget. EvolveLevel's level timer brackets stop before the
    recursive EvolveLevel(level+1) call, so each Level_NN is the
    EXCLUSIVE per-level work (verified: their sum is ~87% of Total, the
    remainder being RebuildHierarchy and I/O, which sit outside the
    brackets). They decompose the same time by level that the named
    sections decompose by phase, so the two must never be added
    together.
  - Some named sections are CROSS-CUTTING: they are entered from inside
    another timed section, so their time is already counted there.
    Verified call sites: SolveForPotential is started inside
    PrepareDensityField.C, and the communication pump
    (CommReceiveHandler / MPIWaitReceive / CommBufferedSend) is entered
    from SetBoundaryConditions, PrepareDensityField, UpdateFromFinerGrids
    and others. These are listed but excluded from the coverage sum.

Coverage = attributed / Total is the number to watch: it says how much
of the run the audit can currently name. Coverage meaningfully above
100% means a section pair overlaps that is not yet in CROSS_CUTTING.
"""

import os
import sys

# Sections entered from within other timed sections; their time overlaps
# whatever called them, so they attribute without partitioning.
CROSS_CUTTING = {
    "CommReceiveHandler",   # invoked from SetBoundaryConditions, PrepareDensityField, ...
    "MPIWaitReceive",       # the pure MPI_Waitsome wait inside the pump
    "CommBufferedSend",     # send side, likewise called from within sections
    "SolveForPotential",    # started inside PrepareDensityField.C:414
}


def read_sections(path):
    """Sum column 1 per section label over every cycle block."""
    totals = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2 or parts[0] == "Cycle_Number":
                continue
            try:
                totals[parts[0]] = totals.get(parts[0], 0.0) + float(parts[1])
            except ValueError:
                continue
    return totals


def report(rundir):
    path = os.path.join(rundir, "performance.out")
    if not os.path.isfile(path):
        print(f"{rundir}: no performance.out")
        return
    totals = read_sections(path)
    total = totals.pop("Total", 0.0)
    if total <= 0:
        print(f"{rundir}: no Total line")
        return

    levels = {k: v for k, v in totals.items() if k.startswith("Level_")}
    named = {k: v for k, v in totals.items() if not k.startswith("Level_")}
    exclusive = {k: v for k, v in named.items() if k not in CROSS_CUTTING}
    crossing = {k: v for k, v in named.items() if k in CROSS_CUTTING}
    attributed = sum(exclusive.values())

    print(f"\n=== {rundir} ===")
    print(f"Total wall: {total:.1f} s")
    print("\n  attributed sections (exclusive):")
    for k, v in sorted(exclusive.items(), key=lambda kv: -kv[1]):
        print(f"    {v:9.1f} s  {100*v/total:5.1f}%  {k}")
    print(f"    {attributed:9.1f} s  {100*attributed/total:5.1f}%  ATTRIBUTED")
    print(f"    {total-attributed:9.1f} s  {100*(total-attributed)/total:5.1f}%  unattributed")

    if crossing:
        print("\n  cross-cutting (already inside the above):")
        for k, v in sorted(crossing.items(), key=lambda kv: -kv[1]):
            print(f"    {v:9.1f} s  {100*v/total:5.1f}%  {k}")

    if levels:
        lsum = sum(levels.values())
        print("\n  per-level decomposition (exclusive of recursion;"
              " a separate axis - do NOT add to the sections above):")
        for k, v in sorted(levels.items()):
            print(f"    {v:9.1f} s  {100*v/total:5.1f}%  {k}")
        print(f"    {lsum:9.1f} s  {100*lsum/total:5.1f}%  sum of levels"
              " (remainder is rebuild + I/O, outside the level brackets)")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for rundir in sys.argv[1:]:
        report(rundir)


if __name__ == "__main__":
    main()
