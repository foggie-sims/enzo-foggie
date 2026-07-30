#!/usr/bin/env python3
"""Attribute a FOGGIE-bench run's wall time across performance.out sections.

    python3 time_budget.py <run_dir> [<run_dir> ...]

Enzo's performance.out carries one block per timed cycle: a `Total` line,
one `Level_NN` line per level, and one line per named TIMER section (see
the instrumentation passes recorded in audit/status.yml T1.8 / T0.6).
Column 1 is the mean-over-ranks time for that cycle; this script sums it
over cycles and reports each section as a share of `Total`.

Two things make naive summing misleading, so they are handled here:

  - Level_NN lines are NESTED: EvolveLevel recurses, so Level_09 time is
    contained within Level_08 and so on. They are reported separately
    from named sections and never added into the attribution total.
  - Some named sections are CROSS-CUTTING: the communication pump
    (CommReceiveHandler, MPIWaitReceive, CommBufferedSend) and
    GravityAccel are entered from inside other timed sections, so their
    time is already counted in those sections. They are listed but
    excluded from the coverage sum, which therefore stays <= 100%.

Coverage = attributed / Total is the number to watch: it says how much
of the run the audit can currently name.
"""

import os
import sys

# Sections entered from within other timed sections; their time overlaps
# whatever called them, so they attribute without partitioning.
CROSS_CUTTING = {
    "CommReceiveHandler",   # invoked from SetBoundaryConditions, PrepareDensityField, ...
    "MPIWaitReceive",       # the pure MPI_Waitsome wait inside the pump
    "CommBufferedSend",     # send side, likewise called from within sections
    "GravityAccel",         # contains the per-grid SolveForPotential calls
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
        print("\n  level loops (nested; for reference only):")
        for k, v in sorted(levels.items()):
            print(f"    {v:9.1f} s  {100*v/total:5.1f}%  {k}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for rundir in sys.argv[1:]:
        report(rundir)


if __name__ == "__main__":
    main()
