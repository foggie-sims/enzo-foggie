#!/usr/bin/env python3
"""Report how evenly an Enzo run is spread across ranks.

Works on ANY Enzo dump written with grouped HDF5 output, including runs
finished long ago, because rank ownership is already recorded: each grid's
data lives in the .cpuNNNN file of the rank that owned it.  Nothing has to
be instrumented, and nothing has to be enabled at run time.

For each dump it reports, per rank and per level, the quantities the
balancer does and does not weigh by:

  cells      what CommunicationLoadBalanceGrids actually equalises
  particles  what it does not - dark matter plus stars
  stars      the subset that drives feedback work
  grids      how many indivisible atoms the rank was handed

and, if a T2.1 work map (gridwork_rank*.txt) is present alongside, the
MEASURED chemistry and particle-gravity time per rank, which is the only
one of these that is actual work rather than a proxy for it.

The headline number is max/mean.  A rank cannot start the next
synchronised phase until the slowest finishes, so max/mean is what the
barriers charge, and 1.00 is perfect.

Usage:
    load_balance.py DD0414                  # one dump
    load_balance.py rundir                  # every dump in a run
    load_balance.py rundir --csv out.csv    # machine-readable, for trends
    load_balance.py rundir --watch          # poll for new dumps as they land

Comparing two runs is just running it on both: the per-level max/mean
columns are directly comparable as long as the rank count matches.
"""

import argparse
import csv
import glob
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("load_balance.py needs h5py")


# --------------------------------------------------------------------------
# hierarchy: which level each grid is on
# --------------------------------------------------------------------------

def grid_levels(hpath):
    """{group_name: level} from a .hierarchy file.

    The HDF5 grid groups carry no attributes in this Enzo build, so level
    has to come from the hierarchy's pointer graph rather than the data.
    """
    text = open(hpath).read()
    nxt = re.findall(r"Pointer:\s*Grid\[(\d+)\]->NextGridNextLevel\s*=\s*(\d+)",
                     text)
    same = re.findall(r"Pointer:\s*Grid\[(\d+)\]->NextGridThisLevel\s*=\s*(\d+)",
                      text)
    level_of = {1: 0}
    links = ([(int(a), int(b), 1) for a, b in nxt if int(b)] +
             [(int(a), int(b), 0) for a, b in same if int(b)])
    changed = True
    while changed:                      # pointers are not in level order
        changed = False
        for a, b, dl in links:
            if a in level_of and b not in level_of:
                level_of[b] = level_of[a] + dl
                changed = True
    return {"Grid%08d" % g: l for g, l in level_of.items()}


# --------------------------------------------------------------------------
# one dump
# --------------------------------------------------------------------------

def scan_dump(dumppath):
    """Per-(rank, level) totals for one dump."""
    base = os.path.basename(dumppath.rstrip("/"))
    hpath = os.path.join(dumppath, base) + ".hierarchy"
    if not os.path.exists(hpath):
        hpath = dumppath + ".hierarchy"
    levels = grid_levels(hpath) if os.path.exists(hpath) else {}

    files = sorted(glob.glob(os.path.join(dumppath, "*.cpu*")))
    if not files:
        files = sorted(glob.glob(dumppath + ".cpu*"))
    if not files:
        return None

    stats = defaultdict(lambda: defaultdict(float))
    for path in files:
        m = re.search(r"\.cpu(\d+)$", path)
        if not m:
            continue
        rank = int(m.group(1))
        with h5py.File(path, "r") as h:
            for gname in h:
                # each file also carries a Metadata group, which is not a grid
                if not re.fullmatch(r"Grid\d{8}", gname):
                    continue
                g = h[gname]
                level = levels.get(gname, -1)
                key = (rank, level)

                # active cells: field arrays hold active zones only
                dens = g.get("Density")
                if dens is not None:
                    stats[key]["cells"] += float(np.prod(dens.shape))

                ptype = g.get("particle_type")
                if ptype is not None:
                    t = ptype[()]
                    stats[key]["particles"] += float(t.size)
                    stats[key]["stars"] += float((t == 2).sum())

                stats[key]["grids"] += 1.0
    return stats


def add_work_map(stats, rundir):
    """Fold in measured work per rank, if a T2.1 work map is present."""
    files = sorted(glob.glob(os.path.join(rundir, "gridwork_rank*.txt")))
    if not files:
        return False
    for path in files:
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # seq cycle level rank cells chem [grav] subcycles ...
            rank, level = int(f[3]), int(f[2])
            stats[(rank, level)]["chem_s"] += float(f[5])
            if len(f) >= 14:                     # gravity column present
                stats[(rank, level)]["grav_s"] += float(f[6])
    return True


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def imbalance(values, nranks):
    """max/mean over ranks, counting idle ranks as the zeros they are."""
    v = np.zeros(nranks)
    v[:len(values)] = values
    return v.max() / v.mean() if v.mean() > 0 else float("nan")


def report(stats, label, nranks=None, quiet=False):
    ranks = sorted(set(r for r, _ in stats))
    levels = sorted(set(l for _, l in stats))
    nranks = nranks or (max(ranks) + 1 if ranks else 0)
    metrics = ["cells", "particles", "stars", "grids"]
    if any("chem_s" in stats[k] for k in stats):
        metrics.append("chem_s")
    if any("grav_s" in stats[k] for k in stats):
        metrics.append("grav_s")

    rows = []
    if not quiet:
        print()
        print("=== %s ===" % label)
        print("    %d ranks; imbalance = max/mean over ranks (1.00 is perfect)"
              % nranks)
        print()
        print("  lvl  grids" + "".join("%12s" % m for m in metrics))
        print("  " + "-" * (11 + 12 * len(metrics)))

    for level in levels + ["ALL"]:
        per = {m: np.zeros(nranks) for m in metrics}
        ngrids = 0.0
        for (r, l), d in stats.items():
            if level != "ALL" and l != level:
                continue
            if r >= nranks:
                continue
            for m in metrics:
                per[m][r] += d.get(m, 0.0)
            ngrids += d.get("grids", 0.0)
        if ngrids == 0:
            continue
        vals = {m: (per[m].max() / per[m].mean() if per[m].mean() > 0
                    else float("nan")) for m in metrics}
        if not quiet:
            # a nan here means the quantity is absent at this level (no
            # stars yet, say), not that the balance is unknown
            cells_ = ["%12.2f" % vals[m] if vals[m] == vals[m] else "%12s" % "-"
                      for m in metrics]
            print("  %3s %6d" % (str(level), ngrids) + "".join(cells_))
        row = dict(label=label, level=level, grids=int(ngrids))
        row.update({("%s_maxmean" % m): vals[m] for m in metrics})
        row.update({("%s_total" % m): float(per[m].sum()) for m in metrics})
        rows.append(row)

    if not quiet:
        # idle ranks are invisible in max/mean but matter a lot
        allc = np.zeros(nranks)
        for (r, l), d in stats.items():
            if r < nranks:
                allc[r] += d.get("cells", 0.0)
        idle = int((allc == 0).sum())
        if idle:
            print()
            print("  NOTE: %d of %d ranks hold no cells at all" % (idle, nranks))
    return rows


# --------------------------------------------------------------------------

def per_rank(stats, metric, nranks):
    v = np.zeros(nranks)
    for (r, l), d in stats.items():
        if r < nranks:
            v[r] += d.get(metric, 0.0)
    return v


def plot(runs, outpath, nranks, sort=False):
    """Per-rank distribution for each metric, one panel each.

    runs is [(label, stats), ...]; several overlay for comparison.
    Ranks are plotted in index order by default, because that is what the
    balancer sees; --sort reorders each series independently by descending
    load, which makes the shape of the imbalance easier to read but
    breaks the correspondence between series at a given x.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [("cells", "active cells"),
               ("particles", "particles (DM + stars)"),
               ("stars", "star particles")]
    if any("chem_s" in s[k] for _, s in runs for k in s):
        metrics.append(("chem_s", "measured chemistry (s)"))
    if any("grav_s" in s[k] for _, s in runs for k in s):
        metrics.append(("grav_s", "measured particle gravity (s)"))

    fig, axes = plt.subplots(len(metrics), 1,
                             figsize=(11, 2.6 * len(metrics)), sharex=True)
    if len(metrics) == 1:
        axes = [axes]
    colours = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]

    def fmt(x):
        if x >= 1e6:
            return "%.2fM" % (x / 1e6)
        if x >= 1e3:
            return "%.1fk" % (x / 1e3)
        return "%.3g" % x

    for ax, (m, title) in zip(axes, metrics):
        notes = []
        for i, (label, stats) in enumerate(runs):
            v = per_rank(stats, m, nranks)
            if v.sum() <= 0:
                continue
            y = np.sort(v)[::-1] if sort else v
            mean, med = v.mean(), np.median(v)
            ratio = v.max() / mean if mean > 0 else float("nan")
            c = colours[i % 4]
            ax.plot(np.arange(nranks), y, lw=1.0, color=c,
                    label="%s  (max/mean %.2f)" % (label, ratio))
            ax.axhline(mean, color=c, ls=":", lw=0.9, alpha=0.8)
            ax.axhline(med, color=c, ls="--", lw=0.9, alpha=0.5)

            # mean vs median says WHERE the imbalance is: mean well above
            # median means a few heavy ranks; below means a few idle ones
            idle = int((v == 0).sum())
            cv = v.std() / mean if mean > 0 else float("nan")
            notes.append(
                "%s\n  mean %s   median %s\n  min %s   max %s\n"
                "  max/mean %.2f   max/median %.2f\n  CV %.2f   idle ranks %d"
                % (label, fmt(mean), fmt(med), fmt(v.min()), fmt(v.max()),
                   ratio, (v.max() / med if med > 0 else float("nan")),
                   cv, idle))

        if notes:
            ax.text(0.005, 0.97, "\n".join(notes), transform=ax.transAxes,
                    fontsize=6.6, va="top", ha="left", family="monospace",
                    bbox=dict(boxstyle="round,pad=0.35", fc="white",
                              ec="0.7", alpha=0.85))
        ax.set_ylabel(title, fontsize=9)
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(alpha=0.25, lw=0.5)
        ax.margins(x=0.01)

    axes[-1].set_xlabel("MPI rank" + (" (sorted by load)" if sort else ""))
    axes[0].set_title("Per-rank load balance   -   dotted line = mean "
                      "(the even share),  dashed = median", fontsize=11)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    print("  wrote %s" % outpath)


def find_dumps(path):
    if os.path.isdir(path) and glob.glob(os.path.join(path, "*.cpu*")):
        return [path]
    out = []
    for pat in ("DD????", "RD????"):
        out += sorted(glob.glob(os.path.join(path, pat)))
    return [d for d in out if glob.glob(os.path.join(d, "*.cpu*"))]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", help="a dump directory, or a run directory")
    ap.add_argument("--csv", help="append per-level rows here, for trends")
    ap.add_argument("--ranks", type=int,
                    help="expected rank count (so idle ranks are counted)")
    ap.add_argument("--plot", metavar="PNG",
                    help="write a per-rank plot here")
    ap.add_argument("--compare", nargs="*", default=[], metavar="DUMP",
                    help="overlay these dumps on the plot, for A/B")
    ap.add_argument("--sort", action="store_true",
                    help="sort ranks by load in the plot (shape, not identity)")
    ap.add_argument("--watch", action="store_true",
                    help="poll for new dumps and report each as it lands")
    ap.add_argument("--interval", type=float, default=60.0,
                    help="seconds between polls with --watch (default 60)")
    args = ap.parse_args()

    if args.plot:
        runs = []
        for p in [args.path] + args.compare:
            for d in find_dumps(p):
                s = scan_dump(d)
                if not s:
                    continue
                add_work_map(s, os.path.dirname(d.rstrip("/")) or ".")
                # label by run directory when comparing arms, since the
                # dump name is usually identical across them
                parent = os.path.basename(os.path.dirname(d.rstrip("/")))
                runs.append(("%s/%s" % (parent, os.path.basename(d.rstrip("/"))), s))
                report(s, runs[-1][0], nranks=args.ranks)
        if not runs:
            sys.exit("no dumps found to plot")
        nranks = args.ranks or max(
            (r for _, s in runs for r, _ in s), default=0) + 1
        plot(runs, args.plot, nranks, sort=args.sort)
        return

    seen = set()
    allrows = []
    while True:
        dumps = [d for d in find_dumps(args.path) if d not in seen]
        for d in dumps:
            seen.add(d)
            stats = scan_dump(d)
            if not stats:
                continue
            add_work_map(stats, os.path.dirname(d.rstrip("/")) or ".")
            allrows += report(stats, os.path.basename(d.rstrip("/")),
                              nranks=args.ranks)
        if args.csv and allrows:
            new = not os.path.exists(args.csv)
            with open(args.csv, "a", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=sorted(allrows[0]))
                if new:
                    w.writeheader()
                w.writerows(allrows)
            print("\n  wrote %d rows to %s" % (len(allrows), args.csv))
            allrows = []
        if not args.watch:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
