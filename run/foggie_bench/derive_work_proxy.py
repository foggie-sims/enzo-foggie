#!/usr/bin/env python3
"""Find a LOCAL predictor of per-grid work (audit T2.1 / work-capped B-R).

A work-capped Berger-Rigoutsos needs to know what a candidate subgrid will
cost at the moment it decides where to cut.  Measured work lives on the
outgoing subgrids, which are scattered across ranks - only 9 of 350 L10
grids sit on the rank owning their own root tile - so using it needs a
rank-agreement Allreduce.  That is affordable but not free: the equivalent
Allreduce in the particle-weighting experiment measured ~2.6% of wall.

A predictor computed from fields the parent already holds needs no
communication at all, and cannot produce rank-divergent grid structure,
which is a hard failure (MPI mismatch) rather than a slow one.

This asks how good such a predictor can be, by regressing MEASURED work
from a T2.1 work map against quantities available locally in the same
run's dump.  Physically the leading candidate is the cooling rate: Grackle
subcycles a cell roughly dt/t_cool times, so sum(1/t_cool) over a grid
should track its chemistry cost far better than the cell count does.

Scored by correlation in the log, because the decision the cap makes is
"is this grid several times more expensive than a fair share", which is a
ratio question, not an absolute one.
"""

import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("needs h5py")


def hierarchy_geometry(hpath):
    """{group_name: (level, left_edge_tuple)} from a .hierarchy file."""
    text = open(hpath).read()
    nxt = re.findall(r"Pointer:\s*Grid\[(\d+)\]->NextGridNextLevel\s*=\s*(\d+)", text)
    same = re.findall(r"Pointer:\s*Grid\[(\d+)\]->NextGridThisLevel\s*=\s*(\d+)", text)
    level_of = {1: 0}
    links = ([(int(a), int(b), 1) for a, b in nxt if int(b)] +
             [(int(a), int(b), 0) for a, b in same if int(b)])
    changed = True
    while changed:
        changed = False
        for a, b, dl in links:
            if a in level_of and b not in level_of:
                level_of[b] = level_of[a] + dl
                changed = True

    out = {}
    for block in re.split(r"\n(?=Grid = )", text):
        m = re.match(r"Grid = (\d+)", block)
        if not m:
            continue
        gid = int(m.group(1))
        mm = re.search(r"GridLeftEdge\s*=\s*([-\d\. eE+]+)", block)
        if not mm:
            continue
        le = tuple(round(float(x), 9) for x in mm.group(1).split()[:3])
        out["Grid%08d" % gid] = (level_of.get(gid, 0), le)
    return out


def load_work(rundir):
    """{(level, left_edge): mean work rate} from the work map."""
    acc = defaultdict(list)
    for path in glob.glob(os.path.join(rundir, "gridwork_rank*.txt")):
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            level = int(f[2])
            chem = float(f[5])
            grav = float(f[6]) if len(f) >= 14 else 0.0
            # 14-col (with gravity): ... chem[5] grav[6] subcycles[7] xl[8]
            # 13-col (chemistry only): ... chem[5] subcycles[6] xl[7]
            sub = max(int(f[7]) if len(f) >= 14 else int(f[6]), 1)
            le = tuple(round(float(x), 9)
                       for x in (f[8:11] if len(f) >= 14 else f[7:10]))
            acc[(level, le)].append((chem + grav) / sub)
    return {k: float(np.mean(v)) for k, v in acc.items()}


def grid_predictors(g):
    """Local, communication-free quantities for one grid."""
    d = g["Density"][()].astype(np.float64)
    cells = float(d.size)
    out = {"cells": cells, "mass": float(d.sum())}

    t = g.get("Temperature")
    ct = g.get("Cooling_Time")
    e = g.get("Electron_Density")

    if ct is not None:
        c = np.abs(ct[()].astype(np.float64))
        c = np.where(c > 0, c, np.inf)
        # Grackle subcycles ~ dt/t_cool, so the cooling RATE summed over
        # cells is the physically motivated cost estimate
        out["coolrate"] = float((1.0 / c).sum())
        out["coolrate_d"] = float((d / c).sum())
    if t is not None:
        tt = t[()].astype(np.float64)
        out["invT"] = float((1.0 / np.maximum(tt, 1.0)).sum())
        out["logT"] = float(np.log10(np.maximum(tt, 1.0)).sum())
    if e is not None:
        out["electron"] = float(e[()].astype(np.float64).sum())

    pt = g.get("particle_type")
    if pt is not None:
        p = pt[()]
        out["particles"] = float(p.size)
        out["stars"] = float((p == 2).sum())
    else:
        out["particles"] = 0.0
        out["stars"] = 0.0
    return out


def main():
    rundir = sys.argv[1] if len(sys.argv) > 1 else "."
    dumps = sorted(glob.glob(os.path.join(rundir, "DD????")))
    if not dumps:
        sys.exit("no dumps in %s" % rundir)
    dump = dumps[-1]
    base = os.path.basename(dump)

    work = load_work(rundir)
    geom = hierarchy_geometry(os.path.join(dump, base + ".hierarchy"))
    print("  work-map entries %d, hierarchy grids %d" % (len(work), len(geom)))

    rows = []
    for path in sorted(glob.glob(os.path.join(dump, "*.cpu*"))):
        with h5py.File(path, "r") as h:
            for gname in h:
                if not re.fullmatch(r"Grid\d{8}", gname):
                    continue
                key = geom.get(gname)
                if key is None or key not in work:
                    continue
                p = grid_predictors(h[gname])
                p["level"] = key[0]
                p["work"] = work[key]
                rows.append(p)

    if len(rows) < 50:
        sys.exit("only %d grids matched work records - cannot regress" % len(rows))

    keys = [k for k in rows[0] if k not in ("work", "level")]
    w = np.array([r["work"] for r in rows])
    lev = np.array([r["level"] for r in rows])
    good = w > 0
    print("  matched %d grids (%d with non-zero work)" % (len(rows), good.sum()))
    print()
    print("  predictor        corr(log)   corr(log, L>=8)   ratio spread")
    print("  " + "-" * 62)

    lw = np.log10(w[good])
    deep = lev[good] >= 8
    scored = []
    for k in keys:
        v = np.array([r.get(k, 0.0) for r in rows])[good]
        if np.all(v <= 0) or np.std(v) == 0:
            continue
        lv = np.log10(np.maximum(v, 1e-300))
        if np.std(lv) == 0:
            continue
        c = float(np.corrcoef(lv, lw)[0, 1])
        cd = (float(np.corrcoef(lv[deep], lw[deep])[0, 1])
              if deep.sum() > 10 and np.std(lv[deep]) > 0 else float("nan"))
        # how far off is work/predictor across grids?  1.0 is a perfect
        # proportional predictor; this is what the cap actually relies on
        # Only grids where the predictor is actually defined: a grid with
        # zero stars says nothing about how well stars predict work, and
        # dividing by it manufactures a meaningless spread.
        nz = v > 0
        if nz.sum() < 20:
            continue
        ratio = w[good][nz] / v[nz]
        spread = float(np.percentile(ratio, 90) / np.percentile(ratio, 10))
        scored.append((abs(c), k, c, cd, spread))

    for _, k, c, cd, spread in sorted(scored, reverse=True):
        print("  %-14s %10.3f %17.3f %14.1fx" % (k, c, cd, spread))

    print()
    print("  'ratio spread' is the 90th/10th percentile of work divided by")
    print("  the predictor: how badly the cap would misjudge a grid.  The")
    print("  cell count is the incumbent - anything that does not beat it")
    print("  clearly is not worth the complexity.")


if __name__ == "__main__":
    main()
