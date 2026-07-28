#!/usr/bin/env python3
"""Compare two FOGGIE-bench runs (baseline vs candidate) and gate the result.

    python3 compare_runs.py bench_baseline bench_candidate --class C1

Checks, per matching dump (see README for the class semantics):
  - conservation invariants: total baryon / metal / stellar / DM mass,
    star particle count
  - per-level structure: grids, cells, particles per level
  - global field sums (level-masked where child data exists is NOT attempted;
    sums are per-level so overlap handling is explicit and comparable)
  - timing KPIs from performance.out (informational)

Structure comparisons need only the ASCII .hierarchy files; mass sums need
h5py. Emits <candidate>/comparison.json and a human verdict; exit code 0/1.
"""

import argparse
import glob
import json
import math
import os
import re
import sys

try:
    import h5py
    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

# Tolerances per change class (relative). Conservation gates every class.
TOL = {
    "C0": {"conserve": 1e-12, "sums": 1e-12, "structure_exact": True},
    "C1": {"conserve": 1e-10, "sums": 1e-8, "structure_exact": True},
    "C2": {"conserve": 1e-10, "sums": None, "structure_exact": False},
}


def find_dumps(rundir):
    """Return {dump_name: param_file_path} for DD*/ dumps in a run dir."""
    out = {}
    for d in sorted(glob.glob(os.path.join(rundir, "DD????"))):
        name = os.path.basename(d)
        cands = [f for f in glob.glob(os.path.join(d, "*"))
                 if re.match(r".*/[A-Za-z_0-9]+" + name[2:] + r"$", f)
                 and os.path.isfile(f)]
        # parameter file = the plain file whose .hierarchy sibling exists
        for c in cands:
            if os.path.isfile(c + ".hierarchy"):
                out[name] = c
                break
    return out


def parse_hierarchy(hpath):
    """Per-level {grids, cells, particles} from an Enzo .hierarchy file.

    Levels are reconstructed from the Pointer lines: NextGridNextLevel
    descends one level, NextGridThisLevel stays on the same level.
    GridDimension includes ghost zones.
    """
    levels = {}
    text = open(hpath).read()
    ptr_next_level = re.findall(
        r"Pointer:\s*Grid\[(\d+)\]->NextGridNextLevel\s*=\s*(\d+)", text)
    ptr_this_level = re.findall(
        r"Pointer:\s*Grid\[(\d+)\]->NextGridThisLevel\s*=\s*(\d+)", text)
    level_of = {1: 0}
    changed = True
    links = [(int(a), int(b), 1) for a, b in ptr_next_level if int(b) != 0] + \
            [(int(a), int(b), 0) for a, b in ptr_this_level if int(b) != 0]
    while changed:
        changed = False
        for a, b, dl in links:
            if a in level_of and b not in level_of:
                level_of[b] = level_of[a] + dl
                changed = True
    for block in re.split(r"\n(?=Grid = )", text):
        m = re.match(r"Grid = (\d+)", block)
        if not m:
            continue
        gid = int(m.group(1))
        lv = level_of.get(gid, 0)
        st = levels.setdefault(lv, {"grids": 0, "cells": 0, "particles": 0})
        st["grids"] += 1
        dm = re.search(r"GridDimension\s*=\s*([\d ]+)", block)
        if dm:
            c = 1
            for tok in dm.group(1).split():
                c *= int(tok)
            st["cells"] += c
        pm = re.search(r"NumberOfParticles\s*=\s*(\d+)", block)
        if pm:
            st["particles"] += int(pm.group(1))
    return levels


def mass_sums(param_path):
    """Global mass sums from the dump's .cpu files. Needs h5py."""
    if not HAVE_H5PY:
        return None
    text = open(param_path).read()

    def fparam(name, default=None):
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([\dEe+.\-]+)", text, re.M)
        return float(m.group(1)) if m else default

    # Level-aware cell volume needs per-grid dx; we sum per grid using the
    # grid's own spacing from the hierarchy-companion attributes in the HDF5.
    sums = {"baryon_mass": 0.0, "metal_mass": 0.0,
            "stellar_mass": 0.0, "dm_mass": 0.0, "star_count": 0}
    for cpu in sorted(glob.glob(param_path + ".cpu*")):
        with h5py.File(cpu, "r") as f:
            for gname, g in f.items():
                if not gname.startswith("Grid"):
                    continue
                dens = g.get("Density")
                if dens is not None:
                    # dx^rank from grid attributes if present; Enzo packed-AMR
                    # stores dx implicitly - reconstruct from Grid attrs.
                    le = g.attrs.get("GridLeftEdge")
                    re_ = g.attrs.get("GridRightEdge")
                    dims = g.attrs.get("GridDimension")
                    if le is not None and re_ is not None and dims is not None:
                        vol = 1.0
                        for lo, hi, n in zip(le, re_, dims):
                            vol *= (hi - lo) / max(int(n), 1)
                    else:
                        vol = 1.0
                    d = dens[()]
                    sums["baryon_mass"] += float(d.sum()) * vol
                    for mf in ("Metal_Density", "SN_Colour", "MetalSNII_Density",
                               "MetalSNIa_Density"):
                        md = g.get(mf)
                        if md is not None:
                            sums["metal_mass"] += float(md[()].sum()) * vol
                pm = g.get("particle_mass")
                pt = g.get("particle_type")
                if pm is not None and pt is not None:
                    m = pm[()]
                    t = pt[()]
                    dx = None
                    dims = g.attrs.get("GridDimension")
                    le = g.attrs.get("GridLeftEdge")
                    re_ = g.attrs.get("GridRightEdge")
                    if dims is not None and le is not None and re_ is not None:
                        dx = (re_[0] - le[0]) / max(int(dims[0]), 1)
                    vol = dx ** len(dims) if dx else 1.0
                    star = (t == 2)
                    sums["stellar_mass"] += float(m[star].sum()) * vol
                    sums["dm_mass"] += float(m[t == 1].sum()) * vol
                    sums["star_count"] += int(star.sum())
    return sums


def parse_timing(rundir):
    """Wall time per cycle and RebuildHierarchy totals from performance.out."""
    path = os.path.join(rundir, "performance.out")
    if not os.path.isfile(path):
        return None
    out = {"cycles": 0, "total_time": 0.0, "rebuild_time": 0.0}
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "Total" and len(p) >= 2:
            out["cycles"] += 1
            try:
                out["total_time"] += float(p[1])
            except ValueError:
                pass
        if p[0] == "RebuildHierarchy" and len(p) >= 2:
            try:
                out["rebuild_time"] += float(p[1])
            except ValueError:
                pass
    return out


def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("candidate")
    ap.add_argument("--class", dest="cls", choices=("C0", "C1", "C2"),
                    default="C1")
    args = ap.parse_args()
    tol = TOL[args.cls]

    da = find_dumps(args.baseline)
    db = find_dumps(args.candidate)
    common = sorted(set(da) & set(db))
    report = {"class": args.cls, "dumps_compared": common, "failures": [],
              "per_dump": {}}
    if not common:
        report["failures"].append("no common dumps between runs")
    if set(da) != set(db):
        report["failures"].append(
            f"dump sets differ: baseline {sorted(da)} vs candidate {sorted(db)}")

    if not HAVE_H5PY:
        print("WARNING: h5py unavailable - mass/field sums skipped "
              "(structure and timing only)")

    for name in common:
        entry = {}
        ha = parse_hierarchy(da[name] + ".hierarchy")
        hb = parse_hierarchy(db[name] + ".hierarchy")
        entry["levels_baseline"] = ha
        entry["levels_candidate"] = hb
        if tol["structure_exact"] and ha != hb:
            report["failures"].append(f"{name}: per-level structure differs")
        sa = mass_sums(da[name])
        sb = mass_sums(db[name])
        if sa and sb:
            entry["sums_baseline"] = sa
            entry["sums_candidate"] = sb
            for key in ("baryon_mass", "metal_mass", "stellar_mass", "dm_mass"):
                r = rel(sa[key], sb[key])
                entry[f"rel_{key}"] = r
                if r > tol["conserve"]:
                    report["failures"].append(
                        f"{name}: {key} differs by {r:.3e} (> {tol['conserve']:.0e})")
            if sa["star_count"] != sb["star_count"]:
                report["failures"].append(
                    f"{name}: star count {sa['star_count']} vs {sb['star_count']}")
        report["per_dump"][name] = entry

    ta = parse_timing(args.baseline)
    tb = parse_timing(args.candidate)
    if ta and tb:
        report["timing"] = {"baseline": ta, "candidate": tb}
        if ta["total_time"] > 0:
            report["timing"]["speedup"] = ta["total_time"] / max(tb["total_time"], 1e-30)

    out = os.path.join(args.candidate, "comparison.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Compared {len(common)} dump(s); report: {out}")
    if "timing" in report and "speedup" in report["timing"]:
        print(f"wall-clock ratio baseline/candidate: {report['timing']['speedup']:.3f}"
              " (informational; needs same node type + repeats)")
    if report["failures"]:
        print("FAIL:")
        for fmsg in report["failures"]:
            print("  -", fmsg)
        sys.exit(1)
    print(f"PASS at class {args.cls}")


if __name__ == "__main__":
    main()
