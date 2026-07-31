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


def grid_geometry(hpath):
    """Per-grid geometry from the .hierarchy file, keyed by HDF5 group name.

    The HDF5 grid groups carry NO attributes in this Enzo build, so the
    cell volume cannot be read from there - it must come from the
    hierarchy.  Returns {group_name: dict(level, dx, vol, left, right,
    start, end, dims)}.
    """
    text = open(hpath).read()
    ptr_next = re.findall(
        r"Pointer:\s*Grid\[(\d+)\]->NextGridNextLevel\s*=\s*(\d+)", text)
    ptr_this = re.findall(
        r"Pointer:\s*Grid\[(\d+)\]->NextGridThisLevel\s*=\s*(\d+)", text)
    level_of = {1: 0}
    links = [(int(a), int(b), 1) for a, b in ptr_next if int(b)] + \
            [(int(a), int(b), 0) for a, b in ptr_this if int(b)]
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

        def vec(name, cast):
            mm = re.search(rf"{name}\s*=\s*([-\d\. eE+]+)", block)
            return [cast(x) for x in mm.group(1).split()] if mm else None

        le = vec("GridLeftEdge", float)
        re_ = vec("GridRightEdge", float)
        st = vec("GridStartIndex", int)
        en = vec("GridEndIndex", int)
        dm = vec("GridDimension", int)
        if not (le and re_ and st and en and dm):
            continue
        nactive = en[0] - st[0] + 1
        dx = (re_[0] - le[0]) / max(nactive, 1)
        out["Grid%08d" % gid] = {
            "level": level_of.get(gid, 0), "dx": dx, "vol": dx ** 3,
            "left": le, "right": re_, "start": st, "end": en, "dims": dm,
        }
    return out


def child_mask(gg, finer):
    """1.0 where this grid's active cells are NOT covered by a finer grid.

    Returned array is indexed [k][j][i] to match the HDF5 field layout.
    """
    import numpy as np
    nz, ny, nx = gg["shape"]
    keep = np.ones((nz, ny, nx), dtype=float)
    if not finer:
        return keep
    le, dx = gg["left"], gg["dx"]
    dims = (nx, ny, nz)
    for _, ch in finer:
        # index range of the child's footprint within this grid
        lo = [int(round((ch["left"][d] - le[d]) / dx)) for d in range(3)]
        hi = [int(round((ch["right"][d] - le[d]) / dx)) for d in range(3)]
        lo = [max(0, v) for v in lo]
        hi = [min(dims[d], hi[d]) for d in range(3)]
        if any(hi[d] <= lo[d] for d in range(3)):
            continue
        keep[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]] = 0.0
    return keep


def mass_sums(param_path):
    """Global mass sums from the dump's .cpu files. Needs h5py.

    Cell volumes come from the .hierarchy geometry (the HDF5 groups have
    no attributes).  Baryon sums are child-masked: cells covered by a
    finer grid are excluded, so the total is a true volume integral
    rather than a sum that double-counts refined regions.  Particle
    masses are stored by Enzo as densities relative to the HOST grid's
    cell volume, so they are multiplied by that grid's dx^3 - without
    this a particle that merely moves to a different level appears to
    change mass by the refinement factor cubed.

    Both corrections matter only when the two runs have DIFFERENT
    refinement structure (a C2-class change).  For same-structure
    comparisons the errors cancel and the older uncorrected sums were
    still valid as change detectors.
    """
    if not HAVE_H5PY:
        return None
    text = open(param_path).read()

    def fparam(name, default=None):
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([\dEe+.\-]+)", text, re.M)
        return float(m.group(1)) if m else default

    geom = grid_geometry(param_path + ".hierarchy")

    # Group grids by level so each grid can be masked against finer grids.
    by_level = {}
    for gname, gg in geom.items():
        by_level.setdefault(gg["level"], []).append((gname, gg))

    sums = {"baryon_mass": 0.0, "metal_mass": 0.0,
            "stellar_mass": 0.0, "dm_mass": 0.0, "star_count": 0,
            "metal_fields": {}}
    for cpu in sorted(glob.glob(param_path + ".cpu*")):
        with h5py.File(cpu, "r") as f:
            for gname, g in f.items():
                if not gname.startswith("Grid"):
                    continue
                gg = geom.get(gname)
                if gg is None:
                    continue
                vol = gg["vol"]
                dens = g.get("Density")
                if dens is not None:
                    d = dens[()]
                    # The HDF5 field arrays hold ACTIVE zones only (no
                    # ghosts), so the whole array is the active region.
                    gg["shape"] = d.shape
                    keep = child_mask(gg, by_level.get(gg["level"]+1, []))
                    sums["baryon_mass"] += float((d*keep).sum()) * vol
                    # Sum every metal-like field individually as well as in
                    # total: code versions name/split these differently
                    # (e.g. split-star-particles vs enzo-performance), and a
                    # single opaque total misreads cross-version comparisons
                    # as enormous metal differences.
                    for mf in g.keys():
                        if ((("etal" in mf or "olou" in mf or "olor" in mf)
                             and "Density" in mf)
                                or mf == "SN_Colour"):
                            v = float((g[mf][()]*keep).sum()) * vol
                            sums["metal_mass"] += v
                            sums["metal_fields"][mf] = \
                                sums["metal_fields"].get(mf, 0.0) + v
                pm = g.get("particle_mass")
                pt = g.get("particle_type")
                if pm is not None and pt is not None:
                    m = pm[()]
                    t = pt[()]
                    # particle_mass is a density w.r.t. THIS grid's cell
                    # volume; without the dx^3 a particle that changes
                    # level appears to change mass by RefineBy^3.
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


def mesh_and_migrations(la, lb):
    """Split a per-level structure diff into (mesh_differs, migrations).

    mesh_differs: grid or cell counts differ on any level - a genuine
    refinement-structure change, gated exactly.  migrations: total
    |particle-count delta| across levels with the mesh identical - single
    particles crossing level boundaries under parallel nondeterminism
    (T1.14), observed between identical-code runs at the production
    anchor, so gated like the star count rather than as structure.
    """
    keys = set(la) | set(lb)
    mesh_differs = any(
        la.get(l, {}).get("grids") != lb.get(l, {}).get("grids") or
        la.get(l, {}).get("cells") != lb.get(l, {}).get("cells")
        for l in keys)
    migrations = sum(
        abs(la.get(l, {}).get("particles", 0) -
            lb.get(l, {}).get("particles", 0)) for l in keys)
    return mesh_differs, migrations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("candidate")
    ap.add_argument("--class", dest="cls", choices=("C0", "C1", "C2"),
                    default="C1")
    ap.add_argument("--noise-floor", metavar="COMPARISON_JSON",
                    help="comparison.json from an identical-code pair (A1 vs "
                         "A2). Gates become max(class tolerance, factor x "
                         "floor) per metric per dump; structure/star-count "
                         "requirements relax from the first dump at which the "
                         "floor pair itself diverged. At production scale "
                         "parallel nondeterminism (ledger T1.14) makes this "
                         "mode mandatory - absolute tolerances only suit "
                         "serial/deterministic runs.")
    ap.add_argument("--noise-factor", type=float, default=10.0,
                    help="multiplier on the measured floor (default 10)")
    args = ap.parse_args()
    tol = TOL[args.cls]

    floor_rel = {}          # dump -> metric -> measured rel diff
    floor_star = {}         # dump -> star count delta
    floor_migr = {}         # dump -> per-level particle migration sum
    floor_struct_break = None  # first dump where the floor MESH diverged
    if args.noise_floor:
        fl = json.load(open(args.noise_floor))
        for dname, entry in fl.get("per_dump", {}).items():
            floor_rel[dname] = {k[len("rel_"):]: v for k, v in entry.items()
                                if k.startswith("rel_")}
            sa, sb = entry.get("sums_baseline"), entry.get("sums_candidate")
            if sa and sb:
                floor_star[dname] = abs(sa["star_count"] - sb["star_count"])
            la, lb = entry.get("levels_baseline"), entry.get("levels_candidate")
            if la is not None and lb is not None:
                m, migr = mesh_and_migrations(la, lb)
                floor_migr[dname] = migr
                if m and (floor_struct_break is None
                          or dname < floor_struct_break):
                    floor_struct_break = dname
            elif "particle_level_migrations" in entry:
                # merged envelope floors (merge_noise_floors.py) carry the
                # classification instead of the raw level tables
                floor_migr[dname] = entry["particle_level_migrations"]
                if entry.get("mesh_differs") and (
                        floor_struct_break is None
                        or dname < floor_struct_break):
                    floor_struct_break = dname
        # Legacy floors (before the mesh/particle split): fall back to the
        # failure messages, which conflated the two.
        if floor_struct_break is None and not floor_migr:
            for fmsg in fl.get("failures", []):
                if "per-level structure differs" in fmsg:
                    d = fmsg.split(":")[0].strip()
                    if floor_struct_break is None or d < floor_struct_break:
                        floor_struct_break = d

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
        entry["structure_differs"] = (ha != hb)
        mesh_differs, migrations = mesh_and_migrations(ha, hb)
        entry["mesh_differs"] = mesh_differs
        entry["particle_level_migrations"] = migrations
        struct_required = tol["structure_exact"] and not (
            floor_struct_break is not None and name >= floor_struct_break)
        if mesh_differs:
            if struct_required:
                report["failures"].append(f"{name}: mesh structure differs "
                                          "(per-level grid/cell counts)")
            else:
                print(f"note: {name}: mesh differs (within floor-pair "
                      f"divergence regime, not gated)")
        migr_gate = 0 if not args.noise_floor else \
            int(args.noise_factor * max(1, floor_migr.get(name, 0)))
        entry["gate_particle_level_migrations"] = migr_gate
        if migrations > migr_gate and tol["structure_exact"]:
            report["failures"].append(
                f"{name}: {migrations} particle level migration(s) "
                f"(> gate {migr_gate})")
        elif migrations:
            print(f"note: {name}: {migrations} particle(s) changed level "
                  f"with identical mesh (gate {migr_gate})")
        sa = mass_sums(da[name])
        sb = mass_sums(db[name])
        if sa and sb:
            entry["sums_baseline"] = sa
            entry["sums_candidate"] = sb
            for key in ("baryon_mass", "metal_mass", "stellar_mass", "dm_mass"):
                r = rel(sa[key], sb[key])
                entry[f"rel_{key}"] = r
                gate = tol["conserve"]
                if gate is not None and args.noise_floor:
                    gate = max(gate,
                               args.noise_factor *
                               floor_rel.get(name, {}).get(key, 0.0))
                entry[f"gate_{key}"] = gate
                if gate is not None and r > gate:
                    report["failures"].append(
                        f"{name}: {key} differs by {r:.3e} (> gate {gate:.3e})")
            dstar = abs(sa["star_count"] - sb["star_count"])
            entry["star_count_diff"] = dstar
            star_gate = 0 if not args.noise_floor else \
                int(args.noise_factor * max(1, floor_star.get(name, 0)))
            if dstar > star_gate:
                report["failures"].append(
                    f"{name}: star count {sa['star_count']} vs {sb['star_count']}"
                    f" (|diff| {dstar} > gate {star_gate})")
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
