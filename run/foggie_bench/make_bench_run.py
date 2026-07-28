#!/usr/bin/env python3
"""Prepare an abbreviated FOGGIE production run from a restart dump.

Creates an output directory containing:
  - an edited copy of the restart parameter file (StopCycle limited,
    per-root-step cycle-based dumps, dt/redshift-based dumps disabled)
  - launch.sh generated from the PBS template
  - bench_manifest.json recording exactly what was prepared

The production restart data is referenced in place and never modified.

Example:
    python3 make_bench_run.py \
        --restart /nobackup/.../DD1234/DD1234 \
        --enzo /home/.../enzo-foggie/src/enzo/enzo.exe \
        --nsteps 5 --ranks 512 --out bench_baseline
"""

import argparse
import json
import os
import re
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "pbs_pleiades.template.sh")

# Parameters we override for the abbreviated run: cycle-driven output at
# every root step, everything time/redshift-driven off, so baseline and
# candidate dump at identical cycles regardless of dt differences.
OVERRIDES_DROP = re.compile(
    r"^\s*(StopCycle|CycleSkipDataDump|dtDataDump|TimeLastDataDump|"
    r"CycleLastDataDump|StopTime|CosmologyOutputRedshift\[\d+\])\s*=")


def read_param(text, name):
    m = re.search(rf"^\s*{re.escape(name)}\s*=\s*(\S+)", text, re.M)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restart", required=True,
                    help="Path to the restart parameter file, e.g. .../DD1234/DD1234")
    ap.add_argument("--enzo", required=True, help="Path to enzo.exe to test")
    ap.add_argument("--nsteps", type=int, default=5,
                    help="Number of root-grid steps to run (default 5)")
    ap.add_argument("--ranks", type=int, default=512, help="MPI ranks")
    ap.add_argument("--out", required=True, help="Run directory to create")
    ap.add_argument("--walltime", default="4:00:00")
    ap.add_argument("--queue", default="normal")
    ap.add_argument("--group", default="s2961")
    ap.add_argument("--model", default="rom_ait",
                    help="PBS node model (e.g. rom_ait for Aitken Rome, bro for Pleiades Broadwell)")
    ap.add_argument("--ncpus", type=int, default=128, help="Cores per node")
    args = ap.parse_args()

    restart = os.path.abspath(args.restart)
    if not os.path.isfile(restart):
        sys.exit(f"error: restart parameter file not found: {restart}")
    enzo = os.path.abspath(args.enzo)
    if not os.path.isfile(enzo):
        sys.exit(f"error: enzo.exe not found: {enzo}")

    out = os.path.abspath(args.out)
    if os.path.exists(out):
        sys.exit(f"error: {out} already exists - refusing to overwrite a run")
    os.makedirs(out)

    text = open(restart).read()
    cycle = read_param(text, "CycleNumber")
    if cycle is None:
        sys.exit("error: CycleNumber not found in restart parameter file")
    stop_cycle = int(cycle) + args.nsteps

    kept = [ln for ln in text.splitlines() if not OVERRIDES_DROP.match(ln)]
    kept += [
        "",
        "# ---- FOGGIE-bench overrides (make_bench_run.py) ----",
        f"StopCycle             = {stop_cycle}",
        "CycleSkipDataDump     = 1",
        "dtDataDump            = 0.0",
        "StopTime              = 1e20",
    ]
    param_copy = os.path.join(out, os.path.basename(restart) + ".bench")
    open(param_copy, "w").write("\n".join(kept) + "\n")

    nodes = max(1, (args.ranks + args.ncpus - 1) // args.ncpus)
    template = open(TEMPLATE).read()
    launch = (template
              .replace("@OUT@", out)
              .replace("@RANKS@", str(args.ranks))
              .replace("@NODES@", str(nodes))
              .replace("@NCPUS@", str(args.ncpus))
              .replace("@MODEL@", args.model)
              .replace("@WALLTIME@", args.walltime)
              .replace("@QUEUE@", args.queue)
              .replace("@GROUP@", args.group)
              .replace("@ENZO@", enzo)
              .replace("@PARAM@", param_copy)
              .replace("@RESTART_DIR@", os.path.dirname(restart)))
    launch_path = os.path.join(out, "launch.sh")
    open(launch_path, "w").write(launch)
    os.chmod(launch_path, os.stat(launch_path).st_mode | stat.S_IXUSR)

    manifest = {
        "restart": restart,
        "restart_cycle": int(cycle),
        "stop_cycle": stop_cycle,
        "nsteps": args.nsteps,
        "enzo": enzo,
        "ranks": args.ranks,
        "nodes": nodes,
        "param_file": param_copy,
    }
    with open(os.path.join(out, "bench_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Prepared {out}")
    print(f"  restart cycle {cycle} -> StopCycle {stop_cycle} ({args.nsteps} root steps)")
    print(f"  {args.ranks} ranks on {nodes} x {args.model} nodes")
    print(f"  submit with: qsub {launch_path}")


if __name__ == "__main__":
    main()
