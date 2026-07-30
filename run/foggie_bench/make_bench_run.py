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
        --nsteps 5 --ranks 128 --out bench_baseline
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
    r"CycleLastDataDump|StopTime|CosmologyOutputRedshift\[\d+\]|"
    r"NumberOfOutputsBeforeExit|TimingCycleSkip)\s*=")


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
    ap.add_argument("--ranks", type=int, default=128,
                    help="MPI ranks (default matches FOGGIE production: "
                         "one node, 128 ranks)")
    ap.add_argument("--out", required=True, help="Run directory to create")
    ap.add_argument("--walltime", default="4:00:00")
    ap.add_argument("--queue", default="normal")
    ap.add_argument("--group", default="s3128")
    ap.add_argument("--mpiexec",
                    default="/u/jtumlins/installs/mpich-4.0.3/usr/local/bin/mpiexec",
                    help="MPI launcher matching the MPI the build links against "
                         "(Make.mach.pleiades-mpich links MPICH 4.0.3, so the "
                         "system MPT mpiexec will NOT work)")
    ap.add_argument("--model", default="mil_ait",
                    help="PBS node model (mil_ait = Aitken Milan, the FOGGIE "
                         "production model; rom_ait = Aitken Rome, bro = "
                         "Pleiades Broadwell)")
    ap.add_argument("--ncpus", type=int, default=128, help="Cores per node")
    args = ap.parse_args()

    restart = os.path.abspath(args.restart)
    if not os.path.isfile(restart):
        sys.exit(f"error: restart parameter file not found: {restart}")
    enzo = os.path.abspath(args.enzo)
    if not os.path.isfile(enzo):
        sys.exit(f"error: enzo.exe not found: {enzo}")

    # The devel queue enforces a 2-hour walltime limit; clamp so the
    # submission is not rejected.
    if args.queue == "devel":
        h = int(args.walltime.split(":")[0])
        if h >= 2 and args.walltime != "2:00:00":
            print(f"note: devel queue limit is 2:00:00 - clamping walltime "
                  f"(was {args.walltime})")
            args.walltime = "2:00:00"

    out = os.path.abspath(args.out)
    if os.path.exists(out):
        sys.exit(f"error: {out} already exists - refusing to overwrite a run")
    os.makedirs(out)

    text = open(restart).read()
    # Enzo dump parameter files record the cycle as InitialCycleNumber
    # (WriteParameterFile.C); accept CycleNumber as a fallback.
    cycle = read_param(text, "InitialCycleNumber")
    if cycle is None:
        cycle = read_param(text, "CycleNumber")
    if cycle is None:
        sys.exit("error: InitialCycleNumber not found in restart parameter file")
    stop_cycle = int(cycle) + args.nsteps

    kept = [ln for ln in text.splitlines() if not OVERRIDES_DROP.match(ln)]
    kept += [
        "",
        "# ---- FOGGIE-bench overrides (make_bench_run.py) ----",
        f"StopCycle             = {stop_cycle}",
        "CycleSkipDataDump     = 1",
        "dtDataDump            = 0.0",
        "StopTime              = 1e20",
        # Production decks use NumberOfOutputsBeforeExit for requeue
        # chaining; with dumps at every root step an inherited value
        # would exit the run at the first output.
        "NumberOfOutputsBeforeExit = 0",
        # Benches need per-root-cycle timing data even if production
        # decks adopt the T0.6 recommendation of TimingCycleSkip = 10.
        "TimingCycleSkip       = 1",
    ]

    # Build a shadow snapshot directory inside the run dir: Enzo snapshots
    # are directories and are restarted as "-r RD0016/RD0016", with the
    # hierarchy's data paths relative to the run root. The shadow directory
    # symlinks every snapshot file in place and substitutes the edited
    # parameter file under its original name, so Enzo sees a normal
    # snapshot while the production dump stays untouched.
    snap = os.path.basename(restart)
    snap_dir = os.path.join(out, snap)
    os.makedirs(snap_dir)
    src_dir = os.path.dirname(restart)
    for entry in sorted(os.listdir(src_dir)):
        if entry == snap:
            continue  # replaced by the edited copy below
        os.symlink(os.path.join(src_dir, entry), os.path.join(snap_dir, entry))
    param_copy = os.path.join(snap_dir, snap)
    open(param_copy, "w").write("\n".join(kept) + "\n")
    restart_rel = os.path.join(snap, snap)

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
              .replace("@MPIEXEC@", args.mpiexec)
              .replace("@RESTART_REL@", restart_rel))
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
