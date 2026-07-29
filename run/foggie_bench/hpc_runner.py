#!/usr/bin/env python3
"""FOGGIE-bench pull-model runner. Invoked by hpc_runner_cron.sh on a NAS
front end; do not run outside the environment that wrapper prepares.

Each tick:
  1. hard-syncs a runner-owned clone of enzo-performance (the queue + code),
  2. advances every queue entry through a small state machine:
       pending -> built -> submitted -> ran -> done | failed
     (building each unique ref once, submitting runs with make_bench_run.py
      + qsub, detecting completion via the bench_exit_status file the PBS
      script writes, then running compare_runs.py),
  3. pushes per-entry state.json / comparison.json / log tails to the
     bench-results branch.

State lives ONLY on bench-results; the queue on enzo-performance is never
written. Everything is idempotent per tick; the cron wrapper holds a lock so
ticks never overlap. Long builds/runs simply span multiple ticks.

Layout under $FOGGIE_BENCH_ROOT:
  repo/     clone of <REPO_URL> at enzo-performance (runner-managed)
  results/  clone at bench-results (runner-managed, pushed)
  builds/<sha>/   one build tree per unique commit
  runs/<id>/      bench run directories
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

try:
    import yaml
except ImportError:
    sys.exit("hpc_runner: PyYAML required (anaconda python3 has it)")

ROOT = os.environ.get("FOGGIE_BENCH_ROOT",
                      os.path.expanduser("~/foggie_bench_root"))
REPO_URL = os.environ.get("FOGGIE_BENCH_REPO",
                          "https://github.com/foggie-sims/enzo-foggie.git")
BRANCH = "enzo-performance"
RESULTS_BRANCH = "bench-results"

REPO = os.path.join(ROOT, "repo")
RESULTS = os.path.join(ROOT, "results")
BUILDS = os.path.join(ROOT, "builds")
RUNS = os.path.join(ROOT, "runs")


def sh(cmd, cwd=None, check=True, capture=True):
    r = subprocess.run(cmd, shell=True, cwd=cwd, text=True,
                       capture_output=capture)
    if check and r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}\n"
                           f"stdout: {(r.stdout or '')[-2000:]}\n"
                           f"stderr: {(r.stderr or '')[-2000:]}")
    return (r.stdout or "").strip()


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def sync_clone(path, branch, create_branch_if_missing=False):
    if not os.path.isdir(os.path.join(path, ".git")):
        sh(f"git clone --branch {BRANCH} {REPO_URL} {path}")
    sh("git fetch origin", cwd=path)
    remote = sh(f"git ls-remote --heads origin {branch}", cwd=path, check=False)
    if remote:
        sh(f"git checkout -q {branch} 2>/dev/null || git checkout -q -b {branch} origin/{branch}",
           cwd=path)
        sh(f"git reset --hard origin/{branch}", cwd=path)
    elif create_branch_if_missing:
        sh(f"git checkout -q --orphan {branch}", cwd=path)
        sh("git rm -rf --quiet . || true", cwd=path)
        open(os.path.join(path, "README.md"), "w").write(
            "FOGGIE-bench results branch - written only by hpc_runner.py\n")
        sh("git add README.md && git commit -m 'Initialize bench-results'",
           cwd=path)
        sh(f"git push -u origin {branch}", cwd=path)
    else:
        raise RuntimeError(f"branch {branch} missing and creation not allowed")


def resolve_sha(ref):
    return sh(f"git rev-parse {ref}^{{commit}}", cwd=REPO)


def build(sha):
    """Build enzo.exe for a commit; returns path. Cached by sha."""
    bdir = os.path.join(BUILDS, sha[:12])
    exe = os.path.join(bdir, "src", "enzo", "enzo.exe")
    if os.path.isfile(exe):
        return exe
    if os.path.isdir(bdir):
        shutil.rmtree(bdir)
    log(f"building {sha[:12]}")
    sh(f"git clone --no-checkout {REPO} {bdir}")
    sh(f"git checkout -q {sha}", cwd=bdir)
    sh("./configure", cwd=bdir)
    ed = os.path.join(bdir, "src", "enzo")
    sh("make machine-pleiades-mpich && make grackle-yes && make opt-high",
       cwd=ed)
    sh("make -j8", cwd=ed, capture=True)
    if not os.path.isfile(exe):
        raise RuntimeError(f"build produced no enzo.exe for {sha[:12]}")
    return exe


def state_path(eid):
    return os.path.join(RESULTS, "results", eid, "state.json")


def load_state(eid):
    p = state_path(eid)
    if os.path.isfile(p):
        return json.load(open(p))
    return {"phase": "pending"}


def save_state(eid, st):
    p = state_path(eid)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    st["updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    json.dump(st, open(p, "w"), indent=2)


def copy_to_results(eid, src, name):
    dst = os.path.join(RESULTS, "results", eid, name)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)


def tail_of(path, n=100):
    if not os.path.isfile(path):
        return "(missing)"
    lines = open(path, errors="replace").readlines()
    return "".join(lines[-n:])


def run_labels(entry):
    same = entry["baseline_ref"] == entry["candidate_ref"]
    return ("A1", "A2") if same else ("A", "B")


def advance(entry):
    eid = entry["id"]
    st = load_state(eid)
    la, lb = run_labels(entry)

    if st["phase"] == "pending":
        st["baseline_sha"] = resolve_sha(entry["baseline_ref"])
        st["candidate_sha"] = resolve_sha(entry["candidate_ref"])
        st["exe_a"] = build(st["baseline_sha"])
        st["exe_b"] = build(st["candidate_sha"])
        st["phase"] = "built"
        save_state(eid, st)

    if st["phase"] == "built":
        rundir = os.path.join(RUNS, eid)
        os.makedirs(rundir, exist_ok=True)
        mk = os.path.join(REPO, "run", "foggie_bench", "make_bench_run.py")
        st["jobs"] = {}
        for label, exe in ((la, st["exe_a"]), (lb, st["exe_b"])):
            out = os.path.join(rundir, f"bench_{label}")
            if not os.path.isdir(out):
                sh(f"python3 {mk} --restart {entry['restart']} --enzo {exe}"
                   f" --nsteps {entry.get('nsteps', 5)}"
                   f" --ranks {entry.get('ranks', 512)}"
                   f" --model {entry.get('model', 'rom_ait')}"
                   f" --ncpus {entry.get('ncpus', 128)}"
                   f" --queue {entry.get('queue', 'devel')}"
                   f" --out {out}")
            jid = sh(f"qsub {os.path.join(out, 'launch.sh')}")
            st["jobs"][label] = {"dir": out, "jobid": jid}
            log(f"{eid}: submitted {label} as {jid}")
        st["phase"] = "submitted"
        save_state(eid, st)
        return  # give the jobs at least one tick

    if st["phase"] == "submitted":
        done, failed = [], []
        for label, j in st["jobs"].items():
            status_file = os.path.join(j["dir"], "bench_exit_status")
            if os.path.isfile(status_file):
                code = open(status_file).read().strip()
                (done if code == "exit=0" else failed).append(label)
                copy_to_results(eid, os.path.join(j["dir"], "enzo_bench.log"),
                                f"enzo_bench_{label}.tail.log") \
                    if os.path.isfile(os.path.join(j["dir"], "enzo_bench.log")) else None
            else:
                q = sh(f"qstat {j['jobid']} 2>/dev/null", check=False)
                if not q:
                    failed.append(label)
                    st.setdefault("errors", []).append(
                        f"{label}: job {j['jobid']} left the queue without "
                        f"writing bench_exit_status")
        if failed:
            st["phase"] = "failed"
            for label in failed:
                logf = os.path.join(st["jobs"][label]["dir"], "enzo_bench.log")
                st.setdefault("errors", []).append(
                    f"{label} log tail:\n{tail_of(logf, 40)}")
            save_state(eid, st)
        elif len(done) == len(st["jobs"]):
            st["phase"] = "ran"
            save_state(eid, st)
        return

    if st["phase"] == "ran":
        cmp_script = os.path.join(REPO, "run", "foggie_bench", "compare_runs.py")
        da = st["jobs"][la]["dir"]
        db = st["jobs"][lb]["dir"]
        r = subprocess.run(
            f"python3 {cmp_script} {da} {db} --class {entry.get('class', 'C1')}",
            shell=True, text=True, capture_output=True)
        st["compare_stdout"] = r.stdout[-4000:]
        st["compare_exit"] = r.returncode
        cj = os.path.join(db, "comparison.json")
        if os.path.isfile(cj):
            copy_to_results(eid, cj, "comparison.json")
        st["phase"] = "done" if r.returncode == 0 else "done-with-failures"
        save_state(eid, st)
        log(f"{eid}: comparison exit {r.returncode}")


def push_results():
    sh("git add -A", cwd=RESULTS)
    if sh("git status --porcelain", cwd=RESULTS):
        sh("git commit -m 'bench runner update'", cwd=RESULTS)
        sh(f"git push origin {RESULTS_BRANCH}", cwd=RESULTS)
        log("pushed results")


def main():
    for d in (ROOT, BUILDS, RUNS):
        os.makedirs(d, exist_ok=True)
    sync_clone(REPO, BRANCH)
    sync_clone(RESULTS, RESULTS_BRANCH, create_branch_if_missing=True)

    qpath = os.path.join(REPO, "run", "foggie_bench", "queue.yml")
    queue = yaml.safe_load(open(qpath)) or {}
    entries = queue.get("entries") or []
    log(f"tick: {len(entries)} queue entrie(s)")

    for entry in entries:
        eid = entry.get("id")
        if not eid or not re.match(r"^[A-Za-z0-9._-]+$", eid):
            log(f"skipping entry with bad id: {entry}")
            continue
        st = load_state(eid)
        if st["phase"] in ("done", "done-with-failures", "failed"):
            continue
        try:
            advance(entry)
        except Exception:
            st = load_state(eid)
            st["phase"] = "failed"
            st.setdefault("errors", []).append(traceback.format_exc()[-4000:])
            save_state(eid, st)
            log(f"{eid}: FAILED (see state.json)")

    push_results()


if __name__ == "__main__":
    main()
