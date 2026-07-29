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

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback


def _self_hash():
    try:
        return hashlib.sha1(open(os.path.abspath(__file__), "rb").read()
                            ).hexdigest()
    except OSError:
        return None


RUNNER_HASH_AT_START = _self_hash()

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


def sh(cmd, cwd=None, check=True, capture=True, timeout=600):
    # python3.6-compatible subprocess invocation: NAS system
    # python may be what cron finds. Every command gets a timeout so a
    # wedged filesystem or PBS call fails the entry loudly instead of
    # freezing the tick under the lock forever.
    pipe = subprocess.PIPE if capture else None
    try:
        r = subprocess.run(cmd, shell=True, cwd=cwd, universal_newlines=True,
                           stdout=pipe, stderr=pipe, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"command timed out after {timeout}s: {cmd}")
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
        return
    if not create_branch_if_missing:
        raise RuntimeError(f"branch {branch} missing and creation not allowed")

    # The branch does not exist on the remote. Creation must be idempotent:
    # an earlier tick may have created the local branch and its initial
    # commit but died at the push (typically the PAT not being set up yet).
    local = sh(f"git rev-parse --verify --quiet refs/heads/{branch}",
               cwd=path, check=False)
    if local:
        sh(f"git checkout -q {branch}", cwd=path)
        sh("git reset --hard && git clean -fdq", cwd=path)
    else:
        sh(f"git checkout -q --orphan {branch}", cwd=path)
        sh("git rm -rf --quiet . || true", cwd=path)
        open(os.path.join(path, "README.md"), "w").write(
            "FOGGIE-bench results branch - written only by hpc_runner.py\n")
        sh("git add README.md && git commit -m 'Initialize bench-results'",
           cwd=path)
    try:
        sh(f"git push -u origin {branch}", cwd=path)
    except RuntimeError as e:
        raise RuntimeError(
            str(e) + "\nhint: pushing bench-results failed - the GitHub PAT "
            "credential is probably not set up for this clone yet "
            "(HPC_SETUP.md step 2). The local branch is intact; rerun the "
            "tick after fixing credentials.")


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
    # Overlay the CURRENT branch tip's machine file: toolchain paths
    # (compiler, MPICH, grackle) must be current even when the code under
    # test is historical, or old commits fail on stale install paths.
    machfile = "Make.mach.nasa-aitken-milan-mpich"
    shutil.copyfile(
        os.path.join(REPO, "src", "enzo", machfile),
        os.path.join(bdir, "src", "enzo", machfile))
    sh("./configure", cwd=bdir)
    ed = os.path.join(bdir, "src", "enzo")
    sh("make machine-nasa-aitken-milan-mpich && make grackle-yes && make opt-high",
       cwd=ed)
    sh("make -j8", cwd=ed, capture=True, timeout=7200)
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
        st.setdefault("jobs", {})
        for label, exe in ((la, st["exe_a"]), (lb, st["exe_b"])):
            if label in st["jobs"]:
                continue
            out = os.path.join(rundir, f"bench_{label}")
            if not os.path.isdir(out):
                log(f"{eid}: preparing {label}")
                sh(f"python3 {mk} --restart {entry['restart']} --enzo {exe}"
                   f" --nsteps {entry.get('nsteps', 5)}"
                   f" --ranks {entry.get('ranks', 512)}"
                   f" --model {entry.get('model', 'rom_ait')}"
                   f" --ncpus {entry.get('ncpus', 128)}"
                   f" --queue {entry.get('queue', 'devel')}"
                   f" --out {out}", timeout=900)
            log(f"{eid}: submitting {label}")
            try:
                jid = sh(f"qsub {os.path.join(out, 'launch.sh')}", timeout=300)
            except RuntimeError as e:
                # Per-user queue limits (devel allows very few concurrent
                # jobs) are transient: keep the jobs already submitted, stay
                # in 'built', and retry the remainder on later ticks.
                if "would exceed" in str(e) or "per-user limit" in str(e):
                    log(f"{eid}: queue limit rejected {label}; retrying on a "
                        f"later tick")
                    break
                raise
            st["jobs"][label] = {"dir": out, "jobid": jid}
            log(f"{eid}: submitted {label} as {jid}")
        if len(st["jobs"]) == 2:
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
        # Optional measured noise floor: entry field 'noise_floor_json' is a
        # path relative to the bench-results clone (e.g.
        # results/t19-manual/noise_floor_A1_vs_A2.json). Mandatory in
        # practice for production-scale gating (ledger T1.14).
        floor_arg = ""
        nf = entry.get("noise_floor_json")
        if nf:
            nf_path = os.path.join(RESULTS, nf)
            if not os.path.isfile(nf_path):
                raise RuntimeError(
                    f"noise_floor_json not found in results clone: {nf}")
            floor_arg = f" --noise-floor {nf_path}"
        r = subprocess.run(
            f"python3 {cmp_script} {da} {db} --class {entry.get('class', 'C1')}"
            + floor_arg,
            shell=True, universal_newlines=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

    # If the sync just updated this very file, the in-memory code is stale
    # (a tick once processed a new queue entry with the previous build
    # logic and burned it). Re-exec so this tick runs the updated runner;
    # the wrapper's flock survives exec.
    if RUNNER_HASH_AT_START and _self_hash() not in (None, RUNNER_HASH_AT_START):
        log("runner source changed during sync; re-executing updated runner")
        sys.stdout.flush()
        os.execv(sys.executable, [sys.executable, os.path.abspath(__file__)])

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
