# One-time setup: FOGGIE-bench pull-model automation on NAS

The runner is a cron job on one front end that pulls the `enzo-performance`
branch, executes `queue.yml` entries (build -> submit -> compare), and pushes
results to the `bench-results` branch. Nothing connects *into* NAS; the
front end only makes outbound HTTPS to GitHub. Queue authorship
(enzo-performance) and result reporting (bench-results) are strictly one-way
in opposite directions, so the runner and the branch's developers never
contend for the same ref.

## 1. Choose a home and clone

    export FOGGIE_BENCH_ROOT=$HOME/foggie_bench_root      # or /nobackup/...
    mkdir -p $FOGGIE_BENCH_ROOT
    git clone -b enzo-performance \
        https://github.com/foggie-sims/enzo-foggie.git $FOGGIE_BENCH_ROOT/repo

(The runner manages `repo/` with hard resets - do not use it as a working
copy for development.)

## 2. Push credential for bench-results

The runner needs push access for the `bench-results` branch only. Create a
**fine-grained GitHub PAT** restricted to `foggie-sims/enzo-foggie` with
Contents: Read and write, then store it on NAS only:

    git config --global credential.helper "store --file $HOME/.foggie_bench_git_creds"
    cd $FOGGIE_BENCH_ROOT/repo && git fetch   # enter username + PAT once when prompted for push later

(Or scope the helper to just these clones with per-repo `git config` if you
prefer.) The token lives only on NAS hardware, is yours, and is revocable
from GitHub at any time.

## 3. Test a tick manually

    sh $FOGGIE_BENCH_ROOT/repo/run/foggie_bench/hpc_runner_cron.sh

With an empty queue this should: sync `repo/`, create and push the
`bench-results` branch (first run only), print "tick: 0 queue entrie(s)",
and exit. Fix any credential prompts now - cron cannot answer them.

## 4. Install the crontab (host-specific!)

NAS crontabs live per-host: pick ONE front end (e.g. `ssh pfe21`) and
always edit the crontab there.

    crontab -e
    */30 * * * * $HOME/foggie_bench_root/repo/run/foggie_bench/hpc_runner_cron.sh >> $HOME/foggie_bench_cron.log 2>&1

The lock in the wrapper makes overlapping ticks impossible; a tick that
starts a 20-minute build simply finishes it while later ticks skip.

## 5. Using it

Enqueue work by committing an entry to `run/foggie_bench/queue.yml` on
`enzo-performance` (schema documented in the file). Watch progress in
`results/<id>/state.json` on the `bench-results` branch; final verdicts are
`results/<id>/comparison.json` plus the comparator stdout in the state file.
Phases: `pending -> built -> submitted -> ran -> done` (or `failed` /
`done-with-failures`, with error tails captured in state.json).

## Operational notes

- **Builds** are cached per commit sha under `builds/`; each unique ref is
  compiled once (~10 min).
- **Disk**: every entry's runs write full dumps under `runs/<id>/`; clean
  old entries' run dirs by hand when done (the results branch keeps only
  small JSON/log files).
- **Stopping**: `crontab -e` and comment the line out. In-flight PBS jobs
  finish on their own; the next enabled tick picks their status up.
- **Untested-until-first-tick**: like the rest of the harness, expect a
  shakedown - run step 3 manually and watch the first real entry end to end
  before trusting the cron cadence.
