#!/bin/sh
# Cron wrapper for the FOGGIE-bench pull-model runner. Prepares the build
# and runtime environment (cron starts with a bare environment), prevents
# overlapping ticks with a lock, and invokes hpc_runner.py.
#
# Install (one specific front end - NAS crontabs are per-host):
#   crontab -e
#   */30 * * * * $HOME/foggie_bench_root/repo/run/foggie_bench/hpc_runner_cron.sh >> $HOME/foggie_bench_cron.log 2>&1
#
# Manual test first:  sh hpc_runner_cron.sh
#
# See HPC_SETUP.md for the one-time setup (clone, credentials).

export FOGGIE_BENCH_ROOT="${FOGGIE_BENCH_ROOT:-$HOME/foggie_bench_root}"

# --- toolchain (mirrors the PBS template and qsub_compile_enzo.sh) ---
if [ -f /usr/share/modules/init/sh ]; then
  . /usr/share/modules/init/sh
fi
module load comp-intel/2020.4.304 2>/dev/null
module load hdf5/1.8.18_serial 2>/dev/null

export LD_LIBRARY_PATH="/u/jtumlins/installs/mpich-4.0.3/usr/local/lib":"/u/jtumlins/installs/mpich-4.0.3/usr/lib":"/u/jtumlins/grackle/grackle-3.3.1-dev/build/lib64":$LD_LIBRARY_PATH
export PATH="/home1/jtumlins/anaconda3/bin:/nobackup/jtumlins/anaconda3/bin:/u/scicon/tools/bin/:/u/jtumlins/installs/mpich-4.0.3/usr/local/bin:$PATH"

# --- single-instance lock: builds and long ticks must not overlap ---
LOCK="$FOGGIE_BENCH_ROOT/.runner.lock"
mkdir -p "$FOGGIE_BENCH_ROOT"
if command -v flock >/dev/null 2>&1; then
  exec flock -n "$LOCK" python3 "$(dirname "$0")/hpc_runner.py"
else
  if mkdir "$LOCK.d" 2>/dev/null; then
    trap 'rmdir "$LOCK.d"' EXIT
    python3 "$(dirname "$0")/hpc_runner.py"
  else
    echo "runner already active; skipping tick"
  fi
fi
