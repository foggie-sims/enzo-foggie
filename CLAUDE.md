# Claude Code guidance for enzo-foggie (enzo-performance branch)

This branch carries the FOGGIE performance audit and its remediation.
Before doing audit work, read:

1. `audit/HANDOFF.md` - project state, standing policies, environment
   facts, and the current in-flight task. Start here.
2. `audit/status.yml` - the remediation ledger (authoritative status of
   every audit item). Update it in the same push as any fix.

Hard rules currently in force:

- **No PRs to main until further notice.** All work lands directly on
  the `enzo-performance` branch.
- Never commit `audit/dashboard.html` by hand - a GitHub Action rebuilds
  and commits it when `audit/status.yml` changes.
- The `bench-results` branch is a one-way results archive; never merge
  it with `enzo-performance`.
- Production A/B comparisons must use the measured noise floor
  (`compare_runs.py --noise-floor ...`); parallel Enzo is not bitwise
  reproducible (see T1.14 in the ledger and `run/foggie_bench/anchor.md`).

Build on Aitken: `make machine-nasa-aitken-milan-mpich && make
grackle-yes && make opt-high && make -j8` in `src/enzo` (toolchain
details and PBS gotchas in `audit/HANDOFF.md`).
