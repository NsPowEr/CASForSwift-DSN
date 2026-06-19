#!/usr/bin/env bash
# pre_commit_hook.sh — Repo-tracked pre-commit gate.
# Install with: ln -sf ../../scripts/pre_commit_hook.sh .git/hooks/pre-commit
#
# This script is a thin orchestrator: the four core debt gates
# (anti-monolith, orphan sources, untracked/stale, HARDCODE-OF-PASSAGE)
# are owned by scripts/debt_gate.sh so that CI (non-bypassable) and
# local hooks (--no-verify bypassable) share the exact same logic.
#
# The disabled-test justification check (REGOLA 0.2) is enforced separately in
# the commit-msg hook (scripts/commit_msg_hook.sh), where the final commit
# message is reliably available — reading COMMIT_EDITMSG here was unreliable.
#
# Override: `git commit --no-verify` (use only for emergency, never for agent commits).

set -u

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT" || exit 2

fail=0
divider="────────────────────────────────────────────────────────────"

# ── Shared gates (also enforced by CI) ───────────────────────────────────────
if ! out=$(bash scripts/debt_gate.sh --staged 2>&1); then
  echo "$out"
  fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "$divider"
  echo "Pre-commit blocked. Resolve the above, or use 'git commit --no-verify' for emergencies."
  exit 1
fi

exit 0
