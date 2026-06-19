#!/usr/bin/env bash
# commit_msg_hook.sh — Repo-tracked commit-msg gate.
# Install with: ln -sf ../../scripts/commit_msg_hook.sh .git/hooks/commit-msg
#
# Enforces CLAUDE.md REGOLA 0.2: a commit that introduces a disabled test
# (DISABLED_/GTEST_SKIP/GTEST_DISABLE) MUST justify it in the commit body with
# a `TEST-DISABLED-JUSTIFY: <reason>` tag.
#
# Why here and not in pre-commit: at pre-commit time the commit message is not
# yet reliably written (esp. for `git commit -m`), so reading COMMIT_EDITMSG
# there false-passed/false-blocked. The commit-msg hook receives the final
# message file as $1 — the only correct place for a message-aware gate.
#
# Override: `git commit --no-verify` (emergencies only, never for agent commits).

set -u

MSG_FILE="${1:-}"
[[ -z "$MSG_FILE" || ! -f "$MSG_FILE" ]] && exit 0

ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0
cd "$ROOT" || exit 0

NEW_DISABLED=$(git diff --cached -U0 -- 'test/**/*.cpp' 'test/**/*.hpp' 2>/dev/null \
              | grep -E '^\+' | grep -E 'DISABLED_|GTEST_SKIP\(|GTEST_DISABLE' || true)

[[ -z "$NEW_DISABLED" ]] && exit 0

if grep -q 'TEST-DISABLED-JUSTIFY:' "$MSG_FILE"; then
  exit 0
fi

divider="────────────────────────────────────────────────────────────"
echo "$divider"
echo "✗ COMMIT-MSG BLOCK: new DISABLED_/GTEST_SKIP introduced without justification"
echo "$NEW_DISABLED" | sed 's/^/    /'
echo "  CLAUDE.md REGOLA 0.2: disabilitare un test richiede una giustificazione esplicita."
echo "  Aggiungi al corpo del commit:  TEST-DISABLED-JUSTIFY: <motivo>"
echo "  (oppure rimuovi il DISABLED_ e risolvi il bug — preferito.)"
echo "$divider"
exit 1
