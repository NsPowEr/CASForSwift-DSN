#!/usr/bin/env bash
# SessionStart hook: inject a compact live-state snapshot so every fresh
# session starts oriented without re-deriving branch/ledger/golden/oracle
# state by hand (token cost paid once here, saved every session).
#
# Output on stdout becomes session context (Claude Code hook protocol).
# Hard budget: keep it under ~30 lines; details live in TASKLIST_MASTER.md.
# Never fails the session: every probe is best-effort.

set -u
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT" || exit 0

echo "=== SNAPSHOT PROGETTO (session_state.sh) ==="

# Git: branch + dirty files + last commit.
BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")
DIRTY=$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')
LAST=$(git log --oneline -1 2>/dev/null || echo "?")
echo "git: branch=${BRANCH} · file modificati/untracked=${DIRTY} · last: ${LAST}"
if [ "${DIRTY}" != "0" ]; then
    echo "in-flight (git status --short, primi 8):"
    git status --short 2>/dev/null | head -8 | sed 's/^/  /'
fi

# Hardcode ledger: one summary line.
if [ -f scripts/ledger_index.py ]; then
    STATS=$(python3 scripts/ledger_index.py stats 2>/dev/null | head -2 | tr '\n' ' ')
    [ -n "$STATS" ] && echo "ledger: ${STATS}"
fi

# Golden ratchet floor (mathematical-truth guardian).
if [ -f scripts/golden_baseline.txt ]; then
    FLOOR=$(grep -E '^(STATUS|PASS_FLOOR|FAIL_CEILING):' scripts/golden_baseline.txt | tr '\n' ' ')
    echo "golden ratchet vs Maxima: ${FLOOR}"
fi

# Reference oracles present?
M_BIN=$(command -v maxima >/dev/null 2>&1 && echo ok || echo MANCANTE)
G_BIN=$(command -v icas >/dev/null 2>&1 && echo ok || echo MANCANTE)
echo "oracoli: maxima=${M_BIN} · giac(icas)=${G_BIN}"

# Open/in-progress tasks from the single source of truth (headers only).
if [ -f TASKLIST_MASTER.md ]; then
    OPEN=$(grep -E '^### ' TASKLIST_MASTER.md | grep -E 'APERT|IN CORSO|QUASI-FATTO' | sed -E 's/^### //; s/ — \[E.*$//' | head -8)
    if [ -n "$OPEN" ]; then
        echo "task aperti/in corso (TASKLIST_MASTER.md):"
        echo "$OPEN" | sed 's/^/  /'
    fi
fi

echo "SoT: TASKLIST_MASTER.md · regole: CLAUDE.md → docs/rules/ · storico: docs/archive/"
echo "============================================"
exit 0
