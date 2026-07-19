#!/usr/bin/env bash
# Stop hook (advisory, never blocking): end-of-turn state reminder.
#
# Fires when the agent ends its turn. If the working tree carries uncommitted
# changes under src/ or include/, emit a one-line reminder so "work declared
# done" and "work actually gated+committed" cannot silently diverge across
# sessions (memory: uncommitted-sweep-audit-2026-07-07).
#
# Advisory only: prints context, never blocks the stop (a blocking Stop hook
# can wedge a session in a loop — explicitly avoided).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT" || exit 0

command -v git >/dev/null 2>&1 || exit 0
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || exit 0

DIRTY=$(git status --porcelain -- src include test CMakeLists.txt 2>/dev/null | head -20 || true)
[[ -z "$DIRTY" ]] && exit 0

# Debounce: emit only when the dirty set CHANGES, not on every stop — a stable
# checkpoint (e.g. another session's work in progress) must not spam each turn.
STATE_DIR="$REPO_ROOT/.cache/stop_state_hook"
mkdir -p "$STATE_DIR" 2>/dev/null || exit 0
CUR_HASH=$(echo "$DIRTY" | shasum -a 256 | awk '{print $1}')
LAST_HASH=$(cat "$STATE_DIR/last_dirty_hash" 2>/dev/null || echo "")
[[ "$CUR_HASH" == "$LAST_HASH" ]] && exit 0
echo "$CUR_HASH" > "$STATE_DIR/last_dirty_hash"

N=$(echo "$DIRTY" | wc -l | tr -d ' ')
FILES=$(echo "$DIRTY" | awk '{print $2}' | head -5 | tr '\n' ' ')

DIRTY_N="$N" DIRTY_FILES="$FILES" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "Stop",
    "additionalContext":
      "stato tree a fine turno: " + os.environ["DIRTY_N"] +
      " file sorgente modificati non committati (" + os.environ["DIRTY_FILES"].strip() +
      " ...). Se il lavoro è dichiarato finito: gate + commit atomico + stato in TASKLIST_MASTER.md. "
      "Se è un checkpoint intenzionale o lavoro di ALTRA sessione: nessuna azione (non committare roba altrui)."
  }
}))'
exit 0
