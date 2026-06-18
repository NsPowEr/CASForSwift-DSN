#!/usr/bin/env bash
# PostToolUse Edit|Write hook: keep graphify knowledge-graph fresh.
#
# Rationale (CLAUDE.md graphify section + scaling to 4x): at >300k LOC raw grep
# stops scaling and `graphify query` becomes the primary navigation tool. The
# graph must therefore track the working tree without manual `graphify update .`.
#
# Design:
#   - Only triggers on source/header edits (.cpp/.hpp/.h/.cc) under src/ or include/.
#   - Debounced: at most one update per DEBOUNCE_SECS window (a burst of edits in
#     one task does not launch N concurrent rebuilds).
#   - Detached + backgrounded: never blocks the edit; AST-only, no API cost.
#   - Single-flight lock: if an update is already running, skip (the next edit
#     after it finishes will pick up all changes).

set -euo pipefail

DEBOUNCE_SECS=45

FILE=$(python3 -c "import json,sys; d=json.load(sys.stdin); ti=d.get('tool_input',{}); print(ti.get('file_path') or ti.get('path') or '')" 2>/dev/null || echo "")
[[ -z "$FILE" ]] && exit 0

case "$FILE" in
    *.cpp|*.hpp|*.h|*.cc) ;;
    *) exit 0 ;;
esac

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REL="${FILE#"$REPO_ROOT"/}"
case "$REL" in
    src/*|include/*) ;;
    *) exit 0 ;;
esac

# Graph not initialised yet → nothing to keep fresh.
[[ -f "$REPO_ROOT/graphify-out/graph.json" ]] || exit 0
command -v graphify >/dev/null 2>&1 || exit 0

STATE_DIR="$REPO_ROOT/.cache/graphify_hook"
mkdir -p "$STATE_DIR"
LOCK="$STATE_DIR/update.lock"
STAMP="$STATE_DIR/last_update"

# Debounce: skip if a successful update ran within the window.
if [[ -f "$STAMP" ]]; then
    now=$(date +%s)
    last=$(cat "$STAMP" 2>/dev/null || echo 0)
    (( now - last < DEBOUNCE_SECS )) && exit 0
fi

# Single-flight: skip if an update is already in flight.
if [[ -f "$LOCK" ]] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
    exit 0
fi

# Detach so the edit returns immediately; AST-only rebuild is cheap but non-zero.
(
    echo $$ > "$LOCK"
    cd "$REPO_ROOT" && graphify update . >/dev/null 2>&1 || true
    date +%s > "$STAMP"
    rm -f "$LOCK"
) </dev/null >/dev/null 2>&1 &
disown 2>/dev/null || true

exit 0
