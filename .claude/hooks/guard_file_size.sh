#!/usr/bin/env bash
# PostToolUse Edit|Write guard: enforce CLAUDE.md anti-monolith on touched file.
#
# Thresholds (with mini-tolerance to avoid churn at 501-550 lines):
#   <= 500  → OK
#   501-550 → WARN (advisory, non-blocking; reminds to plan a split)
#   >  550  → BLOCK (hard fail per CLAUDE.md §ANTI-MONOLITH)
#
# Whitelist (scripts/file_size_whitelist.txt) is honored — exempt files skipped.
# Only scans .cpp/.hpp under src/ or include/.

set -euo pipefail

FILE=$(python3 -c "import json,sys; d=json.load(sys.stdin); ti=d.get('tool_input',{}); print(ti.get('file_path') or ti.get('path') or '')" 2>/dev/null || echo "")

if [[ -z "$FILE" || ! -f "$FILE" ]]; then
    exit 0
fi

# Reject paths with newlines: defends downstream consumers and keeps the path a
# single shell token. (Values are passed to python via argv below, never via
# source interpolation, so quote/metachar injection is not possible.)
if [[ "$FILE" == *$'\n'* ]]; then
    exit 0
fi

case "$FILE" in
    *.cpp|*.hpp) ;;
    *) exit 0 ;;
esac

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REL="${FILE#$REPO_ROOT/}"

case "$REL" in
    src/*|include/*) ;;
    *) exit 0 ;;
esac

WHITELIST="$REPO_ROOT/scripts/file_size_whitelist.txt"
if [[ -f "$WHITELIST" ]] && grep -qE "^[[:space:]]*${REL}([[:space:]]|#|$)" "$WHITELIST"; then
    exit 0
fi

LINES=$(wc -l < "$FILE" | tr -d ' ')
SOFT=500
HARD=550

if (( LINES <= SOFT )); then
    exit 0
fi

if (( LINES <= HARD )); then
    # Advisory warn — surface to model via additionalContext, do not block.
    # Values passed as argv (never interpolated into python source).
    python3 - "$REL" "$LINES" "$SOFT" "$HARD" <<'PY'
import json, sys
rel, lines, soft, hard = sys.argv[1:5]
msg = (f"anti-monolith warning: {rel} = {lines} lines "
       f"(soft limit {soft}, hard {hard}). Pianifica split a breve; "
       f"nessuna azione richiesta ora.")
print(json.dumps({
    "hookSpecificOutput": {
        "hookEventName": "PostToolUse",
        "additionalContext": msg,
    }
}))
PY
    exit 0
fi

# Hard block. Values passed as argv (never interpolated into python source).
python3 - "$REL" "$LINES" "$HARD" <<'PY'
import json, sys
rel, lines, hard = sys.argv[1:4]
msg = (f"anti-monolith BLOCK: {rel} = {lines} lines > {hard} "
       f"(CLAUDE.md §ANTI-MONOLITH, max 500 + 10% tolleranza). "
       f"Split obbligatorio o whitelist in scripts/file_size_whitelist.txt con ticket.")
print(json.dumps({"decision": "block", "reason": msg}))
PY
exit 0
