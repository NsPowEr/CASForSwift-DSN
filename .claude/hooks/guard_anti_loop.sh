#!/usr/bin/env bash
# PreToolUse Bash guard: mechanize the CLAUDE.md "PROTOCOLLO ANTI-LOOP".
#
# The constitution says: after 3 consecutive failed attempts on the same error
# the agent must STOP, restore a safe state (git stash), and write a Stallo
# report. That rule was prose-only and unenforced. This hook gives it teeth by
# detecting blind repetition of *heavy* commands (build/test) within a time
# window.
#
# Signals (only "heavy" cmds counted — ninja/cmake/ctest/test binary/test_quick):
#   - same normalized command >= WARN times in window  -> non-blocking reminder
#   - same normalized command >= HARD times in window  -> DENY (force a stop)
#   - <= 2 distinct heavy cmds across >= OSC_TOTAL runs -> oscillation reminder
#
# Why "normalized": trivial wrappers (timeout/gtimeout/ASAN env) are stripped so
# the same logical run matches, but the --gtest_filter is KEPT: changing the
# filter means the agent is narrowing the problem (productive), not looping.
#
# Time window (WINDOW_SECS) makes the state self-healing: it decays on its own,
# so the hook can never wedge the agent permanently. A genuinely different
# strategy (new command / new filter) also clears the count immediately.
#
# Reads tool_input JSON from stdin (Claude Code hook protocol). Never throws.

set -euo pipefail

WINDOW_SECS="${ANTILOOP_WINDOW_SECS:-1200}"   # 20 min rolling window
WARN="${ANTILOOP_WARN:-3}"                    # nudge at 3rd identical run
HARD="${ANTILOOP_HARD:-6}"                    # block at 6th identical run
OSC_DISTINCT="${ANTILOOP_OSC_DISTINCT:-2}"    # oscillation: <= N distinct cmds
OSC_TOTAL="${ANTILOOP_OSC_TOTAL:-6}"          # ... across >= M total heavy runs

CMD=$(python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('tool_input',{}).get('command',''))" 2>/dev/null || echo "")
[[ -z "$CMD" ]] && exit 0

# Only track heavy build/test commands; everything else is noise for loop detection.
if ! echo "$CMD" | grep -qE '(ninja|cmake|ctest|cas_foundation_tests|test_quick\.sh|cas_tests|benchmark_tests)'; then
    exit 0
fi

# Normalize: drop leading env assignments and timeout wrappers, collapse spaces.
# Keep flags/filters intact so distinct attempts stay distinct.
KEY=$(echo "$CMD" \
    | sed -E 's/^([[:space:]]*[A-Za-z_][A-Za-z0-9_]*=[^[:space:]]+[[:space:]]+)+//' \
    | sed -E 's#(/usr/bin/)?g?timeout[[:space:]]+[0-9]+[[:space:]]+##g' \
    | tr -s '[:space:]' ' ' \
    | sed -E 's/^ //; s/ $//')
[[ -z "$KEY" ]] && exit 0

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STATE_DIR="$REPO_ROOT/.cache/anti_loop"
mkdir -p "$STATE_DIR"
HIST="$STATE_DIR/history.tsv"
touch "$HIST"

NOW=$(date +%s)
CUTOFF=$(( NOW - WINDOW_SECS ))

# Count occurrences of KEY and distinct heavy cmds within the window, then prune
# stale lines and append the current attempt. KEY passed via env (no awk -v
# backslash mangling). Output: "<count_of_key> <distinct> <total>".
read -r COUNT DISTINCT TOTAL < <(
    KEY="$KEY" CUTOFF="$CUTOFF" awk -F '\t' '
        BEGIN { key=ENVIRON["KEY"]; cutoff=ENVIRON["CUTOFF"]+0 }
        ($1+0) >= cutoff {
            total++
            seen[$2]=1
            if ($2 == key) c++
        }
        END {
            d=0; for (k in seen) d++
            # +1 to include the current (not-yet-written) attempt
            printf "%d %d %d\n", c+1, (d + ((key in seen)?0:1)), total+1
        }
    ' "$HIST"
)

# Rewrite history: keep in-window lines, append current attempt.
TMP="$STATE_DIR/.history.tmp.$$"
KEY="$KEY" CUTOFF="$CUTOFF" awk -F '\t' '($1+0) >= (ENVIRON["CUTOFF"]+0)' "$HIST" > "$TMP" 2>/dev/null || true
printf '%s\t%s\n' "$NOW" "$KEY" >> "$TMP"
mv -f "$TMP" "$HIST" 2>/dev/null || true

emit_deny() {
    KEY="$KEY" COUNT="$COUNT" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason":
      "PROTOCOLLO ANTI-LOOP (CLAUDE.md): comando build/test ripetuto identico "
      + os.environ["COUNT"] + " volte in <20min senza progresso:\n  "
      + os.environ["KEY"][:160] + "\n"
      "STOP ai fix ciechi. Esegui il protocollo: (1) fermati, (2) `git stash push` "
      "per stato sicuro (MAI reset --hard), (3) scrivi un Report di Stallo "
      "(natura errore + 3 strategie fallite + causa sistemica sospetta), "
      "(4) attendi intervento umano. Per ritentare DAVVERO cambia strategia "
      "(comando/filtro diversi) o attendi il decadimento della finestra."
  }
}))'
}

emit_warn() {
    local msg="$1"
    MSG="$msg" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "additionalContext": os.environ["MSG"]
  }
}))'
}

if (( COUNT >= HARD )); then
    emit_deny
    exit 0
fi

if (( COUNT >= WARN )); then
    emit_warn "anti-loop: hai lanciato lo STESSO comando build/test ${COUNT}× in <20min. Se l'errore non cambia, NON ripetere alla cieca: cambia ipotesi, restringi il filtro, o instrumenta con std::cerr. Al ${HARD}° identico il comando verrà bloccato (protocollo Stallo)."
    exit 0
fi

if (( DISTINCT <= OSC_DISTINCT && TOTAL >= OSC_TOTAL )); then
    emit_warn "anti-loop: oscillazione rilevata — ${TOTAL} run build/test ma solo ${DISTINCT} comandi distinti (pattern 'fix A rompe B, fix B rompe A'). Fermati a ragionare sulla causa sistemica condivisa prima del prossimo tentativo."
    exit 0
fi

exit 0
