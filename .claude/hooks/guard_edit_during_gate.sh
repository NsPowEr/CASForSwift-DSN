#!/usr/bin/env bash
# PreToolUse Edit|Write guard: no source edits while a gate is running.
#
# Mechanizes the memory rule "no edit sorgenti durante gate bg"
# (ledger-classifier-traps): editing src/include/test mid-gate poisons the
# incremental build the suite is running against, and editing a script that
# bash is currently EXECUTING corrupts it (bash reads the file incrementally).
#
# Behaviour:
#   - No gate running                        -> allow silently.
#   - Target under src/ include/ test/ or a  -> DENY with explanation.
#     top-level CMakeLists.txt
#   - Target under scripts/ AND that script  -> DENY (running-script corruption).
#     appears in the running gate cmdline
#   - Everything else (docs/, .claude/, scratchpad, memoria, ...) -> allow.
#
# Reads tool_input JSON from stdin (Claude Code hook protocol). Never throws.

set -euo pipefail

FILE=$(python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('tool_input',{}).get('file_path',''))" 2>/dev/null || echo "")
[[ -z "$FILE" ]] && exit 0

GATE_PATTERN='test_quick\.sh|debt_gate\.sh|run_golden_measurement\.sh|check_golden_ratchet\.sh|benchmark\.sh|cas_foundation_tests'
RUNNING=$(pgrep -lf "$GATE_PATTERN" 2>/dev/null | grep -v "pgrep" || true)
[[ -z "$RUNNING" ]] && exit 0

ROOT="${CLAUDE_PROJECT_DIR:-$(pwd)}"
REL="${FILE#"$ROOT"/}"

BLOCK=0
case "$REL" in
    src/*|include/*|test/*|CMakeLists.txt) BLOCK=1 ;;
    scripts/*)
        if echo "$RUNNING" | grep -qF "$(basename "$REL")"; then BLOCK=1; fi ;;
esac
[[ "$BLOCK" != "1" ]] && exit 0

SUMMARY=$(echo "$RUNNING" | head -3 | awk '{$1=""; print substr($0,2,120)}' | tr '\n' ';')

RUNNING_SUMMARY="$SUMMARY" TARGET="$REL" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason":
      "EDIT-DURANTE-GATE (memoria ledger-classifier-traps): un gate/suite è in esecuzione:\n  "
      + os.environ["RUNNING_SUMMARY"] + "\n"
      "Modificare ora `" + os.environ["TARGET"] + "` avvelena la build incrementale "
      "su cui il gate sta girando (o corrompe uno script in esecuzione). Attendi la "
      "fine del gate (pgrep -lf test_quick ; pgrep -lf run_golden_measurement), poi riprova. "
      "Edit liberi nel frattempo: docs/, .claude/, scratchpad, memoria."
  }
}))'
exit 0
