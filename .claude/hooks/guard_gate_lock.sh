#!/usr/bin/env bash
# PreToolUse Bash guard: no heavy build/test while a gate is already running.
#
# Mechanizes the memory rule "no-concurrent-gates-load-timeout": launching a
# build or a test binary while test_quick.sh / debt_gate.sh / a benchmark run
# is in flight causes wall-clock timeout false positives (the quick suite runs
# close to its 1200s cap even on an idle machine) and poisons any measurement.
#
# Detection is process-based (pgrep) — zero coupling with the gate scripts, so
# it also protects against gates launched by OTHER sessions on this machine.
#
# Behaviour:
#   - Incoming command is heavy (ninja/cmake/ctest/test binary/benchmark) AND a
#     gate process is already running elsewhere -> DENY with explanation.
#   - Everything else -> allow silently.
#
# Reads tool_input JSON from stdin (Claude Code hook protocol). Never throws.

set -euo pipefail

CMD=$(python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('tool_input',{}).get('command',''))" 2>/dev/null || echo "")
[[ -z "$CMD" ]] && exit 0

# Read-only inspection commands are never gated, even when they MENTION a
# heavy binary (pgrep/ps/grep on process names, tail on logs, ...): the guard
# must not deny the very commands used to check whether a gate is running.
FIRST_TOKEN=$(echo "$CMD" | sed -E 's/^[[:space:]]*//' | awk '{print $1}' | xargs -I{} basename {} 2>/dev/null || echo "")
case "$FIRST_TOKEN" in
    pgrep|pkill|ps|grep|rg|cat|tail|head|wc|ls|echo|find|stat|file) exit 0 ;;
esac

# Only heavy commands are gated.
if ! echo "$CMD" | grep -qE '(ninja|cmake --build|ctest|cas_foundation_tests|test_quick\.sh|benchmark\.sh|benchmark_tests|run_golden_measurement\.sh)'; then
    exit 0
fi

# Gate/suite processes that must not run concurrently with new heavy work.
# pgrep -f matches full command lines of live processes.
GATE_PATTERN='test_quick\.sh|debt_gate\.sh|run_golden_measurement\.sh|benchmark\.sh|cas_foundation_tests'

RUNNING=$(pgrep -lf "$GATE_PATTERN" 2>/dev/null | grep -v "pgrep" || true)
[[ -z "$RUNNING" ]] && exit 0

SUMMARY=$(echo "$RUNNING" | head -3 | awk '{$1=""; print substr($0,2,120)}' | tr '\n' ';')

RUNNING_SUMMARY="$SUMMARY" python3 -c '
import json, os
print(json.dumps({
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason":
      "GATE-LOCK (memoria no-concurrent-gates): un gate/suite è GIÀ in esecuzione su questa macchina:\n  "
      + os.environ["RUNNING_SUMMARY"] + "\n"
      "Lanciare build/test ora causa timeout falso-positivi (quick suite vicina al cap 1200s) "
      "e misure inquinate. Attendi che il gate finisca (ricontrolla con pgrep separati: "
      "pgrep -lf test_quick ; pgrep -lf debt_gate ; pgrep -lf cas_foundation_tests) "
      "senza toccare il tree, poi rilancia."
  }
}))'
exit 0
