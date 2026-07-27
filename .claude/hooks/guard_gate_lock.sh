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

# A command is heavy only if a heavy program is the COMMAND WORD of one of its
# segments — never because a heavy name appears as text (grep/sed/echo that
# merely mention a gate script must pass: 2026-07-27 a read-only multi-line
# inspection was denied because the old first-token allowlist read only line 1
# and the fallback regex matched script names cited as arguments).
HEAVY=$(CMD="$CMD" python3 - <<'PYEOF' 2>/dev/null || echo "OK"
import os, re, shlex
cmd = os.environ.get("CMD", "")

def split_segments(text):
    # Quote-aware split on unquoted ; | & ( ) and newline — a naive regex split
    # would break quoted strings (e.g. the | alternation inside a pgrep
    # pattern) and turn their pieces into phantom command words.
    segs, cur, state, i, n = [], [], None, 0, len(text)
    while i < n:
        c = text[i]
        if state == "'":
            cur.append(c)
            if c == "'":
                state = None
        elif state == '"':
            cur.append(c)
            if c == "\\" and i + 1 < n:
                cur.append(text[i + 1]); i += 1
            elif c == '"':
                state = None
        elif c == "\\" and i + 1 < n:
            cur.append(c); cur.append(text[i + 1]); i += 1
        elif c in ("'", '"'):
            cur.append(c); state = c
        elif c in ";|&()\n":
            segs.append("".join(cur)); cur = []
        else:
            cur.append(c)
        i += 1
    segs.append("".join(cur))
    return [s for s in segs if s.strip()]

segments = split_segments(cmd)
heavy_base = re.compile(
    r'^(ninja|ctest)$'
    r'|^(cas_foundation_tests|benchmark_tests)$'
    r'|^(test_quick|benchmark|run_golden_measurement|check_golden_ratchet|debt_gate)\.sh$')
wrappers = {"env", "nice", "time", "caffeinate", "stdbuf", "xvfb-run"}
def is_heavy(seg):
    try:
        toks = shlex.split(seg.strip())
    except ValueError:
        toks = seg.strip().split()
    while toks and (re.match(r'^[A-Za-z_][A-Za-z0-9_]*=', toks[0]) or toks[0] in wrappers):
        toks = toks[1:]
    if toks and toks[0] == "timeout":
        toks = toks[2:]
    if toks and toks[0] in ("bash", "sh", "zsh"):
        toks = toks[1:]
    if not toks:
        return False
    base = toks[0].rsplit("/", 1)[-1]
    if base == "cmake":
        return "--build" in toks
    return bool(heavy_base.match(base))
print("HEAVY" if any(is_heavy(s) for s in segments) else "OK")
PYEOF
)
[[ "$HEAVY" != "HEAVY" ]] && exit 0

# Gate/suite processes that must not run concurrently with new heavy work.
# pgrep -f matches full command lines of live processes.
GATE_PATTERN='test_quick\.sh|debt_gate\.sh|run_golden_measurement\.sh|check_golden_ratchet\.sh|benchmark\.sh|cas_foundation_tests'

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
