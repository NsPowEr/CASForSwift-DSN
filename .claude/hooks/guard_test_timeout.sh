#!/usr/bin/env bash
# PreToolUse Bash guard: enforce CLAUDE.md "REGOLA TIMEOUT TEST".
#
# Blocks raw invocations of cas_foundation_tests (or any GoogleTest binary in
# build*/ trees) that lack BOTH --gtest_filter AND an explicit shell timeout.
#
# Allowed forms:
#   - bash scripts/test_quick.sh [--slow]
#   - timeout <N> ... cas_foundation_tests --gtest_filter=...
#   - ctest --timeout <N>
#
# Reads tool_input JSON from stdin (Claude Code hook protocol).

set -euo pipefail

CMD=$(python3 -c "import json,sys; d=json.load(sys.stdin); print(d.get('tool_input',{}).get('command',''))" 2>/dev/null || echo "")

# Only inspect commands that mention the test binary directly.
if ! echo "$CMD" | grep -qE '(cas_foundation_tests|ctest)'; then
    exit 0
fi

# Allow the canonical wrapper.
if echo "$CMD" | grep -qE 'scripts/test_quick\.sh'; then
    exit 0
fi

has_filter=0
has_timeout=0
echo "$CMD" | grep -qE -- '--gtest_filter=' && has_filter=1
echo "$CMD" | grep -qE '(^|[[:space:]])timeout[[:space:]]+[0-9]+' && has_timeout=1
echo "$CMD" | grep -qE 'ctest.*--timeout[[:space:]]+[0-9]+' && has_timeout=1

if [[ $has_filter -eq 1 && $has_timeout -eq 1 ]]; then
    exit 0
fi

REASON="REGOLA TIMEOUT TEST (CLAUDE.md): cas_foundation_tests richiede --gtest_filter + timeout esplicito,"
REASON+=" oppure usare \`bash scripts/test_quick.sh [--slow]\`."
REASON+=" Mancano: "
[[ $has_filter -eq 0 ]] && REASON+="--gtest_filter "
[[ $has_timeout -eq 0 ]] && REASON+="timeout(N) "

python3 -c "
import json, sys
print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'PreToolUse',
        'permissionDecision': 'deny',
        'permissionDecisionReason': '''$REASON'''
    }
}))
"
exit 0
