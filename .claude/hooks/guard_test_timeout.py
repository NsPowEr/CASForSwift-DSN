#!/usr/bin/env python3
# Logic for the PreToolUse "REGOLA TIMEOUT TEST" guard. Invoked by
# guard_test_timeout.sh so the hook's stdin (the tool_input JSON) reaches us
# intact — a heredoc (`python3 - <<PY`) would consume stdin with the source.
#
# Fires only when a project test/benchmark binary is the *executed* command of
# a pipeline segment (after stripping env assignments + a timeout wrapper), so
# commands that merely mention a binary (echo/grep/ls/gh) are not blocked.
import json, sys, os, re, shlex

try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(0)

cmd = (data.get("tool_input", {}) or {}).get("command", "") or ""
if not cmd.strip():
    sys.exit(0)

# Canonical self-capping wrapper — always allowed.
if "scripts/test_quick.sh" in cmd:
    sys.exit(0)

GTEST_BINS = {
    "cas_foundation_tests", "assumptions_stability_test", "cas_property_tests",
    "cas_gui_vm_tests", "cas_gui_qml_smoke_tests", "cas_tests",
}
BENCH_BINS = {"benchmark_tests", "cas_benchmarks"}

segments = re.split(r"&&|\|\||;|\||\n", cmd)

def tokenize(seg):
    try:
        return shlex.split(seg, comments=False, posix=True)
    except ValueError:
        return seg.split()

def first_exec(seg):
    toks = tokenize(seg)
    i = 0
    while i < len(toks) and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", toks[i]):
        i += 1  # leading env assignment
    if i < len(toks) and os.path.basename(toks[i]) in ("timeout", "gtimeout"):
        i += 1
        while i < len(toks) and (toks[i].startswith("-") or re.match(r"^[0-9]+[smhd]?$", toks[i])):
            i += 1  # timeout's options + duration
    if i < len(toks):
        return os.path.basename(toks[i]), toks[i:]
    return None, []

def has_timeout_wrapper(seg):
    toks = tokenize(seg)
    for j, t in enumerate(toks):
        if os.path.basename(t) in ("timeout", "gtimeout"):
            if j + 1 < len(toks) and re.match(r"^-|^[0-9]", toks[j + 1]):
                return True
    return False

def real_filter(seg):
    m = re.search(r"--gtest_filter=(\S*)", seg)
    if not m:
        return False
    return m.group(1).strip("'\"") not in ("", "*")

def deny(reason):
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": reason,
        }
    }))
    sys.exit(0)

for seg in segments:
    base, rest = first_exec(seg)
    if base is None:
        continue
    rest_str = " ".join(rest)

    if base in GTEST_BINS:
        if "--gtest_list_tests" in rest_str:
            continue
        missing = []
        if not has_timeout_wrapper(seg):
            missing.append("timeout(N)")
        if not real_filter(seg):
            missing.append("--gtest_filter=<reale, non vuoto e non '*'>")
        if missing:
            deny("REGOLA TIMEOUT TEST (CLAUDE.md): '%s' richiede %s, "
                 "oppure usa `bash scripts/test_quick.sh [--slow]`." %
                 (base, " + ".join(missing)))

    elif base in BENCH_BINS:
        if not has_timeout_wrapper(seg):
            deny("REGOLA TIMEOUT TEST (CLAUDE.md): il benchmark '%s' va lanciato "
                 "sotto `timeout <N>` per evitare hang silenziosi." % base)

    elif base == "ctest":
        if not re.search(r"--timeout\s+[0-9]+", seg) and not re.search(r"(^|\s)-R(\s|=)", seg):
            deny("REGOLA TIMEOUT TEST (CLAUDE.md): ctest richiede `--timeout <N>` "
                 "oppure `-R <regex>` per restringere/limitare l'esecuzione.")

sys.exit(0)
