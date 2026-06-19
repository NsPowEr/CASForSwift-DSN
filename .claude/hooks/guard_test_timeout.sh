#!/usr/bin/env bash
# PreToolUse Bash guard: enforce CLAUDE.md "REGOLA TIMEOUT TEST".
#
# Thin wrapper: the logic lives in guard_test_timeout.py so the hook's stdin
# (the tool_input JSON) reaches the parser intact. v2 detects the *executed*
# command of each pipeline segment (after stripping env + timeout wrapper) so
# commands that merely mention a test binary (echo/grep/ls/gh) are NOT blocked,
# while a real un-capped / un-filtered test run still is. Covers all project
# test + benchmark binaries, and rejects the fake filter `--gtest_filter=*`.
set -euo pipefail
exec python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/guard_test_timeout.py"
