#!/usr/bin/env bash
# PreToolUse Bash guard: mechanize the CLAUDE.md "REGOLA EVIDENCE-FIRST" for
# file deletion. Thin wrapper (same idiom as guard_test_timeout.sh) so the
# tool_input JSON on stdin reaches the Python parser intact.
#
# Denies `rm` / `find -delete` / `find -exec rm` against: git-tracked paths,
# protected source/infra directories, and absolute paths outside the repo
# (except temp/attic). Regenerable build output stays deletable.
set -euo pipefail
exec python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/guard_rm_safety.py"
