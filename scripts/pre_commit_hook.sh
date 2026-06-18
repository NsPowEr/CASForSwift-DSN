#!/usr/bin/env bash
# pre_commit_hook.sh — Repo-tracked pre-commit gate.
# Install with: ln -sf ../../scripts/pre_commit_hook.sh .git/hooks/pre-commit
# Blocks commits that introduce invisible debt:
#   1. Untracked .cpp/.hpp/.h under src/ (orphans without a tracking decision)
#   2. .cpp under src/ not referenced in CMakeLists.txt (and not in allowlist)
#   3. New file >500 LOC under src/ not added to file_size_whitelist.txt
#   4. HARDCODE-OF-PASSAGE comments added to staged files without a matching
#      HARDCODE_LEDGER.md change in the same commit
# Override: `git commit --no-verify` (use only for emergency, never for agent commits).

set -u

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT" || exit 2

fail=0
divider="────────────────────────────────────────────────────────────"

# ── 1. Untracked source files in src/ ────────────────────────────────────────
UNTRACKED=$(git ls-files --others --exclude-standard -- 'src/**/*.cpp' 'src/**/*.hpp' 'src/**/*.h' 2>/dev/null || true)
if [[ -n "$UNTRACKED" ]]; then
  echo "$divider"
  echo "✗ PRE-COMMIT: untracked source(s) in src/ — bonifica required"
  echo "$UNTRACKED" | sed 's/^/    /'
  echo "  Action: add to git, archive (mv to .archive/), or .gitignore with justification."
  fail=1
fi

# ── 2. Orphan sources (not in CMakeLists.txt) ────────────────────────────────
if [[ -x scripts/check_orphan_sources.sh ]]; then
  if ! out=$(bash scripts/check_orphan_sources.sh 2>&1); then
    echo "$divider"
    echo "✗ PRE-COMMIT: orphan sources under src/ (not wired to CMakeLists.txt)"
    echo "$out" | sed 's/^/    /'
    fail=1
  fi
fi

# ── 3. Anti-monolith (delegates to existing script) ──────────────────────────
if [[ -x scripts/check_file_size.sh ]]; then
  if ! out=$(bash scripts/check_file_size.sh 2>&1); then
    echo "$divider"
    echo "✗ PRE-COMMIT: anti-monolith violations"
    echo "$out" | sed 's/^/    /'
    fail=1
  fi
fi

# ── 4. HARDCODE-OF-PASSAGE in staged diff must touch HARDCODE_LEDGER.md ──────
STAGED=$(git diff --cached --name-only --diff-filter=ACMR)
NEW_HC=$(git diff --cached -U0 -- 'src/**/*.cpp' 'src/**/*.hpp' 'src/**/*.h' 2>/dev/null \
        | grep -E '^\+' | grep -F 'HARDCODE-OF-PASSAGE' || true)
if [[ -n "$NEW_HC" ]]; then
  if ! echo "$STAGED" | grep -q 'HARDCODE_LEDGER.md'; then
    echo "$divider"
    echo "✗ PRE-COMMIT: new HARDCODE-OF-PASSAGE comment(s) without HARDCODE_LEDGER.md update"
    echo "$NEW_HC" | sed 's/^/    /'
    echo "  Action: add a corresponding entry to HARDCODE_LEDGER.md in this commit."
    fail=1
  fi
fi

# ── 5. Disabled tests should not be introduced silently ──────────────────────
NEW_DISABLED=$(git diff --cached -U0 -- 'test/**/*.cpp' 'test/**/*.hpp' 2>/dev/null \
              | grep -E '^\+' | grep -E 'DISABLED_|GTEST_SKIP\(|GTEST_DISABLE' || true)
if [[ -n "$NEW_DISABLED" ]]; then
  echo "$divider"
  echo "⚠ PRE-COMMIT WARNING: new DISABLED_/GTEST_SKIP introduced"
  echo "$NEW_DISABLED" | sed 's/^/    /'
  echo "  CLAUDE.md REGOLA 0.2: disabling tests requires explicit justification in commit body."
  echo "  Add 'TEST-DISABLED-JUSTIFY: <reason>' to commit message to acknowledge."
  COMMIT_MSG_FILE="$(git rev-parse --git-path COMMIT_EDITMSG)"
  if [[ -f "$COMMIT_MSG_FILE" ]] && ! grep -q "TEST-DISABLED-JUSTIFY:" "$COMMIT_MSG_FILE"; then
    fail=1
  fi
fi

if [[ $fail -ne 0 ]]; then
  echo "$divider"
  echo "Pre-commit blocked. Resolve the above, or use 'git commit --no-verify' for emergencies."
  exit 1
fi

exit 0
