#!/usr/bin/env bash
# debt_gate.sh — Single entry point for the repo's debt gates.
#
# Aggregates the four machine-verifiable gates that block invisible debt:
#   1. Anti-monolith        (check_file_size.sh)        — files > 500 LOC
#   2. Orphan sources       (check_orphan_sources.sh)   — src/ not wired in CMakeLists.txt
#   3. Untracked / stale    (debt_report.sh --strict)   — untracked .cpp + stale whitelist
#   4. HARDCODE-OF-PASSAGE  (inline grep)               — new HCOP comments without ledger entry
#
# Designed to be invoked from BOTH:
#   - .github/workflows/ci.yml     (non-bypassable wall on every PR / push)
#   - scripts/pre_commit_hook.sh   (local fast feedback)
#
# Usage:
#   bash scripts/debt_gate.sh                # scan tracked files only (CI default)
#   bash scripts/debt_gate.sh --staged       # scan staged diff (pre-commit semantics)
#
# Exit codes:
#   0 — all gates green
#   1 — at least one gate failed (block commit / fail PR)

set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2

MODE="tracked"
if [[ "${1:-}" == "--staged" ]]; then MODE="staged"; fi

fail=0
divider="═══════════════════════════════════════════════════════════"

# ── 1. Anti-monolith ─────────────────────────────────────────────────────────
echo "$divider"
echo "▸ Gate 1/4: anti-monolith (file size)"
echo "$divider"
if ! out=$(bash scripts/check_file_size.sh --verbose 2>&1); then
  echo "$out"
  echo "✗ Gate 1/4 FAILED: anti-monolith violation(s)."
  fail=1
else
  echo "$out"
  echo "✓ Gate 1/4: anti-monolith clean."
fi
echo

# ── 2. Orphan sources ────────────────────────────────────────────────────────
echo "$divider"
echo "▸ Gate 2/4: orphan sources (CMakeLists.txt wiring)"
echo "$divider"
if ! out=$(bash scripts/check_orphan_sources.sh 2>&1); then
  echo "$out"
  echo "✗ Gate 2/4 FAILED: orphan source(s) not wired in CMakeLists.txt."
  fail=1
else
  echo "✓ Gate 2/4: all sources wired in CMakeLists.txt."
fi
echo

# ── 3. Untracked + stale whitelist ───────────────────────────────────────────
echo "$divider"
echo "▸ Gate 3/4: untracked sources + stale whitelist (debt_report --strict)"
echo "$divider"
if ! out=$(bash scripts/debt_report.sh --strict 2>&1); then
  echo "$out"
  echo "✗ Gate 3/4 FAILED: untracked sources or stale whitelist entries."
  fail=1
else
  # debt_report --strict exits 0 when clean but still prints the full table
  # — only show the last summary block to keep CI logs compact.
  echo "$out" | tail -20
  echo "✓ Gate 3/4: no untracked sources, no stale whitelist."
fi
echo

# ── 4. HARDCODE-OF-PASSAGE in new diffs ──────────────────────────────────────
echo "$divider"
echo "▸ Gate 4/4: HARDCODE-OF-PASSAGE without ledger entry"
echo "$divider"

if [[ "$MODE" == "--staged" ]]; then
  DIFF_CMD=(git diff --cached -U0)
  STAGED=$(git diff --cached --name-only --diff-filter=ACMR)
else
  # For CI: scan every change introduced by this PR / push vs base.
  # Falls back to HEAD~1 if no base ref (e.g. push to first commit, or local run).
  # Refuse untrusted / malformed refs to keep the shell-quoted expansion safe.
  BASE="origin/${GITHUB_BASE_REF:-main}"
  if ! [[ "$BASE" =~ ^[A-Za-z0-9._/~-]+$ ]] || ! git rev-parse --verify "$BASE" >/dev/null 2>&1; then
    BASE="HEAD~1"
  fi
  DIFF_CMD=(git diff "$BASE"...HEAD -U0)
  STAGED=$(git diff "$BASE"...HEAD --name-only --diff-filter=ACMR 2>/dev/null || true)
fi

HCOP=$("${DIFF_CMD[@]}" -- 'src/**/*.cpp' 'src/**/*.hpp' 'src/**/*.h' 2>/dev/null \
       | grep -E '^\+' | grep -F 'HARDCODE-OF-PASSAGE' || true)

if [[ -n "$HCOP" ]]; then
  if echo "$STAGED" | grep -q 'HARDCODE_LEDGER.md'; then
    echo "✓ Gate 4/4: HARDCODE-OF-PASSAGE added with matching ledger entry."
  else
    echo "$HCOP" | sed 's/^/    /'
    echo "✗ Gate 4/4 FAILED: HARDCODE-OF-PASSAGE comment(s) introduced without"
    echo "  a HARDCODE_LEDGER.md update in the same commit/PR."
    echo "  See CLAUDE.md REGOLA ZERO and HARDCODE_LEDGER.md for the contract."
    fail=1
  fi
else
  echo "✓ Gate 4/4: no new HARDCODE-OF-PASSAGE comments."
fi
echo

# ── Verdict ──────────────────────────────────────────────────────────────────
echo "$divider"
if [[ $fail -ne 0 ]]; then
  echo "DEBT-GATE: BLOCKED — resolve the failing gate(s) above."
  echo "$divider"
  exit 1
fi
echo "DEBT-GATE: GREEN — all 4 gates passed."
echo "$divider"
exit 0
