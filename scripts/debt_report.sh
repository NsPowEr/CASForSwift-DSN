#!/usr/bin/env bash
# debt_report.sh — Weekly debt snapshot for the CAS Engine.
# Counts machine-verifiable signals of repo health and prints a table.
# Exit code: always 0 (informational). Use --strict to fail when over thresholds.

set -u
shopt -s nullglob

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

STRICT=0
if [[ "${1:-}" == "--strict" ]]; then STRICT=1; fi

# 1. Files > 500 LOC under src/
LARGE_FILES=$(find src -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.h" \) \
  -exec wc -l {} \; 2>/dev/null | awk '$1 > 500 {print}' | sort -rn | head -50)
LARGE_COUNT=$(echo "$LARGE_FILES" | grep -c . || true)

# 2. Whitelist entries (active uncommented src/ lines)
if [[ -f scripts/file_size_whitelist.txt ]]; then
  WL_COUNT=$(grep -E "^src/" scripts/file_size_whitelist.txt | wc -l | tr -d ' ')
else
  WL_COUNT=0
fi

# 3. Untracked .cpp/.hpp in src/  (must be 0 — every orphan is hidden debt)
UNTRACKED=$(git ls-files --others --exclude-standard -- 'src/**/*.cpp' 'src/**/*.hpp' 'src/**/*.h' 2>/dev/null || true)
UNTRACKED_COUNT=$(echo "$UNTRACKED" | grep -c . || true)

# 4. Sources in src/ NOT referenced in CMakeLists.txt
ALL_SRC=$(find src -type f -name "*.cpp" | sed 's|^|/|' | sort)
ORPHAN_SRC=""
while IFS= read -r f; do
  rel="${f#/}"
  if ! grep -qF "$rel" CMakeLists.txt 2>/dev/null; then
    ORPHAN_SRC="$ORPHAN_SRC$rel"$'\n'
  fi
done <<< "$ALL_SRC"
ORPHAN_COUNT=$(echo "$ORPHAN_SRC" | grep -c . || true)

# 5. Unimplemented / Unsupported tags in src/
UNIMPL=$(grep -rcE "CASErrorKind::Unimplemented|make_unimplemented" src/ 2>/dev/null | awk -F: '{s+=$2} END{print s+0}')

# 6. HARDCODE-OF-PASSAGE comments
HCPASS=$(grep -rcE "HARDCODE-OF-PASSAGE" src/ 2>/dev/null | awk -F: '{s+=$2} END{print s+0}')

# 7. TODO / FIXME comments
TODOS=$(grep -rcE "//\s*(TODO|FIXME|XXX|HACK)" src/ 2>/dev/null | awk -F: '{s+=$2} END{print s+0}')

# 8. HARDCODE_LEDGER open entries (status OPEN / PENDING / PARTIAL)
if [[ -f HARDCODE_LEDGER.md ]]; then
  LEDGER_OPEN=$(grep -cE "^\s*-\s*\*\*Status:\*\*\s*(OPEN|PENDING|PARTIAL|MITIGATED)" HARDCODE_LEDGER.md 2>/dev/null || echo 0)
else
  LEDGER_OPEN=0
fi

# 9. Disabled tests (DISABLED_, GTEST_SKIP without TODO)
DISABLED=$(grep -rcE "DISABLED_|GTEST_SKIP\(" test/ 2>/dev/null | awk -F: '{s+=$2} END{print s+0}')

# 10. Whitelist entries that no longer exceed 500 LOC (stale whitelist)
STALE_WL=0
if [[ -f scripts/file_size_whitelist.txt ]]; then
  while IFS= read -r path; do
    [[ -z "$path" || "$path" == \#* ]] && continue
    if [[ -f "$path" ]]; then
      loc=$(wc -l < "$path")
      if [[ $loc -le 500 ]]; then
        STALE_WL=$((STALE_WL+1))
      fi
    fi
  done < scripts/file_size_whitelist.txt
fi

STAMP=$(date +"%Y-%m-%d %H:%M:%S")
GIT_HEAD=$(git rev-parse --short HEAD 2>/dev/null || echo "?")

echo "═══════════════════════════════════════════════════════════"
echo "  CAS Engine — Debt Snapshot"
echo "  $STAMP  @ $GIT_HEAD"
echo "═══════════════════════════════════════════════════════════"
printf "%-44s %s\n" "files >500 LOC in src/"              "$LARGE_COUNT"
printf "%-44s %s\n" "  whitelisted entries"               "$WL_COUNT"
printf "%-44s %s\n" "  stale whitelist (now ≤500)"        "$STALE_WL"
printf "%-44s %s\n" "untracked sources in src/ ⚠"         "$UNTRACKED_COUNT"
printf "%-44s %s\n" "sources NOT referenced in CMakeLists" "$ORPHAN_COUNT"
printf "%-44s %s\n" "Unimplemented call sites"            "$UNIMPL"
printf "%-44s %s\n" "HARDCODE-OF-PASSAGE comments"        "$HCPASS"
printf "%-44s %s\n" "TODO/FIXME/XXX/HACK markers"         "$TODOS"
printf "%-44s %s\n" "HARDCODE_LEDGER open/pending"        "$LEDGER_OPEN"
printf "%-44s %s\n" "DISABLED_ tests / GTEST_SKIP"        "$DISABLED"
echo "═══════════════════════════════════════════════════════════"

if [[ "$UNTRACKED_COUNT" != "0" ]]; then
  echo
  echo "Untracked sources in src/ — bonifica obbligatoria:"
  echo "$UNTRACKED" | sed 's/^/  • /'
fi

if [[ "$ORPHAN_COUNT" != "0" ]]; then
  echo
  echo "Sources NOT wired in CMakeLists.txt:"
  echo "$ORPHAN_SRC" | sed 's/^/  • /'
fi

if [[ "$STRICT" == "1" ]]; then
  fail=0
  [[ "$UNTRACKED_COUNT" != "0" ]] && fail=1
  [[ "$ORPHAN_COUNT" != "0" ]] && fail=1
  [[ "$STALE_WL" != "0" ]] && fail=1
  exit $fail
fi

exit 0
