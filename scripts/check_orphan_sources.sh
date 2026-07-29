#!/usr/bin/env bash
# check_orphan_sources.sh — Fail if any .cpp in src/ is not referenced from CMakeLists.txt.
# Allowlist via scripts/orphan_sources_allow.txt (paths relative to repo root, one per line, # for comments).
# Used by the pre-commit hook to block invisible debt.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

ALLOW_FILE="scripts/orphan_sources_allow.txt"
ALLOW_LIST=""
if [[ -f "$ALLOW_FILE" ]]; then
  ALLOW_LIST=$(grep -vE '^\s*(#|$)' "$ALLOW_FILE")
fi

is_allowed() {
  local needle="$1"
  while IFS= read -r entry; do
    [[ -z "$entry" ]] && continue
    if [[ "$entry" == "$needle" ]]; then return 0; fi
  done <<< "$ALLOW_LIST"
  return 1
}

violations=0
while IFS= read -r f; do
  rel="${f#./}"
  if is_allowed "$rel"; then continue; fi
  # Reference check — exact path match in CMakeLists.txt.
  if ! grep -qF "$rel" CMakeLists.txt 2>/dev/null; then
    echo "ORPHAN: $rel — present in src/ but not in CMakeLists.txt"
    violations=$((violations+1))
  fi
done < <(find src -type f -name "*.cpp")

if [[ $violations -gt 0 ]]; then
  echo
  echo "Found $violations orphan source file(s)."
  echo "Either wire them in CMakeLists.txt, archive them, or add to $ALLOW_FILE with justification."
  exit 1
fi

exit 0
