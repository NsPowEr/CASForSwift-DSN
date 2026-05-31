#!/usr/bin/env bash
# regenerate_maxima_manifest.sh
#
# Rigenera scripts/maxima_5.49.0_manifest.sha256 da una installazione clean.
# Uso: SOLO dopo brew install / brew reinstall maxima su sistema verificato.
# Vedi CLAUDE.md Regola 6.

set -euo pipefail

readonly EXPECTED_VERSION="5.49.0"
readonly EXPECTED_PREFIX="/opt/homebrew/Cellar/maxima/${EXPECTED_VERSION}_5"
readonly OUT="$(cd "$(dirname "$0")" && pwd)/maxima_5.49.0_manifest.sha256"

if [[ ! -d "$EXPECTED_PREFIX" ]]; then
  echo "ERROR: $EXPECTED_PREFIX not found. brew install maxima first." >&2
  exit 1
fi

cd "$EXPECTED_PREFIX"
{
  shasum -a 256 bin/maxima
  find share/maxima/${EXPECTED_VERSION}/src -name "*.lisp" -type f -print0 \
    | sort -z \
    | xargs -0 shasum -a 256
} > "$OUT"

echo "Wrote $(wc -l < "$OUT") entries to $OUT"
echo "Commit this manifest as the integrity baseline."
