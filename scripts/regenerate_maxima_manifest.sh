#!/usr/bin/env bash
# regenerate_maxima_manifest.sh
#
# Rigenera scripts/maxima_5.49.0_manifest.sha256 da una installazione clean.
# Uso: SOLO dopo brew install / brew reinstall maxima su sistema verificato.
# Vedi CLAUDE.md Regola 6.
#
# I due file col prefix Cellar embedded (bin/maxima, autoconf-variables.lisp)
# sono hashati sul contenuto NORMALIZZATO (prefix → @MAXIMA_CELLAR@) così il
# manifest sopravvive ai bump di revision del bottle Homebrew (_5 → _6, …)
# che riscrivono solo quel path. La stessa normalizzazione è applicata in
# verify_maxima_integrity.sh.

set -euo pipefail

readonly EXPECTED_VERSION="5.49.0"
readonly CELLAR_GLOB="/opt/homebrew/Cellar/maxima/${EXPECTED_VERSION}"
readonly OUT="$(cd "$(dirname "$0")" && pwd)/maxima_5.49.0_manifest.sha256"
readonly PATH_EMBEDDED_FILES=(
  "bin/maxima"
  "share/maxima/${EXPECTED_VERSION}/src/autoconf-variables.lisp"
)

shopt -s nullglob
prefixes=( "$CELLAR_GLOB" "$CELLAR_GLOB"_* )
shopt -u nullglob
existing=()
for p in ${prefixes[@]+"${prefixes[@]}"}; do
  [[ -d "$p" ]] && existing+=("$p")
done
if [[ ${#existing[@]} -ne 1 ]]; then
  echo "ERROR: need exactly one Cellar prefix matching ${CELLAR_GLOB}[_*], found: ${existing[*]:-none}" >&2
  exit 1
fi
readonly PREFIX="${existing[0]}"

cd "$PREFIX"
{
  # Path-embedded files: normalized-content hash.
  for f in "${PATH_EMBEDDED_FILES[@]}"; do
    h="$(sed "s|$PREFIX|@MAXIMA_CELLAR@|g" "$f" | shasum -a 256 | awk '{print $1}')"
    printf '%s  %s\n' "$h" "$f"
  done
  # Everything else: literal content.
  find share/maxima/${EXPECTED_VERSION}/src -name "*.lisp" -type f -print0 \
    | sort -z \
    | xargs -0 shasum -a 256 \
    | grep -v "share/maxima/${EXPECTED_VERSION}/src/autoconf-variables\\.lisp"
} > "$OUT"

echo "Wrote $(wc -l < "$OUT") entries to $OUT (prefix: $PREFIX)"
echo "Commit this manifest as the integrity baseline."
