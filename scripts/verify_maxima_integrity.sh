#!/usr/bin/env bash
# verify_maxima_integrity.sh
#
# Verifica che l'installazione di Maxima sia intatta (sorgenti NON modificati).
# Vedi CLAUDE.md Regola 6: Maxima è oracolo immutabile (GPL-2.0-only).
#
# Exit code:
#   0 = manifesto OK, integrità preservata
#   1 = mismatch (sorgenti modificati o manifesto stale)
#   2 = maxima non trovato o non eseguibile
#
# Output: log su stderr; OK/MISMATCH su stdout (CI-friendly).

set -euo pipefail

readonly EXPECTED_VERSION="5.49.0"
readonly EXPECTED_PREFIX="/opt/homebrew/Cellar/maxima/${EXPECTED_VERSION}_5"
readonly MANIFEST="$(cd "$(dirname "$0")" && pwd)/maxima_5.49.0_manifest.sha256"

log() { printf "[verify_maxima] %s\n" "$*" >&2; }

if ! command -v maxima >/dev/null 2>&1; then
  log "ERROR: maxima not found in PATH"
  exit 2
fi

ACTUAL_VERSION=$(maxima --version 2>&1 | head -1 | awk '{print $2}')
if [[ "$ACTUAL_VERSION" != "$EXPECTED_VERSION" ]]; then
  log "ERROR: Maxima version mismatch (expected $EXPECTED_VERSION, got $ACTUAL_VERSION)"
  log "PLAN_HP_PRIME_PARITY.md pins $EXPECTED_VERSION. Aggiornare manifesto o re-install."
  exit 1
fi

if [[ ! -d "$EXPECTED_PREFIX" ]]; then
  log "ERROR: Maxima Cellar prefix not found: $EXPECTED_PREFIX"
  exit 2
fi

if [[ ! -f "$MANIFEST" ]]; then
  log "ERROR: Manifest not found: $MANIFEST"
  log "Run scripts/regenerate_maxima_manifest.sh on a clean install to bootstrap."
  exit 2
fi

log "Verifying Maxima integrity against pinned manifest ($MANIFEST)..."
cd "$EXPECTED_PREFIX"
if shasum -a 256 -c "$MANIFEST" --status; then
  log "OK: Maxima ${EXPECTED_VERSION} integrity verified ($(wc -l < "$MANIFEST") entries)"
  printf "OK\n"
  exit 0
else
  log "MISMATCH: Maxima sources or binary differ from pinned manifest!"
  log "This indicates either (a) tampered installation, (b) drift via brew upgrade,"
  log "or (c) accidental modification. Per CLAUDE.md Regola 6, modification is FORBIDDEN."
  shasum -a 256 -c "$MANIFEST" 2>&1 | grep -v ": OK$" | head -20 >&2
  printf "MISMATCH\n"
  exit 1
fi
