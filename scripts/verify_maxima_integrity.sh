#!/usr/bin/env bash
# verify_maxima_integrity.sh
#
# Verifica che l'installazione di Maxima sia intatta (sorgenti NON modificati).
# Vedi CLAUDE.md Regola 6: Maxima è oracolo immutabile (GPL-2.0-only).
#
# Il pin normativo è sulla VERSIONE upstream (5.49.0). La revision del bottle
# Homebrew (_5, _6, …) è packaging: un bump rilinka le dipendenze e riscrive
# il prefix Cellar embedded in due file (bin/maxima, autoconf-variables.lisp)
# senza toccare alcun sorgente. Quei due file sono quindi verificati sul
# contenuto NORMALIZZATO (prefix → @MAXIMA_CELLAR@); tutti gli altri sul
# contenuto letterale. Verificato empiricamente al bump _5→_6 (2026-07-16):
# 308/310 hash identici, 2 differenti solo per la stringa del prefix.
#
# Exit code:
#   0 = manifesto OK, integrità preservata
#   1 = mismatch (sorgenti modificati o manifesto stale)
#   2 = maxima non trovato / prefix ambiguo o assente
#
# Output: log su stderr; OK/MISMATCH su stdout (CI-friendly).

set -euo pipefail

readonly EXPECTED_VERSION="5.49.0"
readonly CELLAR_GLOB="/opt/homebrew/Cellar/maxima/${EXPECTED_VERSION}"
readonly MANIFEST="$(cd "$(dirname "$0")" && pwd)/maxima_5.49.0_manifest.sha256"
# File noti col prefix Cellar embedded nel contenuto (set chiuso deliberato:
# lo shell-wrapper e le variabili autoconf; nessun sorgente algoritmico).
readonly PATH_EMBEDDED_FILES=(
  "bin/maxima"
  "share/maxima/${EXPECTED_VERSION}/src/autoconf-variables.lisp"
)

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

# Risolvi la revision bottle corrente (5.49.0, 5.49.0_5, 5.49.0_6, …).
# Più prefix installati = ambiguità → errore esplicito, mai scelta silenziosa.
shopt -s nullglob
prefixes=( "$CELLAR_GLOB" "$CELLAR_GLOB"_* )
shopt -u nullglob
existing=()
for p in ${prefixes[@]+"${prefixes[@]}"}; do
  [[ -d "$p" ]] && existing+=("$p")
done
if [[ ${#existing[@]} -eq 0 ]]; then
  log "ERROR: no Maxima Cellar prefix found matching ${CELLAR_GLOB}[_*]"
  exit 2
fi
if [[ ${#existing[@]} -gt 1 ]]; then
  log "ERROR: multiple Maxima Cellar prefixes found: ${existing[*]}"
  log "Ambiguous install; remove stale revisions (brew cleanup maxima) and retry."
  exit 2
fi
readonly PREFIX="${existing[0]}"
log "Cellar prefix: $PREFIX"

if [[ ! -f "$MANIFEST" ]]; then
  log "ERROR: Manifest not found: $MANIFEST"
  log "Run scripts/regenerate_maxima_manifest.sh on a clean install to bootstrap."
  exit 2
fi

is_path_embedded() {
  local f="$1"
  for pe in "${PATH_EMBEDDED_FILES[@]}"; do
    [[ "$f" == "$pe" ]] && return 0
  done
  return 1
}

log "Verifying Maxima integrity against pinned manifest ($MANIFEST)..."
cd "$PREFIX"

status=0

# 1. File a contenuto letterale: shasum -c sul manifest depurato dei path-embedded.
literal_manifest="$(mktemp)"
trap 'rm -f "$literal_manifest"' EXIT
while IFS= read -r line; do
  f="${line#*  }"
  is_path_embedded "$f" || printf '%s\n' "$line"
done < "$MANIFEST" > "$literal_manifest"
if ! shasum -a 256 -c "$literal_manifest" --status; then
  log "MISMATCH in literal-content files:"
  shasum -a 256 -c "$literal_manifest" 2>&1 | grep -v ": OK$" | head -20 >&2
  status=1
fi

# 2. File path-embedded: hash del contenuto normalizzato (prefix → placeholder).
for f in "${PATH_EMBEDDED_FILES[@]}"; do
  expected="$(awk -v f="$f" '$0 ~ "  "f"$" {print $1}' "$MANIFEST")"
  if [[ -z "$expected" ]]; then
    log "MISMATCH: manifest has no entry for path-embedded file $f"
    status=1
    continue
  fi
  if [[ ! -f "$f" ]]; then
    log "MISMATCH: path-embedded file missing on disk: $f"
    status=1
    continue
  fi
  actual="$(sed "s|$PREFIX|@MAXIMA_CELLAR@|g" "$f" | shasum -a 256 | awk '{print $1}')"
  if [[ "$actual" != "$expected" ]]; then
    log "MISMATCH (normalized content): $f"
    status=1
  fi
done

if [[ $status -eq 0 ]]; then
  log "OK: Maxima ${EXPECTED_VERSION} integrity verified ($(wc -l < "$MANIFEST") entries)"
  printf "OK\n"
  exit 0
fi

log "MISMATCH: Maxima sources or binary differ from pinned manifest!"
log "This indicates either (a) tampered installation, (b) drift via brew upgrade,"
log "or (c) accidental modification. Per CLAUDE.md Regola 6, modification is FORBIDDEN."
printf "MISMATCH\n"
exit 1
