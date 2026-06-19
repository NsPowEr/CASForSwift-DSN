#!/usr/bin/env bash
# pre_push_hook.sh — Tier-2 local correctness gate.
# Install: ln -sf ../../scripts/pre_push_hook.sh .git/hooks/pre-push
#
# PERCHÉ ESISTE: per un repo orientato all'autonomia la rete di sicurezza
# PORTANTE deve essere LOCALE, non la CI — che qui è gated sul billing e comunque
# lenta. I controlli strutturali economici stanno nel pre-commit; il segnale di
# correttezza più pesante (build + quick suite) sta QUI, al push (meno frequente
# del commit). Così l'agente resta onesto anche con la CI spenta.
#
# Stadi (fail-fast, dal più economico):
#   A. debt-gate (se presente)         — riprende il debito sfuggito a --no-verify
#   B. build incrementale di build/    — solo se build/ esiste; altrimenti skip
#   C. quick suite (quarantine-aware)  — il vero segnale di correttezza
#   D. golden ratchet (opt-in)         — CAS_PREPUSH_GOLDEN=1, richiede Maxima
#
# Bypass (emergenze, MAI per push agente):   git push --no-verify
# Salta solo questo gate:                     CAS_SKIP_PREPUSH=1 git push
# Anteprima senza buildare/testare:           CAS_PREPUSH_DRYRUN=1 ...

set -u
ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" || exit 0
cd "$ROOT" || exit 0

if [[ "${CAS_SKIP_PREPUSH:-0}" == "1" ]]; then
  echo "[pre-push] CAS_SKIP_PREPUSH=1 → gate saltato (assicurati sia voluto)."
  exit 0
fi
DRY="${CAS_PREPUSH_DRYRUN:-0}"

div="════════════════════════════════════════════════════════════"
fail=0
echo "$div"
echo "  pre-push gate (Tier-2 locale — la rete quando la CI è giù)"
[[ "$DRY" == "1" ]] && echo "  [DRY-RUN: nessuna build/test reale, solo orchestrazione]"
echo "$div"

run_or_dry() { # <label> <logfile> <cmd...>
  local label="$1" log="$2"; shift 2
  if [[ "$DRY" == "1" ]]; then echo "  → [dry] $label: $*"; return 0; fi
  "$@" >"$log" 2>&1
}

# ── A. debt-gate ─────────────────────────────────────────────────────────────
if [[ -x scripts/debt_gate.sh ]]; then
  if run_or_dry "debt-gate" /tmp/prepush_debt.log bash scripts/debt_gate.sh; then
    echo "  ✓ A debt-gate"
  else
    echo "  ✗ A debt-gate — vedi /tmp/prepush_debt.log"; fail=1
  fi
else
  echo "  ⊘ A debt-gate (scripts/debt_gate.sh assente)"
fi

# ── B. build incrementale + C. quick suite ───────────────────────────────────
if [[ -f build/build.ninja || -f build/Makefile ]]; then
  if run_or_dry "build cas_foundation_tests" /tmp/prepush_build.log cmake --build build --target cas_foundation_tests; then
    echo "  ✓ B build incrementale (build/)"
    if run_or_dry "quick suite" /tmp/prepush_test.log bash scripts/test_quick.sh; then
      echo "  ✓ C quick suite (quarantine-aware)"
    else
      echo "  ✗ C quick suite — vedi /tmp/prepush_test.log"; fail=1
    fi
  else
    echo "  ✗ B build FALLITO — vedi /tmp/prepush_build.log"; fail=1
  fi
else
  echo "  ⊘ B/C build+test saltati: build/ assente (cmake -B build -G Ninja per attivarli)"
fi

# ── D. golden ratchet (opt-in, richiede Maxima) ──────────────────────────────
if [[ "${CAS_PREPUSH_GOLDEN:-0}" == "1" ]]; then
  if [[ -x scripts/check_golden_ratchet.sh ]]; then
    if run_or_dry "golden ratchet" /tmp/prepush_golden.log bash scripts/check_golden_ratchet.sh; then
      echo "  ✓ D golden ratchet"
    else
      echo "  ✗ D golden ratchet — vedi /tmp/prepush_golden.log"; fail=1
    fi
  else
    echo "  ⊘ D golden ratchet (scripts/check_golden_ratchet.sh non ancora presente)"
  fi
fi

echo "$div"
if [[ $fail -ne 0 ]]; then
  echo "  PRE-PUSH BLOCCATO. Risolvi sopra, o (emergenza) git push --no-verify."
  echo "$div"
  exit 1
fi
echo "  PRE-PUSH OK — stato locale sano."
echo "$div"
exit 0
