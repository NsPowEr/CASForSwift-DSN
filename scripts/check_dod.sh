#!/usr/bin/env bash
# check_dod.sh — Mechanizable Definition-of-Done gates for the CAS Engine.
# Runs the 🤖 gates from DEFINITION_OF_DONE.md (G1,G3,G4,G7,G8,G9,G10).
# Human/agent gates (G2 generality, G5 test quality, G6 golden oracle) are
# reported as reminders — they require judgment and live in the commit/PR.
#
# Usage:
#   bash scripts/check_dod.sh            # quick gates only (skips slow G7/G8)
#   bash scripts/check_dod.sh --full     # includes test_quick.sh + benchmark.sh
#
# Exit 0 if all run gates pass, 1 otherwise.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

FULL=0
[[ "${1:-}" == "--full" ]] && FULL=1

pass=0; fail=0; skip=0
div="────────────────────────────────────────────────────────────"

ok()   { printf "  ✓ %s\n" "$1"; pass=$((pass+1)); }
ko()   { printf "  ✗ %s\n" "$1"; fail=$((fail+1)); }
note() { printf "  • %s\n" "$1"; }
skipg(){ printf "  ⊘ %s\n" "$1"; skip=$((skip+1)); }

echo "$div"
echo "  CAS Engine — Definition of Done gate check"
echo "  $(date +%H:%M:%S)  @ $(git rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "$div"

# ── G3 — every HARDCODE-OF-PASSAGE id in src/ must exist in the ledger ────────
# Real gate: extract the <id> from each `HARDCODE-OF-PASSAGE: <id>` comment and
# verify it is present in HARDCODE_LEDGER.md. (The old check only failed when the
# ledger had *zero* entries — a tautology, since it has hundreds.)
echo "G3  Hardcode discipline"
HC_IDS=$(grep -rhoE 'HARDCODE-OF-PASSAGE:[[:space:]]*[A-Za-z0-9._-]+' src/ 2>/dev/null \
         | sed -E 's/.*:[[:space:]]*//' | sort -u)
HC_UNLABELED=$(grep -rE 'HARDCODE-OF-PASSAGE' src/ 2>/dev/null \
               | grep -vE 'HARDCODE-OF-PASSAGE:[[:space:]]*[A-Za-z0-9._-]+' | wc -l | tr -d ' ')
if [[ -z "$HC_IDS" ]]; then
  ok "no labeled HARDCODE-OF-PASSAGE ids in src/"
else
  MISSING=""
  while IFS= read -r id; do
    [[ -z "$id" ]] && continue
    grep -qF "$id" HARDCODE_LEDGER.md 2>/dev/null || MISSING+=" $id"
  done <<< "$HC_IDS"
  if [[ -n "$MISSING" ]]; then
    ko "HARDCODE-OF-PASSAGE id(s) in src/ NOT found in HARDCODE_LEDGER.md:$MISSING"
  else
    ok "all $(echo "$HC_IDS" | grep -c .) HARDCODE-OF-PASSAGE id(s) are ledgered"
  fi
fi
[[ "$HC_UNLABELED" -gt 0 ]] && note "HARDCODE-OF-PASSAGE comments without a parseable id: $HC_UNLABELED (assegna un id ledgered)"

# ── G4 — no int64_t/double in symbolic core (ADVISORY, non-blocking) ──────────
# Heuristic only: many hits are legitimate numeric bridges. Reported as advisory
# so it never inflates the pass count nor implies a guarantee it cannot give.
echo "G4  Exact arithmetic in core (advisory)"
RAW_DBL=$(grep -rnE '\b(double|int64_t|long long)\b' src/symbolic src/algebra src/calculus 2>/dev/null \
          | grep -vE '//|/\*|numeric|NumericEval|\.benchmark|to_double|as_double|approx' | wc -l | tr -d ' ')
if [[ "$RAW_DBL" -gt 0 ]]; then
  note "candidate double/int64 in symbolic dirs: $RAW_DBL (advisory — revisione manuale, molti sono ponti numerici legittimi)"
else
  note "no raw double/int64 in symbolic core dirs (advisory)"
fi

# ── G9 — build hygiene: anti-monolith + orphan sources ───────────────────────
echo "G9  Build hygiene"
if [[ -x scripts/check_file_size.sh ]]; then
  if out=$(bash scripts/check_file_size.sh 2>&1); then ok "anti-monolith: no new violations"
  else ko "anti-monolith violations"; echo "$out" | sed 's/^/      /'; fi
else skipg "check_file_size.sh not found"; fi

# ── G10 — tracker hygiene: orphans + unknown statuses ────────────────────────
echo "G10 Tracker hygiene"
if [[ -x scripts/check_orphan_sources.sh ]]; then
  if out=$(bash scripts/check_orphan_sources.sh 2>&1); then ok "no orphan sources outside allowlist"
  else ko "orphan sources present"; echo "$out" | sed 's/^/      /'; fi
else skipg "check_orphan_sources.sh not found"; fi

UNTRACKED=$(git ls-files --others --exclude-standard -- 'src/**/*.cpp' 'src/**/*.hpp' 'src/**/*.h' 2>/dev/null | grep -c . || true)
if [[ "$UNTRACKED" -eq 0 ]]; then ok "no untracked sources in src/"
else ko "untracked sources in src/: $UNTRACKED"; fi

if [[ -x scripts/tasks_audit.sh ]]; then
  UNK=$(bash scripts/tasks_audit.sh --tsv 2>/dev/null | awk -F'\t' '$4=="UNKNOWN"' | wc -l | tr -d ' ')
  note "task ledger UNKNOWN statuses: $UNK (bonificare quando tocchi il file)"
fi

# ── G7 — regression suite (slow, opt-in) ─────────────────────────────────────
echo "G7  Regression suite"
if [[ "$FULL" -eq 1 ]]; then
  if [[ -x scripts/test_quick.sh ]]; then
    if bash scripts/test_quick.sh >/tmp/dod_testquick.log 2>&1; then ok "test_quick.sh green"
    else ko "test_quick.sh failed — see /tmp/dod_testquick.log"; fi
  else skipg "test_quick.sh not found"; fi
else
  skipg "G7 skipped (run with --full)"
fi

# ── G8 — benchmark gate (slow, opt-in) ───────────────────────────────────────
echo "G8  Benchmark gate"
if [[ "$FULL" -eq 1 ]]; then
  if [[ -x scripts/check_bench_regression.sh ]]; then
    if bash scripts/check_bench_regression.sh >/tmp/dod_bench.log 2>&1; then ok "no benchmark regression"
    else ko "benchmark regression — see /tmp/dod_bench.log"; fi
  else skipg "check_bench_regression.sh not found"; fi
else
  skipg "G8 skipped (run with --full)"
fi

# ── Judgment gates reminder ──────────────────────────────────────────────────
echo "Manual gates (document in commit/PR):"
note "G1  Spec read (MISSING_FEATURES_SPECS/<task>.md)"
note "G2  General algorithm, no shortcut (REGOLA ZERO)"
note "G5  Tests: distinct vars + equivalent forms + anti-hardcode"
note "G6  Golden diff vs Maxima 5.49.0 green on declared domain"

echo "$div"
printf "  pass: %d   fail: %d   skipped: %d\n" "$pass" "$fail" "$skip"
echo "$div"

[[ "$fail" -eq 0 ]] && exit 0 || exit 1
