#!/usr/bin/env bash
# check_golden_ratchet.sh — Ratchet sulla correttezza matematica (corpus golden).
#
# CONSUMA un report.json prodotto da cas_golden_runner (NON misura: la misura è il
# passo pesante separato che richiede Maxima). Confronta i PASS/FAIL aggregati con
# il baseline (scripts/golden_baseline.txt) e fallisce su regressione.
#
# Separazione netta (per non shippare orchestrazione non validata):
#   MISURA   → run_golden_maxima.sh + cas_golden_runner --json report.json   (pesante)
#   CONFRONTO→ questo script                                                   (leggero)
#
# Uso:
#   bash scripts/check_golden_ratchet.sh [--report <path>] [--baseline <path>]
#   bash scripts/check_golden_ratchet.sh --report report.json --update-baseline
#
# Default report: build-golden/report.json (lo scrivono CI / job periodico / misura
# manuale). Report assente → SKIP (non si può far rispettare ciò che non è misurato).
#
# Exit: 0 ok/skip/report-only, 1 regressione golden.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2

REPORT="build-golden/report.json"
BASELINE="scripts/golden_baseline.txt"
UPDATE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --report)         REPORT="$2"; shift 2 ;;
    --baseline)       BASELINE="$2"; shift 2 ;;
    --update-baseline) UPDATE=1; shift ;;
    *) echo "arg sconosciuto: $1" >&2; exit 2 ;;
  esac
done

if [[ ! -f "$REPORT" ]]; then
  echo "golden ratchet: report assente ($REPORT) → SKIP."
  echo "  Misura con: run_golden_maxima.sh + build/cas_golden_runner ... --json $REPORT"
  exit 0
fi

# Somma PASS/FAIL aggregati dal report.json (parse robusto via python).
read -r CUR_PASS CUR_FAIL CUR_UNDECIDED < <(python3 - "$REPORT" <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except Exception as e:
    print("ERR", "ERR", "ERR"); sys.exit(0)
p = sum(int(a.get("pass", 0)) for a in d.values() if isinstance(a, dict))
f = sum(int(a.get("fail", 0)) for a in d.values() if isinstance(a, dict))
# A51: entry NON DECISE = skip + over_budget. Nessuna delle due e' un
# verdetto sulla matematica (l'una perche' il confronto non e' concludente,
# l'altra perche' il motore e' stato troncato), e la RIPARTIZIONE fra le due
# dipende dal cap per-entry: misurato, due entry (~72 s) si spostano da OVER
# a SKIP passando da 60 s a 90 s. La SOMMA invece e' invariante (43 in
# entrambe le misure), ed e' percio' l'unica grandezza su cui un gate puo'
# poggiare senza tornare a dipendere dalla macchina.
# Va comunque vincolata: se cresce, il gate perde potere discriminante in
# silenzio — e chi volesse "passare" il ratchet potrebbe far scadere le entry
# scomode invece di ripararle.
o = sum(int(a.get("skip", 0)) + int(a.get("over_budget", 0))
        for a in d.values() if isinstance(a, dict))
print(p, f, o)
PY
)
if [[ "$CUR_PASS" == "ERR" ]]; then
  echo "✗ golden ratchet: $REPORT non è un JSON valido di cas_golden_runner."
  exit 1
fi
CUR_RATE="n/a"
(( CUR_PASS + CUR_FAIL > 0 )) && CUR_RATE=$(awk "BEGIN{printf \"%.1f%%\", 100*$CUR_PASS/($CUR_PASS+$CUR_FAIL)}")
echo "golden ratchet: misura corrente = $CUR_PASS pass / $CUR_FAIL fail / $CUR_UNDECIDED non decise ($CUR_RATE)"

get() { grep -E "^$1:" "$BASELINE" 2>/dev/null | head -1 | awk '{print $2}'; }
STATUS=$(get STATUS); FLOOR=$(get PASS_FLOOR); CEIL=$(get FAIL_CEILING); TOL=$(get TOLERANCE)
UNDEC_CEIL=$(get UNDECIDED_CEILING)
: "${STATUS:=UNSET}"; : "${FLOOR:=0}"; : "${CEIL:=0}"; : "${TOL:=0}"; : "${UNDEC_CEIL:=0}"

# --- aggiornamento baseline (ratchet tightening) ---
if [[ "$UPDATE" == "1" ]]; then
  tmp="$(mktemp)"
  awk -v p="$CUR_PASS" -v f="$CUR_FAIL" -v u="$CUR_UNDECIDED" '
    /^STATUS:/              {print "STATUS: SET"; next}
    /^PASS_FLOOR:/          {print "PASS_FLOOR: " p; next}
    /^FAIL_CEILING:/        {print "FAIL_CEILING: " f; next}
    /^UNDECIDED_CEILING:/   {print "UNDECIDED_CEILING: " u; seen=1; next}
    {print}
    END { if (!seen) print "UNDECIDED_CEILING: " u }
  ' "$BASELINE" > "$tmp" && mv "$tmp" "$BASELINE"
  echo "✓ baseline aggiornato: PASS_FLOOR=$CUR_PASS FAIL_CEILING=$CUR_FAIL UNDECIDED_CEILING=$CUR_UNDECIDED (STATUS=SET)."
  exit 0
fi

# --- report-only finché non misurato ---
if [[ "$STATUS" != "SET" ]]; then
  echo "  • STATUS=$STATUS → gate in SOLA LETTURA (baseline non ancora misurato)."
  echo "  • per attivare l'enforcement: --update-baseline da una misura reale."
  exit 0
fi

# --- enforcement ---
fail=0
if (( CUR_PASS < FLOOR - TOL )); then
  echo "✗ REGRESSIONE golden: PASS $CUR_PASS < floor $FLOOR (tol $TOL)."
  fail=1
fi
if (( CUR_FAIL > CEIL + TOL )); then
  echo "✗ REGRESSIONE golden: FAIL $CUR_FAIL > ceiling $CEIL (tol $TOL)."
  fail=1
fi
if (( CUR_UNDECIDED > UNDEC_CEIL + TOL )); then
  echo "✗ REGRESSIONE golden: NON DECISE $CUR_UNDECIDED > ceiling $UNDEC_CEIL (tol $TOL)."
  echo "  Entry che il motore non porta piu' a un verdetto (skip + over-budget):"
  echo "  e' una regressione, non un dato neutro. docs/rules/perf-root-cause.md."
  fail=1
fi
if (( fail == 0 )); then
  echo "  ✓ nessuna regressione (floor $FLOOR / ceiling $CEIL / non decise $UNDEC_CEIL, tol $TOL)."
  if (( CUR_PASS > FLOOR || CUR_FAIL < CEIL || CUR_UNDECIDED < UNDEC_CEIL )); then
    echo "  • miglioramento rilevato → consolida con --update-baseline."
  fi
fi
exit $fail
