#!/usr/bin/env bash
# check_quarantine_ratchet.sh — Ratchet gate sul debito test "quarantena".
#
# Garantisce che il numero di test noti-rossi (scripts/test_quarantine.txt) non
# CRESCA in silenzio: aggiungere un test richiede alzare deliberatamente CEILING.
# È il meccanismo che impedisce alla quarantena di diventare la prossima "lista
# di esclusione che cresce e nasconde il rosso".
#
# Regole:
#   - FAIL se (entry reali) > CEILING.
#   - FAIL se una entry ha id ledger vuoto (ogni rosso DEVE essere tracciato).
#   - ADVISORY se un id non è trovato in HARDCODE_LEDGER.md (può vivere in un PLAN).
#   - ADVISORY se (entry reali) < CEILING → abbassa CEILING per bloccare il progresso.
#
# Da eseguire a tempo di commit (idealmente cablato in debt_gate.sh / CI) e in
# check_dod.sh. Exit 0 = ok, 1 = ratchet violato.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2
FILE="${1:-scripts/test_quarantine.txt}"   # path override per test/flessibilità
LEDGER="HARDCODE_LEDGER.md"

[[ -f "$FILE" ]] || { echo "quarantine: $FILE assente — nessuna quarantena, ok."; exit 0; }

CEILING=$(grep -E '^CEILING:' "$FILE" | head -1 | awk '{print $2}')
if ! [[ "$CEILING" =~ ^[0-9]+$ ]]; then
    echo "✗ quarantine: CEILING mancante o non numerico in $FILE"
    exit 1
fi

count=0
empty_id=0
unledgered=""
while IFS= read -r line; do
    case "$line" in '#'*|'') continue ;; CEILING:*) continue ;; esac
    name=$(awk '{print $1}' <<< "$line")
    id=$(awk '{print $2}' <<< "$line")
    [[ -z "$name" ]] && continue
    count=$((count+1))
    if [[ -z "$id" || "$id" == "$name" ]]; then
        echo "  ✗ entry senza id ledger: $name"
        empty_id=$((empty_id+1))
        continue
    fi
    grep -qF "$id" "$LEDGER" 2>/dev/null || unledgered+=" $id"
done < "$FILE"

echo "quarantine ratchet: $count entry / CEILING $CEILING"

fail=0
if (( count > CEILING )); then
    echo "✗ RATCHET VIOLATO: $count > CEILING $CEILING."
    echo "  Per aggiungere un test alza CEILING in $FILE (atto deliberato e reviewabile),"
    echo "  oppure — preferito — fixa il test invece di nasconderlo (REGOLA 0.2)."
    fail=1
fi
if (( empty_id > 0 )); then
    echo "✗ $empty_id entry senza id ledger (ogni rosso quarantenato DEVE essere tracciato)."
    fail=1
fi
if [[ -n "$unledgered" ]]; then
    echo "  • advisory: id non trovati in $LEDGER (verifica che vivano in un PLAN):$unledgered"
fi
if (( count < CEILING && fail == 0 )); then
    echo "  • advisory: $count < CEILING $CEILING → abbassa CEILING a $count per bloccare il progresso."
fi

exit $fail
