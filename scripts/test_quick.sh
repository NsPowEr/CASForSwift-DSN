#!/usr/bin/env bash
# test_quick.sh — esegue la suite "quick" con esclusioni governate e protezione
# anti-hang via timeout. Tre categorie di esclusione, semanticamente distinte:
#
#   STRUCTURAL  *Stress* / *Fuzz* / *Disabled*  → esclusi sempre (eseguiti
#               altrove: Stress dal regression-guard, Fuzz via libFuzzer,
#               Disabled è convenzione gtest). Legittimo.
#   SLOW_OK     test che PASSANO ma sono lenti (>30s). Esclusi dal quick;
#               INCLUSI in --slow (cap 1800s). Nessun debito: solo lenti.
#   QUARANTINE  test NOTI-ROSSI (fail/hang) → scripts/test_quarantine.txt.
#               Esclusi da quick E slow (altrimenti il loro hang rende il gate
#               inusabile). NON sono nascosti: banner ad ogni run + ledger id +
#               ratchet (check_quarantine_ratchet.sh). Debito visibile e capped.
#
# Differenza chiave vs la versione precedente: prima quarantena e lenti-legittimi
# erano mescolati nella stessa lista, e --slow li re-includeva → il gate "serio"
# andava in hang e nessuno lo lanciava, mentre il quick dava "verde" nascondendo
# regressioni. Ora il verde del quick è onesto e --slow è di nuovo eseguibile.
#
# Uso:
#   bash scripts/test_quick.sh                 # quick, cap 1200s
#   bash scripts/test_quick.sh --slow          # + SLOW_OK, cap 1800s (gate pre-push)
#   bash scripts/test_quick.sh --filter X      # filtro positivo aggiuntivo
#   bash scripts/test_quick.sh --quarantine    # esegue SOLO i quarantenati, uno per
#                                              #   uno con timeout stretto: l'hang
#                                              #   diventa un FAIL veloce e visibile
#   bash scripts/test_quick.sh --print-filter  # dry-run: stampa filtro+banner, no run
#   bash scripts/test_quick.sh --print-ctest-regex  # stessa esclusione, sintassi
#                                              #   regex per `ctest -E` (CI: la CI
#                                              #   gira via ctest, non via binario
#                                              #   diretto — un solo elenco condiviso,
#                                              #   niente lista duplicata in ci.yml)
#
# Se introduci una regressione di complessità che porta un test sopra il cap, NON
# aggiungerlo qui: indaga la causa (probabile O(2^n) accidentale) e fixa il codice.
# Se è davvero rosso e non fixabile subito, va in test_quarantine.txt con un id di
# ledger e alzando CEILING — un atto deliberato, non un'esclusione silenziosa.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${REPO_ROOT}/build/cas_foundation_tests"
QUAR_FILE="${REPO_ROOT}/scripts/test_quarantine.txt"

# ── Categorie ────────────────────────────────────────────────────────────────
STRUCTURAL=( '*Stress*' '*Fuzz*' '*Disabled*' )

# Lenti ma verdi (eseguiti in --slow). Tempi su macOS arm64.
SLOW_OK=(
    VanHoeijDirect.Deg16_EightQuadratics_FindsRealFactor                  # 181s
    FactorizationTowerNTest.SplitsX2Minus5_Over_Q_Sqrt2_Sqrt3_Sqrt5       # >240s
    PrimitiveElementTest.RedundantMixedTower_Sqrt2_Sqrt3_Sqrt5_Sqrt6      # >60s (oggi DISABLED_)
    PrimitiveElementTest.DetectNLevelTower_MultiBetaNested                # ~30s (macOS ASan)
    GaloisDeg5Test.C5_RealCyclotomic11                                    # ~360s (macOS ASan)
    PermMaximalTest.CrossCheckDegree7BothAmbients                         # ~45s (2× lattice n=7, macOS ASan)
    GaloisStauduharTest.CrossCheckDegree6                                 # ~45s (oracolo = pipeline deg6, macOS ASan)
    GaloisStauduharTest.CrossCheckDegree7                                 # ~57s (oracolo = pipeline deg7, macOS ASan)
    GaloisStauduharTest.IdentifyDegree8X8Minus2Order16                    # ~25s isolato (catena deg-8 completa, macOS ASan) — borderline sotto carico
    GaloisStauduharTest.IdentifyDegree10FullWreathViaStructuralRoute      # ~17s isolato (x¹⁰−2 → ordine 40, route strutturale Brick 3.75, macOS ASan)
    GaloisWreathMaximalTest.ScottLemmaGroundTruthOnA5xA5                  # ~20s isolato (ground truth Scott su A₅×A₅ ordine 3600, macOS ASan)
    GaloisDeg8E2E.Deg8Irreducible                                        # ~42s (corpus Brick 4: S8 + 5 ciclotomiche + x⁸−2, macOS ASan)
    GaloisDeg8E2E.Deg9And10Irreducible                                   # ~26s (S9 + x¹⁰ full-wreath via galois_group pubblico, macOS ASan)
    GaloisDeg8E2E.ReducibleRecursionDeg8                                 # ~8s (ricorsione fattori deg-8, macOS ASan)
    GaloisDeg8E2E.ReducibleRecursionDeg9                                 # ~21s (ricorsione fattori deg-9, macOS ASan)
    GaloisDeg8E2E.ReducibleRecursionDeg10                                # ~22s (ricorsione fattori deg-10, macOS ASan)
    # 2026-07-16: i 4 test sotto (>100s ciascuno, verdi) assorbivano il 49%
    # del quick (574s/1176s) portandolo al 98% del cap 1200s. Tempi verificati
    # ISOLATI (WeierstrassDirectProbe 196s da solo → intrinseco, non carico).
    # WeierstrassDirectProbe è un probe diagnostico senza asserzioni: il suo
    # costo vive in weierstrass_zero_diff/simplify → candidato indagine perf.
    IntegrateInverseTrigTest.WeierstrassDirectProbe                       # ~196-222s (weierstrass_zero_diff su diff di antiderivate)
    IntegrationAdvancedTest.RischDE_ExpQuadraticArg_HigherDegree          # ~142s (RischDE grado alto, macOS ASan)
    DefiniteSummationTest.QuadraticIrreducibleLinearNum_ViaRootOfDigamma  # ~108s (RootOf digamma, macOS ASan)
    DefiniteSummationTest.QuadraticIrreducibleDenom_ViaRootOfDigamma      # ~101s (RootOf digamma, macOS ASan)
)

# Quarantena (noti-rossi) — caricata dal file governato.
QUARANTINE=()
if [[ -f "$QUAR_FILE" ]]; then
    while IFS= read -r line; do
        case "$line" in
            '#'*|'') continue ;;
            CEILING:*) continue ;;
        esac
        name=$(awk '{print $1}' <<< "$line")
        [[ -n "$name" ]] && QUARANTINE+=("$name")
    done < "$QUAR_FILE"
fi

join_colon() { local IFS=':'; echo "$*"; }

# ── Parsing argomenti ────────────────────────────────────────────────────────
# Due dimensioni ORTOGONALI (tenerle separate evita che --slow e --print-filter
# si sovrascrivano): SUITE = cosa escludere; ACTION = cosa fare.
POSITIVE_FILTER=''
# Misura reale suite quick: ~1178s @2729 test (2026-07-16, isolata) PRIMA
# dello spostamento dei 4 test >100s in SLOW_OK → attesa ~600s dopo. Cap =
# misura + ~50% margine anti-flake da carico; se la suite sfora QUESTO cap,
# è una regressione di complessità da indagare, non un numero da alzare
# alla cieca. (Storico: ~803s @2565 test, 2026-07-09.)
CAP=1200
SUITE='quick'          # quick | slow
ACTION='run'           # run | print | quarantine
QCAP=60                # per-test cap in --quarantine

while [[ $# -gt 0 ]]; do
    case "$1" in
        --filter)        POSITIVE_FILTER="$2"; shift 2 ;;
        --slow)          SUITE='slow'; CAP=1800; shift ;;
        --quarantine)    ACTION='quarantine'; shift ;;
        --print-filter|--dry-run) ACTION='print'; shift ;;
        --print-ctest-regex) ACTION='print-ctest-regex'; shift ;;
        --cap)           CAP="$2"; shift 2 ;;
        --qcap)          QCAP="$2"; shift 2 ;;
        -h|--help)       sed -n '1,45p' "$0"; exit 0 ;;
        *) echo "[test_quick] unknown arg: $1" >&2; exit 2 ;;
    esac
done

# ── Banner quarantena (onestà: il verde è "verde ECCETTO questi rossi noti") ──
print_quarantine_banner() {
    local n="${#QUARANTINE[@]}"
    if (( n > 0 )); then
        echo "[test_quick] ⚠ $n test in QUARANTENA esclusi (NON coperti da questo run):"
        printf '            - %s\n' "${QUARANTINE[@]}"
        echo "            → noti-rossi tracciati in scripts/test_quarantine.txt (+ ledger)."
        echo "            → stato reale: bash scripts/test_quick.sh --quarantine"
    fi
}

# ── Azione quarantine: ogni test isolato con timeout stretto ──────────────────
if [[ "$ACTION" == "quarantine" ]]; then
    [[ ! -x "$BIN" ]] && { echo "[test_quick] binario non buildato: $BIN" >&2; exit 2; }
    echo "[test_quick] --quarantine: ${#QUARANTINE[@]} test, per-test cap ${QCAP}s"
    pass=0; fail=0; tmo=0
    for name in ${QUARANTINE[@]+"${QUARANTINE[@]}"}; do
        if timeout "$QCAP" "$BIN" --gtest_filter="$name" >/dev/null 2>&1; then
            echo "  PASS    $name   (potenziale candidato a de-quarantena!)"
            pass=$((pass+1))
        else
            rc=$?
            if [[ $rc -eq 124 ]]; then echo "  TIMEOUT $name   (hang > ${QCAP}s)"; tmo=$((tmo+1))
            else echo "  FAIL    $name   (exit $rc)"; fail=$((fail+1)); fi
        fi
    done
    echo "[test_quick] quarantena: pass=$pass fail=$fail timeout=$tmo"
    # I PASS qui sono test pronti a uscire dalla quarantena (abbassa CEILING).
    exit 0
fi

# ── Costruzione filtro quick / slow ──────────────────────────────────────────
if [[ "$SUITE" == "slow" ]]; then
    # ${arr[@]+...}: bash 3.2 + set -u tratta l'array vuoto come unbound
    NEG=( "${STRUCTURAL[@]}" ${QUARANTINE[@]+"${QUARANTINE[@]}"} )   # SLOW_OK INCLUSI
else
    NEG=( "${STRUCTURAL[@]}" "${SLOW_OK[@]}" ${QUARANTINE[@]+"${QUARANTINE[@]}"} )
fi
EXCLUDE="-$(join_colon "${NEG[@]}")"

if [[ -n "$POSITIVE_FILTER" ]]; then
    FILTER="${POSITIVE_FILTER}${EXCLUDE}"
else
    FILTER="$EXCLUDE"
fi

if [[ "$ACTION" == "print" ]]; then
    print_quarantine_banner
    echo "[test_quick] suite=${SUITE} action=print cap=${CAP}s"
    echo "[test_quick] filter='${FILTER}'"
    exit 0
fi

# gtest_discover_tests registra ogni case come test ctest a sé (nome
# "Suite.Case"): stessa NEG-list del filtro binario, tradotta in ERE per
# `ctest -E`. Un solo elenco sorgente — CI (che gira ctest, non il binario
# diretto) ed esecuzione locale non possono divergere.
if [[ "$ACTION" == "print-ctest-regex" ]]; then
    parts=()
    for pat in "${NEG[@]}"; do
        escaped="${pat//./\\.}"
        escaped="${escaped//\*/.*}"
        parts+=("^${escaped}\$")
    done
    IFS='|'
    echo "${parts[*]}"
    exit 0
fi

# ── Esecuzione ───────────────────────────────────────────────────────────────
if [[ ! -x "$BIN" ]]; then
    echo "[test_quick] cas_foundation_tests non buildato. Esegui: cmake --build build --target cas_foundation_tests" >&2
    exit 2
fi

print_quarantine_banner
echo "[test_quick] suite=${SUITE} cap=${CAP}s"
echo "[test_quick] filter='${FILTER}'"
exec timeout "${CAP}" "$BIN" --gtest_filter="$FILTER"
