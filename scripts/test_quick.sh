#!/usr/bin/env bash
# test_quick.sh — esegue la suite di test "quick" con esclusione canonica
# dei test pesanti (>30s) e protezione anti-hang via timeout.
#
# Uso:
#   bash scripts/test_quick.sh              # default 600s wall-clock cap (quick suite)
#   bash scripts/test_quick.sh --filter X   # filtro positivo aggiuntivo
#   bash scripts/test_quick.sh --slow       # include anche i test slow (cap 1800s)
#
# Per CLAUDE.md "REGOLA TIMEOUT TEST": questa è l'unica via raccomandata per
# eseguire la suite durante lo sviluppo iterativo. La suite completa con tutti
# i test slow va eseguita solo come gate pre-commit (`--slow`).
#
# Lista di esclusione (test slow noti, pre-esistenti, >30s individualmente):
#   - *Stress*              — esclusi da policy CLAUDE.md (cas-regression-guard).
#   - *Fuzz*                — esclusi: invocati separatamente via libFuzzer.
#   - *Disabled*            — gtest DISABLED_ convention.
#   - VanHoeijDirect.Deg16_EightQuadratics_FindsRealFactor (181s legittimo).
#   - FactorizationTowerNTest.SplitsX2Minus5_Over_Q_Sqrt2_Sqrt3_Sqrt5 (>240s).
#   - PrimitiveElementTest.RedundantMixedTower_Sqrt2_Sqrt3_Sqrt5_Sqrt6 (>60s).
#   - PrimitiveElementTest.DetectTowerNLevel_SqrtTwoSqrtThreeSqrtFive (31s).
#
# I tempi sopra riflettono baseline su macOS arm64. Se introduci una regressione
# di complessità che porta un test sopra il cap, NON aggiungerlo qui: indaga la
# causa (probabilmente bug O(2^n) accidentale) e fixa il codice.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${REPO_ROOT}/build/cas_foundation_tests"
if [[ ! -x "$BIN" ]]; then
    echo "[test_quick] cas_foundation_tests not built. Run: cmake --build build --target cas_foundation_tests" >&2
    exit 2
fi

# Canonical exclusion patterns (single '-' prefix; ':' separates negative globs).
EXCLUDE='-*Stress*:*Fuzz*:*Disabled*'
EXCLUDE+=':VanHoeijDirect.Deg16_EightQuadratics_FindsRealFactor'
EXCLUDE+=':FactorizationTowerNTest.SplitsX2Minus5_Over_Q_Sqrt2_Sqrt3_Sqrt5'
EXCLUDE+=':PrimitiveElementTest.RedundantMixedTower_Sqrt2_Sqrt3_Sqrt5_Sqrt6'
EXCLUDE+=':PrimitiveElementTest.DetectTowerNLevel_SqrtTwoSqrtThreeSqrtFive'
EXCLUDE+=':VanHoeijFactorTest.AcceptanceGate_AG2_SwinnertonDyer_SD3_Irreducible'

POSITIVE_FILTER=''
CAP=600
INCLUDE_SLOW=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --filter)
            POSITIVE_FILTER="$2"
            shift 2
            ;;
        --slow)
            INCLUDE_SLOW=1
            CAP=1800
            EXCLUDE='-*Stress*:*Fuzz*:*Disabled*'
            shift
            ;;
        --cap)
            CAP="$2"
            shift 2
            ;;
        -h|--help)
            sed -n '1,40p' "$0"
            exit 0
            ;;
        *)
            echo "[test_quick] unknown arg: $1" >&2
            exit 2
            ;;
    esac
done

if [[ -n "$POSITIVE_FILTER" ]]; then
    FILTER="${POSITIVE_FILTER}${EXCLUDE}"
else
    FILTER="$EXCLUDE"
fi

echo "[test_quick] cap=${CAP}s slow=${INCLUDE_SLOW}"
echo "[test_quick] filter='${FILTER}'"
exec timeout "${CAP}" "$BIN" --gtest_filter="$FILTER"
