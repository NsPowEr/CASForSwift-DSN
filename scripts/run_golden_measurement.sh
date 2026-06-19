#!/usr/bin/env bash
# run_golden_measurement.sh — Orchestrator: Maxima refs + golden_runner + merge + ratchet.
#
# Usage:
#   bash scripts/run_golden_measurement.sh                    # full measurement
#   bash scripts/run_golden_measurement.sh --update-baseline  # + lock baseline
#   bash scripts/run_golden_measurement.sh --skip-maxima      # reuse existing refs
#   bash scripts/run_golden_measurement.sh --area simplify    # single area only
#
# Output: build-golden/report.json  (merged, consumable by check_golden_ratchet.sh)
#
# Prerequisiti:
#   - build/cas_golden_runner buildato (cmake --build build --target cas_golden_runner)
#   - maxima 5.49.0 nel PATH (/opt/homebrew/bin/maxima)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

GOLDEN_BIN="build/cas_golden_runner"
CORPUS_DIR="test/golden/corpus"
OUT_DIR="build-golden"
MAXIMA_OUT="${OUT_DIR}/maxima_out"
REPORT="${OUT_DIR}/report.json"
PER_ENTRY_TIMEOUT=30

SKIP_MAXIMA=0
UPDATE_BASELINE=0
SINGLE_AREA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-maxima)      SKIP_MAXIMA=1; shift ;;
        --update-baseline)  UPDATE_BASELINE=1; shift ;;
        --area)             SINGLE_AREA="$2"; shift 2 ;;
        --per-entry-timeout) PER_ENTRY_TIMEOUT="$2"; shift 2 ;;
        -h|--help)
            sed -n '1,14p' "$0"
            exit 0
            ;;
        *) echo "arg sconosciuto: $1" >&2; exit 2 ;;
    esac
done

if [[ ! -x "$GOLDEN_BIN" ]]; then
    echo "ERROR: $GOLDEN_BIN non trovato. Builda con:" >&2
    echo "  cmake --build build --target cas_golden_runner" >&2
    exit 1
fi

AREAS=()
if [[ -n "$SINGLE_AREA" ]]; then
    AREAS=("$SINGLE_AREA")
else
    for d in "$CORPUS_DIR"/*/; do
        area=$(basename "$d")
        AREAS+=("$area")
    done
fi

corpus_file_for() {
    local area="$1"
    if [[ "$area" == "bronstein" ]]; then
        echo "${CORPUS_DIR}/bronstein/integrals.jsonl"
    else
        echo "${CORPUS_DIR}/${area}/basic.jsonl"
    fi
}

mkdir -p "$OUT_DIR"

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Golden Measurement — CAS vs Maxima 5.49.0         ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Aree:    ${AREAS[*]}"
echo "Output:  $REPORT"
echo ""

total_pass=0
total_fail=0
total_skip=0

for area in "${AREAS[@]}"; do
    corpus=$(corpus_file_for "$area")
    if [[ ! -f "$corpus" ]]; then
        echo "WARN: corpus assente per '$area': $corpus — skip" >&2
        continue
    fi

    area_maxima_dir="${MAXIMA_OUT}/${area}"
    area_json="${OUT_DIR}/golden_${area}.json"
    lines=$(grep -c '{' "$corpus" 2>/dev/null || echo "0")

    echo "────────────────────────────────────────────────────────────────"
    echo "  Area: $area ($lines entries)"
    echo "────────────────────────────────────────────────────────────────"

    # Step 1: Maxima reference generation
    if [[ "$SKIP_MAXIMA" == "0" ]]; then
        echo "  [1/2] Maxima refs → ${area_maxima_dir}/"
        mkdir -p "$area_maxima_dir"
        if ! timeout 600 bash scripts/run_golden_maxima.sh "$corpus" "$area_maxima_dir" 2>&1 | \
             grep -E 'Done\.|ERROR|TIMEOUT|Total'; then
            echo "  WARN: run_golden_maxima.sh returned non-zero for $area"
        fi
    else
        echo "  [1/2] Maxima refs → SKIP (--skip-maxima, riuso esistenti)"
        if [[ ! -d "$area_maxima_dir" ]]; then
            echo "  WARN: $area_maxima_dir non esiste — area skippata" >&2
            continue
        fi
    fi

    # Step 2: golden_runner comparison
    echo "  [2/2] golden_runner → ${area_json}"
    if ! timeout 300 "$GOLDEN_BIN" "$corpus" "$area_maxima_dir" \
         --json "$area_json" --per-entry-timeout "$PER_ENTRY_TIMEOUT" 2>&1 | \
         grep -E '^=|^-|^TOTAL|Area|PASS%'; then
        echo "  WARN: golden_runner returned non-zero for $area"
    fi

    # Extract pass/fail from the per-area JSON
    if [[ -f "$area_json" ]]; then
        read -r ap af as <<< "$(python3 -c "
import json, sys
d = json.load(open('$area_json'))
p = sum(v.get('pass',0) for v in d.values() if isinstance(v, dict))
f = sum(v.get('fail',0) for v in d.values() if isinstance(v, dict))
s = sum(v.get('skip',0) for v in d.values() if isinstance(v, dict))
print(p, f, s)
" 2>/dev/null || echo "0 0 0")"
        total_pass=$((total_pass + ap))
        total_fail=$((total_fail + af))
        total_skip=$((total_skip + as))
        echo "  → $area: pass=$ap fail=$af skip=$as"
    fi
    echo ""
done

# Step 3: Merge per-area JSONs into one report.json
echo "────────────────────────────────────────────────────────────────"
echo "  Merging → $REPORT"
echo "────────────────────────────────────────────────────────────────"

python3 - "$OUT_DIR" "$REPORT" <<'PY'
import json, sys, os, glob

out_dir = sys.argv[1]
report_path = sys.argv[2]
merged = {}

for f in sorted(glob.glob(os.path.join(out_dir, "golden_*.json"))):
    try:
        d = json.load(open(f))
        for area, stats in d.items():
            if isinstance(stats, dict):
                if area in merged:
                    for k in ("pass", "fail", "skip"):
                        merged[area][k] = merged[area].get(k, 0) + stats.get(k, 0)
                    ex = merged[area].get("examples_fail", [])
                    ex.extend(stats.get("examples_fail", []))
                    merged[area]["examples_fail"] = ex[:10]
                else:
                    merged[area] = dict(stats)
    except Exception as e:
        print(f"WARN: skip {f}: {e}", file=sys.stderr)

with open(report_path, "w") as fp:
    json.dump(merged, fp, indent=2)

total_p = sum(v.get("pass", 0) for v in merged.values())
total_f = sum(v.get("fail", 0) for v in merged.values())
total_s = sum(v.get("skip", 0) for v in merged.values())
rate = 100 * total_p / (total_p + total_f) if (total_p + total_f) > 0 else 0
print(f"  Merged: {len(merged)} areas, {total_p} pass / {total_f} fail / {total_s} skip ({rate:.1f}%)")
PY

echo ""

# Step 4: Summary
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                    MEASUREMENT SUMMARY                     ║"
echo "╚══════════════════════════════════════════════════════════════╝"
if (( total_pass + total_fail > 0 )); then
    rate=$(awk "BEGIN{printf \"%.1f\", 100*$total_pass/($total_pass+$total_fail)}")
else
    rate="n/a"
fi
echo "  PASS:  $total_pass"
echo "  FAIL:  $total_fail"
echo "  SKIP:  $total_skip"
echo "  RATE:  ${rate}%"
echo ""

# Step 5: Ratchet check
echo "  Running ratchet check..."
bash scripts/check_golden_ratchet.sh --report "$REPORT"
echo ""

# Step 6: Optional baseline update
if [[ "$UPDATE_BASELINE" == "1" ]]; then
    echo "  Updating baseline..."
    bash scripts/check_golden_ratchet.sh --report "$REPORT" --update-baseline
fi

echo "Done."
