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

STALE_AREAS=()
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

    # Step 2: golden_runner comparison.
    # The full per-entry stdout (PASS/FAIL/SKIP + skip reason) is kept in
    # log_<area>.txt: the JSON only carries aggregates, so that log is the only
    # source for per-entry skip triage (scripts/golden_skip_triage.py).
    area_log="${OUT_DIR}/log_${area}.txt"
    # A37: the area budget must dominate the worst case (every entry burning its
    # own per-entry timeout), otherwise `timeout` kills the runner before it
    # writes the JSON and the area silently keeps its PREVIOUS measurement — the
    # merged report then mixes fresh and months-old numbers (that is exactly how
    # bronstein carried a 35-day-old result into the ratchet). A fixed 300s cap
    # was below the worst case of every area with >10 entries.
    area_budget=$(( lines * PER_ENTRY_TIMEOUT + 120 ))
    # Move the previous measurement aside: the run is fresh iff it recreates
    # the file. (A timestamp comparison is not usable here — `-nt` has
    # one-second granularity and reports a same-second rewrite as stale.)
    prev_json="${area_json}.prev"
    rm -f "$prev_json"
    [[ -f "$area_json" ]] && mv "$area_json" "$prev_json"
    echo "  [2/2] golden_runner → ${area_json}  (log: ${area_log}, budget ${area_budget}s)"
    if ! timeout "$area_budget" "$GOLDEN_BIN" "$corpus" "$area_maxima_dir" \
         --json "$area_json" --per-entry-timeout "$PER_ENTRY_TIMEOUT" 2>&1 | \
         tee "$area_log" | grep -E '^=|^-|^TOTAL|Area|PASS%'; then
        echo "  WARN: golden_runner returned non-zero for $area"
    fi

    # Staleness guard: no new JSON means the runner was killed or crashed.
    # Reporting the previous numbers would be a measurement lie, so restore
    # them for inspection but refuse to aggregate.
    if [[ ! -f "$area_json" ]]; then
        [[ -f "$prev_json" ]] && mv "$prev_json" "$area_json"
        echo "  ERROR: $area_json non rigenerato (runner ucciso o crashato):" >&2
        echo "         il dato precedente NON viene aggregato. Log: $area_log" >&2
        STALE_AREAS+=("$area")
        echo ""
        continue
    fi
    rm -f "$prev_json"

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

# Refuse to merge a mixed-freshness report: the merge globs every
# golden_*.json, so a stale one would be aggregated as if it were current and
# the ratchet would gate on numbers nobody measured (A37).
if [[ ${#STALE_AREAS[@]} -gt 0 ]]; then
    echo "────────────────────────────────────────────────────────────────" >&2
    echo "ERRORE: aree non rigenerate: ${STALE_AREAS[*]}" >&2
    echo "Il merge e il ratchet sono ABORTITI: report.json terrebbe dati stantii." >&2
    echo "Rilancia le aree fallite (--area <a> --skip-maxima) e controlla il log." >&2
    exit 3
fi

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
