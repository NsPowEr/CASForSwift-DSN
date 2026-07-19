#!/usr/bin/env bash
# run_golden_giac.sh — run a JSONL corpus file through Giac 2.0.0 (icas) and
# save one .giac.out file per input.
#
# Usage: bash scripts/run_golden_giac.sh <corpus.jsonl> <output_dir> [--per-entry-timeout N]
#
# Each line in <corpus.jsonl> must be a JSON object with at least:
#   {"input": "<expr>", "area": "<area>", "ref": "<ref>"}
#
# For each input line the script:
#   1. Converts the CAS-style input to Giac syntax (see translate_to_giac()).
#   2. Invokes icas '<expr>' as a separate process (argv batch mode).
#   3. Writes the cleaned text output to <output_dir>/<line_index>.giac.out
#      (first line = verdict tag: ANSWERED | UNEVALUATED | TIMEOUT | ERROR,
#       second line onward = raw giac result).
#   4. Writes a manifest JSON to <output_dir>/manifest.json summarising results.
#
# GPL-3.0-or-later compliance (CLAUDE.md §6): giac is invoked ONLY as an
# external subprocess (fork/exec) and its textual output parsed. Its sources
# are never consulted, modified, or embedded.

set -euo pipefail

GIAC_BIN="${GIAC_BIN:-/opt/homebrew/bin/icas}"
PER_ENTRY_TIMEOUT="${GIAC_PER_ENTRY_TIMEOUT:-30}"
CORPUS_FILE=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 <corpus.jsonl> <output_dir> [--per-entry-timeout N]" >&2
    exit 1
}

if [[ $# -lt 2 ]]; then usage; fi
CORPUS_FILE="$1"
OUTPUT_DIR="$2"
shift 2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --per-entry-timeout)
            PER_ENTRY_TIMEOUT="${2:?--per-entry-timeout requires a value}"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            ;;
    esac
done

if [[ ! -f "$CORPUS_FILE" ]]; then
    echo "ERROR: corpus file not found: $CORPUS_FILE" >&2
    exit 1
fi

if [[ ! -x "$GIAC_BIN" ]]; then
    echo "ERROR: giac binary not found or not executable: $GIAC_BIN" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# translate_to_giac: convert our CAS syntax to Giac syntax.
# Giac's surface syntax is very close to ours; only the deltas are mapped.
# Anything giac cannot parse will surface as ERROR and be classified as such
# (never silently dropped).
# ---------------------------------------------------------------------------
translate_to_giac() {
    local out="$1"

    # limit(f, x, a, plus)  -> limit(f, x, a, 1)
    # limit(f, x, a, minus) -> limit(f, x, a, -1)
    out=$(echo "$out" | sed 's/, *plus)/,1)/g; s/, *minus)/,-1)/g')

    # series(f, x, a, n) -> series(f, x=a, n)   (giac wants x=a form)
    out=$(echo "$out" | sed -E 's/series\(([^,]+), *([A-Za-z_][A-Za-z0-9_]*) *, *([^,]+), *([^)]+)\)/series(\1,\2=\3,\4)/g')

    # inverse( -> inv(   (giac matrix inverse)
    out=$(echo "$out" | sed 's/inverse(/inv(/g')

    # eigenvalues( -> eigenvals(
    out=$(echo "$out" | sed 's/eigenvalues(/eigenvals(/g')

    # log( in our corpus is natural log; giac's log() is natural log too — same.
    # pi, e, i, inf, [[..],[..]] matrix literals, det(, trace(, transpose(,
    # rank(, factorial(, binomial(, gamma(, erf(, sqrt(, abs( — identical.

    echo "$out"
}

# ---------------------------------------------------------------------------
# classify: derive the verdict tag from giac's cleaned output.
#   $1 = cleaned result text, $2 = translated input expression
# ---------------------------------------------------------------------------
classify() {
    local result="$1"
    local expr="$2"

    if [[ -z "$result" ]]; then
        echo "ERROR"
        return
    fi
    if echo "$result" | grep -qiE 'error|not implemented|unable to'; then
        echo "ERROR"
        return
    fi
    # Unevaluated echo: giac returns the same top-level call it was given
    # (integrate(...), limit(...), sum(...)) when it cannot solve it.
    local head_fn
    head_fn=$(echo "$expr" | grep -oE '^[A-Za-z_][A-Za-z0-9_]*' || true)
    if [[ -n "$head_fn" ]] && echo "$result" | grep -qE "(^|[^A-Za-z0-9_])(${head_fn}|int|integration|limit|sum|product)\(" ; then
        echo "UNEVALUATED"
        return
    fi
    echo "ANSWERED"
}

# ---------------------------------------------------------------------------
# run_one: run a single input through giac and write the output file.
# ---------------------------------------------------------------------------
run_one() {
    local idx="$1"
    local input_expr="$2"
    local out_file="${OUTPUT_DIR}/${idx}.giac.out"

    local giac_expr
    giac_expr=$(translate_to_giac "$input_expr")

    local raw verdict result
    if raw=$(timeout "${PER_ENTRY_TIMEOUT}s" "$GIAC_BIN" "$giac_expr" 2>/dev/null); then
        # Strip banner/comment noise; keep the final expression line(s).
        result=$(echo "$raw" \
            | grep -vE '^//|^Welcome|^Added [0-9]+ synonyms|^Unable to (find keyword|open HTML)|^\*\*\*|^Press CTRL|^Type \?|^See http|^Released|^May contain|^Homepage|^\(c\)|^-+$' \
            | sed '/^[[:space:]]*$/d' \
            | tail -5)
        verdict=$(classify "$result" "$giac_expr")
        { echo "$verdict"; echo "$result"; } > "$out_file"
    else
        local exit_code=$?
        if [[ $exit_code -eq 124 ]]; then
            echo "TIMEOUT" > "$out_file"
        else
            echo "ERROR" > "$out_file"
            echo "exit=$exit_code" >> "$out_file"
        fi
    fi
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
total=0
idx=0
declare -a INPUTS=()
while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    input=$(echo "$line" | python3 -c "import json,sys; print(json.load(sys.stdin).get('input',''))" 2>/dev/null || echo "")
    [[ -z "$input" ]] && continue
    INPUTS[$idx]="$input"
    idx=$((idx + 1))
done < "$CORPUS_FILE"
total=$idx

echo "run_golden_giac: ${total} entries from ${CORPUS_FILE} -> ${OUTPUT_DIR} (timeout ${PER_ENTRY_TIMEOUT}s/entry)"

answered=0; uneval=0; timeouts=0; errors=0
for ((i = 0; i < total; i++)); do
    run_one "$i" "${INPUTS[$i]}"
    v=$(head -1 "${OUTPUT_DIR}/${i}.giac.out")
    case "$v" in
        ANSWERED)    answered=$((answered + 1)) ;;
        UNEVALUATED) uneval=$((uneval + 1)) ;;
        TIMEOUT)     timeouts=$((timeouts + 1)) ;;
        *)           errors=$((errors + 1)) ;;
    esac
done

# Manifest
python3 - "$OUTPUT_DIR" "$CORPUS_FILE" "$total" "$answered" "$uneval" "$timeouts" "$errors" <<'PY'
import json, sys, datetime
out_dir, corpus, total, answered, uneval, timeouts, errors = sys.argv[1:8]
manifest = {
    "corpus": corpus,
    "oracle": "giac-2.0.0 (icas, fork/exec)",
    "generated": datetime.datetime.now().isoformat(timespec="seconds"),
    "total": int(total),
    "answered": int(answered),
    "unevaluated": int(uneval),
    "timeout": int(timeouts),
    "error": int(errors),
}
with open(out_dir + "/manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)
PY

echo "Done. answered=${answered} unevaluated=${uneval} timeout=${timeouts} error=${errors} (total ${total})"
