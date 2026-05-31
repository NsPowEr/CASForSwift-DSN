#!/usr/bin/env bash
# run_golden_maxima.sh — run a JSONL corpus file through Maxima 5.49.0 and
# save one .maxima.out file per input.
#
# Usage: bash scripts/run_golden_maxima.sh <corpus.jsonl> <output_dir> [--parallel N]
#
# Each line in <corpus.jsonl> must be a JSON object with at least:
#   {"input": "<expr>", "area": "<area>", "ref": "<ref>"}
#
# For each input line the script:
#   1. Converts the CAS-style input to a Maxima batch string (see translate_to_maxima()).
#   2. Invokes maxima --very-quiet --batch-string="..."
#   3. Writes the raw text output to <output_dir>/<line_index>.maxima.out
#   4. Writes a manifest JSON to <output_dir>/manifest.json summarising results.
#
# The script does NOT modify Maxima in any way (GPL-2.0-only compliance).
# Maxima is invoked only as an external subprocess.

set -euo pipefail

MAXIMA_BIN="${MAXIMA_BIN:-/opt/homebrew/bin/maxima}"
PARALLEL_N=1
CORPUS_FILE=""
OUTPUT_DIR=""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
    echo "Usage: $0 <corpus.jsonl> <output_dir> [--parallel N]" >&2
    exit 1
}

if [[ $# -lt 2 ]]; then usage; fi
CORPUS_FILE="$1"
OUTPUT_DIR="$2"
shift 2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            PARALLEL_N="${2:-4}"
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

if [[ ! -x "$MAXIMA_BIN" ]]; then
    echo "ERROR: maxima binary not found or not executable: $MAXIMA_BIN" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ---------------------------------------------------------------------------
# translate_to_maxima: convert our CAS syntax to Maxima batch syntax.
# Rules documented here cover the differences between our parser and Maxima.
# ---------------------------------------------------------------------------
translate_to_maxima() {
    local expr="$1"

    # Wrap top-level calls that are not already Maxima functions.
    # Our syntax: integrate(f, x)   -> Maxima: integrate(f, x)  (same)
    # Our syntax: diff(f, x)        -> Maxima: diff(f, x)        (same)
    # Our syntax: diff(f, x, n)     -> Maxima: diff(f, x, n)     (same)
    # Our syntax: limit(f, x, a)    -> Maxima: limit(f, x, a)    (same)
    # Our syntax: limit(f, x, a, plus/minus) -> Maxima: limit(f,x,a,plus/minus)
    # Our syntax: solve(f, x)       -> Maxima: solve(f, x) but needs [x]
    # Our syntax: factor(p)         -> Maxima: factor(p)          (same)
    # Our syntax: gcd(p, q)         -> Maxima: gcd(p, q)          (same, returns poly)
    # Our syntax: series(f,x,a,n)   -> Maxima: taylor(f, x, a, n)
    # Our syntax: det([[...]])       -> Maxima: determinant(matrix(...))
    # Our syntax: trace([[...]])     -> Maxima: mattrace(matrix(...))
    # Our syntax: transpose([[...]]) -> Maxima: transpose(matrix(...))
    # Our syntax: inverse([[...]])   -> Maxima: invert(matrix(...))
    # Our syntax: rank([[...]])      -> Maxima: rank(matrix(...))
    # Our syntax: eigenvalues([[...]])-> Maxima: eigenvalues(matrix(...))
    # Our syntax: pi                 -> Maxima: %pi
    # Our syntax: e (as constant)    -> Maxima: %e
    # Our syntax: i (imaginary unit) -> Maxima: %i
    # Our syntax: inf                -> Maxima: inf
    # Our syntax: sqrt(x)            -> Maxima: sqrt(x)           (same)
    # Our syntax: log(x)             -> Maxima: log(x)            (same, natural log)
    # Our syntax: abs(x)             -> Maxima: abs(x)            (same)
    # Our syntax: atan(x)            -> Maxima: atan(x)           (same)
    # Our syntax: asin(x)            -> Maxima: asin(x)           (same)
    # Our syntax: acos(x)            -> Maxima: acos(x)           (same)
    # Our syntax: gamma(x)           -> Maxima: gamma(x)          (same)
    # Our syntax: factorial(n)       -> Maxima: n!                (or factorial(n))
    # Our syntax: binomial(n,k)      -> Maxima: binomial(n, k)    (same)
    # Our syntax: zeta(s)            -> Maxima: zeta(s)           (same)
    # Our syntax: erf(x)             -> Maxima: erf(x)            (same)
    # Our syntax: BesselJ(n, x)      -> Maxima: bessel_j(n, x)
    # Our syntax: LegendreP(n, x)    -> Maxima: legendre_p(n, x)
    # Our syntax: sinh(x)/cosh(x)    -> Maxima: same
    # Our syntax: LambertW(x)        -> Maxima: lambert_w(x)
    # Our syntax: LaguerreL(n,x)     -> Maxima: laguerre(x, n) (note arg order!)

    local out="$expr"

    # series -> taylor  (BSD sed on macOS: use [[:<:]] word boundary or just literal match)
    out=$(echo "$out" | sed 's/series(/taylor(/g')

    # BesselJ -> bessel_j
    out=$(echo "$out" | sed 's/BesselJ(/bessel_j(/g')

    # LegendreP -> legendre_p
    out=$(echo "$out" | sed 's/LegendreP(/legendre_p(/g')

    # LambertW -> lambert_w
    out=$(echo "$out" | sed 's/LambertW(/lambert_w(/g')

    # LaguerreL(n, x) -> laguerre(x, n)  - handled specially, skip sed for now
    # (The runner will mark these SKIP if they fail parse)

    # matrix literal [[a,b],[c,d]] -> matrix([a,b],[c,d])
    # This sed handles one level of nesting
    out=$(echo "$out" | sed 's/\[\[/matrix([/g; s/\],\[/],[/g; s/\]\]/])/g')

    # det( -> determinant(
    out=$(echo "$out" | sed 's/det(/determinant(/g')

    # trace( -> mattrace(
    out=$(echo "$out" | sed 's/trace(/mattrace(/g')

    # inverse( -> invert(
    out=$(echo "$out" | sed 's/inverse(/invert(/g')

    # eigenvalues( stays eigenvalues in Maxima

    # rank( stays rank in Maxima

    # solve(f, x) -> solve([f], [x])
    # BSD sed: use basic regex without \b
    out=$(echo "$out" | sed 's/solve(\([^,]*\), *\([^)]*\))/solve([\1], [\2])/g')

    # gcd in Maxima uses gcd(poly, poly, var) for polynomial gcd,
    # but the simple gcd(a,b) works for integers.
    # Our corpus uses gcd(expr, expr) — pass through as-is; Maxima's gcd handles it.

    # pi -> %pi  (BSD sed: match standalone 'pi' using bracket workaround)
    # Replace pi that is not preceded/followed by alphanumeric
    out=$(echo "$out" | sed 's/pi/%pi/g')

    # inf -> inf (already correct in Maxima)
    # plus/minus direction in limit: already correct keyword strings

    echo "$out"
}

# ---------------------------------------------------------------------------
# run_one: run a single input through Maxima and write output file.
# Args: $1=index(0-based)  $2=input_expr  $3=output_dir
# ---------------------------------------------------------------------------
run_one() {
    local idx="$1"
    local input_expr="$2"
    local out_dir="$3"
    local out_file="${out_dir}/${idx}.maxima.out"

    local maxima_expr
    maxima_expr=$(translate_to_maxima "$input_expr")

    # Build batch string: display2d:false to get 1-line output, then the expr.
    local batch="display2d:false$ ${maxima_expr};"

    # Run maxima; capture stdout+stderr; timeout 30s per input.
    local result
    if result=$(timeout 30s "$MAXIMA_BIN" --very-quiet --batch-string="$batch" 2>&1); then
        # Strip the input echo line that Maxima sometimes emits
        # Keep only lines that start with (%o or are a pure expression
        echo "$result" > "$out_file"
    else
        local exit_code=$?
        if [[ $exit_code -eq 124 ]]; then
            echo "TIMEOUT" > "$out_file"
        else
            echo "ERROR:exit=$exit_code" > "$out_file"
            echo "$result" >> "$out_file"
        fi
    fi
}

export -f run_one translate_to_maxima
export MAXIMA_BIN

# ---------------------------------------------------------------------------
# Main: read corpus line by line, dispatch to run_one
# ---------------------------------------------------------------------------
MANIFEST_FILE="${OUTPUT_DIR}/manifest.json"
TOTAL=0
ERRORS=0

echo "Running corpus: $CORPUS_FILE"
echo "Output dir:     $OUTPUT_DIR"
echo "Maxima binary:  $MAXIMA_BIN"
echo "Parallelism:    $PARALLEL_N"
echo ""

# Collect jobs for optional parallel execution
declare -a JOB_INDICES=()
declare -a JOB_INPUTS=()

IDX=0
while IFS= read -r line; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" == \#* ]] && continue

    # Extract "input" field using basic parameter expansion (no jq dependency)
    input_expr=$(echo "$line" | python3 -c "import sys,json; d=json.loads(sys.stdin.read()); print(d.get('input',''))" 2>/dev/null || true)
    if [[ -z "$input_expr" ]]; then
        echo "WARN: could not parse line $IDX: $line" >&2
        IDX=$((IDX + 1))
        continue
    fi

    JOB_INDICES+=("$IDX")
    JOB_INPUTS+=("$input_expr")
    IDX=$((IDX + 1))
    TOTAL=$IDX
done < "$CORPUS_FILE"

echo "Total inputs: $TOTAL"
echo ""

# Run jobs (sequential or parallel)
if [[ "$PARALLEL_N" -le 1 ]]; then
    for i in "${!JOB_INDICES[@]}"; do
        idx="${JOB_INDICES[$i]}"
        inp="${JOB_INPUTS[$i]}"
        printf "  [%3d/%3d] %s ... " "$((idx+1))" "$TOTAL" "$inp"
        run_one "$idx" "$inp" "$OUTPUT_DIR"
        out_file="${OUTPUT_DIR}/${idx}.maxima.out"
        if grep -q "^TIMEOUT\|^ERROR:" "$out_file" 2>/dev/null; then
            status=$(head -1 "$out_file")
            echo "$status"
            ERRORS=$((ERRORS + 1))
        else
            echo "OK"
        fi
    done
else
    # GNU parallel or xargs -P fallback
    if command -v parallel &>/dev/null; then
        printf '%s\n' "${!JOB_INDICES[@]}" | parallel -j"$PARALLEL_N" \
            run_one "${JOB_INDICES[{1}]}" "${JOB_INPUTS[{1}]}" "$OUTPUT_DIR"
    else
        # xargs fallback: write index:input pairs to temp file
        TMPF=$(mktemp)
        for i in "${!JOB_INDICES[@]}"; do
            echo "${JOB_INDICES[$i]}|${JOB_INPUTS[$i]}"
        done > "$TMPF"
        cat "$TMPF" | xargs -P"$PARALLEL_N" -I{} bash -c '
            IFS="|" read -r idx inp <<< "{}"
            run_one "$idx" "$inp" "'"$OUTPUT_DIR"'"
        '
        rm -f "$TMPF"
    fi
fi

# ---------------------------------------------------------------------------
# Write manifest
# ---------------------------------------------------------------------------
{
    echo "{"
    echo "  \"corpus\": \"$CORPUS_FILE\","
    echo "  \"total\": $TOTAL,"
    echo "  \"errors\": $ERRORS,"
    echo "  \"maxima_bin\": \"$MAXIMA_BIN\","
    echo "  \"timestamp\": \"$(date -u +%Y-%m-%dT%H:%M:%SZ)\""
    echo "}"
} > "$MANIFEST_FILE"

echo ""
echo "Done. $TOTAL inputs processed, $ERRORS errors."
echo "Manifest: $MANIFEST_FILE"
