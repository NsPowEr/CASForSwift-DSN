#!/usr/bin/env bash
# F0.6 — Benchmark regression gate.
#
# Compares a current benchmark output file against a baseline.
# Exits 1 if any metric regresses by more than REGRESSION_THRESHOLD_PCT.
#
# Usage:
#   scripts/check_bench_regression.sh \
#     --current  <path-to-current.txt> \
#     --baseline <path-to-baseline.txt> \
#     [--threshold 10] [--verbose]
#
# File format (one metric per line):
#   <name> <value_ms>
#   parse_implicit_mul_ms 0.009

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

CURRENT_FILE=""
BASELINE_FILE="$repo_root/test/benchmarks/baseline_release.txt"
THRESHOLD_PCT=10
VERBOSE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --current)   CURRENT_FILE="$2";  shift 2 ;;
        --baseline)  BASELINE_FILE="$2"; shift 2 ;;
        --threshold) THRESHOLD_PCT="$2"; shift 2 ;;
        --verbose)   VERBOSE=1;           shift   ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$CURRENT_FILE" ]]; then
    echo "Error: --current <file> is required." >&2
    exit 1
fi
if [[ ! -f "$CURRENT_FILE" ]]; then
    echo "Error: current file not found: $CURRENT_FILE" >&2
    exit 1
fi
if [[ ! -f "$BASELINE_FILE" ]]; then
    echo "Error: baseline file not found: $BASELINE_FILE" >&2
    echo "  Run 'bash scripts/benchmark.sh --update-baseline' to create it." >&2
    exit 1
fi

# Delegate comparison to Python3 (bash 3.2 on macOS lacks associative arrays)
python3 - "$BASELINE_FILE" "$CURRENT_FILE" "$THRESHOLD_PCT" "$VERBOSE" << 'PYEOF'
import sys

baseline_file  = sys.argv[1]
current_file   = sys.argv[2]
threshold_pct  = float(sys.argv[3])
verbose        = sys.argv[4] == "1"

def parse_bench(path):
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                vals[parts[0]] = float(parts[1])
    return vals

baseline = parse_bench(baseline_file)
current  = parse_bench(current_file)

regressions = 0
total = 0
passed = 0

for metric, bval in baseline.items():
    cval = current.get(metric)
    if cval is None:
        print(f"WARN  {metric}: not found in current output — skipping.")
        continue
    total += 1
    if bval == 0:
        pct = 0.0
    else:
        pct = (cval - bval) / bval * 100.0

    if pct > threshold_pct:
        print(f"FAIL  {metric}: {bval:.3f}ms -> {cval:.3f}ms  (+{pct:.1f}% > {threshold_pct:.0f}% threshold)")
        regressions += 1
    else:
        if verbose:
            print(f"PASS  {metric}: {bval:.3f}ms -> {cval:.3f}ms  ({pct:+.1f}%)")
        passed += 1

print(f"---")
print(f"Benchmark gate: {passed}/{total} passed, {regressions} regression(s) (threshold {threshold_pct:.0f}%).")

if regressions > 0:
    print(f"FAIL: {regressions} metric(s) regressed beyond threshold.", file=sys.stderr)
    sys.exit(1)

print("OK: no regressions detected.")
sys.exit(0)
PYEOF
