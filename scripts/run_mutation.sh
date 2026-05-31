#!/usr/bin/env bash
# F0.6 Mutation testing via mull-runner.
#
# Requires: mull-runner installed (brew install mull on macOS, or system package).
# Optional CI job (weekly schedule, non-blocking on PR).
#
# Usage:
#   scripts/run_mutation.sh [--build-dir build] [--module cas_foundation]
#   scripts/run_mutation.sh --all-modules
#
# Produces:
#   mutation-report/index.html  — per-module mutation score
#   mutation-report/summary.txt — console summary

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

BUILD_DIR="$repo_root/build-mutation"
REPORT_DIR="$repo_root/mutation-report"
MODULE="cas_foundation_tests"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-dir) BUILD_DIR="$2"; shift 2 ;;
        --module)    MODULE="$2";    shift 2 ;;
        --all-modules) MODULE="all"; shift   ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

# --- Check mull-runner is available ---
if ! command -v mull-runner &>/dev/null; then
    cat >&2 <<'EOF'
mull-runner not found. To install on macOS:
  brew install mull

On Linux (Ubuntu):
  apt-get install mull

This script is optional (weekly CI, non-blocking). Skipping.
EOF
    exit 0  # Exit 0: non-blocking
fi

MULL_VERSION=$(mull-runner --version 2>&1 | head -1 || echo "unknown")
echo "mull-runner: $MULL_VERSION"

# --- Configure Release build with bitcode for mull ---
echo "==> Configuring mutation build in $BUILD_DIR ..."
cmake -S "$repo_root" -B "$BUILD_DIR" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-fembed-bitcode-marker -g" \
    -DCAS_ENABLE_SANITIZERS=OFF \
    -DCAS_ENABLE_TESTING=ON

cmake --build "$BUILD_DIR" --target cas_foundation_tests -- -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"

mkdir -p "$REPORT_DIR"

TEST_BINARY="$BUILD_DIR/cas_foundation_tests"
if [[ ! -f "$TEST_BINARY" ]]; then
    echo "Error: test binary not found at $TEST_BINARY" >&2
    exit 1
fi

echo "==> Running mull-runner on $TEST_BINARY ..."
mull-runner \
    --reporters HTML \
    --reporters SQLite \
    --report-dir "$REPORT_DIR" \
    --timeout 30 \
    "$TEST_BINARY" \
    || true  # non-blocking: mutation failures don't fail the script

# --- Summary ---
SUMMARY="$REPORT_DIR/summary.txt"
{
    echo "Mutation testing summary — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Module: $MODULE"
    echo "Binary: $TEST_BINARY"
    echo "Report: $REPORT_DIR/index.html"
    echo ""
    echo "Interpretation:"
    echo "  Mutation score >= 70%: acceptable (F0.6 gate)"
    echo "  Mutation score < 70%:  indicates undertested code paths"
    echo ""
    echo "For CI integration: add a weekly scheduled job (non-blocking)."
    echo "See .github/workflows/ci.yml job 'mutation-weekly'."
} | tee "$SUMMARY"

echo ""
echo "Done. Report at: $REPORT_DIR/index.html"
