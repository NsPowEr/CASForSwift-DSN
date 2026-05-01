#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <milestone>" >&2
  exit 1
fi

milestone="$1"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_root="$(mktemp -d "${TMPDIR:-/tmp}/cas-$(printf '%s' "$milestone" | tr '[:upper:]' '[:lower:]').XXXXXX")"
trap 'rm -rf "$build_root"' EXIT

run_debug_suite() {
  local build_dir="$build_root/debug"
  cmake -S "$repo_root" -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug
  cmake --build "$build_dir"
  ctest --test-dir "$build_dir" --output-on-failure
}

run_ast_smoke() {
  local smoke_object="$build_root/cas_expr_smoke.o"
  cat <<'EOF' | c++ -std=c++20 -x c++ -I "$repo_root/include" - -c -o "$smoke_object"
#include "cas/ast.hpp"
int main() { cas::AstArena arena; return static_cast<int>(arena.size() != 0U); }
EOF
  rm -f "$smoke_object"
}

run_benchmark_gate() {
  local benchmark_build_dir="$build_root/bench"
  bash "$repo_root/scripts/benchmark.sh" \
    --build-dir "$benchmark_build_dir" \
    --check
}

case "$milestone" in
  M0)
    run_debug_suite
    run_ast_smoke
    ;;
  M1)
    run_debug_suite
    ;;
  M1b|M2|M2b|M3|M3b|M4|M5|M6|M7|M9)
    run_debug_suite
    run_benchmark_gate
    ;;
  *)
    echo "unsupported milestone: $milestone" >&2
    exit 1
    ;;
esac
