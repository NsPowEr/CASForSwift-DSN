#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
build_dir="$repo_root/build-bench"
policy_file="$repo_root/test/benchmarks/policy.json"
profile="release"
baseline_file=""
mode="run"
output_file=""
report_json=""
strict_metrics=0

resolve_path() {
  local raw_path="$1"
  if [[ "$raw_path" = /* ]]; then
    printf '%s\n' "$raw_path"
  else
    printf '%s\n' "$repo_root/$raw_path"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      build_dir="$2"
      shift 2
      ;;
    --baseline)
      baseline_file="$2"
      shift 2
      ;;
    --policy)
      policy_file="$2"
      shift 2
      ;;
    --profile)
      profile="$2"
      shift 2
      ;;
    --output)
      output_file="$2"
      shift 2
      ;;
    --report-json)
      report_json="$2"
      shift 2
      ;;
    --check)
      mode="check"
      shift
      ;;
    --update-baseline)
      mode="update"
      shift
      ;;
    --strict-metrics)
      strict_metrics=1
      shift
      ;;
    *)
      echo "argomento non supportato: $1" >&2
      exit 1
      ;;
  esac
done

build_dir="$(resolve_path "$build_dir")"
policy_file="$(resolve_path "$policy_file")"
if [[ -n "$baseline_file" ]]; then
  baseline_file="$(resolve_path "$baseline_file")"
fi
if [[ -n "$output_file" ]]; then
  output_file="$(resolve_path "$output_file")"
fi
if [[ -n "$report_json" ]]; then
  report_json="$(resolve_path "$report_json")"
fi

cmake -S "$repo_root" -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCAS_ENABLE_TESTING=OFF -DCAS_ENABLE_SANITIZERS=OFF
cmake --build "$build_dir" --target cas_benchmarks

tmp_output="$(mktemp "${TMPDIR:-/tmp}/cas-bench-output.XXXXXX")"
trap 'rm -f "$tmp_output"' EXIT
"$build_dir/cas_benchmarks" | tee "$tmp_output"

if [[ -n "$output_file" ]]; then
  mkdir -p "$(dirname "$output_file")"
  cp "$tmp_output" "$output_file"
fi

case "$mode" in
  run)
    if [[ -n "$report_json" ]]; then
      compare_args=(
        --policy "$policy_file"
        --profile "$profile"
        --current "$tmp_output"
        --report-json "$report_json"
      )
      if [[ -n "$baseline_file" ]]; then
        compare_args+=(--baseline "$baseline_file")
      fi
      if [[ "$strict_metrics" -eq 1 ]]; then
        compare_args+=(--strict-metrics)
      fi
      if ! python3 "$script_dir/benchmark_report.py" "${compare_args[@]}"; then
        echo "benchmark_report stato=warning motivo=regressione_rilevata_in_modalita_run" >&2
      fi
    fi
    ;;
  update)
    if [[ -z "$baseline_file" ]]; then
      baseline_file="$(python3 - "$policy_file" "$profile" <<'PY'
import json
import pathlib
import sys

policy_path = pathlib.Path(sys.argv[1]).resolve()
profile = sys.argv[2]
payload = json.loads(policy_path.read_text(encoding="utf-8"))
baseline = payload["profiles"][profile]["baseline"]
print((policy_path.parent.parent.parent / baseline).resolve())
PY
)"
    fi
    mkdir -p "$(dirname "$baseline_file")"
    cp "$tmp_output" "$baseline_file"
    ;;
  check)
    compare_args=(
      --policy "$policy_file"
      --profile "$profile"
      --current "$tmp_output"
    )
    if [[ -n "$baseline_file" ]]; then
      compare_args+=(--baseline "$baseline_file")
    fi
    if [[ -n "$report_json" ]]; then
      compare_args+=(--report-json "$report_json")
    fi
    if [[ "$strict_metrics" -eq 1 ]]; then
      compare_args+=(--strict-metrics)
    fi
    python3 "$script_dir/benchmark_report.py" "${compare_args[@]}"
    ;;
esac
