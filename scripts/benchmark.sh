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
# Su questa macchina (laptop, no core pinning) la varianza run-a-run sulle
# metriche sub-millisecondo arriva a 2× (misurato 2026-07-16): un singolo
# campione non è un segnale. Default: mediana per-metrica di 5 run.
runs=5

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
    --runs)
      runs="$2"
      shift 2
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
tmp_runs="$(mktemp "${TMPDIR:-/tmp}/cas-bench-runs.XXXXXX")"
trap 'rm -f "$tmp_output" "$tmp_runs"' EXIT

if ! [[ "$runs" =~ ^[0-9]+$ ]] || [[ "$runs" -lt 1 ]]; then
  echo "--runs richiede un intero >= 1 (ricevuto: $runs)" >&2
  exit 1
fi

# Misure sub-millisecondo sono garbage sotto carico (2026-07-16: load 33+ da
# processi utente ha raddoppiato poly_gcd in 2 minuti a codice identico).
# Warning esplicito, non blocco: chi misura deve SAPERE, il numero finisce
# comunque nel log.
load_1min="$(sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}' || uptime | awk -F'averages?: ' '{print $2}' | awk '{print $1}')"
ncpu="$(sysctl -n hw.ncpu 2>/dev/null || echo 8)"
if python3 -c "import sys; sys.exit(0 if float('$load_1min'.replace(',', '.')) > 0.5 * int('$ncpu') else 1)" 2>/dev/null; then
  echo "⚠ ATTENZIONE: load average 1-min = $load_1min con $ncpu CPU — la macchina è" >&2
  echo "  sotto carico; le misure saranno gonfiate e ad alta varianza. Per un" >&2
  echo "  baseline o un gate affidabile, rieseguire a macchina scarica." >&2
fi

for ((i = 1; i <= runs; i++)); do
  echo "── run $i/$runs ──" >&2
  timeout 300 "$build_dir/cas_benchmarks" | tee -a "$tmp_runs" >&2
done

# Mediana per metrica sui run raccolti (robusta agli outlier di scheduling).
python3 - "$tmp_runs" <<'PY' | tee "$tmp_output"
import statistics
import sys

samples: dict[str, list[float]] = {}
order: list[str] = []
with open(sys.argv[1], encoding="utf-8") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) != 2:
            continue
        name, value = parts
        try:
            v = float(value)
        except ValueError:
            continue
        if name not in samples:
            samples[name] = []
            order.append(name)
        samples[name].append(v)

for name in order:
    print(f"{name} {statistics.median(samples[name]):.3f}")
PY

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
