---
name: benchmark-gate
description: Esegue scripts/benchmark.sh (profilo release) e confronta i risultati contro la baseline (test/benchmarks/baseline_release.txt) tramite scripts/check_bench_regression.sh. Fallisce su regressione oltre soglia configurabile. Gate OBBLIGATORIO pre-merge per CLAUDE.md §WORKFLOW.
disable-model-invocation: true
---

# benchmark-gate

Gate prestazionale pre-merge. Lancia la suite benchmark del CAS e blocca se una metrica regredisce oltre soglia.

## Quando usare

- Prima di un commit/PR che tocca core simbolico, simplifier, BigInt, polynomial, Groebner, integrazione, parser/printer.
- Dopo refactor di un modulo "hot path".
- Per stabilire una nuova baseline dopo ottimizzazione verificata.

## Protocollo

1. **Build release** (se non già aggiornata):
   ```bash
   cmake --build build-bench --config Release -j
   ```
2. **Esegui benchmark**:
   ```bash
   bash scripts/benchmark.sh
   ```
   Output canonico → `build-bench/benchmark_results.txt` (verifica path effettivo nello script).
3. **Diff vs baseline**:
   ```bash
   bash scripts/check_bench_regression.sh \
     --current build-bench/benchmark_results.txt \
     --baseline test/benchmarks/baseline_release.txt
   ```
   Exit 1 se una metrica supera `REGRESSION_THRESHOLD_PCT` (default 5%).
4. **Report**: incolla diff per metrica regredita; identifica modulo responsabile.

## Vincoli

- **NO** modifica baseline senza giustificazione misurata (commit dedicato `perf: update baseline — <reason>`).
- **NO** disabilitazione benchmark per "far passare" il gate (CLAUDE.md §REGOLA 0.2 — divieto disabilitare test).
- Se la regressione è strutturale e accettata (es. correzione algoritmica che migliora correttezza a costo perf), documenta in `HARDCODE_LEDGER.md` o `PERF_LEDGER.md` con motivazione.

## Output atteso

- PASS → log + "OK: no regression vs baseline".
- FAIL → elenco metriche regredite con `current`/`baseline`/`delta%`; investigare prima del merge.
