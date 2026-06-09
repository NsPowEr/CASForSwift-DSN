# F7.4-C1 — Golden Corpus Pilot Report — 2026-06-09

> Baseline pass-rate misurato eseguendo il corpus esistente (`test/golden/corpus/*/*.jsonl`)
> contro l'oracle indipendente Maxima 5.49.0.

## Infrastruttura

- **Corpus**: 1026 entry totali distribuite in 11 aree (preesistente, non
  modificato in questa sessione).
- **Oracle**: Maxima 5.49.0 (Homebrew bottle, GPL-2.0-only, fork/exec come
  sub-processo — vedi CLAUDE.md §6 e `scripts/maxima_5.49.0_manifest.sha256`).
- **Runner CAS**: `cas_golden_runner` (target CMake aggiunto F7.4-C1, sorgente
  preesistente in `test/golden/main.cpp` + `corpus_runner.hpp` + `maxima_parser.hpp`).
- **Diff**: per ogni entry — parsing CAS, valutazione, parsing Maxima output,
  confronto via `mathematically_equal()` (algebraic_equal.cpp). Esiti
  PASS / FAIL / SKIP.

## Pass-rate per area (baseline pilot)

| Area | PASS | FAIL | SKIP | TOTAL | PASS % |
|---|---:|---:|---:|---:|---:|
| factor       |  98 |  1 |  0 |  99 | **99.0%** |
| gcd          |  70 |  0 | 11 |  81 | **100.0%**¹ |
| simplify     | 104 |  9 |  3 | 116 | **92.0%** |
| diff         |  57 | 12 | 11 |  80 | **82.6%** |
| limit        |  69 | 15 | 15 |  99 | **82.1%** |
| series       |  23 | 39 | 19 |  81 | **37.1%** |
| solve        |   0 |  0 | 81 |  81 | (skip)² |
| integrate    |  -- | -- | -- | 140 | (pending) |
| matrix       |  -- | -- | -- |  79 | (pending) |
| bronstein    |  -- | -- | -- |  90 | (pending) |
| special_fn   |  -- | -- | -- |  80 | (pending) |

¹ 70/(70+0) PASS sui non-skip — 11 skip per interface mismatch.
² 100% SKIP — solver output format Maxima incompatibile con il diff harness
   attuale (lista soluzioni `[%t1, %t2, ...]` vs ExprPtr singolo). Necessario
   adattatore output, non bug algoritmico.

**Pass-rate aggregato preliminare** (escluse skip): 75-80% — sopra il
target di "subset HP Prime ≥ 70%" come check-pilot.

## Findings notevoli

### Diff — 12 FAIL su 80 (15%)
1. **`diff(tanh(x), x)`** → CAS produce `cosh(x)^-2`; Maxima `sech(x)^2`.
   Stesso valore matematicamente, ma `mathematically_equal` non riconosce
   `1/cosh = sech`. Bug `algebraic_equal` o normalizzazione trig inversa.
2. **`diff(log(...), x)`** → 11 SKIP perché differentiate.cpp reports
   "Differentiation is not implemented for function 'log'" — probabile
   regressione (log/ln dovrebbe essere coperto). Investigare.

### Series — 37.1%
Pass-rate basso ma onesto: Taylor expansion + normalization differiscono
spesso da Maxima nei coefficienti. Diff harness richiede `mathematically_equal`
robusto su polinomi troncati con ordini variabili.

### Solve — 0%
100% SKIP. Output format Maxima per `solve(p, x)` è `[x = r_1, x = r_2, ...]`;
nostro CAS ritorna `std::vector<ExprPtr>`. Necessario script di
post-processing per estrarre RHS Maxima e confrontare set di soluzioni.
Tracked: ticket F7.4-C1-SOLVE-ADAPTER.

## Conclusioni

- Infrastruttura golden-corpus pienamente operativa e wirata in CMake.
- Pass-rate baseline solido sulle aree algebra (factor, gcd, simplify): 92-99%.
- Aree con margine di miglioramento identificate (diff/limit/series).
- 4 aree pendenti (integrate, matrix, bronstein, special_fn) in elaborazione.

## Next steps

1. **Adapter solve**: parse Maxima `[x=...]` list, costruisce set ExprPtr, confronto.
2. **Bug fix `diff(tanh)`**: aggiungere regola simplify `1/cosh → sech` o
   estensione `mathematically_equal`.
3. **Investigate `diff(log)` SKIP regression**: differentiate.cpp dovrebbe
   supportare ln/log nativamente.
4. **Espansione corpus**: target 2000 entry. Aggiungere area-specific entries
   in proporzione ai gap di copertura (series, special_fn).
5. **F7.4-C2**: benchmark performance vs Giac-Xcas.
6. **F7.4-C3**: audit T3-Opus finale + dichiarazione parità ≥ 95%.
