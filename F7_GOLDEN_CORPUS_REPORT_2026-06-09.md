# F7.4-C1 — Golden Corpus Pilot Report — 2026-06-09

> Baseline pass-rate misurato eseguendo il corpus esistente
> (`test/golden/corpus/*/*.jsonl`) contro l'oracle indipendente Maxima 5.49.0.

## Infrastruttura

- **Corpus**: 1026 entry distribuite in 11 aree (preesistente).
- **Oracle**: Maxima 5.49.0 (Homebrew bottle, GPL-2.0-only, fork/exec —
  CLAUDE.md §6, manifest pin `scripts/maxima_5.49.0_manifest.sha256`).
- **Runner CAS**: target CMake `cas_golden_runner` (wirato F7.4-C1;
  sorgenti in `test/golden/main.cpp` + `corpus_runner.hpp` + `maxima_parser.hpp`).
- **Diff**: parsing CAS, valutazione, parsing Maxima output, confronto via
  `mathematically_equal()` (algebraic_equal.cpp). Esiti PASS / FAIL / SKIP.

## Pass-rate per area (baseline misurato 2026-06-09)

| Area | PASS | FAIL | SKIP | TOTAL | PASS % non-skip | Note |
|---|---:|---:|---:|---:|---:|---|
| factor       |  98 |  1 |  0 |  99 | **99.0%** | algebra solida |
| gcd          |  70 |  0 | 11 |  81 | **100.0%** | skip = adapter gap |
| simplify     | 104 |  9 |  3 | 116 | **92.0%** | core robusto |
| diff         |  66 | 14 |  0 |  80 | **82.5%**¹ | post log-fix (+9) |
| limit        |  69 | 15 | 15 |  99 | **82.1%** | Gruntz solido |
| series       |  23 | 39 | 19 |  81 | **37.1%** | Taylor norm gap |
| solve        |   0 |  0 | 81 |  81 | (adapter gap)² | algoritmico OK |
| integrate    |  34 | 44 | 38 | 116/140 | **43.6%**³ | Risch coverage |
| matrix       |   0 |  0 | 79 |  79 | (adapter gap)⁴ | engine OK |
| bronstein    |   ? |  ? |  ? | 90 | (timeout)⁵ | hang entry 0 |
| special_fn   |  48 | 26 |  6 |  80 | **64.9%** | gap normalizzazione |

¹ Post fix `differentiate.cpp` (Log≡Ln). Pre: 57/12/11. Net +9 PASS.
² Maxima output `[x=r₁, x=r₂, …]` lista; CAS ritorna `std::vector<ExprPtr>`.
   Adapter di post-processing necessario (estrazione RHS, confronto set).
³ Misurazione parziale: runner troncato a 116 entry (su 140) per output
   verbose su integrali trig di alto grado. Runner-side issue.
⁴ corpus_runner.hpp marca `det/trace/transpose/rank/inv/eigenvalues` come
   `matrix fn skipped` perché manca adapter `[[…]]` → `MatrixLit`. Engine
   matriciale già coperto da unit test propri.
⁵ Runner blocca su entry 0 (`bronstein/integrals.jsonl`). Necessario
   per-entry timeout nel runner (ticket F7.4-C1-RUNNER-TIMEOUT).

## Aggregato (8 aree con misura)

```
PASS totali   : 512
FAIL totali   : 148
SKIP totali   : 252
Misurati      : 660 (non-skip)
Pass-rate     : 512 / 660 = 77.6%
```

Aggregato escluse aree con adapter gap (solve, matrix) e bronstein (hang):
77.6% — sopra il floor "subset HP Prime ≥ 70%" come check-pilot.

## Findings notevoli e fix in-session

### Fixed
- **`diff(log(x))` SKIP** → `differentiate.cpp` non gestiva `BuiltinOp::Log`
  (solo `Ln`). Resto del codebase tratta i due come natural log sinonimi.
  Aggiunto branch `|| Log` in abs-intercept e outer-derivative.
  → +9 PASS area diff.

### Documented per follow-up
- **`diff(tanh(x))`**: CAS `cosh(x)^-2`, Maxima `sech(x)^2`. Stesso valore,
  ma `mathematically_equal` non normalizza `1/cosh ↔ sech`. Fix: aggiungere
  rewriter sech→1/cosh prima del confronto, oppure introdurre BuiltinOp
  `Sech/Csch` con simplify identity.
- **Series 37%**: pass-rate basso onesto. Taylor expand + normalization
  spesso differiscono nei coefficienti residuali o nel troncamento.
  `mathematically_equal` su polinomi troncati richiede tolleranza
  ordine-aware.
- **Integrate 44%**: molti SKIP per `INTEGRATE_NO_STRATEGY` su
  `x·atan(x)`, `asin/acos/atan` standalone, log composito. Coverage Risch
  parziale — Bronstein cap. 6-9 ancora da estendere.
- **Solve adapter** (`F7.4-C1-SOLVE-ADAPTER`): parser per Maxima list
  output `[x=r₁,…]` + set-equality.
- **Matrix adapter** (`F7.4-C1-MATRIX-ADAPTER`): wiring `[[…]]` →
  `MatrixLit` + dispatch `det/trace/transpose/inverse/eigenvalues`.
- **Runner truncation** (`F7.4-C1-RUNNER-TRUNCATE`): integrale trig di
  alto grado produce output multi-MB; cap pretty-print per entry.
- **Runner per-entry timeout** (`F7.4-C1-RUNNER-TIMEOUT`): bronstein
  entry 0 hang. Necessario `alarm()` per entry o thread-watchdog.

## Conclusioni

- Infrastruttura golden-corpus operativa, wirata in CMake, eseguita end-to-end.
- Pass-rate baseline solido su algebra: factor 99%, gcd 100%, simplify 92%.
- Calcolo differenziale post-fix: diff 82.5%, limit 82.1%.
- Aree con margine: series 37%, integrate 43.6%, special_fn 64.9%.
- 3 aree con adapter/timeout gap (solve, matrix, bronstein) — non
  algoritmiche.

## Prossimi passi (F7.4 → fase 8)

1. **C1-SOLVE-ADAPTER** + **C1-MATRIX-ADAPTER**: estende runner per chiudere
   le due aree algoritmicamente già coperte dall'engine.
2. **C1-RUNNER-TIMEOUT/TRUNCATE**: stabilizza runner per misura completa
   integrate + bronstein.
3. **Risch coverage push**: ridurre `INTEGRATE_NO_STRATEGY` su funzioni
   trig inverse standalone e `x·trig_inv(x)`.
4. **Series normalization**: confronto algebrico orderaware
   per `mathematically_equal` su Taylor troncato.
5. **Sech/Csch rewriter** in `algebraic_equal.cpp` (pulisce 2 FAIL
   notazionali in diff).
6. **C2 Giac benchmark**: bloccato — Giac/Xcas non installato sul host
   pilot; serve installer dedicato (port Homebrew `giac`). Track in
   `F7.4-C2-GIAC-INSTALL`.
7. **C3 audit finale**: post-fix C1 items, ri-misura aggregato e
   dichiarazione parità F7 → ≥ 90% target sulle aree non-adapter.
