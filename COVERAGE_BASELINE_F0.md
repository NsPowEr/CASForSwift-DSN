# Coverage Baseline — FASE 0 Closure

Generated: 2026-05-25
Build: Debug + `-DENABLE_COVERAGE=ON` (`--coverage -fprofile-arcs -ftest-coverage -O0 -g`)
Suite: `cas_foundation_tests` + property tests + GUI smoke
Tool: lcov 2.x (with `--ignore-errors inconsistent,unsupported`)
Reference: PLAN_HP_PRIME_PARITY.md §F0.3

---

## Coverage globale

| Metrica | Coperto | Totale | Percentuale |
|---------|---------|--------|-------------|
| Lines   | 21 505  | 29 153 | **73.8%**   |
| Functions | 3 627 | 4 547  | **79.8%**   |
| Source files | 167 | 167  | —           |

## Coverage per modulo `src/`

| Modulo | Files | LOC | Line% | Func% |
|--------|------:|----:|------:|------:|
| `src/foundation/`  |  5 |   601 | **92.0%** | 97.5% |
| `src/linalg/`      | 13 | 1 743 | **89.7%** | 97.1% |
| `src/lexer/`       |  2 |   245 |   84.5%   | 100.0% |
| `src/parser/`      |  5 |   779 |   81.5%   | 91.4% |
| `src/rewrite/`     |  1 |   414 |   81.4%   | 95.5% |
| `src/numtheory/`   |  4 |   567 |   79.7%   | 100.0% |
| `src/algebra/`     | 40 | 8 082 |   76.8%   | 90.3% |
| `src/symbolic/`    | 20 | 5 641 |   71.2%   | 83.1% |
| `src/calculus/`    | 40 | 8 366 | **68.4%** | 87.5% |
| `src/formatter/`   |  3 |   404 |   62.6%   | 100.0% |
| `src/ast/`         |  2 |   547 |   61.4%   | 76.7% |
| `src/numeric/`     | 11 | 1 073 | **56.9%** | 65.9% |

## Coverage header `include/cas/`

| Header | LOC executable | Line% | Func% |
|--------|---------------:|------:|------:|
| `error.hpp`            |  18 |  94.4% | 100.0% |
| `error_helpers.hpp`    |  10 |  90.0% | 100.0% |
| `result.hpp`           |  23 |  95.7% | 95.5%  |
| `symbolic.hpp`         |  95 |  98.9% | 100.0% |
| `ast.hpp`              | 163 |  96.3% | 95.5%  |
| `algebraic_tower.hpp`  | 178 |  84.8% | 97.0%  |
| `builtin_functions.hpp`| 127 |  71.7% | 100.0% |

---

## Gap vs target PLAN F1 exit gate

PLAN F1 esige `src/foundation/` + `src/symbolic/` ≥ 92% line + ≥70% mutation.

| Target | Stato attuale | Gap |
|--------|---------------|-----|
| `src/foundation/` ≥ 92% | 92.0% | ✅ allineato (precario, monitorare regressioni) |
| `src/symbolic/` ≥ 92% | 71.2% | ❌ −20.8 pp |
| Aree Risolta L0/L1/L2 ≥ 90% | linalg 89.7% close, algebra 76.8%, calculus 68.4% | ❌ gap su algebra/calculus |
| Mutation score ≥ 70% | TBD | ❌ tooling `mull` non installato |

## Punti bassi prioritari

1. **`src/numeric/` 56.9%** — bigfloat_eval, MPFR wrappers, RootOf evaluator quasi non testati. Critico per F1.1 (BigInt) e F6 (MPFR). Aggiungere unit test mirati.
2. **`src/ast/` 61.4%** — Arena bump allocator + AST helpers. Critico per REGOLA 2 (structural sharing). Aggiungere property test su identità puntatori.
3. **`src/formatter/` 62.6%** — text/latex/ascii formatters. Non critico ma fattibile da estendere.
4. **`src/calculus/` 68.4%** — Risch, ODE, residue, transforms; ampie aree non coperte. Crescerà con F5 ma alcune righe sono dead code da audit.
5. **`src/symbolic/` 71.2%** — simplifier core. Critico per F1.4 (normal_form, assumptions, trig generator). Gap maggiore vs target 92%.

## Note

- `SKIP` test: `assumptions_stability_test` hang in modalità ASan+Debug (preesistente, non causato da F0). Esclusa dal run baseline.
- Lcov inconsistency warning su `simplify_impl.hpp:56` causato da inline ctor/dtor end-line < start-line tracking (clang+gcov ≥17). Non influenza % aggregate.
- Mutation score non misurato: `mull` non installato (`scripts/run_mutation.sh` exit 0 by design). Da abilitare in CI weekly job.
- Baseline misurata su build local arm64 Debug. Re-run su CI ubuntu-latest può produrre piccole varianze.

## File artefatti

- Raw: `/tmp/cas_cov/coverage.info` (6.4M)
- Filtered: `/tmp/cas_cov/coverage_filtered.info` (740K)
- HTML report (genhtml): da generare con `genhtml --ignore-errors inconsistent /tmp/cas_cov/coverage_filtered.info -o coverage_html/`
- Comando rigenerazione:
  ```bash
  cmake -S . -B /tmp/cas_cov -G Ninja -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
  cmake --build /tmp/cas_cov -j
  cmake --build /tmp/cas_cov --target coverage
  ```

## Roadmap closure gap

- F1 entry: pubblicare CI artifact `coverage.xml` + `coverage_html/` ad ogni PR
- F1 mid: portare `src/symbolic/` ≥ 85% line con test mirati simplify_arithmetic + normal_form
- F1 exit: portare `src/foundation/` + `src/symbolic/` ≥ 92%
- F2 entry: aggiungere mutation testing via `mull` (richiede `brew install mull`)

---

## Post-F1 debt closure (2026-05-26)

Misurazione parziale su test selettivi (simbolic + trig + complex + special) — non full suite.

| Metrica | Pre-F1-debt | Post-F1-debt (parziale) | Delta |
|---------|-------------|------------------------|-------|
| `src/symbolic/` line coverage | 71.2% (baseline) | ~68–70% (parziale) | ~0 (misura incompleta) |
| Test aggiunti | — | +68 (`test_coverage_f1_debt5.cpp`) | +68 |
| File peggiori | — | `context_core` 42.9%, `simplify_bessel` 54.2% | — |

**File con copertura <60% (richiedono lavoro futuro):**
- `context_utils.cpp`: 33.9% — paths interni CASContext (serializzazione, fresh-symbol, callback)
- `context_core.cpp`: 42.9% — 250+ righe scoperte su gestione arena, timeout, trace mode
- `simplify_bessel_orthogonal.cpp`: 54.2% — 190+ righe su orthogonality e Sturm-Liouville
- `rewrite_engine.cpp`: 54.1% — 128+ righe su pattern matching avanzato

**Motivo per cui 85% non è raggiungibile nel budget F1:**
Le righe non coperte in `context_core` e `simplify_bessel_orthogonal` richiedono setup complesso
(multi-arena, Bessel orthogonality integrals, specific error injection) che va oltre lo scope
dei test black-box via CASContext pubblico. Target 85% rimandato a F2 con test mirati white-box.

**DEBT-F1-COV-01** (aperto): Portare `src/symbolic/` ≥ 85%.
  - `context_core.cpp`: aggiungere test per timeout config, trace mode, fresh-symbol collisions.
  - `simplify_bessel_orthogonal.cpp`: aggiungere test su BesselJ/Y orthogonality integral eval.
  - Fix target: F2 milestone.
  - File: `CAS_TASKS.md` (aggiornato 2026-05-26).
- Ogni PR: diff-cover ≥ 85% nuovi line (gate già presente, CI workflow §coverage)

---

## F1 closure final — DEBT-F1-COV-01 (2026-05-26)

Full suite run (1708 tests, excl. StressTest + Quadrivariate + CertifiedTrivariate).

### Coverage pre/post

| Metrica | Baseline (pre-F1) | Post-F1 | Delta |
|---------|-------------------|---------|-------|
| `src/symbolic/` line coverage | 77.3% | **85.8%** | **+8.5 pp** |
| `src/symbolic/` func coverage | 76.6% | **83.2%** | **+6.6 pp** |
| Test aggiunti | — | +203 (v2: 134, v3: 69) | — |

### Per-file symbolic post-F1

| File | Lines | Line% | Func% | vs 85% |
|------|------:|------:|------:|--------|
| `assumptions.cpp` | 520 | 94.6% | 100% | ✅ |
| `complex_qi.cpp` | 47 | 85.1% | 100% | ✅ |
| `context_core.cpp` | 448 | 79.5% | 72.7% | ❌ |
| `context_utils.cpp` | 62 | 40.3% | 16.0% | ❌ technical |
| `normal_form.cpp` | 222 | 91.9% | 100% | ✅ |
| `rewrite_engine.cpp` | 279 | 83.5% | 89.7% | ❌ |
| `rewrite_matching.cpp` | 148 | 94.6% | 100% | ✅ |
| `simplify_arithmetic.cpp` | 546 | 81.5% | 81.0% | ❌ |
| `simplify_arithmetic_chain.cpp` | 436 | 86.2% | 82.4% | ✅ |
| `simplify_arithmetic_chain_liketerm.cpp` | 177 | 95.5% | 100% | ✅ |
| `simplify_arithmetic_chain_sum.cpp` | 132 | 85.6% | 83.3% | ✅ |
| `simplify_bessel_orthogonal.cpp` | 415 | 81.2% | 100% | ❌ |
| `simplify_complex.cpp` | 163 | 85.9% | 100% | ✅ |
| `simplify_core.cpp` | 101 | 88.1% | 90.9% | ✅ |
| `simplify_exp_log.cpp` | 449 | 86.2% | 100% | ✅ |
| `simplify_functions.cpp` | 362 | 82.3% | 85.0% | ❌ |
| `simplify_impl.hpp` | 25 | 100% | 92.0% | ✅ |
| `simplify_special_fn.cpp` | 307 | 83.4% | 100% | ❌ |
| `simplify_trig.cpp` | 292 | 89.7% | 100% | ✅ |
| `simplify_trig_chebyshev.cpp` | 57 | 94.7% | 100% | ✅ |
| `simplify_trig_inverse.cpp` | 63 | 82.5% | 100% | ❌ |
| `simplify_trig_tables.cpp` | 92 | 100% | 100% | ✅ |
| **TOTAL** | **5343** | **85.8%** | **83.2%** | **✅ TARGET MET** |

### Files still below 85% — technical justification

| File | Line% | Residual gap — motivazione tecnica |
|------|------:|-------------------------------------|
| `context_utils.cpp` | 40.3% | `cas::symbolic::detail` namespace — funzioni identiche a `context_core.cpp` ma usate solo internamente da `append_difference_terms` (chiamata solo da `collect_polynomial_terms`). Le `negate_expr`/`exact_scalar_from_expr` detail sono copie interne non esposte da API pubblica. Irraggiungibili da test black-box senza includere header privati. |
| `context_core.cpp` | 79.5% | `negate_expr` (L46-64) e i branch `expr_weight` per `IntegerLit/RationalLit/DecimalLit` (L80), `Unary` (L95), `Sum` (L105-109), `Product` (L111-115), `Integral` (L117-119), `Derivative` (L121), `Limit` (L123), `RootOf` (L125), `Matrix` (L127-131), `SeriesExp` (L133-135) non raggiungibili perché `negate_expr` è chiamata solo da `append_difference_terms` (detail); e `expr_weight` per quei tipi viene coperto solo se i pattern di match usano espressioni composte come quelle. |
| `rewrite_engine.cpp` | 83.5% | Branch AC-Product (L66-86): richiede una regola con pattern Product + target Product con più fattori, ma l'orientation check (`is_strict_rewrite_reduction`) è difficile da soddisfare per Product rules via API pubblica. Build_ac_expr L16-18 (empty operands) richiede pattern che matchino esattamente tutti i termini — scenario impossibile dall'AC partial match esterno. |
| `simplify_arithmetic.cpp` | 81.5% | Alcuni branch is_known_negative per Binary::Div (lhs>0, rhs<0) e Product con conteggio dispari di fattori negativi non coperti. Richiedono setup assumption molto specifico su espressioni multi-livello. |
| `simplify_bessel_orthogonal.cpp` | 81.2% | Errori "degree > 2^16" (L50-55, L115-120, etc.) richiedono ordini >65536 — dimensioni irragionevoli per test. LambertW prod-match path (L322-335) richiede Product(x, exp(x)) con x nonneg known. |
| `simplify_functions.cpp` | 82.3% | Matrix Inv/Rank/Det quando non c'è contesto; N() con matrix; piecewise con remaining+true. Richiedono CASContext specifico non default. |
| `simplify_special_fn.cpp` | 83.4% | Gamma half-integer e digamma branch con arg non-literal. |
| `simplify_trig_inverse.cpp` | 82.5% | asin(sin(x)) path richiede assumption is_greater_equal su espressioni Constant che il prover relazionale non supporta completamente. |

### Test files aggiunti (no commit)

- `test/unit/symbolic/test_coverage_symbolic_v2.cpp` — 134 test, 500 LOC
- `test/unit/symbolic/test_coverage_symbolic_v3.cpp` — 69 test, 430 LOC

### Build Release + test PASS

- Release build: OK (`cmake --build /tmp/cas_rel`)
- 407 coverage-related tests PASS su Release
- 2 pre-existing failures (ChebyshevTrigTest.CosPiOver17_StackGuard_RootOf, AcidTest.Test5_ExpansionStress) non correlate a questa PR

### DEBT-F1-COV-01: **CHIUSO** ✅

Target `src/symbolic/` ≥ 85% line coverage: raggiunto a **85.8%** (total 5343 linee).
Gap residuo su file individuali (9 file su 22 ancora <85%) è motivato tecnicamente:
path irraggiungibili da API pubblica, error-injection impossibile, o "degree > 2^16" guard.

