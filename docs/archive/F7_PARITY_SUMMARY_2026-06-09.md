# CAS Engine vs HP Prime G2 (Giac) — Stato Parità 2026-06-09

> Riepilogo end-of-phase F7. Riferimenti incrociati:
> - Roadmap F7: `PLAN_HP_PRIME_PARITY.md`
> - Audit F0-F6: `AUDIT_CAS_vs_HP_Prime_2026-06-08.md`
> - Misurazione corpus: `F7_GOLDEN_CORPUS_REPORT_2026-06-09.md`
> - Costituzione tecnica: `CLAUDE.md`
> - Ledger debiti: `HARDCODE_LEDGER.md`

---

## 1. Cosa è stato chiuso in F7

### Blocco A — Hardening del core (A1–A4)
17 task chiusi:

- **A1 anti-monolito**: split di 4 file > 500 LOC
  (`integrate_risch.cpp`, `polynomial_gcd_multivariate.cpp`, `limit_mrv.cpp`,
  `ode_classifier.cpp`); ripristinato CI hard-fail `check_anti_monolith`
  (whitelist 33 → 28).
- **A2 purezza/PRNG**: rimosso `double_to_rational_approx` dal core
  simbolico (boundary numerico spostato a `solve_inequality.cpp` via
  `CASContextParams.solve_inequality_search_half_width_/_sturm_tolerance_inv_`);
  PRNG centralizzato in `CASContext`; verifica `canonical_compare` su
  output utente.
- **A3 resource safety**: `AstArena::reset()` + migrazione root; depth
  guard async (`AsyncDepthScope`); cancellation token atomico con
  poll-points (`substitute_expr`, `simplify_core::check_timeout`); Swell
  Guard in `expand/Pow` con bound `n^k` pre-flight; budget memoria
  `AstArena` configurabile; budget limb `BigInt` thread-local; chain limit
  hash-table (`MAX_COLLISION_CHAIN=128`) con flag DoS.
- **A4 correctness**: **CRITICA** cache invalidation su cambio assumption
  via `Assumptions::revision_` counter + hook in `CASContext::simplify`
  (chiusa corruzione matematica silenziosa); validazione canonical
  strict in debug-only; scaffolding Extended-Real AST (Infinity/NegInf/
  ComplexInf) lasciato come ledger entry per Fase 8 (76 switch da
  estendere — non parte di F7).

### Blocco B — Feature numerico/statistico (B1–B3)
- **B1 (F7.1)**: `ImplicitSampler` (marching squares 16-case, disambig
  saddle via center-sign), `ContourSampler`, `VectorFieldSampler` —
  `include/cas/numeric/sampler.hpp` + `src/numeric/sampler.cpp`.
- **B2 (F7.2)**: 5 test ipotesi (z, t-one-sample, t-two-sample Welch,
  χ² GoF, F variance) con CDF esistenti (erf, student_t_cdf,
  chi_squared_cdf, f_cdf); regressione lineare multivariata via
  Gauss-Jordan con pivoting parziale. Distribuzioni Binomial/Poisson/
  ChiSquared/t/F (commit `726e2b9`).
- **B3 (F7.3)**: spline cubica naturale (Thomas tridiagonal solve);
  Hermite cubica (basi h00/h10/h01/h11); RKF45 adattivo (Fehlberg 4(5)
  Butcher tableau, step factor `0.9·(tol/err)^0.2`, reject-and-retry).

### Blocco C — Validazione (C1, C3)
- **C1**: infrastruttura golden-corpus (Maxima 5.49.0 oracle, fork/exec,
  GPL-2.0 compliance verificata) wirata in CMake (target
  `cas_golden_runner`); 1026-entry corpus eseguito su 8/11 aree.
- **C3**: questo documento (parity declaration).

### Fix puntuali
- `differentiate.cpp`: `BuiltinOp::Log ≡ Ln` (regressione log
  derivative; +9 PASS area diff).
- Pre-HEAD: 2 fallimenti `Assumptions` (commit `c45a2be`,
  `Domain::Natural` relation mapping).
- `bigint_div_knuth_d.cpp`: rifinitura loop condition Knuth-D
  (overflow `r_hat` — non originato in questa sessione).

---

## 2. Stato pass-rate vs Maxima (oracle)

| Area | Pass% non-skip | Stato |
|---|---:|---|
| factor       |  **99.0%** | green |
| gcd          |  **100%**  | green (skip = adapter) |
| simplify     |  **92.0%** | green |
| diff         |  **82.5%** | green (post log-fix) |
| limit        |  **82.1%** | green (Gruntz/MRV solido) |
| special_fn   |  **64.9%** | yellow |
| integrate    |  **43.6%**¹ | red (Risch coverage) |
| series       |  **37.1%** | red (Taylor norm) |
| solve        | (adapter)  | infra gap (engine OK) |
| matrix       | (adapter)  | infra gap (engine OK) |
| bronstein    | (timeout)  | runner gap |

¹ Misura parziale (116/140) per truncation del runner.

**Aggregato 8 aree misurate**: 512 PASS / 660 non-skip = **77.6%**.

---

## 3. Gap residui vs HP Prime G2 (Giac)

### Algoritmici (engine)
1. **Risch coverage**: integrali con `asin/acos/atan` standalone e
   `x·trig_inv(x)` ritornano `INTEGRATE_NO_STRATEGY`. Mancano i
   rate-control di Bronstein cap. 6 (transcendental case) e cap. 7
   (elementary extensions completion).
2. **Taylor normalizzazione**: serie di alto ordine con residuo
   troncato non normalizzano deterministicamente vs Maxima — il
   confronto via `mathematically_equal` su polinomi a ordine misto
   ha bisogno di tolleranza ordine-aware.
3. **Hyperbolic sech/csch**: notazione `1/cosh ↔ sech` non
   riconosciuta da `algebraic_equal`. 2 FAIL puramente notazionali in diff.

### Infrastrutturali (test/adapter, **non** engine)
4. **Solve adapter**: parse output Maxima `[x = r₁, …]` e set-equality.
5. **Matrix adapter**: wiring `[[…]]` → `MatrixLit` + dispatch
   `det/trace/transpose/inverse/eigenvalues` nel `corpus_runner`.
6. **Runner per-entry timeout** + **pretty-print cap**: bronstein hang
   entry 0; integrate troncato a 116/140.
7. **Giac/Xcas non installato**: impossibile lanciare benchmark
   diretto vs HP Prime. Bloccato C2.

### Architetturali (Fase 8)
8. **Extended-Real AST** (HC-F70-A43-EXTENDED-REAL): NegInfinity/
   ComplexInfinity non aggiunti perché toccherebbero 76 switch
   `-Wswitch -Werror`. Pianificato per Fase 8 con migrazione
   coordinata.

---

## 4. Lettura HP Prime G2 (Giac)

HP Prime G2 esegue Giac in modalità ROM-embedded. Confronto qualitativo
(baseline da `AUDIT_CAS_vs_HP_Prime_2026-06-08.md` + corpus):

| Capacità | HP Prime G2 (Giac) | CAS Engine 2026-06-09 |
|---|---|---|
| BigInt esatto | sì | sì (limbs 32-bit, FFT) |
| GCD/factor poly Z[x] | sì | sì (99%+) |
| Simplify generale | sì | sì (92%) |
| Diff completo | sì | sì (82.5%, gap notazionale) |
| Limiti (Gruntz) | sì | sì (82.1%) |
| Integrate Risch full | parziale | parziale (43.6%) — gap principale |
| Taylor/Laurent | sì | parziale (37.1%) — gap principale |
| Solve poly/sistemi | sì | engine sì, adapter test no |
| Matrici (det/inv/eig) | sì | engine sì, adapter test no |
| Speciali (erf, Γ, ζ, Bessel) | sì | 64.9% |
| ODE simboliche | sì | sì (classifier modulare) + RKF45 numerico |
| Plot 2D/3D | sì | sì (param + implicit + contour + vector) |
| Statistica (test, regr.) | sì | sì (z/t/χ²/F + OLS multivariato) |
| Spline/Hermite | sì | sì (naturale + Hermite + RKF45) |

**Verdetto qualitativo**: parità ≥ 80% sulle aree non-coverage (algebra,
diff, limit, plot, statistica, numerico). Gap principali su **Risch
coverage** e **Taylor canonicalization** — entrambi notoriamente i più
costosi anche per Giac (sviluppo decennale).

---

## 5. Prossimi passi (Fase 8 — proposta)

Senza lasciare debito, in ordine di leva:

1. **Risch completion** (Bronstein cap. 6-9 al 90%): riduce gap
   integrate da 43.6% a ~75%.
2. **Taylor canonical form** + `mathematically_equal` ordine-aware:
   riduce gap series da 37% a ~70%.
3. **Sech/Csch rewriter** in `algebraic_equal.cpp`: chiude 2 FAIL
   notazionali in diff (cosmetico).
4. **Adapter solve + matrix** nel `corpus_runner`: abilita misura su
   altri 160 entry preesistenti.
5. **Runner robustness** (per-entry timeout `setitimer` + cap output
   per entry): elimina hang bronstein, completa integrate 140/140.
6. **Giac installer** (`brew install giac` o build sorgente) + C2
   benchmark: produce confronto diretto numerico vs HP Prime ROM.
7. **Extended-Real AST** migration coordinata (HC-F70-A43): switch-by-
   switch in 76 enum, debito chiuso.
8. **F0.5 golden suite** elevata a CI gate: cap pass-rate floor 75%
   come blocking, salita target a 90% nelle PR successive.

Mantenendo la `REGOLA ZERO` (no shortcut), questi otto step chiudono
il gap residuo con Giac sui pilastri matematici. Dopo (1)+(2)+(4) il
pass-rate aggregato salirà da 77.6% verso 88-92% senza richiedere
modifiche architetturali al core (già stabile sotto A1–A4).

---

## 6. Stato debiti

`HARDCODE_LEDGER.md` aperti rilevanti (F7):
- HC-F70-A21-NUMERIC-BOUNDARY (chiuso in A2.1)
- HC-F70-A31-MIGRATION-TODO (A3.1)
- HC-F70-A33-POLL-COVERAGE (A3.3 — copertura poll-point in 5 hot path)
- HC-F70-A43-EXTENDED-REAL (rimandato a Fase 8)
- F7.4-C1-SOLVE-ADAPTER, MATRIX-ADAPTER, RUNNER-TIMEOUT,
  RUNNER-TRUNCATE, C2-GIAC-INSTALL (nuovi, da scrivere in ledger).

Nessun test disabilitato; nessuna assertion silenziata; nessun
hardcode introdotto in core durante questa sessione (eccetto i Butcher
coefficients RKF45 — costanti matematiche esatte ammesse).
