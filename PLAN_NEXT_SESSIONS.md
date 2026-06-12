# Plan Next Sessions — Post-2026-06-12

> Documento operativo per Opus autonomy. Letto in `SESSION_BOOTSTRAP_OPUS.md §2`. Ogni sessione esegue il blocco corrente; al termine sposta items "DONE" sotto §Log e promuove la prossima sessione a "current".

## Vincolo invariante per ogni sessione

1. Spec read `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<task>.md` PRIMA del codice (CLAUDE.md REGOLA 0.1).
2. Ogni soglia/budget computazionale → ctx param (REGOLA ZERO).
3. Test mirato `--gtest_filter` + timeout 60-180s.
4. Suite quick gate post-modifica: `bash scripts/test_quick.sh` (≤700s).
5. Ledger update atomico per ogni task chiuso.
6. 1 task = 1 commit, conventional format, body multi-riga con citation.
7. NO disabilitazione test (REGOLA 0.2).
8. Anti-monolith scan 0 violation.
9. 3 fallimenti consecutivi → stop + report (anti-loop).

---

## Stato debiti aperti (snapshot 2026-06-12)

| ID | Descrizione | Effort | Blocking |
|----|-------------|--------|----------|
| Task 3.3 res / HC-KV-02 | Kovacic Laurent √r poli ordine pari ≥ 4 | 3-5 gg | nessuno (5.4 chiuso) |
| Task 3.3 res / HC-KV-03 | Kovacic Case 2/3 SL(2,C) | 1-2 wk | AlgebraicNumber tower estesa |
| Task 4.2 | Schönhage-Strassen FFT BigInt mul | 2-3 wk | — |
| Task 5.1 | Resultanti Modulari CRT | 1-2 wk | Task 4.2 raccomandato |
| Task 5.2 | Stauduhar Galois ≥ 6 | 3-4 wk | — |
| Task 5.5 | Collins CAD lifting | 2-3 wk | Task 5.1 |
| Task 7.1 | Slater pFq trasformazioni | 1-2 wk | — |
| Task 7.2 | Meijer G fallback integrator | 2 wk | Task 7.1 |
| HC-F8-PENDING-10 | Wang EEZ Kronecker | 2 wk | — |
| HC-F8-PENDING-11 | Zippel sparse GCD | 2 wk | — |
| HC-F8-PENDING-17 res | Risch parametric df>0 RP-2 Hermite | 1-2 wk | — |
| HC-F8-PENDING-20 res | Branch-cut BC-1..BC-5 | 1 wk | UnwindingNumber chiuso ✓ |
| HC-F8-QR-HOUSEHOLDER-BAILOUT | Soglie complessità symbolic QR | 3-5 gg | Branch-cut sqrt rule |
| HC-F8-FLAKY-COS-7PI-16 | Flaky order-dependent test | 1-3 gg | — |
| HC-F8-PENDING-25 | Monolith split 28 file >500 LOC | 1-2 wk | — (mechanical, ultima) |

## Ordine di priorità (math impact + dipendenze)

1. **Sessione 1** (current): Branch-cut BC-1..BC-5 → flaky bisect → QR bailout removal
2. **Sessione 2**: Wang EEZ (Task 10) + Zippel (Task 11)
3. **Sessione 3**: Risch parametric RP-2 + Kovacic Case 2 start
4. **Sessione 4**: HC-KV-02 Laurent √r + Kovacic Case 3
5. **Sessione 5+**: Stauduhar GA-1..GA-5 (Task 9)
6. **Sessione 6+**: Schönhage-Strassen SS-1..SS-5 (Task 4)
7. **Sessione 7**: CRT Resultanti (Task 5.1) — gated su SSA
8. **Sessione 8**: Slater + Meijer G (Task 7.1, 7.2)
9. **Sessione 9**: Collins CAD lifting (Task 5.5) — gated su CRT
10. **Sessione 10+**: Monolith split (1 file/sessione)

---

## Sessione 1 (CURRENT) — Branch-cut + flaky + QR debt

**Goal**: chiudere Task 20 res + ridurre baseline rosso + sbloccare HC-F8-QR-HOUSEHOLDER-BAILOUT.

**Budget**: 4-6 h.

### Step 1.1 — Branch-cut propagation (Task 20 res)

Spec: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Branch_Cut_Propagation.md`

- [ ] **BC-1** — Tabella propagazione completa
  - File: `src/symbolic/simplify_branch_cut.cpp` (nuovo, ≤ 500 LOC)
  - Header: `include/cas/simplify_branch_cut.hpp` (nuovo)
  - Regole: sqrt, pow rationale, log, arctan2 con UnwindingNumber K(z)
  - Wired in `simplify_impl.cpp` quando `ctx.strict_branch_cuts`
- [ ] **BC-2** — `(z^a)^b → z^(ab) · e^(2πi b K(a ln z))`
  - Wiring in `src/symbolic/simplify_exp_log.cpp` Pow handler
- [ ] **BC-3** — `ln(z1·z2)` / `ln(z1/z2)` strict gating
  - Gating quando `strict_branch_cuts` e non `all_pos` assumptions
- [ ] **BC-4** — Branch-cut direction-limit table
  - `sqrt(x+iε)` per ε → 0± su x < 0
  - File: `src/symbolic/simplify_branch_cut.cpp` (estensione BC-1)
- [ ] **BC-5** — Tests
  - File: `test/unit/symbolic/test_branch_cuts_global.cpp` (estensione)
  - Coverage: ogni regola BC-1..BC-4 con almeno 2 test (positive + negative path)
- [ ] **Ledger**: chiudere `HC-F8-PENDING-20` PARTIAL → DONE; aggiornare `TODO_PH8.md` §FASE 6 Task 6.2.

### Step 1.2 — Flaky cos7π/16 bisect (HC-F8-FLAKY-COS-7PI-16)

- [ ] **F-1** — Riproduzione deterministica
  - `--gtest_shuffle --gtest_random_seed=N` bisect su seed che riproduce fail
  - Salvare seed in ledger entry
- [ ] **F-2** — Strumentazione Chebyshev pipeline
  - Dump in `src/symbolic/simplify_special_fn.cpp` su cos(7π/16) attivata da env var
- [ ] **F-3** — Root cause fix
  - Reset esplicito di global statics in `SpecialFunctionsTest::SetUp` SE confermato
  - OR fix sorgente statica identificata
- [ ] **Ledger**: chiudere `HC-F8-FLAKY-COS-7PI-16` con seed + root cause.

### Step 1.3 — Symbolic QR bailout removal (HC-F8-QR-HOUSEHOLDER-BAILOUT)

Prerequisito: BC-1 con regola `sqrt(p)·sqrt(p) → p` quando `p ≥ 0` via assumptions.

- [ ] **HH-5** — Expose `ctx.symbolic_qr_max_norm_complexity()` ctx param (default 4 o non-bool gating).
- [ ] **HH-6** — Rimuovere `total_degree > 0` bail-out in `matrix_qr.cpp:173`.
- [ ] **HH-7** — Riabilitare `QRTest.SymbolicQR_DefaultSignConvention_2x2`.
- [ ] **Ledger**: chiudere `HC-F8-QR-HOUSEHOLDER-BAILOUT`; aggiornare `HC-F8-PENDING-12` PARTIAL → DONE.

### Output sessione 1

- 3 commit minimi (BC, flaky fix, QR bailout)
- Suite quick: 2358 → ≥ 2358 test, 2 FAIL → 1 FAIL (solo `F2GateBenchmark` baseline)
- Ledger: 3 entry chiuse
- `TODO_PH8.md` aggiornato in §Log

---

## Sessione 2 — Wang EEZ + Zippel

**Goal**: factorization coverage 65→90% + rimozione hardcode `+8` campioni.

**Budget**: 4-6 h.

### Step 2.1 — Wang EEZ Kronecker (Task 10)

Spec: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Hensel_Lifting.md`

- [ ] **WE-1** — File skeleton + ctx params
  - `src/algebra/factorization_wang_eez.cpp` (nuovo)
  - `ctx.max_hensel_lift_attempts()` default 8
  - `ctx.kronecker_max_degree()` default 8
- [ ] **WE-2** — Dispatcher Hensel-attempts con counter bad-prime-rate
- [ ] **WE-3** — Kronecker per deg ≤ 8
  - `src/algebra/factorization_kronecker.cpp` (nuovo)
  - Knuth TAOCP §4.6.2: valutazione in n+1 punti, ricostruzione interpolatoria
- [ ] **WE-4** — Tests `test/unit/algebra/test_wang_eez_kronecker.cpp`
- [ ] **WE-5** — Coverage report Wang + ledger.

### Step 2.2 — Zippel sparse GCD (Task 11)

Spec: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Zippel_Sparse_Interpolation.md`

- [ ] **ZP-1** — Ctx params
  - `ctx.zippel_confidence_samples()` (calcolato)
  - `ctx.zippel_error_probability()` (default 1e-6)
  - Formula: `ceil(log δ / log(1 - p_hit))` sostituisce hardcode `+8`
- [ ] **ZP-2** — Sparse multivariate interpolation
  - `src/algebra/polynomial_zippel_sparse.cpp` (nuovo)
  - Zippel 1979: skeleton via probing, sparse Newton
- [ ] **ZP-3** — Wiring in `polynomial_gcd_multivariate.cpp`
- [ ] **ZP-4** — Tests `test/unit/algebra/test_zippel_sparse_gcd.cpp`
- [ ] **ZP-5** — Rimuovere `+8` hardcode da `polynomial_gcd_multivariate.cpp:741`.

### Output sessione 2

- 2 commit (Wang, Zippel)
- Ledger: 2 entry chiuse (HC-F8-PENDING-10, HC-F8-PENDING-11)

---

## Sessione 3 — Risch parametric RP-2 + Kovacic Case 2 start

### Step 3.1 — Risch RP-2 Hermite parametrica (HC-F8-PENDING-17 res)

Spec: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Transcendental_Cap8.md`

- [ ] **RP-2** — Hermite reduction parametrica
  - Risoluzione coefficienti via sistemi lineari parametrici
  - File: `src/calculus/risch_rde_bronstein_hermite.cpp` (nuovo)
- [ ] **RP-3** — Wiring + corpus Bronstein cap 6-8 (coverage 0% → ≥60%)
- [ ] **RP-4** — Tests + ledger.

### Step 3.2 — Kovacic Case 2 (Task 3.3 res / HC-KV-03 partial)

- [ ] **K2-1** — Estendere `case2_omega` usando RootOf+IsolatingBound (5.4 chiuso ✓)
  - File: `src/calculus/ode_kovacic_case2.cpp` (nuovo)
- [ ] **K2-2** — Tests Bessel/Airy classici
- [ ] **K2-3** — Diagnostico Case 3 (SL(2,C) finiti subgroups: A₄, S₄, A₅).

---

## Sessione 4 — Kovacic Laurent + Case 3

### Step 4.1 — HC-KV-02 Laurent √r poli ordine pari ≥ 4

- [ ] **KV2-1** — File `src/calculus/ode_kovacic_puiseux.cpp` (nuovo)
- [ ] **KV2-2** — Laurent expansion di √r per poli ordine 2k
- [ ] **KV2-3** — Wiring in `compute_r` + `case1_omega`
- [ ] **KV2-4** — Tests.

### Step 4.2 — Kovacic Case 3

- [ ] **K3-1** — Tabella sottogruppi finiti SL(2,C) (A₄, S₄, A₅)
- [ ] **K3-2** — `case3_omega` con polinomio minimo del campo algebrico
- [ ] **K3-3** — Tests Lamé/Heun equations.

---

## Sessioni 5-10+ — Vedi §"Ordine di priorità" sopra

Sessioni successive: Stauduhar Galois → Schönhage-Strassen → CRT Resultanti → Slater+Meijer G → Collins CAD → Monolith split.

Ogni sessione segue lo stesso protocollo: spec read → audit → ctx params → impl → test mirato → quick suite → ledger atomico → commit.

---

## Log sessioni (più recente in alto)

### 2026-06-12 (sessione audit + commit WIP enorme)
- 14 commit landed (`6bff716..895840a`).
- Closed: AST F8.0 split, Sturm-bigfloat MPFR polish, Q(x)[α] arithmetic, Kovacic Case 1, Risch tower descent, MGS→Householder QR, extended_real dedup, cycle-guard ctx params, RootOf bound serialisation, UnwindingNumber registration, SimplifyHints + Sum degree-order + differentiator visitor.
- Ledger: HPP-F4.1-QR-HOUSEHOLDER riaperto PARTIAL, HC-F8-PENDING-12 partial, HC-F8-QR-HOUSEHOLDER-BAILOUT nuovo, HC-F8-FLAKY-COS-7PI-16 nuovo.
- Suite quick: 2358 test, 2354 PASS, 2 SKIPPED, 2 FAIL (1 baseline + 1 flaky).
- **Next session pickup**: Sessione 1 sopra — Branch-cut BC-1..BC-5.
