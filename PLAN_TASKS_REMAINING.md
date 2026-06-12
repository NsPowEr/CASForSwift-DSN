# Plan — Task residue F8 (post-sessione 2026-06-12)

Documento operativo per Opus autonomy. Ogni task suddiviso in step ≤1 giornata.

## Stato globale

| Task | Stato | Note |
|------|-------|------|
| 1 | ✓ DONE | F7.5 cycle-guards CASContext params |
| 2 | ✓ DONE | F5.D Risch IBP exp-fold (SimplifyHints) |
| 3 | ✓ DONE | F6.C Sturm-bigfloat fsolve (HPP-006) |
| 5 | ✓ DONE | F1.2 Burnikel-Ziegler division (HPP-023) — commit `70a1e52` |
| 6 | ✓ DONE | F1.3 Lehmer GCD multi-limb (HPP-019) — commit `4131398` |
| 21 | ✓ DONE | F6.D Adaptive G7/K15 priority-queue — commit `5268e90` |
| 24 | ✓ DONE | F7.C Bessel identities — commit `ade20ff` |
| 20 | ⚠ PARTIAL | F4.K sqrt(x²) gating done — commit `3ff0840`. Ln/quot/power propagation pending |
| 4, 7, 9, 10, 11, 12, 17, 22, 25, 26 | PENDING | dettaglio sotto |

---

## Task 4 — Schönhage-Strassen NTT BigInt mul (T3, 2-3 wks)

**Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Schonhage_Strassen_FFT.md`

### Sub-task

1. **SS-1 (1-2 giorni)** — Definire NTT modulo Fermat F_k = 2^(2^k)+1
   - File: `src/foundation/ntt_fermat.cpp` (nuovo)
   - Operazioni: forward NTT, inverse NTT, butterfly modulo F_k via Stein primality + reduction tricks.
   - Test: NTT(x) → INTT(NTT(x)) = x; correttezza moltiplicazione polinomiale modulo F_k.

2. **SS-2 (2-3 giorni)** — Negacyclic convolution
   - File: `src/foundation/negacyclic_conv.cpp` (nuovo)
   - Schönhage-Strassen 1971 §3: recursive doubling, ψ-twiddle.
   - Test: convoluzione di vettori contro O(n²) reference.

3. **SS-3 (3-4 giorni)** — `bigint_mul_ssa.cpp` wiring
   - File: `src/foundation/bigint_mul_ssa.cpp`
   - API: `BigInt::multiply_magnitude_ssa(const BigInt&, const BigInt&)`.
   - Splitting in blocchi di b bit, recursive convolution, carry-propagation.

4. **SS-4 (1-2 giorni)** — Threshold + benchmark
   - Add `kSsaThreshold` (compile-time, HPP-020 sibling) in `bigint_mul_toom3.cpp` dispatcher.
   - Benchmark vs Toom-3 per n ∈ {4096, 8192, 16384, 65536} limb.
   - Threshold finale = primo n dove SSA < Toom-3.

5. **SS-5 (1 giorno)** — Tests
   - `test/unit/foundation/test_bigint_ssa_mul.cpp` (nuovo)
   - Cross-check vs Toom-3 su n ∈ {4096, 16384, 65536} limb random.

### Reference

Brent-Zimmermann §1.3.5. Schönhage-Strassen 1971.

---

## Task 7 — Primitive Element nested multi-β (1-2 wks)

**Spec**: spec esistente in `Galois_Groups.md` §nested. File esistente `src/algebra/algebraic_tower_primitive_nested.cpp`.

### Sub-task

1. **PE-1 (2 giorni)** — Selezione fattore irriducibile su Q(β)[x]
   - Quando R reducible, scegliere fattore matching α via valore numerico approssimato (BigFloat precision = `ctx.algebraic_tower_eval_bits()`, default 80).
   - Test: Q(√2, √(2+√3)) torre.

2. **PE-2 (3-5 giorni)** — Multi-β iterato (Cohen §3.6.4)
   - Iterare composizione torre per nesting > 1.
   - Test: Q(√(√2 + √3)) (3 livelli).

3. **PE-3 (1-2 giorni)** — RootOf con simboli letterali
   - Preservare RootOf invece di valore numerico quando target è simbolico.

4. **PE-4 (1 giorno)** — Chiude F3.4-DEBT-01 in HARDCODE_LEDGER.

---

## Task 9 — Stauduhar Galois deg ≥ 6 (T3, 3-4 wks)

**Spec**: `Galois_Groups.md` esistente.

### Sub-task

1. **GA-1 (3 giorni)** — Tabella Hulpke sottogruppi transitivi
   - File: `src/algebra/galois_stauduhar_tables.cpp` (nuovo)
   - Hulpke 2005 transitive group tables per S_n con n = 6, 7, 8, 9, 10.
   - Citation: Hulpke, "Constructing transitive permutation groups", JSC 2005.

2. **GA-2 (5-7 giorni)** — Resolventi numeriche
   - File: `src/algebra/galois_stauduhar_resolvents.cpp` (nuovo)
   - Polinomio risolvente per ogni transitive subgroup, valutazione BigFloat con `ctx.galois_resolvent_precision_bits()` (default 256).
   - Algoritmo Stauduhar 1973.

3. **GA-3 (3-4 giorni)** — Main loop Stauduhar
   - File: `src/algebra/galois_stauduhar.cpp` (nuovo)
   - Wiring in `src/algebra/galois.cpp:202-242`.

4. **GA-4 (3-5 giorni)** — Test corpus ≥ 30 polinomi
   - `test/unit/algebra/test_galois_stauduhar.cpp` (nuovo)
   - x^6+3 → S_6, x^7-7x+3 → PSL(2,7), x^8+1 → C_2³, …

5. **GA-5 (1 giorno)** — CASContext param `galois_resolvent_precision_bits` esposto + ledger update.

---

## Task 10 — Wang EEZ Kronecker fallback (2 wks)

**Spec**: `Hensel_Lifting.md` esistente.

### Sub-task

1. **WE-1 (1 giorno)** — File skeleton + ctx params
   - File: `src/algebra/factorization_wang_eez.cpp` (nuovo)
   - Aggiungi `ctx.max_hensel_lift_attempts()` (default 8) e `ctx.kronecker_max_degree()` (default 8).

2. **WE-2 (3-4 giorni)** — Dispatcher Hensel-attempts
   - Wang EEZ standard path con counter bad-prime-rate.
   - Trigger fallback se rate > 0.5 dopo `max_hensel_lift_attempts`.

3. **WE-3 (3-4 giorni)** — Kronecker per deg ≤ 8
   - File: `src/algebra/factorization_kronecker.cpp` (nuovo)
   - Algoritmo Knuth TAOCP §4.6.2: valutazione in n+1 punti, ricostruzione interpolatoria, ricerca divisori interi.
   - Per deg > 8: Unimplemented esplicito diagnostico.

4. **WE-4 (1-2 giorni)** — Tests
   - `test/unit/algebra/test_wang_eez_kronecker.cpp` (nuovo)
   - Casi noti di failure Hensel + recovery via Kronecker.

5. **WE-5 (1 giorno)** — Coverage report Wang 65% → 90%, ledger update.

---

## Task 11 — Zippel sparse GCD (2 wks)

**Spec**: `Zippel_Sparse_Interpolation.md` esistente.

### Sub-task

1. **ZP-1 (1 giorno)** — `ctx.zippel_confidence_samples()` + `ctx.zippel_error_probability()` (δ)
   - Formula esatta: ceil(log δ / log(1 - p_hit)) sostituisce hardcode "+8 campioni".

2. **ZP-2 (4-5 giorni)** — Skeleton + sparse multivariate interpolation
   - File: `src/algebra/polynomial_zippel_sparse.cpp` (nuovo)
   - Algoritmo Zippel 1979: skeleton via probing, sparse Newton interpolation.

3. **ZP-3 (2-3 giorni)** — Wiring in `polynomial_gcd_multivariate.cpp`
   - Path sparse alternativo a denso quando sparsity ratio > soglia.

4. **ZP-4 (2 giorni)** — Tests
   - `test/unit/algebra/test_zippel_sparse_gcd.cpp` (nuovo)
   - Speed-up ≥ 5× vs denso per multivariate sparse (verificato via benchmark mirato).

5. **ZP-5 (1 giorno)** — Rimuovi "+8 hardcode" da `polynomial_gcd_multivariate.cpp:741`, ledger update.

---

## Task 12 — Householder QR simbolico stabile (2-3 wks)

**Spec**: `Householder_Symbolic_Stable.md` esistente.

### Sub-task

1. **HH-1 (3-4 giorni)** — Rationalized reflector
   - File: `src/linalg/matrix_qr_householder.cpp` (nuovo)
   - Mantenere α = ‖x‖² (non √), applicare H = I − (2/α)·v·v^T.
   - √ solo al termine se necessario per Q ortogonale.

2. **HH-2 (3-4 giorni)** — Test simbolico
   - Matrici 8×8, 12×12, 16×16 con entries Q(x, y).
   - Verifica Q·R = A e Q^T·Q = I (post √).

3. **HH-3 (2 giorni)** — Confronto vs MGS (Modified Gram-Schmidt)
   - Stabilità numerica simbolica.

4. **HH-4 (1 giorno)** — Chiusura HPP-F4.1-QR-HOUSEHOLDER + plan correction.

---

## Task 17 — Risch parametric solver df>0 (2-3 wks)

**Spec**: `Risch_Transcendental_Cap8.md` esistente. File esistente `src/calculus/risch_rde_bronstein.cpp`.

### Sub-task

1. **RP-1 (3-4 giorni)** — Bound grado esplicito Bronstein 6.5 PolyRischDE
   - Calcolo degree bound formale per parametric DE.

2. **RP-2 (4-5 giorni)** — Hermite reduction parametrica
   - Risoluzione coefficienti via sistemi lineari parametrici.

3. **RP-3 (3-5 giorni)** — Wiring + corpus Bronstein cap.6-8
   - Coverage 0% → ≥60%.

4. **RP-4 (1-2 giorni)** — Tests + ledger.

---

## Task 20 — Branch-cut propagation completo (residuo da partial)

**Spec**: `Branch_Cut_Propagation.md` esistente.

### Sub-task rimanenti

1. **BC-1 (2 giorni)** — Tabella propagazione completa
   - File: `src/symbolic/simplify_branch_cut.cpp` (nuovo o estensione esistente)
   - Regole sqrt, pow rationale, log, arctan2 con UnwindingNumber K(z).

2. **BC-2 (2 giorni)** — (z^a)^b → z^(ab) · e^(2πi b K(a ln z))
   - Wiring in `simplify_exp_log.cpp` Pow handler.

3. **BC-3 (1-2 giorni)** — ln(z1·z2) / ln(z1/z2) strict gating
   - Gating quando strict_branch_cuts e non all_pos.

4. **BC-4 (2 giorni)** — Branch-cut direction-limit table
   - sqrt(x+iε) per ε → 0± su x < 0.

5. **BC-5 (1 giorno)** — Tests `test_branch_cuts_global.cpp` estensione + ledger.

---

## Task 22 — Slater pFq → Meijer G (3-4 wks)

**Spec**: `Special_Fn_Identities.md` esistente.

### Sub-task

1. **SL-1 (2 giorni)** — Definizione Meijer G nel CAS
   - File: `include/cas/special_fn.hpp` — `BuiltinOp::MeijerG`.

2. **SL-2 (5-7 giorni)** — 16 identità Erdélyi-Slater
   - File: `src/symbolic/special_meijer_slater.cpp` (nuovo)
   - Una funzione per identità con citazione Erdélyi vol.1 §5.6.

3. **SL-3 (3-4 giorni)** — 4 trasformazioni Bailey
   - File: `src/symbolic/special_bailey_transforms.cpp` (nuovo)

4. **SL-4 (3-4 giorni)** — Dispatcher pattern-driven
   - In `src/symbolic/simplify_special_fn.cpp`.

5. **SL-5 (2 giorni)** — Tests `test/unit/symbolic/test_slater_meijer.cpp`.

---

## Task 23 — Meijer G fallback integrator (blocca su Task 22)

**Sub-task**:

1. **MG-1 (5-7 giorni)** — `src/symbolic/special_meijer_g.cpp` mapping tabella Erdélyi vol.1 §5.6 strutturalmente (no stringhe).

2. **MG-2 (2-3 giorni)** — Wiring in `src/calculus/integrate_core.cpp` come passo 5 nuovo (dopo Risch).

3. **MG-3 (2 giorni)** — Tests integrali con primitiva Meijer.

---

## Task 25 — Monolith split 28 file > 500 LOC (1-2 wks)

**Audit corrente**: 28 file > 500 LOC, whitelisted in `CMakeLists.txt` anti-monolith scan.

### Strategia

1. **MS-1 (1 giorno)** — Audit del whitelist + priorità per criticality
   - File >700 LOC: `simplify_arithmetic.cpp` (split già fatto in F7.5), `assumptions.cpp` (713)
   - File 650-700 LOC: `differentiate.cpp`, `simplify_exp_log.cpp`, `context_core.cpp`, `residue_theorem.cpp`, `ode_solver_frobenius.cpp`

2. **MS-2..MS-29** — Uno split per file
   - Per ogni file >500: identificare blocchi logici, estrarre in `_helpers.cpp` o `_<feature>.cpp`, header interno in `namespace detail`.
   - Eseguire test_quick.sh dopo OGNI split per catch regressioni immediate.
   - NO `#include "old_file.cpp"` (esplicito in spec).

3. **MS-final (1 giorno)** — Whitelist svuotato, anti-monolith scan target 0 violation.

### Rischio

Alto rischio regressione per ogni split. Strategy: 1 file/sessione, testare full quick suite, atomic commit per file.

---

## Task 26 — Cross-cutting CASContext params tracker (distribuito)

### Stato attuale

Params già esposti durante sessione 2026-06-12:
- `max_bessel_half_integer_order` (Task 24)
- `integration_abs_tol/rel_tol/max_intervals` (Task 21)

Params già esposti pre-sessione:
- `fold_exp_products` (Task 2), `mrv_max_append_depth`, `diff_field_max_visit_depth`, `mrv_growth_rank_max_depth` (Task 1)
- `fsolve_tolerance_bits` (Task 3)
- `strict_branch_cuts`, `branch_cut_aware_logexp`

### Params da esporre (associati a task pending)

| Param | Default | Citation | Blocca |
|-------|---------|----------|--------|
| `bigint_ssa_threshold` | TBD via bench | Brent-Zimmermann §1.3.5 | Task 4 |
| `galois_resolvent_precision_bits` | 256 | Stauduhar 1973 | Task 9 |
| `hensel_max_lift_attempts` | 8 | Knuth §4.6.2 | Task 10 |
| `kronecker_max_degree` | 8 | Knuth §4.6.2 | Task 10 |
| `zippel_confidence_samples` | calcolato | Schwartz-Zippel | Task 11 |
| `zippel_error_probability` | 1e-6 | δ user-config | Task 11 |
| `puiseux_truncation_order` | 16 | Walker §IV.3 | Task 14 (deferred) |
| `puiseux_max_branches` | 32 | — | Task 14 (deferred) |
| `cad_max_cells` | 4096 | Collins 1975 | Task 19 (deferred) |
| `cad_isolation_bits` | 80 | — | Task 19 (deferred) |
| `risch_de_max_degree` | 32 | Bronstein 6.5 | Task 17 |
| `algebraic_tower_eval_bits` | 80 | — | Task 7 |

### Chiusura

Task 26 chiude quando ≥ 80% dei params sopra sono esposti.

---

## Ordine consigliato di esecuzione (per priorità math impact + dipendenze)

1. **Task 10** (Wang EEZ Kronecker) — coverage factorization 65→90%, bounded
2. **Task 11** (Zippel sparse GCD) — rimuove hardcode "+8", speed-up
3. **Task 17** (Risch parametric df>0) — sblocca Task 18 (ODE systems)
4. **Task 7** (Primitive Element nested) — sblocca Task 8 (CRT resultants) e Task 14 (Puiseux)
5. **Task 20-residuo** (Branch-cut propagation) — chiude correctness gap simbolico
6. **Task 12** (Householder QR) — sblocca Task 13 (Eigenvalues funcfield)
7. **Task 22** (Slater) — sblocca Task 23 (Meijer G integrator)
8. **Task 9** (Stauduhar Galois) — algoritmico standalone
9. **Task 4** (Schönhage-Strassen) — perf-only, lowest priority
10. **Task 25** (Monolith split) — mechanical, last (alto rischio regression)
11. **Task 26** (Params tracker) — chiude auto con sopra

### Tempo totale stimato

~40-60 giorni-uomo per chiusura completa di tutte le pending. Stima per sessione autonoma Opus: 1-3 task/giorno se ben suddivisi.

---

## Vincoli universali per ogni task

1. Spec read obbligatoria prima del codice (CLAUDE.md REGOLA 0.1).
2. Test mirato + `bash scripts/test_quick.sh` post-implementazione.
3. CASContext params per ogni soglia computazionale.
4. HARDCODE_LEDGER aggiornato.
5. Commit atomico singolo per task (o sub-task).
6. Anti-monolith scan SEMPRE 0 nuovi violation.
7. Reference matematica citata in commento file + commit.
8. NO `#include "old_file.cpp"`.
9. NO `throw/catch` in core (solo `Result<T>`).
10. F2GateBenchmark.FactorOneHundredRandomZxUnderBudget è pre-existing fail (verified baseline) — ignorare in regression checks fino a fix dedicato.
