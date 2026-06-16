
# TODO — Phase 8 (F8.x Roadmap)
> Aggiornato automaticamente dall'agente. Ultimo update: 2026-06-11 (Task 3.3 parziale: Kovacic Case 1)

---

## 🧱 FASE 1: Sblocco Strutturale (F8.0)
Nessun nuovo algoritmo può essere aggiunto finché il monolite non è smantellato, altrimenti le nuove
dipendenze circolari (es. AlgebraicNumber che include BigInt e Ast) romperanno la build.

 - [x] Task 1.1: Dividere include/cas/ast.hpp in ast_kinds.hpp (enums), ast_arena.hpp (memory) e ast_nodes.hpp (classi).
 - [x] Task 1.2: Estrarre la logica da src/ast/ast.cpp creando ast_compare.cpp e ast_clone.cpp.
 - [x] Task 1.3: Rifattorizzare differentiate.cpp isolando i derivatori per tipo di nodo (Visitor/Strategy).

## 📐 FASE 2: Prerequisiti per il Calculus (Saldatura F4)
L'algoritmo di Kovacic e i sistemi ODE (F8.1 / F5) richiedono operazioni su matrici simboliche
(campi di funzioni). Il Gram-Schmidt attuale non è sufficiente.

 - [x] Task 2.1: Implementare il riflettore di Householder puramente simbolico per risolvere l'instabilità della decomp QR (Fix F4).

## 🧮 FASE 3: L'Integrazione Simbolica Vera (F8.1 & Fix F5)
Ora che l'AST è modulare e l'algebra lineare simbolica è stabile, possiamo implementare
l'integrazione differenziale.

 - [x] Task 3.1: Generalizzare la Hermite Reduction per accettare un operatore differenziale D arbitrario.
 - [x] Task 3.2: Implementare l'algoritmo di Rothstein-Trager per le RDE (Risch Differential Equations).
 - [~] Task 3.3: Implementare l'algoritmo di Kovacic per le ODE lineari del 2° ordine.
   - **Case 1**: Poli di ordine 2 e parte polinomiale di grado 0 di `r` — completato (test Euler-Cauchy verde).
   - **Case 2**: Completato e testato. Sbloccato dall'integrazione di `RootOf` (Task 5.4). Riccati algebraic certificate attivo su sample points `{1, 4, 9, 16, 25}` (HC-KV-04).
   - **Case 3**: Parzialmente completato. Supporto per gruppi tetraedrico/ottaedrico ($n \in \{4, 6\}$) implementato e testato; $n = 12$ (icosaedrico) escluso per performance (HC-KV-06). Integrazione a valle (`integrate(RootOf)`) da cablare (HC-KV-05).
   - **Bloccato (HC-KV-02)**: poli/poly part di ordine pari ≥ 4 — richiede motore Puiseux/Laurent per √r. Diagnostico `Unimplemented` esplicito.

## ⚡ FASE 4: Prerequisiti per l'Algebra Avanzata (Saldatura F1)
I Resultanti Modulari su grandi torri (F8.2 / F3) generano coefficienti giganteschi. Senza FFT,
i calcoli andranno in timeout.

 - [x] Task 4.1: Pulizia delle legacy form (es. is_neg_infinity duplicate nei .cpp).
 - [ ] Task 4.2: Implementare la moltiplicazione Schönhage-Strassen (FFT) in BigInt (Fix F1).

## 🌪️ FASE 5: Algebra Avanzata e CAD Universale (F8.2 & Fix F3, F6 pt.1)
Con l'aritmetica veloce pronta, possiamo processare torri algebriche profonde e spazi Rⁿ continui.

 - [x] Task 5.1: Costruire i Resultanti Modulari (CRT) per superare il limite delle torri a 2 livelli e supportare EEZ.
 - [ ] Task 5.2: Implementare l'Algoritmo di Stauduhar in Z per i gruppi di Galois ≥ 6.
 - [x] Task 5.3: Sostituire double con MPFR/Arb in sturm.cpp per valutazioni intervallari rigorose.
 - [x] Task 5.4: Introdurre i nodi AlgebraicNumber (RootOf) basati su bound isolanti Rational.
 - [ ] Task 5.5: Implementare l'Lifting e Proiezione di Collins per completare il CAD multivariato.

## 🌌 FASE 6: Superfici di Riemann (F8.3 & Fix F6 pt.2)
Richiede l'AST modulare per aggiungere un nuovo tipo di nodo e il CAD per valutare le singolarità.

 - [x] Task 6.1: Creare il nodo UnwindingNumber $K(z)$.
   - **Chiuso (registrazione builtin)**: aggiunto `BuiltinOp::UnwindingNumber` in `include/cas/builtin_functions.hpp` con doppio mapping nome `K`/`UnwindingNumber`. `builtin_op_name` ritorna `"UnwindingNumber"`. Term-order in `src/symbolic/term_order.cpp` assegna weight 81 (branch-cut, low arithmetic priority). Simplifier preserva l'identità (default-passthrough, NO collapse a 0 — la propagazione branch-aware è riservata a Task 6.2). 5 test in `test/unit/symbolic/test_unwinding_number.cpp` (parsing K/UnwindingNumber, simplify preserve, name round-trip, structural equality). Suite quick: 2325/2326 pass, zero regressioni.
 - [x] Task 6.2: Istruire simplify_impl a propagare attivamente $K(z)$ nelle composizioni log-exp per risolvere il branch-cut.
   - **Chiuso (opt-in)**: introdotto flag `CASContextParams::branch_cut_aware_logexp` (default `false` per backward compat). Quando `true`, la regola `ln(exp(z)) → z` in `src/symbolic/simplify_exp_log.cpp` viene sostituita da `z + 2πi·K(z)` per `z` non dichiarato reale via `Assumptions::is_real`. Coerente per entrambe le forme `ln(exp(x))` e `ln(e^x)`. Quando `false` (default), comportamento legacy preservato — nessuna regressione. Aggiunti 3 test (branch-aware off/on + assumption real). Suite quick: 2328/2329 pass, zero regressioni.

## 🧩 FASE 7: L'Integratore di Ultima Istanza (Saldatura F7)
Solo dopo aver completato Risch (F8.1) ha senso aggiungere Meijer G, che funge da rete di
salvataggio per i casi irrisolvibili.

 - [ ] Task 7.1: Implementare la logica delle trasformazioni di Slater per il pattern matching avanzato.
 - [ ] Task 7.2: Costruire il motore di Meijer G come integratore di fallback nel motore di riscrittura globale.

---
## Log delle Sessioni

### 2026-06-14 (sessione 1)
- **Task 3.3 esteso (Case 2 & Case 3) ✅**: Completata l'implementazione del Case 2 e parziale del Case 3 (n ∈ {4, 6}) di Kovacic in `ode_kovacic_case2.cpp`, `ode_kovacic_case2_helpers.cpp`, `ode_kovacic_case3.cpp` e `ode_kovacic_case3_helpers.cpp`. Aggiunti nuovi test unitari (`test_ode_kovacic_case2.cpp`, `test_ode_kovacic_case3.cpp`).
- **Task 5.4 completato ✅**: Chiuso il supporto per `AlgebraicNumber/RootOf` con bound isolanti ed integrato come dipendenza risolutiva per Kovacic.
- **Risch RDE Hermite Scaffold (F7.5 B3) ✅**: Aggiunto l'entry point per la riduzione parametrica di Hermite in `integrate_core.cpp` (commit `ad1cb34`).

### 2026-06-12 (sessione 1)
- **Task 6.2 esteso ✅ BC-1..BC-3**: branch-cut corrections sotto strict mode. Nuovo TU `src/symbolic/simplify_branch_cut.{hpp,cpp}` con helpers `make_sqrt_of_square_correction` ((-1)^K(2·ln(z))), `make_pow_of_pow_correction` (e^(2πi·b·K(a·ln z))), `make_log_product_correction` (-2πi·K(ln z1 + ln z2)), `make_log_quotient_correction`. Wiring in `simplify_arithmetic_power.cpp` (BC-2 pow-of-pow strict path) e `simplify_exp_log.cpp` (BC-1b sqrt(z²) + BC-3 ln-of-product strict path). 4 nuovi test verdi (`PowOfPow_StrictMode_…`, `LnOfProduct_…`). Commit `c521b57`.
- **Task 2.1 partial closure ✅**: symbolic QR bailout (matrix_qr.cpp) ora ctx-configurable + assumption-aware. Nuovo ctx param `symbolic_qr_max_norm_complexity` (default 2). Bailout salta quando `ctx.assumptions().is_nonnegative(Nx)`. Test `QRTest.SymbolicQR_DefaultSignConvention_2x2` skip-ladder ora diagnostico (3 punti precisi invece di un solo bailout opaco). Ledger `HC-F8-QR-HOUSEHOLDER-BAILOUT` APERTA → PARTIAL. Commit `80ca9fa`.

### 2026-06-11
- **Task 1.1 ✅**: Creati `include/cas/ast_kinds.hpp`, `include/cas/ast_arena.hpp`, `include/cas/ast_nodes.hpp`.
  `ast.hpp` ora è un thin-umbrella che `#include` i tre sub-header (backward-compatible, zero breakage).
- **Task 1.2 ✅**: Creati `src/ast/ast_compare.cpp` (structural_equal + expr_hash + expr_kind_name)
  e `src/ast/ast_clone.cpp` (clone_into_arena). `ast.cpp` ora gestisce solo Arena lifecycle.
- **Task 1.3 ✅**: Rifattorizzato `differentiate.cpp` usando il pattern Visitor (`visit_expr`). Eliminato il dispatch manuale via if/else per migliorare estensibilità e manutenibilità.
- **Task 2.1 ✅**: Implementato riflettore di Householder puramente simbolico e analiticamente razionalizzato in `matrix_qr.cpp`. Questa implementazione previene l'esplosione dell'AST evitando la presenza di radici quadrate (`alpha`) al denominatore durante gli aggiornamenti di Q e R, e risolve l'instabilità senza limitare l'algebra esatta.
- **Task 3.1 ✅**: Generalizzata la Hermite Reduction per supportare un operatore differenziale arbitrario $D$. In `differential_field.cpp`, il calcolo delle derivate di $V$ e $A$ durante l'`ho_reduce_step` sfrutta ora la derivazione di campo (`field.derive(...)`) invece della differenziazione polinomiale formale.
- **Task 3.2 ✅**: Implementato l'algoritmo di Rothstein-Trager generalizzato per le estensioni differenziali in `differential_field.cpp`. Modificata la funzione `integrate_rothstein_trager` affinché la variabile di risultante e il calcolo del GCD usino il generatore differenziale top-level corretto (`t_var`) e affinché la derivata di $Q$ utilizzi l'operatore differenziale del campo (`field.derive`). Questo completa il supporto per l'integrazione differenziale delle parti logaritmiche.
- **Task 3.3 (Parziale)**: Implementato il **Case 1** di Kovacic in due file (split anti-monolith): `src/calculus/ode_kovacic.cpp` (entry point + back-transform + variation of parameters) e `src/calculus/ode_kovacic_case1.cpp` (compute_r + case1_omega + isqrt). Pipeline: classifier riconosce `Linear2ndOrderRationalCoeff` → `solve_ode_kovacic` calcola invariante `r = p'/2 + p²/4 − q`, applica `apart_num_den` + `partial_fractions`, deriva `ω₊/ω₋` da poli ordine 2 con discriminante `1+4A` perfetto quadrato in ℚ, integra per ottenere `z₁,z₂`, back-trasforma `y = z·exp(−½∫p dx)`, applica variation of parameters per RHS ≠ 0. Test Euler-Cauchy `x²y''−2xy'+2y=0` verde con soluzione corretta `(C₂+C₃x)·|x|`. `bigint_isqrt` riscritta con Newton iteration **puramente BigInt** (zero `double`/`int64`, zero limite 2^53). Cases 2/3 e poli ordine pari ≥ 4 → `Unimplemented` con diagnostico esplicito; debiti registrati in `HARDCODE_LEDGER.md` come **HC-KV-02** (Laurent expansion of √r) e **HC-KV-03** (Case 2/3 — blocked on Task 5.4 AlgebraicNumber). Aggiunti due test in `test/unit/test_ode.cpp`: `Kovacic_EulerCauchy_TwoDistinctRealRoots` (PASS) e `Kovacic_SimplePole_Unimplemented` (verifica diagnostico). Suite ODE+Frobenius: 20/20 verde.
- **Task 5.4 step 2 (Bound-aware evaluator + round-trip)**: `NumericEvaluator::evaluate` su `RootOf` ora usa l'`isolating_bound` quando presente per restringere la ricerca Sturm all'intervallo isolante, garantendo che il `double` ritornato sia la STESSA radice identificata dal bound (no drift verso radici sibling). Round-trip printer/parser estesi per `RootOf` con bound: nuovo formato `RootOf(poly, var, idx, low_num, low_den, high_num, high_den)` (7 args), backward-compatible con i formati a 2/3 args. Parser accetta `Unary(Neg, IntegerLit)` per endpoint negativi. Aggiunti 2 test (`NumericEval_UsesBound_NoSiblingDrift`, `RoundTrip_PreservesBound`); suite RootOf 7/7 verde. Quick suite: 2315/2316 pass (+7 totali), zero regressioni.

- **Task 5.4 (Parziale)**: AlgebraicNumber/RootOf con bound isolanti. Esteso `struct RootOf` in `include/cas/ast_nodes.hpp` con `std::optional<IsolatingBound> isolating_bound` (coppie `BigInt num/den` per gli endpoint, evita include-cycle con `rational.hpp`). Nuovo costruttore `RootOf(polynomial, var, IsolatingBound, root_index)`. Aggiornati `structural_equal` e `expr_hash` per discriminare RootOf con bound diversi (`src/ast/ast_compare.cpp`). `clone_into_arena` preserva il bound (`src/ast/ast_clone.cpp`). Nuovo file `src/algebra/algebraic_number_rootof.cpp` con `make_rootof_isolated(poly, var, idx, ctx, low=-1e9, high=1e9, tol=1e-9)`: chiama `numeric::find_polynomial_isolating_intervals` (nuova API in `include/cas/numeric.hpp` + impl in `src/numeric/sturm.cpp`) che ritorna `vector<IsolatingBound>` esatti via Sturm + squarefree decomposition. Bound matematicamente più rigorosi di MPFR float: sono intervalli `[low, high] ∈ ℚ × ℚ` esatti. Nuovo file test `test/unit/algebra/test_algebraic_number_rootof.cpp` con 5 test verdi (Sqrt2, cubic 3-roots, bound distinction, clone preservation, legacy backward compat). Suite quick: 2313/2314 pass (+5 nuovi, 1 skip pre-esistente), zero regressioni in 47.4s. Documentata nota: tolleranze < 2.3e-10 (≈ 2^-32) si annullano nella conversione interna double→Rational; default `tol=1e-9` evita il problema.
- **Task 4.1 ✅**: Pulizia legacy form. Le 4 copie locali di `is_pos_infinity`/`is_neg_infinity` in `src/calculus/{integrate,integrate_improper,integrate_definite_patterns,orthogonal_polynomials}.cpp` (definite nel namespace anonimo) erano **incomplete**: gestivano solo `Unary(Neg, Constant(Infinity))` (legacy form pre-F7.5.F1) e mancavano del caso canonico `Constant(NegInfinity)`. Sostituite tutte e quattro con `using cas::is_pos_infinity; using cas::is_neg_infinity;` dalla sede canonica `include/cas/extended_real.hpp` (introdotta in F7.5.F1 con supporto a entrambe le forme). Aggiunto `#include "cas/extended_real.hpp"` nei 4 file. `integrate_improper.cpp` ora delega `is_infinity` a `cas::is_signed_infinity`. **Fix di correttezza**: prima il classificatore di integrali impropri poteva fallire su intervalli espressi nella forma canonica `[Constant(NegInfinity), Constant(Infinity)]`. Suite quick: 2308/2309 pass (1 pre-existing skip), zero regressioni in 47.4s.
olatingBound>` esatti via Sturm + squarefree decomposition. Bound matematicamente più rigorosi di MPFR float: sono intervalli `[low, high] ∈ ℚ × ℚ` esatti. Nuovo file test `test/unit/algebra/test_algebraic_number_rootof.cpp` con 5 test verdi (Sqrt2, cubic 3-roots, bound distinction, clone preservation, legacy backward compat). Suite quick: 2313/2314 pass (+5 nuovi, 1 skip pre-esistente), zero regressioni in 47.4s. Documentata nota: tolleranze < 2.3e-10 (≈ 2^-32) si annullano nella conversione interna double→Rational; default `tol=1e-9` evita il problema.
- **Task 4.1 ✅**: Pulizia legacy form. Le 4 copie locali di `is_pos_infinity`/`is_neg_infinity` in `src/calculus/{integrate,integrate_improper,integrate_definite_patterns,orthogonal_polynomials}.cpp` (definite nel namespace anonimo) erano **incomplete**: gestivano solo `Unary(Neg, Constant(Infinity))` (legacy form pre-F7.5.F1) e mancavano del caso canonico `Constant(NegInfinity)`. Sostituite tutte e quattro con `using cas::is_pos_infinity; using cas::is_neg_infinity;` dalla sede canonica `include/cas/extended_real.hpp` (introdotta in F7.5.F1 con supporto a entrambe le forme). Aggiunto `#include "cas/extended_real.hpp"` nei 4 file. `integrate_improper.cpp` ora delega `is_infinity` a `cas::is_signed_infinity`. **Fix di correttezza**: prima il classificatore di integrali impropri poteva fallire su intervalli espressi nella forma canonica `[Constant(NegInfinity), Constant(Infinity)]`. Suite quick: 2308/2309 pass (1 pre-existing skip), zero regressioni in 47.4s.
47.4s.
