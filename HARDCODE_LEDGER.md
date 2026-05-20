# HARDCODE LEDGER — Registro Hardcode di Passaggio

> **Scopo**: tracciare ogni hardcode-of-passage introdotto nel codice in modo che
> nessuno rimanga dimenticato. Ogni voce deve essere risolta prima della prossima
> release. Una voce risolta passa allo storico in fondo al file.
>
> **Protocollo (cf. CLAUDE.md REGOLA ZERO)**:
> 1. Quando un hardcode-of-passage è inevitabile (prerequisito mancante,
>    infrastruttura non pronta), aggiungere `// HARDCODE-OF-PASSAGE: <id>`
>    sopra la riga incriminata, dove `<id>` è la chiave usata qui sotto.
> 2. Aprire una voce in questo ledger con tutti i campi obbligatori.
> 3. Il commit message DEVE menzionare l'id dell'hardcode introdotto.
> 4. Al momento del fix, spostare la voce nella sezione "Storico (risolti)"
>    aggiungendo data + commit di risoluzione.

---

## Voci aperte

### KNOWN-DEBT-001 — `-Werror` disabilitato (CMakeLists.txt:21)
- **File**: `CMakeLists.txt:21` (non-MSVC branch).
- **Categoria**: build hygiene; CLAUDE.md § REGOLA ZERO richiede zero
  warning compilation.
- **Stato**: ~6 pre-existing warnings impediscono restore immediato:
  - `src/rewrite/builtin_rewrite.cpp:179` — unused-function
  - `src/algebra/polynomial_cyclotomic.cpp:11,101` — unused-function /
    -parameter
  - `src/algebra/polynomial_conversions.cpp:170` — non-exhaustive switch
  - `libcas_symbolic.a(functions_special.cpp.o)` — empty TU
- **Fix**: pulire warnings + restore `-Werror`. Effort stimato ~1h.

### KNOWN-DEBT-002 — Test coverage ratio 0.64 (89 test / 139 src)
- **Categoria**: testing.
- **Stato**: Moduli grandi senza unit test dedicato:
  - `src/algebra/factorization_polynomials.cpp` (40.8K)
  - `src/algebra/polynomial_gcd_multivariate.cpp` (34.4K)
  - `src/calculus/integrate_risch.cpp` (36.8K)
  - `src/calculus/differentiate.cpp` (24.1K)
  - `src/calculus/limit.cpp` (25.1K)
- **Fix**: stub happy-path test per ciascuno. Effort ~3h totale.

### KNOWN-DEBT-004 — Risch `∫ 1/(x·ln(x)) dx` produces wrong result — RISOLTO 2026-05-20
- **File**: `src/calculus/integrate_risch.cpp` (path that handles
  `1 / (x * ln(x))`).
- **Categoria**: math correctness — silent wrong answer.
- **Discovered**: 2026-05-20 audit probe
  `test/unit/calculus/test_risch_logarithmic_probe.cpp`.
- **Symptom (pre-fix)**: input `1/(x·ln(x))` → engine returned
  `ln(x)^(-1) · ln(abs(x))`, wrong (not the antiderivative).
- **Fix**: logarithmic-derivative recognizer (Risch structure
  theorem) added at the entry of `integrate_risch` BEFORE the
  Hermite/Rothstein-Trager pipeline. For each extension generator
  `t` (log or exp) the code builds candidate `F = ln(g)` where
  `g = ln(u)` for log ext or `g = exp(u)` for exp ext, then
  searches a small set of trial constants `c ∈ {±1, ±1/2, ±2}`
  for which `c · D(F) ≡ integrand` holds after `together()` +
  `simplify()`. The match returns `c · F` as the closed form.
  Verification by structural differentiation guarantees the
  result is an exact antiderivative — no shape inspection, no
  fragile pattern matching.
- **Side fix**: IBP (`src/calculus/integrate_parts.cpp`) now
  verifies its own result via D(F)=integrand using `together()`
  normalization, rejecting partial/cyclic returns that previously
  leaked the wrong `ln(x)^-1·ln|x|` form.
- **Probe**: `test_risch_logarithmic_probe.cpp::IntegralOfReciprocalOfXLnX`
  asserts the diff-inverse invariant holds.

### KNOWN-DEBT-003 — Test DISABLED senza task aperto
- `test_residue_theorem.cpp:138`: double-pole Q(α) GTEST_SKIP
- `test_equivalence_subset.cpp:105`: DISABLED_ExpOfLogSumWithoutPositivity...
- `test_factorization_tower.cpp:202,217`: DISABLED_SplitsProductOfQuadratics...
- `test_factorization_trager.cpp:281,291,310,328`: 4× GTEST_SKIP Galois
- **Fix**: collegare ciascuno a task aperto in CAS_TASKS.md.

---

## Storico (risolti)

### HC-001 — Bridge depth limit — RISOLTO 2026-05-16
- Fix: `ctx.max_q_alpha_bridge_depth()` esposto in `CASContext`
  (default 256, min 8). `express_recursive` in
  `src/algebra/algebraic_number_bridge.cpp` consuma il valore dal contesto.
  Costante `kMaxBridgeDepth` rimossa.

### HC-002 — Gamma half-integer recursion bound — RISOLTO 2026-05-16
- Fix: `ctx.max_gamma_recursion()` esposto (default 1024, min 16). Blocco
  `Gamma half-integer` in `src/symbolic/simplify_functions.cpp` usa il
  valore via `context_->max_gamma_recursion()`.

### HC-004 — Frobenius/Trager/MRV fresh symbols — RISOLTO 2026-05-16
- Fix: `CASContext::make_fresh_symbol(prefix)` implementato in
  `src/symbolic/context_core.cpp`. Counter monotono +
  probe contro `variables_` garantiscono unicità rispetto a tutte le
  chiamate precedenti e ai simboli definiti dall'utente.
- Siti aggiornati:
  - `src/calculus/ode_solver_frobenius.cpp`: `_C1_`/`_C2_` → fresh symbol
    con prefisso `C`.
  - `src/calculus/differential_field.cpp`: Trager `_rt_t_` → fresh symbol
    `rt_t_<n>`.
  - `src/calculus/limit_mrv.cpp`: MRV `__mrv_w` → fresh symbol `mrv_w_<n>`
    con re-roll se collisione strutturale.
- Test anti-collisione: `MakeFreshSymbol_UniqueAndAvoidsUserScope` in
  `test/unit/symbolic/test_special_functions.cpp`.

### HC-005 — Residue-theorem polynomial degree heuristic — RISOLTO 2026-05-16
- Fix: `poly_degree_rational` in `src/calculus/residue_theorem.cpp` ora usa
  `algebra::parse_polynomial(expr, var, ctx)` direttamente (header privato
  `../algebra/polynomial_internal.hpp` incluso). Loop Taylor + early-exit
  "8 trailing zeros" rimosso.

### HC-006 — Improper convergence scan window — RISOLTO 2026-05-16
- Fix: `ctx.improper_leading_order_scan()` esposto (default 8, min 1).
  `effective_leading_order` in
  `src/calculus/integrate_improper.cpp` legge il valore dal contesto
  quando il parametro chiamante è 0.

### HC-007 — Mignotte GCDHEU `B<1000` override — RISOLTO 2026-05-19
- **File**: `src/algebra/polynomial_gcd_heuristic.cpp:268` (era).
- **Categoria CLAUDE.md**: Cat 2 — costante magica in algoritmo algebrico.
- **Fix applicato (commit `cce829b`)**: rimpiazzato `if (B<1000) B=1000` con
  bound Mignotte rigoroso `B = 2·(max_coeff+1)·2^(D+1)+1`, dominante per
  ogni (D, max_coeff). Floor strutturale `B≥3` per signed-digit
  reconstruction. Test `KroneckerAtRigorousMignotteBoundReturnsTrueGcd`
  aggiornato (era basato su B=1000 spurious).

### HC-008 — LPO precedence non esaustiva — RISOLTO 2026-05-19
- **File**: `src/symbolic/term_order.cpp:68-85`.
- **Categoria CLAUDE.md**: Cat 5 — ordinamento non configurabile,
  default=50 collisione su 40+ BuiltinOp.
- **Fix applicato (commit `ed22283`)**: switch esaustivo su tutti i
  BuiltinOp con precedenze distinte derivate dai vincoli rewrite rules
  R1-R12 (exp(ln x)→x, sqrt(x²)→|x|, tan→sin/cos, ecc.). Default branch
  ritorna 0 per nuovi builtin senza UB.

### HC-009 — Bareiss pivot score magici — RISOLTO 2026-05-19
- **File**: `src/linalg/matrix_bareiss.cpp:110-118`.
- **Categoria CLAUDE.md**: Cat 2 — costanti `1000/900/800/500-cplx`.
- **Fix applicato (commit `2892492`)**: `PivotScore` lessicografico
  `(certainty, -total_degree, -complexity)`. `certainty` ∈ {0..3}
  derivato da assumptions; `total_degree` recursive AST walk;
  `complexity` AST-size tiebreaker.

### HC-010 — Limit depth=16 fisso — RISOLTO 2026-05-19
- **File**: `src/calculus/limit.cpp:118` (era `if (depth >= 16U)`).
- **Categoria CLAUDE.md**: Cat 1 — budget computazionale non
  configurabile.
- **Fix applicato (commit `a8d3e75`)**: bound dinamico
  `max(8, 2·tower_height + 4)` con `tower_height` calcolato via
  `transcendental_tower_depth()` helper (Gruntz §3.5). Cap scala con
  altezza tower MRV.

### HC-011 — Buchberger kMaxBuchbergerPairs+kMaxBasisSize — RISOLTO 2026-05-19
- **File**: `src/algebra/polynomial_groebner_f4_buchberger.cpp:159-160`
  (era).
- **Categoria CLAUDE.md**: Cat 1 — budget non configurabile.
- **Fix applicato (commit `43dd1fb`)**: rimossi entrambi i cap. Sugar
  selection strategy (Giovini-Mora-Niesi-Robbiano 1991) + Gebauer-Moeller
  pruning garantiscono basis minimale per Hilbert basis theorem. `Pair`
  struct esteso con `sugar` field; `basis_sugar` parallel vector;
  `select_pair` ora lex `(sugar, total_lcm_degree)`.

### HC-012 — fsolve num_samples=400 (poly path) — RISOLTO 2026-05-19
- **File**: `src/algebra/fsolve.cpp:76` (era).
- **Categoria CLAUDE.md**: Cat 2 — costante magica.
- **Fix applicato (commit `1108880`)**: Sturm sequence 1829 in nuovo
  `src/numeric/sturm.cpp`. Conta esattamente radici reali distinte in
  intervallo; squarefree via `gcd(f, f')`; bisection rationale + Newton
  polish. Path polinomiale ora esatto, grid scan resta solo per
  transcendental (sostituito in HC-014).

### HC-013 — F4 kMaxF4Batches=2048 — RISOLTO 2026-05-19
- **File**: `src/algebra/polynomial_groebner_f4.cpp:200` (era).
- **Categoria CLAUDE.md**: Cat 1 — batch counter arbitrario.
- **Fix applicato (commit `fb23498`)**: rimosso. Hilbert basis theorem
  + Buchberger termination + Sugar selection garantiscono P.empty() in
  tempo finito. Memory guards `kMaxMacaulay*` mantenuti (vedi HC-016).

### HC-014 — fsolve tol/max_iter (transc path) — RISOLTO 2026-05-19
- **File**: `src/algebra/fsolve.cpp:78-79` (era).
- **Categoria CLAUDE.md**: Cat 2 — costanti precision/iter magiche.
- **Fix applicato (commit `5f5e068`)**: Lipschitz dyadic refinement
  (Hansen-style) in nuovo `src/numeric/lipschitz.cpp`. Esclusione
  intervalli senza root via stima Lipschitz 3-punto; descent dyadic;
  Newton polish (30 iter da quadratic convergence). `sin(50x)` su [0,1]
  trova ~16 root, era miss-dipendenti su grid.

### HC-015 — Recombination kMaxSubsets=32768 — RISOLTO 2026-05-19
- **File**: `src/algebra/factorization_recombination.cpp:105` (era).
- **Categoria CLAUDE.md**: Cat 1 — cap subset enumeration.
- **Fix applicato (commit `6c809cd`)**: rimosso. Landau-Mignotte
  coefficient bound pruning (TAOCP §4.6.2 Thm F): `||h||_inf ≤ 2^d·||f||_inf`
  per ogni fattore h. `lift_and_check_subset` rifiuta candidati che
  eccedono il bound prima della divisione costosa. Pruning polynomial
  in pratica; pathological Swinnerton-Dyer resta esponenziale → Step 8.

### HC-016 — F4 Macaulay caps configurable — RISOLTO 2026-05-19
- **File**: `src/algebra/polynomial_groebner_f4.cpp:201-203` (era
  `constexpr`).
- **Categoria CLAUDE.md**: Cat 1 — memory safety hardcoded.
- **Fix applicato (commit precedente in questa sessione)**:
  `ctx.f4_max_macaulay_rows()`, `f4_max_macaulay_monomials()`,
  `f4_max_pending_monomials()` esposti in `CASContext`. Default
  preserva valori storici (512/512/1024). Esposizione, non upgrade
  algoritmico — Block F4 vero (Faugère 2002 §4.3) resta follow-up.

### HC-003 — Zeta closed-form lookup table — RISOLTO 2026-05-16

### HC-003 — Zeta closed-form lookup table — RISOLTO 2026-05-16

- **File**: `src/symbolic/simplify_functions.cpp` (zeta block) + nuovo
  `src/numtheory/bernoulli.cpp` + `include/cas/numtheory.hpp`.
- **Categoria CLAUDE.md**: Cat 3 — REGOLA ZERO violation.
- **Fix applicato**:
  1. Esposto `cas::numtheory::bernoulli_numbers(unsigned)` e
     `bernoulli_number(unsigned)` come API pubblica in
     `include/cas/numtheory.hpp`. Implementazione spostata in
     `src/numtheory/bernoulli.cpp` (Akiyama–Tanigawa).
  2. `summation.cpp` ora consuma `cas::numtheory::bernoulli_numbers`
     (rimossa la funzione `static` locale).
  3. `cas_symbolic` linka `cas_numtheory` (CMake aggiornato).
  4. `simplify_functions.cpp` zeta block: lookup `{2,4,6,8,10,12}` rimpiazzato
     con formula generale `ζ(2k) = (-1)^(k+1) · 2^(2k-1) · π^(2k) · B_{2k}/(2k)!`
     per qualsiasi `2k`. Negativi dispari `ζ(-(2k-1)) = -B_{2k}/(2k)` idem.
  5. Test anti-hardcode aggiunti in
     `test/unit/symbolic/test_special_functions.cpp`:
     `ZetaFourteenViaBernoulli_AntiHardcode`,
     `ZetaSixteenViaBernoulli_AntiHardcode`,
     `ZetaNegativeNineViaBernoulli_AntiHardcode`,
     `ZetaNegativeElevenViaBernoulli_AntiHardcode`. Tutti i test zeta
     esistenti continuano a passare (10/10 verde).

---

## Algorithm Fidelity Audit (2026-05-19)

Onestà rispetto a citazioni letterarie nei commit message:

### Implementati al massimo (5)
- **Sturm 1829** (HC-012, commit `1108880`): sequence + V(a)-V(b) +
  bisection rationale + Newton polish + squarefree decomp completi.
- **Knuth TAOCP §4.6.2 Mignotte** (HC-007, commit `cce829b` e
  HC-015, commit `6c809cd`): formula applicata correttamente.
- **GMNR Sugar 1991 variante inhomogena** (HC-011, commit `43dd1fb`):
  formula S-poly sugar + selezione lex implementate. Variante
  omogena richiederebbe `homogenize()` infrastructure (assente).
- **Hilbert basis theorem** (HC-013, commit `fb23498`): invocazione
  corretta come giustificazione, niente algoritmo da implementare.
- **Buchberger 1985 on-fly tail-reduction** (commit successivo a
  questa audit): variante sicura aggiunta ex-post. NB: piena
  minimization durante run NON sicura — minimization vera resta
  post-loop in `inter_reduce`.

### Implementati semplificati / approssimati (3)
- **Gruntz 1996 §3.5** (HC-010, commit `a8d3e75`): mio contributo =
  helper `transcendental_tower_depth` + bound adattivo. L'algoritmo
  Gruntz MRV vero PREESISTEVA in `limit_mrv.cpp`. Ho sbloccato la
  profondità di ricorsione, non implementato Gruntz da zero.
- **Hansen 1992 interval Newton** (HC-014, commit `5f5e068`): ho
  implementato float-Lipschitz refinement ("Hansen-style"), NON
  interval arithmetic rigoroso. Cap Lipschitz da 3-point sample
  estimate. Vero Hansen richiede MPFR interval lib (futuro L3-01).
- **Lecerf 2007 §3 pruning** (HC-015, commit `6c809cd`): Mignotte
  bound pruning implementato. Variante polynomial-time Lecerf con
  LLL NON implementata.

### Nominali (3) — citati nei plan, NON nel codice
- **Faugère 2002 Block F4 §4.3**: documentato FE-001. NON
  implementato. Solo termination cap rimosso (commit `fb23498`).
- **van Hoeij 2002 Theorem 4.2**: documentato FE-002. Tentativo
  BFS-by-size 2026-05-19 → regressione perf → rolled back.
- **Belabas-vH-Klüners-Steel 2004**: citato nei piani, NON nel codice.

---

## Future Enhancements (NON debts — algoritmo corretto, solo performance)

Queste voci NON sono hardcode né debiti. La correttezza è garantita su
tutti gli input; questi sono upgrade prestazionali per casi specifici.

### FE-001 — Block F4 (Faugère 2002 §4.3)
- Stato: caps memoria F4 ora configurabili via `ctx.f4_max_macaulay_*`
  (HC-016 chiuso 2026-05-19). Block decomposition vero — partition
  monomials in LM-divisibility cones, run Gauss per blocco, merge —
  abiliterebbe input n=6+ var, deg=6+ senza fallback a Buchberger.
- Effort stimato: ~12h (partition helper + block dispatch).

### FE-002 — van Hoeij knapsack lattice (van Hoeij 2002 Thm 4.2)
- Stato: subset enumeration ora con Mignotte pruning (HC-015 chiuso
  2026-05-19). Pathological Swinnerton-Dyer r≥25 resta esponenziale.
- Trace-polynomial lattice via LLL replacerebbe enumeration con
  ricerca polinomiale O(r⁴·N·log p^a). Infrastruttura LLL già presente
  in `src/algebra/lattice_lll.cpp`; servono trace polys via Newton's
  identities + lattice construction.
- Effort stimato: ~14h. Tentativo BFS-by-size 2026-05-19 ha mostrato
  regressione perf su input tipici (6s→28s) → richiede approccio
  lattice puro, non solo riordino enumeration.

### FE-003 — `MAX_BIGINT_LIMBS=10000` — ACCETTATO permanentemente
- CLAUDE.md REGOLA ZERO Eccezione 4: hardware OOM safety.
- Produce `Unimplemented` esplicito a oltrepassare bound (no silent
  wrong result).
- Non rimuovere senza una vera streaming-arithmetic implementation.

---

## Note operative

- **Cadenza revisione**: ad ogni nuova sessione, leggere questo file per primo
  e valutare se almeno una voce può essere chiusa nel ciclo di lavoro corrente.
- **Politica zero crescita**: prima di aprire una voce N+1 in questo ledger,
  verificare se almeno una voce esistente può essere chiusa nello stesso turno.
- **Tag commit**: ogni commit che apre un HC scrive `(introduce HC-NNN)` nel
  body; ogni commit che chiude un HC scrive `(closes HC-NNN)`.
