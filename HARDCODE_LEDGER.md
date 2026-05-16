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

_(nessuna — ledger ripulito 2026-05-16)_

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

## Note operative

- **Cadenza revisione**: ad ogni nuova sessione, leggere questo file per primo
  e valutare se almeno una voce può essere chiusa nel ciclo di lavoro corrente.
- **Politica zero crescita**: prima di aprire una voce N+1 in questo ledger,
  verificare se almeno una voce esistente può essere chiusa nello stesso turno.
- **Tag commit**: ogni commit che apre un HC scrive `(introduce HC-NNN)` nel
  body; ogni commit che chiude un HC scrive `(closes HC-NNN)`.
