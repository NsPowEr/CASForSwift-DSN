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

### HC-001 — Bridge depth limit

- **File**: `src/algebra/algebraic_number_bridge.cpp`
- **Riga**: definizione `kMaxBridgeDepth = 256U`
- **Categoria CLAUDE.md**: Cat 1 (Budget computazionale non configurabile)
- **Descrizione**: profondità massima ricorsione di `express_recursive` durante
  conversione `ExprPtr` → `AlgebraicNumber`. Limite difensivo per alberi adversariali.
- **Fix corretto**: esporre `ctx.max_q_alpha_bridge_depth()` in `CASContext`
  (default 256, configurabile).
- **Blocking dependency**: nessuna (può essere fixato subito).
- **Introdotto**: commit 42c5d62 (2026-05-15)
- **Stato**: aperto

### HC-002 — Gamma half-integer recursion bound

- **File**: `src/symbolic/simplify_functions.cpp`
- **Riga**: blocco `Gamma half-integer`, costante `safety_max = 1024`
- **Categoria CLAUDE.md**: Cat 1 (Budget computazionale non configurabile)
- **Descrizione**: bound iterazioni nella ricorsione `Gamma(z) ↔ Gamma(z±1)` per
  ridurre argomento half-integer arbitrario a `Gamma(1/2)`. Evita loop infinito
  su input malformati.
- **Fix corretto**: esporre `ctx.max_gamma_recursion()` (default 1024).
- **Blocking dependency**: nessuna.
- **Introdotto**: commit 42c5d62 (2026-05-15)
- **Stato**: aperto

### HC-004 — Frobenius integration constants

- **File**: `src/calculus/ode_solver_frobenius.cpp`
- **Riga**: nomi letterali `"_C1_"` e `"_C2_"` per costanti d'integrazione.
- **Categoria CLAUDE.md**: Cat 7 (Nomi variabili interni hardcoded)
- **Descrizione**: la soluzione generale `C₁·y₁(x) + C₂·y₂(x)` usa nomi fissi
  perché manca un generatore di simboli freschi sul `CASContext`. Rischio
  collisione se utente ha `_C1_` o `_C2_` nel proprio scope.
- **Fix corretto**: implementare `CASContext::make_fresh_symbol(std::string prefix)`
  che genera nomi univoci tracciati da counter interno. Aggiornare Frobenius +
  ogni altro sito che oggi usa nomi fissi (vedi anche `_rt_t_` in Trager,
  `__mrv_w` in limit_mrv).
- **Blocking dependency**: introduzione `make_fresh_symbol` (~50 righe in
  `context_core.cpp` + dichiarazione in `symbolic.hpp`). È un gap globale già
  noto e citato in CLAUDE.md Categoria 7.
- **Introdotto**: commit 42c5d62 (2026-05-15)
- **Stato**: aperto — **gap globale**

### HC-005 — Residue-theorem polynomial degree heuristic

- **File**: `src/calculus/residue_theorem.cpp`
- **Riga**: estrazione gradi N/D usa loop con early-exit "8 trailing zeros".
- **Categoria CLAUDE.md**: Cat 1/3
- **Descrizione**: l'agent ha usato uno scan con 8 zeri consecutivi come
  euristica per terminare la rilevazione del grado. Su polinomi con strutture
  particolari può sotto/sovra-stimare.
- **Fix corretto**: utilizzare `algebra::parse_polynomial(expr, var, ctx)`
  pubblica per ottenere direttamente la struttura polinomiale (size = deg + 1).
- **Blocking dependency**: nessuna.
- **Introdotto**: commit 42c5d62 (2026-05-15)
- **Stato**: aperto

### HC-006 — Improper convergence scan window

- **File**: `src/calculus/integrate_improper.cpp`
- **Riga**: `scan_window = 8` in `effective_leading_order`.
- **Categoria CLAUDE.md**: Cat 1 (Budget computazionale non configurabile)
- **Descrizione**: finestra fissa per scansione coefficienti Laurent oltre il
  `leading_order` riportato; serve a riconoscere il "primo coefficiente non
  zero effettivo" quando il numeratore si annulla allo stesso punto del
  denominatore.
- **Fix corretto**: esporre `ctx.improper_leading_order_scan()` con default 8.
- **Blocking dependency**: nessuna.
- **Introdotto**: commit 42c5d62 (2026-05-15)
- **Stato**: aperto

---

## Storico (risolti)

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
