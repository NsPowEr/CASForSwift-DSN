# TASKLIST MASTER — CAS Engine (Single Source of Truth)

> **Unico documento operativo dei task.** Sostituisce e consolida: `CAS_TASKS.md`,
> `PLAN_HP_PRIME_PARITY.md`, `PLAN_F3_F8_GAP_CLOSURE.md`, `PLAN_NEXT_SESSIONS.md`,
> `PLAN_TASKS_REMAINING.md`, `TODO_PH8.md`, `TODO.md`, `HANDOFF_T055.md`, `STATE.md`
> (tutti marcati `⚠️ SUPERSEDED`, non aggiornare).
>
> **Generato 2026-06-26** da audit code-level: ogni voce qui sotto è stata verificata
> **a codice** (file:riga, `git log`, symbol grep), NON copiata dallo stato dichiarato nei
> vecchi file. I vecchi tracker **sovrastimavano** il lavoro aperto e marcavano "RISOLTO"
> codice con gap reali o stub (vedi memoria `tasklist-ledger-stale-2026-06-19`,
> `plan-lies-2026-06-11`). Regola di manutenzione: prima di toccare un task, **ri-verifica
> a codice** — gli status flag di `ledger_index.py` restano inaffidabili.
>
> Registro hardcode autoritativo: `HARDCODE_LEDGER.md` (via `python3 scripts/ledger_index.py`).
> Criteri di chiusura: `DEFINITION_OF_DONE.md`.

---

## Legenda codici `[E·C·S·R]`

- **E**ffort: E1 <1g · E2 1-3g · E3 1-2sett · E4 2-4sett · E5 >1mese
- **C**omplessità: C1 meccanico · C2 moderato · C3 algoritmico · C4 research-grade
- **S**everità: S1 nice-to-have · S2 debito qualità · S3 gap funzionale · S4 bloccante · S5 critico (correttezza/hang)
- **R**ischio regressione: R1 basso (nuovo file) · R2 moderato (file condivisi+guardia) · R3 alto (hot-path: simplifier/arena/BigInt)

Colonna **Refs** = ID nei vecchi tracker (tracciabilità). `git`-anchor = evidenza verifica.

---

## A — TASK APERTI ACTIONABLE (backlog attivo)

Ordinati per severità/impatto decrescente. Ogni voce verificata aperta a codice 2026-06-26.

### A1 · Risch RP-2 Hermite reduction parametrica (df>0) — AUDIT 2026-06-27 (premessa corretta) — `[E4·C4·S3·R2]`
- **AUDIT 2026-06-27** (spec `Risch_Transcendental_Cap8.md` + `Risch_Hermite_Cap5.md` lette): la premessa "scaffold da sostituire" era **fuorviante**. Stato reale:
  - `risch_rde_bronstein_hermite.cpp::risch_rde_hermite_parametric_stub` è **DEAD CODE**: compilato (CMakeLists:220) ma **mai chiamato** da nessuno (grep verificato). NON è il punto attivo.
  - Il solver parametrico **reale** `solve_risch_de_parametric_field` (`risch_rde_bronstein.cpp:249`) è **già implementato**: denominator bound (Q,B,D), trasformazione `f_new=f−D'/D`, `g_new=g·D`, bound di grado N (Log/Exp, §6.4), ricorsione tower `solve_recursive` con correzioni exp/log. Trial constants **già rimossi** (grep vuoto).
  - **GAP REALE A1** = ramo **`df > 0`** a `risch_rde_bronstein.cpp:362-368`: `return Unimplemented("df > 0 not implemented in parametric field solver")`. Tutto `df ≤ 0` funziona; manca il caso *non-cancellation* in cui `deg_t(f) > 0` (il termine `f·y` domina): `deg(y)=dg−df`, `lc(y)=lc(g)/lc(f)`, poi ridurre+ricorrere (Bronstein 6.5/7.4/8.4 PolyRischDE).
- **Next step (definito)**: implementare il ramo df>0 a riga 362 seguendo Bronstein 6.5 (non-cancellation case) nel contesto parametrico (m vettori g simultanei, soluzioni con costanti `c_i`); cross-validare con Maxima. Valutare se rimuovere il dead-code scaffold + relativa riga CMakeLists. Rischio silent-wrong alto → ogni step deve citare il teorema; su indecidibile `Unimplemented` esplicito (mai best-guess).
- **Spec**: `MISSING_FEATURES_SPECS/Risch_Transcendental_Cap8.md`, `Risch_Hermite_Cap5.md`
- **Ledger**: HC-F8-PENDING-17 · **Refs**: T-007, Task 17, F5.1, F7.5.B3, CAS-L1-02

### A2 · Hard-timeout cancellation token su tutti i path integrazione — RISOLTO 2026-06-26 — `[E2·C2·S5·R2]`
- **Stato a codice**: ✅ chiuso. L'ultimo loop combinatorio non interrompibile dell'integrazione (le due eliminazioni di Gauss pure-Rational in `integrate_risch_rde.cpp`) ora polla `check_interrupt()` per colonna pivot; rimosso lo stub morto del poll mai cablato; dispatcher `solve_risch_de_q` reso fedele alla cancellazione (non ingoia più il Timeout di simplify, non maschera l'interrupt come Unimplemented). +2 test in `test_integrate_interrupt.cpp` (RDE solver level). 147 test Risch/integrate verdi.
- **Residuo non bloccante**: fedeltà error-kind su interrupt asincrono *durante* `parse_polynomial` interno a un ramo (rietichettato Unimplemented, finestra ~µs, nessun hang) — vedi ledger HC-F70-A33.
- **Ledger**: HC-F75-A3-HARD-TIMEOUT (chiuso), HC-F70-A33-POLL-COVERAGE (aggiornato 2026-06-26) · **Refs**: T-015, T-024

### A3 · Perf-hang fattorizzazione torri / VanHoeij — `[E3·C4·S2·R2]`
- **Stato a codice**: hang >400-500s su `Q(√2,√3,√5)` factorization e VanHoeij SD3 Swinnerton-Dyer.
- **Ledger**: HC-F8-FACTORIZATIONTOWER-PERF, HC-F8-FACTORIZATIONTOWER-AntiHardcode-X2Minus2-Sqrt3Sqrt5, HC-F8-SD3-VANHOEIJ-SLOW · **Refs**: BUG-HANG-002

### A4 · Smith Normal Form su Q[x] PID generale — ✅ FATTO (verificato 2026-06-26, stale) — `[E3·C3·S3·R2]`
- **Stato a codice**: già implementato in `src/linalg/matrix_smith_qx.cpp` (`smith_normal_form_qx`): algoritmo PID Bezout-based completo (row/col ops con `polynomial_bezout` + `polynomial_exact_divide`, catena divisibilità, normalizzazione monica). Dispatch automatico da `matrix_smith.cpp:163` quando entrate single-var poly. I citati `:166/:216` sono il bail multivar (correttamente NON-PID) e il guard Z-path — non gap. Q[x] è Euclideo ⇒ no Storjohann/LLL necessario. Coverage rafforzata 2026-06-26: +test non-diagonale `[[x,1],[−1,x]]` → `diag(1, x²+1)` con certificato `U·A·V==S`. 4/4 SmithQxTest.
- **Refs**: F4.2, CAS-L2-04

### A5 · ODE Frobenius: log-term + resonance (radici differing-by-integer) — PARZIALE 2026-06-26 — `[E3·C3·S3·R1]`
- **Stato a codice**: ✅ gap principale chiuso: **radice indiciale doppia** (gap N=0) ora produce la seconda soluzione logaritmica via costruzione parametro-derivata `y₂ = ln(x)·y₁ + x^{r1}·Σ aₙ'(r1) xⁿ` (`build_double_root_log_branch`, recurrence con root simbolico ρ + `d/dρ` + sub ρ=r1). Prima `series_solutions.size()==1` restituiva solo `C1·y₁` (soluzione generale incompleta). Verificato: Euler `x²y''−xy'+y=0` → `C1·x + C2·x·ln x`; Bessel₀ `x²y''+xy'+x²y=0` → `J₀·ln x + x²/4 − 3x⁴/128…`. +2 test, 6/6 FrobeniusTest. Gap intero positivo (`build_log_branch`) + resonance detection già presenti.
- **Residuo**: risonanza secondaria multi-livello durante costruzione bₙ (log branch, `ode_solver_frobenius_series.cpp:201`) → `Unimplemented` esplicito diagnostico; punto singolare regolare a x≠0 (solo x=0). **Refs**: F5.3, CAS-L2-01

### A6 · Stauduhar Galois group deg ≥ 6 — `[E4·C4·S3·R1]`
- **Stato a codice**: solo `src/algebra/galois_deg5.cpp` (deg ≤5: S5/A5/D5/F20/C5). Nessun descent Stauduhar.
- **Cosa**: Stauduhar 1973 + tabelle transitive-subgroups Hulpke (n=6..10) + resolventi numeriche BigFloat.
- **Spec**: `MISSING_FEATURES_SPECS/Galois_Groups.md` · **Refs**: T-010, Task 9, F3.6, CAS-L3-18

### A7 · Slater pFq → Meijer G fallback integrator — `[E4·C4·S3·R1]`
- **Stato a codice**: `BuiltinOp::MeijerG` **non esiste**. Nessun file Slater/Meijer.
- **Cosa**: definire `MeijerG`, identità Erdélyi-Slater, trasformazioni Bailey, dispatcher, wiring come passo finale in `integrate_core`. Rete di salvataggio post-Risch.
- **Spec**: `MISSING_FEATURES_SPECS/Special_Fn_Identities.md` · **Refs**: T-013, Task 22-23, F7.A-B

### A8 · Kovacic Case 3 completo (icosaedrico n=12) — `[E3·C4·S1·R1]`
- **Stato a codice**: `ode_kovacic_case3.cpp` ha tetraedrico/ottaedrico; n=12 = solo diagnostico (recurrence blow-up).
- **Cosa**: ottimizzazione recurrence per A₅ n=12. **Ledger**: HC-KV-06 · **Refs**: T-028, T-009 (parziale), F5.C

### A9 · Tower algebrico ≥3 livelli / multi-β nested non-squarefree — `[E3·C4·S2·R2]`
- **Stato a codice**: `algebraic_tower_primitive_nested.cpp` esiste (F3.4-DEBT-01) ma → `Unimplemented` se la min-poll liftata non è squarefree. Bridge solo 2-level (`detect_two_level_tower`).
- **Cosa**: nesting multi-β >1 + factor Q(β)[x] reducible. **Refs**: F3.5, F3.D, Task 7, CAS-L3-06

### A10 · Eigenvalues n>3 / catene di Jordan sotto RootOf — PROB. STALE 2026-06-26 — `[E3·C3·S3·R2]`
- **Stato a codice**: routing RootOf GIÀ presente (`matrix_jordan.cpp:150` → `null_space_over_extension` quando autovalore è `RootOf`). 8/8 test Jordan verdi: `RootOf_Eigenvalues_2x2_Sqrt2`, **`RootOf_Multiplicity2_CompanionDeg4`** (4×4, n>3, catena mult-2 sotto RootOf), `NonDiagonalizable_2x2_Double`. La descrizione "non affidabili" pare obsoleta per i casi testati.
- **Da fare (verifica)**: probe mirato n≥5 con autovalore RootOf molteplicità ≥3 (catena lunghezza ≥3) per confermare/escludere gap residuo prima di chiudere. **Refs**: CAS-L2-02, F4.4

### A11 · Gruntz MRV growth-rank dinamico — 🧊 FROZEN 2026-06-27 (copertura ottima, residuo C4 deferito) — `[E3·C4·S2·R2]`
- **Stato a codice**: rank statico GIÀ rimosso (`limit_infinite.cpp` — confronto via `compare_growth` dinamico ricorsivo, Cat-10 chiuso). Torri log annidate (F7.5.D1+D2) funzionanti: 11/11 (incluso C6 `(log x+log log x)/log x=1`, C7).
- **FATTO 2026-06-26/27** (spec `Gruntz_Nested_Log.md` letta):
  - potenze frazionarie vs log: `log(x)/sqrt(x)→0` (era `∞`, silent-wrong). Cause: predicati `is_positive_power_growth`/`is_reciprocal_positive_power_growth` richiedevano esponente **intero** → ora qualsiasi razionale (classe di crescita = segno dell'esponente) + riconoscimento `sqrt`; `is_logarithmic_in_var` solo `Ln` → ora anche `Log`.
  - **sum-ratio sign** (A11-minor, commit 3067740): `(x+log x)/(x-log x)→1` (era `-1`), `x/(x-log x)→1` (era `0`). Cause: FuncCall handler in `try_infinite_limit` riconosceva solo `Ln` (→ `log(x)` non valutato a `+∞`); mancava riduzione leading-term per quoziente di Sum ∞/∞ → aggiunto `dominant_growth_term` + ratio via `extract_quotient_view`.
  - Test: `LogOverSqrtX`, `FractionalPowerBeatsLog`, `SumRatioLeadingTermSign`. 14/14 gruntz, 75/75 limit/mrv, 180/180 integrate/series/summation/residue.
- **🧊 RESIDUO C4 CONGELATO (deferito, non bloccante)**: track del **coefficiente lentamente variabile** Gruntz §3.5. Sintomo: `x^(3/2)/(x·log x)` → `ComplexRational division by zero` in `leading_power_w` (`limit_mrv_leading.cpp`). Il coefficiente del termine leader non è costante ma funzione del livello MRV successivo (il `log` è "slowly varying"); `leading_power_w` lo tratta come costante e divide per un coefficiente nullo. Fix corretto = ricorsione Gruntz completa sul coefficiente (Gruntz §3.5, `gruntz.py` come riferimento algoritmico). Effort multi-sessione. **NON ridurre il budget né aggiungere hardcode per mascherare**: quando irrisolvibile per ora, `leading_power_w` deve ritornare `Unimplemented` (non valore errato). **Refs**: CAS-L1-01, F5.2, F7.5.D1

### A12 · RootOf come operatore algebrico semplificabile/valutabile — ✅ FATTO (verificato 2026-06-26) — `[E3·C3·S3·R2]`
- **Stato a codice**: RootOf è operatore completo, non più "semi-inerte":
  (1) **simplify** — `simplify_node(RootOf)` risolve esplicitamente se `deg ≤ ctx.max_rootof_explicit_degree()` (simplify_functions.cpp:241);
  (2) **riduzione di potenza** — `Rⁿ → (xⁿ mod P)(R)` via divisione polinomiale (simplify_arithmetic_power.cpp:122), verificato su quintica irriducibile `x⁵−x−1`: `R⁵→R+1`, `R⁶→R²+R`;
  (3) **eval numerica robusta** — Sturm isolation + isolating-bound del nodo (evaluator.cpp:161), con `root_index`;
  (4) **normalizzazione Q(α)** on-demand via `simplify_in_q_alpha` (usata in residue/eigen);
  (5) parser+printer round-trip con isolating bounds.
  +test `RootOfQuinticPowerReduction`. La forma inerte di default per `1/R` è scelta standard (come Maxima/SymPy: ratsimp esplicito); non è un debito. **Refs**: CAS-L1-05

### A13 · Residue theorem grado arbitrario al denominatore — QUASI-FATTO 2026-06-26 — `[E3·C3·S2·R1]`
- **Stato a codice**: ✅ grado arbitrario già coperto: quadratici + biquadratici esatti (Q(α)) **a molteplicità qualsiasi**, fattori grado ≥3/quartiche generali via fallback numerico Aberth (`numeric_residue_contribution`). Poli di ordine >1 risolti dalla ricorrenza di Laurent in `residue()`; l'assembly Q(α)→ℝ è un funzionale lineare indipendente dalla molteplicità. **Rimosso 2026-06-26** il guard `multiplicity > 1` sul ramo biquadratico (era falso limite). +4 test: `1/(x²+1)³`=3π/8, `1/(x⁴+1)²`=3π/2^(5/2), `1/(x⁴+x²+1)²`=2π/3^(3/2) (oracolo Maxima), `1/(x²+1)²`=π/2 (skip morto → assert duro). 14/14 ResidueTheoremTest, 51/51 area residue+laplace.
- **Residuo (deliberato, fuori scope)**: poli reali (fattori lineari) → rifiutati; Cauchy Principal Value + Jordan lemma per integrandi oscillatori (e^{iax}) non implementati. **Refs**: CAS-L2-22, F5.6

### A14 · Equivalenza trascendente `are_equal` (assumption propagation completa) — ✅ FATTO/SOUND (verificato 2026-06-26) — `[E3·C3·S2·R2]`
- **Stato a codice**: positivity inference NON debole — `Assumptions::is_positive/is_nonnegative` copre Exp/Abs/Sqrt/Product/linear/Pow(even)/Sum + graph proof (`prove_relation`). Verificato che `mathematically_equal` con propagazione dominio:
  (a) **dimostra** sotto ipotesi: `ln(x²)=2ln(x)`, `exp(ln x)=x`, `ln(ab)=ln a+ln b`, `sqrt(x²)=x` con `x,a,b>0`;
  (b) **NON sovra-dichiara** senza ipotesi: `ln(x²)=2ln(x)`→false, `sqrt(x²)=x`→false (no silent-wrong, REGOLA ZERO);
  (c) identità incondizionate: `sqrt(x²)=|x|`→true.
  +2 test `A14TranscendentalEqual`.
- **Nota teorica**: la versione *"completa"* (zero-equivalence elementare totale) è **Richardson-indecidibile** — il contratto corretto è soundness (mai falso-positivo), non completezza. Questa è raggiunta. **Refs**: CAS-L2-19

### A15 · Units SI system — ✅ FATTO (verificato+completato 2026-06-26, era stale) — `[E3·C2·S3·R1]`
- **Stato a codice**: NON "mai iniziato" — sistema completo:
  - AST node `Quantity{value, SIDimensions}` (7 dimensioni base m·kg·s·A·K·mol·cd).
  - Registro `src/symbolic/units.cpp`: base SI + derivati (Hz,N,J,W,Pa,C,V,Ohm) + non-SI (ft,in,mi,lb,cal,eV) con scale Rational esatte.
  - `make_quantity_from_unit` / `convert_quantity` (con check mismatch dimensionale).
  - **Dimensional analysis** in simplify: prodotto combina dimensioni (chain.cpp:325); somma raggruppa per dimensione e **rifiuta** dimensioni incompatibili `1·m+1·s` (chain_sum.cpp:75, no silent-wrong); potenza scala dim·n.
  - **Completato 2026-06-26**: prefissi SI algoritmici (decomposizione `GHz→giga·Hz`, `um→micro·m`, set CGPM esaustivo applicato ai simboli prefissabili, exact-match prioritario) + **costanti fisiche esatte 2019-SI** (`make_physical_constant`: speed_of_light, planck_constant, elementary_charge, boltzmann_constant, avogadro_constant, standard_gravity — tutte Rational esatte).
  - 18 UnitsTest + 10 QuantityTest verdi. **Refs**: F6.6

### A16 · solve_inequality boundary `double` → `Rational` — `[E1·C2·S3·R2]`
- **Stato a codice**: `src/algebra/solve_inequality.cpp:98-100` usa `double low/high/tol` (la logica Sturm è OK, solo il boundary numerico no).
- **Ledger**: HC-F70-A21-NUMERIC-BOUNDARY · **Refs**: T-018

### A17 · IBP doppia applicazione su `Product(Log,…)` — ✅ FATTO (verificato 2026-06-26, stale) — `[E1·C2·S3·R2]`
- **Stato a codice**: ledger HC-F75-B1-IBP-DOUBLE-APPLY CHIUSO (fix 2026-06-10 simplify(vdu) pre-ricorsione + T-016 ILATE-class-da-base per `Pow`). Round-trip `IntegrateByPartsTest.{XLogX,XLogXSquared,LogXCubed}` 3/3 verdi 2026-06-26. Voce era stale nel tasklist.
- **Ledger**: HC-F75-B1-IBP-DOUBLE-APPLY (chiuso) · **Refs**: T-016

### A18 · GCD multivariato budget magico `*16U` → bound Mignotte — PARZIALE 2026-06-26 — `[E1·C2·S2·R2]`
- **Stato a codice**: ✅ galois `*8U+16U` rimosso → bound derivato `prime_budget + Σ bit_length(lc,disc_num,disc_den)` (galois_deg5.cpp:187, HPP-003B). 🟡 Residuo: pad `+16U` additivo in `divides_sparse_z` (brown_helpers/lc_scaling/fp_helpers) — diagnosi: il bound moltiplicativo sottostante è euristico (non provato), il budget→`return false` è un falso-negativo Cat-4 latente pre-esistente. Fix corretto specificato in HPP-003C (rely on well-ordering termination, cap provato `∏(deg_i+1)`, ritorna `Unimplemented` non `false`). R2/R3 hot-path, deferito con ledger.
- **Ledger**: HPP-003B (chiuso), HPP-003C (documentato) · **Refs**: T-020

### A19 · AstArena reset con root migration — `[E2·C3·S2·R3]`
- **Cosa**: migrare radici quando l'arena è resettata con nodi ancora referenziati. ⚠ tocca Arena core.
- **Ledger**: HC-F70-A31-MIGRATION-TODO · **Refs**: T-023

### A20 · Cycle detection simplifier/evaluator — `[E2·C3·S2·R3]`
- **Stato a codice**: parziale. **Cosa**: rilevare ricorsione infinita via depth max + fingerprint stack. **Refs**: CAS-L0-10

### A21 · Pivot Bareiss euristica contestuale — `[E2·C2·S2·R2]`
- **Cosa**: PivotScore contestuale (oggi: first non-zero). **Ledger**: HPP-025 · **Refs**: CAS-L1-17

### A22 · Padé coeff generici Q(π,e,√) — `[E2·C3·S1·R1]`
- **Cosa**: estendere Padé a campo algebrico generato dinamicamente (oggi vincolo Q), riusare LLL. **Refs**: F7.5.C2

### A23 · Structured Unimplemented diagnostic — `[E2·C1·S1·R1]`
- **Cosa**: payload strutturato `{module, fn, input_shape, reason, suggestion, ticket}` su ogni `Unimplemented`.
- **Ledger**: F0.8-ERROR-STRUCT · **Refs**: T-026

### A24 · Runner matrix·scalar / matrix±matrix — ✅ FATTO (verificato 2026-06-26, stale) — `[E1·C2·S2·R1]`
- **Stato a codice**: entrambi ledger CHIUSI. `matrix_adapter.hpp::evaluate_matrix_expression` (evaluatore top-level precedence-aware: scalar·matrix, matrix·matrix, matrix±matrix, matrix/scalar, unary −matrix; `Unimplemented` esplicito su scalar±matrix e div-by-matrix) + `try_evaluate_mattrace_wrapper`. Corpus matrix 79/79 = 100% (era 70.9%). `MatrixAdapterD2Test` 19/19 verdi 2026-06-26.
- **Ledger**: HC-F75-A2-MATRIX-SCALAR-OP, HC-F75-A2-MAXIMA-MATTRACE (entrambi chiusi)

---

## B — RESEARCH / WON'T-DO-NOW (aperta-permanente by-design)

Fuori scope deliberato (vedi `PLAN_HP_PRIME` §STATUS PERMANENTI). Parità HP Prime ≥95% **non** li richiede. Confermati assenti a codice = coerente con la scelta. **Non** sono lavoro dimenticato.

| Voce | Motivazione | Fallback attuale a codice | Ledger |
|---|---|---|---|
| Schönhage-Strassen NTT BigInt mul | SS NTT = mesi effort, GMP 30 anni tuning | Karatsuba + Toom-3 + Burnikel-Ziegler (tutti presenti) | HPP-F1.1-MUL |
| Collins CAD McCallum multivar | complessità doppia esponenziale | solo disequazioni 1-var via Sturm (presente) | F6.A / Fase 8 |
| Risch structure theorem full (Bronstein cap 9) | 600pp, multi-mese math expert | Liouville+Hermite+Trager+log-deriv (presenti) | — |
| Hypergeometric ₚFq recognition completo | tabelle Wilf-Zeilberger vaste | Gauss ₂F₁, Saalschütz, casi noti | — |
| Multi-sheet Riemann surface | AST single-valued incompatibile | branch principale (presente) | — |

---

## C — FATTO E VERIFICATO (non rifare)

Confermato a codice 2026-06-26. Elencato per evitare ri-lavoro (i vecchi tracker davano molti di questi come aperti).

**Foundation/BigInt**: Toom-3 (`bigint_mul_toom3.cpp`), Burnikel-Ziegler (`bigint_div_burnikel_ziegler.cpp`), Lehmer GCD, Miller-Rabin basi fisse. *(commit 70a1e52, 4131398)*

**Algebra univar**: Half-GCD (`polynomial_half_gcd.cpp`), Berlekamp, sparse interpolation, cyclotomic detection >724 *(commit dd27a51, T-022)*.

**Algebra multivar/ext**: Wang EEZ + Kronecker fallback, Zippel sparse GCD CRT+Farey *(commit ad023bf, T-006)*, **CRT modular resultant** (`polynomial_resultant_crt.cpp`, Collins 1971 — *era erroneamente dato aperto come T-012*), FGLM + F5 Groebner, tower 2-level + nested squarefree.

**LinAlg**: **Householder QR simbolico razionalizzato** (`matrix_qr.cpp` — *PLAN_HP_PRIME lo dava aperta-perm/MGS, è stale: codice ha risolto l'AST-explosion via reflector `I−2vvᵀ/N` senza √ al denominatore*), Cholesky LDLᵀ, fresh-symbol pervasivo, adaptive G7/K15 quadrature *(commit 5268e90)*.

**Calculus**: Kovacic Case 1 (Laurent √r poli≥4, *commit 709a437, T-008*) + Case 2; Risch trial-constants rimossi (`integrate_risch.cpp:69`), Rothstein-Trager generalizzato, Hermite con D arbitrario; Zeilberger/Gosper/Petkovšek higher-order *(commit 8440ba5, 3addb54, T-027)*; ∫xᵏ/√(c−dx²) *(commit 830f80e, T-055)*, ∫1/√(Ax²+Bx+C) *(6c5ae48)*, inverse-trig IBP, ∫log(x+√(x²+a)) *(ec9e065)*.

**Complex/branch-cut**: BC-1..BC-3 + BC-1b corrections, ln(a+bi) esatto *(commit 89eae9e, T-017)*, cyclotomic↔exp equal *(4597332, T-025)*, extended-real arithmetic *(b145cb9)*, UnwindingNumber builtin.

**Numerica**: fsolve tolerance ctx-param (`fsolve_tolerance_bits`, T-019), solve_inequality via Sturm, RootOf isolation bounds, Bessel identities *(ade20ff)*, statistics package (`src/statistics/`).

**Infra**: **Phase-7 anti-monolith DONE** — tutti i 23 file ex >500 LOC ora <500 *(commit 7aa5f2e + T-046..050, T-029..051)*; flaky cos(7π/16) fixed via co-function *(d03ff6d, T-004)*; rapidcheck property tests; golden aggregate **94.5%** (F7.5).

---

## D — Grafo dipendenze (task aperti)

```
A1 (Risch RP-2 Hermite) ──> A7 (Slater/Meijer G fallback, raccomandato)
A9 (tower ≥3) ──> A4 (Smith Q[x]) [parziale], A10 (eigen/Jordan RootOf)
A11 (Gruntz dinamico) ──> A13 (residue) [indipendenti ma stesso dominio limit]
A2 (hard-timeout) — indipendente, ALTA priorità (S5 hang)
A16..A24 (debiti minori) — indipendenti, inframmezzabili
```

## E — Ordine raccomandato

1. **A2** (hard-timeout S5 — rischio hang produzione) + **A3** (perf-hang)
2. Debiti rapidi E1: **A16, A17, A18, A23, A24**
3. **A5** (Frobenius) · **A13** (residue) · **A4** (Smith Q[x]) — gap funzionali medi
4. **A1** (Risch RP-2) → **A7** (Slater/Meijer) — research calculus, alto impatto integrate
5. **A6** (Stauduhar) · **A9** (tower≥3) · **A10/A11/A12/A14** — research algebra/limit
6. **A15** (Units SI) · **A8** (Kovacic n=12) · **A19-A22** — bassa priorità

---

## Vincoli invariabili (CLAUDE.md)

BigInt only (no int64/double nel core) · Structural sharing · AstArena bump · LPO orientation
· `Result<T>` (no throw/catch) · 500 LOC/file (hard 550) · fresh-symbol pervasivo · REGOLA ZERO
(no via facile, no hardcode non-ledgered) · REGOLA 0.1 (spec read first) · REGOLA 0.2 (no test disable)
· timeout test espliciti (`scripts/test_quick.sh`).
