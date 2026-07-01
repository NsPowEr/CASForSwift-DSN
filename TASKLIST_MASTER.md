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

### A1 · Risch RP-2 Hermite reduction parametrica (df>0) — 🚧 BLOCCATO (deferito; vedi A7 per la research attiva) — `[E4·C4·S3·R2]`
> **DEFERRAL 2026-06-28**: NON è la prossima research. Il ramo df>0 è irraggiungibile dal corpus (verificato instrumentando) → impl. alla cieca = codice non verificabile su hot-path Risch (silent-wrong > Unimplemented, REGOLA ZERO). **Sblocco = primo step sotto** (costruire una torre differenziale con `deg_t(f)>0`); solo allora implementare Bronstein §6.5 con cross-check Maxima. Direzione research scelta ora = **[[A7]]** (unblocked, additiva, R1). Riprendere A1 solo dopo aver costruito il caso-test che esercita il sito.
- **AUDIT 2026-06-27** (spec `Risch_Transcendental_Cap8.md` + `Risch_Hermite_Cap5.md` lette): la premessa "scaffold da sostituire" era **fuorviante**. Stato reale:
  - `risch_rde_bronstein_hermite.cpp::risch_rde_hermite_parametric_stub` è **DEAD CODE**: compilato (CMakeLists:220) ma **mai chiamato** da nessuno (grep verificato). NON è il punto attivo.
  - Il solver parametrico **reale** `solve_risch_de_parametric_field` (`risch_rde_bronstein.cpp:249`) è **già implementato**: denominator bound (Q,B,D), trasformazione `f_new=f−D'/D`, `g_new=g·D`, bound di grado N (Log/Exp, §6.4), ricorsione tower `solve_recursive` con correzioni exp/log. Trial constants **già rimossi** (grep vuoto).
  - **GAP REALE A1** = ramo **`df > 0`** a `risch_rde_bronstein.cpp:362-368`: `return Unimplemented("df > 0 not implemented in parametric field solver")`. Tutto `df ≤ 0` funziona; manca il caso *non-cancellation* in cui `deg_t(f) > 0` (il termine `f·y` domina): `deg(y)=dg−df`, `lc(y)=lc(g)/lc(f)`, poi ridurre+ricorrere (Bronstein 6.5/7.4/8.4 PolyRischDE).
- **Cleanup 2026-06-27**: dead-code scaffold `risch_rde_bronstein_hermite.cpp` **rimosso** (file + CMakeLists); diagnostic del ramo df>0 reso strutturato (Bronstein §6.5/§7.4/§8.4). 69/69 Risch verdi.
- **🚧 BLOCCO alla chiusura (verità da REGOLA ZERO)**: il ramo df>0 è **irraggiungibile** dall'intero corpus + 16 integrandi-torre costruiti a mano (verificato instrumentando il sito). Senza un caso-test che lo eserciti non c'è verifica possibile → implementarlo alla cieca = codice non verificabile sul hot-path Risch (silent-wrong peggio di Unimplemented).
- **🔬 CAUSA RADICE TROVATA 2026-06-29**: l'irraggiungibilità è **strutturale**, non di corpus. `solve_risch_de_parametric_field` (host del sito df>0) è **dead code: zero caller esterni**, solo auto-ricorsione (grep tree-wide + graph). Il df>0 non sarà mai raggiunto finché la funzione non viene **wired-in** dal path integrate reale. Estratto come task prerequisito **A26** (`[E4·C4·S3·R2]`, deciso 2026-06-29). A1 resta BLOCCATO **fino al completamento di A26**; solo allora df>0 ha un caso-test verificabile (Bronstein §6.5 + cross-check Maxima).
- **Spec (verificata 2026-07-01)**: spec algoritmica **self-contained** = `MISSING_FEATURES_SPECS/Symbolic_Integration_I.md` (testo Bronstein completo con dimostrazioni): **§6.4 SPDE** (Rothstein, algo `SPDE`), **§6.5 The Non-Cancellation Cases** (`PolyRischDENoCancel1/2/3`, Lemma 6.5.1 — dà `deg(q)=deg(c)−deg(b)`, `lc(q)=lc(c)/lc(b)` = esatto per df>0), **§7.1 The Parametric Risch Differential Equation** (`ParamRischDE` + *Parametric SPDE* + *Parametric Poly Risch d.e. – no cancellation*, righe ~7078/7190) = il target reale del df>0 **parametrico**. `Risch_Transcendental_Cap8.md` + `Risch_Hermite_Cap5.md` = solo riassunti high-level (puntano al libro, NON self-contained). ⚠ **Correzione citazione**: "§6.5/§7.4/§8.4" nei commenti storici è errata — **§7.4 non esiste** (cap.7 = 7.1/7.2/7.3) e **§8.4 = Hypertangent Case** (solo tan/complessi, fuori scope torri log/exp). Ref corretto = **§6.5 + §7.1**.
- **Ledger**: HC-F8-PENDING-17 · **Refs**: T-007, Task 17, F5.1, F7.5.B3, CAS-L1-02

### A2 · Hard-timeout cancellation token su tutti i path integrazione — RISOLTO 2026-06-26 — `[E2·C2·S5·R2]`
- **Stato a codice**: ✅ chiuso. L'ultimo loop combinatorio non interrompibile dell'integrazione (le due eliminazioni di Gauss pure-Rational in `integrate_risch_rde.cpp`) ora polla `check_interrupt()` per colonna pivot; rimosso lo stub morto del poll mai cablato; dispatcher `solve_risch_de_q` reso fedele alla cancellazione (non ingoia più il Timeout di simplify, non maschera l'interrupt come Unimplemented). +2 test in `test_integrate_interrupt.cpp` (RDE solver level). 147 test Risch/integrate verdi.
- **Residuo non bloccante**: fedeltà error-kind su interrupt asincrono *durante* `parse_polynomial` interno a un ramo (rietichettato Unimplemented, finestra ~µs, nessun hang) — vedi ledger HC-F70-A33.
- **Ledger**: HC-F75-A3-HARD-TIMEOUT (chiuso), HC-F70-A33-POLL-COVERAGE (aggiornato 2026-06-26) · **Refs**: T-015, T-024

### A3 · Perf-hang fattorizzazione torri / VanHoeij — 🔨 QUASI-FATTO 2026-06-29 (solo deg-16 Hensel residuo) — `[E3·C4·S2·R2]`
- **FATTO 2026-06-29** (analisi `sample` + instrumentazione `[WANG]`, spec `Hensel_Lifting.md` letta): la diagnosi vecchia ("suspect 5c72bc0 / together() GCD") era **ERRATA**. Tre fix reali, zero debito, zero hardcode:
  - **(1) LLL intero fraction-free** (`lattice_lll.cpp`, Cohen 2.6.7 / de Weger): la Gram-Schmidt `Rational` rifatta per-swap esplodeva in frazioni sui reticoli van Hoeij. Ora d[i]/λ[i][j] interi Hadamard-bounded; fallback Rational solo per basi rank-deficient. 25 test LLL/VanHoeij correttezza verdi.
  - **(2) Kronecker cancellabile** (`factorize_kronecker.cpp`): per norma torre perfect-power il recombination cade su ricerca combinatoria esponenziale che non pollava cancellazione → hang. Ora polla `ctx.past_hard_deadline()` (deadline **opt-in** del chiamante; default illimitato per non rompere factor lente-legittime, es. Trager Q(∛2) 27s) + propaga Timeout. `factor_polynomial_tower[_n]` pubblica `set_hard_deadline(now()+ctx.timeout())`.
  - **(3) Verdetto irriducibilità lucky-prime** (`factorization_wang_eez.cpp`, HC-F8-SD3-VANHOEIJ-SLOW ✅): su primo **lucky** (`f` squarefree mod p, check `gcd(f,f') mod p` costante) la recombination **Hensel esaustiva** è completa → "nessun fattore" = prova di irriducibilità → ritorna `{f}` invece di cadere su Kronecker illimitato. Si fida **solo** del verdetto esaustivo (NON di van Hoeij, incompleto) e **solo su lucky prime**. Un primo tentativo "trust van Hoeij" era stato revertito (silent-wrong su RedundantGenerator) — il bug era saltare il fallback esaustivo; risolto.
  - **Esito**: quarantena 7→4. SD3 **219ms** (era hang >400s, rimosso da SLOW_OK). AntiHardcode 3.2s, PreservesLeadingCoeff 5.8s, IrreducibleX2Minus7 3.5s. **243** factor/poly/LLL/Galois/Trager + **352** consumer (integrate/simplify/solve/residue) verdi, **0 regressioni**.
- **🔧 RESIDUO**: estratto in task separato **A25** (deg-16 Hensel-lift mod-p^k, perf-only `[E3·C3·S2·R3]`, ⏸️ in pausa). La correttezza di A3 è **CHIUSA**; ciò che resta è puro debito performance, non più tracciato qui.
- **Ledger**: HC-F8-FACTORIZATIONTOWER-PERF (PARZIALE: solo deg-16 Hensel → A25), HC-F8-FACTORIZATIONTOWER-AntiHardcode-X2Minus2-Sqrt3Sqrt5 (✅ RISOLTO), HC-F8-SD3-VANHOEIJ-SLOW (✅ RISOLTO) · **Refs**: A25, BUG-HANG-002

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

### A7 · Slater pFq → Meijer G fallback integrator — 🚧 SPEC-BLOCKED (nessuna spec formale esiste) — `[E4·C4·S3·R1]`
- **Stato a codice**: `BuiltinOp::MeijerG` **non esiste**. Nessun file Slater/Meijer.
- **⚠️ CORREZIONE 2026-06-28**: la spec citata in precedenza (`Special_Fn_Identities.md`) **NON è di Meijer-G** — è **F7.5.E1** (identità Γ/B/ζ/erf), task **diverso e di fatto FATTO** (special_fn golden = **100%**, vedi §C). Quindi **A7 non ha spec formale**: scrivere Meijer-G ora violerebbe REGOLA 0.1 ("codice senza spec = INVALIDO") e REGOLA ZERO (formule da memoria = allucinazione). **Bloccante reale**: prima va prodotta/approvata una spec `MISSING_FEATURES_SPECS/Meijer_G_Slater.md` (definizione G, identità Erdélyi-Slater, formula integrale, trasformazioni Bailey, scope) — decisione architetturale, non auto-derivabile.
- **Cosa (sketch ingegneristico, valido SOLO dopo la spec)**: definire `MeijerG`, dispatcher, wiring come passo finale post-Risch.
  - Sequenza incrementale (ogni step = commit + gate verde): (1) enum `BuiltinOp::MeijerG` in `include/cas/builtin_functions.hpp` + round-trip print/parse; (2) rappresentazione `MeijerG[[a_p],[b_q]](z)` param esatti + riduzioni degeneri come simplify rules; (3) bridge `pFq→MeijerG` (Slater) con cross-check Maxima (`hgfred`/`makegamma`); (4) `∫ MeijerG dz` (integrale di G = G con indici shiftati); (5) wiring **ultimo** in `src/calculus/integrate_core.cpp` dopo il fallimento Risch (top-level `src/calculus/integrate.cpp:421`), mai prima.
  - **File chiave**: `include/cas/builtin_functions.hpp`, `src/calculus/integrate_core.cpp`, `src/calculus/integrate.cpp:421`.
- **Spec**: ❌ DA CREARE (`Meijer_G_Slater.md`) — la vecchia ref era errata · **Refs**: T-013, Task 22-23, F7.A-B

### A8 · Kovacic Case 3 completo (icosaedrico n=12) — `[E3·C4·S1·R1]`
- **Stato a codice**: `ode_kovacic_case3.cpp` ha tetraedrico/ottaedrico; n=12 = solo diagnostico (recurrence blow-up).
- **Cosa**: ottimizzazione recurrence per A₅ n=12. **Ledger**: HC-KV-06 · **Refs**: T-028, T-009 (parziale), F5.C

### A9 · Tower algebrico ≥3 livelli / multi-β nested non-squarefree — `[E3·C4·S2·R2]`
- **Stato a codice**: `algebraic_tower_primitive_nested.cpp` esiste (F3.4-DEBT-01) ma → `Unimplemented` se la min-poll liftata non è squarefree. Bridge solo 2-level (`detect_two_level_tower`).
- **Cosa**: nesting multi-β >1 + factor Q(β)[x] reducible. **Refs**: F3.5, F3.D, Task 7, CAS-L3-06

### A10 · Eigenvalues n>3 / catene di Jordan sotto RootOf — ✅ FATTO (probe chiuso 2026-06-27) — `[E3·C3·S3·R2]`
- **Stato a codice**: routing RootOf GIÀ presente (`matrix_jordan.cpp:150` → `null_space_over_extension` quando autovalore è `RootOf`). Algoritmo generale (Filippov `n_k=2d_k−d_{k−1}−d_{k+1}`, catene top-down) senza cap su `n`/molteplicità.
- **CHIUSO 2026-06-27**: il probe richiesto è verde. `RootOf_Multiplicity3_CompanionDeg6`: companion di `(x²−2)³ = x⁶−6x⁴+12x²−8` (non-derogatoria ⇒ blocco Jordan singolo size 3 per ciascun coniugato ±√2), **n=6>3, autovalore RootOf molteplicità 3, catena lunghezza 3**, certificato `P·J·P⁻¹==A` esatto su Q(√2) (1.3s). Esercita `null_space_over_extension` su `(A−√2·I)^k` fino a k=3. La descrizione storica "non affidabili" era stale. 5/5 JordanCertTest verdi.
- **Refs**: CAS-L2-02, F4.4

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

### A16 · solve_inequality boundary `double` → `Rational` — 🧊 DEFERITO post-parità (ledgered, era mis-scoped E1) — `[E3·C3·S2·R2]`
- **AUDIT 2026-06-27** (ledger HC-F70-A21 autoritativo + codice letti): la voce "E1, solo `low/high/tol` double→Rational" era **fuorviante**. Stato reale:
  - `low/high/tol` (`solve_inequality.cpp:98-100`) sono il **confine numerico Cat-4 deliberato**: bound/tolerance stoccati `long long` in `CASContextParams`, convertiti a `double` SOLO al call-site `find_polynomial_roots_sturm` (layer numeric). REGOLA 1 vieta `double` nel *core simbolico*, lo **ammette ai confini numerici** → non è un debito, è by-design documentato.
  - Debito reale (più profondo, **non** quello del titolo originale): (a) output radici via `double_to_rational_approx` (riga 115) = boundary **approssimati**, e (b) `InequalityInterval.lower/upper` è `optional<Rational>` → **non può rappresentare radici irrazionali esatte** (√2 ecc.) neppure col fix puro-Rational.
  - **Fix completo (ledgered, deferito Fase 8 post-parità)**: `find_polynomial_roots_sturm_rational` (bisection puro Rational/BigFloat) per (a); per (b) servirebbe `InequalityInterval` boundary `ExprPtr`+`RootOf` via `find_polynomial_isolating_intervals` (esiste, F8.0-5.4). Re-stima realistica **E3·C3**, non E1.
  - **Decisione**: NON toccare ora — versione cosmetica `low/high/tol`→Rational = pigra/no-op (REGOLA ZERO), versione vera contraddice il deferral post-parità deciso (goal corrente = parità HP Prime ≥95%, §B). Blast radius contenuto (3 file, nessun chiamante esterno) quando si riprenderà in Fase 8.
- **Ledger**: HC-F70-A21-NUMERIC-BOUNDARY (PARTIAL, deferral documentato) · **Refs**: T-018

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

### A25 · deg-16 Hensel-lift mod-p^k perf (norma Trager torri 2-livello) — ⏸️ IN PAUSA (perf-only, isolato da A3 2026-06-29) — `[E3·C3·S2·R3]`
- **Origine**: estratto da A3 dopo che LLL fraction-free + Kronecker deadline + verdetto lucky-prime hanno risolto il collo di bottiglia correttezza. **Nessun problema di correttezza**: le fattorizzazioni/antiderivate prodotte sono giuste, il debito è **puramente performance**.
- **Sintomo a codice**: `FactorizationTowerTest.SplitsProductOfQuadraticsOverQSqrt2Sqrt3` e `FactorizationTowerTest.SplitsX4Minus10X2Plus1OverQSqrt2Sqrt3` >200s. Lo splitter deg-4 su torre `Q(√2,√3)` produce una norma di Trager **deg-16 a coefficienti enormi**; col LLL e Kronecker già risolti il profilo (`sample`) si concentra su `hensel_lift` / `quadratic_step` / `poly_*_mod` — il lift mod-p^k della norma deg-16.
- **Fix corretto**: ottimizzare l'Hensel-lift mod-p^k (quadratic Hensel / aritmetica modulare più veloce / scelta primo+precisione mirata sui coeff-bound Mignotte della deg-16), NON euristiche di forma. Spec `Hensel_Lifting.md` (REGOLA 0.1) da rileggere prima dell'impl.
- **Quarantena**: i 2 test restano in `scripts/test_quarantine.txt` (CEILING 4), abilitati via `--slow`/filtro esplicito. `PrimitiveElementTest.SqrtTwoSqrtThreeSqrtFive` (3-livello primitive element) è facet correlato ma distinto.
- **Ledger**: HC-F8-FACTORIZATIONTOWER-PERF (residuo deg-16 Hensel) · **Refs**: A3 (chiuso), BUG-HANG-002

### A26 · Wire-in solver parametrico Risch DE a livello tower (prerequisito A1 df>0) — ✅ CHIUSO per la parte concludibile senza A1 (f=0 + f≠0 razionale + WeakNormalizer chiusi; residuo wiring gated A1) — `[E4·C4·S3·R2]`
- **FATTO 2026-06-29** (validazione-prima-del-wiring, spec Cap8/Hermite/Risch_Algorithm lette): la funzione orfana è ora **raggiungibile + testata** via harness diretto `test/unit/calculus/test_risch_parametric_tower.cpp` (back-substitution SOUND `D(y)+f·y≡Σc_i g_i` con `field.derive`).
  - **✅ memory-safety**: heap-buffer-overflow ASan a `risch_rde_bronstein.cpp:384/425` (`g_polys[s][idx]` con `idx>deg` → `operator[]` non-const non-checked) corretto via `const PolyExpr&` (overload bounds-checked, OOB→0). Era una mina: wiring alla cieca nel hot-path avrebbe introdotto corruzione memoria.
  - **✅ ramo esponenziale validato sound** (`Exp_F0_GExp_AllSound`).
  - **✅ df>0 RIPRODOTTO deterministicamente** (`Log_FEqualsTheta_DfPositive_ReachesNonCancellation`): `f=t` su `Q(x,log x)` raggiunge il ramo non-cancellation → **sblocco concreto di A1** (16 integrandi a mano non ci riuscivano). Quando A1 implementerà Bronstein §6.5, quel test passa da Unimplemented a verify_field_de.
- **✅ FASE 2 — caso f=0 RISOLTO 2026-06-30** (commit `7187c77`): il ramo log (primitivo) genera dalla correzione `i·y·θ'` (θ'=1/x) una forzante **razionale** che `solve_risch_de_parametric_q` (solo `Q[x]`) rifiutava. Nuovo `solve_param_limited_integration_rational_q` (`src/calculus/risch_parametric_rational.cpp`) = **rational limited integration su Q(x)** (Bronstein §7.2/§7.3): integra ogni `g_i`, split `R_i` razionale + atomi log/arctan, null space matrice atomi → `y=Σc_i R_i`, **verifica back-substitution** (sound-by-construction, mai silent-wrong). Wired nel base case `solve_risch_de_parametric_field` (ext_idx==0, fallback f==0). Flippa `Log_F0_G1_PrimitiveDescent_SolvedSound` (era Unimplemented → `y=x,c=1` verificata) + 3 test diretti. 81 risch/parametric verdi, 0 regressioni.
- **✅ FASE 3 — caso f≠0 razionale RISOLTO 2026-06-30** (commit `7725c25`): nuovo `solve_param_risch_de_rational_q` = **ParamRischDE generale su Q(x)** (Bronstein §6.1/§7.1). Ansatz `y=P/D`, `D=lcm(den f,den g_i)`; identità `D·P'+(Fn−D')·P−Σc_i(Gn_i·D)=0` → sistema omogeneo → null space. **Degree bound cancellation** (`deg_H==deg_D−1` → grado omogeneo `−lc(H)/lc(D)`; trova `y=x²` di `y'−2y/x=0`). **Doppia verifica**: residuo Q esatto (interno, robusto) + `verify_field_de` simbolico (test) — entrambi verdi. Output integer-cleared. Wired (ext_idx==0, fallback f≠0). +4 test. 85 verdi, 0 regressioni. *(Sospetto bug-simplify su `(−x/2)/x²` investigato = FALSO ALLARME; engine corretto.)*
- **✅ FASE 4 — WeakNormalizer RISOLTO 2026-07-01** (`inflate_denominator`, spec Cap5/Cap8 rilette): **denominator inflation Bronstein 6.1.1**. Correzione diagnosi: nella convenzione `y'+f·y=g` è il residuo intero **POSITIVO** (non negativo) al polo semplice a creare il polo in `y_h=(x−α)^{−n}` oltre `lcm`; il negativo dà `y_h` polinomiale (già coperto). Fix: residuo `n` ai radici di `gcd(fn−n·fd', s)` (Rothstein-Trager) per `n∈[2,cap]` → `D *= s_n^{n−1}`. **SOUND** (back-sub scarta non-soluzioni). Cap = `ctx.max_risch_rational_ansatz_degree()` (configurabile, zero magic). +3 test critici `WeakNormalizer_*` con **probe load-bearing** (disabilitato→3 fail, riabilitato→pass). 13 ParametricTower + 151 risch/integrate verdi, 0 regressioni.
- **🔧 RESIDUO (non-A26, incompletezza MAI sbagliata)**: poli di f a **molteplicità ≥2**, residuo intero **oltre cap**, residui **algebrici** → Unimplemented pulito; inflate del solver **non-parametrico** `solve_risch_de_rational_q` (R3 hot-path, task separato); wiring in `integrate()` reale subordinato a **df>0 (A1)**.
- **Finding verificato 2026-06-29** (grep tree-wide + graph): `solve_risch_de_parametric_field` (`risch_rde_bronstein.cpp:249`) — la funzione che l'audit 2026-06-27 chiamava "il solver parametrico reale già implementato" — è **dead code**: **zero caller esterni**, raggiunta solo dalla propria auto-ricorsione `ext_idx-1` a L399 (che non parte mai). L'audit aveva verificato che il codice è *scritto* (trial-constants rimossi, logica df≤0 presente) ma **non** che è *raggiunto*. Il path Risch attivo è `solve_risch_de_general` → `solve_risch_de_field` (non-parametrico, single-`g`); il parametrico base-case `solve_risch_de_parametric_q` (`risch_parametric.cpp`) esiste ma non è invocato dal field solver.
- **Conseguenza su A1**: il sito `df > 0` (`risch_rde_bronstein.cpp:362`) è dentro questa funzione orfana ⇒ **strutturalmente irraggiungibile**, non solo "irraggiungibile dal corpus". Implementare df>0 lì = math non verificabile su funzione mai chiamata = **REGOLA ZERO violata**.
- **Scope A26** (vero unblock di A1, deciso 2026-06-29): costruire l'integrazione parametrica **a livello tower** che invoca `solve_risch_de_parametric_field` dal path integrate reale (es. limited-integration / coupled-system / parametri), così che la RDE parametrica diventi raggiungibile **e testabile** con cross-check Maxima. Solo a quel punto il ramo df>0 (Bronstein §6.5/§7.4/§8.4) ha un caso-test che lo esercita ⇒ A1 implementabile in sicurezza. Spec `Risch_Transcendental_Cap8.md` + `Risch_Hermite_Cap5.md` (REGOLA 0.1) da rileggere.
- **Ledger**: HC-F8-PENDING-17 (condiviso con A1) · **Refs**: A1 (sbloccato da qui), A7, T-007

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

### B.1 · Golden `simplify`-fail = soundness deliberata (NON fixare alla cieca)
Audit golden 2026-06-28: i 5 fail `simplify` vs Maxima (`exp(log(x))→x`, `log(x^n)→n·log(x)`, `log(1/x)→−log(x)`, `abs(x^2)→x^2`, `0^x→0`) **non sono bug**: valgono solo sotto ipotesi `x>0` / dominio reale. Maxima sovra-semplifica assumendo dominio reale; il CAS è **complex-aware conservativo** (A14 soundness, `is_known_nonnegative` senza regola even-power perché `x=i⇒x²=−1`). Allinearsi a Maxima = introdurre falsi-positivi (REGOLA ZERO). Il fix corretto sarebbe un **default real-domain configurabile** (decisione architetturale cross-cutting per parità HP Prime, fuori scope sessione singola), non hack per-regola. **Golden aggregate 2026-06-28: ~97.7%** (corpus integrate basic.jsonl = 140 entry, 116/6/18 fresh; ⚠ `build-golden/report.json` *merged view* è stale/inaffidabile — fidarsi del per-area output, non del totale 230≠140).

### B.2 · Golden `integrate`/`diff`-fail = gap canonicalizzazione trig (risposte già corrette) — 🔨 PARZIALE 2026-06-28 — `[E3·C4·S2·R3]`
- **FATTO 2026-06-28** (commit `3b15be1`): aggiunto building-block **sound** `trig_exponential_zero_diff` (`algebraic_equal_trig_exp.cpp`) — forma canonica Laurent in `z=e^{ix}` per **polinomi trig** (e razionali via cross-multiplication), wired in `mathematically_equal` dopo `weierstrass_zero_diff`. Prova product-to-sum freq-diverse + Fourier (`cos(x)cos(2x)=(cos x+cos3x)/2`, `sin³x cos²x` linearizzato) dove il path `t=tan(x/2)` cede. SOUND (solo identità esatte, bail su non-trig-poly → mai falso-positivo; A14 verde). +4 test, suite quick 2489 verde, zero regressione. Anti-monolith: split `algebraic_equal_subset_risch.cpp` (file era 504 LOC).
- **✅ SUB-GAP abs/half-angle CHIUSO 2026-06-28** (commit `ec651d5`): esteso `trig_exponential_zero_diff` da polinomi a **RAZIONALI trig** — la forma esponenziale ora modella `N(z)/D(z)` (Laurent num/den), gestendo **potenze intere negative** e **trig al denominatore** via cross-mult esatta di mappe Laurent. `N/D≡0 ⟺ N≡0` (D≢0), sound perché un Laurent che si annulla su un arco del cerchio unitario è il poly nullo. **Diagnosi corretta**: l'`abs` NON era il problema — `diff` già applica `d/dx ln|u|=u'/u` e droppa l'abs; il blocker reale erano le **potenze negative di half-angle** (`cot(x/2)`, `cos(x/2)⁻²`) nel derivato. **Golden flip verificato (sound)**: integrate `116→118` (famiglia `1/sin x`), diff `78→79` (`sin/cos`). +3 test rational-equiv + 1 rational-soundness-reject (`TrigIdentitiesTest` 21/21), A14 verde, file 247 LOC.
- **✅ FLIP `∫x²/(x²−1)` 2026-06-28** (commit `7ee0935`): root cause era un **gap in `diff`**, non equivalenza — `differentiate_visitor.cpp` faceva bail `Unimplemented` su `ComplexLit` (`"Complex literals differentiation not implemented"`). Ma `ComplexLit` è una **costante numerica** → `d/dx≡0` come `IntegerLit`/`Constant`. L'integratore razionale emette `∫x²/(x²−1)=x−i·arctan(−i·x)` (con `ComplexLit(−i)`, non `Symbol`); `antiderivative_equivalent` non riusciva a derivarla. Fix: `ComplexLit → 0`. Ora `diff(x−i·arctan(−i·x))=1+1/(x²−1)` matcha Maxima. **Golden integrate 118→119**. Sound (costante→0, nessuna ipotesi di dominio). +2 test, quick suite 2495 verde, zero regressione.
- **✅ FLIP `∫√(x²−1)` 2026-06-29** (commit `c27f1ee` + fix simplify): erano DUE bug. (1) **`unknown(x)`**: NON era un bug d'integrazione — l'antiderivata corretta è `½x√(x²−1)−½·acosh(x)`, ma `simplify_node` (`simplify_functions.cpp`) ricostruiva ogni `FuncCall` da `node.func_id`; `acosh/asinh/atanh` NON sono enumerati in `BuiltinOp` (→ `Unknown`), e `builtin_name(Unknown)="unknown"` → il nome veniva corrotto in `unknown(x)`. Fix: ricostruzione name-preserving quando `func_id==Unknown` (generale, ogni funzione non-enumerata). (2) **Equivalenza radicale**: nuovo building-block sound `radical_zero_diff` (`algebraic_equal_radical.cpp`) — estensione singola radice `Q(x)(√p)`, rappresenta il diff come `(A+B·s)/(C+D·s)`, `A,B,C,D∈Q[x]`, `s²→p`; `diff≡0 ⟺ A≡0 ∧ B≡0` (sound incondizionato). Wired in `mathematically_equal` **prima** del dispatch RootOf (un `√quadratica` è RootOf di `t²−p`, altrimenti short-circuit a false), con `operation_active_` azzerato temporaneamente (il `polynomial_normal_form`→`expand` non è idempotente sotto operazione attiva → coordinate non ridotte → falso negativo). **Golden integrate 119→120**. +5 test, quick suite verde, diff golden 79/1 invariato.
- **🔧 RESIDUO (real-domain `abs`, NON forzare — vd. B.1)**: restano 2 fail integrate **inverse-hyp/algebrici a dominio reale** — `∫1/sqrt(x²+1)` (CAS `ln|x+√(x²+1)|` vs Maxima `asinh x`: bridge `asinh→ln` già in `hyperbolic_normalize`, blocca solo l'**`abs`** — `x+√(x²+1)>0` vale solo su ℝ), `∫1/sqrt(a²−x²)` (`arcsin(x/a)` vs `arcsin(x/|a|)`, `abs` su parametro, vale solo `a>0`); e `diff(x/sin x)` = **polinomio×trig misto** (`x` Symbol nudo × trig → bail sound). I casi `abs` sono la tensione real-vs-complex di B.1: fix corretto = **default real-domain configurabile** (cross-cutting), non hack per-caso.
- **🔬 INDAGATO 2026-06-28, NON gating (no fix)**: il `diff` timeout su `arctan(tan(x/2)/√3)` **non blocca alcun golden** — `∫1/(2+cos x)` già PASSA (confronto antiderivate diretto, non dipende da quel `diff`). Root cause accertata (diag): `integrate` ritorna l'antiderivata **non-canonica** (`materialize_expr` senza simplify → scaffolding `1·…`, termine spurio `−(0·x)` dalla back-sub Weierstrass); `diff` su quel wrapper sfora il budget simplify 1000ms (`diff(raw F)` ~1334ms→timeout vs `diff(simplify(F))` ~1044ms OK). **Fix tentato e REVERTITO**: simplify dell'output di `integrate` (con fallback su errore) → **0 flip golden + 3 regressioni unit** (`XTimesAsinhX`: il `diff` della forma asinh semplificata va a sua volta in timeout; `IntegratesClassicalTrigSubstitutionRadicals`, `IntegralOfErfTimesExp`). Conclusione: il timeout è **debito perf R3 di `diff`/`simplify`** su `arctan(g)` con `g` a potenze negative, non risolvibile via canonicalizzazione dell'output senza regredire altre famiglie. Eventuale fix corretto = riduzione mirata del quoziente `u'/(1+u²)` in `diff` (multi-sessione), oppure cleanup *leggero* (prune 0-termini/flatten `1·`) localizzato in `integrate_weierstrass` — ma **valore nullo** (nessun golden gated). Tutti C4·R3, sessioni a sé.
Audit 2026-06-28: i fail `integrate` (~6 su basic.jsonl: `1/sin(x)`, `1/(1+sin x)`, `1/(2+cos x)`, `sin³cos²`, `1/(sin·cos²)`, `1/(sin·cos)` …) e `diff` (2: `x/sin(x)`, `sin/cos`) **non sono bug di correttezza**: le antiderivate CAS sono **verificate corrette** dai unit test via `D(F)==integrand` (`test_integrate_inverse_trig.cpp::expect_integral_correct` passa). Falliscono nel golden perché `antiderivative_equivalent` (`test/golden/integrate_equiv.hpp`) confronta `D(cas)` vs `D(maxima)` e **`simplify` non riesce a provare `D(cas)−D(maxima)=0`** per forme trig equivalenti (Weierstrass `tan(x/2)=cos(x/2)⁻¹·sin(x/2)`, half-angle `cos(2x)`). Es. `∫1/sin(x)`: CAS `ln|tan(x/2)|` (forma da manuale) = Maxima `½log((cos x−1)/(cos x+1))`.
- **Root cause reale**: il simplifier manca di una **forma canonica per razionali trigonometrici** (power-reduction/double-angle folding orientato, es. `½cos(2x)+½ ↔ cos²x`, Weierstrass-collapse). Migliorarlo flipperebbe questi entry E rafforzerebbe `mathematically_equal`/diff/integrate-verify in generale. Research-grade, hot-path simplifier (R3).
- **NON fare ora alla cieca**: rischio regressione su hot-path; e modificare l'oracolo `antiderivative_equivalent` per aggirare (es. confronto `D(cas)==integrand`) è sound ma marginale (testato 2026-06-28: +1 entry, bloccato dallo stesso limite simplify) e tocca un benchmark → cautela REGOLA 0.2 (revertito). **Il vero lavoro è la canonicalizzazione trig nel simplifier.** **Refs**: F7.5 trig, golden integrate/diff corpus

---

## C — FATTO E VERIFICATO (non rifare)

Confermato a codice 2026-06-26. Elencato per evitare ri-lavoro (i vecchi tracker davano molti di questi come aperti).

**Foundation/BigInt**: Toom-3 (`bigint_mul_toom3.cpp`), Burnikel-Ziegler (`bigint_div_burnikel_ziegler.cpp`), Lehmer GCD, Miller-Rabin basi fisse. *(commit 70a1e52, 4131398)*

**Algebra univar**: Half-GCD (`polynomial_half_gcd.cpp`), Berlekamp, sparse interpolation, cyclotomic detection >724 *(commit dd27a51, T-022)*.

**Algebra multivar/ext**: Wang EEZ + Kronecker fallback, Zippel sparse GCD CRT+Farey *(commit ad023bf, T-006)*, **CRT modular resultant** (`polynomial_resultant_crt.cpp`, Collins 1971 — *era erroneamente dato aperto come T-012*), FGLM + F5 Groebner, tower 2-level + nested squarefree.

**LinAlg**: **Householder QR simbolico razionalizzato** (`matrix_qr.cpp` — *PLAN_HP_PRIME lo dava aperta-perm/MGS, è stale: codice ha risolto l'AST-explosion via reflector `I−2vvᵀ/N` senza √ al denominatore*), Cholesky LDLᵀ, fresh-symbol pervasivo, adaptive G7/K15 quadrature *(commit 5268e90)*.

**Calculus**: Kovacic Case 1 (Laurent √r poli≥4, *commit 709a437, T-008*) + Case 2; Risch trial-constants rimossi (`integrate_risch.cpp:69`), Rothstein-Trager generalizzato, Hermite con D arbitrario; Zeilberger/Gosper/Petkovšek higher-order *(commit 8440ba5, 3addb54, T-027)*; ∫xᵏ/√(c−dx²) *(commit 830f80e, T-055)*, ∫1/√(Ax²+Bx+C) *(6c5ae48)*, inverse-trig IBP, ∫log(x+√(x²+a)) *(ec9e065)*; **limit signed-infinity per poli non-polinomiali a punto finito** (`tan(x)/(x−π/2)→−∞` via sign del reciproco, `try_signed_pole_via_reciprocal`, *commit fc537e5*).

**Special functions (F7.5.E1 — Γ/B/ζ/erf identities)**: ✅ **DONE** (verificato golden 2026-06-28: area `special_fn` = **75/0/5 = 100%**, spec target era ≥82%). ζ pari/dispari via **Bernoulli** (no tabella, REGOLA ZERO, `simplify_special_fn.cpp:380-413`), Γ half-integer + functional-eq + poli→ComplexInfinity, Beta→Γ, Digamma/Polygamma/Pochhammer, erf/erfc. *La spec `Special_Fn_Identities.md` era erroneamente citata da A7 (Meijer-G) — sono task diversi.*

**Complex/branch-cut**: BC-1..BC-3 + BC-1b corrections, ln(a+bi) esatto *(commit 89eae9e, T-017)*, cyclotomic↔exp equal *(4597332, T-025)*, extended-real arithmetic *(b145cb9)*, UnwindingNumber builtin.

**Numerica**: fsolve tolerance ctx-param (`fsolve_tolerance_bits`, T-019), solve_inequality via Sturm, RootOf isolation bounds, Bessel identities *(ade20ff)*, statistics package (`src/statistics/`).

**Infra**: **Phase-7 anti-monolith DONE** — tutti i 23 file ex >500 LOC ora <500 *(commit 7aa5f2e + T-046..050, T-029..051)*; flaky cos(7π/16) fixed via co-function *(d03ff6d, T-004)*; rapidcheck property tests; golden aggregate **94.5%** (F7.5); **A2 hard-timeout RDE solver interrompibile** *(commit 388f6f5)*; **untracked top-level `_deps/` cruft + gitlink rotto rimossi, `_deps/` gitignored** — il dep tree reale è `build/_deps` da FetchContent *(2026-06-28)*; `ledger_index.py` 3-tier classifier + `doctor` *(03928f7)*.

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

> **Avanzamento 2026-06-28**: A2 ✅ (388f6f5) · A10 ✅ chiuso con probe (156b18a) · A16 chiarito = deferral ledgered post-parità (non E1) · A4/A12/A14/A15/A17/A24 = già FATTO/stale. **Quick-win esauriti** (done-stale o ledger-deferred). **Prossima direzione = research C4**: entry-point **A7** (unblocked, vedi RESTART GUIDE in A7); A1 resta BLOCCATO (richiede prima la torre-test deg_t(f)>0).
>
> ⬇️ ordine storico (pre-sessione 2026-06-28), conservato per contesto:

1. ~~**A2**~~ ✅ + **A3** (perf-hang C4, hard)
2. Debiti rapidi E1: ~~A16~~ (deferral) · ~~A17/A24~~ (done) · A18 (ledger-deferred) · **A23** (structured Unimplemented, S1, broad)
3. **A5** (Frobenius, residuo by-design) · **A13** (residue, residuo by-design) · ~~A4~~ (done)
4. **A7** (Slater/Meijer) ← **prossima research, entry-point** · A1 (Risch RP-2) BLOCCATO
5. **A6** (Stauduhar) · **A9** (tower≥3) · ~~A10~~ ✅ · A11 (FROZEN) / A12 (done) / A14 (done)
6. **A15** (done) · **A8** (Kovacic n=12) · **A19-A22** — bassa priorità



---

## F — ESTENSIONI ARCHITETTURALI XCAS / GIAC (lungo termine)

> Macro-aree per il raggiungimento della parità applicativa completa con XCAS/Giac.
> Ogni task richiede la redazione e approvazione di una spec formale prima dell'implementazione (REGOLA 0.1).

| ID | Area | Titolo / Descrizione Algoritmica | Codici | Dipendenze | Stato |
|---|---|---|---|---|---|
| XCAS-EXT-01 | Scripting | **Interprete e Control Flow**: execution engine su AST per blocchi `if/then/else`, cicli `for`/`while`, scope locale/globale e definizioni `proc(args)`. | `[E4·C3·S3·R1]` | Lexer/Parser | Pianificato |
| XCAS-EXT-02 | Calculus | **Trasformate Fourier & Z**: integratore di kernel trascendente ($e^{-i\omega t}$ / $z^{-n}$), tabelle identità, proprietà di convoluzione e traslazione. | `[E3·C3·S3·R1]` | Calculus core | Pianificato |
| XCAS-EXT-03 | Tensor | **Calcolo Tensoriale & Geometria Diff**: oggetti n-dimensionali, indici covarianti/controvarianti, contrazione di Einstein, simboli di Christoffel e operatori vettoriali ($\nabla \times, \nabla \cdot$). | `[E4·C4·S3·R1]` | LinAlg/Symbolic | Pianificato |
| XCAS-EXT-04 | Optimization | **Ottimizzazione & Simplesso**: algoritmo del Simplesso esatto (su Rational) e numerico, solutore moltiplicatori di Lagrange per minimi/massimi vincolati. | `[E3·C3·S2·R1]` | LinAlg esatta | Pianificato |
| XCAS-EXT-05 | Graphics | **Motore Plotting Adattivo**: campionamento adattivo di funzioni cartesiane/parametriche/polari 2D e superfici 3D, generatore di maglie (mesh) agnostico da UI. | `[E4·C3·S2·R1]` | Evaluator numerico | Pianificato |
| XCAS-EXT-06 | Geometry | **Geometria Dinamica & Vincoli**: primitive algebriche (rette, coniche), intersezioni esatte e risolutore di vincoli (incidenza, tangenza) via Basi di Gröbner. | `[E5·C4·S2·R1]` | Gröbner/Resultanti | Pianificato |

---

## Vincoli invariabili (CLAUDE.md)

BigInt only (no int64/double nel core) · Structural sharing · AstArena bump · LPO orientation
· `Result<T>` (no throw/catch) · 500 LOC/file (hard 550) · fresh-symbol pervasivo · REGOLA ZERO
(no via facile, no hardcode non-ledgered) · REGOLA 0.1 (spec read first) · REGOLA 0.2 (no test disable)
· timeout test espliciti (`scripts/test_quick.sh`).
