# PIANO COMPLETO CHIUSURA GAP F3 — F8

**Data**: 2026-06-12
**Scopo**: piano dettagliato, matematicamente corretto, **senza MVP / senza debiti / senza hardcode**, per portare ogni modulo del CAS Engine al traguardo di copertura ≥ 80 % (target progetto). Ogni task riportato qui rispetta REGOLA ZERO + Categorie 1-10 HARDCODE di `CLAUDE.md` e la procedura REGOLA 0.1 di lettura spec.

Ordine di esecuzione: per **dipendenza algoritmica**, non per facilità. Un task la cui dipendenza interna non è ancora pronta **NON va iniziato finché la dipendenza non è chiusa** (REGOLA ZERO: costruire prima il prerequisito).

---

## Stato baseline reale (2026-06-12)

| Modulo | Copertura reale | Bugie smentite | Gap restante |
|---|---|---|---|
| F1 Foundation | 95 % | HPP-F1.1-MUL (Schönhage-Strassen) e HPP-023 (Burnikel-Ziegler) sono APERTE PERMANENTI, non DONE. | SS FFT + Burnikel-Ziegler |
| F2 Univariate | 98 % | Nessuna | — |
| F3 Multivar/Ext | ~75 % | F3.2 / F3.3 / F3.4 / F3.5 chiusi 2026-05-29 → 31. F3.5-DEBT-01 chiuso 2026-05-31. F3.4-DEBT-01 RESIDUO aperto su Q(β)[x] reducible + multi-β nesting > 1. | Stauduhar ≥6, EEZ Wang fallback Kronecker, Zippel sparse interp, F3.4-DEBT-01 residuo |
| F4 LinAlg | ~85 % | HPP-F4.1-QR-HOUSEHOLDER APERTA PERMANENTE (codice usa MGS, piano diceva Householder). | Householder simbolico stabile, eigenvalues/SVD su campi di funzioni |
| F5 Calculus | ~50 % | HC-KV-02 Laurent √r APERTA, HC-KV-03 Case 3 INCOMPLETO (Case 2 chiuso), HPP-007 trial constants risolto. Risch RDE Bronstein framework presente ma 1 test fallisce (erf×exp). | Puiseux/Laurent generico, Kovacic Case 3, ODE sistemi su campi di funzioni, exp-fold simplifier rule (debito sessione 2026-06-12) |
| F6 Numeric/CAD | ~35 % | Collins CAD assente, branch-cut parziale (`branch_cut_aware_logexp` opt-in solo per `ln(exp(z))`). | Collins CAD completo (proiezione McCallum + lifting), branch-cut esteso |
| F7 Acceptance | ~70 % | Slater / Meijer G assenti. | Slater transforms, Meijer G fallback |
| F8 Modularità | ~90 % | HC-F8-MONOLITH-WAIVER 28 file > 500 LOC. | Split monoliti rimasti |

---

## Principi metodologici applicabili a ogni task

1. **REGOLA 0.1 lettura spec**: per ogni task elencato qui esiste già il file `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<Spec>.md` oppure ne va scritto uno **prima** di toccare il codice. Vedi colonna "Spec source" per ogni riga.
2. **Anti-hardcode test bench**: per ogni costante numerica introdotta, eseguire le 4 domande di self-check (CLAUDE.md §Regola di Self-Check) e iscriverla in `HARDCODE_LEDGER.md` se non passa.
3. **Aritmetica esatta only**: zero `int64_t` / `double` nel core simbolico. Solo `BigInt`, `Rational`, `ComplexRational`, `AlgebraicNumber*`.
4. **Structural Sharing**: ogni funzione deve restituire l'`ExprPtr` originale se non modifica. Vietato clone non necessario.
5. **Anti-monolite**: ≤ 500 LOC/file. Split in moduli specializzati prima di crescere.
6. **TDD**: test prima di codice (RED→GREEN→REFACTOR). Coverage ≥ 80 % per riga e per branch sul nuovo codice.
7. **Verification gate**: `bash scripts/test_quick.sh` (≤ 600 s) PASS + zero regressioni dopo ogni task. `--slow` (≤ 1800 s) come gate pre-PR.
8. **Spec-driven docstring**: ogni nuova funzione cita teorema/algoritmo + sezione della reference (Bronstein, Cohen, Trager, Stauduhar, Schönhage-Strassen, Wang, Collins, Bronstein 2005).

---

# F1 — Chiusura debiti permanenti (residuo 5 %)

## F1.1 — Schönhage-Strassen NTT multiplication (HPP-F1.1-MUL)

| Campo | Valore |
|---|---|
| File da creare | `src/foundation/bigint_mul_ssa.cpp`, `src/foundation/bigint_mul_ssa_internal.hpp`, `test/unit/foundation/test_bigint_mul_ssa.cpp` |
| File da modificare | `src/foundation/bigint_arithmetic.cpp:231` (switch a SSA per `n_limbs ≥ kSSAThreshold`), `include/cas/bigint.hpp` (eventuale config) |
| Spec source | DA SCRIVERE → `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Schonhage_Strassen_FFT.md` (Brent-Zimmermann §1.3.5; Schönhage-Strassen 1971; Crandall-Pomerance §9.5.7) |
| Algoritmo | NTT modulo primi di Fermat F_k = 2^(2^k)+1, k ∈ {5,6,7,8}. Radix decimation-in-frequency, Schönhage-Strassen recursive doubling per gestire la "negacyclic convolution". Modular reduction esatta via shift+sub (Solinas trick). |
| Costanti vietate | Nessun `kFFTThreshold` magico: deve essere derivato da bound asintotico O(n log n log log n) vs O(n^1.585). Calcolo empirico **una sola volta** in CI per benchmark, poi `ctx.bigint_ssa_threshold()` configurabile con default = risultato del benchmark. |
| Test obbligatori | Cross-check vs Karatsuba/Toom-3 per n ∈ {4096, 8192, 16384, 32768, 65536} limbs. Identità su 1000 inputs random. Stress 2^20 limbs sotto ASan. |
| Acceptance | Suite quick PASS, benchmark `scripts/benchmark.sh` su `MulHuge` mostra speed-up > 2× rispetto a Karatsuba per n ≥ 8192 limbs. Ledger HPP-F1.1-MUL chiuso. |
| Dipendenze | Nessuna (autonomo). |
| Effort stimato | 2-3 settimane T3-Opus (ricerca + impl + tuning). |

## F1.2 — Burnikel-Ziegler quotient (HPP-023)

| Campo | Valore |
|---|---|
| File da creare | `src/foundation/bigint_div_burnikel.cpp` |
| Spec source | DA SCRIVERE → `Burnikel_Ziegler_Division.md` (Burnikel-Ziegler "Fast Recursive Division" 1998 §3-5; Brent-Zimmermann §1.4.3.2) |
| Algoritmo | Divisione D&C: `div_3by2` recursive sub-routine. Threshold derivato dal cross-over asintotico O(M(n) log n) vs O(n²). |
| Costanti vietate | `kBZThreshold` configurabile via `CASContext`. |
| Test | Identità `q*v + r = u` per random inputs su 4096-65536 limbs. Edge cases: divisore monolimb, divisore = u. |
| Acceptance | Quoziente bit-identico a Knuth Algorithm D. Speed-up > 3× per ≥ 8192 limbs. Ledger HPP-023 chiuso. |
| Dipendenze | F1.1 utile ma non bloccante (BZ funziona anche con Karatsuba). |
| Effort | 1-2 settimane. |

## F1.3 — Lehmer GCD multi-limb (HPP-019)

| Campo | Valore |
|---|---|
| File | `src/foundation/bigint_gcd_lehmer.cpp` (estendere) |
| Spec source | DA SCRIVERE → `Lehmer_GCD_Multilimb.md` (Knuth TAOCP vol. 2 §4.5.2, Algoritmo L) |
| Stato attuale | Single-limb surrogate, breaks su `na != nb`. |
| Fix | Estensione multi-limb con accurate cofactor recovery via matrice 2×2. |
| Acceptance | GCD identico a Euclid per random pairs fino a 65536 limbs. Ledger HPP-019 chiuso. |
| Effort | 3-5 giorni. |

---

# F3 — Chiusura ~25 % gap residuo

## F3.A — Stauduhar Galois group identification ≥ degree 6

| Campo | Valore |
|---|---|
| Spec source | ESISTE → `Galois_Groups.md` (11.8 KB, contiene tabella resolventi per S_n / A_n / sottogruppi transitivi fino n=12) |
| File da creare | `src/algebra/galois_stauduhar.cpp`, `src/algebra/galois_stauduhar_resolvents.cpp` (split anti-monolite), `src/algebra/galois_stauduhar_internal.hpp`, test in `test/unit/algebra/test_galois_stauduhar.cpp` |
| File da modificare | `src/algebra/galois.cpp:202-242` (delega a Stauduhar per deg ≥ 6) |
| Algoritmo | Stauduhar 1973 (Math. Comp. 27): resolventi numeriche per ogni sottogruppo H ⊂ S_n, isolamento radici via Sturm-bigfloat (F5.3 DONE), test sull'integralità della resolvente. Per ogni livello di lattice, scelta del successivo H da `transitive_subgroups[n]` (table from Hulpke-Hulpke 1997). |
| Costanti vietate | Tabella `transitive_subgroups` derivata da paper Conway-Hulpke-McKay (referenza: LMFDB) — NON è hardcode perché è un dato matematico immutabile (gruppi finiti sono unici). Esposta come `const`-data file. Precisione bigfloat per la resolvente = `ctx.galois_resolvent_precision_bits()` (default 256), backoff exponential se la resolvente cade su intero ambiguamente vicino a 0. |
| Casi test (≥ 30) | Polinomi noti: x^6 + 3 (S_6), x^6 + x^5 + x^4 + x^3 + x^2 + x + 1 (C_6), x^7 − 7x + 3 (PSL(2,7)), x^8 + 1 (C_2³), x^10 + … (Mathieu M_10). Vedi `Galois_Groups.md §Esempi di validazione`. |
| Anti-hardcode | Tabella resolventi cita Hulpke 2005, NON inserita come `if (poly==...) return "C2";`. Risolutori istanziati via builder generico. |
| Acceptance | ≥ 95 % accuracy su benchmark Hulpke (50 polinomi degree 6-10). Suite ≥ baseline. |
| Dipendenze | F5.3 Sturm-bigfloat (CHIUSO). |
| Effort | 3-4 settimane T3. |

## F3.B — Wang EEZ Hensel-fallback Kronecker per low-degree

| Campo | Valore |
|---|---|
| Spec source | ESISTE → `Hensel_Lifting.md` |
| File da modificare | `src/algebra/factorization_wang_eez.cpp` (entry-point), nuovo `src/algebra/factorization_kronecker_fallback.cpp` |
| Algoritmo | Dispatcher: se Hensel lifting fails dopo `ctx.max_hensel_lift_attempts()` con `bad-prime-rate > 0.5`, fall through a Kronecker-substitution (Knuth TAOCP §4.6.2) per deg ≤ 8. Per deg > 8 Kronecker è O(n^n) → bail-out `Unimplemented` esplicito. |
| Costanti vietate | Threshold `kronecker_max_degree` configurabile (default 8, derivato da Kronecker complexity 2^deg interpolations). |
| Test | Wang corner cases che fanno fallire Hensel (lc divisibile da tutti primi small): `(2x²+1)(3x²+1)` mod p per p∈{2,3,5,7,11} corner. |
| Acceptance | Coverage Wang da 65 % → 90 %. Suite ≥ baseline. |
| Dipendenze | Nessuna. |
| Effort | 2 settimane. |

## F3.C — Zippel sparse interpolation per GCD multivariate

| Campo | Valore |
|---|---|
| Spec source | ESISTE → `Zippel_Sparse_Interpolation.md` (12.4 KB) |
| File da creare | `src/algebra/polynomial_zippel_sparse.cpp` |
| Algoritmo | Zippel 1979: skeleton + sparse multivariate interpolation. Schwartz-Zippel error probability `δ = (deg(f) / |F|)^k` con `k = ctx.zippel_confidence_samples()`. |
| Costanti vietate | `+8` campioni extra (vedi HARDCODE_LEDGER multivariate GCD): RIMUOVERE, sostituire con `ceil(log(δ)/log(1-p_hit))` come da Categoria 2 esempio. |
| Acceptance | GCD multivariato 3+ variabili sparse: speed-up ≥ 5× vs dense. Suite ≥ baseline. |
| Effort | 2 settimane. |

## F3.D — F3.4-DEBT-01 residuo (multi-β nested + reducible Q(β)[x])

| Campo | Valore |
|---|---|
| File | `src/algebra/algebraic_tower_primitive_nested.cpp` (estendere) |
| Spec source | Cohen §3.6.4 (multi-β iteration), GTM 138 |
| Algoritmo | Iterazione su tower di β annidati: per ogni livello, ricerca shift di Trager con squarefree-test. Selezione del fattore irriducibile su Q(β)[x] via test conjugate-realization. |
| Test | Q(√2, √(2+√3)) — caso classico nested. |
| Acceptance | Residuo OPEN chiuso. Test PrimitiveElement nested multi-livello PASS. |
| Effort | 1-2 settimane. |

## F3.E — Modular Resultants CRT (Task 5.1 piano user)

| Campo | Valore |
|---|---|
| Spec source | DA SCRIVERE → `Modular_Resultants_CRT.md` (basato su Collins 1971, Hadamard bound) |
| File | `src/algebra/polynomial_resultant_modular.cpp` |
| Algoritmo | Multi-modulo, ogni primo single-precision, recombination via CRT incrementale, terminazione su Hadamard bound. Già infrastruttura `polynomial_modular.cpp` disponibile. |
| Costanti vietate | Bound CRT calcolato dinamicamente da Hadamard, non hardcoded. Primi sorgenti pool di Mersenne+Solinas (no `p=13` fallback). |
| Acceptance | Resultanti su torri ≥ 3 livelli `Q(α₁, α₂, α₃)[x]` chiudibili in tempo polinomiale. |
| Dipendenze | F3.D utile. |
| Effort | 2 settimane. |

---

# F4 — Chiusura ~15 % gap (Householder + advanced linalg)

## F4.A — Householder QR simbolico stabile (HPP-F4.1-QR-HOUSEHOLDER)

| Campo | Valore |
|---|---|
| File | `src/linalg/matrix_qr_householder.cpp` (nuovo, separato da MGS attuale) |
| Spec source | ESISTE → `Matrix_Adapter.md` (estendere) o DA SCRIVERE `Householder_Symbolic_Stable.md` |
| Problema noto | Householder simbolico esplode l'AST per matrici ≥ 8×8 a causa di radici quadrate al denominatore. |
| Algoritmo corretto | Householder con **rationalized reflector** (lavoro F4 Task 2.1 già fatto per piccole matrici, generalizzare): mantenere `α = ‖x‖²` invece di `√‖x‖²` finché possibile, applicare il riflettore in forma `H = I − (2/α)·v·vᵀ` dove `α = vᵀv`. Niente `√` nel denominatore durante updates di Q e R. Solo al termine, se richiesto Q ortogonale esatto, applicare la `√` (oppure restituire QR fattorizzato senza estrarre `√α`). |
| Anti-hardcode | Nessun threshold size — algoritmo deve scalare su qualsiasi n. |
| Test | Householder vs MGS su matrici simboliche 8×8, 12×12, 16×16 con elementi Q(x,y). Verifica `Q·R = A` e `QᵀQ = I` post-`√` extraction. |
| Acceptance | Householder converge senza esplosione AST per n ≤ 16 con elementi razionali in 2 variabili. Ledger HPP-F4.1 chiuso. |
| Effort | 2-3 settimane. |

## F4.B — Eigenvalues / Jordan canonical form su campi di funzioni

| Campo | Valore |
|---|---|
| File | `src/linalg/matrix_eigenvalues_funcfield.cpp` |
| Algoritmo | Polinomio caratteristico via Faddeev-LeVerrier (già parziale), fattorizzazione su Q(x,...) via `factor_multivariate` (F3 ready), RootOf nodes (F5.4 DONE), Jordan blocks calcolati simbolicamente. |
| Spec source | DA SCRIVERE → `Eigenvalues_FuncField.md` (Gantmacher Matrix Theory vol.2 §4-6) |
| Dipendenze | F3.D, F5.4 (entrambi DONE/in-progress). |
| Effort | 2 settimane. |

## F4.C — SVD esatto / Schur decomposition su campi numerici esatti

| Campo | Valore |
|---|---|
| Stato | Out of scope F8 (richiede sviluppo nucleo nuovo). Iscrivere come **Task F9** futuro. |

---

# F5 — Chiusura ~30 % gap (Calculus completion)

## F5.A — Motore Puiseux/Laurent series generico (sblocca HC-KV-02)

| Campo | Valore |
|---|---|
| File da creare | `src/series/puiseux_engine.cpp`, `src/series/puiseux_arithmetic.cpp`, `src/series/puiseux_expand_algebraic.cpp`, `include/cas/puiseux.hpp` (nuovo namespace `cas::series`). Test in `test/unit/series/`. |
| Spec source | DA SCRIVERE → `Puiseux_Series_Engine.md` (Walker "Algebraic Curves" §IV.3, Brieskorn-Knörrer §8, Bronstein "Symbolic Integration I" §3.4 per applicazione a √r). |
| Algoritmo | Newton polygon per fattorizzazione di Puiseux su Q[[x]][y]. Coefficiente leader via radice n-esima di un razionale. Espansione locale `f(x) = Σ c_k · x^(k/n)` con `n = ramification index` e `c_k ∈ Q̄`. Supporto a poli (ordine negativo). |
| Anti-hardcode | Truncation order `ctx.puiseux_truncation_order()` configurabile (default = `2·v` dove v = ordine del polo, parametrizzato sul caso d'uso). |
| Test | Sviluppo locale di √(x³ − x), Bessel ODE invariant (radice quadrata di razionale con polo di ordine 4), funzione ellittica. |
| Acceptance | API `puiseux_expand(f, x0, n_terms)` con accuratezza esatta su 30+ casi noti. |
| Dipendenze | F3.D (algebraic numbers Q(α)). |
| Effort | 4-6 settimane T3 (algoritmo centrale, ricco di edge cases). |

## F5.B — Kovacic Case 1 completo (chiude HC-KV-02)

| Campo | Valore |
|---|---|
| File | `src/calculus/ode_kovacic_case1.cpp` (estendere `case1_omega`) + `src/calculus/ode_kovacic_laurent.cpp` (nuovo, glue verso F5.A) |
| Algoritmo | Kovacic 1986 §3 Thm 1.2: per polo c di ordine 2v, espansione Laurent di √r = Σ s_k·(x−c)^(k−v) con `s_k = (u_k − Σ_{i=0..k-1} s_i·s_{k−i}) / (2·s_0)`. Calcolo dei termini polari α_c± dal coefficiente s_v. Identica logica per parte polinomiale al ∞ (degree pari ≥ 2). |
| Dipendenze | F5.A obbligatoria. |
| Acceptance | Bessel ODE deg 1 risolta in closed form, Weber/parabolic-cylinder casi, Riemann curve quartica. HC-KV-02 chiuso. |
| Effort | 1 settimana (dopo F5.A). |

## F5.C — Kovacic Case 3 (chiude HC-KV-03)

| Campo | Valore |
|---|---|
| File | `src/calculus/ode_kovacic_case3.cpp` (nuovo) |
| Spec source | DA SCRIVERE → `Kovacic_Case3_SL2C.md` (Kovacic 1986 §4-5; classificazione sottogruppi finiti di SL(2,ℂ): A_4, S_4, A_5) |
| Algoritmo | Per ogni gruppo finito di SL(2,ℂ) (3 classi), test sui poli di r e sul comportamento all'∞ con tabelle di Klein. Ricostruzione della invariante via fattorizzazione di un polinomio di grado 4, 6 o 12 con coefficienti razionali. |
| Anti-hardcode | Tabelle di Klein per i sottogruppi: dati matematici fissi (rappresentanti minimal). Iscritti come `const` con citazione paper. NON hardcode (Categoria "eccezione legittima costanti matematiche esatte"). |
| Test | ODE Schwarz triangle equation per ogni classe (A_4, S_4, A_5). |
| Acceptance | HC-KV-03 chiuso. |
| Dipendenze | F5.A, F5.B, F3.A (Galois). |
| Effort | 3-4 settimane. |

## F5.D — Risch DE Bronstein chiusura `IntegralOfErfTimesExp`

| Campo | Valore |
|---|---|
| File | `src/symbolic/simplify_arithmetic_chain.cpp` (Product simplifier) **+** `include/cas/cas_context_params.hpp` (nuovo flag) |
| Problema | IBP verifier non chiude perché Product simplifier non folda `exp(a)·exp(b) → exp(a+b)`. Tentativo sessione 2026-06-12: fix simplifier rompe MRV / Gruntz / Weierstrass. |
| Strategia corretta | Introdurre `ctx.fold_exp_products` (default `false`). IBP verifier lo abilita temporaneamente con scoped guard. MRV / Gruntz mantengono `false`. Necessità: refactor `CASContext::simplify` per accettare `SimplifyHints { bool fold_exp_products; }`. |
| Anti-hardcode | Flag esposto come parametro contestuale, no costanti. |
| Acceptance | `IntegrateRischExpMixTest.IntegralOfErfTimesExp` PASS. Suite ≥ baseline (nessuna regressione MRV/Gruntz/Weierstrass). |
| Dipendenze | Nessuna. |
| Effort | 3-5 giorni. |

## F5.E — Risch transcendental Cap.8 corpus 60 % → 75 %

| Campo | Valore |
|---|---|
| File | `src/calculus/risch_rde_bronstein.cpp` (estendere `solve_risch_de_field` per df > 0 nel parametric case) |
| Spec source | ESISTE → `Risch_Transcendental_Cap8.md` |
| Stato attuale | Framework completo, `solve_risch_de_parametric_field` ritorna `Unimplemented` per df > 0. |
| Algoritmo | Bronstein 6.5 PolyRischDE: bound grado esplicito per df > 0; risoluzione coefficienti via Hermite reduction parametrica. |
| Acceptance | ≥ 20 casi corpus Bronstein cap. 6-8 PASS. Trial constants già rimossi. |
| Dipendenze | F5.D utile. |
| Effort | 2-3 settimane. |

## F5.F — ODE sistemi lineari su campi di funzioni Q(x)

| Campo | Valore |
|---|---|
| File | `src/calculus/ode_system_linear_funcfield.cpp` |
| Spec source | DA SCRIVERE → `ODE_System_Linear_FuncField.md` (Coddington-Levinson "Theory of ODEs" §III, Putzer algorithm for matrix exponentials) |
| Algoritmo | y' = A·y con A ∈ M_n(Q(x)). Jordan canonical form di A (F4.B), back-substitution termine per termine via solver Risch-DE scalare. |
| Anti-hardcode | Nessun threshold sulle dim della matrice (deve scalare). |
| Acceptance | ODE 2×2 e 3×3 con coefficienti razionali risolte. |
| Dipendenze | F4.B + F5.E. |
| Effort | 2 settimane. |

---

# F6 — Chiusura ~50 % gap (Numeric + CAD)

## F6.A — Collins CAD: proiezione di McCallum

| Campo | Valore |
|---|---|
| File | `src/algebra/cad_collins_projection.cpp`, `src/algebra/cad_collins_lifting.cpp`, `src/algebra/cad_collins_internal.hpp` |
| Spec source | DA SCRIVERE → `Collins_CAD.md` (Collins 1975 LNCS 33, McCallum 1998 JSC 24, Caviness-Johnson textbook §6) |
| Algoritmo | Projection chain: per ogni livello k → k-1, calcolo del set di proiezione (coefficienti leader, discriminanti, resultants). Lifting: per ogni cell della partizione (k-1), root isolation in `Q[α_1,...,α_{k-1}][x_k]` via Sturm-bigfloat. |
| Anti-hardcode | Numero di cell `ctx.cad_max_cells()` configurabile. Precisione root isolation = `ctx.cad_isolation_bits()` (default 256). |
| Test | Decomposizione cilindrica di `x²+y²−1 < 0`, `xy − 1 > 0`, sistema polinomiale 3-variabili. Confronto contro QEPCAD se installato. |
| Acceptance | CAD per polinomi grado ≤ 4 in 3 variabili risolto correttamente con cell count esatto. |
| Dipendenze | F3.E (resultants), F5.A (root isolation). |
| Effort | 6-8 settimane T3. |

## F6.B — Branch-cut completion: `sqrt(x²)=|x|`, `pow(x,1/2)·pow(x,1/2)=x`

| Campo | Valore |
|---|---|
| File | `src/symbolic/simplify_branch_cut.cpp` (estendere) |
| Spec source | DA SCRIVERE → `Branch_Cut_Propagation.md` (Corless-Davenport-Jeffrey "According to Abramowitz and Stegun" 2000) |
| Algoritmo | Tabella completa di propagazione branch-cut per le funzioni multivalore (sqrt, pow rationale, log, arctan2). Uso di `UnwindingNumber` (F6 Task 6.1 DONE) per discriminare il foglio. |
| Acceptance | `sqrt((-1)²)` → `1` solo se x è reale; rimane `sqrt(x²)` con annotazione altrimenti. `pow(x,1/2)·pow(x,1/2)` → `|x|` per x reale, `x·K(x)` altrimenti. |
| Dipendenze | Task 6.2 DONE. |
| Effort | 2-3 settimane. |

## F6.C — Sturm-bigfloat sostituzione completa double in fsolve

| Campo | Valore |
|---|---|
| File | `src/algebra/fsolve.cpp` (HPP-006 in ledger) |
| Algoritmo | Migrare Newton polishing da `double` a `BigFloat` con `ctx.fsolve_tolerance_bits()`. |
| Anti-hardcode | Tolerance già parzialmente configurabile, rendere `bits` non `double`. |
| Acceptance | HPP-006 chiuso. |
| Effort | 3-5 giorni. |

## F6.D — Adaptive numerical integration robusta (esiste spec)

| Spec | `Adaptive_Numerical_Integration.md` (7.8 KB) |
| File | `src/numeric/adaptive_integration.cpp` |
| Algoritmo | Gauss-Kronrod 7/15 con error estimate, recursive subdivision, Wynn epsilon acceleration. |
| Acceptance | Convergenza su integrali oscillanti / quasi-singolari. |
| Effort | 2 settimane. |

---

# F7 — Chiusura ~10 % gap (Acceptance)

## F7.A — Slater transforms per pFq → Meijer G

| Spec | ESISTE → `Special_Fn_Identities.md` (5.2 KB) |
| File | `src/symbolic/special_meijer_slater.cpp` |
| Algoritmo | Erdélyi-Slater tabella trasformazioni: 16 identità canoniche tra pFq di parametri specifici e Meijer G. Implementazione dispatcher pattern-driven con citazione paper per ogni regola. |
| Anti-hardcode | Tabella delle trasformazioni = dati matematici. Esposta come `const` con citazione (NON Categoria 8 perché ogni riga è un teorema dimostrato). |
| Test | Tutti 16 cases Erdélyi-Slater + 4 corner di Bailey. |
| Effort | 3-4 settimane. |

## F7.B — Meijer G fallback integrator

| File | `src/symbolic/special_meijer_g.cpp`, wiring in `src/calculus/integrate_core.cpp` (passo 5 nuovo) |
| Algoritmo | Wiring nel motore di integrazione globale come fallback dopo Risch. Tabella Erdélyi vol.1 §5.6 mappata strutturalmente, NON via stringhe. |
| Effort | 3 settimane. |

## F7.C — Bessel identities (esiste spec)

| Spec | `Bessel_Identities.md` (3.9 KB) |
| File | `src/symbolic/simplify_bessel_orthogonal.cpp` (estendere) |
| Effort | 1-2 settimane. |

## F7.D — Probability/Statistics module (esiste spec)

| Spec | `Probability_Statistics.md` (9.4 KB) |
| File | `src/symbolic/probability_basics.cpp` + `src/symbolic/probability_distributions.cpp` |
| Effort | 4-6 settimane. Iscrivere come parte F9 se priorità F1-F7 prevale. |

---

# F8 — Modularizzazione finale (chiusura HC-F8-MONOLITH-WAIVER)

| Campo | Valore |
|---|---|
| Stato | 28 file > 500 LOC waiver attivo. |
| Lista file critici | `grep -E ":[0-9]+:" anti-monolith report` |
| Strategia | Split per area logica matematica (es. `simplify_arithmetic_chain.cpp` → `_product.cpp` + `_step5.cpp` + `_step8.cpp`). |
| Anti-debito | NIENTE `#include "old_file.cpp"`. Solo split modulare con header internal limitato al namespace `detail`. |
| Acceptance | 0 file > 500 LOC (no whitelist). Ledger HC-F8-MONOLITH-WAIVER chiuso. |
| Effort | 1-2 settimane (pure refactor + test invariati). |

---

# Cycle-guard / Audit residuali (HPP-F75-AUDIT-CYCLE-GUARD-1/2/3)

Tutti e tre sono soglie `kMaxAppendDepth = 1024`, `kVisitRecursiveMaxDepth = 4096`, `kGrowthRankMaxDepth = 1024`. Sono **safety net** per loop infiniti nel motore MRV / differential field / Gruntz comparator. CLAUDE.md Categoria 1: vanno esposte via `CASContext` con default identici.

| File | Action |
|---|---|
| `src/calculus/limit_mrv_exp.cpp:38` | `ctx.mrv_max_append_depth()` |
| `src/calculus/differential_field.cpp:21` | `ctx.diff_field_max_visit_depth()` |
| `src/calculus/limit_mrv_compare.cpp:100` | `ctx.mrv_growth_rank_max_depth()` |
| Effort | 1 giorno totale (refactor meccanico). |

---

# Cross-cutting: nuovi parametri CASContext da introdurre

Tutti questi parametri DEVONO essere esposti in `include/cas/cas_context_params.hpp` con default documentato e citazione bibliografica:

```
bigint_ssa_threshold              (default = benchmark-derived)
bigint_bz_threshold               (default = benchmark-derived)
galois_resolvent_precision_bits   (default = 256)
hensel_max_attempts               (default = 4)
kronecker_max_degree              (default = 8)
zippel_confidence_samples         (default ceil(log δ / log(1-p_hit)))
puiseux_truncation_order          (default = 2*v)
cad_max_cells                     (default = 0 = unlimited; raise OOM-guard)
cad_isolation_bits                (default = 256)
fold_exp_products                 (default = false)
mrv_max_append_depth              (default = 1024)
diff_field_max_visit_depth        (default = 4096)
mrv_growth_rank_max_depth         (default = 1024)
fsolve_tolerance_bits             (default = 80 = 1e-24)
puiseux_max_branches              (default = computed from degree)
risch_de_max_degree               (default = computed from Bronstein 6.5)
```

---

# Ordine di esecuzione raccomandato (dipendenze critiche)

```
F1.1 SSA  ─┬─> F1.2 BZ ─> F1.3 Lehmer
           └─> (sblocca speed)
F3.D Nested ────────┐
F3.E Modular Resul ─┴───> F3.A Stauduhar
F3.B EEZ-Kronecker  ─────> (autonomo)
F3.C Zippel         ─────> (autonomo)
F5.A Puiseux  ───────────> F5.B Kovacic Case 1 ──> F5.C Kovacic Case 3
                                                    F3.A Galois ─┘ (Case 3 dep)
F5.D Risch IBP ───────────> F5.E Risch parametric df>0
F4.A Householder ────────> F4.B Eigen funcfield ──> F5.F ODE systems
F6.A Collins CAD (deps F3.E + F5.A)
F6.B Branch-cut
F6.C Sturm-bigfloat fsolve
F6.D Adaptive
F7.A Slater ─> F7.B Meijer G fallback
F7.C Bessel  (autonomo)
F8 split (parallelo a tutto)
HPP-AUDIT-CYCLE-GUARD-1/2/3  (1 giorno, fare subito)
```

**Path critico (catena più lunga)**: F3.D → F3.E → F6.A (≈ 10-12 settimane).
**Path indipendente più veloce**: F5.D + cycle-guards + F6.C (≈ 2 settimane).

---

# Effort totale stimato

| Fase | Effort settimane T3-Opus |
|---|---|
| F1 (SSA + BZ + Lehmer) | 5-7 |
| F3 (Stauduhar + EEZ + Zippel + nested + CRT) | 9-11 |
| F4 (Householder + Eigen funcfield) | 4-5 |
| F5 (Puiseux + Kovacic Case 1+3 + Risch parametric + IBP exp-fold + ODE systems) | 13-18 |
| F6 (Collins CAD + branch-cut + fsolve + adaptive) | 11-14 |
| F7 (Slater + Meijer + Bessel) | 7-9 |
| F8 monolith split | 1-2 |
| Cycle-guards | 1 giorno |
| **Totale** | **50-66 settimane T3 cumulative** |

Con parallelismo 2-3 worker T3 indipendenti su path indipendenti: **22-30 settimane wall-clock**.

---

# Gating procedurale per ogni task

Prima di iniziare:
1. `read .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<Spec>.md` (o crearlo) + dichiarazione formale "Ho letto la specifica formale".
2. `bash scripts/test_quick.sh` baseline snapshot.
3. Iscrivere task in `TODO_PH8.md` come `in_progress` con timestamp.
4. Per ogni costante: verifica self-check Categorie 1-10. Iscrivere in `HARDCODE_LEDGER.md` se non passa.

Durante:
5. TDD strict: nuovo test PRIMA di nuovo codice.
6. Spot-check anti-monolite ogni 200 LOC aggiunte.
7. Spot-check anti-double / anti-int64 nel core simbolico.

Prima di dichiarare DONE:
8. `bash scripts/test_quick.sh` ≥ baseline pass count.
9. `bash scripts/test_quick.sh --slow` PASS.
10. `bash scripts/benchmark.sh` zero regressioni performance.
11. `ninja -C build 2>&1 | grep -E "error:|warning:"` zero.
12. ASan + UBSan run su nuovi test.
13. Aggiornare `TODO_PH8.md` con stato onesto.
14. Chiudere voci `HARDCODE_LEDGER.md` risolte con riferimento commit.
15. Aggiornare questo `PLAN_F3_F8_GAP_CLOSURE.md` con avanzamento.

---

# Bugie potenziali da evitare (lezione apprese 2026-06-11)

Una bugia di piano è: dichiarare un task DONE quando l'algoritmo dichiarato nel piano NON corrisponde a quello in codice. Esempio storico: F1.1 piano diceva "Schönhage-Strassen exit gate" ma codice usava Karatsuba fallback. Fix: ogni `[DONE]` mark va validato con `grep` + `Read` che codice corrisponde target piano (memory `plan-lies-2026-06-11`).

Ogni `[DONE]` apposto a un task qui sopra:
- DEVE essere accompagnato da SHA commit specifico.
- DEVE elencare ≥ 1 file:line di riferimento del codice rilevante.
- DEVE confermare passing tests rilevanti con pass-count esatto.

---

# Risk register

| Rischio | Probabilità | Mitigazione |
|---|---|---|
| Collins CAD esplode su input adversariali (cell count exponential) | Alta | OOM-guard cap configurabile + diagnostica esplicita |
| Schönhage-Strassen instabile per limbs intermedi | Media | Cross-check vs Karatsuba in CI per ogni n |
| Puiseux: divergenza con coefficienti algebrici annidati > 2 livelli | Alta | Test edge cases + bail-out diagnostico |
| Kovacic Case 3 Klein tables miscoded | Media | Spec scritta verbatim da paper + property tests |
| Refactor `CASContext::simplify` con `SimplifyHints` rompe API esistente | Media | API esistente preservata via overload, hint default `=SimplifyHints{}` |
| Stauduhar resolvent precision insufficiente per radici quasi-degenerate | Media | Backoff exponential precision + diagnostica |

---

**Fine PIANO F3_F8.**
