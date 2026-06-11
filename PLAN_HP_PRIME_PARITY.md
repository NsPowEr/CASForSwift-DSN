# PIANO IMPLEMENTAZIONE — Parità HP Prime G2

> Obiettivo: portare ogni engine matematico al 100% del dominio CAS production-grade (riferimento: HP Prime G2 / Giac-Xcas baseline).
> Metodo: chiusura stratificata bottom-up. Ogni strato chiuso 100% prima di toccare lo strato superiore. Ogni algoritmo deve avere algoritmo canonico citato, certificatore indipendente, test anti-hardcode, e validation gate.
> Anti-pattern vietati: REGOLA ZERO (via facile), via-passaggio non marcata, "Risolta" su subset non dichiarato, hardcode non in `HARDCODE_LEDGER.md`.
> Terminologia: NO "MVP". Solo `Risolta` (dominio dichiarato 100% chiuso e certificato), `Parziale avanzata` (algoritmo vero ma dominio strict subset esplicito), `Parziale` (componenti mancanti: certificatore/test/fallback), `Aperta` (placeholder o non implementato).

---

## QA REVIEW v2 — Aggiunte / correzioni post-audit qualità

Il piano v1 aveva ottimismi residui e shortcut taciti. Correzioni applicate:

1. **Rimossa mitigazione "wrap GMP"**: viola REGOLA ZERO se non marcata. Performance gate BigInt ricalibrato (vedi F1.1).
2. **Eliminato concetto "MVP"**: ogni status deve essere `Risolta` su dominio dichiarato OPPURE `Parziale avanzata` con dominio strict subset esplicito (nome task dichiara il subset, es. "GCD multivariato lineare certificato").
3. **CAD spostato a Fase 8 post-parità**: complessità doppia esponenziale + scope onesto incompatibile con strato L3. F6 ora copre disequazioni 1-var via Sturm/risultanti come strato vero, no CAD subset.
4. **Aggiunto F1.6 Complex foundation**: i complessi sono trasversali (Q[i], C[x], eigenvalues, branch cuts). Devono stare in L0, non L3.
5. **Aggiunto F0.6 Build infrastructure CI**: sanitizer + bench gate + anti-monolith check + diff-coverage automated.
6. **Aggiunto F0.7 Doc-per-algorithm template**: ogni algoritmo Risolta deve avere `.APROJECT_REFERENCES/algorithms/<nome>.md` con: algoritmo canonico, complessità, dominio, limiti, riferimenti bibliografici, esempi.
7. **Aggiunta API stability freeze**: a chiusura ogni fase, public API congelata. Strato successivo non rompe ABI.
8. **Aggiunta cross-phase regression policy**: se Fn rivela bug strutturale in F_{n-k} foundation, ROLLBACK a F_{n-k} prima di proseguire. No layering di shortcut.
9. **Aggiunto pre-review interlock**: prima di audit T3-Opus, un subagent T2-Sonnet-thinking pre-review il blocco appena chiuso, cattura ~50% problemi triviali, riduce token T3-Opus.
10. **Aggiunte voci `Aperta` esplicite permanenti** per: Risch structure theorem full Bronstein cap 9, Galois ≥6 generale, FFT BigInt competitiva GMP, Hypergeometric `_pF_q` recognition completo, CAD generale. Documentate in F7 acceptance: parità HP Prime ≥95% NON significa 100% di questi target research.
11. **Golden corpus aumentato da 500 a 2000 input** (50/area era debole; HP Prime test suite reale è 5000+; 2000 è bilanciato).
12. **Performance gate BigInt ricalibrato**: target onesto **30-50% performance GMP** (non 80%). GMP = 30 anni tuning; 80% richiederebbe via-facile wrap.
13. **Aggiunto mutation testing** come quality gate alternativo a line coverage (line coverage = misura sintetica facile gameable).
14. **Aggiunto anti-monolith CI**: script CMake check 500-line limit per file, fail build se violato.
15. **Aggiunta error diagnostic framework**: ogni `Unimplemented` deve riportare: input shape, modulo, motivo dominio escluso, suggerimento utente, ticket ID. No silenzio.
16. **Aggiunto incremental build cache strategy** (ccache + Ninja) per ridurre rebuild time durante refactor pervasivo.
17. **Aggiunto pin versioni**: Maxima 5.49.0 (Homebrew bottle, GPL-2.0-only, sorgenti immutabili — vedi CLAUDE.md Regola 6), MPFR 4.2.1, GMP 6.3.0 (se usata in tests), C++20 con compilers documentati (Clang 17+, GCC 13+).

---

---

## FASE 0 — Sanitizzazione tracker (PREREQUISITO)

Prima di qualsiasi nuova implementazione, sanare i debiti documentati come "Risolta" ma non chiusi.

### F0.1 — Retroclassificazione CAS_TASKS.md
Retro-marcare a `Parziale` o `Parziale avanzata` con evidenza file:riga le voci sovrastimate. Lista minima:
- L1-01 Gruntz MRV (rank statico, comparazione fragile)
- L1-02 Risch (log-extension generale assente, trial constants hardcoded)
- L1-08 GCD multivariato (Brown stub)
- L1-17 Pivot Bareiss (magic 1000/500/400 vivi in `matrix_ops.cpp:243`)
- L2-04 Smith (solo Z, non PID)
- L2-06 fsolve (tolerance hardcoded, complessi spuri)
- L2-08 Polar/log (ln complesso parziale)
- L2-19 Equivalenza Risch (positivity inference debole)
- L2-22 Residue theorem (solo quadratici+biquadratici)
- L3-04 Funzioni speciali (Bessel/Chebyshev incompleti)
- L3-06 Trager tower (solo 2 livelli)
- L3-18 Galois (deg ≤4)

### F0.2 — HARDCODE_LEDGER.md sanato
Aggiungere voci mancanti:
- `polynomial_hensel.cpp:165` linear-lifting "for simplicity"
- `polynomial_groebner_f4_solver.cpp:17` Lex hardcoded
- `polynomial_gcd_multivariate.cpp:536` magic `*16`
- `integrate.cpp:40-67` double + to_u64() in core simbolico (REGOLA 1 violation)
- `fsolve.cpp:88` kTolerance constexpr 1e-10
- `integrate_risch.cpp:623` trial constants {±1, ±1/2, ±2}
- `ode_solver_advanced.cpp:227` "C1"/"C2" literal
- `ode_solver_1st_order.cpp:37` "C1" literal
- `matrix_solve.cpp:192` "t"+i naming
- `matrix_eigenvalues.cpp:263` "_lambda_" literal
- `matrix_jordan.cpp:123` "_lambda_" literal
- `evaluator.cpp:142-171` RootOf seed scheme deterministic
- `simplify_trig.cpp:279-282` kBaseAngles set chiuso 15 angoli
- `simplify_special_fn.cpp:111` bit_length>16 bail-out

### F0.3 — Test coverage baseline
- Installare `gcov`/`lcov` in CMake (`-fprofile-arcs -ftest-coverage`)
- Eseguire suite completa, generare `coverage/index.html`
- Documentare baseline: copertura per modulo. Target finale: ≥90% line coverage per area Risolta.

### F0.4 — Property-based test framework
- Aggiungere `rapidcheck` C++ a `test/`
- Property minime obbligatorie da subito:
  - `gcd(a,b)·lcm(a,b) ≡ a·b` su Z[x] random deg ≤30
  - `D(integrate(f)) ≡ f` su corpus 200 integrali (verifica via `simplify()`)
  - `f(roots(f)) ≡ 0` per ogni radice numerica/RootOf
  - `A·inv(A) ≡ I` su matrici random rationals 5×5
  - `factor(p) → ∏ fattori ≡ p` su Z[x] random
  - `partial_fractions(f) → riassemblata ≡ f`

### F0.5 — Golden test vs reference
- Corpus **2000 input** matematici (~200 per area)
- Reference primaria: **Maxima 5.49.0** (pin versione, Homebrew bottle GPL-2.0-only, sorgenti **immutabili** per integrità oracolo — vedi CLAUDE.md Regola 6) tramite `maxima --very-quiet --batch-string`. Hash binario verificato via `scripts/verify_maxima_integrity.sh`
- Reference secondaria: SymPy 1.13 per cross-check edge cases
- Confronto AST equivalenza (no toString)
- Pass-rate dashboard per modulo + delta vs commit precedente
- Bronstein "Symbolic Integration I" book esercizi numerati come sub-corpus calculus

### F0.6 — Build infrastructure CI
- Sanitizer (ASan + UBSan + TSan) integrati in CMake target `cas_tests_san`
- Benchmark gate: `cas_bench` con baseline pinned in `benchmarks/baseline_release.json`. Regression ≥10% blocca merge
- Anti-monolith CI: script `scripts/check_file_size.sh` fail se qualsiasi `.cpp` o `.hpp` >500 righe (eccezione documentata in `ANTI_MONOLITHIC_REPORT.md`)
- Diff-coverage CI: `diff-cover` o `codecov` con gate ≥85% sui nuovi line per PR
- Build cache: ccache + Ninja default per ridurre rebuild dopo refactor pervasivo
- Mutation testing: setup `mutest` o `mull` per quality gate alternativo (line coverage non sufficiente da solo)

### F0.7 — Doc-per-algorithm template
- Creare `.APROJECT_REFERENCES/algorithms/_TEMPLATE.md` con sezioni: Algoritmo (nome canonico), Riferimento (citazione bibliografica con pagina), Dominio input/output, Complessità tempo/spazio, Limiti dichiarati (cosa NON copre), Esempi accettati, Esempi rifiutati (con diagnostico), Test certificatore, Performance baseline
- Ogni algoritmo che entra in `Risolta` deve avere doc compilato in `.APROJECT_REFERENCES/algorithms/<nome>.md`
- PR senza doc → blocco merge

### F0.8 — Error diagnostic framework
- Estendere `CASErrorKind::Unimplemented` con payload strutturato: `{module, function, input_shape_descriptor, reason_code, suggestion, ticket_id}`
- Helper `make_unimplemented(module, fn, input, reason, suggestion, ticket)` obbligatorio in tutto src/
- Audit grep `return Unimplemented(` per assicurare payload strutturato ovunque
- Output utente: messaggio multi-linea con causa, suggerimento, link al ticket

**Exit gate F0**:
- Tracker sanato (12 voci retroclassificate)
- Ledger completo (15 voci aggiunte)
- Coverage baseline misurata + dashboard pubblicata
- Property tests verdi su ≥50 property
- Golden suite 2000 input eseguibile, pass-rate baseline misurato per modulo
- CI sanitizer + bench + anti-monolith + diff-coverage attivi
- Template doc disponibile
- Diagnostic framework rolled out

---

## FASE 1 — STRATO L0 FOUNDATION (chiusura 100%)

Senza foundation solido tutto sopra crolla. Nessun engine sopra L0 si tocca finché L0 chiuso.

### F1.1 — BigInt production-grade
**Algoritmo**: limbs `uint32_t` array (già OK). Riferimento: GMP `mpn_*`, Knuth TAOCP Vol 2, Brent-Zimmermann "Modern Computer Arithmetic".
- **Multiply**: schoolbook (n<32), Karatsuba (n<2048), Toom-3 (n<4096). Per n≥4096 fallback Karatsuba — Schönhage-Strassen FFT è **APERTA PERMANENTE** (HPP-F1.1-MUL in HARDCODE_LEDGER.md). Citare Granlund-Möller-Möller per Toom-3.
- **Divide**: Knuth Algorithm D (long division per limb). Burnikel-Ziegler divide-and-conquer per n>1000 — APERTA PERMANENTE (HPP-023 in HARDCODE_LEDGER.md). Knuth D corretto per tutti n; BZ è ottimizzazione O(M(n)log n) vs O(n²) per n grande.
- **GCD**: binary GCD (Stein), Lehmer GCD per integers large. No Euclidean naive.
- **Modexp**: Montgomery reduction o Barrett.
- **Primality**: Miller-Rabin deterministico 12 basi {2,3,5,7,11,13,17,19,23,29,31,37} per n<3.18×10²³ > 2^64 (Sorenson-Webster 2015, Math. Comp. 84, Table 2). Probabilistico (k=40 basi) per n>2^64. Nota: la precedente dichiarazione "BPSW" era imprecisa — l'implementazione usa MR a basi fisse, non BPSW (che richiederebbe MR base-2 + strong Lucas-Selfridge).
- **Integer factorization**: trial small primes + Pollard rho + Pollard p-1 + ECM (Lenstra) + Quadratic Sieve (per n>10^25).

**Certificatori**:
- `gcd(a,b)·lcm(a,b) ≡ a·b`
- `(a·b)/b = a` se b ≠ 0
- `modexp(a,b,m) ≡ a^b mod m` validato con piccoli b
- `is_prime(n) ⇒ Miller-Rabin pass su 40 basi`
- `factor(n) → ∏ p_i^{e_i} ≡ n`

**Anti-hardcode tests**: mul/div su n ∈ {32, 128, 512, 2048, 8192, 32768} limbs random. Fuzz 1000 input per ciascuno.

**Exit (ricalibrato, post-audit 2026-06-11)**: mul ≥10× più veloce sul caso 2048 limbs rispetto a baseline pre-Toom-3. Performance vs GMP su 1024-4096 limbs: target onesto **30-50%** (non 80%). **Schönhage-Strassen NON è exit-gate F1**: è `Aperta permanente` (HPP-F1.1-MUL); pipeline F1 chiude su Karatsuba + Toom-3 + fallback Karatsuba per n≥4096. FFT BigInt resta `Aperta` permanente come research-grade.

### F1.2 — Rational
- Normalizzazione canonica obbligatoria post-ogni-op (già OK)
- `to_continued_fraction(n_max)` per L3 numerica
- Optional: dyadic rationals fast-path se denom = 2^k

### F1.3 — Arena + hashconsing
- Mantenere bump allocator esistente (OK)
- Rimuovere `mutex_` o sostituire con shard locks per build multi-thread
- `make_fresh_symbol(prefix)` deve essere funzione pubblica garantita unica usata ovunque

### F1.4 — Simplifier core canonico
- `simplify_arithmetic`: completare like-term collection su coeff simbolici (`a·x + b·x → (a+b)·x` con a, b simbolici)
- `simplify_trig`: sostituire `kBaseAngles` set chiuso con generatore via Chebyshev minimal polynomial. Per `cos(p·π/q)` qualsiasi q ≤ `ctx.max_trig_denom()`:
  - Costruire polinomio minimo di `cos(π/q)` via Chebyshev `T_q(x) = -1` o ciclotomico Φ_{2q}
  - Per q non costruibile (q=7, q=11, ...) emettere `RootOf` con polinomio minimo
  - Per q costruibile (Fermat primes 17, 257, 65537 + 2^k · prodotti) emettere espressione in radicali via Gauss period
- `simplify_exp_log`: denesting Landau generale (`sqrt(p + q·sqrt(c) + r·sqrt(d) + s·sqrt(cd))`) + Cardano cube root denesting
- `normal_form` trascendente assumption-aware:
  - `ln(x·y) → ln(x)+ln(y)` solo se `x>0 ∧ y>0` (gate strict_branch_cuts)
  - `exp(ln(x)) → x` solo se `x>0`
  - Direzione contraria (contraction): `ln(x)+ln(y) → ln(x·y)` su positivi
- `assumptions` inferenza algebrica:
  - `x>0 ∧ y>0 ⊢ x+y>0` (sum positive)
  - `x>0 ∧ y>0 ⊢ x·y>0` (prod positive)
  - `x≠0 ⊢ x²>0` (square nonzero positive)
  - `x>a ∧ a>b ⊢ x>b` (transitivity già presente — verificare 3-hop)
  - `x∈[a,b] ∧ y∈[c,d] ⊢ x+y∈[a+c,b+d]` (interval arith)
  - `assume(x²>0) ⊢ assume(x≠0)` (square positive implies nonzero)

**Certificatori**:
- Property `simplify(simplify(e)) ≡ simplify(e)` (idempotenza)
- Property `simplify(e1) = simplify(e2) ⇔ e1 ≡ e2` su corpus equivalenze note
- Cycle detection in rewrite engine (già presente, verificare invariant)

### F1.5 — Pattern matcher e rewrite engine
- Discrimination net per indexing regole per fast match O(log n_rules)
- Pattern types tipati: `x_Integer`, `x_Positive`, `x_Symbol`, repeated `x__`, `x___`
- AC-matching con canonical sort pre-match per evitare backtracking esponenziale

### F1.6 — Complex foundation (trasversale L0)
I numeri complessi sottostanno polynomi su C, eigenvalues, branch cuts, residue. Devono essere strato L0, non L3.
- `GaussianInt` Z[i] arithmetic completo (presente parziale `gaussian_int.cpp` 53 LOC) → estendere a factor Z[i] via Gaussian primes
- `ComplexLit(re, im)` come AST first-class con re, im Rational (non DecimalLit)
- Arithmetic complex su Q[i]: add/sub/mul/div esatti
- Conversion `RootOf(min_poly) ↔ AlgebraicComplex` quando radici complesse
- `i` come costante simbolica, no Symbol("i") generic
- `abs(z) = sqrt(re²+im²)` con denesting per RootOf
- `arg(z)` su branch principale (-π, π] (presente parziale)
- Branch cut policy globale via `ctx.complex_branch_policy` enum: PrincipalCut (default), MultiSheet (research)

**Certificatori**:
- `(a+bi)·(a-bi) ≡ a²+b²` (sempre, no assumptions)
- `factor_gaussian(5) ≡ (2+i)(2-i)` (Gauss prime decomposition)
- `factor_gaussian(p) per p prime con p ≡ 1 mod 4 → 2 fattori non triviali`
- `abs(re+im·i)² ≡ re²+im²` su Q[i]

**Exit gate F1**:
- BigInt benchmark ≥ 30-50% performance GMP su input random fino 10^4 cifre (ricalibrato onesto, pipeline Karatsuba+Toom-3 — **SS escluso da F1 exit**, Aperta perm HPP-F1.1-MUL)
- Tutti property test foundation passano
- `assume()` inferenza algebrica copre 50 implicazioni canoniche
- Hashconsing thread-safe sotto 4 thread benchmark senza contesa visibile
- Complex Q[i] arithmetic completo + factor Gaussian primes funzionante
- Coverage `src/foundation/` + `src/symbolic/` ≥ 92% (line) + mutation score ≥70%
- 0 hardcode non-ledger nei file toccati (verifica via grep + CI gate)
- Public API congelata e documentata in `include/cas/*.hpp` con doxygen

---

## FASE 2 — STRATO L1 POLINOMI UNIVARIATI COMPLETI

### F2.1 — GCD univariato production
- Subresultant PRS (Brown-Collins) OK presente
- **Half-GCD (Knuth-Schönhage) INTEGRATO** (F2 Block A, R2+R4 remediation 2026-05-27):
  - `polynomial_half_gcd.cpp`: bounds derivati da invariant matematici (non magic constants).
  - `gcd_integer_poly_dispatch` ora usato da tutti i call-site (factorization_integers.cpp,
    factorization_poly_gcd.cpp, polynomial_square_free.cpp, builtin_rewrite.cpp).
  - Certificatori CERT1-CERT5 (deg 200-400, sparse) PASS.
- **Modular GCD CRT multi-prime** con Chinese remainder accumulato fino a bound di Mignotte

### F2.2 — Resultant production
- Subresultant PRS (presente)
- **Modular resultant CRT** per coefficienti grandi
- **Bivariate resultant** ottimizzato (gate per csolve 2-var)

### F2.3 — Factor Z[x] production
- Squarefree Yun (corretto)
- Mod-p factor via Cantor-Zassenhaus (presente)
- Hensel lifting **quadratico** (non solo lineare) — eliminare `polynomial_hensel.cpp:165` "for simplicity"
- **vanHoeij combinatorial knapsack** sui digit-vectors dei fattori modulo p (sostituisce recombination subset O(2^n))
- Lemma di Hadamard per bound coefficienti
- Test stress: fattorizzare polinomio random Z[x] deg 100 con 30+ fattori mod p

### F2.4 — Factor Fp[x] production
- **Berlekamp INTEGRATO** (F2 Block A, R4 remediation 2026-05-27):
  - `factorization_berlekamp.cpp` ora chiamato dal dispatcher in `factorization_integers.cpp:329`
    quando `deg(f)*p ≤ ctx.max_berlekamp_matrix_size()` (default 1024).
  - Cantor-Zassenhaus usato automaticamente come fallback per input più grandi.
- EDF p=2 (presente)
- Distinct degree factorization (presente)

### F2.5 — Partial fractions production
- Bronstein LRT con RootSum (presente)
- Verificare RootSum su radici irrazionali grado ≥3

### F2.6 — Cyclotomic detection
- Möbius construction (presente, OK)
- **Silent-truncation bug CHIUSO** (F2 Block A, R1 remediation 2026-05-27):
  - `is_cyclotomic()` ora usa bound matematicamente provato: `max_n = max(6, 2*d²)`
    (da φ(n) ≥ √(n/2) → n ≤ 2φ(n)² per n > 6, Rosser-Schoenfeld).
  - Fix copre tutti n con φ(n) = deg(poly): es. φ(18)=6 era mancato con la vecchia formula.
  - Per deg > 724 (n > 2^20), `is_cyclotomic` segnala esplicitamente nullopt (A5-LARGECYCLO),
    non restituisce risultato ambiguo vuoto.
- Rimuovere bound max_cyclotomic_n: generatore on-demand per qualsiasi n (A5-LARGECYCLO)

**Certificatori**:
- `factor(p) → ∏ fattori ≡ p` con `assert(simplify(prod - p) == 0)`
- `gcd(a,b)` divide entrambi: `divrem(a, gcd) == (q, 0)`
- `resultant(p, p') == 0` per p con radice multipla
- `partial_fractions(f) → riassemblata ≡ f`

**Anti-hardcode tests minimi per ogni algoritmo**:
- 3 nominali
- 3 equivalenti AST diversi
- 2 variabili rinominate
- 2 coefficienti razionali grandi (>10^6)
- 2 gradi alti (>50)
- 2 degeneri (deg 0, 1)
- 2 fuori dominio (segnale Unimplemented esplicito)
- 1 property metamorphic
- 1 regressione anti-hardcode storico

**Exit gate F2**:
- Fattorizzare 100 polinomi random Z[x] deg ≤100 in <30s totali
- Brown GCD multivariato spunta da `polynomial_gcd_modular.cpp:44` (preparazione per F3)
- 0 Unimplemented in algebra univariata su Z/Q
- Coverage `src/algebra/polynomial_*.cpp` ≥ 90%

---

## FASE 3 — STRATO L2 POLINOMI MULTIVARIATI + ESTENSIONI ALGEBRICHE

### F3.1 — GCD multivariato production (priorità #1)
- **Brown's modular GCD** (sostituisce stub `polynomial_gcd_modular.cpp:44`)
- **Zippel's sparse interpolation** per polinomi sparsi multivariati (IMPLEMENTATO 2026-06-02 S1/A1 in `src/algebra/polynomial_sparse_interpolation.cpp`; n-variate generale con distinct-prime guarantee + retry-on-singular configurable via `ctx.sparse_interp_max_retries()`)
- **EZ-GCD** (Wang) per cofactors
- Bound Mignotte multivariato

### F3.2 — Factor multivariato Wang
- Wang multivariate factorization (riferimento: Wang 1978)
- Hensel lifting multivariato (presupposto F3.1)
- Leading coefficient determination via evaluation
- Multivariate squarefree

### F3.3 — Groebner basis production
- F4 (presente) + sugar + GM (presente)
- Aggiungere **F5** (Faugère) per signature-based redundancy elimination
- **FGLM**: conversione tra ordini monomiali (GRevLex → Lex per solve)
- Rimuovere `polynomial_groebner_f4_solver.cpp:17` Lex hardcoded → usare FGLM
- Reduced Groebner basis garantita

### F3.4 — Algebraic numbers tower production
- Q(α) arith (presente, OK)
- Q(α_1, ..., α_n): **primitive element theorem** (Trager): trovare θ = ∑ s_i·α_i con s_i piccoli interi tale che Q(θ) = Q(α_1, ..., α_n)
- Risolvere minimo polinomio di θ via resultant
- Espressione di ogni α_i come polinomio in θ
- Esportare via `algebraic_tower_bridge.cpp` come API stabile

### F3.5 — Factor Q(α) production
- Trager su Q(α) (presente per 2-level)
- Estensione a tower ≥3 livelli (gate F3.4)

### F3.6 — Galois group ≥5
- Discriminante computation (presente)
- Soubin-Stauduhar deg-5 (S5, A5, D5, F20, C5)
- Algoritmi per deg 6,7,8 (riferimento: Cohen "A course in computational algebraic number theory" cap. 6)

**Certificatori**:
- `factor_multivariate(p) → ∏ ≡ p` (assert simplify)
- Groebner: `f mod G ≡ 0 ⇔ f ∈ ideal(G)` (verifica via membership query)
- Primitive element: `Q(θ) ⊇ Q(α_i)` verificato (α_i polinomio in θ)
- Galois group cardinality divide n!

**Exit gate F3**:
- GCD multivariato 3 var deg 20 cofattori coprimi: <5s
- factor su Q(√2, √3, √5) attivato e testato
- F4 + FGLM su sistema 3×3 random deg 4: <60s
- 0 Unimplemented "primitive element" nel tower bridge

---

## FASE 4 — STRATO L2 LINEAR ALGEBRA COMPLETA

### F4.1 — Decomposizioni production
- LU con **partial pivoting numerico** vero (oggi: first non-zero)
- LU con **PivotScore** (oggi: header mente, codice no)
- **QR via Modified Gram-Schmidt** (Trefethen-Bau §8) — sostituisce GS classico instabile. **Householder QR è APERTA PERMANENTE** (HPP-F4.1-QR-HOUSEHOLDER): AST explosion su Q simbolico per matrici ≥8×8 (vedi HC-F4-QR-SYMBOLIC-TIMEOUT chiuso 2026-05-XX ledger). MGS mantiene rational updates → no timeout. Householder require Trefethen-Bau §10 con riflettori `I - 2vv^T/v^Tv` che generano coefficienti razionali esplosivi in simbolico esatto.
- **Cholesky LDL^T** per matrici simmetriche
- Implementare Bareiss in UN file `matrix_determinant.cpp`; eliminare duplicato in `matrix_ops.cpp:243-249`

### F4.2 — Forme canoniche production
- Jordan canonical form: routing autovalori RootOf → `null_space_over_extension` (presente). Oggi `matrix_jordan.cpp` non lo chiama
- Smith Normal Form generalizzata: **Q[x] PID** (oggi solo Z). Riferimento: Storjohann LLL-based Smith
- **Hermite Normal Form** per Z e Q[x]
- **Companion matrix** API pubblica

### F4.3 — Determinanti speciali
- Tridiagonale (presente)
- **Vandermonde** formula chiusa ∏(x_j-x_i)
- **Toeplitz** via Levinson
- **Circulant** via DFT
- **Banda** generale

### F4.4 — Eigen completo
- Eigenvalues delega solve_polynomial (OK)
- Eigenvectors: rimuovere kernel error swallowing (`matrix_eigenvalues.cpp:324`)
- Multipli autovalori: catene di Jordan affidabili sotto RootOf

### F4.5 — Fresh-symbol pervasivo
- Sostituire `"t"+i` in `matrix_solve.cpp:192` con `ctx.make_fresh_symbol("t")`
- Sostituire `"_lambda_"` in `matrix_eigenvalues.cpp:263` e `matrix_jordan.cpp:123`

**Certificatori**:
- `A·inv(A) ≡ I` per A invertibile random 8×8 rationals
- `LU(A) → L·U ≡ P·A`
- `QR(A) → Q·R ≡ A ∧ Q^T·Q ≡ I` (verifica su MGS; Householder = Aperta perm)
- `det(A) ≡ ∏ eigenvalues(A)` (counting molteplicità)
- Jordan: `A ≡ P·J·P^{-1}` con J Jordan canonical
- Smith: `U·A·V ≡ D` con D diagonale di invariant factors

**Exit gate F4**:
- Matrici 10×10 random Q: LU/QR/inv/det <100ms
- Jordan su companion matrix random deg 5: <2s
- Smith su Z^{6×6} <500ms
- 0 hardcode literal `"C1"/"_lambda_"/"t1"`

---

## FASE 5 — STRATO L2 CALCULUS COMPLETO

### F5.1 — Differential field + Risch production
- Bronstein "Symbolic Integration I" book completo
  - Capitolo 5: rational integration (presente, Hermite)
  - Capitolo 6: integration of transcendental functions (parziale)
  - Capitolo 7: parametric problems (assente)
  - Capitolo 8: Risch differential equation generale (parziale)
  - Capitolo 9: structure theorem completo (assente — critico)
- Sostituire trial constants `{±1, ±1/2, ±2}` in `integrate_risch.cpp:623` con risoluzione formale del coefficiente costante via residue field equation
- Log-extension polynomial integration generale (oggi: `:815` Unimplemented)
- Liouville theorem: decidere se F è elementare esprimibile

### F5.2 — Limit MRV/Gruntz production
- Algoritmo Gruntz §3.5 completo: comparazione asintotica via Mrv set ricorsivo + serie Taylor leading
- Sostituire growth rank statico con confronto dinamico
- Cancellation tower con livelli arbitrari
- Test corpus: limiti da Gruntz tesi 1996 (100+ casi)

### F5.3 — ODE production
- ODE classifier emettere TUTTI i tipi: separable, linear, exact, Bernoulli, homogeneous, Riccati, Clairaut, d'Alembert
- ODE 1st order: solver per ogni tipo classificato (oggi solo Linear1st reale)
- ODE 2nd order: const-coeff (presente) + Euler equation + Cauchy-Euler + variazione parametri (presente parziale)
- ODE Frobenius: caso radici differing-by-integer con log term (oggi Unimplemented `:201-205`)
- ODE Frobenius resonance generale
- ODE Laplace transform method (presente)
- Sostituire literal `"C1"/"C2"` con `ctx.make_fresh_symbol("C")` ovunque

### F5.4 — Integrate definite + improper
- Hadamard finite part per poli ordine ≥2 generale (oggi parziale)
- Improper su singolarità trascendenti
- Riemann-Liouville fractional integral
- Sostituire `integrate.cpp:40-67` `double + to_u64()` con bound BigInt esatto

### F5.5 — Series production
- Taylor generale via diff iterata (presente, OK costoso)
- Laurent: pole order arbitrary (oggi limite 4)
- **Puiseux series** per branch point algebrici
- Padé: rimuovere vincolo coeff Q (oggi fallisce su Taylor con π, e, sqrt)

### F5.6 — Residue theorem completo
- Espandere `residue_theorem.cpp` oltre quadratici+biquadratici
- Polinomi qualsiasi grado al denominatore via fattorizzazione + residui ai poli
- Contour integration generale (Jordan lemma per integrandi exp(iax)/Q(x))
- Mellin transform via contour

### F5.7 — Summation simbolica
- **Gosper's algorithm** per hypergeometric closed form (IMPLEMENTATO 2026-06-02 S2/A2 in `src/symbolic/summation_gosper.cpp`; chiusura Petkovšek pulita su term=1, term=k, term=1/(k(k+1)), term=1/(k²+1)→nullopt. Closes HC-F4-GOSPER-CONSTANT-HANG)
- **Zeilberger's algorithm** (creative telescoping)
- **WZ pair method** (Petkovsek-Wilf-Zeilberger)
- **Abramov's algorithm** per rational function summation
- Bernoulli numbers (presente)

### F5.8 — Transforms
- Laplace transform: estendere tabella + algoritmo Bronstein per inverse Laplace via residue
- Fourier transform: aggiungere
- Mellin transform: aggiungere
- Z-transform: aggiungere

### F5.9 — Funzioni speciali
- Gamma, Bessel, Erf, Zeta (parziali)
- Hypergeometric `_pF_q(a;b;z)` come prima-class builtin
- Contiguous relations
- Casi closed-form noti (Gauss 2F1(a,b;c;1), Saalschütz)
- Jacobi P_n^{(α,β)}
- Elliptic integrals K, E, Π, F
- Rimuovere bail-out `bit_length>16` in `simplify_special_fn.cpp:111`

**Certificatori**:
- `D(integrate(f)) ≡ f` (verifica via simplify) — già presente, estendere
- `taylor(f, x, x0, n) → polinomio` con `D(F)(x0) = (1/k!)·D^k(f)(x0)`
- `laplace(L^{-1}(F)) ≡ F`
- ODE: sostituire soluzione → 0
- Residue theorem: confronto contro numerical contour integration
- Summation closed-form: somma finita coincide su 10 valori N

**Exit gate F5**:
- Bronstein book 50 integrali test corpus: 90% PASS
- Gruntz tesi 100 limiti: 95% PASS
- ODE Bernoulli/Exact/Riccati riconosciute e risolte
- Residue theorem su deg 5,6,7,8 polinomi PASS
- Zeilberger su corpus 20 identità ipergeometriche PASS

---

## FASE 6 — STRATO L3 NUMERICA + COMPLEX + UNITS

### F6.1 — MPFR integrazione completa
- BigFloat (presente, OK)
- bigfloat_eval unificare env con BigFloat (oggi env `map<string,double>` perde precisione)
- Adaptive precision: aumentare prec se cancellazione catastrofica detected
- Complex MPFR (MPC library)

### F6.2 — Complex completo
- AlgebraicNumber tower (gate F3.4)
- ln(z) generale `ln|z| + i·arg(z)` su branch principale (presente parziale L2-17)
- Branch cuts globali con `strict_branch_cuts` enforced
- Multi-valued representation come Riemann sheet (research target)
- Z[i] Gaussian factorization (oggi solo gcd)
- Estensioni Q(α) complesse

### F6.3 — Equazioni trascendenti
- fsolve: rimuovere `kTolerance = 1e-10` hardcoded → ctx
- Filter reale vs complesso garantito (oggi: complessi spuri tra reali)
- Newton-Raphson per multi-equazione sistema trascendente
- Continuation methods per sistema non lineare

### F6.4 — Disequazioni 1-var via Sturm + risultanti (NO CAD in F6)
CAD generale spostato a Fase 8 (post-parità). Complessità doppia esponenziale + scope corretto incompatibile con strato L3.
- `solve_inequality(poly, var, ctx)` per polinomi univariati via Sturm sequence: isolate radici reali, segno tra radici, intervalli soluzione
- `solve_inequality(poly_system, vars, ctx)` multi-var deferito a Fase 8
- Disequazioni con trascendenti (es. `sin(x) > 1/2`): periodicità + Sturm su un periodo
- Fase 8 (post-parità): CAD McCallum completo come strato L4 research

### F6.5 — RootOf evaluator
- Fix seed scheme bug in `evaluator.cpp:142-171`: usare Sturm isolation interval come bracketing per Newton, garantire unicità per root_index

### F6.6 — Unit system (SI completo)
- Quantità con unità: `1.5*meter`
- Conversioni automatiche
- Dimensional analysis check su +,-,=
- Base SI + derived + prefissi (k, M, m, μ, n, ...)
- Costanti fisiche (c, h, k_B, ...)

**Certificatori**:
- MPFR: numeric eval di π@1000 digits confronta cifra-per-cifra con MPFR diretto
- Branch cuts: `ln(-1) ≡ i·π`, `sqrt(-1) ≡ i`
- CAD: regione output verificata su 100 punti campione
- Units: `(1*meter)+(1*second)` errore dimensionale

**Exit gate F6**:
- MPFR 10000 digit eval di expressions complesse <1s
- fsolve restituisce solo reali su input reali (filtraggio garantito)
- Disequazioni 1-var Sturm-based: 50 test PASS
- Unit system 50 conversioni test PASS
- CAD multivar marker permanente `Aperta` (Fase 8 post-parità)

---

## FASE 7 — INTEGRAZIONE FINALE + PLOTTING + PARITÀ HP PRIME G2

### F7.1 — Plotting 2D/3D
- Adaptive sampler (presente parziale)
- Implicit plot
- Parametric plot
- Contour plot
- Vector field plot

### F7.2 — Statistics package
- Distribuzioni: Normal, Binomial, Poisson, χ², t, F
- pdf, cdf, quantile, sampling
- Hypothesis testing
- Linear regression

### F7.3 — Numerical methods package
- ODE numerica: RK4/RKF45 (presente, fix error-gate bug)
- Adaptive quadrature: Gauss-Kronrod (oggi solo Simpson)
- Interpolation: spline, Hermite, Lagrange

### F7.4 — Parità HP Prime G2 final acceptance
- Eseguire **TUTTO** il test corpus HP Prime G2 reference (vedere `AUDIT_CAS_vs_HP_Prime_2026-05-04.md`)
- Pass-rate target: ≥95%
- Performance: ≥80% di Giac-Xcas su benchmark identici
- Documentation: API reference completo

**Exit gate F7 = parità HP Prime G2 dichiarata**.

---

## FASE 7.5 — CLOSURE GAP PASS-RATE (post-C1, target aggregato 85-88%)

> Inserita 2026-06-09 post-misurazione C1 pilot (vedi
> `F7_GOLDEN_CORPUS_REPORT_2026-06-09.md` + `F7_PARITY_SUMMARY_2026-06-09.md`).
> Aggregato corrente: **77.6%** (512 PASS / 660 non-skip su 8 aree).
> Obiettivo F7.5 onesto: **86-88%**. Non 95%: la chiusura di Risch
> structure theorem full + hypergeometric `_pF_q` completo + Galois ≥5
> resta `Aperta` permanente come dichiarato in QA REVIEW v2 §10.
>
> **Regola operativa F7.5**: ogni task qui sotto è bound da
> `CLAUDE.md` REGOLA ZERO (niente shortcut), REGOLA 0.1 (spec read
> first), REGOLA 0.2 (zero test disable), 500 LOC anti-monolito, BigInt
> only, no `throw/catch`, `Result<T>` ovunque. Ogni hardcode-of-passage
> → `HARDCODE_LEDGER.md`. Pass-rate goals per area sono **vincolanti**
> per chiusura: se la misura non raggiunge il floor, il task resta
> in_progress (no `completed` ottimistico).

### Pass-rate target per area (binding)

| Area | Oggi | Target F7.5 | Leva principale |
|---|---:|---:|---|
| factor       |  99.0% |  99.0% | freeze (no regression) |
| gcd          | 100.0% | 100.0% | freeze |
| simplify     |  92.0% |  95.0% | identità trig/log mancanti |
| diff         |  82.5% |  92.0% | sech rewriter + asin/acos chain |
| limit        |  82.1% |  88.0% | Gruntz §3.5 + nested log |
| special_fn   |  64.9% |  82.0% | identità Γ/erf/Bessel canonicalization |
| integrate    |  43.6% |  75.0% | Risch log-ext + transcendental case |
| series       |  37.1% |  75.0% | Taylor canonical + ordine-aware diff |
| solve        |  (gap) |  92.0% | adapter test + set-equality |
| matrix       |  (gap) |  92.0% | adapter test |
| bronstein    | (hang) |  70.0% | runner robustness + Risch overlap |

**Aggregato pesato atteso**: ~86-88% (512+X / 1026, X = ~340 nuovi PASS).

---

### F7.5.A — Test infrastructure completion

Lavoro su `test/golden/` e `corpus_runner`. Nessuna modifica al core
simbolico. Effort: ~3 giorni T1-Sonnet.

#### F7.5.A1 — Solve adapter (Maxima list → set ExprPtr)
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Solve_Adapter.md` (DA SCRIVERE PRIMA — REGOLA 0.1).
- **File**: `test/golden/maxima_parser.hpp` (estendere), nuovo
  `test/golden/solve_set_equal.hpp`.
- **Algoritmo**: parse `[x = r₁, x = r₂, …]` → `std::vector<ExprPtr>`
  RHS; confronto via `mathematically_equal` su ogni coppia con
  bijection check (Hopcroft-Karp su grafo bipartito di compatibilità).
  Nessuna heuristic on cardinality matching.
- **Acceptance**: corpus area `solve` (81 entry) → ≥ 90% PASS.
- **No-shortcut**: vietato confrontare solo prima soluzione, ignorare
  soluzioni complesse, troncare set ≥ N.

#### F7.5.A2 — Matrix adapter (`[[…]]` → `MatrixLit` dispatch)
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Matrix_Adapter.md`.
- **File**: `test/golden/corpus_runner.hpp` (rimuovere skip lista
  `det/trace/transpose/rank/inv/eigenvalues`), nuovo
  `test/golden/matrix_parser.hpp`.
- **Algoritmo**: parser `[[a,b],[c,d]]` → `MatrixLit`; dispatch via
  `cas::linalg::determinant/trace/transpose/rank/inverse/eigenvalues`
  già esistenti.
- **Acceptance**: area `matrix` (79 entry) → ≥ 90% PASS. Eigenvalues
  con caratteristico Bareiss su Z[x] esatto.

#### F7.5.A3 — Runner robustness (per-entry timeout + output cap)
- **File**: `test/golden/main.cpp`.
- **Algoritmo**: `setitimer(ITIMER_REAL, …)` per-entry con SIGALRM
  → entry SKIP con motivo "TIMEOUT N s"; cap output `format_expr` a
  4 KB con sentinel `<truncated N bytes>`; CTest timeout 1800 s
  globale.
- **Acceptance**: bronstein corpus (90 entry) completato end-to-end;
  integrate corpus 140/140 (non più 116/140).

#### F7.5.A4 — Sech/Csch normalization in `algebraic_equal`
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Sech_Csch_Identity.md`.
- **File**: `src/algebra/algebraic_equal.cpp` (helper privato),
  nuovo unit test `test/unit/algebra/test_sech_csch_normal.cpp`.
- **Algoritmo**: traversal pre-confronto che riscrive `FuncCall("sech",
  x) → 1/cosh(x)`, `FuncCall("csch", x) → 1/sinh(x)`, idem `coth →
  cosh/sinh`, `tanh → sinh/cosh`. Applicato simmetricamente su lhs +
  rhs prima di `polynomial_normal_form`.
- **Acceptance**: corpus diff 80/80 entry produce 0 FAIL su mismatch
  notazionale sech/csch/coth (oggi 2 FAIL). diff pass-rate → ≥ 90%.
- **No-shortcut**: non aggiungere `BuiltinOp::Sech/Csch` (toccherebbe
  76 switch enum — fuori scope F7.5; lavoro per Fase 8 con
  Extended-Real migration).

---

### F7.5.B — Risch completion (integrate gap)

Target integrate: 43.6% → 75%. Lavoro su Bronstein cap. 5-8 reali, no
pattern matching tabellare. Effort: ~3-4 settimane T3-Opus.

#### F7.5.B1 — `arctan/asin/acos` standalone integrali
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Inverse_Trig.md`.
- **File**: nuovo `src/calculus/integrate_inverse_trig.cpp`
  (split anti-monolito da `integrate_core.cpp` se necessario).
- **Algoritmo**: integration by parts come dispatcher generale (non
  caso speciale), `∫f·g' = f·g - ∫f'·g`, con scelta `f = trig_inv(x)`,
  `g' = poly(x)`. Coverage: `atan`, `asin`, `acos`, `atan2`, e
  composizioni `x^n·trig_inv(x)` per n ≥ 1.
- **Acceptance**: ≥ 30 entry corpus integrate che oggi falliscono
  `INTEGRATE_NO_STRATEGY` su trig_inv passano. Pass-rate integrate
  area → ≥ 60%.
- **No-shortcut**: vietato pattern `if (is_atan(arg)) return x·atan(x) -
  ½ln(1+x²)`; deve essere by-parts strutturale.

#### F7.5.B2 — Bronstein cap. 5 (Hermite reduction completa)
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Hermite_Cap5.md`.
- **File**: `src/calculus/integrate_hermite.cpp` (già splittato in A1.1).
- **Algoritmo**: completare Hermite reduction su denominatore
  square-full (oggi parziale); LRT (Lazard-Rioboo-Trager) con RootSum
  formale (presente, da estendere a campo `Q(α)` con α algebrico di
  grado ≥ 3).
- **Acceptance**: corpus `bronstein/integrals.jsonl` (90 entry) → ≥ 70%
  PASS (era 0% per hang). 90 entry sono Bronstein book esempi
  numerati, copertura misurabile entry-per-entry.

#### F7.5.B3 — Risch transcendental case (Bronstein cap. 8)
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Risch_Transcendental_Cap8.md`.
- **File**: `src/calculus/integrate_rde.cpp` (presente, parziale),
  `src/calculus/differential_field.cpp`.
- **Algoritmo**: Risch differential equation generale per estensioni
  `θ = log(u)` e `θ = exp(u)` con `u ∈ Q(x, …)`. Sostituire trial
  constants hardcoded `{±1, ±1/2, ±2}` (ledger esistente, non chiuso)
  con risoluzione formale residue field equation.
- **Acceptance**: ≥ 20 integrali corpus che combinano `exp + log + poly`
  passano. Pass-rate integrate → ≥ 75%.
- **No-shortcut**: rimozione trial constants è obbligatoria (è già nel
  ledger debiti).

---

### F7.5.C — Taylor canonical form (series gap)

Target series: 37.1% → 75%. Effort: ~2 settimane T2-Sonnet-thinking.

#### F7.5.C1 — Canonical truncated polynomial form
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Taylor_Canonical.md`.
- **File**: nuovo `src/calculus/series_canonical.cpp`,
  `src/algebra/algebraic_equal.cpp` estensione.
- **Algoritmo**: rappresentazione `(coefficients[0..n], order_remainder)`
  per Taylor troncato; confronto pairwise coefficient + verifica
  order_remainder compatibile (`O(x^min(n_lhs, n_rhs))`). Nessun
  confronto via differenza simplify (causa dei FAIL attuali).
- **Acceptance**: `mathematically_equal` su `(1 + x + x²/2 + O(x³),
  1 + x + x²/2 + x³/6 + O(x⁴))` ritorna `true` (oggi `false`).
  series pass-rate → ≥ 70%.

#### F7.5.C2 — Padé completion (coeff Q senza vincolo)
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Pade_Generic_Coeff.md`.
- **File**: `src/calculus/pade.cpp` (presente, vincolo Q).
- **Algoritmo**: estensione a `Q(π, e, sqrt(2), …)` via campo
  algebrico generato dinamicamente. Riusare LLL già presente per
  basis minimization.
- **Acceptance**: 5 esempi Padé corpus con coefficienti irrazionali
  passano.

---

### F7.5.D — Limit gap closure

Target limit: 82.1% → 88%. Effort: ~1 settimana T2-Sonnet-thinking.

#### F7.5.D1 — Gruntz §3.5 nested log tower
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Gruntz_Nested_Log.md`.
- **File**: `src/calculus/limit_mrv.cpp` (split A1.3 — moduli MRV).
- **Algoritmo**: estensione Mrv set per torri `log(log(log(x)))` e
  combinazioni `exp(log(x)^k)`. Algorithm Gruntz §3.5 reale, no
  pattern matching su forma.
- **Acceptance**: 8 entry corpus limit con torri nested log passano.
  limit pass-rate → ≥ 88%.

---

### F7.5.E — Special functions canonicalization

Target special_fn: 64.9% → 82%. Effort: ~1.5 settimane T2-Sonnet-thinking.

#### F7.5.E1 — Identità Γ/B/ζ/erf canonical
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Special_Fn_Identities.md`.
- **File**: nuovo `src/symbolic/simplify_special_fn.cpp`,
  estensione `src/algebra/algebraic_equal.cpp`.
- **Algoritmo**: applicazione identità note (`Γ(n+1) = n!`,
  `Γ(x)·Γ(1-x) = π/sin(πx)`, `B(x,y) = Γ(x)Γ(y)/Γ(x+y)`, `erf(0) = 0`,
  `erf(-x) = -erf(x)`, `ζ(2) = π²/6`) come rewrite rules orientate
  LPO. Nessun cataloghi chiusi.
- **Acceptance**: special_fn pass-rate → ≥ 82%.

#### F7.5.E2 — Bessel `J_n / Y_n / I_n / K_n` parser + simplify
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Bessel_Identities.md`.
- **File**: `include/cas/builtin_functions.hpp` (+ entries),
  `src/symbolic/simplify_special_fn.cpp`.
- **Algoritmo**: ricorrenze Bessel (`J_{n-1}(x) + J_{n+1}(x) =
  (2n/x)·J_n(x)`) + valori speciali (`J_0(0) = 1`, `J_n(0) = 0` per
  `n ≥ 1`).

---

### F7.5.F — Extended-Real AST migration (HC-F70-A43 close)

Effort: ~1 settimana T1-Sonnet (lavoro meccanico, no creatività).

#### F7.5.F1 — `MathConstant::NegInfinity / ComplexInfinity` propagation
- **Spec**: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md`.
- **File**: `include/cas/ast.hpp` (enum extension), tutti i 76 switch
  in `src/symbolic/*.cpp`, `src/algebra/*.cpp`, `src/calculus/*.cpp`,
  `src/numeric/*.cpp` (script di scan automatico tramite
  `scripts/find_constant_switch.sh` da scrivere).
- **Acceptance**: build pulita con `-Wswitch -Werror`; nuovi unit test
  in `test/unit/symbolic/test_extended_real.cpp` (≥ 15 casi:
  `lim x→∞ x = +∞`, `1/0⁺ = +∞`, `log(0⁺) = -∞`, `(-∞) + ∞ =
  indeterminate`, `0·∞ = indeterminate`, ecc.).
- **Ledger**: chiude HC-F70-A43-EXTENDED-REAL.
- **No-shortcut**: vietato `default: return Unimplemented` silenzioso;
  ogni switch deve gestire NegInf/ComplexInf esplicitamente o
  ritornare `Indeterminate` diagnostico.

---

### F7.5.G — Giac install + C2 benchmark

Effort: ~2 giorni T1-Sonnet.

#### F7.5.G1 — Build Giac da sorgente o Homebrew bottle
- **File**: `scripts/install_giac.sh`, `scripts/giac_integrity.sh`,
  `scripts/run_corpus_giac.sh` (analogo a `run_golden_maxima.sh`).
- **Algoritmo**: priorità (a) `brew install giac`, fallback (b) build
  sorgente con pin versione + sha256 manifest (analogo Maxima).
- **Acceptance**: `giac --version` printa, `scripts/giac_integrity.sh`
  verde, hash binario in `scripts/giac_X.Y.Z_manifest.sha256`.

#### F7.5.G2 — Benchmark vs Giac come secondo oracolo
- **File**: `test/golden/giac_parser.hpp` (analogo `maxima_parser.hpp`),
  estensione `test/golden/main.cpp` con flag `--oracle giac|maxima`.
- **Algoritmo**: identico runner, normalizzazione output Giac (`x^2`,
  `*`, etc. — Giac syntax ≈ Maxima).
- **Acceptance**: corpus 1026 entry eseguito contro Giac → tabella
  pass-rate per area. Confronto delta con Maxima (entry su cui
  divergono fra loro → SKIP scientifico, non FAIL).

---

### F7.5.H — Final re-measure + F7 declaration

#### F7.5.H1 — Re-measure aggregato
- Esecuzione `cas_golden_runner` su tutte le 11 aree post-A/B/C/D/E/F.
- Output: `F7.5_FINAL_PARITY_REPORT_YYYY-MM-DD.md`.
- **Acceptance**: pass-rate aggregato non-skip ≥ 86%. Se < 86%,
  identificare aree sotto target binding (tabella sopra) e iterare il
  task pertinente fino a raggiungere il floor. No `completed`
  ottimistico.

#### F7.5.H2 — Audit indipendente
- Re-run audit equivalente a `AUDIT_CAS_vs_HP_Prime_2026-06-08.md`,
  questa volta numerico (Maxima + Giac).
- Sign-off utente come gate Fase 8.

---

### Vincoli F7.5 (binding)

1. **Anti-monolito**: ogni nuovo file rispetta 500 LOC. Split
   pre-emptive se l'algoritmo lo richiede.
2. **No test disable** (REGOLA 0.2). Se un entry corpus passa pre-F7.5
   e fallisce post-F7.5, regressione = blocking commit.
3. **No hardcode-of-passage non documentato** (REGOLA ZERO). Ogni
   costante numerica nuova → giustificazione matematica nel commit
   message o entry ledger.
4. **Spec read first** (REGOLA 0.1). Ogni task F7.5.X richiede il
   file `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<Spec>.md`
   scritto **prima** del codice. Plan mode reply deve citare la spec.
5. **Result<T>**, no `throw/catch`. **BigInt only**. **Structural
   sharing**. **AstArena**.
6. **Per-task ledger entries** se introducono `Unimplemented` per casi
   fuori scope (es. inverse trig su campo non-reale).

### Stima totale F7.5

- T1-Sonnet (meccanico): ~6 giorni (A1, A2, A3, F1, G1, G2)
- T2-Sonnet-thinking: ~4.5 settimane (A4, C1, C2, D1, E1, E2)
- T3-Opus (research-grade): ~3-4 settimane (B1, B2, B3)
- Totale calendario: **6-8 settimane** se sequenziale,
  **4-5 settimane** con 2-3 task paralleli (dove indipendenti).

### Exit gate F7.5

- Aggregato corpus ≥ 86%
- Bronstein book corpus 90 entry ≥ 70%
- HC-F70-A43-EXTENDED-REAL chiuso
- Solve + Matrix area ≥ 90%
- Audit indipendente firmato
- Nessun test disabilitato, nessun hardcode non ledgered, nessun file
  > 500 LOC

Post-F7.5 → Fase 8 (research target permanenti `Aperta`: Risch
structure theorem full, Galois ≥6, CAD McCallum, hypergeometric `_pF_q`).

---

## MODEL ROUTING — Orchestrator Opus + Subagent Tier

**Setup obbligatorio**:
- Utente seleziona **Opus 4.7** come modello sessione via `/model`
- Opus 4.7 = **orchestrator permanente**. Mai esegue lavoro coding pesante diretto
- Sottagenti spawnati via `Agent({subagent_type, model, prompt})` con `model: "sonnet"|"opus"` override per-call
- Haiku **escluso** per scelta utente (tier T0 mechanical → Sonnet)

### Tier subagent

| Tier | Modello | Ruolo | Quando usare |
|---|---|---|---|
| **T1-Sonnet** | `model: "sonnet"` (Sonnet 4.6 normal) | Coding canonico + meccanico | Refactor pervasivo, fresh-symbol, ledger append, test boilerplate, algoritmi well-documented con pseudo-codice (Karatsuba, Berlekamp, Hensel, Householder QR, Bareiss, Cholesky), build/test/fix loop su errori chiari, plot/stats setup |
| **T2-Sonnet-thinking** | `model: "sonnet"` con prompt che invoca extended reasoning | Coding con invariant tracking | Brown GCD, Zippel sparse, Trager primitive element, ODE classifier completo, Frobenius log-term, residue grado arbitrario, transforms, branch cuts, RootOf seed fix, simplifier positivity-aware. Prompt forza "ragiona step-by-step, identifica invariant, considera edge cases prima di scrivere" |
| **T3-Opus** | `model: "opus"` (Opus 4.7 extended thinking implicito) | Research-grade + audit critici | Risch structure theorem, Gruntz §3.5 generale, CAD McCallum, Wang multivariate, Galois ≥5, FFT Schönhage-Strassen, audit bi-settimanali, code review architetturale, decisione "Risolta vs Parziale", cas-regression-guard giudizio finale |

### Regole di routing (orchestrator Opus segue queste)

1. **Default subagent**: T1-Sonnet. Escalate a T2-Sonnet-thinking solo se task richiede invariant multi-livello. Escalate a T3-Opus solo per research-grade o audit.
2. **Mai T3-Opus per refactor meccanico** — spreco senza benefit.
3. **Mai T1-Sonnet per task con multi-invariant subtle** — replica problema tracker ottimista.
4. **Audit periodici bi-settimanali**: sempre T3-Opus + prompt esplicito "no ottimismo, cita file:riga, percentuali oneste".
5. **cas-regression-guard finale**: sempre T3-Opus per giudizio pass/fail.

### Ottimizzazione: batching + subagent long-lived + pre-review interlock

**Anti-pattern da evitare**: 100 subagent spawn/dispose. Ogni spawn costa context-loading (~5-15k token re-read file rilevanti). Batching aggressivo riduce overhead 5-10×.

**Pattern raccomandato per ogni fase**:

```
Per fase F_n:
  Opus orchestrator:
    1. Spawn 1 subagent per macro-area indipendente (max 3-5 paralleli)
    2. Ogni subagent riceve BLOCCO di task correlati (5-15 task), non singolo task
    3. Subagent itera sequenziale internamente (build/test/fix loop)
    4. Subagent ritorna report compatto: cosa fatto, cosa skip, file:riga modificati
    5. **Pre-review interlock**: spawn T2-Sonnet-thinking peer-reviewer su diff del subagent
       - Cattura ~50% problemi triviali (hardcode dimenticati, naming, fresh-symbol)
       - Riduce token consumati da T3-Opus audit successivo
       - Riporta a Opus orchestrator solo issue non-banali
    6. Opus valida, decide next batch o escalation a T3-Opus audit per giudizio finale
    7. Riusa stesso subagent via SendMessage per follow-up correlati (no re-spawn)
```

**Cross-phase regression policy**:
- Se subagent in F_n scopre bug strutturale in F_{n-k} foundation (es. F5 integrate scopre `simplify` regression introdotta in F1):
  - **STOP** F_n immediatamente
  - Apri ticket regression
  - ROLLBACK a F_{n-k}, apply fix come PR separata
  - Re-run audit T3-Opus su F_{n-k}
  - Riprende F_n solo dopo F_{n-k} re-certificato

**API stability freeze**:
- A chiusura di ogni fase, public API congelata in `include/cas/*.hpp`
- Modifiche public API in fase successiva → richiedono semver bump + migration note + deprecation 1 release
- Implementazione interna libera; signature pubbliche stabili

**Esempio concreto F4 (linalg)**:

```
Opus spawn T1-Sonnet "F4-block-A":
  prompt: "Refactor Bareiss collapse (matrix_ops.cpp:243 → matrix_bareiss.cpp).
           Fresh-symbol replace _lambda_/t<i>/C<i> in matrix_solve, matrix_eigenvalues,
           matrix_jordan, ode_solver_1st_order, ode_solver_advanced.
           Routing matrix_jordan kernel → null_space_over_extension per autovalori RootOf.
           Fix eigenvectors kernel error swallowing matrix_eigenvalues.cpp:324.
           Build incrementale dopo ogni edit. Test mirati. Report compatto."

Opus spawn T2-Sonnet-thinking "F4-block-B" parallelo:
  prompt: "Implementa Modified Gram-Schmidt QR (sostituisce GS classico matrix_qr.cpp; Householder = Aperta perm HPP-F4.1-QR-HOUSEHOLDER, AST explosion).
           Implementa Cholesky LDL^T (nuovo matrix_cholesky.cpp).
           Implementa Vandermonde + Toeplitz Levinson in matrix_structured_determinant.cpp.
           Implementa Smith Q[x] generalization (matrix_smith.cpp:125 Unimplemented).
           Per ognuno: certificatore + 10 test anti-hardcode + property metamorphic.
           Build + ctest. Report compatto."

Opus attende entrambi. Audit T3-Opus su risultati. Promote a Risolta solo se gates passati.
```

**Risparmio stimato vs 1-task-1-agent**: ~70% token, ~60% wall time, ~80% spawn overhead.

### Parallelismo controllato

- Max **3-5 subagent paralleli** per fase (oltre = merge complexity esplode)
- Subagent paralleli SOLO su file/aree disgiunte (no race condition su stesso modulo)
- Audit sempre **sequenziale** dopo subagent paralleli (Opus aspetta tutti, poi audit)

### Token budget orientativi per fase

| Fase | Spawn count tipico | Token totali stimati | Modello mix |
|---|---|---|---|
| F0 | 1-2 subagent T1-Sonnet | ~100k | 100% Sonnet |
| F1 | 4-5 subagent (2 T1 + 2 T2 + 1 T3 audit) | ~800k | 60% Sonnet, 40% Opus |
| F2 | 3 subagent (1 T1 + 2 T2) + 1 T3 audit | ~500k | 70% Sonnet, 30% Opus |
| F3 | 5 subagent (1 T1 + 2 T2 + 2 T3 per Wang+Galois) + 2 T3 audit | ~1.5M | 40% Sonnet, 60% Opus |
| F4 | 2-3 subagent T1+T2 paralleli + 1 audit | ~400k | 80% Sonnet, 20% Opus (Householder escluso — Aperta perm) |
| F5 | 6-8 subagent (Risch+Gruntz=T3, ODE+residue+transforms=T2, refactor=T1) + 2 audit | ~2M | 30% Sonnet, 70% Opus |
| F6 | 3-4 subagent (CAD=T3, MPFR+branch=T2, units=T1) + 1 audit | ~700k | 50% Sonnet, 50% Opus |
| F7 | 3 subagent T1+T2 + 1 audit final T3 | ~500k | 70% Sonnet, 30% Opus |

**Totale stimato**: ~6.5M token. Costo orientativo: ~$300-500 con tier-routing, vs ~$2000-3000 Opus-everything.

### Gates utente obbligatori per fase

Dopo ogni fase, Opus orchestrator presenta a utente:
1. Riepilogo cosa fatto + cosa skip con motivazione
2. Audit T3-Opus report (file:riga, percentuali oneste, REGOLA ZERO check)
3. Coverage report delta
4. Property/golden suite pass-rate
5. Lista nuovi Unimplemented + ledger updates
6. Richiesta esplicita "approva passaggio a F_{n+1}?"

Senza approvazione utente, Opus NON procede a fase successiva. Anti-drift a "scaffold finto-completo".

---

## CRONOLOGIA STIMATA — Executor = AI (Claude Code-class)

Esecuzione delegata a AI orchestrata da utente con gates di validazione. Timeline ricontestualizzata.

### Modello di esecuzione AI

**Vantaggi AI vs human dev**:
- Coding throughput 5-10× su task ben-specificati (Knuth Algorithm D, Karatsuba, Berlekamp, Bareiss, ecc.: algoritmi canonici con pseudo-codice in letteratura)
- Parallelismo: agent paralleli per task indipendenti (es. MGS QR + Smith Q[x] + Hermite NF in parallelo)
- Test/property generation rapida
- Refactor pervasivo (fresh-symbol, Bareiss collapse) in minuti
- Documentazione/commit message zero-friction

**Limiti AI**:
- Debug loop iterativo (build → test fail → fix → rebuild): cap fisico ~10-15 cicli/ora
- Context window: ogni task richiede ri-leggere file rilevanti; PR grandi vanno spezzate
- Algoritmi research-grade (Risch structure theorem, CAD McCallum, Wang multivariate, Gruntz §3.5 generale): AI tende a produrre MVP convincenti ma incompleti — serve revisione umana spec-vs-impl
- Performance tuning (FFT BigInt, Schönhage-Strassen): correttezza facile, costanti competitive con GMP difficili
- Test e2e contro CAS reference: setup veloce, ma analisi failure cases richiede iterazione umana
- Session quota / cost windows: lavoro a burst, non continuo 24/7
- Drift di astrazione: senza checkpoint utente AI converge verso "scaffold + Unimplemented" se non vincolato (vedi audit attuale)

**Acceleratore richiesto**: gates utente per ogni F-fase. Senza checkpoint umani il rischio è ri-creare lo stesso problema attuale (tracker ottimista, percentuali finte).

### Stima AI-driven con gates utente

| Fase | Effort AI realistico | Note |
|---|---|---|
| F0 sanitizzazione | **1-2 giorni** | T1-Sonnet (1 subagent). Lavoro meccanico: edit tracker, append ledger, setup gcov/rapidcheck/Maxima golden |
| F1 L0 foundation | **2-3 settimane** | T1-Sonnet (BigInt Karatsuba/Toom/Knuth-D); T2-Sonnet-thinking (simplifier_trig generatore + assumption inference + normal_form positivity); T3-Opus (audit). SS FFT escluso — Aperta perm HPP-F1.1-MUL |
| F2 L1 poly univariati | **1-2 settimane** | T1-Sonnet (half-GCD, Berlekamp, Hensel quadratic); T2-Sonnet-thinking (vanHoeij knapsack); T3-Opus audit finale |
| F3 L2 multivar + alg ext | **3-4 settimane** | T2-Sonnet-thinking (Brown GCD, Zippel, primitive element Trager, F5+FGLM); T3-Opus (Wang multivariate, Galois ≥5); T3-Opus audit |
| F4 L2 linalg | **5-7 giorni** | T1-Sonnet (refactor pervasivo + Bareiss collapse + fresh-symbol); T2-Sonnet-thinking (MGS QR, Cholesky, Smith Q[x], Jordan routing); T3-Opus audit. Householder QR escluso — Aperta perm HPP-F4.1-QR-HOUSEHOLDER |
| F5 L2 calculus | **4-6 settimane** | T2-Sonnet-thinking (ODE classifier+Frobenius+Padé+transforms+Laurent); T3-Opus (Risch structure theorem, Gruntz §3.5, residue grado arbitrario, Zeilberger, hypergeometric); T3-Opus audit ogni 2 sett |
| F6 L3 numerica+complex+CAD+units | **2-3 settimane** | T1-Sonnet (units SI, RootOf seed fix); T2-Sonnet-thinking (MPFR unified, branch cuts); T3-Opus (CAD McCallum); T3-Opus audit |
| F7 plotting + stats + acceptance | **1-2 settimane** | T1-Sonnet (plot adaptive, stats distributions); T2-Sonnet-thinking (Gauss-Kronrod, spline); T3-Opus acceptance HP Prime corpus + bugfix |

**Totale AI-driven**: **~3-4 mesi calendario reali** con gates utente, build/test loop, debug iterativo.

### Distribuzione effort

- ~50% del tempo AI = build/test/fix loop (non coding pure)
- ~20% = analisi failure cases + debug edge case
- ~15% = coding nuovo
- ~10% = refactor + cleanup + ledger
- ~5% = documentazione

### Sliding parallelism

Più agent in parallelo per task indipendenti accorcia clock-time ma aumenta token cost e merge complexity. Esempio: F4 (linalg) interamente parallelizzabile su 4 agent → 2 giorni invece di 5-7. F5 (calculus) parallelizzabile a fasi: Risch, ODE, residue, summation indipendenti → potenziale -40% clock time.

Realistico con parallelismo aggressivo: **2-3 mesi calendario**.

### Punti critici dove AI fallisce regolarmente (richiedono escalation)

1. **Risch structure theorem completo (Bronstein cap 9)**: spec ambigua in letteratura, AI converge su MVP. Mitigazione: spezzare in 6-8 PR incrementali con validazione utente ogni step.
2. **CAD McCallum**: complessità doppia esponenziale, AI tende a stub silenziosi. Mitigazione: ridurre scope a 2-var deg ≤4 confermato dal utente come MVP onesto.
3. **Wang multivariate factorization**: leading coefficient determination + bad zero detection sono trappole classiche. Mitigazione: corpus 100 input random pre-marker Risolta, no shortcut.
4. **FFT BigInt competitiva con GMP**: AI può scrivere FFT corretto ma 5-10× più lento. Mitigazione: accettare fallback `--use-gmp` come opzione runtime se performance non target.
5. **Hypergeometric `_pF_q` closed-form recognition**: tabelle Wilf-Zeilberger sono vaste, AI non le sa tutte. Mitigazione: codificare solo Gauss 2F1(a,b;c;1) + Saalschütz + casi ovvi, marker research target per resto.

### Quality gates obbligatori (anti-regressione AI)

Per evitare ripetizione del problema attuale (tracker ottimista, codice MVP):

- Ogni "Risolta" claim AI deve essere validata da **cas-regression-guard** + property tests + golden Maxima ≥20 input
- Ogni `Unimplemented` aggiunto deve iscriversi in `HARDCODE_LEDGER.md` con motivazione esplicita
- Ogni PR deve avere **diff coverage ≥85%** dei nuovi line
- Ogni 2 settimane: audit critico utente vs tracker (come quello fatto oggi)
- Sanity check finale F7: re-run audit indipendente; pass-rate HP Prime corpus ≥95% obbligatorio

---

## REGOLE DI INGAGGIO PER OGNI TASK

Per ogni task in ogni fase, applicare il protocollo standard:

0. **LETTURA SPECIFICA OBBLIGATORIA (CLAUDE.md REGOLA 0.1)**: Verificare e leggere il file `.md` corrispondente in `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/` prima di scrivere qualsiasi codice.
1. **Contratto matematico**: dominio, output, casi esclusi, invarianti, errori
2. **Algoritmo canonico** con citazione bibliografica (autore, anno, pagina)
3. **Certificatore indipendente** dal path principale
4. **Anti-hardcode tests** (vedere F2 minimi)
5. **Coverage gate**: ≥90% line coverage + ≥70% mutation score per il modulo
6. **Property-based**: almeno 1 property metamorphic
7. **Golden vs Maxima**: ≥20 input confrontati (≥50 per algoritmi top-tier)
8. **Performance gate**: benchmark vs baseline (regression <10%) o vs Giac/Maxima
9. **0 nuovi hardcode** non in HARDCODE_LEDGER.md
10. **Tracker update** con stato onesto e evidenza file:riga
11. **Doc algoritmo** compilato in `.APROJECT_REFERENCES/algorithms/<nome>.md`
12. **Error diagnostic** strutturato su ogni Unimplemented path
13. **Anti-monolith**: nessun file >500 righe (CI gate)
14. **Public API stability** rispettata (no breakage senza semver bump)
15. **Status onesto**: `Risolta` solo se dominio dichiarato 100% chiuso; altrimenti `Parziale avanzata` con subset esplicito nel nome task

Una pull request che viola anche solo 1 di questi 15 punti viene rifiutata.

---

## VINCOLI INVARIABILI

- BigInt only, no `int64_t`/`double` nel core simbolico (REGOLA 1)
- Structural sharing (REGOLA 2)
- Memory Arena bump (REGOLA 3)
- LPO orientation rewrite rules
- Result<T> monadico, no throw/catch
- 500 righe max per file sorgente
- Fresh-symbol pervasivo, mai literal names

---

## DELIVERABLE PER FASE

Ogni fase chiude con:
1. PR singola (o serie) referenziata in CAS_TASKS.md
2. Aggiornamento HARDCODE_LEDGER.md
3. Bench report numerico vs baseline
4. Coverage report HTML
5. Golden suite pass-rate
6. Test count delta
7. Documento `PHASE_N_CLOSURE.md` con: cosa fatto, cosa NON fatto, evidence file:riga

---

## RISCHI E MITIGAZIONI

| Rischio | Mitigazione |
|---|---|
| FFT BigInt non competitivo GMP | Accettare 30-50% GMP onesto come `Risolta`; gap permanente marcato `Aperta` "FFT vs GMP competitivo" research; NO wrap GMP (REGOLA ZERO) |
| Risch structure theorem grande effort | Spezzare in incrementi: log-extension prima, exp dopo, mix dopo. Ogni step richiede dominio dichiarato esplicito + certificatore D(F)=f |
| Brown GCD ricorsione esplosione | Test su corpus 1000 input random pre-marker Risolta. Property `gcd·lcm=a·b` su 10000 random Z[x,y,z] |
| CAD doppia esponenziale | Spostata fuori scope parità HP Prime: Fase 8 post-parità come strato L4. F6 implementa solo disequazioni 1-var via Sturm (algoritmo intero, no subset) |
| Galois deg ≥5 effort enorme | F3 chiude deg 5 (S5/A5/D5/F20/C5) come `Risolta`. Deg ≥6 marker `Aperta` research permanente |
| Performance regressione cumulativa | Bench gate obbligatorio ogni PR + dashboard tracking |
| Cross-phase regression | Policy ROLLBACK obbligatoria se F_n rivela bug in F_{n-k}; no shortcut layering |
| Subagent scaffold finto-completo | Pre-review interlock T2-Sonnet-thinking + audit T3-Opus + gates utente obbligatori |
| API breakage cross-fase | API freeze a chiusura ogni fase; semver bump + deprecation per modifiche |
| Build time esplosione su refactor | ccache + Ninja + split file >500 righe via CI gate |

---

## STATUS PERMANENTI `Aperta` POST-PARITÀ

Onestamente accettati come research target oltre F7. Parità HP Prime ≥95% non richiede chiusura di questi:

| Voce | Motivazione | Strategia |
|---|---|---|
| Risch structure theorem full Bronstein cap 9 | 600pp libro, multi-mese effort math expert | F5.1 chiude subset Liouville+Hermite+Trager esplicito; structure theorem completo = Fase 8 |
| Galois group deg ≥6 generale | Algoritmi research-grade (Cohen cap 6 + Magma papers) | F3.6 chiude deg ≤5; deg ≥6 = Fase 8 |
| CAD McCallum generale multivar | Complessità doppia esponenziale, scope onesto incompatibile L3 | F6.4 chiude solo 1-var Sturm; CAD = Fase 8 |
| FFT BigInt competitivo GMP (>50%) | GMP = 30 anni tuning low-level assembly | F1.1 chiude 30-50% GMP onesto; tuning oltre = Fase 8 |
| Schönhage-Strassen BigInt FFT (HPP-F1.1-MUL) | n≥4096 limbs fallback Karatsuba; SS NTT richiede modular FFT prime-finding + bit-reversal stable in BigInt esatto, multi-mese effort | F1 chiude su Karatsuba+Toom-3; SS = Fase 8 research |
| Householder QR simbolico (HPP-F4.1-QR-HOUSEHOLDER) | Riflettori `I - 2vv^T/v^Tv` esplodono coefficienti razionali AST per matrici ≥8×8 simboliche; Modified Gram-Schmidt è la scelta corretta per CAS esatto | F4 chiude su MGS (Trefethen-Bau §8); Householder numerico = Fase 8 con AlgebraicNumber tower stabile |
| Hypergeometric `_pF_q` recognition completo | Tabelle Wilf-Zeilberger vaste, knowledge-heavy | F5.9 chiude Gauss 2F1, Saalschütz, casi closed-form noti; recognition completo = Fase 8 |
| Multi-sheet Riemann surface complex | Modello AST single-valued incompatibile | F6.2 chiude branch principale; multi-sheet = Fase 9 ricerca |

Tutte iscritte in CAS_TASKS.md come `Aperta` permanente con motivazione esplicita. Nessuna pretesa di chiusura.

---

## CONCLUSIONE

Stato attuale (audit critico): ~40% di CAS production-grade.
Stato target: 100% parità HP Prime G2 (su domini dichiarati).
Strategia: bottom-up stratificato con gates rigorosi, certificatori indipendenti, anti-hardcode pervasivo, tracker onesto.

**Tempo AI-driven realistico**: **3-4 mesi** calendario seriale, **2-3 mesi** con parallelismo aggressivo (max 3-5 subagent), condizionato a:
- Gates utente ogni fase (no AI auto-promotion a "Risolta")
- Build/test/fix loop iterativo onesto
- Refusal di shortcut MVP non documentati
- Audit indipendente bi-settimanale T3-Opus (replica di quello fatto oggi)
- Accettazione che 2-3 algoritmi top-tier (Risch full, CAD, Galois ≥6) restino research target oltre F7

Senza gates utente: AI converge a "scaffold finto-completo" in 2-3 settimane, replicando il problema attuale a scala più grande.

### Setup operativo utente

1. **Sessione**: `/model opus` (Opus 4.7 orchestrator permanente)
2. **Apri questo file**: orchestrator legge `PLAN_HP_PRIME_PARITY.md` come spec
3. **Conferma fase di partenza**: utente dice "esegui F0" o "esegui F1 block A"
4. **Orchestrator Opus**: spawna subagent con tier corretto (T1/T2/T3) seguendo MODEL ROUTING
5. **Gates**: dopo ogni fase, Opus presenta report + chiede approvazione esplicita per fase successiva
6. **Audit bi-settimanale**: utente o orchestrator triggera audit T3-Opus indipendente, confronta con tracker, retro-classifica voci ottimistiche

### Costo orientativo totale

- Tier-routing (questo piano): **~$300-500** in 3-4 mesi
- Confronto Opus-everywhere: ~$2000-3000 stesso periodo
- Risparmio: **~80%** senza perdita qualità (tier-routing è ottimale, non compromesso)
