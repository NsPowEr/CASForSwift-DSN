# DIVIETO HARDCODE — Catalogo Anti-Pattern (dettaglio)

> Dettaglio normativo delle 10 categorie di hardcode vietate. L'**indice** vive
> in `CLAUDE.md` (sezione "DIVIETO HARDCODE"); questo file contiene tabelle ed
> esempi completi. Carica questo file quando lavori sulla categoria pertinente.
>
> **Regola assoluta**: un valore hardcoded in un algoritmo matematico è un bug
> latente. L'unica eccezione sono le costanti matematiche esatte (π, e, φ) e i
> default configurabili tramite `CASContext`.

## Principio Fondamentale

Un CAS professionale differisce da un prototipo didattico perché applica **algoritmi universali** invece di **tavole lookup o soglie empiriche**. Ogni costante che non può essere derivata matematicamente dalla struttura del problema è un hardcode vietato. Se il valore "funziona" solo per i test presenti ma fallirebbe su input più grandi/complessi, è un hardcode.

**Test di autodichiarazione obbligatorio prima di ogni commit:**
```
- [ ] Ogni costante numerica nel codice ha una giustificazione matematica formale (bound, formula, teorema).
- [ ] Ogni limite computazionale è configurabile via CASContext (non solo tramite ricompilazione).
- [ ] Nessun set fisso di "casi speciali" che esclude silenziosamente input validi.
- [ ] Nessun fallback a Unimplemented causato da un parametro fisso superabile con più iterazioni/campioni.
```

---

## Categoria 1: Budget Computazionali Non Configurabili

**VIETATO**: Limiti di profondità, iterazioni o ricorsione hardcoded che troncano calcoli matematicamente validi.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `if (depth >= 16U) return Unimplemented` | `integrate_core.cpp:17` | Blocca `∫x^n*e^x dx` con n>14 | `ctx.max_integration_depth` configurabile, default 32 |
| `constexpr int MAX_SIMPLIFICATION_DEPTH = 300` | `simplify_impl.hpp:17` | Blocca det(A) simbolico 5×5 | Configurabile + distinzione ciclo/profondità legittima |
| `if (subset_count >= 32768) return nullopt` | `factorization_recombination.cpp:105` | Timeout-count, non timeout-tempo | Timeout temporale configurabile O non conta subset |
| `max_steps` Hensel calcolato ma con bail-out fisso | `polynomial_gcd_multivariate.cpp:539` | GCD multivariato complesso → Unimplemented | Esporre `max_hensel_steps` in CASContext |

**Regola**: Ogni limite computazionale deve essere: (a) configurabile via `CASContext`, (b) accompagnato da un `Unimplemented` **esplicito e diagnostico** quando raggiunto (non un risultato silenziosamente sbagliato).

---

## Categoria 2: Costanti Magiche in Algoritmi Algebrici

**VIETATO**: Costanti nei calcoli algebrici senza derivazione matematica formale.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `BigInt B = 2*max_coeff + 100; B *= 1000` | `polynomial_gcd_heuristic.cpp:146-151` | Per coeff > 10^6 la collisione è probabile | Bound di Mignotte: `B ≥ 2^deg * max_coeff` |
| `score = 1000` (integer), `500 - min(400, cplx)` (simbolico) | `matrix_bareiss.cpp:110-115` | Ignora assumptions su segno/dominio | Score da `is_known_nonzero(assumptions)` |
| `max_samples = required_samples + 8U` | `polynomial_gcd_multivariate.cpp:741` | "+8" arbitrario, nessun livello di confidenza | `ceil(log(δ)/log(1-p_hit))` con δ configurabile |
| `delta_val = 0.75` come unico default LLL | `polynomial_internal.hpp:201` | Non ottimale per tutti i reticoli | Già configurabile; non ripristinare default hardcoded |

**Regola**: Le costanti in algoritmi probabilistici o euristici devono derivare da un **bound matematico** (Mignotte, Hadamard, Schwartz-Zippel) o da un **parametro di confidenza** espresso esplicitamente.

---

## Categoria 3: Set e Range di Ricerca Fissi

**VIETATO**: Insiemi finiti che escludono silenziosamente input matematicamente validi oltre la soglia.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `kCandidates[] = {13,17,...,47}` con fallback a 13 | `factorization_integers.cpp:232` | Fallback deterministico compromette randomizzazione | Hash del polinomio come seed per selezione ciclica |
| `for (int n = 1; n <= 100; ++n)` ciclotomici | `polynomial_cyclotomic.cpp:87` | `x^101 - 1 = 0` non riconosciuto come ciclotomico | Limite come parametro; generazione on-demand per n arbitrario |
| Tabella hardcoded `sin/cos/exp` in serie di Taylor | `limit_series.cpp:190-237` | `tan(x)`, `W(x)`, composizioni → Unimplemented | Generatore sistematico via derivate successive |
| Riduzione trig solo per `k*π/12` | `simplify_functions.cpp:231` | `sin(π/5)`, `sin(π/17)` ignorati | Polinomi di Chebyshev o algoritmo di riduzione generale |

**Regola**: Se un set fisso `{a₁, a₂, ..., aₙ}` esiste nel codice, deve (a) essere esaustivo per tutti i casi matematicamente possibili, oppure (b) avere un generatore algoritmico per i casi rimanenti. Non è mai accettabile il silenzio (fallire senza segnalare che il caso è fuori range).

---

## Categoria 4: Bail-out su Tipo di Dato

**VIETATO**: Rifiuto di input valido basato sul tipo invece che sul dominio matematico.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `if (!expr_cast<IntegerLit>(value)) return Unimplemented` | `algebra_core.cpp:311` | Eval su Q o RootOf è matematicamente valida | Estendere a `IntegerLit \| RationalLit \| RootOf` |
| `if (expr_is<DecimalLit>(expr)) return Unimplemented` in diff/integrate | `differentiate.cpp:114`, `integrate_core.cpp:142` | `0.5*x` è differenziabile | Conversione DecimalLit→Rational al parser (non nel core) |
| `return Unimplemented("solo sqrt(n) per ora")` | `factorization_polynomials.cpp:649` | Fattorizzazione su `Q(∛2)` è valida | Dispatcher su tipo estensione; stub onesto per non-sqrt |
| Uso di `int64_t` o `double` nel core simbolico | (vietato da Regola 1) | Overflow silenzioso o errore numerici | Solo `BigInt` e `Rational` |

**Regola**: Il rifiuto di un tipo deve essere **esplicito, diagnostico, e temporaneo** (marcato con un task aperto in CAS_TASKS.md). Non è mai accettabile un Unimplemented silenzioso che non spiega cosa manca.

---

## Categoria 5: Ordinamenti e Strutture Fisse Non Configurabili

**VIETATO**: Scelte strutturali che bloccano algoritmi alternativi corretti e più efficienti.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `MonomialLexComparator` forzato in F4 | `polynomial_groebner_f4.cpp` | `grevlex` è drasticamente più efficiente per Groebner | Ordinamento come parametro via `MonomialOrder` enum (→ L3-20) |
| Unico modulo di Mersenne (2³¹-1) senza accumulazione | Già corretto con CRT | — | Già implementato correttamente — non regredire |
| RREF con pivot fisso senza selezione contestuale | `matrix_bareiss.cpp` | Pivoting non-numerico su simbolici richiede euristica domain-aware | Scoring con assumptions (→ L1-17) |

---

## Categoria 6: Seed, Randomness e Parametri Probabilistici

**VIETATO**: Sorgenti di randomness deterministica non derivata dall'input o non configurabile.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| Seed fisso `42` o simili in Cantor-Zassenhaus | `factorization_integers.cpp` | Stessa sequenza → stesse collisioni su stessi input | Hash del polinomio come seed; seed esplicito in CASContext |
| Fallback fisso a primo `p=13` | `factorization_integers.cpp:239` | Se lc divisibile da 13, degrado algoritmo | Selezione ciclica hash-based su pool esteso |
| `+8` campioni extra fissi in GCD | `polynomial_gcd_multivariate.cpp:741` | Nessuna garanzia probabilistica formale | `ceil(log(δ)/log(1-p_hit))` con δ = `ctx.gcd_error_probability` |

---

## Categoria 7: Nomi di Variabili Interni Hardcoded

**VIETATO**: Nomi "magici" per variabili ausiliarie interne che possono collidere con simboli utente.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `std::string w_name = "__mrv_w"` | `limit_mrv.cpp:540` | Anche con collision check, è fragile; utente potrebbe usare `__mrv_w` | Generatore di simboli freschi via `ctx.make_fresh_symbol("mrv")` che garantisce unicità |
| Nomi fissi `C1`, `C2`, ... per costanti ODE | `ode_solver_advanced.cpp` | Collisione se utente ha già `C1` in scope | `ctx.make_fresh_symbol("C")` → `C_1`, `C_2` garantiti unici |
| Variabili `__cas_internal_*` non generate freschie | qualsiasi | Un simbolo fisso può apparire in output utente | Sempre `make_fresh_symbol` con prefisso; mai nomi letterali |

---

## Categoria 8: Pattern Matching a Tabella Chiusa

**VIETATO**: Dispatch su forma dell'espressione tramite lista fissa di pattern che non scala.

| Pattern vietato | Perché sbagliato | Alternativa |
|---|---|---|
| `if (is_exp_x) return exp_x; if (is_ln_x) return ln_x; ...` per Risch | Copertura < 20% dell'algoritmo; ogni caso aggiunto manualmente | Implementazione formale del differential field e riduzione algoritmica |
| Integrazione via lista fissa di forme riconosciute | Su input non in lista → Unimplemented, anche se integrabile | Risch + Hermite + LRT come pipeline algoritmica completa |
| Serie Taylor via formule separate per sin/cos/exp | `tan(x)` fuori lista → Unimplemented | Generatore sistematico via derivazione ripetuta |
| ODE: `if (type == Type1) solver1; if (type == Type2) solver2` | Fallisce su ODE che non rientrano nei tipi previsti | Classificatore extensible con fallback Unimplemented diagnostico |

**Regola**: Il pattern matching su forma è accettabile **solo** come ottimizzazione (fast-path) per casi comuni, **mai** come unico algoritmo. Deve sempre esistere un path algoritmico generale.

---

## Categoria 9: Intervalli di Controllo e Polling Fissi

**VIETATO**: Intervalli fissi per operazioni che devono adattarsi al costo computazionale reale.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `kTimeoutCheckInterval = 1024U` fisso | `symbolic_internal.hpp:49` | Su Groebner pesante, 1024 ops = secondi → timeout non reattivo | Configurabile in CASContext; adattivo per algoritmo pesante (64-256) |
| Check timeout ogni N operazioni senza stima costo | qualsiasi | "Operazione" ha costo variabile (1 vs milioni di cicli) | Campionamento temporale: check ogni `min(1024, N)` dove N stima costo |

---

## Categoria 10: Gerarchie di Crescita e Rank Statici

**VIETATO**: Classificazioni statiche che assegnano lo stesso rank a oggetti matematicamente distinti.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `GrowthRank` 0/1/2/3 assegnato staticamente | `limit_infinite.cpp` | `e^x` e `e^(e^x)` hanno stesso rank 3 → Cancellation Tower errata | Rank calcolato dinamicamente via confronto asintotico ricorsivo (algoritmo Gruntz) |
| Comparazione `compare_growth()` su struttura senza coefficienti leader | `limit_mrv.cpp:115-160` | `2*e^x` vs `e^x` indistinguibili se si ignora il coefficiente | Calcolo esplicito del coefficiente leader per comparabili stesso ordine |

---

## Eccezioni Legittime (Hardcode Accettabili)

Questi hardcode **sono permessi** perché matematicamente fondati o architetturalmente necessari:

1. **Costanti matematiche esatte**: `π`, `e`, `φ = (1+√5)/2`, `γ` (Eulero-Mascheroni) — sono oggetti matematici, non parametri.
2. **Casi speciali del simplifier** derivati da identità provate: `sin(0)=0`, `exp(0)=1`, `log(1)=0` — non sono "tabella" ma assiomi.
3. **Default configurabili in CASContext** con documentazione esplicita: `default_integration_depth = 32` è accettabile se esposto come `ctx.max_integration_depth` e documentato.
4. **Limiti di sicurezza hardware**: `MAX_BIGINT_LIMBS = 10000` per evitare OOM — accettabile se causa `Unimplemented` esplicito con messaggio diagnostico.
5. **Bound CRT**: il primo primo `2^31-1` come punto di partenza del CRT multi-prime è accettabile perché viene accumulato dinamicamente fino al bound di Hadamard.
