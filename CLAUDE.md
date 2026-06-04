# CAS Engine — Costituzione Tecnica per Agenti e Sviluppatori

> Questo file è la **Legge Suprema** del progetto. Ogni agente DEVE leggerlo e rispettarlo.
> Ignorare queste regole porterà al rigetto immediato del codice (Architectural Fail).

---

## REGOLA ZERO — DIVIETO DELLA VIA FACILE

> **MAI E POI MAI** scegliere la via più facile o più veloce solo per comodità.

Ogni decisione di design (algoritmo, struttura dati, condizione di terminazione,
gestione del tipo, scope dei casi) deve essere presa **scegliendo la soluzione
matematicamente corretta e generale**, anche se richiede più codice, più tempo,
o un investimento prerequisito. Se la soluzione completa richiede un'infrastruttura
non ancora presente, due opzioni sono ammissibili:

1. **Costruire prima il prerequisito**, poi tornare al problema (preferita).
2. **Implementare una versione di passaggio chiaramente segnalata** come tale,
   con:
   - Commento `// HARDCODE-OF-PASSAGE: <descrizione>` nel codice.
   - Iscrizione obbligatoria nel file `HARDCODE_LEDGER.md` con: id, file, riga,
     descrizione, categoria CLAUDE.md, fix corretto, blocking dependency.
   - Diagnostico `Unimplemented` esplicito per casi fuori scope (mai silenzio).

**Vietate** in particolare:
- Lookup tables con set chiusi quando esiste un algoritmo generale.
- Scorciatoie "funzionano per il test presente" senza copertura del caso generale.
- Pattern matching su forma quando esiste un algoritmo strutturale.
- Bail-out su tipo o range invece di estendere il dispatcher.
- Costanti magiche al posto di parametri configurabili in `CASContext`.

Il commit messaggio deve dichiarare ogni hardcode-of-passage introdotto.
Una pull request che lascia hardcode non documentato in `HARDCODE_LEDGER.md`
viene rifiutata.

---

## DIVIETO HARDCODE — Catalogo Anti-Pattern

> **Regola assoluta**: Un valore hardcoded in un algoritmo matematico è un bug latente.
> L'unica eccezione sono le costanti matematiche esatte (π, e, φ) e i default configurabili tramite `CASContext`.

### Principio Fondamentale

Un CAS professionale differisce da un prototipo didattico perché applica **algoritmi universali** invece di **tavole lookup o soglie empiriche**. Ogni costante che non può essere derivata matematicamente dalla struttura del problema è un hardcode vietato. Se il valore "funziona" solo per i test presenti ma fallirebbe su input più grandi/complessi, è un hardcode.

**Test di autodichiarazione obbligatorio prima di ogni commit:**
```
- [ ] Ogni costante numerica nel codice ha una giustificazione matematica formale (bound, formula, teorema).
- [ ] Ogni limite computazionale è configurabile via CASContext (non solo tramite ricompilazione).
- [ ] Nessun set fisso di "casi speciali" che esclude silenziosamente input validi.
- [ ] Nessun fallback a Unimplemented causato da un parametro fisso superabile con più iterazioni/campioni.
```

---

### Categoria 1: Budget Computazionali Non Configurabili

**VIETATO**: Limiti di profondità, iterazioni o ricorsione hardcoded che troncano calcoli matematicamente validi.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `if (depth >= 16U) return Unimplemented` | `integrate_core.cpp:17` | Blocca `∫x^n*e^x dx` con n>14 | `ctx.max_integration_depth` configurabile, default 32 |
| `constexpr int MAX_SIMPLIFICATION_DEPTH = 300` | `simplify_impl.hpp:17` | Blocca det(A) simbolico 5×5 | Configurabile + distinzione ciclo/profondità legittima |
| `if (subset_count >= 32768) return nullopt` | `factorization_recombination.cpp:105` | Timeout-count, non timeout-tempo | Timeout temporale configurabile O non conta subset |
| `max_steps` Hensel calcolato ma con bail-out fisso | `polynomial_gcd_multivariate.cpp:539` | GCD multivariato complesso → Unimplemented | Esporre `max_hensel_steps` in CASContext |

**Regola**: Ogni limite computazionale deve essere: (a) configurabile via `CASContext`, (b) accompagnato da un `Unimplemented` **esplicito e diagnostico** quando raggiunto (non un risultato silenziosamente sbagliato).

---

### Categoria 2: Costanti Magiche in Algoritmi Algebrici

**VIETATO**: Costanti nei calcoli algebrici senza derivazione matematica formale.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `BigInt B = 2*max_coeff + 100; B *= 1000` | `polynomial_gcd_heuristic.cpp:146-151` | Per coeff > 10^6 la collisione è probabile | Bound di Mignotte: `B ≥ 2^deg * max_coeff` |
| `score = 1000` (integer), `500 - min(400, cplx)` (simbolico) | `matrix_bareiss.cpp:110-115` | Ignora assumptions su segno/dominio | Score da `is_known_nonzero(assumptions)` |
| `max_samples = required_samples + 8U` | `polynomial_gcd_multivariate.cpp:741` | "+8" arbitrario, nessun livello di confidenza | `ceil(log(δ)/log(1-p_hit))` con δ configurabile |
| `delta_val = 0.75` come unico default LLL | `polynomial_internal.hpp:201` | Non ottimale per tutti i reticoli | Già configurabile; non ripristinare default hardcoded |

**Regola**: Le costanti in algoritmi probabilistici o euristici devono derivare da un **bound matematico** (Mignotte, Hadamard, Schwartz-Zippel) o da un **parametro di confidenza** espresso esplicitamente.

---

### Categoria 3: Set e Range di Ricerca Fissi

**VIETATO**: Insiemi finiti che escludono silenziosamente input matematicamente validi oltre la soglia.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `kCandidates[] = {13,17,...,47}` con fallback a 13 | `factorization_integers.cpp:232` | Fallback deterministico compromette randomizzazione | Hash del polinomio come seed per selezione ciclica |
| `for (int n = 1; n <= 100; ++n)` ciclotomici | `polynomial_cyclotomic.cpp:87` | `x^101 - 1 = 0` non riconosciuto come ciclotomico | Limite come parametro; generazione on-demand per n arbitrario |
| Tabella hardcoded `sin/cos/exp` in serie di Taylor | `limit_series.cpp:190-237` | `tan(x)`, `W(x)`, composizioni → Unimplemented | Generatore sistematico via derivate successive |
| Riduzione trig solo per `k*π/12` | `simplify_functions.cpp:231` | `sin(π/5)`, `sin(π/17)` ignorati | Polinomi di Chebyshev o algoritmo di riduzione generale |

**Regola**: Se un set fisso `{a₁, a₂, ..., aₙ}` esiste nel codice, deve (a) essere esaustivo per tutti i casi matematicamente possibili, oppure (b) avere un generatore algoritmico per i casi rimanenti. Non è mai accettabile il silenzio (fallire senza segnalare che il caso è fuori range).

---

### Categoria 4: Bail-out su Tipo di Dato

**VIETATO**: Rifiuto di input valido basato sul tipo invece che sul dominio matematico.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `if (!expr_cast<IntegerLit>(value)) return Unimplemented` | `algebra_core.cpp:311` | Eval su Q o RootOf è matematicamente valida | Estendere a `IntegerLit \| RationalLit \| RootOf` |
| `if (expr_is<DecimalLit>(expr)) return Unimplemented` in diff/integrate | `differentiate.cpp:114`, `integrate_core.cpp:142` | `0.5*x` è differenziabile | Conversione DecimalLit→Rational al parser (non nel core) |
| `return Unimplemented("solo sqrt(n) per ora")` | `factorization_polynomials.cpp:649` | Fattorizzazione su `Q(∛2)` è valida | Dispatcher su tipo estensione; stub onesto per non-sqrt |
| Uso di `int64_t` o `double` nel core simbolico | (vietato da Regola 1) | Overflow silenzioso o errore numerici | Solo `BigInt` e `Rational` |

**Regola**: Il rifiuto di un tipo deve essere **esplicito, diagnostico, e temporaneo** (marcato con un task aperto in CAS_TASKS.md). Non è mai accettabile un Unimplemented silenzioso che non spiega cosa manca.

---

### Categoria 5: Ordinamenti e Strutture Fisse Non Configurabili

**VIETATO**: Scelte strutturali che bloccano algoritmi alternativi corretti e più efficienti.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `MonomialLexComparator` forzato in F4 | `polynomial_groebner_f4.cpp` | `grevlex` è drasticamente più efficiente per Groebner | Ordinamento come parametro via `MonomialOrder` enum (→ L3-20) |
| Unico modulo di Mersenne (2³¹-1) senza accumulazione | Già corretto con CRT | — | Già implementato correttamente — non regredire |
| RREF con pivot fisso senza selezione contestuale | `matrix_bareiss.cpp` | Pivoting non-numerico su simbolici richiede euristica domain-aware | Scoring con assumptions (→ L1-17) |

---

### Categoria 6: Seed, Randomness e Parametri Probabilistici

**VIETATO**: Sorgenti di randomness deterministica non derivata dall'input o non configurabile.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| Seed fisso `42` o simili in Cantor-Zassenhaus | `factorization_integers.cpp` | Stessa sequenza → stesse collisioni su stessi input | Hash del polinomio come seed; seed esplicito in CASContext |
| Fallback fisso a primo `p=13` | `factorization_integers.cpp:239` | Se lc divisibile da 13, degrado algoritmo | Selezione ciclica hash-based su pool esteso |
| `+8` campioni extra fissi in GCD | `polynomial_gcd_multivariate.cpp:741` | Nessuna garanzia probabilistica formale | `ceil(log(δ)/log(1-p_hit))` con δ = `ctx.gcd_error_probability` |

---

### Categoria 7: Nomi di Variabili Interni Hardcoded

**VIETATO**: Nomi "magici" per variabili ausiliarie interne che possono collidere con simboli utente.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `std::string w_name = "__mrv_w"` | `limit_mrv.cpp:540` | Anche con collision check, è fragile; utente potrebbe usare `__mrv_w` | Generatore di simboli freschi via `ctx.make_fresh_symbol("mrv")` che garantisce unicità |
| Nomi fissi `C1`, `C2`, ... per costanti ODE | `ode_solver_advanced.cpp` | Collisione se utente ha già `C1` in scope | `ctx.make_fresh_symbol("C")` → `C_1`, `C_2` garantiti unici |
| Variabili `__cas_internal_*` non generate freschie | qualsiasi | Un simbolo fisso può apparire in output utente | Sempre `make_fresh_symbol` con prefisso; mai nomi letterali |

---

### Categoria 8: Pattern Matching a Tabella Chiusa

**VIETATO**: Dispatch su forma dell'espressione tramite lista fissa di pattern che non scala.

| Pattern vietato | Perché sbagliato | Alternativa |
|---|---|---|
| `if (is_exp_x) return exp_x; if (is_ln_x) return ln_x; ...` per Risch | Copertura < 20% dell'algoritmo; ogni caso aggiunto manualmente | Implementazione formale del differential field e riduzione algoritmica |
| Integrazione via lista fissa di forme riconosciute | Su input non in lista → Unimplemented, anche se integrabile | Risch + Hermite + LRT come pipeline algoritmica completa |
| Serie Taylor via formule separate per sin/cos/exp | `tan(x)` fuori lista → Unimplemented | Generatore sistematico via derivazione ripetuta |
| ODE: `if (type == Type1) solver1; if (type == Type2) solver2` | Fallisce su ODE che non rientrano nei tipi previsti | Classificatore extensible con fallback Unimplemented diagnostico |

**Regola**: Il pattern matching su forma è accettabile **solo** come ottimizzazione (fast-path) per casi comuni, **mai** come unico algoritmo. Deve sempre esistere un path algoritmico generale.

---

### Categoria 9: Intervalli di Controllo e Polling Fissi

**VIETATO**: Intervalli fissi per operazioni che devono adattarsi al costo computazionale reale.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `kTimeoutCheckInterval = 1024U` fisso | `symbolic_internal.hpp:49` | Su Groebner pesante, 1024 ops = secondi → timeout non reattivo | Configurabile in CASContext; adattivo per algoritmo pesante (64-256) |
| Check timeout ogni N operazioni senza stima costo | qualsiasi | "Operazione" ha costo variabile (1 vs milioni di cicli) | Campionamento temporale: check ogni `min(1024, N)` dove N stima costo |

---

### Categoria 10: Gerarchie di Crescita e Rank Statici

**VIETATO**: Classificazioni statiche che assegnano lo stesso rank a oggetti matematicamente distinti.

| Pattern vietato | File (esempio) | Perché sbagliato | Alternativa |
|---|---|---|---|
| `GrowthRank` 0/1/2/3 assegnato staticamente | `limit_infinite.cpp` | `e^x` e `e^(e^x)` hanno stesso rank 3 → Cancellation Tower errata | Rank calcolato dinamicamente via confronto asintotico ricorsivo (algoritmo Gruntz) |
| Comparazione `compare_growth()` su struttura senza coefficienti leader | `limit_mrv.cpp:115-160` | `2*e^x` vs `e^x` indistinguibili se si ignora il coefficiente | Calcolo esplicito del coefficiente leader per comparabili stesso ordine |

---

### Eccezioni Legittime (Hardcode Accettabili)

Questi hardcode **sono permessi** perché matematicamente fondati o architetturalmente necessari:

1. **Costanti matematiche esatte**: `π`, `e`, `φ = (1+√5)/2`, `γ` (Eulero-Mascheroni) — sono oggetti matematici, non parametri.
2. **Casi speciali del simplifier** derivati da identità provate: `sin(0)=0`, `exp(0)=1`, `log(1)=0` — non sono "tabella" ma assiomi.
3. **Default configurabili in CASContext** con documentazione esplicita: `default_integration_depth = 32` è accettabile se esposto come `ctx.max_integration_depth` e documentato.
4. **Limiti di sicurezza hardware**: `MAX_BIGINT_LIMBS = 10000` per evitare OOM — accettabile se causa `Unimplemented` esplicito con messaggio diagnostico.
5. **Bound CRT**: il primo primo `2^31-1` come punto di partenza del CRT multi-prime è accettabile perché viene accumulato dinamicamente fino al bound di Hadamard.

---

### Sicurezza Git (Git Safety Rule)

Mai dare agente accesso a `git reset --hard` o `git restore` without explicit safe pattern.
**NO `git reset --hard`, NO `git restore --source`**, solo backup via `git stash push`.

---

### Regola di Self-Check Prima di Ogni Commit

Prima di introdurre qualsiasi costante numerica `N` nel codice, rispondere a:

1. **"Cosa succede se l'input è 10× più grande?"** — Se la risposta è "il codice fallisce", la costante è un hardcode vietato.
2. **"Questa costante ha un nome in letteratura matematica?"** — Se no, probabilmente è arbitraria.
3. **"L'utente può cambiarla senza ricompilare?"** — Se no, deve diventare un campo di `CASContext`.
4. **"Questa costante causa un risultato matematicamente sbagliato (silenzioso) o solo un Unimplemented esplicito?"** — Il silenzio è sempre vietato; Unimplemented esplicito con diagnostica è accettabile come placeholder.

---

## REGOLE ARCHITETTURALI NON NEGOZIABILI

### 1. Aritmetica Esatta e BigInt a Limbs
- È **vietato** l'uso di `int64_t` o `double` per il calcolo simbolico.
- Il tipo `BigInt` deve usare un'implementazione a **Limbs** (`std::vector<uint32_t>` o `uint64_t`) con aritmetica bit-a-bit. 
- È **vietata** l'aritmetica basata su stringhe decimali (troppo lenta, $O(N^2)$).

### 2. Structural Sharing (Regola d'Oro delle Performance)
- L'AST è immutabile. Se una funzione (es: `simplify`, `substitute`) non modifica un nodo, **DEVE** restituire il puntatore originale (`ExprPtr`).
- È **severamente vietato** effettuare "Deep Copy" o clonazioni ricorsive dell'albero se non strettamente necessario (es. cambio di Arena).
- L'identità del puntatore è il metodo primario per verificare l'uguaglianza in $O(1)$.

### 3. Vera Memory Arena (Bump Allocator)
- I nodi AST risiedono in una `AstArena` che alloca **blocchi contigui di memoria** (es. 64KB).
- È **vietato** usare `std::make_unique` o `new` per ogni singolo nodo (distrugge la cache L1/L2).

### 4. Moltiplicazione Implicita
- Il Parser **deve** iniettare implicitamente il token `*` tra Primary Expressions adiacenti (es: `2x` -> `2*x`).

### 5. Gestione DecimalLit (Confine Simbolico/Numerico)
- I letterali decimali (`DecimalLit`) sono ammessi solo per preservare l'input utente.
- Il Core Simbolico deve restituire `CASErrorKind::Unimplemented` se un'operazione algebrica (es. espansione) richiede di operare su un `DecimalLit`. Il calcolo numerico è un fallback esplicito.

### 6. Maxima Reference Oracle — Sorgente NON Modificabile (Licenza GPL-2.0-only)

Maxima 5.49.0 (Homebrew bottle, `/opt/homebrew/Cellar/maxima/5.49.0/`) è il **reference oracle** primario della golden test suite (F0.5).

**Divieto assoluto**:
- È **vietato** modificare per qualsiasi ragione i sorgenti, i binari, i file `.lisp`, `.mac`, `.fas`, `.dem` o qualunque altro artefatto contenuto in `/opt/homebrew/Cellar/maxima/`, `/opt/homebrew/share/maxima/`, o qualunque altra installazione di Maxima sul sistema.
- È **vietato** patchare, ricompilare con flag personalizzati, o re-distribuire Maxima alterato.
- È **vietato** ridistribuire output di Maxima embedded nei nostri binari senza rispetto della GPL-2.0-only.

**Motivazione**: Maxima è rilasciato sotto **GPL-2.0-only**. Qualsiasi modifica al sorgente trasformerebbe il codice CAS Engine in un *derivative work* soggetto a copyleft, invalidando la nostra licenza proprietaria e invalidando l'uso di Maxima come oracolo indipendente nei test (oracolo modificato = non più indipendente, perde validità scientifica).

**Uso ammesso**:
- Invocazione via `maxima --very-quiet --batch-string="..."` come processo separato (fork/exec).
- Parsing dell'output testuale di Maxima per confronto AST.
- Pin esatto della versione (`5.49.0`) documentato in `PLAN_HP_PRIME_PARITY.md` e CI.
- Citazione di Maxima come reference nei doc, con link alla licenza GPL-2.0-only.

**Verifica integrità**:
- Lo script `scripts/verify_maxima_integrity.sh` deve calcolare lo SHA-256 dell'eseguibile Maxima e dei file `*.lisp` core, e fallire la build se non corrisponde al manifesto pinned `scripts/maxima_5.49.0_manifest.sha256`.
- Ogni golden run logga `maxima --version` + hash binario; mismatch → build fail.

**Conseguenze violazione**: rigetto immediato della PR + audit completo per identificare contaminazione GPL nel sorgente CAS Engine.

---

## REGOLA TIMEOUT TEST (CRITICA — Anti-Hang)

**VIETATO** lanciare l'eseguibile dei test senza un timeout esplicito. Se l'agente introduce un loop infinito (o un O(2^n) in un percorso non testato in precedenza), `cas_foundation_tests` può bloccarsi **per ore** consumando CPU senza alcun output utile.

### Regole operative obbligatorie

1. **Ogni invocazione di `cas_foundation_tests` o di qualsiasi binario di test DEVE essere lanciata con:**
   - `--gtest_filter=` esplicito sulla famiglia di test pertinente al cambio (mai eseguire l'intera suite senza filtro durante lo sviluppo iterativo);
   - timeout di shell hard via Bash tool (`timeout` parameter ≤ 120 s per filtro mirato, ≤ 600 s per la suite completa pre-commit).
2. **Mai usare `run_in_background=true` per la suite intera durante lo sviluppo iterativo**: se il processo si blocca senza output, non c'è notifica e l'agente attende silenzioso. Usare invece foreground con timeout esplicito.
3. **Suite completa solo come gate pre-commit**: dopo che tutti i test mirati passano e il codice è considerato pronto, eseguire `cas_foundation_tests` con timeout 600s. Se supera 600s, abortire e investigare hang/loop prima di riprovare.
4. **Se un test si blocca**: STOP immediato del processo (`TaskStop` o `kill`); aggiungere `std::cerr` mirato per individuare il punto di stallo; **NON ripetere** la stessa invocazione senza prima isolare la causa.
5. **In CMakeLists.txt aggiungere `set_tests_properties(... PROPERTIES TIMEOUT 60)`** per ogni test gtest individuale registrato — protezione di rete via CTest.

### Anti-pattern vietati

- `./build/cas_foundation_tests` senza filtro né timeout.
- `cmake --build && ./build/cas_foundation_tests` come "verifica finale" senza considerare possibili regressioni di complessità.
- Background bash con `tail -5` e attendere senza scadenza: se la suite hang prima del primo output, l'agente non riceve nulla.

### Risposta a un hang

Se l'agente lancia un test e non riceve output entro il timeout fissato:
1. Killare il processo (mai aspettarlo "ancora un po'").
2. Restringere il filtro al sottoset più piccolo che produce hang.
3. Aggiungere instrumentazione (`std::cerr`) nel codice modificato.
4. Identificare il loop/ricorsione; ripristinare lo stato sicuro via `git stash` se il codice introdotto è il colpevole sospetto.

---

## WORKFLOW: PRIMA DI DICHIARARE UN TASK COMPLETATO

1. **Integrità Matematica**: Tutti i test unitari e di integrazione passano al 100%.
2. **Benchmark Gate (OBBLIGATORIO)**: 
   - Esegui `bash scripts/benchmark.sh`. 
   - Se le prestazioni degradano rispetto alla baseline (`baseline_release.txt`), il codice va ottimizzato prima del merge.
3. **Sanitizers**: Zero errori in ASan e UBSan.
4. **Zero Warning**: Compilazione pulita con `-Wall -Wextra -Wpedantic -Werror`.
5. **Trace Validation**: Se il task tocca il `Simplifier`, verifica che il `ComputationTrace` sia accurato e non ridondante.

---

## PROTOCOLLO ANTI-LOOP E GESTIONE DEGLI ERRORI (Agent Panic Stop)

**VIETATO**: Entrare in loop infiniti di correzione errori (es: fix test A → rompe test B → fix test B → rompe test A) alterando ciecamente il codice o stravolgendo l'architettura nel panico da debug.

**Regola Operativa Assoluta**: Se l'agente tenta di risolvere un errore di compilazione, un test fallito o un bug e **fallisce per 3 tentativi consecutivi**, DEVE:
1. **Fermarsi immediatamente**. È severamente vietato tentare "rewrite" totali del file, aggiungere patch casuali o disabilitare test/warning per aggirare il blocco.
2. **Ripristinare lo stato di sicurezza**. Effettuare un `git stash` o scartare le modifiche non funzionanti per riportare la working tree all'ultimo stato noto stabile.
3. **Produrre un Report di Stallo**. Segnalare esplicitamente all'operatore umano: *"Ho raggiunto il limite di 3 tentativi. Ecco la natura dell'errore, le 3 strategie fallite e la potenziale causa sistemica."*
4. Attendere un intervento umano o un chiarimento architetturale. Nessun ulteriore tentativo di esecuzione codice è permesso senza autorizzazione.

---

## STANDARD TECNICI E ANTI-DEBITO

- **Linguaggio**: C++20 rigoroso.
- **Anti-Monolito**: Limite massimo di **500 righe per file sorgente**. Oltre tale soglia, lo split in moduli specializzati è **obbligatorio** (es. `simplify_arithmetic.cpp`, `simplify_functions.cpp`).
- **Error Handling**: Divieto assoluto di `throw/catch` nel core. Solo `Result<T>` monadico. Ogni path di errore deve essere tracciabile via `CASError`.
- **Gestione Memoria**: 
  - Solo `AstArena`. Vietato l'uso di `std::shared_ptr` o `std::unique_ptr` per i nodi AST.
  - **Dangling Protection**: Vietato passare `ExprPtr` tra arene diverse a meno che non sia una copia esplicita (re-interning). Nei test, usa sempre l'arena del `CASContext` per il parsing.
- **Integrità Matematica**: 
  - Ogni regola di semplificazione deve essere orientata via **LPO (Lexicographic Path Ordering)** per prevenire cicli infiniti.
  - Niente "cerotti" (if-else hardcoded per casi base). Usa il sistema di rewrite universale.
- **Testing**: Solo confronto strutturale o equivalenza matematica. **Mai** usare `toString()` per validare la logica.
- **UI/UX Desktop**: Il tool `cas_ui` è un visualizzatore professionale basato su **ImGui**. Supporta campionamento adattivo e plotting 2D/3D (Fase 10).

---

## DOCUMENTAZIONE DI RIFERIMENTO
1. Regole Architetturali: `.APROJECT_REFERENCES/01_PROJECT_GOALS/02_architectural_rules.md`
2. Roadmap (FASE 1b): `.APROJECT_REFERENCES/02_ROADMAP/01_full_roadmap.md`
3. Plotting & Sampling: `.APROJECT_REFERENCES/03_ENGINE_MODULES/11_plotting_sampling.md`

---

## PROTOCOLLO DI AUTO-EVOLUZIONE (Self-Update)

L'IA ha il mandato di proporre aggiornamenti a questo file nei seguenti momenti chiave della timeline:

### 1. Transizione di Fase (Context Shift)
All'inizio di ogni nuova Fase della Roadmap (es. passaggio dalla Fase 1b alla Fase 2), l'IA deve:
- Aggiornare i link alla documentazione di riferimento nel paragrafo sopra.
- Aggiungere eventuali "Vincoli di Dominio" specifici della nuova fase (es. regole per la derivazione simbolica nella Fase 4).

### 2. Post-Mortem di un Bug Sistemico
Se durante lo sviluppo emerge un errore strutturale o un pattern che causa regressioni di performance:
- L'IA deve aggiungere una riga nella sezione "REGOLE ARCHITETTURALI" per proibire esplicitamente quel pattern (es. "Vietato usare X perché causa Y").

### 3. Ottimizzazione del Workflow
Se uno script di automazione (es. `benchmark.sh`) viene modificato o ne viene aggiunto uno nuovo:
- L'IA deve aggiornare la sezione "WORKFLOW" per includere il nuovo comando obbligatorio.

### Vincoli di Modifica
- **MAI** rimuovere o indebolire le Regole 1, 2 e 3 (BigInt, Structural Sharing, Arena). Sono le fondamenta immutabili.
- Ogni modifica al `CLAUDE.md` deve essere comunicata allo sviluppatore con la dicitura: *"Aggiornamento Costituzione Tecnica: [MOTIVAZIONE]"*.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
