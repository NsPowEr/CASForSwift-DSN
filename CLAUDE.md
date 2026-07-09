# CAS Engine — Costituzione Tecnica (indice)

> **Legge Suprema** del progetto. Ogni agente DEVE rispettarla. Violazione → rigetto del codice (Architectural Fail).
>
> Questo file è un **indice operativo**. Il dettaglio normativo vive in `docs/rules/`
> e si carica on-demand quando pertinente — così il context resta libero a scala 4x.
> Stato corrente del progetto: **`STATE.md`** (non i vecchi `*_SESSION_*.md`, ora in `docs/archive/`).

---

## REGOLA ZERO — DIVIETO DELLA VIA FACILE

**MAI** scegliere la via più facile/veloce per comodità. Ogni decisione di design va presa scegliendo la soluzione **matematicamente corretta e generale**, anche se costa più codice/tempo. Se manca l'infrastruttura:

1. **Costruisci prima il prerequisito**, poi torna al problema (preferito); oppure
2. **Versione di passaggio segnalata** con: commento `// HARDCODE-OF-PASSAGE: <id>`, voce in `HARDCODE_LEDGER.md` (id, file, riga, descrizione, categoria, fix corretto, blocking dep), `Unimplemented` esplicito (mai silenzio).

Il commit message DEVE dichiarare ogni hardcode-of-passage. PR con hardcode non ledgered → rifiutata.
**Vietati**: lookup table chiuse quando esiste algoritmo generale; scorciatoie "passa il test attuale"; pattern matching su forma al posto di algoritmo strutturale; bail-out su tipo/range; costanti magiche al posto di parametri `CASContext`.

### 0.1 — Mandatory Specification Check (anti-allucinazione)
Prima di scrivere codice per un task di `PLAN_HP_PRIME_PARITY.md`/`CAS_TASKS.md`: `ls -lR .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/`, leggi la spec `.md` pertinente, e dichiara *"Ho letto la specifica formale [NomeFile.md] e implementerò le formule e i vincoli ivi contenuti"*. Subagent spawnato → l'orchestratore gli passa il path assoluto della spec. Codice senza spec letta = INVALIDO.

### 0.2 — Divieto assoluto di disabilitare test
**MAI** disabilitare test, commentare asserzioni, alterare benchmark per nascondere fallimenti. Test rosso = bug matematico/architetturale da risolvere nel motore. Vietati `#if 0`, macro `DISABLE`, commenti elusivi → rigetto immediato.

---

## DIVIETO HARDCODE — indice 10 categorie

Un valore hardcoded in un algoritmo matematico è un bug latente. Eccezioni: costanti matematiche esatte (π, e, φ, γ) e default configurabili via `CASContext`.
**Dettaglio completo + tabelle + esempi → [`docs/rules/hardcode-catalog.md`](docs/rules/hardcode-catalog.md)** (carica quando lavori sulla categoria).

| # | Categoria | Regola in una riga |
|---|---|---|
| 1 | Budget computazionali | Limiti depth/iter/ricorsione → configurabili via `CASContext` + `Unimplemented` diagnostico |
| 2 | Costanti magiche algebriche | Derivare da bound (Mignotte/Hadamard/Schwartz-Zippel) o parametro di confidenza |
| 3 | Set/range fissi | Set chiuso → esaustivo, oppure generatore algoritmico. Mai silenzio fuori range |
| 4 | Bail-out su tipo | Rifiuto per tipo → esplicito, diagnostico, temporaneo (task aperto). Mai silenzioso |
| 5 | Ordinamenti/strutture fisse | Scelte strutturali (MonomialOrder, pivot) parametrizzabili |
| 6 | Seed/randomness | Seed da hash dell'input o esplicito in `CASContext`; mai seed fisso |
| 7 | Nomi variabili interni | Sempre `ctx.make_fresh_symbol(prefix)`; mai nomi letterali (`__mrv_w`, `C1`) |
| 8 | Pattern matching a tabella | Solo come fast-path; deve esistere sempre un path algoritmico generale |
| 9 | Intervalli polling fissi | Check timeout adattivo al costo reale, configurabile |
| 10 | Rank di crescita statici | Rank dinamico via confronto asintotico (Gruntz) + coefficiente leader |

**Self-check prima di ogni costante `N`**: (1) input 10× più grande → fallisce? = hardcode. (2) Ha un nome in letteratura? No → arbitraria. (3) Cambiabile senza ricompilare? No → campo `CASContext`. (4) Causa risultato sbagliato silenzioso? Vietato; `Unimplemented` diagnostico ok.

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
Il Parser **deve** iniettare implicitamente `*` tra Primary Expressions adiacenti (`2x` → `2*x`).

### 5. Gestione DecimalLit (confine simbolico/numerico)
`DecimalLit` ammessi solo per preservare l'input utente. Il core simbolico ritorna `CASErrorKind::Unimplemented` se un'operazione algebrica richiede di operare su `DecimalLit`. Il calcolo numerico è fallback esplicito.

### 6. Maxima Reference Oracle — GPL-2.0-only, NON modificabile
Maxima 5.49.0 (`/opt/homebrew/Cellar/maxima/5.49.0/`) è l'oracolo della golden suite. **VIETATO** modificare/patchare/ricompilare qualsiasi sorgente/binario/`.lisp`/`.mac`/`.fas` Maxima, o embeddarne output nei binari. Uso ammesso: solo fork/exec (`maxima --very-quiet --batch-string=...`) + parsing testuale. Integrità via `scripts/verify_maxima_integrity.sh` (SHA-256 pinned). Motivo: modifica → derivative work copyleft + oracolo non più indipendente. **Dettaglio → [`docs/rules/maxima-oracle.md`](docs/rules/maxima-oracle.md)**.

---

## REGOLA TIMEOUT TEST (CRITICA — anti-hang)

**VIETATO** lanciare i test senza timeout esplicito: un loop introdotto può bloccare `cas_foundation_tests` per ore senza output.
- **Modo preferito**: `bash scripts/test_quick.sh` (cap 1200s — misura reale ~803s @2565 test ASan/M1 Pro 2026-07-09 — esclude slow noti) per sviluppo; `--slow` (cap 1800s) come gate pre-commit. Non modificare la lista esclusioni senza indagare la regressione.
- Invocazione diretta sempre con `--gtest_filter=` mirato + timeout shell (≤120s mirato, ≤1200s quick, ≤1800s slow).
- **Mai** `run_in_background` per la suite intera in sviluppo (hang silenzioso). Foreground + timeout.
- Test bloccato → STOP immediato (kill), restringi filtro, instrumenta con `std::cerr`, isola la causa. **NON** ripetere la stessa invocazione alla cieca.
- CMakeLists: `set_tests_properties(... PROPERTIES TIMEOUT 60)` per ogni gtest registrato.

---

## PROTOCOLLO ANTI-LOOP (Agent Panic Stop)

**VIETATO** loop di fix ciechi (fix A rompe B, fix B rompe A). Dopo **3 tentativi consecutivi falliti** su un errore l'agente DEVE:
1. Fermarsi (no rewrite totali, no patch casuali, no disabilitare test/warning).
2. Ripristinare stato sicuro (`git stash` — **mai** `git reset --hard` / `git restore --source`).
3. Produrre Report di Stallo: natura dell'errore, 3 strategie fallite, causa sistemica sospetta.
4. Attendere intervento umano. Nessuna ulteriore esecuzione senza autorizzazione.

**Git safety**: solo `git stash push` per backup. NO `git reset --hard`, NO `git restore --source`.

---

## WORKFLOW — prima di dichiarare un task completato

1. **Integrità matematica**: tutti i test unit + integrazione passano al 100%.
2. **Benchmark gate** (OBBLIGATORIO): `bash scripts/benchmark.sh`; regressione vs `baseline_release.txt` → ottimizza prima del merge.
3. **Sanitizers**: zero errori ASan + UBSan.
4. **Zero warning**: `-Wall -Wextra -Wpedantic -Werror` pulito.
5. **Trace validation**: se tocchi il `Simplifier`, verifica `ComputationTrace` accurato e non ridondante.

---

## STANDARD TECNICI E ANTI-DEBITO

- **Linguaggio**: C++20 rigoroso.
- **Anti-Monolito**: max **500 righe/file** (soft), hard block a 550 (hook `guard_file_size.sh`). Oltre → split obbligatorio o whitelist con ticket.
- **Error Handling**: niente `throw/catch` nel core. Solo `Result<T>` monadico, ogni errore tracciabile via `CASError`.
- **Memoria**: solo `AstArena`. Vietati `shared_ptr`/`unique_ptr` per nodi AST. Vietato passare `ExprPtr` tra arene diverse senza re-interning; nei test usa l'arena del `CASContext`.
- **Integrità**: ogni regola di rewrite orientata via **LPO** (no cicli). Niente "cerotti" if-else hardcoded: usa il rewrite universale.
- **Testing**: solo confronto strutturale o equivalenza matematica. **Mai** `toString()` per validare la logica.
- **UI**: `cas_ui` è visualizzatore ImGui (campionamento adattivo, plot 2D/3D, Fase 10).

---

## NAVIGAZIONE & TOOLING (scala 4x)

- **Codice**: `graphify query "<domanda>"` PRIMA di grep raw (subgraph scoped, molto più piccolo). `graphify path/explain` per relazioni/concetti. `graphify update .` dopo modifiche (auto via hook `graphify_autoupdate.sh`).
- **Ledger/task**: NON aprire i markdown da 100-180 KB per una voce. Usa `python3 scripts/ledger_index.py {hc <id> | task <id> | search <testo> | open | stats}`.
- **Stato**: `STATE.md`. **Storico**: `docs/archive/`. **Regole dettaglio**: `docs/rules/`.

---

## DOCUMENTAZIONE DI RIFERIMENTO
1. Regole architetturali: `.APROJECT_REFERENCES/01_PROJECT_GOALS/02_architectural_rules.md`
2. Roadmap: `.APROJECT_REFERENCES/02_ROADMAP/01_full_roadmap.md`
3. Plotting & Sampling: `.APROJECT_REFERENCES/03_ENGINE_MODULES/11_plotting_sampling.md`
4. Catalogo hardcode (dettaglio): `docs/rules/hardcode-catalog.md`

---

## AUTO-EVOLUZIONE (Self-Update)
L'IA propone aggiornamenti a questo file: (1) a ogni transizione di Fase (aggiorna link + vincoli di dominio), (2) post-mortem di bug sistemico (aggiungi divieto esplicito del pattern), (3) nuovo script di automazione (aggiorna WORKFLOW).
**Vincoli**: MAI rimuovere o indebolire le Regole 1/2/3 (BigInt, Structural Sharing, Arena) — fondamenta immutabili. Ogni modifica a `CLAUDE.md` comunicata con *"Aggiornamento Costituzione Tecnica: [MOTIVAZIONE]"*.
