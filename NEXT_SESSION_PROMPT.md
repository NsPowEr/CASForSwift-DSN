# Sessione CAS Engine — Continuazione Sviluppo Post-F4

## Stato attuale (consolidato)

**Fasi chiuse**: F0 → F3 + F4 (linalg core: Bareiss, LU, Cholesky LDL^T,
QR Householder, Hermite, Smith Z+Q[x], Jordan multi-blocchi Filippov,
Companion, Vandermonde, Tridiagonal, Circulant n≤4).

**Suite**: 1914+ test PASS. Stress F4 pass tranne 2 noti:
  - `F4StressTest.Householder_QR_8x8_RandomQ` → timeout (QR-symbolic-simplify)
  - `StressTest.MatrixDeterminant100x100GiantInteger` → `stod out_of_range`

**HARDCODE_LEDGER aperti** (zero crescita, chiudere prima di aprirne nuovi):
  - HC-F43-CIRCULANT-GT4  (formula chiusa via Q(ω_n))
  - HC-F43-TOEPLITZ        (Trench/Levinson simbolico)
  - HC-F43-BANDED          (LU banded O(n·k²))
  - + voci pre-F4 elencate in `HARDCODE_LEDGER.md`

## Obiettivi prossima sessione

1. **Chiudere DUE debt residui pre-esistenti**:
   a. **QR-symbolic-simplify**: 8x8 random Q timeout. Audit:
      - simplify ricorsivo su sqrt(Σ x_i²) genera AST esponenziale?
      - `together` mancante post-riflessione?
      - alternativa: Gram-Schmidt modificato con normalizzazione lazy +
        memoizzazione di Σ x_i² per colonna.
      Soluzione corretta: bound matematico sul numero di simplify per
      colonna (≤ k+1 per Householder), no fallback rimosso.

   b. **stod out_of_range 100x100**: parser numerico chiama `std::stod`
      su literali troppo grandi. Sostituire con `BigInt` o `Decimal`
      arbitrary-precision per tutti i numeric literals di matrici stress.

2. **Avviare F5 (Calculus avanzato)** — Risch + LRT completo:
   - Implementare differential field algebraicamente, non pattern-matching
     (vietata Cat. 8 CLAUDE.md).
   - Hermite-Ostrogradsky reduction su Q(x).
   - Lazard-Rioboo-Trager (LRT) per parte logaritmica completa.
   - Risch principal step per trascendentali elementari (exp, log).
   - Algoritmo Bronstein (cap. 5) come riferimento — niente scorciatoie.

## Regole operative non negoziabili

### REGOLA ZERO — divieto via facile
Riferimento: `CLAUDE.md`. Ogni decisione di design deve scegliere la
soluzione **matematicamente corretta e generale**, anche se costa più
tempo e codice. Se prerequisito infrastrutturale manca:
  - opzione A (preferita): costruire il prerequisito prima;
  - opzione B: hardcode-of-passage con commento `HARDCODE-OF-PASSAGE:`
    + entry in `HARDCODE_LEDGER.md` + diagnostic `Unimplemented`
    esplicito per casi fuori scope.

### 10 categorie hardcode vietati
Vietato:
  1. budget computazionali non configurabili (uso `CASContext` getter/setter);
  2. costanti magiche in algoritmi algebrici (bound Mignotte/Hadamard/Schwartz-Zippel);
  3. set chiusi che escludono input validi (generatori algoritmici);
  4. bail-out su tipo (estendere dispatcher);
  5. ordinamenti/strutture fisse (parametrizzare via enum);
  6. seed e randomness deterministica non hash-derivata;
  7. nomi interni hardcoded (usare `ctx.make_fresh_symbol`);
  8. pattern matching a tabella chiusa (impl algoritmica generale);
  9. intervalli polling fissi (configurabili + adattivi);
  10. gerarchie crescita statiche (Gruntz dinamico).

### Workflow per ogni task
1. Leggere `CLAUDE.md` (Legge Suprema), `HARDCODE_LEDGER.md`, `MEMORY.md`.
2. Lettura targeted di file/grep — NIENTE Explore agent se filename noto.
3. Progettare la soluzione matematica corretta prima di scrivere codice.
4. Se proposta soluzione "facile" affiora, fermarsi e ragionare se esiste
   path generale; spiegare scelta nel codice/commit.
5. TDD: scrivere prima i test cert simbolici (non solo numeric razionali);
   includere edge case RootOf, simboli con assumptions positive/negative,
   matrici simboliche miste.
6. Implementare la soluzione corretta.
7. Build (`cmake --build build-release` con Ninja) + gtest filter mirato.
8. AcidTest + SupremeTest (subagent `cas-regression-guard`).
9. Hardcode audit (skill `hardcode-audit`) sui file toccati.
10. Benchmark gate (`bash scripts/benchmark.sh` vs baseline).
11. Commit con messaggio che dichiara ogni HC introdotto/chiuso.

### Test simbolico obbligatorio
Per ogni algoritmo linalg/algebra/calculus:
  - test su input numerico razionale (sanity);
  - test su input con almeno 2 simboli + assumptions positive/negative;
  - test su input RootOf (quando supportato);
  - test certificato (Q·R == A, A·A⁻¹ == I, ∂(∫f) == f, ecc.);
  - test edge case (matrice singolare, polinomio costante, ecc.).

### Oracolo Maxima
Riferimento `Maxima 5.49.0` (GPL-2.0-only, NON modificabile). Per ogni
nuovo algoritmo del core mat:
  - golden test che invoca Maxima via `maxima --very-quiet --batch-string`;
  - confronto AST normalizzato (non `toString()`);
  - verifica hash binario via `scripts/verify_maxima_integrity.sh`.

### Anti-loop debug
Se 3 tentativi consecutivi di fix falliscono:
  - STOP;
  - `git stash` per ripristinare working tree;
  - report stallo all'operatore: natura errore, 3 strategie fallite,
    ipotesi causa sistemica;
  - attendere intervento umano.

### Standard tecnici
  - C++20 rigoroso; `-Wall -Wextra -Wpedantic -Werror`;
  - solo `BigInt`/`Rational` nel core simbolico (vietati `int64_t`,
    `double`);
  - solo `Result<T>` monadico (vietato `throw/catch`);
  - solo `AstArena` (vietati `shared_ptr`/`unique_ptr` per nodi AST);
  - structural sharing: se non modifica, ritorna `ExprPtr` originale;
  - file <500 LOC, split obbligatorio oltre soglia;
  - LPO orientation per ogni rewrite rule;
  - zero `throw/catch` nel core;
  - zero warning su -Wall.

### Git safety
  - NO `git reset --hard`, NO `git restore --source`;
  - solo `git stash push` per backup;
  - mai bypass hook (`--no-verify`) senza richiesta esplicita;
  - mai amend; sempre nuovo commit;
  - mai force push.

## Priorità ordinata

P0 (blocker prima di F5):
  1. QR-symbolic-simplify (bound matematico, no fallback)
  2. stod out_of_range parser (sostituzione con arbitrary-precision)

P1 (apertura F5 — Risch foundation):
  3. Differential field algebraico (Bronstein cap. 3)
  4. Hermite-Ostrogradsky reduction
  5. LRT (Lazard-Rioboo-Trager) completo
  6. Risch principal step exp/log

P2 (chiusura debt F4 strutturali):
  7. HC-F43-BANDED → LU banded O(n·k²)
  8. HC-F43-TOEPLITZ → Trench/Levinson simbolico
  9. HC-F43-CIRCULANT-GT4 → aritmetica Q(ω_n)

## Cosa NON fare
  - Non aprire nuovi HC senza chiuderne uno esistente.
  - Non implementare via lookup table o pattern matching chiuso.
  - Non bypass `is_known_nonzero` con check strutturale.
  - Non aggiungere `int64_t`/`double` per "comodità".
  - Non scrivere file >500 LOC senza split.
  - Non claim "production-ready" senza test cert simbolico completo.
  - Non sopprimere warning, never `-Wno-*` ad-hoc.
  - Non dichiarare task completo se anche un solo test fallisce
    (anche pre-esistente — va ledgered esplicitamente).

## Riferimenti documentali obbligatori (leggere a inizio sessione)
  - `CLAUDE.md` (Legge Suprema)
  - `HARDCODE_LEDGER.md` (debt registry)
  - `.APROJECT_REFERENCES/01_PROJECT_GOALS/02_architectural_rules.md`
  - `.APROJECT_REFERENCES/02_ROADMAP/01_full_roadmap.md`
  - `memory/MEMORY.md` (memoria auto-gestita)

## Filosofia operativa
> Cerchiamo la perfezione, non scrittura di codice tanto per scrivere.
> La via più facile è quasi sempre la via sbagliata. Ogni scorciatoia
> oggi è un debito domani che blocca un input legittimo del CAS.
> Tempo investito in soluzione corretta = anni di non-regressioni.

Ogni decisione architetturale deve essere giustificata con riferimento
a un teorema, un bound matematico, o una scelta di design consapevole
documentata. Il CAS è un sistema matematico, non un toolkit di euristiche.
