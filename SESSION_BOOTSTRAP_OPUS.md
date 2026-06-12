# Session Bootstrap — Opus autonomy CAS Engine F8

Documento di ingresso obbligatorio. Leggere PRIMA di ogni azione.

---

## 0. Mission

Chiusura completa task pending elencate in `PLAN_TASKS_REMAINING.md` rispettando vincoli `CLAUDE.md` REGOLA ZERO. Output finale atteso: 0 task pending, 0 nuove voci HARDCODE_LEDGER, suite quick verde, ledger riconciliato.

---

## 1. Pre-flight checklist (eseguire SEMPRE all'avvio sessione)

```bash
# 1.1 Stato repo
git status -s
git log --oneline -10

# 1.2 Plan + ledger
cat PLAN_TASKS_REMAINING.md     # roadmap task residue
cat HARDCODE_LEDGER.md | head -60   # debiti aperti

# 1.3 Suite baseline
bash scripts/test_quick.sh > /tmp/baseline.log 2>&1
grep -E "PASSED|FAILED" /tmp/baseline.log | tail -5
```

**Stato baseline atteso (commit `a5c9ee9` o successivi)**:
- 2346+ test PASS
- 1 FAIL pre-esistente: `F2GateBenchmark.FactorOneHundredRandomZxUnderBudget` — NON regressione, baseline-confirmed
- 2 SKIPPED noti

Se baseline diverge: investigare PRIMA di toccare nuovo codice.

---

## 2. Task selection algorithm

**Sorgente autoritativa**: `PLAN_NEXT_SESSIONS.md` §"Sessione N (CURRENT)".

Procedura:
1. Apri `PLAN_NEXT_SESSIONS.md`, identifica la sessione marcata `CURRENT`.
2. Esegui gli step in ordine (Step 1.1, 1.2, 1.3, ...).
3. Ogni step ha checklist `- [ ]` → marca `- [x]` al completamento.
4. Al termine sessione: sposta blocco intero sotto §"Log sessioni" (in alto), promuovi la sessione successiva a `CURRENT`.

**Fallback**: se `PLAN_NEXT_SESSIONS.md` non esiste o sessione corrente è completata senza nuova promozione, consultare `PLAN_TASKS_REMAINING.md` §ordering come backup.

```
For each step in PLAN_NEXT_SESSIONS.md §"Sessione CURRENT":
  if step.status == [x]: skip
  if step.blockedBy contains non-completed step: skip
  if step.effort > budget rimasto: defer + nuovo sub-step in PLAN
  else: select step → proceed
```

**Budget sessione = 4-6 ore di lavoro produttivo**. Se task richiede > 1 giorno-uomo, suddividere in sub-task per sessione corrente. Non iniziare ciò che non si finisce nella sessione.

---

## 3. Working style (per ogni task)

### 3.1 Sequenza obbligatoria

```
1. Lettura spec   → .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<TaskName>.md
                    + dichiarazione esplicita nel primo messaggio:
                    "Ho letto la specifica formale [NomeFile.md] e implementerò
                     le formule esatte e i vincoli ivi contenuti"
2. Audit corrente → grep esistente impl, leggere file rilevanti
3. Plan ctx params → quali nuovi ctx params occorrono (REGOLA ZERO)
4. Implement      → file nuovi/modificati, citation in commento header
5. Test mirato    → --gtest_filter='<TestSuite>.*' + timeout 60-180s
6. Quick suite    → bash scripts/test_quick.sh (capo timeout 700s)
7. Ledger update  → aggiungere/chiudere voci HPP/HC corrispondenti
8. Graph update   → graphify update . (background)
9. Commit atomic  → conventional commit, HEREDOC per body multi-riga
10. Mark done     → /Users/davidesaba/.claude/tasks/<dir>/<id>.json status=completed
```

### 3.2 Anti-pattern (rifiutare immediatamente)

| Tentazione | Azione corretta |
|------------|-----------------|
| "Aggiungo `int64_t` per velocità" | Solo `BigInt`/`Rational` (REGOLA 1) |
| "Hardcodo 16 come soglia" | Ctx param + citazione matematica |
| "Lascio Unimplemented silenzioso" | Diagnostico esplicito con error::reason_codes |
| "Disabilito test che fallisce" | Fix bug nel motore (REGOLA 0.2) |
| "Provo 4 fix successivi al volo" | Stop a 3 tentativi → report umano (protocollo anti-loop) |
| "Touch 28 file in un commit" | 1 file/commit, test quick dopo ogni split (Task 25) |
| "Cito 'Optimal'/'Faster'/'Best'" | Cita paper specifico + sezione |
| "git reset --hard per pulire" | Solo `git stash push` |
| "Bench con 1 sample" | Bench multi-size + tabella regression |

### 3.3 Decisione defer-vs-push

```
if (spec esiste AND infra base presente AND tests scritti AND < 1 giorno):
    push (implement + test + commit)
elif (spec esiste AND ≥ 80% gia` fatto):
    verify pre-existing impl + reconcile ledger + mark done
elif (spec esiste AND infra base assente):
    defer: aggiungere sub-task ≤1day in PLAN_TASKS_REMAINING.md
elif (spec mancante):
    spec write first (Knuth/Bronstein/Cohen/Stauduhar reference)
elif (effort > sessione):
    suddividere in sub-task atomiche, eseguire primo, defer resto
```

---

## 4. Reasoning patterns

### 4.1 Per algoritmi matematici

Domande da rispondere PRIMA di codice:
1. Qual è il **bound matematico** della soglia/iterazione? (Mignotte, Hadamard, Schwartz-Zippel)
2. Qual è il **caso peggiore** (worst-case complexity)?
3. Cosa succede se input è **10× più grande**? (test scalabilità)
4. Il risultato è **strutturalmente esatto** o solo numericamente?
5. Esiste un **oracle indipendente** (Maxima, SymPy) per cross-check?

### 4.2 Per scelte architetturali

```
Domanda: "X o Y?"
  Risposta = (more_general(X, Y), citation_strength(X, Y),
              ctx_configurability(X, Y), AST_immutability_preserved(X, Y))
```

NON pesare velocità di codifica. Pesare correttezza matematica + estensibilità.

### 4.3 Per debugging

Protocollo:
1. Isola test che riproduce bug a granularità minima.
2. Aggiungi `std::cerr` mirato (non printf su tutto).
3. Verifica stato AST con `debug_print(ExprPtr)`.
4. Cross-check con Maxima se identità algebrica.
5. 3 tentativi max → stop + report (REGOLA ANTI-LOOP).

---

## 5. Skills/strumenti disponibili (ordine di preferenza)

### 5.1 Per ricerca codice/spec
- **Read** path absoluto (precision)
- **Bash grep -n** con pattern stringa (NON regex con caratteri speciali, ripgrep escape buggy)
- **graphify query "<question>"** se grafo aggiornato

### 5.2 Per oracle indipendente
- **Maxima 5.49.0** (`/opt/homebrew/Cellar/maxima/5.49.0/`) — fork/exec only
- Mai modificare sorgenti Maxima (GPL-2.0-only, REGOLA 6)

### 5.3 Per test scaling
- Quick: `bash scripts/test_quick.sh` (≤600s)
- Slow gate: `bash scripts/test_quick.sh --slow` (≤1800s)
- Mirato: `timeout 60 ./build/cas_foundation_tests --gtest_filter='Suite.*'`

### 5.4 Per benchmark
- `bash scripts/benchmark.sh` vs baseline `baseline_release.txt`

### 5.5 Subagent (uso parsimonioso)
- `hardcode-auditor` — scansione anti-pattern (NO modifica codice)
- `maxima-golden-diff` — cross-check espressioni vs Maxima
- `cas-regression-guard` — AcidTest + SupremeTest post-modifica core
- `Explore` — broad code search (NO open-ended analysis)
- `general-purpose` — multi-step research quando confidente di trovare risposta lontano dal cursore

Avoid: dispatching subagent per task >50 token di lavoro recuperabile localmente.

---

## 6. Commit policy

### 6.1 Conventional format

```
<type>(<scope>): <subject-line-≤72-chars>

<body multi-paragraph: cosa, perché, dove, test, refs>

References:
- Paper author year title section
- spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<File>.md
```

Types ammessi: `feat`, `fix`, `refactor`, `docs`, `test`, `perf`, `chore`, `ci`.

### 6.2 Atomicità

1 commit = 1 task chiuso O 1 sub-task chiuso. Mai mischiare task indipendenti.

### 6.3 Body obbligatorio per

- Algoritmi matematici (citation + complexity + correctness argument)
- Closure di voci ledger
- Modifiche header pubblici
- Splits di file monolith (motivazione del cut)

### 6.4 NO `--no-verify`, NO `--amend` (sempre new commit), NO `git push --force`

---

## 7. Quick reference — file critici

```
.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/   # spec formali (read first)
PLAN_NEXT_SESSIONS.md                          # **AUTORITATIVO** — sessione corrente
PLAN_TASKS_REMAINING.md                        # roadmap pending (legacy, fallback)
TODO_PH8.md                                    # Phase 8 tracker (auto-update)
HARDCODE_LEDGER.md                             # debiti tracciati
CLAUDE.md                                      # legge suprema
SESSION_BOOTSTRAP_OPUS.md                      # questo file
include/cas/cas_context_params.hpp             # ctx params
include/cas/symbolic.hpp                       # CASContext API
include/cas/bigint.hpp                         # BigInt API
src/symbolic/simplify_impl.hpp                 # Simplifier class
scripts/test_quick.sh                          # gate quick (≤600s)
```

---

## 8. Stato commit attuale (sessione 2026-06-12)

```
a5c9ee9 docs(F4.1 / Task 12): close HPP-F4.1-QR-HOUSEHOLDER ledger entry
6b6671a docs(F8): comprehensive plan for remaining tasks + ledger entries
3ff0840 feat(F4.K partial): sqrt(x^2) branch-cut gating
5268e90 feat(F6.D): adaptive Gauss-Kronrod G7/K15 with priority-queue panels
70a1e52 feat(F1.2 / HPP-023): Burnikel-Ziegler recursive division
ade20ff feat(F7.5.E2 / Bessel): half-integer recurrence + I/K integer expansion
4131398 feat(F1.3 / HPP-019): double-digit Lehmer GCD (Knuth L + Jebelean)
```

Task done: 1, 2, 3, 5, 6, 11, 12, 21, 24 (9 totali).
Task partial: 20.
Task pending: 4, 7, 9, 10, 17, 22, 25, 26.

---

## 9. Failure recovery

Se sessione precedente ha lasciato:
- Working tree sporca: `git status` → decidere `git stash push` (NON `restore --source`) o continuare.
- Test rosso non pre-existing: investigare PRIMA di proseguire.
- Build rotto: NON pushare. Fix prima.
- Ledger inconsistente con codice: riconciliare ledger (come fatto per Task 12).

---

## 10. Termination

Sessione termina quando:
1. Task target chiuso + test + commit + ledger update + task json updated, OPPURE
2. Task target richiede > budget rimasto → defer in plan + stop, OPPURE
3. 3 fallimenti consecutivi (REGOLA ANTI-LOOP) → stop + report umano.

Report finale obbligatorio:
- Task chiuse (con commit hash)
- Task partial (con scope completato/residuo)
- Task defer (con motivazione)
- Suite quick: pass/fail count vs baseline
- Next session pickup point

---

## 11. Caveman mode

Tutti i message user-facing in caveman ultra (drop articoli/filler/pleasantries; fragments OK; abbreviazioni: DB/auth/config/req/res/fn/impl; arrows X → Y per causalità).

**Eccezioni** (write normal): commit messages, PR descriptions, code comments, security warnings, irreversible action confirmations, multi-step sequences ad alto rischio misread.

Persistente per tutta sessione, anche dopo molti turni. Off solo su `stop caveman` esplicito.
