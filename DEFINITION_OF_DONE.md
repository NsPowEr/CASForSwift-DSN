# Definition of Done (DoD) — CAS Engine

> Documento canonico. Consolida i criteri "Risolta" sparsi in `CAS_TASKS.md`,
> `PLAN_HP_PRIME_PARITY.md` e la CHECKLIST ANTI-FURBIZIA. Ogni task di
> implementazione DEVE soddisfare questa DoD prima di essere marcato `DONE`.
> Subordinata a `CLAUDE.md` (Legge Suprema): in caso di conflitto vince CLAUDE.md.

---

## Vocabolario di stato (unico, machine-parsable)

`scripts/tasks_audit.sh` normalizza ogni ledger a questi 5 stati. Significato
**vincolante**:

| Stato | Significato | Quando usarlo |
|---|---|---|
| **DONE** (`Risolta`) | Dominio dichiarato chiuso al 100% + tutta la DoD soddisfatta + certificatore indipendente verde | Solo quando OGNI gate sotto passa |
| **PARTIAL** (`Parziale`/`Parziale avanzata`) | Algoritmo vero ma dominio è un **subset esplicito dichiarato nel nome del task** | Subset reale, mai per "MVP" o "funziona sui test presenti" |
| **OPEN** (`Aperta`) | Placeholder, stub, o non implementato | Niente codice reale, oppure `Unimplemented` diagnostico |
| **PENDING** | In coda, non iniziato | Backlog |
| **UNKNOWN** | Stato non riconosciuto dal parser | Da bonificare: riscrivere lo stato |

> **Vietato il termine "MVP".** Un task è `PARTIAL` con subset esplicito, oppure `DONE`.

---

## I 10 gate della DoD

Un task passa a **DONE** solo se TUTTI i gate applicabili sono verdi.
I gate marcati 🤖 sono verificabili da `scripts/check_dod.sh`.

### G1 — Spec letta (REGOLA 0.1) 🤖
Se esiste `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<task>.md`, la spec è stata
letta e il commit/PR dichiara: *"Ho letto la specifica formale [NomeFile.md]"*.

### G2 — Algoritmo generale, non scorciatoia (REGOLA ZERO)
Nessuna lookup table a set chiuso dove esiste algoritmo generale. Nessun
pattern-matching su forma come unico path. Ogni fast-path ha un fallback algoritmico.

### G3 — Zero hardcode vietato (Categorie 1-10) 🤖
Nessuna costante numerica priva di giustificazione matematica formale. Ogni budget
computazionale è configurabile via `CASContext`. Eventuali hardcode-of-passage sono
marcati `// HARDCODE-OF-PASSAGE:` **e** iscritti in `HARDCODE_LEDGER.md`.
Verifica: `Skill hardcode-audit` o agent `hardcode-auditor`.

### G4 — Aritmetica esatta (REGOLA 1) 🤖
Solo `BigInt` e `Rational` nel core simbolico. Nessun `int64_t`/`double` per calcolo
simbolico (eccetto numeric engine esplicito).

### G5 — Test robusti + anti-hardcode
- Unit test con **variabili diverse** (non solo `x`).
- Test su **forme sintattiche equivalenti**.
- Almeno un **test anti-hardcode** (input 10× più grande del caso base).
- Validazione per **equivalenza matematica o struttura**, mai via `toString()`.

### G6 — Certificatore indipendente (golden oracle)
Per capacità matematiche: golden diff vs **Maxima 5.49.0** (e/o SymPy) verde sul
dominio dichiarato. Agent `maxima-golden-diff`. Maxima fork/exec only (CLAUDE.md §6).

### G7 — Zero regressioni 🤖
`bash scripts/test_quick.sh` verde (≤600s). Nessun test prima verde ora rosso.
Agent `cas-regression-guard` (AcidTest + SupremeTest).

### G8 — Benchmark gate 🤖
`bash scripts/benchmark.sh` non degrada vs `baseline_release.txt`. Se degrada,
ottimizzare prima del merge (REGOLA 0.2: vietato alterare il benchmark).

### G9 — Build pulito 🤖
Compilazione `-Wall -Wextra -Wpedantic -Werror` zero warning. Zero errori ASan/UBSan.
Anti-monolith: nessun file >500 LOC fuori whitelist (`scripts/check_file_size.sh`).

### G10 — Tracker aggiornato + limiti dichiarati 🤖
Stato aggiornato nel ledger di competenza. Limiti residui dichiarati esplicitamente.
Nessun sorgente orfano (`scripts/check_orphan_sources.sh`). Nessun `Unimplemented`
silenzioso: ogni bail-out riporta input shape + modulo + motivo + ticket ID.

---

## Checklist pre-commit (copia in ogni PR)

```
DoD — Definition of Done
[ ] G1  Spec letta (se esiste in MISSING_FEATURES_SPECS/)
[ ] G2  Algoritmo generale, nessuna scorciatoia (REGOLA ZERO)
[ ] G3  Zero hardcode vietato; eventuali passaggi in HARDCODE_LEDGER.md
[ ] G4  Solo BigInt/Rational nel core
[ ] G5  Test: variabili diverse + forme equivalenti + anti-hardcode
[ ] G6  Golden diff vs Maxima verde sul dominio dichiarato
[ ] G7  scripts/test_quick.sh verde, zero regressioni
[ ] G8  scripts/benchmark.sh non degrada vs baseline
[ ] G9  -Werror zero warning, ASan/UBSan clean, anti-monolith ok
[ ] G10 Tracker aggiornato, limiti dichiarati, zero orfani, zero Unimplemented muto
```

Stato finale dichiarato: **DONE** solo con tutti i box ✓ sul **dominio dichiarato**.
Altrimenti **PARTIAL** col subset esplicito nel nome, o **OPEN**.

---

## Applicazione retroattiva ai task pending

Tutti i task `PENDING`/`PARTIAL`/`OPEN` nei ledger sono soggetti a questa DoD al
momento della lavorazione. La promozione a `DONE` di un task storico richiede il
re-pass completo dei 10 gate, non la sola esistenza di codice.

`scripts/check_dod.sh` automatizza G1/G3/G4/G7/G8/G9/G10 (gate 🤖); G2/G5/G6
restano verifiche di giudizio umano/agent documentate nel commit.
