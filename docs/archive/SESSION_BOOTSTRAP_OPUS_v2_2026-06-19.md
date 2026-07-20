# SESSION_BOOTSTRAP_OPUS — entry-point sessione autonoma

> Punto di ingresso per una sessione Opus autonoma sul CAS Engine. Leggi questo
> file per primo, **non** caricare i markdown da 100-180 KB se non serve.
> A scala 4x il context è una risorsa scarsa: naviga con gli indici, non con i dump.

## §1 — Carica solo questo, in quest'ordine

1. `CLAUDE.md` — Costituzione (indice). Il dettaglio normativo è in `docs/rules/`,
   caricalo on-demand quando lavori sulla categoria pertinente.
2. `STATE.md` — stato vivo. Fasi chiuse, debiti aperti, corpus golden, prossimi passi.
3. `PLAN_NEXT_SESSIONS.md` — la sessione marcata `(CURRENT)`. Esegui quel blocco.
4. La spec formale del task: `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/<task>.md`
   (REGOLA 0.1 — obbligatoria PRIMA di scrivere codice).

Indici (preferiscili all'apertura dei markdown grandi):
- Task / hardcode: `python3 scripts/ledger_index.py {task <id> | hc <id> | search <txt> | stats}`
- Codice: `graphify query "<domanda>"` (NON grep raw). `graphify path/explain` per relazioni.

## §2 — Vincoli invarianti per ogni sessione (gate di chiusura task)

1. **Spec letta** prima del codice (REGOLA 0.1) — dichiara la frase di conferma.
2. **No via facile** (REGOLA ZERO): ogni soglia/budget → parametro `CASContext`,
   mai costante magica. Hardcode di passaggio solo se ledgered + commento + `Unimplemented`.
3. **Test mirato**: `--gtest_filter` reale (mai `*`) + `timeout 60-180s`. Mai la suite
   nuda. Per la suite usa `bash scripts/test_quick.sh` (cap 600s) — `--slow` come gate.
4. **No disabilitazione test** (REGOLA 0.2). Test rosso = bug nel motore, non nel test.
   Un `DISABLED_` richiede `TEST-DISABLED-JUSTIFY:` nel commit (hook commit-msg).
5. **Anti-monolito**: file ≤500 LOC (hard block 550). Split o whitelist con ticket.
6. **Debt-gate verde**: `bash scripts/debt_gate.sh --staged` prima del commit.
7. **Definition of Done**: `bash scripts/check_dod.sh [--full]` — vedi `DEFINITION_OF_DONE.md`.
8. **Ledger atomico**: ogni task chiuso → update di `HARDCODE_LEDGER.md` / tracker nello
   stesso commit. 1 task = 1 commit, conventional format, body con citazione spec.

## §3 — Protocollo anti-loop (NON negoziabile)

Dopo **3 fallimenti consecutivi** sullo stesso errore: **STOP**. Non riscrivere tutto,
non patchare a caso, non disabilitare test/warning. Esegui: (1) fermati, (2) `git stash push`
(MAI `reset --hard` / `restore --source`), (3) Report di Stallo (natura errore + 3 strategie
fallite + causa sistemica sospetta), (4) attendi intervento umano.
È anche meccanizzato: l'hook `guard_anti_loop.sh` avvisa al 3° comando build/test identico
e blocca al 6°. Se vieni bloccato, **non** aggirare — è il segnale di fare lo stallo.

## §4 — Stato degli strumenti di verifica (LEGGI prima di fidarti del "verde")

- ⚠️ La **CI GitHub** può essere ferma per billing: se i job non partono, i gate
  CI (sanitizers, benchmark, debt-gate, property) **non sono attivi** — non assumere
  copertura. Verifica con `gh run list` a inizio sessione.
- ⚠️ `scripts/test_quick.sh` esclude una famiglia di test con **regressioni perf reali**
  e un baseline che fallisce: un "verde" lì **non** certifica quelle aree. Vedi le note
  inline nello script e i ledger `HC-F8-FACTORIZATIONTOWER-PERF` / `F2GateBenchmark`.
- Oracolo forte: golden diff vs Maxima 5.49.0 (`scripts/run_golden_maxima.sh`,
  fork/exec only — CLAUDE.md §6). Usalo per validare equivalenza, non `toString()`.

## §5 — A fine sessione

Aggiorna `STATE.md` e sposta in `PLAN_NEXT_SESSIONS.md` gli item "DONE" sotto §Log,
promuovendo la sessione successiva a `(CURRENT)`. Non creare nuovi `*_SESSION_SUMMARY_*.md`
(lo storico vive in `docs/archive/`).
