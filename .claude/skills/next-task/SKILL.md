---
name: next-task
description: Playbook del ciclo di sviluppo autonomo del CAS Engine — seleziona la prossima task da TASKLIST_MASTER.md (guidato dal gap PARITY_GIAC), esegue spec→piano→implementazione→verifica→gate→commit, aggiorna lo stato SOLO in TASKLIST_MASTER.md. Progettato per /loop e job in background. Include stop-conditions obbligatorie.
---

# next-task — ciclo canonico di sviluppo autonomo

Un'iterazione = UNA task portata a uno stato verificabile (fatta, o avanzata a
un checkpoint committabile, o bloccata con report). Mai due task in parallelo
nella stessa iterazione.

## 0. Pre-flight (ogni iterazione, sempre)

1. `git status` — tree sporco NON tuo? → l'iterazione DIVENTA recupero: skill
   `stale-work-recovery` (ricostruzione forense + completamento definitivo del
   lavoro ereditato), MAI costruire sopra senza audit (memoria
   `uncommitted-sweep-audit`) — eccetto se l'autore è ancora attivo (sessione
   parallela viva: coordinarsi, non toccare). Tree sporco tuo di
   un'iterazione precedente → riprendi quella task, non sceglierne una nuova.
2. Gate/build di altre sessioni in corso (`pgrep -f "test_quick|debt_gate|ninja"`)
   → STOP e attendi: no lavori concorrenti ai gate (memoria
   `no-concurrent-gates-load-timeout`).
3. `python3 scripts/ledger_index.py stats` + `TASKLIST_MASTER.md` sezione A per
   il backlog attuale.

## 1. Selezione task

Priorità di scelta (in ordine):
1. Task già IN CORSO in `TASKLIST_MASTER.md` (finire prima di aprire).
2. Severità S5/S4 aperte (correttezza/hang battono le feature).
3. Task che chiudono il gap più largo in `PARITY_GIAC.md` (se lo scoreboard è
   fresco; altrimenti rigenerarlo via skill `giac-parity-scan` conta come task).
4. Ordine raccomandato della tasklist.

Dichiarare SEMPRE nel log dell'iterazione: task scelta + perché.

## 2. Spec (REGOLA 0.1 — bloccante)

- `spec-fetcher` con l'id task. Risposta `SPEC MISSING` → spawna
  `spec-researcher` per produrre la BOZZA di spec, poi STOP iterazione:
  la bozza va validata dall'umano prima di implementare.
- Frase obbligatoria di conferma lettura spec nel log.

## 3. Piano minimo

- `module-locator` per la mappa file:riga (non esplorare a mano a scala 4x).
- Task perf-class (sintomo lento/timeout/budget raggiunto) → carica
  `docs/rules/perf-root-cause.md` PRIMA del piano: protocollo 4 passi
  obbligatorio, "alza il budget" è ultima risorsa con giustificazione.
- Definire: file da toccare, test nuovi, criterio di verifica matematica
  (equivalenza strutturale o certificato numerico — skill `numeric-certify`).

## 4. Implementazione

- Regole costituzione: zero hardcode (10 categorie), Result<T>, arena, LPO,
  max 500 righe/file. Ogni costante N → self-check delle 4 domande.
- Test SEMPRE con timeout + filtro (il hook lo impone comunque).

## 5. Verifica (prima di dichiarare qualsiasi cosa)

Chiusura = SINTOMO originale rimisurato prima/dopo, non suite verde —
protocollo e trappole note in `docs/rules/verification-discipline.md`.
In parallelo dove possibile:
- `cas-regression-guard` (AcidTest+SupremeTest).
- `hardcode-auditor` sul diff.
- Test mirati nuovi + `numeric-certify` sui risultati matematici nuovi.
- Se tocchi aree golden: `giac-golden-diff` / `maxima-golden-diff` su 2-3 input
  rappresentativi.

## 6. Gate + commit

- `bash scripts/test_quick.sh` (cap 1200s) — MAI in background, MAI durante
  altri gate.
- Benchmark gate se hot-path (skill `benchmark-gate`).
- Commit atomico SOLO dei file della task (mai sweep del tree); messaggio
  dichiara eventuali HARDCODE-OF-PASSAGE.
- Se hai toccato integrate/simplify: `bash scripts/check_golden_ratchet.sh`
  OLTRE alla quick (il ratchet golden NON è nella quick — lezione 2026-07-19).
- Aggiorna `TASKLIST_MASTER.md`: stato task (✅ FATTO / 🚧 checkpoint / ⛔
  BLOCCATA + causa). NESSUN altro file di stato: TASKLIST_MASTER.md è l'unico
  tracker.
- **Promozione lezioni** (CLAUDE.md §AUTO-EVOLUZIONE): c'è una lezione
  generalizzabile in questa iterazione (diagnosi smentita, trappola, pattern
  ricorrente)? → 2-5 righe nel file `docs/rules/` pertinente, con riferimento
  al caso (`A<N>`/commit). La memoria privata NON basta: i subagent e le
  sessioni future non la vedono.
- **Raccomandazione routing (OBBLIGATORIA, ultima riga dell'iterazione)** —
  `docs/rules/model-routing.md`:
  `PROSSIMA SESSIONE → task: A<N> · modello: <Sonnet|Opus|Fable> · effort: <low|medium|high|max> · motivo: <classe C + modificatori>`

## Stop-conditions (interrompono l'iterazione, non si aggirano)

- 3 tentativi falliti sullo stesso errore → PROTOCOLLO ANTI-LOOP (Report di
  Stallo, attendi umano).
- Decisione architetturale non coperta da costituzione/spec → formulare la
  domanda per l'umano e fermarsi.
- Spec mancante (dopo bozza spec-researcher) → fermarsi.
- Gate rosso non causato dal tuo diff → segnalare, non "sistemare" test altrui.
- Budget iterazione esaurito senza checkpoint committabile → `git stash push`
  con nome descrittivo + nota in TASKLIST_MASTER.md.

## Uso con /loop

`/loop /next-task` (o job schedulato): ogni firing esegue UNA iterazione del
ciclo sopra. Il pacing lo decide il loop; questa skill definisce solo il
contenuto dell'iterazione e le stop-conditions.
