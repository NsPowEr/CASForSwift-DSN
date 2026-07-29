---
name: stale-work-recovery
description: Presa in carico di lavoro fermo o a metà — di altre sessioni AI, dell'utente, o interrotto — trovato nel tree/stash/branch. Ricostruisce il ragionamento sottostante da evidenze (diff, spec, task, memoria), decide la strada corretta, e COMPLETA definitivamente quel codice al 100% senza lasciare debiti, con gate pieni. Da invocare quando il pre-flight di next-task trova un tree sporco non proprio, uno stash orfano, o codice con TODO/stub recenti.
---

# stale-work-recovery — completare il lavoro degli altri, senza debiti

Lavoro a metà è un debito: blocca i passi successivi e marcisce (API sotto di
lui cambiano, il contesto si perde). Questa skill lo trasforma in lavoro FINITO
o, quando completarlo è impossibile, in stato esplicito e recuperabile. Mai
lasciarlo com'è, mai cancellarlo.

## Trigger

- Pre-flight `next-task` trova tree sporco non prodotto da questa sessione.
- `git stash list` con stash non documentati.
- Branch avanti rispetto a main con lavoro non mergiato e nessuna sessione attiva.
- Hook Stop segnala file sorgente sporchi stabili da più turni e nessun'altra
  sessione risulta attiva.
- L'utente indica esplicitamente lavoro suo o di altra AI da finire.
- Codice incontrato con `Unimplemented`/stub/`HARDCODE-OF-PASSAGE` che BLOCCA
  il task corrente.

## Fase 1 — Forense (ricostruire il ragionamento, non indovinarlo)

Ordine delle fonti, dalla più affidabile:

1. **Il diff stesso**: `git diff` + `git diff --stat` file per file. Leggere
   TUTTO il codice toccato, non un campione. Per ogni hunk: cosa fa, cosa
   presuppone, cosa manca.
2. **Traiettoria**: `git log --oneline -15` sul branch — i commit precedenti
   dello stesso filone dicono dove il lavoro stava andando (es. "Step 1-4
   committati" ⇒ il diff è lo Step 5).
3. **Tracker**: task pertinente in `TASKLIST_MASTER.md` (grep per file/modulo);
   `python3 scripts/ledger_index.py search "<modulo>"` per debiti collegati.
4. **Spec**: la spec formale in `.APROJECT_REFERENCES/MISSING_FEATURES_SPECS/`
   riferita dal task (REGOLA 0.1 vale anche per lavoro ereditato: leggerla
   PRIMA di completare).
5. **Memoria di progetto**: note su quel filone (lezioni, dead-end già battuti
   — non ripercorrerli).
6. **Test**: test nuovi/modificati nel diff = specifica ESEGUIBILE dell'intento.

Output fase 1 (obbligatorio, nel log): **Ricostruzione** = intento, stato reale
(cosa è finito / cosa è a metà / cosa è rotto), strada che l'autore stava
seguendo, e verdetto: la strada è quella corretta secondo spec+costituzione?

## Fase 2 — Decisione (con evidenze, non opinioni)

- **Strada corretta, lavoro incompleto** → caso normale: completare (Fase 3).
- **Strada corretta ma con bug** → completare + fixare; il bug va capito
  (root cause), non cerottato.
- **Strada sbagliata** (viola spec/costituzione/architettura) → NON buttare in
  silenzio: documentare perché con file:riga ed evidenze, salvare il lavoro
  (`git stash push -m "descrizione"` o attic), proporre la strada giusta.
  Se il rework è grosso → task in TASKLIST_MASTER.md e conferma utente.
- **Impossibile ricostruire l'intento** → NON completare alla cieca (rischio
  silent-wrong): checkpoint sicuro + domanda mirata all'utente con ciò che si
  è capito e le 2-3 interpretazioni possibili.

## Fase 3 — Completamento (standard pieno, nessuno sconto)

Il codice ereditato diventa TUO a tutti gli effetti:

1. Completare seguendo spec + costituzione (zero hardcode, Result<T>, arena,
   500 righe/file). "Era già così" NON giustifica un pattern vietato: se il
   lavoro ereditato contiene hardcode non ledgered, o si sistema o si ledgera.
2. Test: quelli lasciati a metà si finiscono; coprire i path nuovi; verifica
   matematica vera (strutturale o `numeric-certify`).
3. Verifica completa: `cas-regression-guard` + `hardcode-auditor` + gate
   (`test_quick.sh`, benchmark se hot-path). Regole anti-collisione valide
   (gate-lock, no concorrenza).
4. Commit atomico che dichiara la presa in carico: cosa è stato ereditato,
   cosa è stato completato/corretto. Mai mischiare col proprio lavoro d'altro
   tipo nello stesso commit.
5. `TASKLIST_MASTER.md`: aggiornare il task esistente (o crearlo se il lavoro
   orfano non era tracciato — evidenza del perché esiste). Stato SOLO lì.

## Divieti

- MAI completare senza Fase 1 completa (completare alla cieca = silent-wrong
  peggio del lavoro fermo).
- MAI `git reset --hard` / cancellare lavoro non capito (REGOLA EVIDENCE-FIRST).
- MAI dichiarare "completato" senza gate verdi: un lavoro ereditato finito a
  metà due volte è debito al quadrato.
- MAI prenderlo in carico mentre l'autore è ANCORA attivo (sessione parallela
  viva sullo stesso filone): coordinarsi, non collidere — verificare con
  pgrep/attività recente sui file (mtime) prima di toccare.
