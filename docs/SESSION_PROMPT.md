# Prompt canonico di apertura sessione

> Copia-incolla per riprendere le lavorazioni senza sprechi. Lo snapshot di
> stato è iniettato automaticamente dal hook `session_state.sh` a ogni
> apertura: NON chiedere "a che punto siamo" — la sessione lo sa già.

## Procedura (3 passi)

1. Apri Claude Code nella root del repo.
2. Imposta **modello + effort** secondo l'ultima raccomandazione dichiarata
   (fine sessione precedente, formato `PROSSIMA SESSIONE → ...`;
   regola: `docs/rules/model-routing.md`). Senza raccomandazione: Sonnet medium.
3. Incolla UNO dei prompt sotto.

## Prompt standard (iterazione singola autonoma)

```
/next-task
Autorizzazioni sessione: procedi in autonomia — implementazione, test, gate,
commit atomici — senza chiedere conferme intermedie. Fermati SOLO alle
stop-conditions della skill (stallo, decisione architetturale, spec mancante,
gate rosso non tuo). Dopo modifiche a integrate/simplify: check_golden_ratchet
oltre alla quick. Chiusura: TASKLIST_MASTER.md aggiornato + raccomandazione
modello/effort per la prossima sessione (docs/rules/model-routing.md).
```

## Prompt con focus (quando vuoi forzare la task)

```
/next-task
Focus: <A38 | A37 | ...> — salta la selezione, esegui questa.
[+ stesse autorizzazioni del prompt standard]
```

## Prompt continuativo (più iterazioni nella stessa sessione)

```
/loop /next-task
[+ stesse autorizzazioni del prompt standard]
```

Il loop si auto-regola (una task per iterazione, stop-conditions attive);
interrompilo quando vuoi. Per le task C4 (research-grade) preferisci la
sessione singola dedicata con Fable high, non il loop.

## Perché non serve altro nel prompt

- Stato: `session_state.sh` (SessionStart) inietta branch/dirty/ledger/ratchet.
- Protocollo: skill `next-task` (pre-flight anti-collisione, spec-check 0.1,
  verifica, gate, stop-conditions) + `stale-work-recovery` per WIP ereditati.
- Enforcement: gli hook negano comunque le azioni vietate (timeout, gate
  concorrenti, git distruttivi, rm, monoliti).
- Priorità: TASKLIST_MASTER.md sezione E + `PARITY_GIAC.md` (gap vs Giac).
