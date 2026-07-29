---
name: timer-post
description: >
  Mette in pausa l'agente AI per un tempo definito dall'utente, restituendo un messaggio
  specifico al termine. L'agente deve attendere l'esecuzione completa del timer prima di
  procedere. Utilizzato tramite il comando `/timer-post messaggio;tempo`.
---

# Skill: Timer Post (Bloccante)

Questa skill permette all'utente di imporre una pausa all'AI prima di procedere con la generazione o con le azioni successive.

## Trigger

L'utente invoca la skill con la sintassi:
`/timer-post <messaggio da restituire> ; <tempo in linguaggio naturale>`

Esempio:
`/timer-post procedi con l'analisi ; 17 minuti`
`/timer-post test superato ; 1 ora e mezza`

## Workflow dell'AI (Istruzioni Rigorose)

Quando questa skill viene invocata, l'agente **deve** seguire esattamente questi passaggi:

1. **Parsing del comando:** Estrai il `<messaggio>` e il `<tempo in linguaggio naturale>`.
2. **Conversione in Secondi:** L'AI deve calcolare internamente quanti secondi corrispondono al tempo descritto (es. "17 minuti" = 1020, "1 ora" = 3600).
3. **Esecuzione dello Script:** Esegui lo script Python dedicato al timer tramite shell command:
   ```bash
   python3 /Users/REAL_CAS_ENGINE_C++/.claude/ai_timer.py "<messaggio>" <secondi>
   ```
4. **Attesa Bloccante:** **DEVI** attendere che il processo della shell termini. Lo script stamperà un piccolo "tick" ogni 60 secondi (consumando un numero minimo di token) per evitare che il tuo sistema vada in timeout. NON interrompere il processo.
5. **Azione Finale:** Quando lo script termina stamperà `[Timer Completato] <messaggio>`. Da quel momento in poi, elabora il `<messaggio>` come se te lo avesse appena scritto l'utente e procedi con il workflow richiesto.

## Regole Aggiuntive
- Non eseguire altri tool in parallelo durante l'attesa.
- Non rispondere all'utente dicendo "Sto aspettando", avvia semplicemente il comando shell, che si occuperà di mostrare il feedback tramite l'output.
