# Model & Effort Routing — regola operativa

> Obiettivo: capacità del modello proporzionata alla task. Mai Fable/Opus su
> lavoro meccanico; mai Sonnet-low su matematica research-grade. Il routing
> usa la legenda `[E·C·S·R]` di `TASKLIST_MASTER.md` (già assegnata per task).

## Tabella di routing (sessione principale)

| Classe | Natura | Modello | Effort | Esempi dal progetto |
|---|---|---|---|---|
| **C1** meccanico | rename, split anti-monolito, wiring noto, doc, whitelist | **Sonnet** | low/medium | A29 split, registrare test in CMake, aggiornare TASKLIST |
| **C2** moderato | fix guidato da spec chiara, estensioni tabellari, test nuovi, triage/misura | **Sonnet** | medium | A36, entry tabella §5, parity scan + triage A37, golden regen |
| **C3** algoritmico | algoritmo nuovo da spec completa, debugging non banale, wiring con semantica | **Opus** medium — oppure Sonnet **high** se R1 + spec molto dettagliata | medium/high | A38 wiring RDE, A13 residui, fold inversi G, A22 Padé |
| **C4** research-grade | matematica nuova, spec da scrivere/verificare, caccia silent-wrong, teoria profonda | **Fable** high (fallback: Opus high) | high/max | Risch algebrico (Bronstein §8), Mellin §6.7, Galois deg>10, A31 fase 3 |

**Modificatori** (un gradino in su di modello O di effort, non entrambi di default):
- **S5** (correttezza/hang) o sospetto silent-wrong → +1.
- **R3** (hot-path: simplifier/arena/BigInt) → effort +1 (il costo di un errore è churn di massa).
- Debug "impossibile" / 2+ sessioni fallite sulla stessa task → passa a Fable high.
- Sessione multi-task → dimensiona sul picco (la task più pesante prevista).
- Opus in **fast mode** (`/fast`): ok per C1–C2 lunghe quando serve throughput, non per C3+.

## Subagenti (già fissati nel frontmatter — NON sovradimensionare)

| Agente | Modello | Perché |
|---|---|---|
| cas-regression-guard | haiku | run-and-report meccanico |
| module-locator, spec-fetcher, hardcode-auditor, maxima/giac-golden-diff | sonnet | lookup/audit con giudizio limitato |
| spec-researcher | opus | ricerca matematica autoritativa + verifica mpmath |
| cavecrew-* | (definiti dal plugin) | già ottimizzati per costo |

Mai spawnare un subagente con modello superiore a quello della tabella per
"sicurezza": se la task del subagente sembra richiedere di più, è la task
principale a essere sotto-classificata — riclassifica quella.

## Obbligo di raccomandazione (chiusura sessione/iterazione)

Ogni fine iterazione `next-task` (e comunque ogni fine sessione) dichiara:

```
PROSSIMA SESSIONE → task: A<N> · modello: <Sonnet|Opus|Fable> · effort: <low|medium|high|max> · motivo: <1 riga con classe C e modificatori>
```

L'utente imposta modello/effort PRIMA di aprire la sessione successiva
(selettore modello + eventuale `/fast`). Se la raccomandazione manca,
default = Sonnet medium e riclassificare alla prima iterazione.

## Anti-pattern vietati

- Fable/Opus per: aggiornare markdown, rigenerare golden, lanciare suite,
  commit, parity scan, split file. (Sono C1/C2: Sonnet.)
- Sonnet low per: scrivere/validare una spec matematica, chiudere task C4,
  decidere semantica branch-cut. (Rischio allucinazione formale.)
- Cambiare modello A METÀ di una task C4 per "risparmiare": il contesto
  perso costa più del modello.
