# DISCIPLINA DI VERIFICA — come si PROVA una chiusura

> Dettaglio normativo su cosa conta come **prova** che una task è chiusa.
> L'indice vive in `CLAUDE.md`. Carica questo file prima di dichiarare
> completata una task, e sempre quando erediti claim da ledger/tasklist.
>
> **Motivazione storica**: chiusure delegate rivelate false a suite verde;
> claim di ledger smentiti alla verifica empirica (A27 "clean bail" che era
> silent-wrong; A8 "perf blow-up" che era shape-bug; "no factor" della
> ricombinazione che non implicava irriducibilità).

## Regola d'oro

Una task è chiusa quando il **SINTOMO originale, rimisurato, è sparito** —
non quando la suite è verde. La suite verde è condizione necessaria, mai
sufficiente: il sintomo può essere un WARN debug-only, una misura di
performance o un comportamento che nessun test copre.

## Protocollo

1. **Sintomo prima/dopo**: riproduci il sintomo su un worktree al commit
   precedente al fix, poi sul fix. Stesso comando, stessi input, stesse
   condizioni di carico. Riporta entrambe le misure.
2. **Equivalenza matematica dubbia** → certificato numerico multi-punto ad
   alta precisione (skill `numeric-certify`, mpmath). **Limite noto**: il
   multi-punto NON prova identità che coinvolgono trascendenti in generale
   (lezione A47) — lì serve un argomento strutturale o un certificato esatto.
3. **Claim ereditati** (ledger, tasklist, commit message: "clean bail",
   "perf blow-up", "no factor", "già coperto") → verificare empiricamente
   PRIMA di costruirci sopra. Leggere sempre `file:riga` reali prima di
   scrivere spec o piani (lezione A31).
4. **Misure golden**: dati freschi obbligatori — `cas_golden_runner` NON si
   ricompila da solo, va rebuildato dopo modifiche al codice o si misura
   stale; ratchet solo in foreground; guard di freschezza attivo.

## Trappole note (repo-specifiche)

- `polynomial_normal_form` NON è idempotente sotto `operation_active_`:
  helper zero-diff dentro `mathematically_equal` danno falsi negativi;
  pulire `operation_active_` al call site.
- Test di no-pollution delle side-conditions: servono DUE emittenti di
  condizioni **diverse** — uno spacer neutro maschera il bug (lezione A31).
- Test flaky sotto carico (timeout wall-clock, caso A30): verificare in
  isolamento prima di attribuire una regressione al proprio diff.
- Mai `toString()` per validare la logica (costituzione §Testing): solo
  confronto strutturale o equivalenza matematica.
- Nessuna verifica (build/test/benchmark propri) mentre un gate gira in
  background: timeout falsi-positivi e risultati non attendibili.
