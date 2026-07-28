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

## Un gate deve misurare il codice, non la macchina (A51)

Un gate il cui esito dipende dal carico è peggio di nessun gate: dà un verde
che non significa niente e un rosso che non si riproduce.

1. **Ripetere la misura non basta come prova di determinismo.** Due run
   consecutive identiche sulla stessa macchina scarica non escludono che il
   verdetto sia deciso dal tempo. La prova forte è l'**invarianza rispetto al
   budget**: stessa misura con cap diversi (30 s vs 60 s vs 300 s) deve dare
   gli stessi verdetti. In A51 non li dava — e a cap diversi cambiavano entry
   *diverse*, cosa che due run identiche non avrebbero mai rivelato.
2. **Un esito troncato non è un verdetto.** Se il motore viene interrotto a
   metà, ciò che il runner registra (`false` da un confronto incompleto,
   `NO_STRATEGY` da una strategia interrotta) parla della macchina, non della
   matematica. Va classificato in una categoria propria (`over_budget`) e
   tenuto fuori da pass/fail — con un suo tetto nel ratchet, o si "passa" il
   gate lasciando scadere le entry scomode.
3. **La soglia si sceglie nel vuoto della distribuzione, misurandola.** Cap
   utile = margine ≥3× fra l'entry decisa più lenta e la soglia, con nessuna
   entry nella fascia intermedia. Serve strumentazione (`--ops-report`:
   ops + ms per entry): senza dati la soglia è un numero preso a intuito e le
   entry al confine restano invisibili finché non oscillano.
4. **Un budget che non morde mai non è un gate.** Il gate deterministico di
   A30 (`max_operation_ops`) non ha mai deciso nulla nel golden runner: il
   wall-clock scattava a ~400k ops contro un tetto di 2M. Verificare sempre
   *quale* dei due limiti taglia per primo, sul lavoro reale.

## Il gate vale solo per il codice su cui è girato

Una suite verde **non si eredita** attraverso modifiche successive. Nel caso
A51 la quick suite era verde, poi sono arrivate altre modifiche al motore, e il
commit è partito sulla fiducia in quel verde: la suite rieseguita dopo ha
trovato uno stack overflow. Regola operativa: l'ultima esecuzione della suite
completa deve essere **posteriore all'ultima modifica** a `src/` o `include/`;
se una modifica arriva dopo, il gate va rifatto, anche se "tocca solo commenti"
(nel caso reale erano commenti *e* uno split di header).

Corollario: i test mirati (`--gtest_filter`) servono a iterare, non a chiudere.
Un fix può essere verde sui test della sua area e rompere un'altra area — e
succede soprattutto quando si stringe un budget condiviso.

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
- **Il NOME di un test può escluderlo dal gate.** `test_quick.sh` filtra via
  `*Stress*`, `*Fuzz*` e `*Disabled*` (convenzione gtest): un test chiamato
  `...DisabledOpsGate...` non sarebbe mai girato, restando verde per non
  esistere (A53, intercettato contando i `[ RUN ]`). Dopo aver aggiunto test:
  contarli nel log del gate, non fidarsi del totale — `2894` invece di `2895`
  è l'unico segnale che ne manca uno.
- **Un "risultato sbagliato" può essere il VERIFICATORE che sbaglia** (A54).
  `D(F) ≠ f` accusa il produttore (Risch), ma l'antiderivata era corretta —
  certificata mpmath a 40 cifre in 7 punti — e a fallire era `is_zero_expr`,
  perché `expand` non garantiva la propria post-condizione (nessun `Sum` sotto
  un `Product`) e una differenza identicamente nulla non si cancellava.
  **Prima di aprire il codice del produttore, certificare numericamente il suo
  output**: se il valore è giusto, il difetto è a valle. Corollario: quando un
  verdetto dipende dalla FORMA dell'input, sospettare per primo chi normalizza,
  non chi decide.
- **Un flag che cambia il risultato di `simplify` deve entrare nella chiave di
  cache, o bypassarla** (A54). `simplify_cache_` è indicizzata sul solo
  `ExprPtr`: un parametro nuovo che altera l'output rende gli hit precedenti
  risposte sbagliate — e il difetto si ripresenta *solo* quando la cache è
  calda, cioè nei test end-to-end e mai nei microtest a contesto fresco.
- **Il benchmark è un gate solo a macchina scarica.** Eseguito nella coda di
  una misura golden (load 30 su 10 CPU) produce numeri gonfiati; lo script
  avvisa, ma l'avviso si perde nel log di build. Attendere `load < 4`, e usare
  `--check --baseline <file>`: senza `--check` lo script stampa le misure e
  NON confronta nulla (A53).
