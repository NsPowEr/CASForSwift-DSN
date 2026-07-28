# DIAGNOSI PERFORMANCE — Protocollo Root-Cause (anti "alza il budget")

> Dettaglio normativo del metodo di diagnosi per problemi di performance
> (lentezza, timeout, budget/depth raggiunti). L'indice vive in `CLAUDE.md`.
> Carica questo file PRIMA di toccare qualsiasi budget, depth, timeout o
> soglia, e PRIMA di ottimizzare dopo una regressione benchmark.
>
> **Motivazione storica**: 3 diagnosi perf registrate nel ledger si sono
> rivelate TUTTE sbagliate alla misura (A6/A25/F2Gate); il "perf blow-up" di
> A8/Kovacic era in realtà un bug di shape-matching (Cat-8); A42/A47 hanno
> mostrato che il collo era la strategia di verifica, non il motore. Il fix
> istintivo "aumenta il budget" cura il sintomo e nasconde il bug.

## Regola d'oro

Un budget raggiunto è un **SINTOMO**, non una causa. È vietato alzarlo prima
di aver completato i 4 passi sotto, con evidenze misurate. Alzarlo è legittimo
SOLO con giustificazione scritta nel commit: *"l'algoritmo è ottimale per
questa classe di input; il costo è intrinseco e misurato; il margine serve per
\<classe di input documentata\>"*.

---

## Passo 1 — Misura fasi + contatori (mai a occhio)

- Instrumenta il percorso sospetto: timer per fase + contatori (numero
  chiamate, nodi visitati, dimensione AST in ingresso/uscita).
- **Vietato** attribuire il costo per intuizione o analogia ("somiglia al
  caso X"): nel progetto 3 diagnosi su 3 sono state smentite dalla misura.
- Output atteso: tabella `fase → ms → #chiamate`. Il collo reale è spesso
  diverso da quello sospettato.
- **Vietato** misurare mentre un gate/build gira in background: la quick suite
  è già vicina al cap pulita, il carico produce timeout falsi-positivi.

## Passo 2 — Raggiungibilità: shape-bug o carico genuino?

Domanda obbligatoria: *il ramo lento (o il bail-out) è raggiunto perché la
**forma** dell'input non è quella che il dispatch si aspetta, o perché il
carico è algoritmicamente reale?*

Casi reali del progetto:
- `Pow(Q,-1)` e `Product(N,Pow(D,-1))` non arrivavano a `integrate_rational`
  → 40s (path generico) vs 0.1s (path razionale). La forma decideva il costo.
- A8/Kovacic: `extract_pole_loc` non vedeva i poli in forma `Sum` canonica →
  diagnosi "perf blow-up" errata; il fix era un parse shape-independent.

Check pratico: stampa la shape AST esatta all'ingresso del ramo; confrontala
col pattern atteso dal dispatch. Gli helper di estrazione basati su
`expr_cast<...>` su forme fisse sono i primi sospetti.

## Passo 3 — Lavoro duplicato: memo mancante, costo dentro loop

- Il costo totale è N volte la stessa operazione? → memoizzazione strutturale
  (l'identità del puntatore è O(1) grazie allo structural sharing) oppure
  hoisting dell'operazione fuori dal loop.
- Casi reali: A47 — `apart_num_den` costava 5.2s per ~6 `simplify` per passo
  dentro la verifica IBP; A42 — `diff()` su un'antiderivata da 526 nodi nella
  verifica trasformava 405ms in 192s. In entrambi il fix vero è cambiare
  **strategia** (es. linearizzazione trigonometrica), non il budget.
- Attenzione: alcuni `simplify` intermedi servono alla LOGICA dell'algoritmo,
  non solo all'estetica (lezione A47) — verificare prima di rimuoverli.

## Passo 4 — Solo ora: complessità intrinseca e budget

- Se i passi 1–3 escludono shape-bug e duplicazione: stima la crescita
  empirica (raddoppia la dimensione dell'input, misura il rapporto).
- Crescita **esponenziale nel budget di ricorsione** (caso A49) = struttura
  ricorsiva sbagliata (stesso sottoproblema riesplorato), NON budget stretto.
- Bump legittimo: sempre via `CASContext` (Cat-1 del catalogo hardcode), mai
  costante ricompilata; commit con classe di input, misura, e perché
  l'algoritmo è già ottimale.

## Passo 5 — Chi possiede il budget? (A51)

Un budget vincola solo l'operazione che lo **apre**. Nel progetto le ops si
contano in due soli punti (`Simplifier` e `Substituter`), e ogni operazione
top-level azzera contatore e timer: un motore che chiama `ctx.simplify()`
migliaia di volte senza aprire un'operazione propria resta quindi **senza
alcun limite complessivo**, perché ogni passo riparte da zero. Misurato su
`calculus::integrate`: consumava per intero qualunque cap gli si desse (30 s,
60 s, 300 s), e il risultato cambiava col budget invece che con l'integranda.

- Domanda obbligatoria prima di attribuire un costo a "l'algoritmo è lento":
  *questo motore apre un'operazione con budget proprio, o eredita quella del
  chiamante?* Se nessuno la apre, il budget che credi attivo non esiste.
- Fix strutturale: `CASContext::OperationScope` (RAII rientrante) all'ingresso
  pubblico del motore. Rientrante è essenziale — una chiamata annidata non
  deve riaprire il budget, o il costo dei rami interni non viene addebitato.
- Verifica quale dei due limiti taglia per primo sul lavoro reale
  (`--ops-report` nel golden runner): un gate ops che non viene mai raggiunto
  prima del wall-clock non è deterministico, è decorativo.

**Stringere un budget ESPONE difetti latenti a valle.** Non è un effetto
collaterale raro, è la norma: il codice a valle è stato scritto assumendo che
`simplify` riesca sempre. Chiudendo il budget di `integrate` (A51/A53) sono
emersi in un colpo solo una ricorsione non terminante (`poly_extended_gcd`
contava sulla riduzione dei coefficienti per far calare il grado → stack
overflow) e un cambio di semantica delle side-conditions (le condizioni dei
rami scartati sopravvivevano nel risultato). Perciò: dopo un cambio di budget
si esegue la suite **completa**, mai i soli test mirati — e ogni difetto
esposto va giudicato per sé, non assorbito allargando di nuovo il budget.

**Una soglia si calibra sul VUOTO della distribuzione, non su un multiplo.**
Misurare il costo reale con il gate spento (`--ops-report --max-ops 0`),
separare le esecuzioni **decise** da quelle troncate, e cercare l'intervallo
dove non cade nulla. Per il budget di `integrate` (A53): legittime fino a
310'184 ops, prima patologica a 700'911, nessuna entry decisa in mezzo →
soglia 500'000 (1.61× sopra il massimo legittimo, 1.40× sotto il primo caso
patologico). È lo stesso criterio con cui A51 scelse il cap per-entry di 60 s.
Un multiplo scelto a intuito sul solo massimo osservato non dice **niente** su
dove comincia la patologia.

**Rimisura il sintomo prima di aprire il codice: la task può descrivere uno
stato non più vero.** A47 era nata su `verify:simplify` a 60 s (cap) e
`verify:together` a 46 s; alla rimisura, tre giorni e due task dopo, erano 70 ms
e 2.9 s — il difetto d'origine era già stato chiuso da A53. Senza quella misura
il lavoro sarebbe stato attribuito alla task sbagliata e diretto al punto
sbagliato. Vale anche il contrario: il costo residuo era reale, ma altrove
(`apart_num_den`, 96% del test).

**Un'ottimizzazione che cambia il RISULTATO non è un'ottimizzazione.** In A47
l'identità esatta `N₁/D + N₂/D = (N₁+N₂)/D` valeva un altro 25%, ma cambiando la
forma passata a valle faceva produrre a Risch un'antiderivata sbagliata. La
regola non è "l'identità era sbagliata" (era esatta): è che un guadagno di costo
non si paga con un silent-wrong. Si rimuove, si apre il difetto esposto come
task propria (A54), e la si riattiva quando quella chiude.

**Costo per livello ≠ numero di livelli.** In una ricorsione, prima di attaccare
la profondità misura il costo del SINGOLO livello: in A49 le tre ipotesi naturali
sono state tutte smentite dai contatori — passi di ricorsione identici (18 contro
18, e il caso peggiore ne faceva 9), espressioni di dimensione costante (6 contro
8 nodi), `expand`/`parse_polynomial` a +10% — mentre il costo reale era che ogni
livello **rilanciava l'intera pipeline** sul sotto-problema (56% dei campioni in
`integrate_risch`). Il campionamento (`sample <pid>`) dice dove si sta, i
contatori dicono quante volte: servono entrambi, e vanno confrontati a parità di
unità (contare le RIGHE di un albero di `sample` non è contare i campioni).

**Il fix migliore è spesso un teorema, non un'ottimizzazione.** Se una strategia
non può avere successo per ragioni matematiche, non ottimizzarla: non chiamarla.
`∫e^{P(x)}·R(x)` con un polo genuino non ha primitiva elementare (Liouville) →
l'intera catena per parti era lavoro garantito sprecato, e il guard che la
esclude vale 27-50× (A49). Prima di ottimizzare un ramo caro, chiedersi se quel
ramo poteva concludere qualcosa.

**Un budget morde solo dove la grandezza si conta.** Le ops si incrementano in
`Simplifier` e `Substituter`: un'integranda il cui costo sta in Risch o
nell'algebra dei polinomi resta limitata dal solo wall-clock, per quanto si
stringa il tetto (misurato in A53: `sin(log x)·cos(log x)/x³` gira per minuti
con un tetto di 20'000 ops, perché quelle ops non le spende mai). Prima di
concludere che "il budget non funziona", verificare con `ops_high_water()`
**se il contatore avanza affatto** su quell'input.

---

## Ordine di preferenza dei fix

1. **Fix del dispatch/shape** — canonicalizzazione all'ingresso, parse
   shape-independent (il più frequente nel progetto).
2. **Eliminazione lavoro duplicato** — memo, hoisting, strategia alternativa
   più economica.
3. **Algoritmo migliore** — complessità inferiore.
4. **Budget più alto via `CASContext`** — ultima risorsa, con giustificazione
   scritta.

## Divieti

- Vietato qualunque "fix" che cambia un numero senza sapere dove va il tempo.
- Vietato dichiarare chiusa una task perf senza rimisurare il SINTOMO
  originale → [`verification-discipline.md`](verification-discipline.md).
- Vietato aggiornare `baseline_release.txt` per assorbire una regressione non
  diagnosticata (skill `benchmark-gate`).
