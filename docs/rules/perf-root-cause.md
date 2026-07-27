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
