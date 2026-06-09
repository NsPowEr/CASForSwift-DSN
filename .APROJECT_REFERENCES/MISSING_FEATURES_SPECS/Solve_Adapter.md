# Solve Adapter — Maxima list output → CAS set equality

> Task F7.5.A1 — chiude gap test-infra area `solve` (oggi 0/0/81 SKIP).
> Engine algoritmico `solve_polynomial` / `solve_equation` già coperto
> da unit test propri; manca solo lato runner.

## Problema

Maxima oracle restituisce le soluzioni come **lista di equazioni**:

```
solve(x^2 - 4 = 0, x) → [x = -2, x = 2]
solve(x^3 - 1 = 0, x) → [x = 1, x = (-1 - %i*sqrt(3))/2, x = (-1 + %i*sqrt(3))/2]
```

Il CAS Engine restituisce `std::vector<ExprPtr>` con i valori `r_k`
(non equazioni). Confronto naïve via `mathematically_equal` su un
unico ExprPtr fallisce.

Il runner `corpus_runner` oggi filtra qualsiasi output Maxima che
inizia per `[` (vedi `maxima_parser.hpp:147-148`) → 100% SKIP su area
solve.

## Specifica formale

### Input

- **Maxima raw**: stringa `"[x = r_1, x = r_2, ..., x = r_n]"`
  - `x` può essere qualunque simbolo (non solo `x`)
  - `r_k` può contenere `%i`, `%pi`, `%e`, `sqrt(...)`, etc.
  - le `r_k` possono ripetersi (radici multiple) — Maxima per default
    NON le ripete (usa `multiplicities` separato)
  - lista vuota `[]` ammessa (nessuna soluzione reale)
- **CAS result**: `std::vector<ExprPtr>` di lunghezza `m`.

### Output

`Result<bool>` — `true` se i due multiset di soluzioni sono uguali
(matching biunivoco via `mathematically_equal`).

### Algoritmo

1. **Parse list Maxima**: tokenizer dedicato che riconosce delimitatori
   `[`, `]`, `,` rispettando parentesi annidate (NON regex naïve, che
   spezza `f(a, b)` su virgola interna). Output: `std::vector<std::string>`
   degli elementi raw `x = r_k`.
2. **Estrai RHS**: per ogni elemento, parser che split su `=`
   bilanciato (l'unico `=` top-level), restituisce `r_k` raw. Validare
   che LHS sia un simbolo (non literal).
3. **Normalizza + parse** ogni `r_k` via `normalize_maxima_output` +
   `parse_maxima_expr` esistenti → `ExprPtr`.
3b. **Dedup modulo `mathematically_equal`** su entrambi i set.
   Maxima `solve()` di default ritorna ogni radice distinta una sola
   volta (la molteplicità è in `multiplicities` separato). Il CAS
   include la molteplicità nel vettore. Confronto come **set** di
   radici uniche.
4. **Build bipartite graph**: vertici lato sinistro = soluzioni CAS,
   lato destro = soluzioni Maxima. Edge `(i, j)` se
   `mathematically_equal(cas[i], maxima[j], ctx)` ritorna `Ok(true)`.
5. **Maximum bipartite matching** (Hopcroft-Karp, `O(E·√V)`).
6. Confronto: **uguali iff** `|matching| == |cas| == |maxima|` (matching
   perfetto). Caso disuguale → riporta diff (mancanti / extra).

### Algoritmo Hopcroft-Karp (deterministico)

Implementazione in `test/golden/bipartite_matching.hpp` (header-only
template). Complessità `O(E·√V)`, sufficiente per `|V| ≤ 20`
(corpus solve ha al massimo 6 soluzioni per entry).

**Vietato**:
- greedy matching (può fallire matching valido se ordine sfortunato)
- confronto solo per cardinalità
- assumere ordine canonical (Maxima ordina per modulo, CAS per
  ordine canonical_compare → ordinamento divergente per radici complesse)

### Edge cases

| Caso | Comportamento |
|---|---|
| CAS `{r, r}` (multiplicity), Maxima `{r}` | PASS (dedup pre-comparison) |
| Maxima `[]`, CAS `{}` | PASS |
| Maxima `[]`, CAS `{r}` | FAIL (CAS spurious root) |
| Maxima `[x=1, x=1]` (multiplicity esplicita) | matching multiset OK |
| CAS produce `RootOf` simbolico | richiama `mathematically_equal` standard |
| Maxima soluzione con `%i` | parser `%i → i` esistente, OK |
| `mathematically_equal` ritorna `Unimplemented` su una coppia | edge contato come non presente; se matching non-perfetto, SKIP (non FAIL) |
| Lista Maxima contiene `false` | NULL — corpus entry SKIP |

## Acceptance criteria

- Corpus area `solve` (81 entry, `test/golden/corpus/solve/basic.jsonl`):
  - pass-rate ≥ **90%** sui non-SKIP
  - SKIP ammessi solo per `Unimplemented` interno (non per output format)
- Unit test dedicato `test/unit/golden/test_solve_set_equal.cpp` (≥ 10
  casi: vuota, singleton, doppia reale, complessa, multiplicity,
  ordering reverse, mismatch cardinalità, mismatch valori, RootOf).
- Nessuna regressione su altre aree (re-run completo).

## File da modificare/creare

- `test/golden/maxima_parser.hpp` — rimuovere filtro `s.front() == '['`
  → null, delegare a adapter.
- `test/golden/solve_set_equal.hpp` (nuovo) — algorithm sopra.
- `test/golden/bipartite_matching.hpp` (nuovo) — Hopcroft-Karp.
- `test/golden/main.cpp` — branch su `area == "solve"` chiama
  `compare_solve_sets()` invece di `mathematically_equal()` standard.
- `test/unit/golden/test_solve_set_equal.cpp` (nuovo).
- `CMakeLists.txt` — registra nuovo unit test.

## Vincoli (CLAUDE.md)

- 500 LOC max per file (header-only template solve_set_equal probabile
  ~250 LOC, bipartite_matching ~150 LOC — within budget).
- BigInt only nel comparison; Result<T> per error reporting.
- Nessun `throw/catch`. Nessun test disabled.
- Anti-monolito: se solve_set_equal cresce > 400 LOC, split
  `solve_set_equal_impl.hpp` per algoritmi interni.
- No hardcode-of-passage. Hopcroft-Karp è algoritmo standard, no
  costanti magiche.

## Note implementative

Hopcroft-Karp BFS layered + DFS augmenting paths classico:
```
while (bfs_finds_augmenting()):
    for each unmatched u on left:
        if dfs_augment(u):
            ++matching_size
return matching_size
```

`mathematically_equal` cache: per ogni edge `(i, j)`, calcoliamo una
sola volta. Memoize in `std::unordered_map<std::pair<int,int>, bool>`
per evitare re-call su DFS retries.
