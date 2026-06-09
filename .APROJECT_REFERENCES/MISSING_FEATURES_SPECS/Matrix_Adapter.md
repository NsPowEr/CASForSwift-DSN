# Matrix Adapter — `[[…]]` corpus syntax → MatrixLit dispatch

> Task F7.5.A2 — chiude gap test-infra area `matrix` (oggi 0/79 PASS).
> Engine linalg (det/trace/transpose/rank/inverse/eigenvalues) già
> coperto da unit test propri; manca solo lato runner.

## Problema

Il corpus area `matrix` usa la sintassi notation
`det([[1,2],[3,4]])`, `transpose([[a,b],[c,d]])`, ecc.

`corpus_runner.hpp:142-145` marca tutte queste come "matrix fn
skipped" prima di provarne il parsing, lasciando 0/79 PASS sul
corpus.

Maxima emette matrici come `matrix(row1, row2, ...)` che il parser
attuale filtra con `s.rfind("matrix(", 0) == 0 → return ""` (linea
143-144 di `maxima_parser.hpp`).

## Specifica formale

### Input corpus

```
det([[1,2],[3,4]])              → -2
det([[a,b],[c,d]])              → a*d - b*c
trace([[1,2],[3,4]])            → 5
transpose([[1,2],[3,4]])        → [[1,3],[2,4]]
rank([[1,2],[2,4]])             → 1
inverse([[1,2],[3,4]])          → [[-2, 1], [3/2, -1/2]]
eigenvalues([[2,0],[0,3]])      → [2, 3]
```

### Parsing `[[r1c1,r1c2,...],[r2c1,r2c2,...],...]`

Parser `parse_matrix_lit(text, arena) → Result<ExprPtr>` che:
1. Riconosce `[` outer, splitta righe su `,` rilevati a depth 0
   bilanciato (riusa logica `split_maxima_list` da F7.5.A1).
2. Per ogni riga, riconosce inner `[…]`, splitta elementi.
3. Parsa ogni elemento via `parse_expr` esistente.
4. Costruisce `MatrixLit` con `rows × cols` (verifica rettangolarità).

Vietato regex naïve (`f(a,b)` matrici interne romperebbero).

### Dispatch funzioni

In `corpus_runner.hpp::evaluate_cas`:

```cpp
if (cmd.fn == "det") {
    auto m = parse_matrix_lit(cmd.arg_strs[0], ctx);
    return cas::linalg::determinant(m, ctx);
}
// idem trace, transpose, rank, inverse, eigenvalues
```

### Maxima output `matrix(row, row, ...)`

Parser dedicato `parse_maxima_matrix(text, arena)`:
1. Strip outer `matrix(` e `)` finale.
2. Splitta righe top-level (depth 0 sui parens).
3. Per ogni riga `[a, b, …]` (Maxima usa già `[…]` interno), riusa
   `split_maxima_list` per elementi.
4. Costruisce `MatrixLit`.

Rimuovere il filtro `s.rfind("matrix(", 0) == 0` da
`maxima_parser.hpp` quando area==matrix.

### eigenvalues output

Maxima `eigenvalues(M)` ritorna `[[λ_1, …, λ_k], [m_1, …, m_k]]`
(autovalori + molteplicità). Adapter prende solo la prima sublista
(autovalori) e applica `compare_solve_sets` (riuso F7.5.A1) per
matching modulo dedup.

### transpose / inverse: ritornano matrici

Maxima `transpose(M)` → `matrix(row_t1, row_t2, …)`. Confronto:
1. Parsa entrambe le matrici (CAS + Maxima).
2. Verifica `rows × cols` uguali.
3. `mathematically_equal` su ogni elemento `[i][j]`.

Funzione helper `compare_matrices(cas, maxima, ctx) → Result<bool>`.

### rank / det / trace: ritornano scalari

Confronto standard via `mathematically_equal`.

## Edge cases

| Caso | Comportamento |
|---|---|
| Matrice 1×1 `[[a]]` | OK, scalar conversion non necessaria |
| Matrice non-quadrata su det | CAS ritorna error → SKIP |
| eigenvalues con simbolici (es. `[[a,1],[0,b]]`) | Maxima output `[[a,b],[1,1]]` → adapter usa prima sublista |
| inverse di singolare | CAS error → SKIP, non FAIL |

## Acceptance criteria

- Corpus `matrix/basic.jsonl` (79 entry): pass-rate ≥ **90%** non-skip.
- SKIP ammessi solo per:
  - eigenvalues con CAS RootOf irriducibile vs Maxima esponenziale
    (HC-F75-CYCLOTOMIC-ROOTOF già aperto).
  - inverse di matrici singolari (sia CAS sia Maxima error → entrambi
    SKIP).
- Unit test `test/unit/golden/test_matrix_adapter.cpp`:
  - parse `[[1,2],[3,4]]` → MatrixLit 2×2.
  - parse `matrix([1,2],[3,4])` → MatrixLit 2×2.
  - compare_matrices identità.
  - dispatch det/trace/transpose/rank/inverse.
- Nessuna regressione su altre aree.

## File da modificare/creare

- `test/golden/matrix_parser.hpp` (nuovo) — `parse_matrix_lit`,
  `parse_maxima_matrix`, `compare_matrices`.
- `test/golden/corpus_runner.hpp` — rimuovere skip list matrix,
  aggiungere dispatch det/trace/...
- `test/golden/maxima_parser.hpp` — bypassare il filtro `matrix(` per
  area==matrix (rimuovere o condizionare).
- `test/golden/main.cpp` — branch `area == "matrix"` analogo a
  `area == "solve"` (route a matrix-specific compare).
- `test/unit/golden/test_matrix_adapter.cpp` (nuovo).
- `CMakeLists.txt` — registrazione unit test.

## Vincoli (CLAUDE.md)

- 500 LOC max. matrix_parser.hpp stimato ~250 LOC.
- BigInt only nelle comparazioni elemento-wise.
- Result<T> ovunque.
- Nessun throw/catch.
- Structural sharing dove possibile.
- No hardcode-of-passage.
- Spec read first (REGOLA 0.1).
