# Sech/Csch/Coth/Tanh Identity Normalization

> Task F7.5.A4 — chiude gap notazionale `diff` (cosh(x)^-2 ↔ sech(x)^2)
> e `simplify`/`limit` su identità iperboliche.

## Problema

Maxima emette forme canoniche con `sech`, `csch`, `coth` come
funzioni native. Il CAS Engine non ha `BuiltinOp::Sech/Csch/Coth`:
parsa `sech(x)` come `FuncCall("sech", x)` con `func_id ==
BuiltinOp::Custom` (opaco). `mathematically_equal` non riconosce
`1/cosh(x) ≡ sech(x)`.

Esempi corpus FAIL:
```
diff(tanh(x), x) → CAS cosh(x)^-2  vs  Maxima sech(x)^2
diff(coth(x), x) → CAS -1/sinh(x)^2  vs  Maxima -csch(x)^2
```

## Decisione architetturale

**Vietato** aggiungere `BuiltinOp::Sech/Csch/Coth` perché toccherebbe
~76 switch enum con `-Wswitch -Werror` (lavoro di Fase 8 con
Extended-Real migration, tracked F7.5.F1).

**Scelto**: rewriter pre-confronto in `algebraic_equal.cpp` che
normalizza `sech/csch/coth/tanh` come `FuncCall` Custom a forme
canoniche standard:

| Da | A |
|---|---|
| `sech(u)`   | `1/cosh(u)`       |
| `csch(u)`   | `1/sinh(u)`       |
| `coth(u)`   | `cosh(u)/sinh(u)` |
| `tanh(u)`   | `sinh(u)/cosh(u)` (opzionale, ridurre FAIL composti) |

Applicato simmetricamente su `lhs` e `rhs` PRIMA di
`polynomial_normal_form` e `split_num_den`.

## Specifica formale

### Algoritmo

`hyperbolic_normalize(expr, arena) → ExprPtr`:

1. Bottom-up traversal AST.
2. Caso `FuncCall(name, [u])` con `name ∈ {"sech","csch","coth","tanh"}`:
   - `u' = hyperbolic_normalize(u, arena)` ricorsivo.
   - Riscrittura via `make_binary(BinaryOp::Div, …)` con `cosh/sinh`
     come `FuncCall(BuiltinOp::Cosh/Sinh, u')`.
3. Casi `Sum/Product/Binary/Unary/Power`: ricorsione sui figli,
   structural sharing se nessun figlio cambia.
4. Casi atomici (`IntegerLit/RationalLit/Symbol/Constant`): identity.

### Integrazione in `mathematically_equal`

```cpp
auto lhs_h = hyperbolic_normalize(lhs_s.value(), ctx.arena());
auto rhs_h = hyperbolic_normalize(rhs_s.value(), ctx.arena());
// continua con structural_equal, normal form, etc. su lhs_h/rhs_h
```

Applicato dopo simplify iniziale, prima di structural_equal. Nessuna
modifica al simplify stesso (per non alterare canonical form).

### Garanzie

- **Idempotenza**: `hyperbolic_normalize(hyperbolic_normalize(x)) ==
  hyperbolic_normalize(x)` structural.
- **Soundness**: trasformazioni preservano valore matematico per ogni
  `u` reale o complesso fuori dai poli (cosh ≠ 0, sinh ≠ 0).
- **Structural sharing**: ritorna `ExprPtr` originale se nessun
  rewrite avviene nel sottoalbero.

## Acceptance criteria

- Corpus `diff` (80 entry): pass-rate ≥ 90% (era 82.5%).
  - 2 FAIL noti `tanh/coth/sech/csch` notazionali → PASS.
- Unit test `test/unit/algebra/test_hyperbolic_normalize.cpp`
  (≥ 8 casi):
  - `mathematically_equal(sech(x)^2, 1/cosh(x)^2) == true`
  - `mathematically_equal(csch(x), 1/sinh(x)) == true`
  - `mathematically_equal(coth(x), cosh(x)/sinh(x)) == true`
  - `mathematically_equal(tanh(x)*cosh(x), sinh(x)) == true`
  - idempotenza
  - structural sharing (puntatore identico se no rewrite)
  - nested composition: `sech(sin(x))` → `1/cosh(sin(x))`
  - no false positive: `mathematically_equal(sinh(x), cosh(x)) ==
    false`
- Nessuna regressione (run corpus completo, no area peggiora).

## File da modificare/creare

- `src/algebra/algebraic_equal.cpp` — aggiungere
  `hyperbolic_normalize` helper privato + call in
  `mathematically_equal`.
- `test/unit/algebra/test_hyperbolic_normalize.cpp` — nuovo unit test.
- `CMakeLists.txt` — registra test.

## Vincoli (CLAUDE.md)

- 500 LOC max. `algebraic_equal.cpp` oggi 286 LOC; +50 stimati.
- Result<T> ovunque. No throw/catch.
- BigInt only (no `int64_t`).
- No hardcode-of-passage (rewriter è universale).
- Structural sharing obbligatorio per identity branches.
- Spec read first (REGOLA 0.1): questo file letto prima del codice.
