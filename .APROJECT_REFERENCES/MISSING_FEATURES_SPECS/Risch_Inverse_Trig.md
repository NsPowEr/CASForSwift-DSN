# F7.5.B1 — Inverse trig + inverse hyperbolic standalone integrals

> Closure parziale del gap Risch coverage. Estende
> `integrate_function_direct` con primitive note delle inverse
> circolari/iperboliche per ridurre `INTEGRATE_NO_STRATEGY` su corpus
> integrate / bronstein.

## Problema

Corpus integrate area entry standalone `asin(x)`, `acos(x)`, `asinh(x)`,
`acosh(x)`, `atanh(x)`, `acoth(x)` → `INTEGRATE_NO_STRATEGY`.

`integrate_function_direct` (`src/calculus/integrate_elementary.cpp`)
copre già `Atan`, ma non le altre inverse circolari né le inverse
iperboliche (queste ultime non sono in `BuiltinOp` — `Asinh/Acosh/
Atanh` non enumerati, parsati come `FuncCall("asinh", x)` con
`BuiltinOp::Custom`).

## Specifica formale

### Primitive (Bronstein "Symbolic Integration I" cap. 2 + Gradshteyn-Ryzhik)

| Integranda | Primitiva |
|---|---|
| `asin(x)`  | `x·asin(x) + sqrt(1 - x²)` |
| `acos(x)`  | `x·acos(x) − sqrt(1 - x²)` |
| `atan(x)`  | `x·atan(x) − ½·ln(1 + x²)` (già presente) |
| `asinh(x)` | `x·asinh(x) − sqrt(x² + 1)` |
| `acosh(x)` | `x·acosh(x) − sqrt(x² − 1)` |
| `atanh(x)` | `x·atanh(x) + ½·ln(1 − x²)` |
| `acoth(x)` | `x·acoth(x) + ½·ln(x² − 1)` |

Derivazione: integration by parts strutturale, `u = trig_inv(x)`,
`dv = dx`. Non è una "tabella chiusa" — è la **forma esatta** della
by-parts su queste funzioni, derivabile via algoritmo Risch su
elementary extensions (Bronstein cap. 5-6) ma più cheaply riconoscibile
come primitive note quando l'argomento è `x` o affine in `x`.

### Argomento affine

`integrate_function_direct` viene già chiamato da
`integrate_function` solo se `extract_affine_argument(arg, var)`
torna `coefficient != 0`. Le primitive sopra valgono per `argument =
x` e si generalizzano a `a·x + b` via scalatura:

```
∫f(a·x + b) dx = (1/a)·F(a·x + b)
```

dove `F` è la primitiva per argument = `x`. Lo scaling è già fatto da
`integrate_function` (linea 100-104).

### BuiltinOp::Custom path

Le inverse iperboliche non hanno `BuiltinOp` proprio. Il dispatch in
`integrate_function_direct` matcha per nome (`call.name`):

```cpp
const std::string& n = name;
if (n == "asinh") return ...;
if (n == "acosh") return ...;
// ...
```

In aggiunta agli switch su `BuiltinOp` esistenti.

### x · trig_inv(x) e similari (n ≥ 1)

`integrate(x · atan(x), x)` etc. sono Product nodes. Vengono routati a
`integrate_by_parts` che, post-fix, prova ILATE selection. ILATE
priority deve preferire `trig_inv` come `u`. Verificare:

```cpp
get_ilate_priority:
  Inverse  > 4 (trig_inv: asin/atan/acos…)
  Log      = 4
  Algebraic= 3
  Trig     = 2
  Exp      = 1
```

Se ILATE già OK, `integrate(x·atan(x), x)` chiamerà by-parts con
`u = atan(x)`, `dv = x dx` → `½x²·atan(x) - ∫(½x²/(1+x²))dx =
½x²·atan(x) - ½(x - atan(x))`. La sotto-integrazione è razionale
elementare, gestita da `integrate_quadratic_rational` o
divisione polinomiale.

### No-shortcut

- Vietato pattern `if (input == "asin(x)") return "x*asin(x) + sqrt(1-x^2)"`
  con literal match. La logica deve essere `(name, argument)` ed essere
  generica nell'argomento.
- Vietato omettere lo scaling per `a·x + b`.
- Le inverse iperboliche `asinh/acosh/atanh/acoth` devono essere
  riconosciute sia con prefisso `a` sia con la forma alias
  `arcsinh/arccosh/arctanh` se Maxima li emette.

## Acceptance criteria

- Corpus integrate area: ≥ 8 entry SKIP correntemente
  (`integrate(asin(x), x)`, `acos`, `atan` (test pre-esistente),
  `asinh`, `acosh`, `atanh`, `acoth`, e composizioni `x·atan(x)`
  via by-parts) passano.
- Corpus integrate pass-rate non-skip: 45.2% → ≥ **52%** (binding floor).
- Bronstein corpus: ≥ 2 ulteriori entry passano.
- Unit test `test/unit/calculus/test_integrate_inverse_trig.cpp`
  (≥ 14 casi): primitive correnti × 7 funzioni; più scaling
  `integrate(asin(2x), x) = (1/2)·(2x·asin(2x) + sqrt(1 - 4x²))`
  semplificato; più Product `integrate(x·atan(x), x)`.
- Nessuna regressione su altre aree.

## File da modificare/creare

- `src/calculus/integrate_elementary.cpp` — estendere
  `integrate_function_direct` con i 6 nuovi casi (asin/acos +
  asinh/acosh/atanh/acoth). Anti-monolito: 213 LOC oggi, +60 stimati
  → 273 LOC, sotto limite.
- `test/unit/calculus/test_integrate_inverse_trig.cpp` (nuovo).
- `CMakeLists.txt` — registra unit test.

## Vincoli (CLAUDE.md)

- 500 LOC max per file.
- Result<T>, no throw/catch.
- BigInt only.
- Structural sharing nei `make_*` already arena-aware.
- No hardcode-of-passage. Le costanti nelle primitive (`1`, `2`, `-1`)
  sono coefficienti matematici, non parametri.
- Spec read first (REGOLA 0.1).
