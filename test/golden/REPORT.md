# Golden Test Report — FASE 0.5

Generated: 2026-05-24 17:35 UTC

## Pass-Rate by Area

| Area | Pass | Fail | Skip | Total | Pass% | Delta |
|------|------|------|------|-------|-------|-------|
| simplify   |   35 |    5 |    0 |    40 |  87.5% |    new |
| factor     |   30 |    0 |    0 |    30 | 100.0% |    new |
| gcd        |   17 |    0 |    3 |    20 | 100.0% |    new |
| integrate  |    8 |   26 |   36 |    70 |  23.5% |    new |
| diff       |   20 |    7 |    3 |    30 |  74.1% |    new |
| limit      |   22 |    3 |    5 |    30 |  88.0% |    new |
| solve      |    0 |    0 |   20 |    20 |   0.0% |    new |
| series     |    0 |   16 |    4 |    20 |   0.0% |    new |
| special_fn |   11 |    8 |    1 |    20 |  57.9% |    new |
| **TOTAL**  |  143 |   65 |   72 |   280 |  68.8% |    new |

## Failure Examples

### simplify

- `0^x | CAS: 0^x | Maxima: 0`
- `log(1) | CAS: log(1) | Maxima: 0`
- `log(exp(x)) | CAS: log(exp(x)) | Maxima: x`
- `exp(log(x)) | CAS: exp(log(x)) | Maxima: x`
- `abs(x^2) | CAS: abs(x^2) | Maxima: x^2`

### integrate

- `integrate(1/(x^2 - 1), x) | CAS: (1/2) * ln(abs(x + -1)) + (-1/2) * ln(abs(x + 1)) | Maxima: (-1/2) * log(x + 1) + (1/2) * log(x + -1)`
- `integrate(1/(x^2 + x), x) | CAS: ln(abs(x)) - ln(abs(x + 1)) | Maxima: -log(x + 1) + log(x)`
- `integrate(1/(x^3 - x), x) | CAS: -ln(abs(x)) + (1/2) * ln(abs(x + -1)) + (1/2) * ln(abs(x + 1)) where { v1 = 1/2 } | Maxima: -log(x) + (1/2) * log(x + -1) + (1/2) * log(x + 1) where { v1 = 1/2 }`
- `integrate(x/(x^2 - 1), x) | CAS: (1/2) * ln(abs(x^2 - 1)) | Maxima: (1/2) * log(x^2 + -1)`
- `integrate(1/(x^2*(x + 1)), x) | CAS: ln(abs(x + 1)) + (-1/1) * x^-1 - ln(abs(x)) | Maxima: -x^-1 - log(x) + log(x + 1)`

### diff

- `diff(exp(x), x) | CAS: exp(x) | Maxima: e^x`
- `diff(tan(x), x) | CAS: ((1/2) * cos(2 * x) + (1/2) where { v1 = 1/2 })^-1 where { v1 = 1/2 } | Maxima: sec(x)^2`
- `diff(x*exp(x), x) | CAS: x * exp(x) + exp(x) where { v1 = exp(x) } | Maxima: x * e^x + e^x where { v1 = e^x }`
- `diff(exp(x^2), x) | CAS: 2 * x * exp(x^2) | Maxima: 2 * x * e^x^2`
- `diff(sin(exp(x)), x) | CAS: cos(exp(x)) * exp(x) where { v1 = exp(x) } | Maxima: e^x * cos(e^x) where { v1 = e^x }`

### limit

- `limit(log(x), x, 1) | CAS: log(1) | Maxima: 0`
- `limit(log(1 + x)/x, x, 0) | CAS: inf | Maxima: 1`
- `limit(log(x)/sqrt(x), x, inf) | CAS: sqrt(inf)^-1 * log(inf) | Maxima: 0`

### series

- `series(exp(x), x, 0, 5) | CAS: (1/120) * x^5 + (1/24) * x^4 + (1/6) * x^3 + (1/2) * x^2 + x + 1 | Maxima: series(e^x, x, 0, 5)`
- `series(sin(x), x, 0, 6) | CAS: (1/120) * x^5 + (-1/6) * x^3 + x | Maxima: series(sin(x), x, 0, 6)`
- `series(cos(x), x, 0, 6) | CAS: (-1/720) * x^6 + (1/24) * x^4 + (-1/2) * x^2 + 1 | Maxima: series(cos(x), x, 0, 6)`
- `series(1/(1 - x), x, 0, 5) | CAS: x^5 + x^4 + x^3 + x^2 + x + 1 | Maxima: unknown((-x + 1)^-1, x, 0, 5)`
- `series(sqrt(1 + x), x, 0, 4) | CAS: (-5/128) * x^4 + (1/16) * x^3 + (-1/8) * x^2 + (1/2) * x + 1 | Maxima: series(sqrt(x + 1), x, 0, 4)`

### special_fn

- `factorial(5) | CAS: factorial(5) | Maxima: 120`
- `factorial(10) | CAS: factorial(10) | Maxima: 3628800`
- `binomial(5, 2) | CAS: binomial(5, 2) | Maxima: 10`
- `binomial(10, 3) | CAS: binomial(10, 3) | Maxima: 120`
- `binomial(n, 0) | CAS: binomial(n, 0) | Maxima: 1`

## Notes

- **SKIP**: our CAS returned `Unimplemented`, Maxima output missing/unparseable,
  or equality check was inconclusive.
- **Pass%** is computed over decided inputs (pass + fail); SKIP excluded.
- Delta compares against previous `REPORT.prev.md` (if available).
- Maxima reference: 5.49.0 (GPL-2.0-only, unmodified).
- Corpus seed: 1026 inputs across 11 areas (expanded from 280 on 2026-05-25).

## Corpus Expansion Summary (2026-05-25)

| Area              | Before | After | Delta |
|-------------------|--------|-------|-------|
| simplify          |     40 |   116 |  +76  |
| factor            |     30 |    99 |  +69  |
| gcd               |     20 |    81 |  +61  |
| integrate         |     30 |   140 | +110  |
| diff              |     30 |    80 |  +50  |
| limit             |     30 |    99 |  +69  |
| solve             |     20 |    81 |  +61  |
| series            |     20 |    81 |  +61  |
| special_fn        |     20 |    80 |  +60  |
| matrix            |     20 |    79 |  +59  |
| bronstein         |     40 |    90 |  +50  |
| **TOTAL**         |**300** |**1026**|**+726**|

All inputs validated: JSON-parseable, no duplicates, `input` + `area` fields present.
