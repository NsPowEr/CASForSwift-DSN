# REAL CAS ENGINE

> An industrial-grade symbolic computer algebra system in modern C++20.

[![Build](https://img.shields.io/badge/build-passing-brightgreen)](#build)
[![Tests](https://img.shields.io/badge/tests-2307%20passing-brightgreen)](#testing)
[![HP Prime Parity](https://img.shields.io/badge/HP--Prime%20parity-94.5%25-blue)](#golden-corpus)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue)](#requirements)

---

## Overview

REAL CAS ENGINE is a from-scratch symbolic math engine targeting feature parity with HP Prime G2 / Giac. It is designed as a **pure C++ library** with a stable C API, built to be embedded in native applications — in particular as the mathematical backend of a Swift iOS/macOS app.

The engine provides exact arithmetic (BigInt / Rational), a hash-consed immutable AST with arena allocation, and a layered architecture covering symbolic simplification, polynomial algebra, calculus, linear algebra, and numerical evaluation — all without floating-point in the symbolic core.

---

## Capabilities

### Foundation
- **Exact arithmetic** — BigInt (Toom-Cook 3, Knuth Algorithm D, Montgomery CIOS modexp, Pollard p-1), Rational, BigFloat (MPFR wrapper, arbitrary precision)
- **Expression tree** — immutable, hash-consed, arena-allocated AST; structural sharing; O(1) pointer-identity equality
- **Parser / Lexer** — full expression grammar with decimal-to-rational conversion at the boundary; round-trip formatting to plain text and LaTeX

### Symbolic Simplification
- Canonical simplification of sums, products, powers, exponentials, logarithms
- Assumptions engine — `assume(x > 0)`, relational graph with 3-hop transitivity inference, domain propagation
- Transcendental normal form, radical denesting (Borodin–Fagin–Hopcroft–Tompa), abs/sign rules
- Expression caching (LRU, configurable eviction), adaptive simplification depth, cycle detection
- Complex arithmetic — exact Q[i] via `ComplexLit`, principal branch policy, Euler formula, exp/log roundtrip

### Polynomial Algebra
- Univariate factorization over Q — LLL (configurable δ), Zassenhaus recombination with budget guard, EDF (p = 2 trace polynomial), Cantor–Zassenhaus
- Multivariate GCD — Zippel sparse interpolation (n-variate), EZ-GCD, Brown modular GCD, GCDHEU with Mignotte rigorous bound
- Gröbner bases — Buchberger with Gebauer–Möller pair pruning + Sugar selection; F4/F5 reduction; reduced bases; FGLM order conversion (GRevLex → Lex); custom monomial orders (Lex, GRevLex, GLex)
- Resultant, discriminant, cyclotomic polynomials (Möbius inversion, completeness certified for deg ≤ 724)
- Sparse multivariate interpolation (Zippel n-variate, configurable retry budget)

### Calculus
- **Differentiation** — full symbolic diff, implicit diff `F(x,y)=0`, Jacobian, Hessian, gradient, partial derivatives, numeric difference formulas (Forward/Central O(h²)/O(h⁴))
- **Integration (indefinite)** — Risch algorithm (Bronstein Ch. 5–9): Hermite reduction, Trager/LRT, log-derivative recognizer, Risch structure theorem (algebraic independence of log/exp extensions), IBP, Weierstrass substitution `t = tan(x/2)`, u-substitution framework
- **Integration (definite)** — Fundamental Theorem of Calculus, singularity detection (rational poles, algebraic, transcendental), Cauchy principal value, Hadamard finite part (poles of order m ≥ 2), residue theorem (quadratic and biquadratic irreducible factors), multiple integrals via Fubini (rectangular domains), orthogonal polynomial patterns (Legendre, Hermite-H/He, Chebyshev T/U)
- **Limits** — Gruntz algorithm with recursive asymptotic comparison (no static rank), MRV set extraction, cancellation tower, ∞/∞ closure, L'Hôpital; lateral limits, signed −∞; asymptotes (vertical, horizontal, oblique)
- **Series** — Taylor/Maclaurin with generalized binomial and systematic generic fallback; Laurent series (rational fast-path + general transcendental); Padé approximants; partial fractions
- **ODEs** — 1st-order (separable, linear); Frobenius method at regular singular points; variation of parameters; Laplace-transform solver for constant-coefficient n-th order
- **Transforms** — Laplace (direct + inverse, linearity, elementary patterns)
- **Summation** — Abramov algorithm; Zeilberger (partial)
- **Improper integrals** — convergence classification, Cauchy PV, Hadamard finite part

### Linear Algebra
- Matrix arithmetic, determinant (Bareiss, PLU-pivoted), rank, RREF, inverse
- LU / PLU decomposition, QR (classical Gram–Schmidt), LDL^T (Cholesky), Smith normal form
- Jordan normal form, eigenvectors over algebraic extensions Q(α) via null-space RREF
- Symbolic eigenvalue pipeline with `RootOf` fallback for n > 3

### Algebraic Structures
- Algebraic numbers — `AlgebraicNumber` (Q(α)), `AlgebraicTowerTwoLevel` (Q(α₁, α₂)), generic n-level recursive tower `AlgebraicElement<C>`
- Factorization over algebraic extensions (Trager shift, composite resultant, 2-level tower)
- Galois groups (deg ≤ 3 via discriminant; deg ≤ 4 in progress)
- Gaussian integers Z[i] — Euclidean algorithm, canonical form, norm, unit detection

### Equation Solving
- Polynomial solver — closed-form deg 1–4; RootOf fallback for deg ≥ 5 (Abel–Ruffini)
- Transcendental / hybrid — Sturm sequence + Lipschitz dyadic refinement (`fsolve`), configurable tolerance
- System solving — Gröbner (F4), resultant-based fallback for 2-variable nonlinear systems

### Special Functions
- Gamma (integer, half-integer, reflection formula Γ(z)Γ(1−z) = π/sin(πz), functional equation)
- Riemann Zeta — exact values ζ(2k) via Bernoulli (closed-form, no table), ζ(−odd), trivial zeros
- Bessel J/Y/I/K (integer-order recurrence, BesselZero), Legendre P_n (Bonnet recurrence), Hermite H/He, Chebyshev T/U
- erf, erfc, exact-value identities; factorial, binomial, exact combinatorial branches
- Hyperbolic functions sinh/cosh/tanh/coth — exact values, parity, identities; inverse trig compositions

### Numerical
- Arbitrary-precision evaluation via MPFR (`eval_mpfr`, configurable precision 6–10 000 digits)
- Interval arithmetic (Moore convention, conservative sin/cos bounds)
- Numeric differentiation (Abramowitz–Stegun formulas)

### Units
- SI base dimensions (m, kg, s, A, K, mol, cd); 30+ named units; exact rational scaling; dimension mismatch detection

---

## Architecture

```
include/cas/          Public C++ headers (engine API)
include/cas_api.h     Public C API (Swift / FFI boundary)
src/
  foundation/         BigInt, Rational, BigFloat
  ast/                AST nodes, AstArena, hash-consing
  lexer/ parser/      Expression parsing, token stream
  symbolic/           Simplifier, assumptions, rewrite engine
  algebra/            Polynomials, GCD, factorization, Gröbner
  calculus/           Diff, integrate, limits, series, ODE
  linalg/             Matrix algorithms
  numeric/            MPFR evaluation, interval arithmetic
  numtheory/          Bernoulli, Euler phi, primality
  formatter/          Plain text and LaTeX output
  c_api/              C API bridge (cas_api.h implementation)
test/                 GoogleTest unit + property-based tests
  golden/             HP Prime G2 golden corpus (1 026 entries, Maxima oracle)
scripts/              CI helpers: benchmark, anti-monolith, milestone
```

**Key invariants enforced at all times:**

| Rule | Detail |
|------|--------|
| No floating-point in symbolic core | `BigInt` / `Rational` only; `double` forbidden |
| Structural sharing | Functions return the original `ExprPtr` when no change occurs |
| Arena allocation | All AST nodes allocated via `AstArena::make`; `new` / `make_unique` forbidden |
| Zero hardcode | Every constant must have a mathematical justification or a `CASContext` knob |
| Zero warnings | `-Wall -Wextra -Wpedantic -Werror` on all targets |
| 500 LOC gate | Each `.cpp` / `.hpp` ≤ 500 lines (active refactor; 28 files pending split) |

---

## C API

The engine exposes a stable C API designed for FFI consumption (Swift, Python, etc.):

```c
#include <cas_api.h>

CASContextRef ctx = cas_context_create();

CASExprRef expr = NULL;
cas_parse(ctx, "diff(sin(x)^2, x)", &expr);

CASExprRef result = NULL;
cas_simplify(ctx, expr, &result);

char* text = NULL;
cas_format_text(ctx, result, &text);
// → "sin(2*x)"

cas_string_destroy(text);
cas_expr_destroy(expr);
cas_expr_destroy(result);
cas_context_destroy(ctx);
```

Full API surface: [`include/cas_api.h`](include/cas_api.h)

---

## Build

**Requirements:** CMake ≥ 3.20, Ninja, Clang or GCC with C++20 support, MPFR.

```bash
# Configure (Debug — ASan + UBSan enabled automatically)
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug

# Build
cmake --build build

# Run test suite
ctest --test-dir build --output-on-failure
```

### Release build

```bash
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

### Build with ccache (recommended for development)

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
```

`brew install ccache` on macOS. Typically 5–20× faster on incremental rebuilds.

---

## Testing

```bash
# Full test suite (2 307 tests)
ctest --test-dir build --output-on-failure

# Fast smoke suite only
bash scripts/test_quick.sh

# Property-based tests (rapidcheck)
ctest --test-dir build -R property --output-on-failure

# Milestone gate
./scripts/verify_milestone.sh M1b
```

### Sanitizers

ASan and UBSan are enabled by default in `Debug` and `RelWithDebInfo` builds.

```bash
cmake -S . -B build-san -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCAS_ENABLE_SANITIZERS=ON
cmake --build build-san
ASAN_OPTIONS=detect_leaks=0 ctest --test-dir build-san --output-on-failure
```

### Coverage

```bash
cmake -S . -B build-cov -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DENABLE_COVERAGE=ON \
  -DCAS_ENABLE_SANITIZERS=OFF
cmake --build build-cov
cmake --build build-cov --target coverage
# Report: build-cov/coverage/index.html
```

Requires `lcov` + `genhtml` (`brew install lcov`).

---

## Golden Corpus

The engine is validated against a Maxima 5.49.0 oracle across **1 026 test cases** spanning 11 mathematical areas.

| Area | Pass | Fail | Skip | % non-skip |
|------|------|------|------|-----------|
| diff | 77 | 3 | 0 | 96.2 % |
| gcd | 70 | 0 | 11 | 100.0 % |
| factor | 98 | 1 | 0 | 99.0 % |
| simplify | 110 | 6 | 0 | 94.8 % |
| series | 70 | 2 | 9 | 97.2 % |
| special_fn | 74 | 1 | 5 | 98.7 % |
| limit | 83 | 4 | 12 | 95.4 % |
| solve | 81 | 0 | 0 | 100.0 % |
| matrix | 79 | 0 | 0 | 100.0 % |
| integrate | 96 | 17 | 27 | 85.0 % |
| bronstein | 35 | 17 | 38 | 67.3 % |
| **Total** | **873** | **51** | **102** | **94.5 %** |

> *94.5 % excludes SKIP entries. Total coverage including SKIP: 85.1 %.*  
> *Bronstein (Risch integration corpus) is the primary open area; Hermite cap. 5 is next.*

---

## Benchmarks

```bash
# Run Release benchmarks
bash scripts/benchmark.sh

# Check against pinned baseline (fails if any metric regresses > 10 %)
bash scripts/benchmark.sh --check

# Update baseline after an approved change
bash scripts/benchmark.sh --update-baseline
```

Baseline and policy: `test/benchmarks/baseline_release.txt` / `test/benchmarks/policy.json`.

---

## Code Quality Gates

### Anti-monolith
Every `.cpp` and `.hpp` under `src/` and `include/` must be ≤ 500 lines. Pre-existing violations are whitelisted in `scripts/file_size_whitelist.txt` and tracked in `HARDCODE_LEDGER.md` with mandatory split tickets.

```bash
bash scripts/check_file_size.sh --verbose
```

This check also runs as a CI job on every push (`anti-monolith` job).

### Benchmark regression

```bash
bash scripts/check_bench_regression.sh \
  --current build-bench/current.txt \
  --threshold 10
```

### Mutation testing (optional, weekly CI)

```bash
brew install mull
bash scripts/run_mutation.sh
# Report: mutation-report/index.html
```

Target mutation score: ≥ 70 % per module.

---

## Development Status

| Milestone | Status |
|-----------|--------|
| F0 — Foundation, build infrastructure, sanitizers | ✅ Complete |
| F1 / F1b — Parser, lexer, round-trip, benchmark discipline | ✅ Complete |
| F2 — Canonical simplification, assumptions | ✅ Complete |
| F3 — Polynomial algebra (Gröbner, GCD, factorization, FGLM) | ✅ Complete |
| F4 — Differentiation, integration (Risch), limits (Gruntz) | ✅ Complete |
| F5 — Special functions, series, ODE, transforms | ✅ Complete |
| F6 — Linear algebra (Jordan, Smith, QR, eigenvectors) | ✅ Complete |
| F7.5 — HP Prime G2 parity corpus (94.5 % non-skip) | ✅ Complete (conditional) |
| F8 — C API hardening, Swift bridge integration | 🔜 Next |

### Open research targets (permanent)

- Risch structure theorem — Hermite reduction over deep exp+log towers (Bronstein Ch. 5, B2/B3)
- Galois groups for deg ≥ 5 (Stauduhar / Soubin)
- CAD (Cylindrical Algebraic Decomposition) for quantifier elimination
- Hypergeometric `_pFq` — Wilf–Zeilberger, Pochhammer, Gauss 2F1
- Schönhage–Strassen FFT multiplication for BigInt (n ≥ 4096 limbs)

---

## Scope Rules

- Pure C++ engine. No Swift, no UI, no ObjC in this repository.
- No floating-point arithmetic in the symbolic core.
- Error handling via `Result<T>` (monadic); no exceptions in engine code.
- Tests validate AST structure or mathematical equivalence — never formatter output.
- Every constant requires a mathematical justification or a `CASContext` configuration knob (zero-hardcode policy).

See [`CLAUDE.md`](CLAUDE.md) for the full architectural constitution and [`CONTRIBUTING.md`](CONTRIBUTING.md) for contributor rules.

---

## Project Layout

| Path | Purpose |
|------|---------|
| `include/cas/` | Public C++ API headers |
| `include/cas_api.h` | Public C API (FFI boundary) |
| `src/` | Engine source modules |
| `test/` | Unit, property, and golden tests |
| `scripts/` | CI and developer tooling |
| `CLAUDE.md` | Architectural constitution |
| `CAS_TASKS.md` | Milestone tracker and gap inventory |
| `HARDCODE_LEDGER.md` | Hardcode audit ledger |
| `.APROJECT_REFERENCES/` | Module specs and architecture docs |
