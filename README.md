# REAL CAS ENGINE

> An industrial-grade symbolic computer algebra system in modern C++20.

[![Build](https://img.shields.io/badge/build-passing-brightgreen)](#build)
[![Tests](https://img.shields.io/badge/tests-2858%20passing-brightgreen)](#testing)
[![Maxima Golden](https://img.shields.io/badge/golden%20vs%20Maxima-99.6%25%20non--skip-blue)](#golden-corpus)
[![Giac Parity](https://img.shields.io/badge/second%20oracle-Giac%202.0.0-lightgrey)](#golden-corpus)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue)](#requirements)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

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
- **Integration (indefinite)** — Risch algorithm (Bronstein Ch. 5–9): Hermite reduction, Trager/LRT, log-derivative recognizer, Risch structure theorem (algebraic independence of log/exp extensions), IBP, Weierstrass substitution `t = tan(x/2)`, u-substitution framework; non-elementary closed forms (Ei, Si, Ci, Shi, Chi, li, Li₂, erfi) via Meijer-G bridge (Slater expansion) once Liouville's theorem proves no elementary antiderivative exists
- **Integration (definite)** — Fundamental Theorem of Calculus, singularity detection (rational poles, algebraic, transcendental), Cauchy principal value, Hadamard finite part (poles of order m ≥ 2), residue theorem (quadratic and biquadratic irreducible factors), multiple integrals via Fubini (rectangular domains), orthogonal polynomial patterns (Legendre, Hermite-H/He, Chebyshev T/U), Mellin-convolution definite integrator
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
- Galois groups — exact deg ≤ 7 (Newton power-sum resolvents), certified Stauduhar descent deg 8–10 (p-adic root lifting, structural naming from certified invariants only, no lookup tables); deg ≥ 11 out of scope
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

### Statistics
- Distributions — Normal (symbolic pdf/cdf via Erf, Newton quantile), Binomial, Poisson, Chi-squared, Student-t, Fisher-Snedecor F (pdf/cdf/quantile)
- Regression — univariate OLS and multivariate OLS (normal equations via Gauss-Jordan with partial pivoting), R², residual sum of squares
- Hypothesis testing — one/two-sample z-test and t-test (Welch), chi-squared goodness-of-fit, F-test for variance ratio

---

## Architecture

```
include/cas/          Public C++ headers (engine API)
include/cas_api.h     Public C API (Swift / FFI boundary)
src/
  foundation/         BigInt, Rational, BigFloat
  ast/                AST nodes, AstArena, hash-consing
  lexer/ parser/      Expression parsing, token stream
  symbolic/           Simplifier, assumptions, side-conditions, Meijer-G
  algebra/            Polynomials, GCD, factorization, Gröbner, Galois groups
  calculus/           Diff, integrate (Risch), limits (Gruntz), series, ODE
  rewrite/            Discrimination-net builtin rewrite rules
  linalg/             Matrix algorithms
  statistics/         Distributions, OLS regression, hypothesis tests
  numeric/            MPFR evaluation, interval arithmetic
  numtheory/          Bernoulli, Euler phi, primality
  formatter/          Plain text and LaTeX output
  c_api/              C API bridge (cas_api.h implementation)
test/                 GoogleTest unit + property-based tests (2 858 tests)
  golden/             HP Prime G2 golden corpus (1 026 entries, dual oracle: Maxima + Giac)
scripts/              CI helpers: benchmark, anti-monolith, golden/parity measurement
```

**Key invariants enforced at all times:**

| Rule | Detail |
|------|--------|
| No floating-point in symbolic core | `BigInt` / `Rational` only; `double` forbidden |
| Structural sharing | Functions return the original `ExprPtr` when no change occurs |
| Arena allocation | All AST nodes allocated via `AstArena::make`; `new` / `make_unique` forbidden |
| Zero hardcode | Every constant must have a mathematical justification or a `CASContext` knob |
| Zero warnings | `-Wall -Wextra -Wpedantic -Werror` on all targets |
| 500 LOC gate | Each `.cpp` / `.hpp` ≤ 500 lines; enforced pre-commit (0 files pending split) |

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
# Full test suite (2 858 tests)
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

The engine is validated against two independent, unmodified reference implementations (fork/exec only — never linked, patched, or embedded; see [Reference Oracles](#reference-oracles)) across **1 026 test cases** spanning 11 mathematical areas (diff, gcd, factor, simplify, series, special_fn, limit, solve, matrix, integrate, bronstein).

**Primary oracle — Maxima 5.49.0** (ratchet-gated, `scripts/golden_baseline.txt`, tolerance 0):

| Metric | Value |
|--------|-------|
| Pass | 980 |
| Fail | 4 |
| Undecided (skip + over-budget) | 42 |
| Non-skip pass rate | 99.6 % |
| Full coverage (incl. undecided) | 95.5 % |

The ratchet is a floor/ceiling gate (`scripts/check_golden_ratchet.sh`): pass cannot drop, fail/undecided cannot rise, without an explicit `--update-baseline` after human review. **Bronstein (Risch integration corpus)** remains the area with the highest residual fail/skip share.

**Second oracle — Giac 2.0.0** (parity target, not a correctness gate): `PARITY_GIAC.md` is a regenerable scoreboard (`giac-parity-scan` skill / `scripts/giac_parity_report.py`) comparing Giac's own closed-form coverage against the CAS pass rate per area, plus a per-entry cross-diff dump for areas with the largest gap. It is a measurement artifact, not a tracker — every gap it surfaces becomes a task in `TASKLIST_MASTER.md`.

```bash
# Regenerate the Maxima ratchet report
bash scripts/run_golden_measurement.sh
bash scripts/check_golden_ratchet.sh --report <report.json>

# Regenerate the Giac parity scoreboard
bash scripts/run_golden_giac.sh
python3 scripts/giac_parity_report.py
```

### Reference Oracles

| Oracle | Version | License | Usage |
|--------|---------|---------|-------|
| Maxima | 5.49.0 | GPL-2.0-only | `maxima --very-quiet --batch-string=...`, textual output parsing only |
| Giac | 2.0.0 | GPL-3.0-or-later | `icas` via stdin, textual output parsing only |

Both are pinned by SHA-256 manifest (`scripts/verify_maxima_integrity.sh`, `scripts/giac_integrity.sh`). Neither is patched, recompiled, or linked into the engine — modifying either would create a copyleft derivative work and invalidate them as independent oracles.

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
Every `.cpp` and `.hpp` under `src/` and `include/` must be ≤ 500 lines. `scripts/file_size_whitelist.txt` currently holds zero active entries (all past violations split below the limit); any new one requires a split ticket in `HARDCODE_LEDGER.md` before it can be added.

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
| F7.2 — Statistics (distributions, OLS, hypothesis tests) | ✅ Complete |
| F7.5 — HP Prime G2 parity corpus (99.6 % non-skip vs Maxima) | ✅ Complete (conditional) |
| F7.5.G — Giac 2.0.0 second-oracle parity infra (per-entry cross-diff, `PARITY_GIAC.md`) | ✅ Active (continuous measurement) |
| F8 — C API hardening, Swift bridge integration | 🔜 Next |

### Open research targets (permanent)

- Bronstein corpus residual — deepest Risch parametric cases still fail/skip against Maxima (primary open area, see [Golden Corpus](#golden-corpus))
- Galois groups for deg ≥ 11 (beyond certified Stauduhar descent, currently exact through deg 10)
- CAD (Cylindrical Algebraic Decomposition) for quantifier elimination
- Full Zeilberger / WZ certificate algorithm (creative telescoping beyond current partial summation support)
- Schönhage–Strassen FFT multiplication for BigInt (n ≥ 4096 limbs)

---

## Scope Rules

- Pure C++ engine. No Swift, no UI, no ObjC in this repository.
- No floating-point arithmetic in the symbolic core.
- Error handling via `Result<T>` (monadic); no exceptions in engine code.
- Tests validate AST structure or mathematical equivalence — never formatter output.
- Every constant requires a mathematical justification or a `CASContext` configuration knob (zero-hardcode policy).

See [`CLAUDE.md`](CLAUDE.md) for the full architectural constitution and [`CONTRIBUTING.md`](CONTRIBUTING.md) for contributor rules. Licensed under [MIT](LICENSE).

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
| `TASKLIST_MASTER.md` | Single source of truth: open/closed task tracker |
| `PARITY_GIAC.md` | Regenerable Giac parity scoreboard (measurement, not a tracker) |
| `HARDCODE_LEDGER.md` | Hardcode audit ledger |
| `docs/rules/` | Detailed operating protocols (perf root-cause, verification, oracles) |
| `docs/archive/` | Superseded trackers (STATE, TODO*, PLAN_*, CAS_TASKS, …) |
| `.APROJECT_REFERENCES/` | Module specs and architecture docs |
