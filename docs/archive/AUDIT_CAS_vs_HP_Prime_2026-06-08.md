# AUDIT — CAS Engine vs HP Prime G2 — 2026-06-08

> Stato di chiusura per fase del piano `PLAN_HP_PRIME_PARITY.md`.
> Quick suite: **2111/2111 PASS** — zero fail residui.
> Disabilitati documentati: 2 (stress test long-running notturni).

---

## Riepilogo per fase

| Fase | Titolo | Stato | Note |
|---|---|---|---|
| F0 | Sanitizzazione tracker + property test + CI | ✅ Chiusa | F0.1-F0.8 tutti integrati (BigInt CI, rapidcheck, sanitizers, anti-monolith 500-line gate) |
| F1 | STRATO L0 Foundation (BigInt, Rational, Arena, Simplifier, Pattern matcher, **Complex F1.6**) | ✅ Chiusa | F1.6 ComplexLit canonical pipeline finalizzato; ln(a+bi) generale via Abs/Arg |
| F2 | STRATO L1 Univariate polys (GCD, Resultant, Factor Z[x]/Fp[x], Partial fractions, Cyclotomic) | ✅ Chiusa | Half-GCD + Berlekamp + Cantor-Zassenhaus + Van Hoeij LLL + Trager Q(α); F2-GATE 100 poly random Z[x] in 2.98s |
| F3 | STRATO L2 Multivariate + Algebraic numbers tower | ✅ Chiusa | Brown GCD, Zippel sparse interp (n-variate), F4 Groebner, primitive element Trager 3+ livelli, Galois deg-5+ |
| F4 | STRATO L2 Linear algebra | ✅ Chiusa | LU PivotScore, MGS QR, Cholesky LDL^T, Jordan RootOf-aware, Vandermonde closed-form, fresh-symbol pervasivo |
| F5 | STRATO L2 Calculus | ✅ Chiusa | Risch cap.7-9 Bronstein, Gruntz dynamic, ODE 1st-order (Riccati/Bernoulli/Clairaut/d'Alembert), Frobenius log+resonance, Pade non-Q, Puiseux, Residue deg ≥5, Zeilberger creative telescoping + Pochhammer fast-path, Abramov, Fourier/Mellin/Z/Laplace transforms, special fns (Hypergeometric, Bessel/Erf/Zeta, Jacobi P, Elliptic K/E/Π/F) |
| F6 | STRATO L3 Numerica + Complex + Units | ✅ Chiusa | fsolve config, Sturm-based 1-var inequality, RootOf Sturm-isolated, SI dimensional check su Sum, ln/exp branch cut, Gaussian Z[i] factorization |
| F7 | Integration finale + Plotting + Statistics + Parity HP Prime | 🟡 In progress | F7.1 ParametricSampler, F7.2 Normal + Linear regression, F7.3 Gauss-Kronrod + Lagrange — F7.4 corpus HP Prime ≥95% pending |

---

## Debiti aperti documentati in `HARDCODE_LEDGER.md`

| ID | Categoria | Bloccante? |
|---|---|---|
| F5.7-ZEIL-HIGHER-ORDER | Cat 3 algoritmo incompleto (Petkovšek higher-order Hyper) | No — Zeilberger base produce le recurrence; solver superiore richiesto solo per casi research |
| F5.7-B6BIS-QUADRATIC-M-GT-1 | Cat 4 bail-out tipo (polygamma alto ordine su Q-irriducibili m>1) | No — frequenza rara; partial_fractions tipicamente produce m=1 |

### Risolti recentemente
- HC-F16-LN-COMPLEX-FULL → CHIUSA via dispatch ln(a+bi) = ln|z| + i·arg(z) (commit 1a36321 + temp commit)
- HC-F16-TRAGER-QI → CHIUSA via expand_expr_impl ComplexLit leaf (commit a797c02)
- HC-F57-ZEIL-GAMMA-RATIO → CHIUSA via Pochhammer fast-path (commit 1a36321)

---

## Conformità Exit Gate per fase

### F1 Foundation
- ✅ BigInt benchmark target onesto (Toom-3, Knuth-D, Lehmer-GCD, FFT >8192 limbs)
- ✅ Property test foundation passano
- ✅ Assume() copertura algebrica ≥50 implicazioni canoniche (AssumptionsCanonical50 suite verde)
- ✅ Hashconsing arena bump allocator operativo
- ✅ Q[i] arithmetic completo + Gaussian primes factor
- ✅ Public API congelata (`include/cas/*.hpp`)

### F2 Univariate L1
- ✅ Bench: 100 poly random Z[x] in 2.98s (target <30s) → F2GateBenchmark
- ✅ Half-GCD + Berlekamp + Cantor-Zassenhaus + Van Hoeij LLL operativi
- ✅ 0 Unimplemented in algebra univariata Z/Q

### F3 Multivariate L2
- ✅ Brown modular GCD + Zippel n-variate operativi
- ✅ Galois deg ≥5 (S5/A5/D5/F20/C5)
- ✅ Trager su tower ≥3 livelli

### F4 Linear Algebra
- ✅ LU partial pivoting con PivotScore vero
- ✅ QR Modified-Gram-Schmidt simbolico stabile
- ✅ Cholesky LDL^T per matrici simmetriche
- ✅ Jordan via null_space_over_extension routing RootOf
- ✅ Smith Z + Hermite + Vandermonde closed-form
- ✅ Fresh-symbol pervasivo (no più "_lambda_", "C1", "t1" hardcoded)

### F5 Calculus
- ✅ Risch Bronstein cap.7-9 (parametric problems, Risch DE generale, structure theorem)
- ✅ Gruntz dynamic growth rank
- ✅ ODE classifier completo + solver per ogni tipo
- ✅ Hadamard finite part poles ≥2
- ✅ Puiseux + Pade non-Q
- ✅ Residue deg ≥5 via Bronstein
- ✅ Gosper + Zeilberger + Abramov + WZ pair
- ✅ Fourier/Mellin/Z + Laplace residue-based inverse
- ✅ Hypergeometric pFq + Bessel/Erf/Zeta + Jacobi P + Elliptic K/E/Π/F

### F6 Numerica + Complex + Units
- ✅ fsolve real-filter + tolerance configurable
- ✅ Solve_inequality 1-var via Sturm (no CAD)
- ✅ RootOf evaluator deterministic via Sturm bracketing
- ✅ SI unit dimensional check su Sum (dim mismatch → Undefined)

### F7 Integration + Plotting + Statistics — work-in-progress
- ✅ F7.1-T1 ParametricSampler adaptive curvature-based
- ✅ F7.2-T1 Normal pdf/cdf/quantile
- ✅ F7.2-T2 Linear regression OLS univariate
- ✅ F7.3-T1 Gauss-Kronrod 15/7 adaptive quadrature
- ✅ F7.3-T2 Lagrange interpolation barycentric
- 🟡 F7.1 Implicit/Contour/Vector field plot — TODO
- 🟡 F7.2 Binomial/Poisson/χ²/t/F distributions — TODO
- 🟡 F7.2 Hypothesis testing + Linear regression multivariata — TODO
- 🟡 F7.3 Spline + Hermite interpolation — TODO
- 🟡 F7.3 RKF45 adaptive ODE — TODO
- 🟡 F7.4 Golden corpus HP Prime ≥95% — TODO (richiede corpus 2000+ entry)

---

## Architettura corrente

- C++20 stretto, BigInt limbs + Toom-3 + FFT, Rational, AstArena bump allocator, structural sharing
- Simplifier modulare: `simplify_arithmetic_chain.cpp` 440 LOC + `_sum` + `_liketerm` + `_gamma` + `_sqrt` (anti-monolith 500-line gate rispettato)
- Pipeline ComplexLit: parse → simplify (Sum/Product/Unary canonicalize coefficients via ComplexRational) → format (canonical "I"/"-I"/"a ± b * I") → expand/substitute → mathematically_equal con bridge cross-form
- Statistics package: `cas::statistics` namespace (`include/cas/statistics.hpp`) — Normal + LinearRegression
- Numeric package: GK15/G7 quadrature + Lagrange + Sturm + Aberth root isolator + Lipschitz fsolve + Adaptive/Parametric sampler + RK4 ODE
- Algebra: solve_polynomial, solve_inequality (Sturm-based), factor_polynomial (Berlekamp/CZ/LLL/Trager), Groebner F4+F5

---

## Prossimi step naturali verso F7.4 parity HP Prime ≥95%

1. **F7.2** completare distribuzioni Binomial/Poisson/χ²/t/F + pdf/cdf/quantile
2. **F7.3** Spline cubica + Hermite interpolation
3. **F7.3** RKF45 adaptive ODE solver
4. **F7.1** Implicit plot via marching squares + Contour + Vector field
5. **F7.4** Costruire golden corpus HP Prime G2 con 2000+ entry (referenza Maxima + HP Prime emulator)
6. **F7.4** Benchmark performance: target ≥80% Giac-Xcas su benchmark identici

---

## Commit chain F7 (questa sessione)

| Commit | Task | Test |
|---|---|---|
| `fb6215e` | F7.3-T1 Gauss-Kronrod 15/7 | 3/3 |
| `bb00c9f` | F7.3-T2 Lagrange interpolation | 4/4 |
| `dfe5ac8` | F7.2-T1 Normal pdf/cdf/quantile | 4/4 |
| `59b8bc5` | F7.2-T2 Linear regression OLS | 4/4 |
| `af8451e` | F7.1-T1 ParametricSampler | 3/3 |

Totale aggiunte test session F7: **18 test, 18 pass**.
