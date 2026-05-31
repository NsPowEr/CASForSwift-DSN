// polynomial_gcd_modular.cpp — Dispatch layer for modular/probabilistic multivariate GCD.
//
// Entry points:
//   gcd_modular(P, Q)       — Brown's evaluation/interpolation/CRT (F3.1).
//   gcd_probabilistic(P, Q) — Las Vegas dispatch: GCDHEU → eval-interp → Brown's.
//
// The actual Brown algorithm lives in polynomial_gcd_brown.cpp.

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/numtheory.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace cas::algebra {

// gcd_modular: Brown's modular GCD, F3.1.
// Pipeline:
//   1. Trivial cases (one operand zero, constants, univariate → handled by Brown's recursion).
//   2. Try GCDHEU heuristic first (very fast for dense polynomials with small coefficients).
//   3. If GCDHEU fails, use Brown's recursive evaluation/interpolation.
//   4. If Brown's fails (unlucky points or budget), fall back to eval-interp (existing robust path).
// Every error path returns an explicit Unimplemented — never a wrong result.
Result<MultivariatePolynomial> gcd_modular(
        const MultivariatePolynomial& P,
        const MultivariatePolynomial& Q) {
    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);

    // Build a minimal CASContext for budget parameters.
    // Brown's needs ctx for gcd_error_probability() and budget knobs.
    // We create a default context here; callers who want custom budgets
    // should use gcd_brown() directly with their own ctx.
    symbolic::CASContext ctx;

    auto vars = P.variables();
    if (vars.empty()) {
        // Both are constants.
        BigInt g(0);
        for (const auto& t : P.terms()) g = (g.is_zero()) ? t.coefficient.abs() : g;
        for (const auto& t : Q.terms()) {
            BigInt tmp = g; BigInt qt = t.coefficient.abs();
            while (!qt.is_zero()) { BigInt r = tmp % qt; tmp = qt; qt = r; }
            g = tmp;
        }
        return ok(MultivariatePolynomial({{ .coefficient = g, .factors = {} }}));
    }

    // Step 1: Try GCDHEU heuristic (fast path).
    auto heu = gcd_heuristic(P, Q);
    if (heu.is_ok()) return heu;

    // Step 2: Brown's recursive evaluation/interpolation.
    auto brown = gcd_brown(P, Q, ctx);
    if (brown.is_ok()) return brown;

    // Step 3: Fall back to the eval-interp multivariate engine (robust certified path).
    auto interp = gcd_multivariate_eval_interp(P, Q, ctx);
    if (interp.is_ok()) return interp;

    // All paths exhausted — report the Brown's failure (most informative).
    return brown;
}

// gcd_probabilistic: Las Vegas dispatch.
// Guaranteed correct when it returns ok(); may return Unimplemented.
Result<MultivariatePolynomial> gcd_probabilistic(
        const MultivariatePolynomial& P,
        const MultivariatePolynomial& Q) {
    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);
    return gcd_modular(P, Q);
}

}  // namespace cas::algebra
