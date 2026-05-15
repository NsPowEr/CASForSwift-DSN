#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include <vector>

namespace cas::calculus {

enum class OdeType {
    Unknown,
    Separable,
    Linear1stOrder,
    Bernoulli,
    Exact,
    Linear2ndOrderConstantCoeff,
    Linear2ndOrderRationalCoeff,
    LinearNthOrderConstantCoeff
};

struct OdeClassification {
    OdeType type;
    ExprPtr equation; // y' = ... o f(x,y,y') = 0
    Symbol y;
    Symbol x;
    std::vector<ExprPtr> components; // P(x), Q(x), etc.

    OdeClassification(OdeType t, ExprPtr eq, Symbol sy, Symbol sx)
        : type(t), equation(eq), y(std::move(sy)), x(std::move(sx)) {}
};

[[nodiscard]] Result<OdeClassification> classify_ode(ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> solve_ode_1st_order(const OdeClassification& classification, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> solve_ode_advanced(const OdeClassification& classification, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> solve_ode(ExprPtr equation, const Symbol& y, const Symbol& x, symbolic::CASContext& ctx);

// Frobenius series solution for a homogeneous linear 2nd-order ODE
//   a_2(x) y'' + a_1(x) y' + a_0(x) y = 0
// expanded around the regular singular point x = 0.
//
// Computes the indicial equation
//   r(r-1) + p_0 r + q_0 = 0,    p_0 = lim x p(x),   q_0 = lim x^2 q(x),
// (with p = a_1/a_2, q = a_0/a_2), and constructs `num_terms` coefficients
// of each Frobenius series solution.
//
// Currently supports the regular case where the two indicial roots differ
// by a non-integer.  When the difference is a non-negative integer and a
// resonance occurs in the recurrence (forcing a logarithmic term), an
// Unimplemented error is returned with diagnostic context.
//
// The result is the general solution  C_1 * y_1(x) + C_2 * y_2(x)
// using freshly named constants (`_C1_`, `_C2_`).
[[nodiscard]] Result<ExprPtr> solve_ode_frobenius_at_zero(
    ExprPtr a_2,
    ExprPtr a_1,
    ExprPtr a_0,
    const Symbol& x,
    unsigned int num_terms,
    symbolic::CASContext& ctx);

} // namespace cas::calculus
