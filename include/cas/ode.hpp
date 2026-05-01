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
    Linear2ndOrderRationalCoeff
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

} // namespace cas::calculus
