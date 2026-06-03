#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "summation_internal.hpp"
#include <cstdint>
#include <optional>

namespace cas::calculus {

[[nodiscard]] static bool is_one(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value == BigInt(1);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator == BigInt(1) && rational->denominator == BigInt(1);
    }
    return false;
}

[[nodiscard]] static bool is_positive_infinity(ExprPtr expr) {
    const auto* constant = expr_cast<Constant>(expr);
    return constant != nullptr && constant->value == MathConstant::Infinity;
}

[[nodiscard]] static BigInt factorial(unsigned int value) {
    BigInt result(1);
    for (unsigned int factor = 2U; factor <= value; ++factor) {
        result *= BigInt(static_cast<std::int64_t>(factor));
    }
    return result;
}

[[nodiscard]] static BigInt pow_bigint_nonnegative(BigInt base, unsigned int exponent) {
    BigInt result(1);
    while (exponent > 0U) {
        if ((exponent % 2U) == 1U) {
            result *= base;
        }
        exponent /= 2U;
        if (exponent > 0U) {
            base *= base;
        }
    }
    return result;
}

[[nodiscard]] static std::optional<unsigned int> reciprocal_power_exponent(
    ExprPtr term, const Symbol& var) {
    const auto* division = expr_cast<Binary>(term);
    if (division == nullptr || division->op != BinaryOp::Div || !is_one(division->left)) {
        return std::nullopt;
    }
    const auto* power = expr_cast<Binary>(division->right);
    if (power == nullptr || power->op != BinaryOp::Pow) {
        return std::nullopt;
    }
    const auto* symbol = expr_cast<Symbol>(power->left);
    if (symbol == nullptr || symbol->name != var.name) {
        return std::nullopt;
    }
    return positive_integer_u32(power->right);
}

[[nodiscard]] static Result<ExprPtr> zeta_even_value(
    unsigned int exponent, symbolic::CASContext& ctx) {
    if (exponent == 0U || (exponent % 2U) != 0U) {
        return make_unimplemented<ExprPtr>(
            "calculus", "zeta_even_value",
            "p-series with odd or zero exponent",
            cas::error::reason_codes::SUMMATION_GENERAL,
            "Implement zeta at odd arguments (open problem) or return unevaluated Zeta(n) form",
            "F0.8");
    }
    const unsigned int m = exponent / 2U;
    const auto bernoulli = cas::numtheory::bernoulli_numbers(exponent);
    Rational coefficient = bernoulli[exponent]
        * Rational(pow_bigint_nonnegative(BigInt(2), exponent - 1U), factorial(exponent));
    if ((m % 2U) == 0U) {
        coefficient = -coefficient;
    }
    AstArena& arena = ctx.arena();
    ExprPtr pi_power = pow_expr(arena, arena.make<Constant>(MathConstant::Pi), exponent);
    ExprPtr coeff    = rational_expr(arena, coefficient);
    if (const auto* integer = expr_cast<IntegerLit>(coeff);
        integer != nullptr && integer->value == BigInt(1)) {
        return ok(pi_power);
    }
    return ctx.simplify(arena.make<Binary>(BinaryOp::Mul, coeff, pi_power));
}

[[nodiscard]] Result<ExprPtr> symbolic_sum(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    if (is_positive_infinity(upper)) {
        if (const auto* l = expr_cast<IntegerLit>(lower);
            l != nullptr && l->value == BigInt(1)) {
            if (auto exponent = reciprocal_power_exponent(term, var);
                exponent.has_value()) {
                return zeta_even_value(*exponent, ctx);
            }
        }
        // Infinite bounds without a recognised p-series form — fall through.
    } else {
        // Definite finite sum: try Gosper → single-atom polygamma →
        // full Abramov (partial fractions + Q-linear + Q-quadratic RootOf paths).
        auto gosper_res = try_gosper_definite(term, var, lower, upper, ctx);
        if (gosper_res.is_ok()) return gosper_res;
        auto poly_res = try_polygamma_definite(term, var, lower, upper, ctx);
        if (poly_res.is_ok()) return poly_res;
        auto abramov_res = try_abramov_definite(term, var, lower, upper, ctx);
        if (abramov_res.is_ok()) return abramov_res;
        auto zeil_res = try_zeilberger_definite(term, var, lower, upper, ctx);
        if (zeil_res.is_ok()) return zeil_res;
    }

    return make_unimplemented<ExprPtr>(
        "calculus", "sum_closed_form",
        "general summand not in closed-form table",
        cas::error::reason_codes::SUMMATION_GENERAL,
        "Implement Petkovšek-WZ / Zeilberger creative telescoping or extend "
        "Abramov for terms outside the Gosper/polygamma/quadratic-RootOf class",
        "F0.8");
}

Result<ExprPtr> sum(
    ExprPtr expr, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    return symbolic_sum(expr, var, lower, upper, ctx);
}

} // namespace cas::calculus
