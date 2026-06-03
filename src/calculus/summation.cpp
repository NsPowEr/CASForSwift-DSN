#include "cas/ast.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "../symbolic/summation_gosper.hpp"
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

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

[[nodiscard]] static std::optional<unsigned int> positive_integer_u32(ExprPtr expr) {
    const auto* integer = expr_cast<IntegerLit>(expr);
    if (integer == nullptr || integer->value <= BigInt(0)) {
        return std::nullopt;
    }
    if (integer->value.bit_length() > std::numeric_limits<unsigned int>::digits) {
        return std::nullopt;
    }
    return static_cast<unsigned int>(integer->value.to_u64());
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

[[nodiscard]] static ExprPtr rational_expr(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] static ExprPtr pow_expr(AstArena& arena, ExprPtr base, unsigned int exponent) {
    if (exponent == 0U) return arena.make<IntegerLit>(BigInt(1));
    if (exponent == 1U) return base;
    return arena.make<Binary>(BinaryOp::Pow, base, arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(exponent))));
}

[[nodiscard]] static std::optional<unsigned int> reciprocal_power_exponent(ExprPtr term, const Symbol& var) {
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

[[nodiscard]] static Result<ExprPtr> zeta_even_value(unsigned int exponent, symbolic::CASContext& ctx) {
    if (exponent == 0U || (exponent % 2U) != 0U) {
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "zeta_even_value",
            "p-series with odd or zero exponent",
            cas::error::reason_codes::SUMMATION_GENERAL,
            "Implement zeta at odd arguments (open problem) or return unevaluated Zeta(n) form",
            "F0.8");
    }

    const unsigned int m = exponent / 2U;
    const auto bernoulli = cas::numtheory::bernoulli_numbers(exponent);
    Rational coefficient = bernoulli[exponent] * Rational(pow_bigint_nonnegative(BigInt(2), exponent - 1U), factorial(exponent));
    if ((m % 2U) == 0U) {
        coefficient = -coefficient;
    }

    AstArena& arena = ctx.arena();
    ExprPtr pi_power = pow_expr(arena, arena.make<Constant>(MathConstant::Pi), exponent);
    ExprPtr coeff = rational_expr(arena, coefficient);
    if (const auto* integer = expr_cast<IntegerLit>(coeff); integer != nullptr && integer->value == BigInt(1)) {
        return ok(pi_power);
    }
    return ctx.simplify(arena.make<Binary>(BinaryOp::Mul, coeff, pi_power));
}

// F5.7 sub-block 0 — definite hypergeometric summation via Gosper.
//
// gosper_sum(t, k) returns (when it succeeds) an antidifference S(k) such
// that S(k+1) − S(k) = t(k).  The Newton-Leibniz analogue for the finite
// calculus then gives  Σ_{k=a}^{b} t(k) = S(b+1) − S(a).
//
// Indeterminate bounds (free symbols, RootOf, …) are passed verbatim into
// the substitution machinery; the simplifier folds them when possible and
// leaves them symbolic otherwise.
[[nodiscard]] static Result<ExprPtr> try_gosper_definite(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff = symbolic::gosper_sum(term, var, ctx);
    if (antidiff.is_error()) return fail<ExprPtr>(antidiff.error());
    if (!antidiff.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Gosper: term is not Gosper-summable",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff.value().value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub, S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

[[nodiscard]] Result<ExprPtr> symbolic_sum(
    ExprPtr term,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    if (is_positive_infinity(upper)) {
        if (const auto* l = expr_cast<IntegerLit>(lower); l != nullptr && l->value == BigInt(1)) {
            if (auto exponent = reciprocal_power_exponent(term, var); exponent.has_value()) {
                return zeta_even_value(*exponent, ctx);
            }
        }
        // Infinite Gosper-summable series:  if S(k) → 0 as k → ∞ we cannot
        // decide this without an asymptotic analysis, so fall through to
        // diagnostic Unimplemented below rather than risk silent error.
    } else {
        // Definite finite-bound sum — try Gosper.
        auto gosper_res = try_gosper_definite(term, var, lower, upper, ctx);
        if (gosper_res.is_ok()) return gosper_res;
        // Non-Gosper-summable: fall through to diagnostic.
    }

    return make_unimplemented<ExprPtr>(
        "calculus", "sum_closed_form",
        "general summand not in closed-form table",
        cas::error::reason_codes::SUMMATION_GENERAL,
        "Implement Petkovšek-WZ / Zeilberger creative telescoping or Abramov "
        "rational summation for terms outside the Gosper hypergeometric class",
        "F0.8");
}

Result<ExprPtr> sum(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    return symbolic_sum(expr, var, lower, upper, ctx);
}

} // namespace cas::calculus
