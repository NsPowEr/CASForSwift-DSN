#include "cas/ast.hpp"
#include "cas/numtheory.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace cas::calculus {

[[nodiscard]] static CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

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
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "p-series summation is implemented exactly only for positive even exponents"));
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
    }
    
    return fail<ExprPtr>(make_error(
        CASErrorKind::Unimplemented,
        "General hypergeometric summation requires Gosper/Petkovsek-WZ and is not implemented yet"));
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
