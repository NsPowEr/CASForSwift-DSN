// Internal header shared between summation.cpp and summation_abramov.cpp.
// Not part of the public CAS API.
#pragma once
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include <cstdint>
#include <limits>
#include <optional>

namespace cas::calculus {

[[nodiscard]] inline std::optional<unsigned int> positive_integer_u32(ExprPtr expr) {
    const auto* integer = expr_cast<IntegerLit>(expr);
    if (integer == nullptr || integer->value <= BigInt(0)) {
        return std::nullopt;
    }
    if (integer->value.bit_length() > std::numeric_limits<unsigned int>::digits) {
        return std::nullopt;
    }
    return static_cast<unsigned int>(integer->value.to_u64());
}

[[nodiscard]] inline ExprPtr rational_expr(AstArena& arena, const Rational& value) {
    if (value.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] inline ExprPtr pow_expr(AstArena& arena, ExprPtr base, unsigned int exponent) {
    if (exponent == 0U) return arena.make<IntegerLit>(BigInt(1));
    if (exponent == 1U) return base;
    return arena.make<Binary>(BinaryOp::Pow, base,
        arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(exponent))));
}

// Implemented in summation_abramov.cpp:

[[nodiscard]] Result<ExprPtr> try_gosper_definite(
    ExprPtr term, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> try_polygamma_definite(
    ExprPtr term, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> try_abramov_definite(
    ExprPtr term, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx);

// Implemented in summation_zeilberger_driver.cpp:

[[nodiscard]] Result<ExprPtr> try_zeilberger_definite(
    ExprPtr term, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx);

} // namespace cas::calculus
