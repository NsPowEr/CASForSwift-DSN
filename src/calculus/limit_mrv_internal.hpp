#pragma once
// limit_mrv_internal.hpp — private helpers shared across limit_mrv_*.cpp units.
// NOT part of the public CAS API; do not include from include/cas/.

#include "calculus_internal.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"
#include "../symbolic/simplify_impl.hpp"
#include "cas/algebra.hpp"

#include <iostream>
#include <algorithm>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::calculus {

// ---------------------------------------------------------------------------
// LeadingPower — result of extracting the dominant term c·w^p from an
// expression rewritten in terms of the MRV substitute w → 0+.
// ---------------------------------------------------------------------------
struct LeadingPower {
    long long power{};
    ExprPtr coefficient{};
};

// ---------------------------------------------------------------------------
// ExponentialTerm — single term in an exponential-product decomposition:
//   exp(exponent)^coefficient  (coefficient is an integer multiplier).
// ---------------------------------------------------------------------------
struct ExponentialTerm {
    ExprPtr exponent{};
    long long coefficient{};
};

// ---------------------------------------------------------------------------
// integer_value — extract a long long from an IntegerLit or unit-denominator
// RationalLit; returns nullopt when the value does not fit or is not exact.
// ---------------------------------------------------------------------------
[[nodiscard]] inline std::optional<long long> integer_value(ExprPtr expr) {
    constexpr std::uint64_t kMaxI64 =
        static_cast<std::uint64_t>(std::numeric_limits<long long>::max());
    constexpr std::uint64_t kMinI64Magnitude = kMaxI64 + 1ULL;

    auto bigint_to_i64 = [&](const BigInt& value) -> std::optional<long long> {
        const std::uint64_t magnitude = value.abs().to_u64();
        if (!value.is_negative()) {
            if (magnitude > kMaxI64) return std::nullopt;
            return static_cast<long long>(magnitude);
        }
        if (magnitude > kMinI64Magnitude) return std::nullopt;
        if (magnitude == kMinI64Magnitude) return std::numeric_limits<long long>::min();
        return -static_cast<long long>(magnitude);
    };

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return bigint_to_i64(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            return bigint_to_i64(rational->numerator);
        }
    }
    return std::nullopt;
}

// ---------------------------------------------------------------------------
// is_exact_zero / is_exact_nonzero — fast structural checks used in leading
// power cancellation detection.
// ---------------------------------------------------------------------------
[[nodiscard]] inline bool is_exact_zero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return rational->numerator.is_zero();
    }
    symbolic::detail::LiteralRational literal;
    auto exact = symbolic::detail::try_get_exact_rational(expr, literal);
    return exact.is_ok() && exact.value() && literal.value.numerator().is_zero();
}

[[nodiscard]] inline bool is_exact_nonzero(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return !integer->value.is_zero();
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return !rational->numerator.is_zero();
    }
    symbolic::detail::LiteralRational literal;
    auto exact = symbolic::detail::try_get_exact_rational(expr, literal);
    return exact.is_ok() && exact.value() && !literal.value.numerator().is_zero();
}

// ---------------------------------------------------------------------------
// simplify_binary / simplify_unary_neg — thin wrappers that forward to
// CASContext::simplify and convert errors to nullopt.
// ---------------------------------------------------------------------------
[[nodiscard]] inline std::optional<ExprPtr> simplify_binary(
    BinaryOp op,
    ExprPtr lhs,
    ExprPtr rhs,
    symbolic::CASContext& ctx) {
    auto expr = ctx.arena().make<Binary>(op, lhs, rhs);
    auto simplified = ctx.simplify(expr);
    if (simplified.is_error()) {
        std::cerr << "simplify_binary error: " << simplified.error().message
                  << " for op " << static_cast<int>(op) << "\n";
        return std::nullopt;
    }
    return simplified.value();
}

[[nodiscard]] inline std::optional<ExprPtr> simplify_unary_neg(
    ExprPtr operand,
    symbolic::CASContext& ctx) {
    auto simplified = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, operand));
    if (simplified.is_error()) return std::nullopt;
    return simplified.value();
}

// ---------------------------------------------------------------------------
// safe_mul_i64 — overflow-safe signed 64-bit multiply used in exponential
// term scaling.  Returns false (and leaves out unchanged) on overflow.
// ---------------------------------------------------------------------------
[[nodiscard]] inline bool safe_mul_i64(long long lhs, long long rhs, long long& out) {
    if (lhs == 0 || rhs == 0) {
        out = 0;
        return true;
    }
    if (lhs == std::numeric_limits<long long>::min() || rhs == std::numeric_limits<long long>::min()) {
        return false;
    }
    const auto abs_lhs = lhs < 0 ? -lhs : lhs;
    const auto abs_rhs = rhs < 0 ? -rhs : rhs;
    if (abs_lhs > std::numeric_limits<long long>::max() / abs_rhs) return false;
    out = lhs * rhs;
    return true;
}

// ---------------------------------------------------------------------------
// exact_sign — returns -1, 0, or +1 for an exact rational expression; returns
// nullopt when the sign cannot be determined structurally.
// ---------------------------------------------------------------------------
[[nodiscard]] inline std::optional<int> exact_sign(ExprPtr expr) {
    symbolic::detail::LiteralRational lr;
    auto exact = symbolic::detail::try_get_exact_rational(expr, lr);
    if (exact.is_error() || !exact.value()) return std::nullopt;
    if (lr.value.numerator().is_zero()) return 0;
    return lr.value.numerator().is_negative() ? -1 : 1;
}

// ---------------------------------------------------------------------------
// is_negative_infinity — structural test: matches -∞ = Unary(Neg, Infinity).
// ---------------------------------------------------------------------------
[[nodiscard]] inline bool is_negative_infinity(ExprPtr point) {
    const auto* unary = expr_cast<Unary>(point);
    return unary != nullptr && unary->op == UnaryOp::Neg && limit_is_infinity(unary->operand);
}

// ---------------------------------------------------------------------------
// Public wrapper declarations — implemented in limit_mrv_leading.cpp and
// limit_mrv_exp.cpp; called from limit_mrv.cpp.
// ---------------------------------------------------------------------------

// Leading-term wrappers (limit_mrv_leading.cpp)
[[nodiscard]] std::optional<Result<ExprPtr>> mrv_try_quotient_valuation_limit(
    ExprPtr expr, const Symbol& w_var, AstArena& arena, symbolic::CASContext& ctx);

[[nodiscard]] std::optional<Result<ExprPtr>> mrv_try_leading_power_limit(
    ExprPtr expr, const Symbol& w_var, AstArena& arena, symbolic::CASContext& ctx);

[[nodiscard]] std::optional<LeadingPower> mrv_leading_power_w(
    ExprPtr expr, const Symbol& w_var, symbolic::CASContext& ctx);

// Exponential-product wrapper (limit_mrv_exp.cpp)
[[nodiscard]] std::optional<Result<ExprPtr>> mrv_try_exponential_product_limit(
    ExprPtr expr, const Symbol& var, AstArena& arena, symbolic::CASContext& ctx);

} // namespace cas::calculus
