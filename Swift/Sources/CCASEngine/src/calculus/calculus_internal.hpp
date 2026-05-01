#pragma once

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <vector>

namespace cas::calculus {

struct QuotientView {
    ExprPtr numerator;
    ExprPtr denominator;
};

/// @brief Differentiates transcendental functions.
/// Handles asin, acos, atan, sinh, cosh, tanh.
[[nodiscard]] Result<ExprPtr> differentiate_transcendental(
    const std::string& name,
    ExprPtr argument,
    const Symbol& var,
    symbolic::CASContext& context);

/// @brief Integration by parts: ∫ u dv = uv - ∫ v du
[[nodiscard]] Result<ExprPtr> integrate_by_parts(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context);

[[nodiscard]] ExprPtr limit_make_integer(AstArena& arena, long long value);
[[nodiscard]] ExprPtr limit_make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool limit_is_zero(ExprPtr expr);
[[nodiscard]] bool limit_is_infinity(ExprPtr expr);
[[nodiscard]] std::optional<QuotientView> extract_quotient_view(ExprPtr expr, AstArena& arena);
[[nodiscard]] std::optional<ExprPtr> cancel_common_linear_factor(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_polynomial_pole_limit(
    ExprPtr numerator,
    ExprPtr denominator,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] std::optional<Result<ExprPtr>> try_logarithmic_root_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_infinite_limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    AstArena& arena);

}  // namespace cas::calculus
