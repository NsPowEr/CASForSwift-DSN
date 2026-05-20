#include "matrix_structured_determinant.hpp"

#include <utility>

namespace cas::linalg::detail {
namespace {

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] Result<ExprPtr> simplify(symbolic::CASContext& ctx, ExprPtr expr) {
    return ctx.simplify(expr);
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* integer_lit = expr_cast<IntegerLit>(expr)) {
        return integer_lit->value.is_zero();
    }
    if (const auto* rational_lit = expr_cast<RationalLit>(expr)) {
        return rational_lit->numerator.is_zero();
    }
    return false;
}

[[nodiscard]] bool is_one_expr(ExprPtr expr) {
    if (const auto* integer_lit = expr_cast<IntegerLit>(expr)) {
        return integer_lit->value == BigInt(1);
    }
    if (const auto* rational_lit = expr_cast<RationalLit>(expr)) {
        return rational_lit->numerator == BigInt(1) && rational_lit->denominator == BigInt(1);
    }
    return false;
}

[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(rhs)) return ok(lhs);
    if (lhs == rhs) return ok(integer(ctx, 0));
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs) || is_zero_expr(rhs)) return ok(integer(ctx, 0));
    if (is_one_expr(lhs)) return ok(rhs);
    if (is_one_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

[[nodiscard]] bool is_tridiagonal(const MatrixExpr& matrix) {
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            const std::size_t distance = (row > col) ? row - col : col - row;
            if (distance > 1U && !is_zero_expr(matrix(row, col))) {
                return false;
            }
        }
    }
    return true;
}

}  // namespace

Result<std::optional<ExprPtr>> determinant_tridiagonal_if_applicable(
    const MatrixExpr& matrix,
    symbolic::CASContext& ctx) {
    if (!is_tridiagonal(matrix)) {
        return ok(std::optional<ExprPtr>{});
    }

    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(std::optional<ExprPtr>{integer(ctx, 1)});
    if (n == 1U) return ok(std::optional<ExprPtr>{matrix(0U, 0U)});

    ExprPtr previous_previous = integer(ctx, 1);
    ExprPtr previous = matrix(0U, 0U);

    for (std::size_t index = 1U; index < n; ++index) {
        auto diagonal_term = mul_expr(ctx, matrix(index, index), previous);
        if (diagonal_term.is_error()) return fail<std::optional<ExprPtr>>(diagonal_term.error());

        auto off_diagonal_product = mul_expr(ctx, matrix(index, index - 1U), matrix(index - 1U, index));
        if (off_diagonal_product.is_error()) {
            return fail<std::optional<ExprPtr>>(off_diagonal_product.error());
        }

        auto correction = mul_expr(ctx, off_diagonal_product.value(), previous_previous);
        if (correction.is_error()) return fail<std::optional<ExprPtr>>(correction.error());

        auto next = sub_expr(ctx, diagonal_term.value(), correction.value());
        if (next.is_error()) return fail<std::optional<ExprPtr>>(next.error());

        previous_previous = previous;
        previous = next.value();
    }

    return ok(std::optional<ExprPtr>{previous});
}

}  // namespace cas::linalg::detail
