#include "cas/linalg/Matrix.hpp"

#include <optional>
#include <string>
#include <utility>

namespace cas::linalg {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] ExprPtr neg(symbolic::CASContext& ctx, ExprPtr operand) {
    return ctx.arena().make<Unary>(UnaryOp::Neg, operand);
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

[[nodiscard]] bool is_structurally_nonzero(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return !i->value.is_zero();
    if (const auto* r = expr_cast<RationalLit>(expr)) return !r->numerator.is_zero();
    if (expr_is<Constant>(expr)) {
        auto c = expr_ref<Constant>(expr).value;
        return c != MathConstant::NaN;
    }
    if (expr_is<Symbol>(expr)) return true;
    return !is_zero_expr(expr);
}

[[nodiscard]] bool is_known_nonzero_expr(ExprPtr expr) {
    return is_structurally_nonzero(expr);
}

[[nodiscard]] Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs)) return ok(integer(ctx, 0));
    if (is_one_expr(rhs)) return ok(lhs);
    if (lhs == rhs) return ok(integer(ctx, 1));
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs));
}

[[nodiscard]] Result<MatrixExpr> minor_matrix(const MatrixExpr& matrix, std::size_t skip_row, std::size_t skip_col) {
    const std::size_t n = matrix.rows();
    MatrixExpr result(n - 1U, n - 1U);
    std::size_t target_row = 0U;
    for (std::size_t row = 0; row < n; ++row) {
        if (row == skip_row) continue;
        std::size_t target_col = 0U;
        for (std::size_t col = 0; col < n; ++col) {
            if (col == skip_col) continue;
            result(target_row, target_col) = matrix(row, col);
            ++target_col;
        }
        ++target_row;
    }
    return ok(std::move(result));
}

}  // namespace

Result<MatrixExpr> inverse(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<MatrixExpr>(make_error(CASErrorKind::InvalidArgument, "Inverse requires a square matrix"));
    }

    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(MatrixExpr(0U, 0U));

    // Small symbolic inverses stay in adjugate form; larger matrices use RREF
    // with delayed normalization to control expression swell.
    if (n <= 3U) {
        auto det = determinant(matrix, ctx);
        if (det.is_error()) return fail<MatrixExpr>(det.error());
        if (is_zero_expr(det.value())) {
            return fail<MatrixExpr>(make_error(CASErrorKind::Undefined, "Singular matrix has no inverse"));
        }

        MatrixExpr result(n, n);
        if (n == 1U) {
            auto val = div_expr(ctx, integer(ctx, 1), matrix(0U, 0U));
            if (val.is_error()) return fail<MatrixExpr>(val.error());
            result(0U, 0U) = val.value();
            return ok(std::move(result));
        }

        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t col = 0; col < n; ++col) {
                auto minor = minor_matrix(matrix, row, col);
                auto cofactor = determinant(minor.value(), ctx);
                if (cofactor.is_error()) return fail<MatrixExpr>(cofactor.error());
                ExprPtr val = cofactor.value();
                if (((row + col) % 2U) != 0U) {
                    auto negated = simplify(ctx, neg(ctx, val));
                    if (negated.is_error()) return fail<MatrixExpr>(negated.error());
                    val = negated.value();
                }
                auto entry = div_expr(ctx, val, det.value());
                if (entry.is_error()) return fail<MatrixExpr>(entry.error());
                result(col, row) = entry.value();
            }
        }
        return ok(std::move(result));
    }

    MatrixExpr augmented(n, 2U * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            augmented(row, col) = matrix(row, col);
        }
        for (std::size_t col = 0; col < n; ++col) {
            augmented(row, n + col) = integer(ctx, (row == col) ? 1 : 0);
        }
    }

    auto reduced_res = rref(augmented, ctx);
    if (reduced_res.is_error()) return fail<MatrixExpr>(reduced_res.error());
    const auto& reduced = reduced_res.value();

    for (std::size_t i = 0; i < n; ++i) {
        if (!is_known_nonzero_expr(reduced(i, i))) {
            return fail<MatrixExpr>(make_error(CASErrorKind::Undefined, "Singular matrix has no inverse"));
        }
    }

    MatrixExpr result(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            result(row, col) = reduced(row, n + col);
        }
    }
    return ok(std::move(result));
}

}  // namespace cas::linalg
