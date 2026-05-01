#include "cas/linalg/Matrix.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] Result<ExprPtr> simplify(symbolic::CASContext& ctx, ExprPtr expr) {
    return ctx.simplify(expr);
}

[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Mul, lhs, rhs));
}

[[nodiscard]] Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs));
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

[[nodiscard]] std::size_t estimate_complexity(ExprPtr expr) {
    if (!expr) return 0;
    switch (expr->kind) {
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Constant:
        return 1;
    case ExprKind::Symbol:
        return 2;
    case ExprKind::Unary:
        return 1 + estimate_complexity(expr_cast<Unary>(expr)->operand);
    case ExprKind::Binary: {
        const auto* b = expr_cast<Binary>(expr);
        return 1 + estimate_complexity(b->left) + estimate_complexity(b->right);
    }
    case ExprKind::Sum: {
        std::size_t c = 1;
        for (auto t : expr_cast<Sum>(expr)->terms) c += estimate_complexity(t);
        return c;
    }
    case ExprKind::Product: {
        std::size_t c = 1;
        for (auto f : expr_cast<Product>(expr)->factors) c += estimate_complexity(f);
        return c;
    }
    case ExprKind::FuncCall: {
        std::size_t c = 1;
        for (auto a : expr_cast<FuncCall>(expr)->args) c += estimate_complexity(a);
        return c;
    }
    default:
        return 5;
    }
}

[[nodiscard]] bool is_structurally_nonzero(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* i = expr_cast<IntegerLit>(expr)) return !i->value.is_zero();
    if (const auto* r = expr_cast<RationalLit>(expr)) return !r->numerator.is_zero();
    if (expr_is<Constant>(expr)) return true;
    if (expr_is<Symbol>(expr)) return true;
    return !is_zero_expr(expr);
}

void swap_rows(MatrixExpr& matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) return;
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        std::swap(matrix(lhs, col), matrix(rhs, col));
    }
}

}  // namespace

Result<MatrixExpr> bareiss(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    MatrixExpr result(matrix.rows(), matrix.cols(), matrix.elements());
    std::size_t r = 0;
    ExprPtr prev_pivot = integer(ctx, 1);

    for (std::size_t c = 0; c < result.cols() && r < result.rows(); ++c) {
        std::size_t pivot_row = result.rows();
        int best_score = -1;

        for (std::size_t i = r; i < result.rows(); ++i) {
            ExprPtr val = result(i, c);
            if (is_zero_expr(val)) continue;

            int score = 0;
            if (expr_is<IntegerLit>(val) || expr_is<RationalLit>(val)) {
                score = 1000;
            } else if (is_structurally_nonzero(val)) {
                score = 500;
                score -= static_cast<int>(std::min<std::size_t>(400, estimate_complexity(val)));
            } else {
                score = 1;
            }

            if (score > best_score) {
                best_score = score;
                pivot_row = i;
            }
            if (score == 1000) break;
        }

        if (pivot_row == result.rows()) {
            continue;
        }

        if (pivot_row != r) {
            swap_rows(result, r, pivot_row);
        }

        ExprPtr current_pivot = result(r, c);

        for (std::size_t i = 0; i < result.rows(); ++i) {
            if (i == r) continue;
            for (std::size_t j = c + 1; j < result.cols(); ++j) {
                auto left_res = mul_expr(ctx, current_pivot, result(i, j));
                if (left_res.is_error()) return fail<MatrixExpr>(left_res.error());
                
                auto right_res = mul_expr(ctx, result(i, c), result(r, j));
                if (right_res.is_error()) return fail<MatrixExpr>(right_res.error());
                
                auto num_res = sub_expr(ctx, left_res.value(), right_res.value());
                if (num_res.is_error()) return fail<MatrixExpr>(num_res.error());
                
                auto val_res = div_expr(ctx, num_res.value(), prev_pivot);
                if (val_res.is_error()) return fail<MatrixExpr>(val_res.error());
                
                result(i, j) = val_res.value();
            }
            result(i, c) = integer(ctx, 0);
        }

        prev_pivot = current_pivot;
        ++r;
    }

    return ok(std::move(result));
}

}  // namespace cas::linalg
