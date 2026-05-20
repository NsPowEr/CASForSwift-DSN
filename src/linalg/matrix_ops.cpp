#include "cas/linalg/Matrix.hpp"
#include "matrix_structured_determinant.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

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

[[nodiscard]] Result<ExprPtr> add_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs)) return ok(rhs);
    if (is_zero_expr(rhs)) return ok(lhs);
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Add, lhs, rhs));
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

[[nodiscard]] Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    if (is_zero_expr(lhs)) return ok(integer(ctx, 0));
    if (is_one_expr(rhs)) return ok(lhs);
    if (lhs == rhs) return ok(integer(ctx, 1));
    return simplify(ctx, ctx.arena().make<Binary>(BinaryOp::Div, lhs, rhs));
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
    if (expr_is<Constant>(expr)) {
        auto c = expr_ref<Constant>(expr).value;
        return c != MathConstant::NaN;
    }
    if (expr_is<Symbol>(expr)) return true;
    return !is_zero_expr(expr);
}

[[nodiscard]] Result<MatrixExpr> require_same_shape(const MatrixExpr& a, const MatrixExpr& b, std::string operation) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        return fail<MatrixExpr>(make_error(
            CASErrorKind::InvalidArgument,
            "Matrix dimensions must match for " + operation));
    }
    return ok(MatrixExpr(0U, 0U));
}

}  // namespace

Result<MatrixExpr> add(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx) {
    if (auto shape = require_same_shape(a, b, "addition"); shape.is_error()) {
        return shape;
    }
    MatrixExpr result(a.rows(), a.cols());
    for (std::size_t row = 0; row < a.rows(); ++row) {
        for (std::size_t col = 0; col < a.cols(); ++col) {
            auto value = add_expr(ctx, a(row, col), b(row, col));
            if (value.is_error()) return fail<MatrixExpr>(value.error());
            result(row, col) = value.value();
        }
    }
    return ok(std::move(result));
}

Result<MatrixExpr> subtract(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx) {
    if (auto shape = require_same_shape(a, b, "subtraction"); shape.is_error()) {
        return shape;
    }
    MatrixExpr result(a.rows(), a.cols());
    for (std::size_t row = 0; row < a.rows(); ++row) {
        for (std::size_t col = 0; col < a.cols(); ++col) {
            auto value = sub_expr(ctx, a(row, col), b(row, col));
            if (value.is_error()) return fail<MatrixExpr>(value.error());
            result(row, col) = value.value();
        }
    }
    return ok(std::move(result));
}

Result<MatrixExpr> multiply(const MatrixExpr& a, const MatrixExpr& b, symbolic::CASContext& ctx) {
    if (a.cols() != b.rows()) {
        return fail<MatrixExpr>(make_error(
            CASErrorKind::InvalidArgument,
            "Matrix dimensions must match for multiplication"));
    }

    MatrixExpr result(a.rows(), b.cols());
    for (std::size_t row = 0; row < a.rows(); ++row) {
        for (std::size_t col = 0; col < b.cols(); ++col) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t k = 0; k < a.cols(); ++k) {
                auto product = mul_expr(ctx, a(row, k), b(k, col));
                if (product.is_error()) return fail<MatrixExpr>(product.error());
                auto next = add_expr(ctx, sum, product.value());
                if (next.is_error()) return fail<MatrixExpr>(next.error());
                sum = next.value();
            }
            result(row, col) = sum;
        }
    }
    return ok(std::move(result));
}

Result<MatrixExpr> transpose(const MatrixExpr& matrix) {
    MatrixExpr result(matrix.cols(), matrix.rows());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            result(col, row) = matrix(row, col);
        }
    }
    return ok(std::move(result));
}

Result<MatrixExpr> identity(std::size_t n, symbolic::CASContext& ctx) {
    MatrixExpr result(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            result(i, j) = integer(ctx, (i == j) ? 1 : 0);
        }
    }
    return ok(std::move(result));
}

Result<ExprPtr> determinant(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Determinant requires a square matrix"));
    }

    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(integer(ctx, 1));
    if (n == 1U) return ok(matrix(0U, 0U));
    auto structured = detail::determinant_tridiagonal_if_applicable(matrix, ctx);
    if (structured.is_error()) return fail<ExprPtr>(structured.error());
    if (structured.value().has_value()) {
        return ok(*structured.value());
    }

    std::vector<std::vector<ExprPtr>> work(n, std::vector<ExprPtr>(n));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            work[row][col] = matrix(row, col);
        }
    }

    ExprPtr previous_pivot = integer(ctx, 1);
    bool sign_flip = false;

    for (std::size_t k = 0; k + 1U < n; ++k) {
        // Dynamic pivot selection (find simplest non-zero expression)
        std::size_t pivot_row = n;
        int best_score = -1;

        for (std::size_t i = k; i < n; ++i) {
            ExprPtr val = work[i][k];
            if (is_zero_expr(val)) continue;

            int score = 0;
            if (expr_is<IntegerLit>(val) || expr_is<RationalLit>(val)) {
                score = 1000; // Prefer numeric constants
            } else if (is_structurally_nonzero(val)) {
                score = 500; // Symbolic but likely non-zero
                // Penalize complexity to keep formulas smaller
                score -= static_cast<int>(std::min<std::size_t>(400, estimate_complexity(val)));
            } else {
                score = 1; // Last resort
            }

            if (score > best_score) {
                best_score = score;
                pivot_row = i;
            }
            if (score == 1000) break; // Found constant, good enough
        }

        if (pivot_row == n) {
            return ok(integer(ctx, 0));
        }

        if (pivot_row != k) {
            std::swap(work[pivot_row], work[k]);
            sign_flip = !sign_flip;
        }

        // Bareiss step: division by previous_pivot is guaranteed exact
        for (std::size_t row = k + 1U; row < n; ++row) {
            for (std::size_t col = k + 1U; col < n; ++col) {
                auto left = mul_expr(ctx, work[k][k], work[row][col]);
                if (left.is_error()) return left;
                auto right = mul_expr(ctx, work[row][k], work[k][col]);
                if (right.is_error()) return right;
                auto numerator = sub_expr(ctx, left.value(), right.value());
                if (numerator.is_error()) return numerator;
                auto value = div_expr(ctx, numerator.value(), previous_pivot);
                if (value.is_error()) return value;
                work[row][col] = value.value();
            }
        }
        previous_pivot = work[k][k];
    }

    ExprPtr det = work[n - 1U][n - 1U];
    if (sign_flip) {
        auto simplified = simplify(ctx, neg(ctx, det));
        if (simplified.is_error()) return simplified;
        det = simplified.value();
    }
    return ok(det);
}

Result<std::size_t> rank(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    auto reduced = rref(matrix, ctx);
    if (reduced.is_error()) return fail<std::size_t>(reduced.error());

    std::size_t r = 0;
    for (std::size_t row = 0; row < reduced.value().rows(); ++row) {
        bool row_has_nonzero = false;
        for (std::size_t col = 0; col < reduced.value().cols(); ++col) {
            if (!is_zero_expr(reduced.value()(row, col))) {
                row_has_nonzero = true;
                break;
            }
        }
        if (row_has_nonzero) {
            ++r;
        }
    }
    return ok(r);
}

Result<ExprPtr> trace(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    ExprPtr tr = integer(ctx, 0);
    const std::size_t n = std::min(matrix.rows(), matrix.cols());
    for (std::size_t i = 0; i < n; ++i) {
        auto next = add_expr(ctx, tr, matrix(i, i));
        if (next.is_error()) return next;
        tr = next.value();
    }
    return ok(tr);
}

}  // namespace cas::linalg
