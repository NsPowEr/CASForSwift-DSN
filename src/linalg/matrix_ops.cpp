#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "matrix_determinant.hpp"
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
            "Inner matrix dimensions must match for multiplication"));
    }
    MatrixExpr result(a.rows(), b.cols());
    for (std::size_t row = 0; row < a.rows(); ++row) {
        for (std::size_t col = 0; col < b.cols(); ++col) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t k = 0; k < a.cols(); ++k) {
                auto prod = mul_expr(ctx, a(row, k), b(k, col));
                if (prod.is_error()) return fail<MatrixExpr>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<MatrixExpr>(next_sum.error());
                sum = next_sum.value();
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
            result(i, j) = (i == j) ? integer(ctx, 1) : integer(ctx, 0);
        }
    }
    return ok(std::move(result));
}

Result<std::size_t> rank(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    auto b_res = bareiss(matrix, ctx);
    if (b_res.is_error()) return fail<std::size_t>(b_res.error());
    const auto& b = b_res.value();
    
    std::size_t r = 0;
    for (std::size_t i = 0; i < b.rows(); ++i) {
        bool row_nonzero = false;
        for (std::size_t j = 0; j < b.cols(); ++j) {
            if (!is_zero_expr(b(i, j))) {
                row_nonzero = true;
                break;
            }
        }
        if (row_nonzero) ++r;
    }
    return ok(r);
}

Result<ExprPtr> trace(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Trace requires a square matrix"));
    }
    ExprPtr result = integer(ctx, 0);
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        auto next = add_expr(ctx, result, matrix(i, i));
        if (next.is_error()) return next;
        result = next.value();
    }
    return ok(result);
}

Result<ExprPtr> determinant(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Determinant requires a square matrix"));
    }
    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(integer(ctx, 1));
    if (n == 1U) return ok(matrix(0U, 0U));

    auto tri = detail::determinant_tridiagonal_if_applicable(matrix, ctx);
    if (tri.is_error()) return fail<ExprPtr>(tri.error());
    if (tri.value().has_value()) return ok(*tri.value());

    auto van = detail::determinant_vandermonde_if_applicable(matrix, ctx);
    if (van.is_error()) return fail<ExprPtr>(van.error());
    if (van.value().has_value()) return ok(*van.value());

    auto circ = detail::determinant_circulant_if_applicable(matrix, ctx);
    if (circ.is_error()) return fail<ExprPtr>(circ.error());
    if (circ.value().has_value()) return ok(*circ.value());

    auto toep = detail::determinant_toeplitz_if_applicable(matrix, ctx);
    if (toep.is_error()) return fail<ExprPtr>(toep.error());
    if (toep.value().has_value()) return ok(*toep.value());

    auto band = detail::determinant_banded_if_applicable(matrix, ctx);
    if (band.is_error()) return fail<ExprPtr>(band.error());
    if (band.value().has_value()) return ok(*band.value());

    return detail::bareiss_determinant(matrix, ctx);
}

}  // namespace cas::linalg
