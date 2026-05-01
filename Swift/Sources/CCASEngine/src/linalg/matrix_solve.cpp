#include "cas/linalg/Matrix.hpp"

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

[[nodiscard]] bool is_known_nonzero_expr(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* integer_lit = expr_cast<IntegerLit>(expr)) {
        return !integer_lit->value.is_zero();
    }
    if (const auto* rational_lit = expr_cast<RationalLit>(expr)) {
        return !rational_lit->numerator.is_zero();
    }
    return !is_zero_expr(expr);
}

void swap_rows(MatrixExpr& matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) return;
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        std::swap(matrix(lhs, col), matrix(rhs, col));
    }
}

[[nodiscard]] Result<void> scale_row(MatrixExpr& matrix, std::size_t row, ExprPtr divisor, symbolic::CASContext& ctx) {
    if (!is_known_nonzero_expr(divisor)) {
        return fail<void>(make_error(CASErrorKind::Unimplemented, "RREF pivot is not decidably nonzero"));
    }
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        auto value = div_expr(ctx, matrix(row, col), divisor);
        if (value.is_error()) return fail<void>(value.error());
        matrix(row, col) = value.value();
    }
    return ok();
}

[[nodiscard]] Result<void> eliminate_row(
    MatrixExpr& matrix,
    std::size_t target_row,
    std::size_t pivot_row,
    std::size_t start_col,
    symbolic::CASContext& ctx) {
    ExprPtr factor = matrix(target_row, start_col);
    if (is_zero_expr(factor)) return ok();

    for (std::size_t col = start_col; col < matrix.cols(); ++col) {
        auto product = mul_expr(ctx, factor, matrix(pivot_row, col));
        if (product.is_error()) return fail<void>(product.error());
        auto value = sub_expr(ctx, matrix(target_row, col), product.value());
        if (value.is_error()) return fail<void>(value.error());
        matrix(target_row, col) = value.value();
    }
    return ok();
}

[[nodiscard]] bool row_has_nonzero_coefficients(const MatrixExpr& matrix, std::size_t row, std::size_t coefficient_cols) {
    for (std::size_t col = 0; col < coefficient_cols; ++col) {
        if (!is_zero_expr(matrix(row, col))) return true;
    }
    return false;
}

[[nodiscard]] std::optional<std::size_t> first_nonzero_column(
    const MatrixExpr& matrix,
    std::size_t row,
    std::size_t coefficient_cols) {
    for (std::size_t col = 0; col < coefficient_cols; ++col) {
        if (!is_zero_expr(matrix(row, col))) return col;
    }
    return std::nullopt;
}

}  // namespace

Result<MatrixExpr> rref(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    MatrixExpr result(matrix.rows(), matrix.cols(), matrix.elements());
    std::size_t lead = 0U;

    for (std::size_t pivot_row = 0; pivot_row < result.rows() && lead < result.cols(); ++pivot_row) {
        std::size_t candidate_row = pivot_row;
        while (candidate_row < result.rows() && is_zero_expr(result(candidate_row, lead))) {
            ++candidate_row;
        }
        while (candidate_row == result.rows()) {
            ++lead;
            if (lead == result.cols()) {
                return ok(std::move(result));
            }
            candidate_row = pivot_row;
            while (candidate_row < result.rows() && is_zero_expr(result(candidate_row, lead))) {
                ++candidate_row;
            }
        }

        swap_rows(result, pivot_row, candidate_row);
        ExprPtr pivot = result(pivot_row, lead);
        auto scaled = scale_row(result, pivot_row, pivot, ctx);
        if (scaled.is_error()) return fail<MatrixExpr>(scaled.error());

        for (std::size_t row = 0; row < result.rows(); ++row) {
            if (row == pivot_row) continue;
            auto eliminated = eliminate_row(result, row, pivot_row, lead, ctx);
            if (eliminated.is_error()) return fail<MatrixExpr>(eliminated.error());
        }
        ++lead;
    }

    return ok(std::move(result));
}

Result<std::vector<ExprPtr>> linsolve(const MatrixExpr& a, const std::vector<ExprPtr>& b, symbolic::CASContext& ctx) {
    if (a.rows() != b.size()) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "Linear system dimensions must satisfy A.rows == b.size"));
    }

    MatrixExpr augmented(a.rows(), a.cols() + 1U);
    for (std::size_t row = 0; row < a.rows(); ++row) {
        for (std::size_t col = 0; col < a.cols(); ++col) {
            augmented(row, col) = a(row, col);
        }
        augmented(row, a.cols()) = b[row];
    }

    auto reduced = rref(augmented, ctx);
    if (reduced.is_error()) return fail<std::vector<ExprPtr>>(reduced.error());

    for (std::size_t row = 0; row < reduced.value().rows(); ++row) {
        if (!row_has_nonzero_coefficients(reduced.value(), row, a.cols()) &&
            !is_zero_expr(reduced.value()(row, a.cols()))) {
            return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::Undefined, "Linear system is inconsistent"));
        }
    }

    std::vector<std::optional<std::size_t>> pivot_row_for_col(a.cols());
    for (std::size_t row = 0; row < reduced.value().rows(); ++row) {
        auto pivot_col = first_nonzero_column(reduced.value(), row, a.cols());
        if (pivot_col.has_value()) {
            pivot_row_for_col[*pivot_col] = row;
        }
    }

    std::vector<ExprPtr> solution(a.cols(), integer(ctx, 0));
    std::vector<ExprPtr> free_parameters(a.cols());
    std::size_t parameter_index = 1U;
    for (std::size_t col = 0; col < a.cols(); ++col) {
        if (!pivot_row_for_col[col].has_value()) {
            free_parameters[col] = ctx.arena().make<Symbol>("t" + std::to_string(parameter_index++));
            solution[col] = free_parameters[col];
        }
    }

    for (std::size_t col = a.cols(); col > 0U; --col) {
        const std::size_t variable_col = col - 1U;
        if (!pivot_row_for_col[variable_col].has_value()) continue;

        const std::size_t row = *pivot_row_for_col[variable_col];
        ExprPtr value = reduced.value()(row, a.cols());
        for (std::size_t free_col = variable_col + 1U; free_col < a.cols(); ++free_col) {
            if (!free_parameters[free_col]) continue;
            ExprPtr coefficient = reduced.value()(row, free_col);
            if (is_zero_expr(coefficient)) continue;
            auto product = mul_expr(ctx, coefficient, solution[free_col]);
            if (product.is_error()) return fail<std::vector<ExprPtr>>(product.error());
            auto next = sub_expr(ctx, value, product.value());
            if (next.is_error()) return fail<std::vector<ExprPtr>>(next.error());
            value = next.value();
        }
        solution[variable_col] = value;
    }

    return ok(std::move(solution));
}

}  // namespace cas::linalg
