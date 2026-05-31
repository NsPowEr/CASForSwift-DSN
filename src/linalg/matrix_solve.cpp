#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/error_helpers.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

void swap_rows(MatrixExpr& matrix, std::size_t lhs, std::size_t rhs) {
    if (lhs == rhs) return;
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        std::swap(matrix(lhs, col), matrix(rhs, col));
    }
}

[[nodiscard]] Result<void> scale_row(MatrixExpr& matrix, std::size_t row, ExprPtr divisor, symbolic::CASContext& ctx) {
    if (is_one_expr(divisor)) return ok();
    if (!is_known_nonzero(divisor, ctx)) {
        return make_unimplemented<void>(
            "linalg", "scale_row",
            "RREF pivot expression cannot be decided nonzero under current assumptions",
            error::reason_codes::LINALG_RREF_UNDECIDABLE_PIVOT,
            "Add assumptions about pivot variables or use Bareiss integer-preserving elimination",
            "L3-13");
    }
    for (std::size_t col = 0; col < matrix.cols(); ++col) {
        auto value = div_expr(ctx, matrix(row, col), divisor);
        if (value.is_error()) return fail<void>(value.error());
        matrix(row, col) = value.value();
    }
    return ok();
}

[[nodiscard]] Result<void> eliminate_row(MatrixExpr& matrix, std::size_t target_row, std::size_t pivot_row, std::size_t col, symbolic::CASContext& ctx) {
    ExprPtr factor = matrix(target_row, col);
    if (is_zero_expr(factor)) return ok();

    for (std::size_t j = col; j < matrix.cols(); ++j) {
        auto product = mul_expr(ctx, factor, matrix(pivot_row, j));
        if (product.is_error()) return fail<void>(product.error());
        auto diff = sub_expr(ctx, matrix(target_row, j), product.value());
        if (diff.is_error()) return fail<void>(diff.error());
        matrix(target_row, j) = diff.value();
    }
    return ok();
}

} // namespace

Result<MatrixExpr> rref(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t rows = matrix.rows();
    const std::size_t cols = matrix.cols();
    MatrixExpr result(rows, cols, matrix.elements());

    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
        // Find pivot
        std::size_t sel_row = rows;
        PivotScore best_score{-1, 0, 0};
        for (std::size_t i = pivot_row; i < rows; ++i) {
            if (is_zero_expr(result(i, col))) continue;
            PivotScore score = make_pivot_score(result(i, col), ctx);
            if (sel_row == rows || score > best_score) {
                best_score = score;
                sel_row = i;
                if (score.certainty == 3) break;
            }
        }

        if (sel_row == rows || best_score.certainty < 0) continue;

        swap_rows(result, pivot_row, sel_row);

        auto scale_res = scale_row(result, pivot_row, result(pivot_row, col), ctx);
        if (scale_res.is_error()) return fail<MatrixExpr>(scale_res.error());

        for (std::size_t i = 0; i < rows; ++i) {
            if (i != pivot_row) {
                auto elim_res = eliminate_row(result, i, pivot_row, col, ctx);
                if (elim_res.is_error()) return fail<MatrixExpr>(elim_res.error());
            }
        }
        ++pivot_row;
    }

    return ok(result);
}

Result<std::vector<ExprPtr>> linsolve(const MatrixExpr& a, const std::vector<ExprPtr>& b, symbolic::CASContext& ctx) {
    if (a.rows() != b.size()) {
        return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "solve_linear_system: matrix rows and vector size mismatch"));
    }

    const std::size_t n = a.rows();
    const std::size_t m = a.cols();

    // Augmented matrix [A|b]
    MatrixExpr augmented(n, m + 1);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            augmented(i, j) = a(i, j);
        }
        augmented(i, m) = b[i];
    }

    auto rref_res = rref(augmented, ctx);
    if (rref_res.is_error()) return fail<std::vector<ExprPtr>>(rref_res.error());
    const auto& rref = rref_res.value();

    std::vector<ExprPtr> solution(m);
    std::vector<bool> pivot_found(m, false);

    for (std::size_t i = 0; i < n; ++i) {
        std::size_t col = 0;
        while (col < m && is_zero_expr(rref(i, col))) col++;

        if (col == m) {
            // Check consistency
            if (!is_zero_expr(rref(i, m))) {
                return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::Undefined, "solve_linear_system: inconsistent system"));
            }
            continue;
        }
        
        pivot_found[col] = true;
    }

    // 2. Assign symbols to free variables
    for (std::size_t j = 0; j < m; ++j) {
        if (!pivot_found[j]) {
            std::string name = "c" + std::to_string(j + 1);
            solution[j] = ctx.arena().make<Symbol>(name);
        }
    }

    // 3. Back-substitution for pivots
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        std::size_t col = 0;
        while (col < m && is_zero_expr(rref(static_cast<std::size_t>(i), col))) col++;

        if (col < m) {
            ExprPtr val = rref(static_cast<std::size_t>(i), m);
            for (std::size_t j = col + 1; j < m; ++j) {
                if (!is_zero_expr(rref(static_cast<std::size_t>(i), j))) {
                    auto term = mul_expr(ctx, rref(static_cast<std::size_t>(i), j), solution[j]);
                    if (term.is_error()) return fail<std::vector<ExprPtr>>(term.error());
                    auto next_val = sub_expr(ctx, val, term.value());
                    if (next_val.is_error()) return fail<std::vector<ExprPtr>>(next_val.error());
                    val = next_val.value();
                }
            }
            solution[col] = val;
        }
    }

    for (ExprPtr& s : solution) {
        auto simp = simplify(ctx, s);
        if (simp.is_ok()) s = simp.value();
    }

    return ok(solution);
}

Result<std::vector<std::vector<ExprPtr>>> null_space(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    auto rref_res = rref(matrix, ctx);
    if (rref_res.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(rref_res.error());
    const auto& rref = rref_res.value();

    const std::size_t n = matrix.rows();
    const std::size_t m = matrix.cols();

    std::vector<int> pivot_col(n, -1);
    std::vector<bool> is_pivot_col(m, false);
    std::size_t r = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            if (!is_zero_expr(rref(i, j))) {
                pivot_col[i] = static_cast<int>(j);
                is_pivot_col[j] = true;
                ++r;
                break;
            }
        }
    }

    std::vector<std::vector<ExprPtr>> basis;
    for (std::size_t j = 0; j < m; ++j) {
        if (!is_pivot_col[j]) {
            std::vector<ExprPtr> v(m, integer(ctx, 0));
            v[j] = integer(ctx, 1);
            for (std::size_t i = 0; i < r; ++i) {
                if (pivot_col[i] != -1) {
                    auto neg_val = simplify(ctx, ctx.arena().make<Unary>(UnaryOp::Neg, rref(i, j)));
                    if (neg_val.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(neg_val.error());
                    v[static_cast<std::size_t>(pivot_col[i])] = neg_val.value();
                }
            }
            basis.push_back(std::move(v));
        }
    }
    return ok(std::move(basis));
}

}  // namespace cas::linalg
