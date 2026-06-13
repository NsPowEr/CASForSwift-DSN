#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

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
        // HC-F70-A33: poll for external interrupt at outer pivot loop entry.
        if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<MatrixExpr>(chk.error());
        std::size_t pivot_row = result.rows();
        std::optional<PivotScore> best_score;

        for (std::size_t i = r; i < result.rows(); ++i) {
            ExprPtr val = result(i, c);
            if (is_zero_expr(val)) continue;
            if (!is_structurally_nonzero(val)) continue;

            PivotScore score = make_pivot_score(val, ctx);
            if (!best_score.has_value() || *best_score < score) {
                best_score = score;
                pivot_row = i;
            }
            // Early exit: a literal-nonzero, degree-0 pivot dominates any
            // other candidate this column can offer. (Higher certainty
            // would be impossible; lower degree would be impossible too.)
            if (score.certainty == 3 && score.neg_total_degree == 0) break;
        }

        if (pivot_row == result.rows()) {
            continue;
        }

        if (pivot_row != r) {
            swap_rows(result, r, pivot_row);
        }

        ExprPtr current_pivot = result(r, c);

        for (std::size_t i = r + 1; i < result.rows(); ++i) {
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
