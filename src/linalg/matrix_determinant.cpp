// CAS-F4.1d — Bareiss fraction-free determinant (consolidato).
//
// Precedentemente: loop Bareiss replicato inline in matrix_ops::determinant.
// Ora: implementazione unica qui, riusata da matrix_ops.cpp via
// detail::bareiss_determinant. matrix_bareiss.cpp resta separato: calcola la
// forma RREF intera (matrice), non lo scalare determinante.
//
// Algoritmo (Bareiss 1968):
//   for k = 0..n-2:
//     pivot row = argmax_{i>=k} pivot_score(A[i][k])
//     swap rows; track sign flip
//     for i = k+1..n-1, j = k+1..n-1:
//       A[i][j] = (A[k][k]·A[i][j] − A[i][k]·A[k][j]) / prev_pivot
//     prev_pivot = A[k][k]
//   det = (±) A[n-1][n-1]
//
// La divisione per prev_pivot è esatta in Z (e nei domini integrali) →
// nessuna esplosione di denominatori, intermediate-expression-swell limitata.

#include "matrix_determinant.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cas::linalg::detail {

namespace {

[[nodiscard]] ExprPtr neg(symbolic::CASContext& ctx, ExprPtr operand) {
    return ctx.arena().make<Unary>(UnaryOp::Neg, operand);
}

}  // namespace

Result<ExprPtr> bareiss_determinant(const MatrixExpr& matrix,
                                      symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::InvalidArgument,
            "bareiss_determinant: matrix must be square", std::nullopt});
    }
    if (n == 0U) return ok(integer(ctx, 1));
    if (n == 1U) return ok(matrix(0U, 0U));

    std::vector<std::vector<ExprPtr>> work(n, std::vector<ExprPtr>(n));
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t col = 0; col < n; ++col)
            work[row][col] = matrix(row, col);

    ExprPtr previous_pivot = integer(ctx, 1);
    bool sign_flip = false;

    for (std::size_t k = 0; k + 1U < n; ++k) {
        std::size_t pivot_row = n;
        PivotScore best_score{-1, 0, 0};
        for (std::size_t i = k; i < n; ++i) {
            ExprPtr val = work[i][k];
            if (is_zero_expr(val)) continue;
            
            PivotScore score = make_pivot_score(val, ctx);
            if (pivot_row == n || score > best_score) {
                best_score = score;
                pivot_row = i;
                if (score.certainty == 3) break;
            }
        }
        if (pivot_row == n) return ok(integer(ctx, 0));
        if (pivot_row != k) {
            std::swap(work[pivot_row], work[k]);
            sign_flip = !sign_flip;
        }

        // Bareiss step: divisione esatta per previous_pivot.
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

}  // namespace cas::linalg::detail
