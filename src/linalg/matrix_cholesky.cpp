// CAS-F4.1c — Decomposizione Cholesky LDL^T senza radici quadrate.
//
// Per A simmetrica n×n calcola L unit-lower triangular e D diagonale t.c.
// A = L · D · L^T.  Più stabile di Cholesky classico L·L^T per matrici
// indefinite e produce factorizzazioni razionali se A è razionale (no sqrt).
//
// Formule (Golub-Van Loan §4.1.2):
//   D[j] = A[j][j] − Σ_{k<j} L[j][k]² · D[k]
//   L[i][j] = (A[i][j] − Σ_{k<j} L[i][k]·D[k]·L[j][k]) / D[j]    (per i > j)
//
// Pre-condizione: A simmetrica strutturalmente (verificata; le entry sono
// confrontate via simplify(A[i][j] − A[j][i])).  Fallisce con
// LINALG_ZERO_PIVOT se qualche D[j] = 0 (manca pivoting Bunch-Kaufman,
// debito esplicito documentato; per matrici positive definite questo caso
// non si verifica).

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <vector>

namespace cas::linalg {

namespace {

[[nodiscard]] bool symmetric(const MatrixExpr& A, symbolic::CASContext& ctx) {
    const std::size_t n = A.rows();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            auto diff = sub_expr(ctx, A(i, j), A(j, i));
            if (diff.is_error() || !is_zero_expr(diff.value())) return false;
        }
    }
    return true;
}

}  // namespace

Result<LDLTDecomposition> cholesky_ldlt(const MatrixExpr& matrix,
                                          symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<LDLTDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "cholesky_ldlt: matrix must be square", std::nullopt});
    }
    if (n == 0U) {
        return fail<LDLTDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "cholesky_ldlt: empty matrix", std::nullopt});
    }
    if (!symmetric(matrix, ctx)) {
        return fail<LDLTDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "cholesky_ldlt: matrix must be symmetric", std::nullopt});
    }

    MatrixExpr L(n, n);
    std::vector<ExprPtr> D(n);

    for (std::size_t j = 0; j < n; ++j) {
        // Compute D[j] = A[j][j] - Σ_{k<j} L[j][k]² · D[k]
        ExprPtr sum = integer(ctx, 0);
        for (std::size_t k = 0; k < j; ++k) {
            auto l2 = mul_expr(ctx, L(j, k), L(j, k));
            if (l2.is_error()) return fail<LDLTDecomposition>(l2.error());
            auto term = mul_expr(ctx, l2.value(), D[k]);
            if (term.is_error()) return fail<LDLTDecomposition>(term.error());
            auto next_sum = add_expr(ctx, sum, term.value());
            if (next_sum.is_error()) return fail<LDLTDecomposition>(next_sum.error());
            sum = next_sum.value();
        }
        auto dj_res = sub_expr(ctx, matrix(j, j), sum);
        if (dj_res.is_error()) return fail<LDLTDecomposition>(dj_res.error());
        D[j] = dj_res.value();

        // M-FIX: check nonzero D[j] safely (assumptions)
        if (!is_known_nonzero(D[j], ctx)) {
            return make_unimplemented<LDLTDecomposition>(
                "linalg", "cholesky_ldlt",
                "zero pivot D[" + std::to_string(j) + "] encountered",
                error::reason_codes::LINALG_ZERO_PIVOT,
                "Cholesky LDLT requires non-zero diagonal entries; check if matrix is positive definite",
                "L3-16");
        }

        L(j, j) = integer(ctx, 1);
        for (std::size_t k = j + 1; k < n; ++k) L(j, k) = integer(ctx, 0);

        // Compute L[i][j] = (A[i][j] - Σ_{k<j} L[i][k]·D[k]·L[j][k]) / D[j]
        for (std::size_t i = j + 1; i < n; ++i) {
            ExprPtr sum_l = integer(ctx, 0);
            for (std::size_t k = 0; k < j; ++k) {
                auto l_ik = L(i, k);
                auto d_k = D[k];
                auto l_jk = L(j, k);
                auto t1 = mul_expr(ctx, l_ik, d_k);
                if (t1.is_error()) return fail<LDLTDecomposition>(t1.error());
                auto t2 = mul_expr(ctx, t1.value(), l_jk);
                if (t2.is_error()) return fail<LDLTDecomposition>(t2.error());
                auto next_sum = add_expr(ctx, sum_l, t2.value());
                if (next_sum.is_error()) return fail<LDLTDecomposition>(next_sum.error());
                sum_l = next_sum.value();
            }
            auto num = sub_expr(ctx, matrix(i, j), sum_l);
            if (num.is_error()) return fail<LDLTDecomposition>(num.error());
            auto lij = div_expr(ctx, num.value(), D[j]);
            if (lij.is_error()) return fail<LDLTDecomposition>(lij.error());
            
            // Normalize with together
            auto tg = algebra::together(lij.value(), ctx);
            L(i, j) = tg.is_ok() ? tg.value() : lij.value();
        }
    }

    return ok(LDLTDecomposition{std::move(L), std::move(D)});
}

}  // namespace cas::linalg
