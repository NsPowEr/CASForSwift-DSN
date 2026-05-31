// CAS-L3-17 — Symbolic LU decomposition (Doolittle algorithm).
//
// Decomposes square matrix A as A = L · U where:
//   L: lower triangular, L[i][i] = 1 (unit diagonal)
//   U: upper triangular
//
// Doolittle classical formulas (Golub-Van Loan §3.2):
//   U[k][j] = A[k][j] - Σ_{s<k} L[k][s] · U[s][j]      for j ≥ k
//   L[i][k] = (A[i][k] - Σ_{s<k} L[i][s] · U[s][k]) / U[k][k]  for i > k
//
// No pivoting; zero pivot returns Unimplemented (caller can fall back
// to bareiss for det/rank, or permuted LU when implemented).

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cas::linalg {

Result<LUDecomposition> lu_decompose(const MatrixExpr& matrix,
                                      symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<LUDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_decompose: matrix must be square", std::nullopt});
    }
    if (n == 0U) {
        return fail<LUDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_decompose: empty matrix", std::nullopt});
    }

    MatrixExpr L(n, n);
    MatrixExpr U(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            L(i, j) = (i == j) ? integer(ctx, 1) : integer(ctx, 0);
            U(i, j) = integer(ctx, 0);
        }
    }

    for (std::size_t k = 0; k < n; ++k) {
        // Compute U[k][j] for j ≥ k.
        for (std::size_t j = k; j < n; ++j) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t s = 0; s < k; ++s) {
                auto prod = mul_expr(ctx, L(k, s), U(s, j));
                if (prod.is_error()) return fail<LUDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<LUDecomposition>(next_sum.error());
                sum = next_sum.value();
            }
            auto u_kj = sub_expr(ctx, matrix(k, j), sum);
            if (u_kj.is_error()) return fail<LUDecomposition>(u_kj.error());
            U(k, j) = u_kj.value();
        }

        // Check pivot U[k][k] ≠ 0.
        if (is_zero_expr(U(k, k))) {
            return make_unimplemented<LUDecomposition>(
                "linalg", "lu_decompose",
                "zero pivot at row " + std::to_string(k) + " in Doolittle LU",
                error::reason_codes::LINALG_ZERO_PIVOT,
                "Use lu_decompose_pivoted for partial pivoting or Bareiss for symbolic det/rank",
                "L3-17");
        }

        // Compute L[i][k] for i > k.
        for (std::size_t i = k + 1; i < n; ++i) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t s = 0; s < k; ++s) {
                auto prod = mul_expr(ctx, L(i, s), U(s, k));
                if (prod.is_error()) return fail<LUDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<LUDecomposition>(next_sum.error());
                sum = next_sum.value();
            }
            auto num = sub_expr(ctx, matrix(i, k), sum);
            if (num.is_error()) return fail<LUDecomposition>(num.error());
            auto quot = div_expr(ctx, num.value(), U(k, k));
            if (quot.is_error()) return fail<LUDecomposition>(quot.error());
            
            auto t = algebra::together(quot.value(), ctx);
            ExprPtr norm = t.is_ok() ? t.value() : quot.value();
            auto simp = ctx.simplify(norm);
            L(i, k) = simp.is_ok() ? simp.value() : norm;
        }
    }

    return ok(LUDecomposition{std::move(L), std::move(U)});
}

Result<PLUDecomposition> lu_decompose_pivoted(const MatrixExpr& matrix,
                                                symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<PLUDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_decompose_pivoted: matrix must be square", std::nullopt});
    }
    if (n == 0U) {
        return fail<PLUDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_decompose_pivoted: empty matrix", std::nullopt});
    }

    // Work on copy of A; permute rows as we go.
    MatrixExpr A_perm(n, n, matrix.elements());
    std::vector<std::size_t> P(n);
    for (std::size_t i = 0; i < n; ++i) P[i] = i;

    MatrixExpr L(n, n);
    MatrixExpr U(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            L(i, j) = (i == j) ? integer(ctx, 1) : integer(ctx, 0);
            U(i, j) = integer(ctx, 0);
        }
    }

    for (std::size_t k = 0; k < n; ++k) {
        // Find best pivot row among rows [k, n) for column k via PivotScore:
        // numerico esatto > simbolico nonzero (penalizzato per complessità) > zero.
        std::size_t pivot_row = n;
        PivotScore best_score{-1, 0, 0};
        for (std::size_t i = k; i < n; ++i) {
            // Compute candidate pivot value = A_perm[i][k] - Σ L[i][s]·U[s][k]
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t s = 0; s < k; ++s) {
                auto prod = mul_expr(ctx, L(i, s), U(s, k));
                if (prod.is_error()) return fail<PLUDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<PLUDecomposition>(next_sum.error());
                sum = next_sum.value();
            }
            auto cand_res = sub_expr(ctx, A_perm(i, k), sum);
            if (cand_res.is_error()) return fail<PLUDecomposition>(cand_res.error());
            ExprPtr v = cand_res.value();
            
            PivotScore score = is_zero_expr(v) ? PivotScore{-1, 0, 0} : make_pivot_score(v, ctx);
            if (pivot_row == n || score > best_score) {
                best_score = score;
                pivot_row = i;
                if (score.certainty == 3) break;  // numerico esatto, ottimo
            }
        }
        if (best_score.certainty < 0) pivot_row = n;
        if (pivot_row == n) {
            return make_unimplemented<PLUDecomposition>(
                "linalg", "lu_decompose_pivoted",
                "singular matrix: all-zero pivot column " + std::to_string(k),
                error::reason_codes::LINALG_SINGULAR_MATRIX,
                "Check matrix rank via Bareiss or use pseudo-inverse for rank-deficient systems",
                "L3-17");
        }
        // Swap rows k and pivot_row in A_perm and in L (already-built columns).
        if (pivot_row != k) {
            std::swap(P[k], P[pivot_row]);
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(A_perm(k, j), A_perm(pivot_row, j));
            }
            // Swap L entries for columns 0..k-1 (already computed).
            for (std::size_t j = 0; j < k; ++j) {
                std::swap(L(k, j), L(pivot_row, j));
            }
        }
        // Now compute U[k][j] for j ≥ k.
        for (std::size_t j = k; j < n; ++j) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t s = 0; s < k; ++s) {
                auto prod = mul_expr(ctx, L(k, s), U(s, j));
                if (prod.is_error()) return fail<PLUDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<PLUDecomposition>(next_sum.error());
                sum = next_sum.value();
            }
            auto u_kj = sub_expr(ctx, A_perm(k, j), sum);
            if (u_kj.is_error()) return fail<PLUDecomposition>(u_kj.error());
            U(k, j) = u_kj.value();
        }
        // L[i][k] per i > k.
        for (std::size_t i = k + 1; i < n; ++i) {
            ExprPtr sum = integer(ctx, 0);
            for (std::size_t s = 0; s < k; ++s) {
                auto prod = mul_expr(ctx, L(i, s), U(s, k));
                if (prod.is_error()) return fail<PLUDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, sum, prod.value());
                if (next_sum.is_error()) return fail<PLUDecomposition>(next_sum.error());
                sum = next_sum.value();
            }
            auto num = sub_expr(ctx, A_perm(i, k), sum);
            if (num.is_error()) return fail<PLUDecomposition>(num.error());
            auto quot = div_expr(ctx, num.value(), U(k, k));
            if (quot.is_error()) return fail<PLUDecomposition>(quot.error());
            
            auto t = algebra::together(quot.value(), ctx);
            ExprPtr norm = t.is_ok() ? t.value() : quot.value();
            auto simp = ctx.simplify(norm);
            L(i, k) = simp.is_ok() ? simp.value() : norm;
        }
    }
    return ok(PLUDecomposition{std::move(P), std::move(L), std::move(U)});
}

Result<std::vector<ExprPtr>> lu_solve(const LUDecomposition& lu,
                                       const std::vector<ExprPtr>& b,
                                       symbolic::CASContext& ctx) {
    const std::size_t n = lu.L.rows();
    if (lu.L.cols() != n || lu.U.rows() != n || lu.U.cols() != n) {
        return fail<std::vector<ExprPtr>>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_solve: non-square factors", std::nullopt});
    }
    if (b.size() != n) {
        return fail<std::vector<ExprPtr>>(CASError{
            CASErrorKind::InvalidArgument,
            "lu_solve: b size mismatch", std::nullopt});
    }

    // Forward substitution: L·y = b. L unit triangular → y[i] = b[i] - Σ L[i][j] y[j]
    std::vector<ExprPtr> y(n);
    for (std::size_t i = 0; i < n; ++i) {
        ExprPtr sum = integer(ctx, 0);
        for (std::size_t j = 0; j < i; ++j) {
            auto prod = mul_expr(ctx, lu.L(i, j), y[j]);
            if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
            auto next_sum = add_expr(ctx, sum, prod.value());
            if (next_sum.is_error()) return fail<std::vector<ExprPtr>>(next_sum.error());
            sum = next_sum.value();
        }
        auto yi = sub_expr(ctx, b[i], sum);
        if (yi.is_error()) return fail<std::vector<ExprPtr>>(yi.error());
        y[i] = yi.value();
    }
    // Back substitution: U·x = y. x[i] = (y[i] - Σ_{j>i} U[i][j] x[j]) / U[i][i]
    std::vector<ExprPtr> x(n);
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
        ExprPtr sum = integer(ctx, 0);
        for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n; ++j) {
            auto prod = mul_expr(ctx, lu.U(i, j), x[j]);
            if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
            auto next_sum = add_expr(ctx, sum, prod.value());
            if (next_sum.is_error()) return fail<std::vector<ExprPtr>>(next_sum.error());
            sum = next_sum.value();
        }
        auto num = sub_expr(ctx, y[i], sum);
        if (num.is_error()) return fail<std::vector<ExprPtr>>(num.error());
        auto quot = div_expr(ctx, num.value(), lu.U(i, i));
        if (quot.is_error()) return fail<std::vector<ExprPtr>>(quot.error());
        
        auto t = algebra::together(quot.value(), ctx);
        ExprPtr norm = t.is_ok() ? t.value() : quot.value();
        auto simp = ctx.simplify(norm);
        x[i] = simp.is_ok() ? simp.value() : norm;
    }
    return ok(x);
}

}  // namespace cas::linalg
