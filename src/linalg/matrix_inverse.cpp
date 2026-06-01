// CAS-F4.6 — Symbolic matrix inverse via fraction-free Gauss-Jordan
// (Bareiss-Edmonds; Geddes/Czapor/Labahn "Algorithms for Computer Algebra"
// §9.5, Algorithm 9.2). For n ≥ 3 the augmented matrix [A | I] is reduced in
// place with the Bareiss step
//
//   M[i][j] ← (M[k][k] · M[i][j] − M[i][k] · M[k][j]) / d_{k-1}
//
// applied to all rows i ≠ k (full Gauss-Jordan), with d_{-1} = 1.  By the
// Sylvester identity the division is exact in the integral closure of the
// original entries, so intermediate expressions stay polynomial — no
// rational blow-up.  At termination the bottom-right pivot M(n-1, n-1)
// equals det(A) (up to row-swap sign) and the right block C is the
// adjugate up to the same sign; hence A^{-1} = C / det(A) is extracted
// with a single division.  The 2×2 fast path is preserved for brevity
// and to avoid the overhead of the augmented sweep on the trivial case.
//
// Canonical extraction requires the simplifier to flatten
// Pow(Product, n_int) into Product(Pow(factor_i, n_int)) so the final
// `C[i][j] / det(A)` cross-cancels Pow terms on shared symbols.  This
// rule is enabled in simplify_arithmetic.cpp::simplify_power for any
// integer n with |n| ≤ 20.  Reference: HC-F4-INV-SYMBOLIC-CANONICAL.

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"

#include <optional>
#include <string>
#include <utility>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] Result<MatrixExpr> inverse_bareiss_jordan(
    const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    const std::size_t total_cols = 2U * n;

    // Augmented matrix M = [A | I] of size n × 2n.
    MatrixExpr M(n, total_cols);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) M(i, j) = matrix(i, j);
        for (std::size_t j = 0; j < n; ++j) {
            M(i, n + j) = (i == j) ? integer(ctx, 1) : integer(ctx, 0);
        }
    }

    ExprPtr d_prev = integer(ctx, 1);

    for (std::size_t k = 0; k < n; ++k) {
        // Pivot selection via PivotScore: literal-nonzero beats
        // structurally-nonzero, lower degree/complexity preferred.
        std::size_t pivot_row = n;
        PivotScore best{-1, 0, 0};
        for (std::size_t i = k; i < n; ++i) {
            if (is_zero_expr(M(i, k))) continue;
            PivotScore s = make_pivot_score(M(i, k), ctx);
            if (pivot_row == n || s > best) {
                best = s;
                pivot_row = i;
                if (s.certainty == 3) break;
            }
        }
        if (pivot_row == n || best.certainty < 0) {
            return fail<MatrixExpr>(make_error(CASErrorKind::Undefined,
                "inverse: matrix is singular (or pivot cannot be decided nonzero)"));
        }

        if (pivot_row != k) {
            for (std::size_t j = 0; j < total_cols; ++j) {
                std::swap(M(k, j), M(pivot_row, j));
            }
            // Row swap negates det(A) and permutes the right block in a
            // matching way, so the relation C · A = M(n-1, n-1) · I (carrying
            // the same sign) is preserved and A^{-1} = C / M(n-1, n-1)
            // remains correct without explicit sign correction.
        }

        ExprPtr pivot = M(k, k);

        // Bareiss step: for rows i ≠ k, columns j > k, update via
        //   M[i][j] = (pivot · M[i][j] − M[i][k] · M[k][j]) / d_prev
        // and zero column k.
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            ExprPtr m_ik = M(i, k);
            const bool m_ik_zero = is_zero_expr(m_ik);
            for (std::size_t j = k + 1U; j < total_cols; ++j) {
                ExprPtr new_val;
                if (m_ik_zero) {
                    auto lhs = mul_expr(ctx, pivot, M(i, j));
                    if (lhs.is_error()) return fail<MatrixExpr>(lhs.error());
                    auto val = div_expr(ctx, lhs.value(), d_prev);
                    if (val.is_error()) return fail<MatrixExpr>(val.error());
                    new_val = val.value();
                } else {
                    auto lhs = mul_expr(ctx, pivot, M(i, j));
                    if (lhs.is_error()) return fail<MatrixExpr>(lhs.error());
                    auto rhs = mul_expr(ctx, m_ik, M(k, j));
                    if (rhs.is_error()) return fail<MatrixExpr>(rhs.error());
                    auto diff = sub_expr(ctx, lhs.value(), rhs.value());
                    if (diff.is_error()) return fail<MatrixExpr>(diff.error());
                    auto val = div_expr(ctx, diff.value(), d_prev);
                    if (val.is_error()) return fail<MatrixExpr>(val.error());
                    new_val = val.value();
                }
                M(i, j) = new_val;
            }
            M(i, k) = integer(ctx, 0);
        }

        d_prev = pivot;
    }

    // A^{-1}[i][j] = M(i, n+j) / det where det = M(n-1, n-1).
    ExprPtr det = M(n - 1U, n - 1U);
    if (is_zero_expr(det)) {
        return fail<MatrixExpr>(make_error(CASErrorKind::Undefined,
            "inverse: matrix is singular"));
    }

    MatrixExpr inv(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            auto q = div_expr(ctx, M(i, n + j), det);
            if (q.is_error()) return fail<MatrixExpr>(q.error());
            auto t = algebra::together(q.value(), ctx);
            ExprPtr val = t.is_ok() ? t.value() : q.value();
            auto s = simplify(ctx, val);
            inv(i, j) = s.is_ok() ? s.value() : val;
        }
    }
    return ok(std::move(inv));
}

} // namespace

Result<MatrixExpr> inverse(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<MatrixExpr>(make_error(CASErrorKind::InvalidArgument, "Inverse requires a square matrix"));
    }
    if (n == 0) {
        return fail<MatrixExpr>(make_error(CASErrorKind::InvalidArgument, "Inverse of empty matrix"));
    }

    if (n == 2) {
        auto det_res = determinant(matrix, ctx);
        if (det_res.is_error()) return fail<MatrixExpr>(det_res.error());
        ExprPtr det = det_res.value();
        if (is_zero_expr(det)) return fail<MatrixExpr>(make_error(CASErrorKind::Undefined, "inverse: matrix is singular"));

        MatrixExpr res(2, 2);
        auto d_res = div_expr(ctx, matrix(1, 1), det);
        auto mb_r = negate_expr(ctx, matrix(0, 1));
        auto mc_r = negate_expr(ctx, matrix(1, 0));
        if (mb_r.is_error() || mc_r.is_error()) return fail<MatrixExpr>(make_error(CASErrorKind::InternalError, "negate failed"));
        
        auto mb_res = div_expr(ctx, mb_r.value(), det);
        auto mc_res = div_expr(ctx, mc_r.value(), det);
        auto a_res = div_expr(ctx, matrix(0, 0), det);
        
        if (d_res.is_error() || mb_res.is_error() || mc_res.is_error() || a_res.is_error())
            return fail<MatrixExpr>(make_error(CASErrorKind::InternalError, "2x2 inverse failed"));
            
        res(0, 0) = d_res.value();
        res(0, 1) = mb_res.value();
        res(1, 0) = mc_res.value();
        res(1, 1) = a_res.value();
        
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                auto s = simplify(ctx, res(i, j));
                if (s.is_ok()) res(i, j) = s.value();
            }
        }
        return ok(std::move(res));
    }

    return inverse_bareiss_jordan(matrix, ctx);
}

} // namespace cas::linalg
