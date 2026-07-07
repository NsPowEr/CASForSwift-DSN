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

// Cofactor adjugate inverse: A⁻¹[i][j] = (-1)^(i+j) · det(M_{j,i}) / det(A)
// where M_{j,i} is the (j,i)-minor (matrix obtained by deleting row j and
// column i of A).  Costs n² calls to `determinant` on (n-1)×(n-1) sub-
// matrices, hence O(n²·T_det(n-1)).  Used as fallback when Bareiss-Jordan
// cannot decide a pivot on symbolic/RootOf inputs (the per-call
// `bareiss_determinant` handles small-n RootOf arithmetic via the
// Sylvester identity without requiring pivot certification at every
// elimination step).  Reference: HC-F4-JORDAN-INVERSE-ROOTOF closure.
[[nodiscard]] Result<MatrixExpr> inverse_via_cofactor(
    const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();

    auto det_res = determinant(matrix, ctx);
    if (det_res.is_error()) return fail<MatrixExpr>(det_res.error());
    ExprPtr det = det_res.value();
    if (is_zero_expr(det)) {
        return fail<MatrixExpr>(make_error(CASErrorKind::Undefined,
            "inverse: matrix is singular"));
    }

    MatrixExpr inv(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            // Build minor M_{j,i}: delete row j and column i.
            MatrixExpr minor(n - 1U, n - 1U);
            std::size_t rr = 0;
            for (std::size_t r = 0; r < n; ++r) {
                if (r == j) continue;
                std::size_t cc = 0;
                for (std::size_t c = 0; c < n; ++c) {
                    if (c == i) continue;
                    minor(rr, cc) = matrix(r, c);
                    ++cc;
                }
                ++rr;
            }
            auto m_det = determinant(minor, ctx);
            if (m_det.is_error()) return fail<MatrixExpr>(m_det.error());
            ExprPtr cof = m_det.value();
            if ((i + j) % 2U != 0U) {
                auto neg = negate_expr(ctx, cof);
                if (neg.is_error()) return fail<MatrixExpr>(neg.error());
                cof = neg.value();
            }
            auto q = div_expr(ctx, cof, det);
            if (q.is_error()) return fail<MatrixExpr>(q.error());
            auto t = algebra::together(q.value(), ctx);
            ExprPtr val = t.is_ok() ? t.value() : q.value();
            auto s = simplify(ctx, val);
            inv(i, j) = s.is_ok() ? s.value() : val;
        }
    }
    return ok(std::move(inv));
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
            }
        }
        if (pivot_row == n || best.certainty < 0) {
            // No pivot certified nonzero in column k. For matrices with
            // RootOf-extension entries this is typically NOT true
            // singularity: the residual M(i, k) is a polynomial in the
            // algebraic generators that needs reduction modulo their
            // minimal polynomials to decide nonzero. Signal Unimplemented
            // with the standard reason code so the outer `inverse` can
            // fall back to the cofactor-adjugate path on small matrices.
            return make_unimplemented<MatrixExpr>(
                "linalg", "inverse",
                "pivot expression cannot be decided nonzero",
                error::reason_codes::LINALG_RREF_UNDECIDABLE_PIVOT,
                "Add assumptions about pivot variables or use the cofactor "
                "fallback (n <= 5)",
                "F4.6");
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

    // Bareiss-Edmonds Gauss-Jordan: optimal for matrices with structurally-
    // certifiable pivots (rational/symbolic entries with `is_known_nonzero`
    // decidable).  Returns Unimplemented(LINALG_RREF_UNDECIDABLE_PIVOT) when
    // a residual pivot cannot be decided non-zero — typical for matrices
    // with RootOf-extension entries whose nonzero status requires reduction
    // modulo the minimal polynomial.  In that case fall back to the
    // cofactor-adjugate algorithm, which routes every minor through
    // `bareiss_determinant` (Sylvester identity) and tolerates undecidable
    // intermediates because the final determinant is computed in one shot.
    //
    // For small n (≤ 5) the cofactor fallback is asymptotically cheap
    // (≤ 25 calls to bareiss(4×4)); for larger n the BJ failure is more
    // likely structural (true singularity) than RootOf-undecidability, so
    // the cofactor pass is skipped to avoid O(n²·n³) = O(n^5) cost on
    // expensive symbolic inputs.
    auto bj_res = inverse_bareiss_jordan(matrix, ctx);
    if (bj_res.is_ok()) return bj_res;
    if (n <= 5U && bj_res.error().kind == CASErrorKind::Unimplemented) {
        return inverse_via_cofactor(matrix, ctx);
    }
    return bj_res;
}

} // namespace cas::linalg
