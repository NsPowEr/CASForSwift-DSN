// CAS-F4.1b — QR decomposition via Householder reflections (symbolic).
//
// A = Q · R, dove Q ha colonne ortonormali (Q^T·Q = I_n) e R è
// triangolare superiore.  Restituisce la factorizzazione "reduced":
//   Q : m×n
//   R : n×n
// con la convenzione segno: alpha_k = sgn(x_{0}^{(k)}) · ‖x^{(k)}‖.
//
// Algoritmo (Golub-Van Loan §5.1.3):
//   per ogni colonna k in [0, min(n, m−1)):
//     x = R_work[k:m, k]
//     alpha = sgn(x_0) · ‖x‖     (sgn=+ se sconosciuto)
//     u = x − alpha · e_1
//     se u = 0 → colonna già canonica, salta riflessione
//     H_k = I − 2 u u^T / (u^T u)  applicato in-place a R_work[k:m, k:n]
//     accumula H_k a destra in Q_full
//
// Singolarità: se durante la riflessione k-esima la norma del residuo è
// zero (colonna dipendente dalle precedenti), restituisce Unimplemented
// con reason code LINALG_LINEAR_DEPENDENT.

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

Result<QRDecomposition> qr_decompose(const MatrixExpr& matrix,
                                       symbolic::CASContext& ctx) {
    const std::size_t m = matrix.rows();
    const std::size_t n = matrix.cols();
    if (m == 0U || n == 0U) {
        return fail<QRDecomposition>(CASError{
            CASErrorKind::InvalidArgument, "qr_decompose: empty matrix", std::nullopt});
    }

    MatrixExpr R_work(m, n, matrix.elements());
    auto Qi = identity(m, ctx);
    if (Qi.is_error()) return fail<QRDecomposition>(Qi.error());
    MatrixExpr Q_full = std::move(Qi.value());

    for (std::size_t k = 0; k < std::min(n, m); ++k) {
        const std::size_t len = m - k;
        if (len == 0U) break;

        std::vector<ExprPtr> x(len);
        for (std::size_t i = 0; i < len; ++i) x[i] = R_work(k + i, k);

        auto norm_sq_res = sym_norm_sq(x, ctx);
        if (norm_sq_res.is_error()) return fail<QRDecomposition>(norm_sq_res.error());
        auto norm_sq_simp = simplify(ctx, norm_sq_res.value());
        ExprPtr norm_sq = norm_sq_simp.is_ok() ? norm_sq_simp.value() : norm_sq_res.value();
if (is_zero_expr(norm_sq)) {
    return make_unimplemented<QRDecomposition>(
        "linalg", "qr_decompose", "matrix column is linearly dependent",
        error::reason_codes::LINALG_LINEAR_DEPENDENT,
        "Matrix must have full column rank for this QR implementation",
        "F4.1b");
}

        auto norm_res = sym_norm(x, ctx);
        if (norm_res.is_error()) return fail<QRDecomposition>(norm_res.error());
        ExprPtr norm = norm_res.value();

        // alpha = sgn(x[0]) * norm.
        ExprPtr alpha = norm;
        if (ctx.assumptions().is_negative(x[0])) {
            auto neg_norm = negate_expr(ctx, norm);
            if (neg_norm.is_error()) return fail<QRDecomposition>(neg_norm.error());
            alpha = neg_norm.value();
        }

        // u = x - alpha * e_1
        std::vector<ExprPtr> u = x;
        auto u0_res = sub_expr(ctx, u[0], alpha);
        if (u0_res.is_error()) return fail<QRDecomposition>(u0_res.error());
        u[0] = u0_res.value();

        // denom = 2 * (norm^2 - alpha * x[0])
        auto ax0 = mul_expr(ctx, alpha, x[0]);
        if (ax0.is_error()) return fail<QRDecomposition>(ax0.error());
        auto diff = sub_expr(ctx, norm_sq, ax0.value());
        if (diff.is_error()) return fail<QRDecomposition>(diff.error());
        auto denom_res = mul_expr(ctx, integer(ctx, 2), diff.value());
        if (denom_res.is_error()) return fail<QRDecomposition>(denom_res.error());
        ExprPtr denom = denom_res.value();

        if (is_zero_expr(denom)) continue;

        // Apply Householder to R_work
        for (std::size_t j = k; j < n; ++j) {
            ExprPtr uv = integer(ctx, 0);
            for (std::size_t i = 0; i < len; ++i) {
                auto prod = mul_expr(ctx, u[i], R_work(k + i, j));
                if (prod.is_error()) return fail<QRDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, uv, prod.value());
                if (next_sum.is_error()) return fail<QRDecomposition>(next_sum.error());
                uv = next_sum.value();
            }
            auto num = mul_expr(ctx, integer(ctx, 2), uv);
            if (num.is_error()) return fail<QRDecomposition>(num.error());
            auto coeff_res = div_expr(ctx, num.value(), denom);
            if (coeff_res.is_error()) return fail<QRDecomposition>(coeff_res.error());
            ExprPtr coeff = coeff_res.value();

            if (is_zero_expr(coeff)) continue;
            for (std::size_t i = 0; i < len; ++i) {
                auto update = mul_expr(ctx, coeff, u[i]);
                if (update.is_error()) return fail<QRDecomposition>(update.error());
                auto next_r = sub_expr(ctx, R_work(k + i, j), update.value());
                if (next_r.is_error()) return fail<QRDecomposition>(next_r.error());
                auto simp_r = simplify(ctx, next_r.value());
                R_work(k + i, j) = simp_r.is_ok() ? simp_r.value() : next_r.value();
            }
        }

        // Apply Householder to Q_full (accumulate from right: Q = H1*H2*...*Hk)
        for (std::size_t j = 0; j < m; ++j) {
            ExprPtr uv = integer(ctx, 0);
            for (std::size_t i = 0; i < len; ++i) {
                auto prod = mul_expr(ctx, u[i], Q_full(j, k + i));
                if (prod.is_error()) return fail<QRDecomposition>(prod.error());
                auto next_sum = add_expr(ctx, uv, prod.value());
                if (next_sum.is_error()) return fail<QRDecomposition>(next_sum.error());
                uv = next_sum.value();
            }
            auto num = mul_expr(ctx, integer(ctx, 2), uv);
            if (num.is_error()) return fail<QRDecomposition>(num.error());
            auto coeff_res = div_expr(ctx, num.value(), denom);
            if (coeff_res.is_error()) return fail<QRDecomposition>(coeff_res.error());
            ExprPtr coeff = coeff_res.value();

            if (is_zero_expr(coeff)) continue;
            for (std::size_t i = 0; i < len; ++i) {
                auto update = mul_expr(ctx, coeff, u[i]);
                if (update.is_error()) return fail<QRDecomposition>(update.error());
                auto next_q = sub_expr(ctx, Q_full(j, k + i), update.value());
                if (next_q.is_error()) return fail<QRDecomposition>(next_q.error());
                auto simp_q = simplify(ctx, next_q.value());
                Q_full(j, k + i) = simp_q.is_ok() ? simp_q.value() : next_q.value();
            }
        }
    }

    // Reduced QR: return first n columns of Q and first n rows of R.
    MatrixExpr Q_red(m, n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) Q_red(i, j) = Q_full(i, j);

    MatrixExpr R_red(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) R_red(i, j) = (i <= j) ? R_work(i, j) : integer(ctx, 0);

    return ok(QRDecomposition{std::move(Q_red), std::move(R_red)});
}

}  // namespace cas::linalg
