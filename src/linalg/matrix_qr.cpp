// CAS-F4.1b — QR decomposition via Modified Gram-Schmidt (symbolic).
//
// Decomposes A (m×n, m ≥ n, full column rank) into Q (m×n) e R (n×n)
// triangolare superiore tali che A = Q·R e Q^T·Q = I.
//
// Algoritmo (Trefethen-Bau "Numerical Linear Algebra" Lecture 8, Modified GS):
//   for k in [0, n):
//       v_k = A[:, k] - Σ_{j<k} r_{jk}·q_j        (rational pipeline)
//       N_k = v_k^T · v_k                          (rational)
//       R[k][k] = sqrt(N_k)                        (irrazionale)
//       Q[:, k] = v_k / sqrt(N_k)
//       for j in (k, n):
//           dot_{k,j} = v_k^T · v_j (correnti)    (rational)
//           R[k][j] = dot_{k,j} / sqrt(N_k)       = dot_{k,j} · sqrt(N_k) / N_k
//           v_j    -= (dot_{k,j} / N_k) · v_k     (rational update, NO sqrt)
//
// Vantaggio rispetto a Householder simbolico: gli aggiornamenti dei vettori
// rimanenti sono *razionali puri* (coeff = dot/N_k, no sqrt nelle entry V),
// quindi l'AST non esplode. Le radici quadrate appaiono solo nelle colonne
// di Q e nella diagonale di R. La cert Q^T·Q == I e Q·R == A semplifica
// correttamente grazie alle regole sqrt(a)·sqrt(a) → a nel simplifier.
//
// Singolarità: se N_k = 0 (colonna linearmente dipendente dalle precedenti),
// restituisce Unimplemented con reason code LINALG_LINEAR_DEPENDENT.

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

    // V is the running matrix of partially-orthogonalized columns (rational
    // pipeline). Initialized to A.
    MatrixExpr V(m, n, matrix.elements());
    MatrixExpr Q(m, n);
    MatrixExpr R(n, n);
    R.fill(integer(ctx, 0));

    for (std::size_t k = 0; k < n; ++k) {
        // Compute N_k = V[:, k]^T · V[:, k] (rational if A rational).
        std::vector<ExprPtr> vk(m);
        for (std::size_t i = 0; i < m; ++i) vk[i] = V(i, k);
        auto norm_sq_res = sym_norm_sq(vk, ctx);
        if (norm_sq_res.is_error()) return fail<QRDecomposition>(norm_sq_res.error());
        ExprPtr norm_sq = norm_sq_res.value();

        if (is_zero_expr(norm_sq)) {
            return make_unimplemented<QRDecomposition>(
                "linalg", "qr_decompose", "matrix column is linearly dependent",
                error::reason_codes::LINALG_LINEAR_DEPENDENT,
                "Matrix must have full column rank for this QR implementation",
                "F4.1b");
        }

        // R[k][k] = sqrt(N_k). Left unevaluated; the simplifier will reduce
        // sqrt(N_k)·sqrt(N_k) to N_k where it appears (orthogonality cert).
        ExprPtr norm = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{norm_sq});
        R(k, k) = norm;

        // Q[:, k] = V[:, k] / sqrt(N_k).
        for (std::size_t i = 0; i < m; ++i) {
            auto q_res = div_expr(ctx, V(i, k), norm);
            if (q_res.is_error()) return fail<QRDecomposition>(q_res.error());
            Q(i, k) = q_res.value();
        }

        // For j > k:
        //   dot_{k,j} = V[:, k]^T · V[:, j]            (rational)
        //   R[k][j]   = dot_{k,j} / sqrt(N_k)
        //   V[:, j]  -= (dot_{k,j} / N_k) · V[:, k]    (rational, NO sqrt)
        for (std::size_t j = k + 1; j < n; ++j) {
            ExprPtr dot = integer(ctx, 0);
            for (std::size_t i = 0; i < m; ++i) {
                auto p = mul_expr(ctx, V(i, k), V(i, j));
                if (p.is_error()) return fail<QRDecomposition>(p.error());
                auto s = add_expr(ctx, dot, p.value());
                if (s.is_error()) return fail<QRDecomposition>(s.error());
                dot = s.value();
            }

            // R[k][j] = dot / sqrt(N_k).
            auto rkj = div_expr(ctx, dot, norm);
            if (rkj.is_error()) return fail<QRDecomposition>(rkj.error());
            R(k, j) = rkj.value();

            // coeff = dot / N_k (rational).
            auto coeff_res = div_expr(ctx, dot, norm_sq);
            if (coeff_res.is_error()) return fail<QRDecomposition>(coeff_res.error());
            ExprPtr coeff = coeff_res.value();
            if (is_zero_expr(coeff)) continue;

            for (std::size_t i = 0; i < m; ++i) {
                auto t = mul_expr(ctx, coeff, V(i, k));
                if (t.is_error()) return fail<QRDecomposition>(t.error());
                auto u = sub_expr(ctx, V(i, j), t.value());
                if (u.is_error()) return fail<QRDecomposition>(u.error());
                V(i, j) = u.value();
            }
        }
    }

    return ok(QRDecomposition{std::move(Q), std::move(R)});
}

}  // namespace cas::linalg
