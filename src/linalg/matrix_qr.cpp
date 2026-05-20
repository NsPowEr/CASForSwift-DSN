// CAS-L3-17 — QR decomposition via classical Gram-Schmidt (symbolic).
//
// A = Q · R, Q orthonormal columns, R upper triangular.
//
// Algoritmo:
//   q_1 = a_1 / ||a_1||
//   r_{1,1} = ||a_1||
//   for k = 2..n:
//     v_k = a_k - Σ_{j<k} (q_j · a_k) · q_j
//     r_{k,k} = ||v_k||
//     q_k = v_k / r_{k,k}
//     r_{j,k} = q_j · a_k  for j < k
//
// Norma simbolica: sqrt(Σ entry²).

#include "cas/linalg/Matrix.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::linalg {

namespace {

[[nodiscard]] ExprPtr zero_e(symbolic::CASContext& ctx) {
    return ctx.arena().make<IntegerLit>(BigInt(0));
}

// Inner product of two column vectors (length n).
[[nodiscard]] ExprPtr inner_product(const std::vector<ExprPtr>& a,
                                     const std::vector<ExprPtr>& b,
                                     symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr sum = zero_e(ctx);
    for (std::size_t i = 0; i < a.size(); ++i) {
        ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, a[i], b[i]);
        sum = arena.make<Binary>(BinaryOp::Add, sum, prod);
    }
    auto s = ctx.simplify(sum);
    return s.is_ok() ? s.value() : sum;
}

// Norm of column = sqrt(<v, v>).
[[nodiscard]] ExprPtr norm(const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    ExprPtr ip = inner_product(v, v, ctx);
    ExprPtr sqrt_ip = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{ip});
    auto s = ctx.simplify(sqrt_ip);
    return s.is_ok() ? s.value() : sqrt_ip;
}

}  // namespace

Result<QRDecomposition> qr_decompose(const MatrixExpr& matrix,
                                       symbolic::CASContext& ctx) {
    const std::size_t m = matrix.rows();
    const std::size_t n = matrix.cols();
    if (m < n) {
        return fail<QRDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "qr_decompose: requires rows ≥ cols (tall or square)", std::nullopt});
    }
    if (n == 0U) {
        return fail<QRDecomposition>(CASError{
            CASErrorKind::InvalidArgument,
            "qr_decompose: empty matrix", std::nullopt});
    }

    AstArena& arena = ctx.arena();
    MatrixExpr Q(m, n);
    MatrixExpr R(n, n);
    // Initialize.
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) Q(i, j) = zero_e(ctx);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) R(i, j) = zero_e(ctx);

    // Extract column k as vector.
    auto column = [&](const MatrixExpr& M, std::size_t k) {
        std::vector<ExprPtr> col(M.rows());
        for (std::size_t i = 0; i < M.rows(); ++i) col[i] = M(i, k);
        return col;
    };

    std::vector<std::vector<ExprPtr>> q_cols;
    for (std::size_t k = 0; k < n; ++k) {
        std::vector<ExprPtr> a_k = column(matrix, k);
        std::vector<ExprPtr> v_k = a_k;
        // Subtract projections onto previous q_j.
        for (std::size_t j = 0; j < k; ++j) {
            ExprPtr proj = inner_product(q_cols[j], a_k, ctx);
            R(j, k) = proj;
            for (std::size_t i = 0; i < m; ++i) {
                ExprPtr scaled = arena.make<Binary>(BinaryOp::Mul, proj, q_cols[j][i]);
                v_k[i] = arena.make<Binary>(BinaryOp::Sub, v_k[i], scaled);
                auto s = ctx.simplify(v_k[i]);
                if (s.is_ok()) v_k[i] = s.value();
            }
        }
        ExprPtr norm_v = norm(v_k, ctx);
        // Check non-zero norm (linear independence).
        if (auto* il = expr_cast<IntegerLit>(norm_v); il && il->value.is_zero()) {
            return fail<QRDecomposition>(CASError{
                CASErrorKind::Unimplemented,
                "qr_decompose: column " + std::to_string(k)
                    + " linearly dependent (zero norm)", std::nullopt});
        }
        R(k, k) = norm_v;
        // q_k = v_k / norm_v
        std::vector<ExprPtr> q_k(m);
        for (std::size_t i = 0; i < m; ++i) {
            ExprPtr q_i = arena.make<Binary>(BinaryOp::Div, v_k[i], norm_v);
            auto t = algebra::together(q_i, ctx);
            ExprPtr norm_q = t.is_ok() ? t.value() : q_i;
            auto s = ctx.simplify(norm_q);
            q_k[i] = s.is_ok() ? s.value() : norm_q;
        }
        q_cols.push_back(q_k);
        for (std::size_t i = 0; i < m; ++i) Q(i, k) = q_k[i];
    }
    return ok(QRDecomposition{std::move(Q), std::move(R)});
}

}  // namespace cas::linalg
