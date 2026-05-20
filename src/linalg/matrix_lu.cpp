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

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::linalg {

namespace {

[[nodiscard]] ExprPtr make_zero(symbolic::CASContext& ctx) {
    return ctx.arena().make<IntegerLit>(BigInt(0));
}
[[nodiscard]] ExprPtr make_one(symbolic::CASContext& ctx) {
    return ctx.arena().make<IntegerLit>(BigInt(1));
}

[[nodiscard]] bool is_zero_expr(ExprPtr e) {
    if (!e) return true;
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

}  // namespace

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
            L(i, j) = (i == j) ? make_one(ctx) : make_zero(ctx);
            U(i, j) = make_zero(ctx);
        }
    }

    AstArena& arena = ctx.arena();
    for (std::size_t k = 0; k < n; ++k) {
        // Compute U[k][j] for j ≥ k.
        for (std::size_t j = k; j < n; ++j) {
            ExprPtr sum = make_zero(ctx);
            for (std::size_t s = 0; s < k; ++s) {
                ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, L(k, s), U(s, j));
                sum = arena.make<Binary>(BinaryOp::Add, sum, prod);
            }
            ExprPtr u_kj = arena.make<Binary>(BinaryOp::Sub, matrix(k, j), sum);
            auto simp = ctx.simplify(u_kj);
            U(k, j) = simp.is_ok() ? simp.value() : u_kj;
        }

        // Check pivot U[k][k] ≠ 0.
        if (is_zero_expr(U(k, k))) {
            return fail<LUDecomposition>(CASError{
                CASErrorKind::Unimplemented,
                "lu_decompose: zero pivot at row " + std::to_string(k)
                    + " — use permuted LU (follow-up)",
                std::nullopt});
        }

        // Compute L[i][k] for i > k.
        for (std::size_t i = k + 1; i < n; ++i) {
            ExprPtr sum = make_zero(ctx);
            for (std::size_t s = 0; s < k; ++s) {
                ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, L(i, s), U(s, k));
                sum = arena.make<Binary>(BinaryOp::Add, sum, prod);
            }
            ExprPtr num = arena.make<Binary>(BinaryOp::Sub, matrix(i, k), sum);
            ExprPtr quot = arena.make<Binary>(BinaryOp::Div, num, U(k, k));
            auto t = algebra::together(quot, ctx);
            ExprPtr norm = t.is_ok() ? t.value() : quot;
            auto simp = ctx.simplify(norm);
            L(i, k) = simp.is_ok() ? simp.value() : norm;
        }
    }

    return ok(LUDecomposition{std::move(L), std::move(U)});
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
    AstArena& arena = ctx.arena();
    // Forward substitution: L·y = b. L unit triangular → y[i] = b[i] - Σ L[i][j] y[j]
    std::vector<ExprPtr> y(n);
    for (std::size_t i = 0; i < n; ++i) {
        ExprPtr sum = arena.make<IntegerLit>(BigInt(0));
        for (std::size_t j = 0; j < i; ++j) {
            ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, lu.L(i, j), y[j]);
            sum = arena.make<Binary>(BinaryOp::Add, sum, prod);
        }
        ExprPtr yi = arena.make<Binary>(BinaryOp::Sub, b[i], sum);
        auto simp = ctx.simplify(yi);
        y[i] = simp.is_ok() ? simp.value() : yi;
    }
    // Back substitution: U·x = y. x[i] = (y[i] - Σ_{j>i} U[i][j] x[j]) / U[i][i]
    std::vector<ExprPtr> x(n);
    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i) {
        ExprPtr sum = arena.make<IntegerLit>(BigInt(0));
        for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n; ++j) {
            ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, lu.U(i, j), x[j]);
            sum = arena.make<Binary>(BinaryOp::Add, sum, prod);
        }
        ExprPtr num = arena.make<Binary>(BinaryOp::Sub, y[i], sum);
        ExprPtr quot = arena.make<Binary>(BinaryOp::Div, num, lu.U(i, i));
        auto t = algebra::together(quot, ctx);
        ExprPtr norm = t.is_ok() ? t.value() : quot;
        auto simp = ctx.simplify(norm);
        x[i] = simp.is_ok() ? simp.value() : norm;
    }
    return ok(x);
}

}  // namespace cas::linalg
