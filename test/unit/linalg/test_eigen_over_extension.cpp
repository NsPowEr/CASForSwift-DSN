// Tests for null_space_over_extension and the RootOf-aware eigenvector path.
// Validates that for matrices whose characteristic polynomial has irrational
// roots, the computed eigenvectors actually satisfy A*v = lambda*v over the
// algebraic extension Q(lambda).

#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast_debug.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <vector>

namespace cas::test {

class EigenOverExtensionTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    [[nodiscard]] ExprPtr I(long long n) { return ctx->arena().make<IntegerLit>(BigInt(n)); }

    // Build the n-by-n companion matrix for x^n - c (rational c).
    [[nodiscard]] linalg::MatrixExpr companion_xn_minus_c(std::size_t n, long long c) {
        linalg::MatrixExpr M(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                M(i, j) = I(0);
        M(0, n - 1U) = I(c);
        for (std::size_t i = 1; i < n; ++i) M(i, i - 1U) = I(1);
        return M;
    }

    // Verify a candidate eigenpair (lambda, v) satisfies A*v - lambda*v == 0
    // componentwise, where the residual is simplified in Q(lambda).
    [[nodiscard]] bool residual_is_zero(
        const linalg::MatrixExpr& A,
        const ExprPtr lambda,
        const std::vector<ExprPtr>& v) {
        const std::size_t n = A.rows();
        for (std::size_t i = 0; i < n; ++i) {
            // (A v)_i = sum_j A(i,j) * v_j
            ExprPtr accum = I(0);
            for (std::size_t j = 0; j < n; ++j) {
                ExprPtr prod = ctx->arena().make<Binary>(BinaryOp::Mul, A(i, j), v[j]);
                accum = ctx->arena().make<Binary>(BinaryOp::Add, accum, prod);
            }
            // residual_i = (A v)_i - lambda * v_i
            ExprPtr lam_vi = ctx->arena().make<Binary>(BinaryOp::Mul, lambda, v[i]);
            ExprPtr residual = ctx->arena().make<Binary>(BinaryOp::Sub, accum, lam_vi);

            auto reduced = cas::algebra::simplify_in_q_alpha(residual, *ctx);
            if (reduced.is_error()) return false;
            const auto* il = expr_cast<IntegerLit>(reduced.value());
            if (!il || !il->value.is_zero()) {
                ADD_FAILURE() << "Component " << i << " residual non-zero: "
                              << debug_print(reduced.value());
                return false;
            }
        }
        return true;
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

// 5x5 companion of x^5 - 2. Char poly = x^5 - 2.
// Degree 5 has no general radical formula, so solve_polynomial returns RootOf
// expressions, which exercises the algebraic-extension eigenvector path.
TEST_F(EigenOverExtensionTest, QuinticCompanionEigenvectorsSatisfyAv) {
    auto A = companion_xn_minus_c(5U, 2);

    auto eigs = linalg::eigenvectors(A, *ctx);
    ASSERT_TRUE(eigs.is_ok()) << "eigenvectors failed: " << eigs.error().message;
    ASSERT_GT(eigs.value().size(), 0U) << "Expected non-empty eigenpair list";

    int verified = 0;
    for (const auto& pair : eigs.value()) {
        EXPECT_TRUE(expr_is<RootOf>(pair.eigenvalue))
            << "Expected RootOf eigenvalue; got: " << debug_print(pair.eigenvalue);
        EXPECT_EQ(pair.eigenvector.size(), 5U);
        EXPECT_TRUE(residual_is_zero(A, pair.eigenvalue, pair.eigenvector));
        ++verified;
    }
    EXPECT_GE(verified, 1) << "At least one eigenpair must verify A*v = lambda*v";
}

// 4x4 companion of x^4 - 2. Same idea, degree 4.
TEST_F(EigenOverExtensionTest, QuarticCompanionEigenvectorsSatisfyAv) {
    auto A = companion_xn_minus_c(4U, 2);

    auto eigs = linalg::eigenvectors(A, *ctx);
    ASSERT_TRUE(eigs.is_ok()) << "eigenvectors failed: " << eigs.error().message;
    ASSERT_GT(eigs.value().size(), 0U) << "Expected non-empty eigenpair list";

    int verified = 0;
    for (const auto& pair : eigs.value()) {
        EXPECT_TRUE(expr_is<RootOf>(pair.eigenvalue))
            << "Expected RootOf eigenvalue; got: " << debug_print(pair.eigenvalue);
        EXPECT_EQ(pair.eigenvector.size(), 4U);
        EXPECT_TRUE(residual_is_zero(A, pair.eigenvalue, pair.eigenvector));
        ++verified;
    }
    EXPECT_GE(verified, 1);
}

// Direct test for null_space_over_extension on a hand-built shifted matrix.
// A = [[0, 0, 2], [1, 0, 0], [0, 1, 0]], alpha = RootOf(x^3 - 2, x, 0).
// shifted = A - alpha*I.  Expected null space dimension over Q(alpha) is 1.
TEST_F(EigenOverExtensionTest, NullSpaceOverExtensionRank2Has1DimensionalKernel) {
    auto A = companion_xn_minus_c(3U, 2);

    // Build alpha = RootOf(x^3 - 2, x, 0).
    ExprPtr x_sym = ctx->arena().make<Symbol>("x");
    ExprPtr x3 = ctx->arena().make<Binary>(BinaryOp::Pow, x_sym, I(3));
    ExprPtr poly = ctx->arena().make<Binary>(BinaryOp::Sub, x3, I(2));
    ExprPtr alpha = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);

    // Build shifted = A - alpha*I.
    linalg::MatrixExpr shifted(3U, 3U);
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            ExprPtr entry = A(r, c);
            if (r == c) {
                entry = ctx->arena().make<Binary>(BinaryOp::Sub, entry, alpha);
            }
            auto simplified = ctx->simplify(entry);
            ASSERT_TRUE(simplified.is_ok());
            shifted(r, c) = simplified.value();
        }
    }

    auto kernel = linalg::null_space_over_extension(shifted, alpha, *ctx);
    ASSERT_TRUE(kernel.is_ok()) << "null_space_over_extension failed: "
                                << kernel.error().message;
    ASSERT_EQ(kernel.value().size(), 1U)
        << "Expected exactly one kernel basis vector for rank-2 shifted matrix";

    // Verify the recovered vector v satisfies shifted * v = 0  (i.e. A v = alpha v).
    EXPECT_TRUE(residual_is_zero(A, alpha, kernel.value()[0]));
}

// Regression: matrix not in Q(alpha) — gracefully falls back to null_space.
TEST_F(EigenOverExtensionTest, NullSpaceOverExtensionFallsBackWhenEntriesNotInQAlpha) {
    // Build a matrix whose entries contain a free symbol y (not in Q(alpha)).
    linalg::MatrixExpr A(2U, 2U);
    A(0, 0) = I(1);
    A(0, 1) = ctx->arena().make<Symbol>("y");
    A(1, 0) = I(0);
    A(1, 1) = I(1);

    ExprPtr x_sym = ctx->arena().make<Symbol>("x");
    ExprPtr x2 = ctx->arena().make<Binary>(BinaryOp::Pow, x_sym, I(2));
    ExprPtr poly = ctx->arena().make<Binary>(BinaryOp::Sub, x2, I(2));
    ExprPtr alpha = ctx->arena().make<RootOf>(poly, Symbol("x"), 0U);

    // Should fall back without erroring out.
    auto kernel = linalg::null_space_over_extension(A, alpha, *ctx);
    ASSERT_TRUE(kernel.is_ok()) << "Expected graceful fallback; got error: "
                                << kernel.error().message;
    // Identity-ish 2x2 matrix has trivial null space.
    EXPECT_EQ(kernel.value().size(), 0U);
}

}  // namespace cas::test
