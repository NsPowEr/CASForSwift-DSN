// CAS-F4.1c — Test Cholesky LDL^T.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class LDLTTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] MatrixExpr from_rows(
        std::initializer_list<std::initializer_list<long long>> rows) {
        std::size_t r = rows.size();
        std::size_t c = rows.begin()->size();
        std::vector<ExprPtr> data;
        for (auto& row : rows) for (auto v : row) data.push_back(lit(v));
        return MatrixExpr(r, c, data);
    }
    [[nodiscard]] bool entries_equal(ExprPtr a, ExprPtr b) {
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
        auto t = algebra::together(delta, ctx);
        auto s = ctx.simplify(t.is_ok() ? t.value() : delta);
        if (s.is_error()) return false;
        if (auto* il = expr_cast<IntegerLit>(s.value())) return il->value.is_zero();
        if (auto* rl = expr_cast<RationalLit>(s.value())) return rl->numerator.is_zero();
        return false;
    }
};

TEST_F(LDLTTest, IdentityIsTrivial) {
    auto I = from_rows({{1, 0}, {0, 1}});
    auto r = cholesky_ldlt(I, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_TRUE(entries_equal(r.value().L(0, 0), lit(1)));
    EXPECT_TRUE(entries_equal(r.value().L(1, 1), lit(1)));
    EXPECT_TRUE(entries_equal(r.value().L(1, 0), lit(0)));
    ASSERT_EQ(r.value().D.size(), 2U);
    EXPECT_TRUE(entries_equal(r.value().D[0], lit(1)));
    EXPECT_TRUE(entries_equal(r.value().D[1], lit(1)));
}

TEST_F(LDLTTest, NonSymmetricRejected) {
    auto A = from_rows({{1, 2}, {3, 4}});
    auto r = cholesky_ldlt(A, ctx);
    EXPECT_TRUE(r.is_error());
}

TEST_F(LDLTTest, Certificator_LDLt_Reconstructs_3x3) {
    // SPD matrix: A = [[4,2,1],[2,5,3],[1,3,6]]
    auto A = from_rows({{4, 2, 1}, {2, 5, 3}, {1, 3, 6}});
    auto r = cholesky_ldlt(A, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;

    // Reconstruct A_recon = L · diag(D) · L^T
    const std::size_t n = 3;
    MatrixExpr D_mat(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            D_mat(i, j) = (i == j) ? r.value().D[i] : lit(0);
    auto LD = multiply(r.value().L, D_mat, ctx);
    ASSERT_TRUE(LD.is_ok());
    auto Lt = transpose(r.value().L);
    ASSERT_TRUE(Lt.is_ok());
    auto recon = multiply(LD.value(), Lt.value(), ctx);
    ASSERT_TRUE(recon.is_ok());

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            EXPECT_TRUE(entries_equal(recon.value()(i, j), A(i, j)))
                << "L·D·L^T mismatch at (" << i << "," << j << ")";
}

TEST_F(LDLTTest, ZeroPivotReported) {
    // Symmetric singular: A[0][0] = 0 → D[0] = 0 → fails
    auto A = from_rows({{0, 1}, {1, 1}});
    auto r = cholesky_ldlt(A, ctx);
    EXPECT_TRUE(r.is_error());
}

}  // namespace
