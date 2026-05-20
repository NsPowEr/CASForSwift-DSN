// CAS-L3-17 — Symbolic LU decomposition tests.
//
// Verifies A = L · U structural invariant + L unit-triangular +
// U upper-triangular shape.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class LUTest : public ::testing::Test {
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
        return false;
    }
    [[nodiscard]] bool matrix_equal(const MatrixExpr& a, const MatrixExpr& b) {
        if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
        for (std::size_t i = 0; i < a.rows(); ++i)
            for (std::size_t j = 0; j < a.cols(); ++j)
                if (!entries_equal(a(i, j), b(i, j))) return false;
        return true;
    }
};

TEST_F(LUTest, IdentityDecomposes) {
    auto I = from_rows({{1, 0}, {0, 1}});
    auto r = lu_decompose(I, ctx);
    ASSERT_TRUE(r.is_ok());
    // L = I, U = I.
    EXPECT_TRUE(matrix_equal(r.value().L, I));
    EXPECT_TRUE(matrix_equal(r.value().U, I));
}

TEST_F(LUTest, ProductReconstructs2x2) {
    // A = [[2, 3], [4, 7]] → L[[1,0],[2,1]], U[[2,3],[0,1]]
    auto A = from_rows({{2, 3}, {4, 7}});
    auto r = lu_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto prod = multiply(r.value().L, r.value().U, ctx);
    ASSERT_TRUE(prod.is_ok());
    EXPECT_TRUE(matrix_equal(prod.value(), A))
        << "L · U should reconstruct A";
}

TEST_F(LUTest, ProductReconstructs3x3) {
    // Use a non-singular 3x3.
    auto A = from_rows({{2, 1, 1}, {4, 3, 3}, {8, 7, 9}});
    auto r = lu_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto prod = multiply(r.value().L, r.value().U, ctx);
    ASSERT_TRUE(prod.is_ok());
    EXPECT_TRUE(matrix_equal(prod.value(), A));
}

TEST_F(LUTest, LIsUnitTriangular) {
    auto A = from_rows({{2, 3}, {4, 7}});
    auto r = lu_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto& L = r.value().L;
    // Diagonal = 1
    for (std::size_t i = 0; i < L.rows(); ++i) {
        auto* il = expr_cast<IntegerLit>(L(i, i));
        ASSERT_NE(il, nullptr) << "L[" << i << "][" << i << "] not int";
        EXPECT_EQ(il->value, BigInt(1));
    }
    // Above diagonal = 0
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t j = i + 1; j < L.cols(); ++j) {
            auto* il = expr_cast<IntegerLit>(L(i, j));
            ASSERT_NE(il, nullptr);
            EXPECT_TRUE(il->value.is_zero());
        }
}

TEST_F(LUTest, UIsUpperTriangular) {
    auto A = from_rows({{2, 3}, {4, 7}});
    auto r = lu_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto& U = r.value().U;
    // Below diagonal = 0
    for (std::size_t i = 1; i < U.rows(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            auto* il = expr_cast<IntegerLit>(U(i, j));
            ASSERT_NE(il, nullptr);
            EXPECT_TRUE(il->value.is_zero());
        }
}

TEST_F(LUTest, ZeroPivotRejected) {
    // A = [[0, 1], [1, 0]] — A[0][0] = 0 → no LU without pivoting.
    auto A = from_rows({{0, 1}, {1, 0}});
    auto r = lu_decompose(A, ctx);
    EXPECT_TRUE(r.is_error()) << "Zero pivot should fail without pivoting";
}

TEST_F(LUTest, AntiHardcodeNonSingular4x4) {
    auto A = from_rows({{1, 2, 3, 4},
                        {2, 5, 7, 8},
                        {3, 7, 11, 13},
                        {4, 8, 13, 19}});
    auto r = lu_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto prod = multiply(r.value().L, r.value().U, ctx);
    ASSERT_TRUE(prod.is_ok());
    EXPECT_TRUE(matrix_equal(prod.value(), A));
}

}  // namespace
