// CAS-F4.3 — Test determinanti speciali (Vandermonde, Toeplitz, Circulant, Banded).

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class SpecialDetTest : public ::testing::Test {
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
        auto eq = symbolic::mathematically_equal(a, b, ctx);
        return eq.is_ok() && eq.value();
    }
};

// Vandermonde 3×3 con x = (1, 2, 3) → det = (2−1)(3−1)(3−2) = 2.
TEST_F(SpecialDetTest, Vandermonde_3x3_Integer) {
    auto V = from_rows({{1, 1, 1}, {1, 2, 4}, {1, 3, 9}});
    auto d = determinant(V, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    EXPECT_TRUE(entries_equal(d.value(), lit(2)));
}

// Vandermonde 4×4 con x = (1, 2, 3, 4) → det = ∏(x_j-x_i) for i<j = 12.
TEST_F(SpecialDetTest, Vandermonde_4x4) {
    auto V = from_rows({{1, 1, 1, 1}, {1, 2, 4, 8}, {1, 3, 9, 27}, {1, 4, 16, 64}});
    auto d = determinant(V, ctx);
    ASSERT_TRUE(d.is_ok());
    EXPECT_TRUE(entries_equal(d.value(), lit(12)));
}

// Circulant n=2: M = [[a,b],[b,a]], det = a²−b² = (a−b)(a+b). Per (3,1) → 8.
TEST_F(SpecialDetTest, Circulant_2x2) {
    auto C = from_rows({{3, 1}, {1, 3}});
    auto d = determinant(C, ctx);
    ASSERT_TRUE(d.is_ok());
    EXPECT_TRUE(entries_equal(d.value(), lit(8)));
}

// Circulant n=3: M = [[1,2,3],[3,1,2],[2,3,1]]. det = 18.
TEST_F(SpecialDetTest, Circulant_3x3) {
    auto C = from_rows({{1, 2, 3}, {3, 1, 2}, {2, 3, 1}});
    auto d = determinant(C, ctx);
    ASSERT_TRUE(d.is_ok());
    EXPECT_TRUE(entries_equal(d.value(), lit(18)));
}

// Circulant n=4: M[0]=[1,2,3,4]. det = (10)(-2)((1-3)²+(2-4)²) = 10·(-2)·8 = -160.
TEST_F(SpecialDetTest, Circulant_4x4) {
    auto C = from_rows({{1, 2, 3, 4}, {4, 1, 2, 3}, {3, 4, 1, 2}, {2, 3, 4, 1}});
    auto d = determinant(C, ctx);
    ASSERT_TRUE(d.is_ok());
    EXPECT_TRUE(entries_equal(d.value(), lit(-160)));
}

// Toeplitz 3×3: M = [[2,3,5],[1,2,3],[4,1,2]]. Toeplitz (diag costanti).
// Bareiss applicato → det = computed directly.
TEST_F(SpecialDetTest, Toeplitz_3x3_DetectedAndCorrect) {
    auto T = from_rows({{2, 3, 5}, {1, 2, 3}, {4, 1, 2}});
    auto d = determinant(T, ctx);
    ASSERT_TRUE(d.is_ok());
    // det = 2(4-3) - 3(2-12) + 5(1-8) = 2 + 30 - 35 = -3.
    EXPECT_TRUE(entries_equal(d.value(), lit(-3)));
}

// Banded 5×5: route a Bareiss generale (la specializzazione "fraction-free
// banded" è follow-up).  Cert: det non-nullo e consistente con cofactor.
TEST_F(SpecialDetTest, Banded_5x5_PentadiagonalCorrect) {
    auto B = from_rows({
        {2, -1,  1,  0,  0},
        {-1, 2, -1,  1,  0},
        { 1, -1, 2, -1,  1},
        { 0,  1, -1, 2, -1},
        { 0,  0,  1, -1, 2},
    });
    auto d = determinant(B, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    // Espansione di Laplace lungo la riga 0 → riduzione a 4 minori 4×4
    // computati come det di sotto-matrici.  Confronto contro questo riferimento
    // ricorsivo (chiamata determinante stessa: il pattern bandato del minore
    // 4×4 lo route a Bareiss general anch'esso → comparabile).
    // Espressione 2·M00 + 1·M01 + 1·M02 (M0i = minor i di riga 0).
    auto minor = [&](std::initializer_list<std::size_t> rows,
                     std::initializer_list<std::size_t> cols) {
        std::vector<ExprPtr> data;
        for (auto r : rows) for (auto c : cols) data.push_back(B(r, c));
        return MatrixExpr(rows.size(), cols.size(), data);
    };
    auto m00 = determinant(minor({1,2,3,4}, {1,2,3,4}), ctx);
    auto m01 = determinant(minor({1,2,3,4}, {0,2,3,4}), ctx);
    auto m02 = determinant(minor({1,2,3,4}, {0,1,3,4}), ctx);
    ASSERT_TRUE(m00.is_ok() && m01.is_ok() && m02.is_ok());
    // det(B) = 2·m00 - (-1)·m01 + 1·m02 = 2·m00 + m01 + m02
    ExprPtr expected = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        ctx.arena().make<Binary>(BinaryOp::Mul, lit(2), m00.value()),
        m01.value(),
        m02.value()
    });
    EXPECT_TRUE(entries_equal(d.value(), expected))
        << "Banded routed to Bareiss general should equal cofactor expansion";
}

}  // namespace
