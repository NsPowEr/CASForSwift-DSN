// CAS-F4.4b — Test certificatori Jordan: A ≡ P·J·P^{-1}.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class JordanCertTest : public ::testing::Test {
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
    [[nodiscard]] Result<MatrixExpr> reconstruct(const JordanDecomposition& jd) {
        auto Pinv = inverse(jd.P, ctx);
        if (Pinv.is_error()) return fail<MatrixExpr>(Pinv.error());
        auto PJ = multiply(jd.P, jd.J, ctx);
        if (PJ.is_error()) return PJ;
        return multiply(PJ.value(), Pinv.value(), ctx);
    }
};

// Diagonalizzabile, eigenvalues rational distinti.
TEST_F(JordanCertTest, Diagonalizable_3x3_Rational) {
    auto A = from_rows({{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}});
    auto jd = jordan_normal_form(A, ctx);
    ASSERT_TRUE(jd.is_ok()) << jd.error().message;
    auto recon = reconstruct(jd.value());
    ASSERT_TRUE(recon.is_ok());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_TRUE(entries_equal(recon.value()(i, j), A(i, j)))
                << "P·J·P^-1 mismatch at (" << i << "," << j << ")";
}

// Eigenvalues RootOf (irrazionali): A = [[0,2],[1,0]], char poly λ²−2 → ±√2.
TEST_F(JordanCertTest, RootOf_Eigenvalues_2x2_Sqrt2) {
    auto A = from_rows({{0, 2}, {1, 0}});
    auto jd = jordan_normal_form(A, ctx);
    ASSERT_TRUE(jd.is_ok()) << jd.error().message;
    auto recon = reconstruct(jd.value());
    ASSERT_TRUE(recon.is_ok()) << recon.error().message;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(recon.value()(i, j), A(i, j)))
                << "P·J·P^-1 mismatch at (" << i << "," << j << ")";
}

// Block 2x2 con autovalore doppio (Jordan non-trivial).
TEST_F(JordanCertTest, NonDiagonalizable_2x2_Double) {
    // A = [[3, 1], [0, 3]] → eigenvalue 3 doppio, Jordan = [[3,1],[0,3]]
    auto A = from_rows({{3, 1}, {0, 3}});
    auto jd = jordan_normal_form(A, ctx);
    ASSERT_TRUE(jd.is_ok()) << jd.error().message;
    auto recon = reconstruct(jd.value());
    ASSERT_TRUE(recon.is_ok());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(recon.value()(i, j), A(i, j)));
}

// F4.4b — Companion matrix di (x²−2)² ha autovalori RootOf ±√2 con
// molteplicità 2 → Jordan deve costruire catene size 2 per ciascuno.
TEST_F(JordanCertTest, RootOf_Multiplicity2_CompanionDeg4) {
    // p(x) = x^4 - 4x^2 + 4 = (x^2 - 2)^2
    auto A = from_rows({{0, 0, 0, -4}, {1, 0, 0, 0}, {0, 1, 0, 4}, {0, 0, 1, 0}});
    auto jd = jordan_normal_form(A, ctx);
    ASSERT_TRUE(jd.is_ok()) << jd.error().message;
    auto recon = reconstruct(jd.value());
    ASSERT_TRUE(recon.is_ok()) << recon.error().message;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            EXPECT_TRUE(entries_equal(recon.value()(i, j), A(i, j)))
                << "P·J·P^-1 mismatch at (" << i << "," << j << ")";
}

}  // namespace
