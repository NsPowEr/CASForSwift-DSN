// CAS-L3-17 — QR decomposition tests.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class QRTest : public ::testing::Test {
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
};

TEST_F(QRTest, IdentityQRIsIdentity) {
    auto I = from_rows({{1, 0}, {0, 1}});
    auto r = qr_decompose(I, ctx);
    ASSERT_TRUE(r.is_ok());
    // Q = I, R = I.
    EXPECT_TRUE(entries_equal(r.value().Q(0, 0), lit(1)));
    EXPECT_TRUE(entries_equal(r.value().Q(1, 1), lit(1)));
    EXPECT_TRUE(entries_equal(r.value().Q(0, 1), lit(0)));
    EXPECT_TRUE(entries_equal(r.value().R(0, 0), lit(1)));
}

TEST_F(QRTest, ProductReconstructs2x2) {
    auto A = from_rows({{1, 2}, {0, 1}});
    auto r = qr_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    auto prod = multiply(r.value().Q, r.value().R, ctx);
    ASSERT_TRUE(prod.is_ok());
    // Q·R must equal A.
    EXPECT_TRUE(entries_equal(prod.value()(0, 0), A(0, 0)));
    EXPECT_TRUE(entries_equal(prod.value()(0, 1), A(0, 1)));
    EXPECT_TRUE(entries_equal(prod.value()(1, 0), A(1, 0)));
    EXPECT_TRUE(entries_equal(prod.value()(1, 1), A(1, 1)));
}

TEST_F(QRTest, RIsUpperTriangular) {
    auto A = from_rows({{2, 1}, {1, 3}});
    auto r = qr_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());
    // R[1][0] must be zero.
    EXPECT_TRUE(entries_equal(r.value().R(1, 0), lit(0)));
}

TEST_F(QRTest, SingularMatrixRejected) {
    // Linearly dependent columns.
    auto A = from_rows({{1, 2}, {2, 4}});
    auto r = qr_decompose(A, ctx);
    EXPECT_TRUE(r.is_error())
        << "Linearly dependent columns should fail";
}

TEST_F(QRTest, SymbolicQR_PositiveUpperTriangular_QR_Reconstructs) {
    // M-FIX QR symbolic coverage: upper-triangular con a,c > 0.
    // ‖(a,0)‖ = a, ‖(c)‖ = c → R rimane uguale ad A, Q = I.
    Symbol sa("a"), sc("c");
    ctx.assumptions().assume_positive(sa);
    ctx.assumptions().assume_positive(sc);
    ExprPtr a = ctx.arena().make<Symbol>("a");
    ExprPtr b = ctx.arena().make<Symbol>("b");
    ExprPtr c = ctx.arena().make<Symbol>("c");

    MatrixExpr A(2, 2);
    A(0, 0) = a; A(0, 1) = b;
    A(1, 0) = lit(0); A(1, 1) = c;

    auto r = qr_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok()) << "QR failed on symbolic upper-triangular";

    auto prod = multiply(r.value().Q, r.value().R, ctx);
    ASSERT_TRUE(prod.is_ok());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(prod.value()(i, j), A(i, j)))
                << "Q·R != A at (" << i << "," << j << ") symbolic";
}

TEST_F(QRTest, SymbolicQR_DefaultSignConvention_2x2) {
    // M-FIX QR sign-default test: per x[0] non noto come negativo, alpha = +‖x‖.
    // Verifica che Q·R = A anche quando primo elemento è simbolico.
    Symbol sx("x"), sy("y");
    ctx.assumptions().assume_positive(sx);
    ctx.assumptions().assume_positive(sy);
    ExprPtr x = ctx.arena().make<Symbol>("x");
    ExprPtr y = ctx.arena().make<Symbol>("y");

    MatrixExpr A(2, 2);
    A(0, 0) = x; A(0, 1) = lit(1);
    A(1, 0) = y; A(1, 1) = lit(2);

    // HC-F8-QR-HOUSEHOLDER-BAILOUT / HC-F8-PENDING-20-RESIDUE (closed 2026-06-13):
    // The Householder pipeline now produces a certifiable Q·R reconstruction
    // because together() reduces (N, D) by polynomial GCD content (see
    // .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Together_Polynomial_GCD_Reduction.md).
    auto r = qr_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok()) << "QR failed on symbolic default-sign 2x2";
    auto prod = multiply(r.value().Q, r.value().R, ctx);
    ASSERT_TRUE(prod.is_ok()) << "Q·R multiplication failed downstream";
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(prod.value()(i, j), A(i, j)))
                << "Q·R != A at (" << i << "," << j << ") symbolic default-sign";
}

TEST_F(QRTest, HouseholderCertificator3x3_QtQ_Identity_And_QR_Eq_A) {
    // F4.1b certificator: Q^T·Q ≡ I_n  AND  Q·R ≡ A.
    auto A = from_rows({{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}});
    auto r = qr_decompose(A, ctx);
    ASSERT_TRUE(r.is_ok());

    // Q·R == A
    auto prod = multiply(r.value().Q, r.value().R, ctx);
    ASSERT_TRUE(prod.is_ok());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_TRUE(entries_equal(prod.value()(i, j), A(i, j)))
                << "Q·R mismatch at (" << i << "," << j << ")";

    // Q^T·Q == I_n
    auto Qt = transpose(r.value().Q);
    ASSERT_TRUE(Qt.is_ok());
    auto qtq = multiply(Qt.value(), r.value().Q, ctx);
    ASSERT_TRUE(qtq.is_ok());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            ExprPtr expected = lit(i == j ? 1 : 0);
            EXPECT_TRUE(entries_equal(qtq.value()(i, j), expected))
                << "Q^T·Q not identity at (" << i << "," << j << ")";
        }
}

}  // namespace
