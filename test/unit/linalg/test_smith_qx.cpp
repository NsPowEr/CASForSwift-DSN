// CAS-F4.2b — Smith Normal Form su Q[x] (PID polinomiale).

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::linalg;

namespace {

class SmithQxTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr lit(long long v) {
        return ctx.arena().make<IntegerLit>(BigInt(v));
    }
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto tokens = Lexer(s).tokenize();
        if (tokens.is_error()) return nullptr;
        Parser p(tokens.value(), ctx.arena());
        auto e = p.parse();
        if (e.is_error()) return nullptr;
        auto simp = ctx.simplify(e.value());
        return simp.is_ok() ? simp.value() : e.value();
    }
    [[nodiscard]] bool entries_equal(ExprPtr a, ExprPtr b) {
        auto eq = symbolic::mathematically_equal(a, b, ctx);
        return eq.is_ok() && eq.value();
    }
};

// Smith su Z (regressione: il path esistente deve continuare a funzionare).
TEST_F(SmithQxTest, Z_path_regression_2x2) {
    MatrixExpr A(2, 2, {lit(2), lit(4), lit(6), lit(10)});
    auto s = smith_normal_form(A, ctx);
    ASSERT_TRUE(s.is_ok());
    // Diagonal entries: invariant factors = (2, 2) (det = 4, gcd = 2 → d1=2, d2=2).
    EXPECT_TRUE(entries_equal(s.value().S(0, 0), lit(2)));
    EXPECT_TRUE(entries_equal(s.value().S(1, 1), lit(2)));
}

// Smith su Q[x]: matrice diagonale [[x, 0], [0, x^2]] → SNF [[x, 0], [0, x^2]].
TEST_F(SmithQxTest, Qx_diagonal_x_x2) {
    auto x = parse("x");
    auto x2 = parse("x^2");
    ASSERT_NE(x, nullptr);
    ASSERT_NE(x2, nullptr);
    MatrixExpr A(2, 2, {x, lit(0), lit(0), x2});
    auto s = smith_normal_form(A, ctx);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    EXPECT_TRUE(entries_equal(s.value().S(0, 0), x));
    EXPECT_TRUE(entries_equal(s.value().S(1, 1), x2));
}

// Smith su Q[x]: matrice [[x+1, 0], [0, x-1]] → invariant factors (1, x^2-1).
TEST_F(SmithQxTest, Qx_diagonal_xp1_xm1) {
    auto xp1 = parse("x + 1");
    auto xm1 = parse("x - 1");
    ASSERT_NE(xp1, nullptr);
    ASSERT_NE(xm1, nullptr);
    MatrixExpr A(2, 2, {xp1, lit(0), lit(0), xm1});
    auto s = smith_normal_form(A, ctx);
    ASSERT_TRUE(s.is_ok()) << s.error().message;
    // d1 = gcd(x+1, x-1) = 1 (Bezout: (x-1) - (x+1) = -2 → gcd = costante 2 → monico 1).
    // d2 = (x+1)(x-1)/d1 = x^2 - 1.
    // Verifica U·A·V == diag(d1, d2): cert structural.
    auto UA = multiply(s.value().U, A, ctx);
    ASSERT_TRUE(UA.is_ok());
    auto UAV = multiply(UA.value(), s.value().V, ctx);
    ASSERT_TRUE(UAV.is_ok());
    // UAV must equal S
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(UAV.value()(i, j), s.value().S(i, j)))
                << "U·A·V != S at (" << i << "," << j << ")";
}

// Smith su Q[x]: matrice NON diagonale x·I − A con A = [[0,−1],[1,0]] (rotazione),
// cioè [[x, 1], [−1, x]].  Esercita davvero il pivoting off-diagonal e la
// Bezout-elimination (il pivot iniziale è la costante 1, grado 0).
//   gcd di tutte le entrate = 1  →  d1 = 1,  d2 = det = x² + 1.
//   SNF = diag(1, x² + 1).
TEST_F(SmithQxTest, Qx_nondiagonal_companion_x2plus1) {
    auto x = parse("x");
    ASSERT_NE(x, nullptr);
    MatrixExpr A(2, 2, {x, lit(1), lit(-1), x});
    auto s = smith_normal_form(A, ctx);
    ASSERT_TRUE(s.is_ok()) << s.error().message;

    // Invariant factors (monic): d1 = 1, d2 = x² + 1.
    EXPECT_TRUE(entries_equal(s.value().S(0, 0), lit(1)))
        << "d1 expected 1";
    EXPECT_TRUE(entries_equal(s.value().S(1, 1), parse("x^2 + 1")))
        << "d2 expected x^2 + 1";
    EXPECT_TRUE(entries_equal(s.value().S(0, 1), lit(0)));
    EXPECT_TRUE(entries_equal(s.value().S(1, 0), lit(0)));

    // Transformation certificate: U·A·V == S.
    auto UA = multiply(s.value().U, A, ctx);
    ASSERT_TRUE(UA.is_ok());
    auto UAV = multiply(UA.value(), s.value().V, ctx);
    ASSERT_TRUE(UAV.is_ok());
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            EXPECT_TRUE(entries_equal(UAV.value()(i, j), s.value().S(i, j)))
                << "U·A·V != S at (" << i << "," << j << ")";
}

}  // namespace
