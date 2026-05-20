// CAS-L3-04 — Jacobi polynomial P_n^{(α,β)}(x) via Bonnet recurrence.
//
// JacobiP(n, α, β, x):
//   P_0 = 1
//   P_1 = (α-β)/2 + (α+β+2)·x/2
//   higher via recurrence.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class JacobiTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
    [[nodiscard]] bool equiv(ExprPtr a, ExprPtr b) {
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, a, b);
        auto tog = algebra::together(delta, ctx);
        auto simp = ctx.simplify(tog.is_ok() ? tog.value() : delta);
        if (simp.is_error()) return false;
        if (auto* il = expr_cast<IntegerLit>(simp.value())) return il->value.is_zero();
        return false;
    }
};

TEST_F(JacobiTest, JacobiZeroEqualsOne) {
    auto e = parse("JacobiP(0, 1, 2, x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto* il = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(il, nullptr);
    EXPECT_EQ(il->value, BigInt(1));
}

TEST_F(JacobiTest, JacobiOneAlphaBetaPattern) {
    // P_1^{(α,β)}(x) = (α-β)/2 + (α+β+2)·x/2
    // For α=1, β=2: P_1(x) = (1-2)/2 + (1+2+2)·x/2 = -1/2 + 5x/2
    auto e = parse("JacobiP(1, 1, 2, x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("-1/2 + (5/2)*x");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(JacobiTest, JacobiOneAlphaEqualsBetaReducesToLegendreP1) {
    // When α=β=0: P_n^{(0,0)}(x) = P_n(x) (Legendre).
    // P_1^{(0,0)}(x) = (0-0)/2 + (0+0+2)·x/2 = x = LegendreP(1, x)
    auto e = parse("JacobiP(1, 0, 0, x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("x");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(JacobiTest, JacobiTwoSpecific) {
    // P_2^{(0,0)}(x) = LegendreP(2, x) = (3x²-1)/2
    auto e = parse("JacobiP(2, 0, 0, x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("(3*x^2 - 1)/2");
    EXPECT_TRUE(equiv(r.value(), expected));
}

TEST_F(JacobiTest, AntiHardcodeHigherDegreeReturnsExpression) {
    // P_5^{(1,1)}(x) — non-trivial expression, must NOT crash.
    auto e = parse("JacobiP(5, 1, 1, x)");
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

}  // namespace
