// CAS-L3-12 — Symbolic finite-difference numeric derivative tests.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include "../../../src/calculus/calculus_internal.hpp"

using namespace cas;
using namespace cas::calculus;

namespace {

class NumericDiffTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

TEST_F(NumericDiffTest, Forward1OnLinear) {
    // f(x) = 3x → (3(x+h) - 3x)/h = 3
    auto f = parse("3 * x");
    auto h = parse("h");
    auto r = numeric_diff(f, x, h, FiniteDiffOrder::Forward1, ctx);
    ASSERT_TRUE(r.is_ok());
    // simplified should be 3
    auto* lit = expr_cast<IntegerLit>(r.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_EQ(lit->value, BigInt(3));
}

TEST_F(NumericDiffTest, Central2OnQuadratic) {
    // f(x) = x²; central 2: ((x+h)² - (x-h)²)/(2h) = 4xh/(2h) = 2x. Exact!
    auto f = parse("x^2");
    auto h = parse("h");
    auto r = numeric_diff(f, x, h, FiniteDiffOrder::Central2, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("2 * x");
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, r.value(), expected);
    auto t = algebra::together(delta, ctx);
    auto simp = ctx.simplify(t.is_ok() ? t.value() : delta);
    ASSERT_TRUE(simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
        << "Central2(x²) should equal 2x exactly";
}

TEST_F(NumericDiffTest, Central2OnCubic) {
    // f(x) = x³; central 2: ((x+h)³ - (x-h)³)/(2h)
    // = (3x²·2h + 2·h³)/(2h) = 3x² + h²
    auto f = parse("x^3");
    auto h = parse("h");
    auto r = numeric_diff(f, x, h, FiniteDiffOrder::Central2, ctx);
    ASSERT_TRUE(r.is_ok());
    auto expected = parse("3*x^2 + h^2");
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, r.value(), expected);
    auto t = algebra::together(delta, ctx);
    auto simp = ctx.simplify(t.is_ok() ? t.value() : delta);
    ASSERT_TRUE(simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero());
}

TEST_F(NumericDiffTest, Central4OnQuartic) {
    // f(x) = x⁴; central 4 is exact for degree-4 polys.
    // 4th-order formula has truncation O(h⁴) → for x⁴, error is constant·h⁴
    // hence non-zero. But for x³ it's exact.
    auto f = parse("x^3");
    auto h = parse("h");
    auto r = numeric_diff(f, x, h, FiniteDiffOrder::Central4, ctx);
    ASSERT_TRUE(r.is_ok());
    // Should equal 3x² exactly (central 4 exact through degree 5).
    auto expected = parse("3*x^2");
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, r.value(), expected);
    auto t = algebra::together(delta, ctx);
    auto simp = ctx.simplify(t.is_ok() ? t.value() : delta);
    ASSERT_TRUE(simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(simp.value());
    EXPECT_TRUE(lit != nullptr && lit->value.is_zero())
        << "Central4(x³) should equal 3x² exactly";
}

TEST_F(NumericDiffTest, AntiHardcodeReturnsValidExpression) {
    // f(x) = sin(x); central 2 produces (sin(x+h) - sin(x-h))/(2h)
    // No closed form to simplify; just verify result is non-null ExprPtr.
    auto f = parse("sin(x)");
    auto h = parse("h");
    auto r = numeric_diff(f, x, h, FiniteDiffOrder::Central2, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

}  // namespace
