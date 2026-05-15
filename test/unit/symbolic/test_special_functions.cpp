// Tests for L3-04 Special Functions identities:
//   - Gamma(1/2) = sqrt(pi).
//   - Gamma(n + 1/2) via descending recursion + sqrt(pi).
//   - Gamma(z + n) -> z*(z+1)*...*(z+n-1) * Gamma(z) functional equation.
//   - erf odd:  erf(-x) = -erf(x).
//   - Derivative formulas continue to apply (regression check).

#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace cas::test {

namespace {
Result<ExprPtr> parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    if (t.is_error()) return fail<ExprPtr>(t.error());
    Parser p(t.value(), arena);
    return p.parse();
}
}  // namespace

class SpecialFunctionsTest : public ::testing::Test {
protected:
    void SetUp() override { ctx = std::make_unique<symbolic::CASContext>(); }

    void expect_simplify_equiv(const std::string& src, const std::string& expected_src) {
        auto e = parse_expr(src, ctx->arena());
        ASSERT_TRUE(e.is_ok()) << src;
        auto exp = parse_expr(expected_src, ctx->arena());
        ASSERT_TRUE(exp.is_ok()) << expected_src;
        auto sa = ctx->simplify(e.value());
        ASSERT_TRUE(sa.is_ok()) << sa.error().message;
        auto sb = ctx->simplify(exp.value());
        ASSERT_TRUE(sb.is_ok());
        auto eq = symbolic::mathematically_equal(sa.value(), sb.value(), *ctx);
        ASSERT_TRUE(eq.is_ok()) << eq.error().message;
        EXPECT_TRUE(eq.value())
            << "src=" << src << "\n"
            << "  got:      " << debug_print(sa.value()) << "\n"
            << "  expected: " << debug_print(sb.value());
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(SpecialFunctionsTest, GammaHalfEqualsSqrtPi) {
    expect_simplify_equiv("gamma(1/2)", "sqrt(pi)");
}

TEST_F(SpecialFunctionsTest, GammaThreeHalvesEqualsHalfSqrtPi) {
    // Gamma(3/2) = (1/2) * Gamma(1/2) = sqrt(pi)/2
    expect_simplify_equiv("gamma(3/2)", "sqrt(pi)/2");
}

TEST_F(SpecialFunctionsTest, GammaFiveHalvesEqualsThreeOverFourSqrtPi) {
    // Gamma(5/2) = (3/2) * Gamma(3/2) = (3/2) * (1/2) * sqrt(pi) = 3*sqrt(pi)/4
    expect_simplify_equiv("gamma(5/2)", "3*sqrt(pi)/4");
}

TEST_F(SpecialFunctionsTest, GammaMinusHalfEqualsMinusTwoSqrtPi) {
    // Gamma(-1/2) = Gamma(1/2)/(-1/2) = -2*sqrt(pi)
    expect_simplify_equiv("gamma(-1/2)", "-2*sqrt(pi)");
}

TEST_F(SpecialFunctionsTest, GammaPositiveIntegerStillFactorial) {
    // Regression: integer arg shortcut.
    expect_simplify_equiv("gamma(5)", "24");
    expect_simplify_equiv("gamma(1)", "1");
}

TEST_F(SpecialFunctionsTest, GammaShiftFunctionalEquationPlus1) {
    // Gamma(z + 1) = z * Gamma(z)
    expect_simplify_equiv("gamma(x + 1)", "x*gamma(x)");
}

TEST_F(SpecialFunctionsTest, GammaShiftFunctionalEquationPlus2) {
    // Gamma(z + 2) = z*(z+1)*Gamma(z)
    expect_simplify_equiv("gamma(x + 2)", "x*(x+1)*gamma(x)");
}

TEST_F(SpecialFunctionsTest, ErfIsOddOnNegSymbol) {
    // erf(-x) -> -erf(x)
    expect_simplify_equiv("erf(-x)", "-erf(x)");
}

TEST_F(SpecialFunctionsTest, ErfAtZeroIsZero) {
    expect_simplify_equiv("erf(0)", "0");
}

// Zeta special values:
TEST_F(SpecialFunctionsTest, ZetaZeroEqualsMinusHalf) {
    expect_simplify_equiv("zeta(0)", "-1/2");
}
TEST_F(SpecialFunctionsTest, ZetaMinusOneEqualsMinusOneTwelfth) {
    expect_simplify_equiv("zeta(-1)", "-1/12");
}
TEST_F(SpecialFunctionsTest, ZetaTwoEqualsBaselValue) {
    expect_simplify_equiv("zeta(2)", "(pi^2)/6");
}
TEST_F(SpecialFunctionsTest, ZetaFourEqualsPi4Over90) {
    expect_simplify_equiv("zeta(4)", "(pi^4)/90");
}
TEST_F(SpecialFunctionsTest, ZetaSixEqualsPi6Over945) {
    expect_simplify_equiv("zeta(6)", "(pi^6)/945");
}
TEST_F(SpecialFunctionsTest, ZetaNegativeEvenIsZero) {
    expect_simplify_equiv("zeta(-2)", "0");
    expect_simplify_equiv("zeta(-4)", "0");
}

TEST_F(SpecialFunctionsTest, DerivativeOfErfPreserved) {
    // d/dx erf(x) = 2/sqrt(pi) * exp(-x^2).  Regression check: ensure the
    // derivative survives (mathematically_equal across Pow(-1)/Div forms is
    // not guaranteed; we compare via the difference normalised to zero).
    auto e = parse_expr("erf(x)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto d = calculus::diff(e.value(), Symbol("x"), 1U, *ctx);
    ASSERT_TRUE(d.is_ok());
    // Note: parser binds Pow tighter than unary Neg, so "exp(-x^2)" would
    // parse as exp((-x)^2); the simplifier doesn't currently fold (-x)^2.
    // Use explicit parens to express -(x^2).
    auto expected = parse_expr("2/sqrt(pi) * exp(-(x^2))", ctx->arena());
    ASSERT_TRUE(expected.is_ok());
    ExprPtr diff_expr = ctx->arena().make<Binary>(BinaryOp::Sub, d.value(), expected.value());
    auto simp = ctx->simplify(diff_expr);
    ASSERT_TRUE(simp.is_ok());
    // Either it simplifies to zero structurally, or it simplifies to an
    // expression mathematically equivalent to zero — both acceptable.
    const auto* il = expr_cast<IntegerLit>(simp.value());
    if (il && il->value.is_zero()) {
        SUCCEED();
        return;
    }
    auto zero = ctx->arena().make<IntegerLit>(BigInt(0));
    auto eq = symbolic::mathematically_equal(simp.value(), zero, *ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value())
        << "Difference did not simplify to zero. Got: " << debug_print(simp.value());
}

}  // namespace cas::test
