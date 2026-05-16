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

// HC-003 anti-hardcode: zeta(2k) must be computed via Bernoulli formula for ANY k,
// not via a fixed lookup table {2,4,6,8,10,12}.  Previously zeta(14) silently
// stayed inert because nu=14 was outside the switch.  After HC-003 fix, the
// closed-form is computed dynamically.
TEST_F(SpecialFunctionsTest, ZetaFourteenViaBernoulli_AntiHardcode) {
    // B_14 = 7/6.  zeta(14) = (-1)^8 · 2^13 · pi^14 · (7/6) / 14!
    //                       = 8192 · 7 / (6 · 87178291200) · pi^14
    //                       = 57344 / 523069747200 · pi^14
    //                       = 2 · pi^14 / 18243225  (after reduction)
    expect_simplify_equiv("zeta(14)", "(2*(pi^14))/18243225");
}

TEST_F(SpecialFunctionsTest, ZetaSixteenViaBernoulli_AntiHardcode) {
    // B_16 = -3617/510.  zeta(16) = 3617 · pi^16 / 325641566250.
    expect_simplify_equiv("zeta(16)", "(3617*(pi^16))/325641566250");
}

TEST_F(SpecialFunctionsTest, ZetaNegativeNineViaBernoulli_AntiHardcode) {
    // zeta(-9) = -B_10 / 10 = -(5/66)/10 = -1/132.
    expect_simplify_equiv("zeta(-9)", "-1/132");
}

TEST_F(SpecialFunctionsTest, ZetaNegativeElevenViaBernoulli_AntiHardcode) {
    // zeta(-11) = -B_12 / 12 = -(-691/2730)/12 = 691/32760.
    expect_simplify_equiv("zeta(-11)", "691/32760");
}

// L3-04: Gamma reflection identity Γ(z)·Γ(1−z) = π / sin(πz).
TEST_F(SpecialFunctionsTest, GammaReflectionSymbolic) {
    // Γ(x) · Γ(1 − x) should reduce to π / sin(π·x).
    expect_simplify_equiv("gamma(x)*gamma(1-x)", "pi/sin(pi*x)");
}

TEST_F(SpecialFunctionsTest, GammaReflectionWithSumPlusInteger) {
    // Γ(x + 2) · Γ(-1 - x).  Sum of args = x+2 + (-1-x) = 1.  Functional
    // equation strips the integer shifts; reflection still applies because
    // Γ(x)·Γ(-x) = -π/(x·sin(πx)) and the (x+1) cancels.  Expected end
    // result is π / sin(πx) (sin is 2π-periodic so this equals
    // π / sin(π(x+2))).
    expect_simplify_equiv("gamma(x+2)*gamma(-1-x)", "pi/sin(pi*x)");
}

// L3-04: Legendre P_n via Bonnet recurrence (no hardcoded table).
TEST_F(SpecialFunctionsTest, LegendreP_0_1_2) {
    expect_simplify_equiv("LegendreP(0, x)", "1");
    expect_simplify_equiv("LegendreP(1, x)", "x");
    // P_2(x) = (3x^2 − 1)/2.
    expect_simplify_equiv("LegendreP(2, x)", "(3*x^2 - 1)/2");
}

TEST_F(SpecialFunctionsTest, LegendreP_3_4) {
    // P_3(x) = (5x^3 − 3x)/2.
    expect_simplify_equiv("LegendreP(3, x)", "(5*x^3 - 3*x)/2");
    // P_4(x) = (35x^4 − 30x^2 + 3)/8.
    expect_simplify_equiv("LegendreP(4, x)", "(35*x^4 - 30*x^2 + 3)/8");
}

TEST_F(SpecialFunctionsTest, LegendreP_AntiHardcode_HighDegree) {
    // P_6(x) = (231x^6 − 315x^4 + 105x^2 − 5)/16.  Anti-hardcode: degree
    // generated by recurrence, not table.
    expect_simplify_equiv("LegendreP(6, x)",
        "(231*x^6 - 315*x^4 + 105*x^2 - 5)/16");
}

TEST_F(SpecialFunctionsTest, LegendreP_AtOne_Equals_One) {
    // Identity: P_n(1) = 1 for any non-negative n.  Verify for n=5.
    expect_simplify_equiv("LegendreP(5, 1)", "1");
}

// HC-004: fresh symbol generator never collides with itself nor with
// user-defined names.
// L3-04: Beta function via Gamma identity.
TEST_F(SpecialFunctionsTest, BetaIntegerArgsReducesToRational) {
    // B(2, 3) = Γ(2)Γ(3)/Γ(5) = 1·2/24 = 1/12.
    expect_simplify_equiv("Beta(2, 3)", "1/12");
}

TEST_F(SpecialFunctionsTest, BetaHalfIntegerProducesPi) {
    // B(1/2, 1/2) = Γ(1/2)^2 / Γ(1) = π.
    expect_simplify_equiv("Beta(1/2, 1/2)", "pi");
}

TEST_F(SpecialFunctionsTest, BetaSymmetryUnderSwap) {
    // B(x, 1-x) ↔ Γ(x)Γ(1-x)/Γ(1) = π/sin(πx) (Gamma reflection composition).
    expect_simplify_equiv("Beta(x, 1-x)", "pi/sin(pi*x)");
}

// L3-04: Pochhammer (rising factorial).
TEST_F(SpecialFunctionsTest, PochhammerIntegerNonNegative) {
    // (x)_0 = 1
    expect_simplify_equiv("Pochhammer(x, 0)", "1");
    // (x)_3 = x*(x+1)*(x+2)
    expect_simplify_equiv("Pochhammer(x, 3)", "x*(x+1)*(x+2)");
    // (1)_n = n!  → (1)_5 = 120
    expect_simplify_equiv("Pochhammer(1, 5)", "120");
}

TEST_F(SpecialFunctionsTest, PochhammerNegativeIndex) {
    // (x)_(-2) = 1 / ((x-1)*(x-2))
    expect_simplify_equiv("Pochhammer(x, -2)", "1/((x-1)*(x-2))");
}

TEST_F(SpecialFunctionsTest, MakeFreshSymbol_UniqueAndAvoidsUserScope) {
    // Pre-populate context with user-defined variables that look like
    // fresh-symbol candidates.
    ctx->define(Symbol("C_1"), ctx->arena().make<IntegerLit>(BigInt(0)));
    ctx->define(Symbol("C_2"), ctx->arena().make<IntegerLit>(BigInt(0)));
    Symbol a = ctx->make_fresh_symbol("C");
    Symbol b = ctx->make_fresh_symbol("C");
    Symbol c = ctx->make_fresh_symbol("C");
    EXPECT_NE(a.name, b.name);
    EXPECT_NE(b.name, c.name);
    EXPECT_NE(a.name, c.name);
    EXPECT_NE(a.name, std::string("C_1"));
    EXPECT_NE(a.name, std::string("C_2"));
    EXPECT_NE(b.name, std::string("C_1"));
    EXPECT_NE(b.name, std::string("C_2"));
    EXPECT_NE(c.name, std::string("C_1"));
    EXPECT_NE(c.name, std::string("C_2"));
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
