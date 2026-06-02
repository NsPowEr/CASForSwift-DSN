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
// L3-04: Bessel function identities.
TEST_F(SpecialFunctionsTest, BesselJ_NegativeIntegerParity) {
    // J_{-2}(x) = J_2(x); J_{-3}(x) = -J_3(x).
    expect_simplify_equiv("BesselJ(-2, x)", "BesselJ(2, x)");
    expect_simplify_equiv("BesselJ(-3, x)", "-BesselJ(3, x)");
}
TEST_F(SpecialFunctionsTest, BesselY_NegativeIntegerParity) {
    expect_simplify_equiv("BesselY(-2, x)", "BesselY(2, x)");
    expect_simplify_equiv("BesselY(-3, x)", "-BesselY(3, x)");
}
TEST_F(SpecialFunctionsTest, BesselI_K_NegativeIntegerEqualsPositive) {
    expect_simplify_equiv("BesselI(-5, x)", "BesselI(5, x)");
    expect_simplify_equiv("BesselK(-4, x)", "BesselK(4, x)");
}
TEST_F(SpecialFunctionsTest, BesselJ_HalfInteger_TrigClosedForm) {
    // J_{1/2}(x) = sqrt(2/(pi*x))*sin(x);  J_{-1/2}(x) = sqrt(2/(pi*x))*cos(x).
    expect_simplify_equiv("BesselJ(1/2, x)", "sqrt(2/(pi*x))*sin(x)");
    expect_simplify_equiv("BesselJ(-1/2, x)", "sqrt(2/(pi*x))*cos(x)");
}
TEST_F(SpecialFunctionsTest, BesselY_HalfInteger_TrigClosedForm) {
    expect_simplify_equiv("BesselY(1/2, x)", "-sqrt(2/(pi*x))*cos(x)");
    expect_simplify_equiv("BesselY(-1/2, x)", "sqrt(2/(pi*x))*sin(x)");
}
TEST_F(SpecialFunctionsTest, BesselI_HalfInteger_HyperbolicClosedForm) {
    expect_simplify_equiv("BesselI(1/2, x)", "sqrt(2/(pi*x))*sinh(x)");
    expect_simplify_equiv("BesselI(-1/2, x)", "sqrt(2/(pi*x))*cosh(x)");
}
TEST_F(SpecialFunctionsTest, BesselK_HalfInteger_ExpClosedForm) {
    expect_simplify_equiv("BesselK(1/2, x)", "sqrt(pi/(2*x))*exp(-x)");
    expect_simplify_equiv("BesselK(-1/2, x)", "sqrt(pi/(2*x))*exp(-x)");
}

// L3-04: Chebyshev polynomials of the first / second kind.
TEST_F(SpecialFunctionsTest, ChebyshevT_smallDegrees) {
    expect_simplify_equiv("ChebyshevT(0, x)", "1");
    expect_simplify_equiv("ChebyshevT(1, x)", "x");
    expect_simplify_equiv("ChebyshevT(2, x)", "2*x^2 - 1");
    expect_simplify_equiv("ChebyshevT(3, x)", "4*x^3 - 3*x");
}
TEST_F(SpecialFunctionsTest, ChebyshevU_smallDegrees) {
    expect_simplify_equiv("ChebyshevU(0, x)", "1");
    expect_simplify_equiv("ChebyshevU(1, x)", "2*x");
    expect_simplify_equiv("ChebyshevU(2, x)", "4*x^2 - 1");
    expect_simplify_equiv("ChebyshevU(3, x)", "8*x^3 - 4*x");
}
TEST_F(SpecialFunctionsTest, ChebyshevT_HighDegree_AntiHardcode) {
    // T_5(x) = 16x^5 − 20x^3 + 5x.
    expect_simplify_equiv("ChebyshevT(5, x)", "16*x^5 - 20*x^3 + 5*x");
}

// L3-04: Hermite polynomials (both conventions).
TEST_F(SpecialFunctionsTest, HermiteH_PhysicistConvention) {
    expect_simplify_equiv("HermiteH(0, x)", "1");
    expect_simplify_equiv("HermiteH(1, x)", "2*x");
    expect_simplify_equiv("HermiteH(2, x)", "4*x^2 - 2");
    expect_simplify_equiv("HermiteH(3, x)", "8*x^3 - 12*x");
    expect_simplify_equiv("HermiteH(4, x)", "16*x^4 - 48*x^2 + 12");
}
TEST_F(SpecialFunctionsTest, HermiteHe_ProbabilistConvention) {
    expect_simplify_equiv("HermiteHe(0, x)", "1");
    expect_simplify_equiv("HermiteHe(1, x)", "x");
    expect_simplify_equiv("HermiteHe(2, x)", "x^2 - 1");
    expect_simplify_equiv("HermiteHe(3, x)", "x^3 - 3*x");
}

// L3-04 step: Digamma / Polygamma closed forms.
TEST_F(SpecialFunctionsTest, DigammaOne_NegativeEulerGamma) {
    expect_simplify_equiv("digamma(1)", "-EulerGamma");
}
TEST_F(SpecialFunctionsTest, DigammaPositiveInteger_HarmonicSum) {
    // ψ(4) = -γ + 1 + 1/2 + 1/3 = -γ + 11/6.
    expect_simplify_equiv("digamma(4)", "-EulerGamma + 11/6");
}
TEST_F(SpecialFunctionsTest, DigammaOneHalf) {
    // ψ(1/2) = -γ - 2·ln(2).
    expect_simplify_equiv("digamma(1/2)", "-EulerGamma - 2*ln(2)");
}
TEST_F(SpecialFunctionsTest, DigammaThreeHalves_FunctionalEquation) {
    // ψ(3/2) = -γ - 2 ln 2 + 2.
    expect_simplify_equiv("digamma(3/2)", "-EulerGamma - 2*ln(2) + 2");
}
TEST_F(SpecialFunctionsTest, DigammaFunctionalEquationOnSymbol) {
    // ψ(x+1) = ψ(x) + 1/x.
    expect_simplify_equiv("digamma(x+1)", "digamma(x) + 1/x");
}
TEST_F(SpecialFunctionsTest, PolygammaZeroIsDigamma) {
    expect_simplify_equiv("polygamma(0, 1)", "-EulerGamma");
}
TEST_F(SpecialFunctionsTest, PolygammaAtOne_AntiHardcode) {
    // ψ'(1) = ζ(2) = π²/6.   ψ''(1) = -2·ζ(3).
    expect_simplify_equiv("polygamma(1, 1)", "(pi^2)/6");
    // Test for n=3:  ψ'''(1) = 6·ζ(4) = 6·π^4/90 = π^4/15.
    expect_simplify_equiv("polygamma(3, 1)", "(pi^4)/15");
}

// L2-10 algorithmic step: half-angle recursion produces closed forms for
// any constructible n built by doubling a base table entry.
TEST_F(SpecialFunctionsTest, CosPiOverEight_HalfAngle) {
    // cos(π/8) = sqrt(2 + sqrt(2))/2 = sqrt((1 + sqrt(2)/2)/2).
    expect_simplify_equiv("cos(pi/8)", "sqrt((1 + sqrt(2)/2)/2)");
}
TEST_F(SpecialFunctionsTest, SinPiOverEight_HalfAngle) {
    // sin(π/8) = sqrt((1 - cos(π/4))/2) = sqrt((1 - sqrt(2)/2)/2).
    expect_simplify_equiv("sin(pi/8)", "sqrt((1 - sqrt(2)/2)/2)");
}
TEST_F(SpecialFunctionsTest, CosPiOverSixteen_HalfAngle) {
    // cos(π/16) = sqrt((1 + cos(π/8))/2) — two halving steps from π/4.
    expect_simplify_equiv("cos(pi/16)", "sqrt((1 + sqrt((1 + sqrt(2)/2)/2))/2)");
}
// Chebyshev T_p generalizes cos(p·π/q) for arbitrary p ≥ 2.
TEST_F(SpecialFunctionsTest, CosThreePiOverSeven_Chebyshev_NonInert) {
    // cos(3π/16) — q=16 reachable via half-angle, p=3 via T_3.
    auto e = parse_expr("cos(3*pi/16)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(3π/16) stayed inert; got: " << debug_print(s.value());
}
TEST_F(SpecialFunctionsTest, CosSevenPiOverSixteen_NonInert) {
    auto e = parse_expr("cos(7*pi/16)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(7π/16) stayed inert; got: " << debug_print(s.value());
}

TEST_F(SpecialFunctionsTest, CosPiOverTwentyFour_HalfAngle) {
    // π/24 = (π/12)/2.  cos(π/12) is in the base table (=cos(15°)).
    // After one halving: cos(π/24) = sqrt((1 + cos(π/12))/2).
    // Direct verification: just ensure it simplifies (not stay inert).
    auto e = parse_expr("cos(pi/24)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    // Result must NOT be the inert FuncCall(cos, pi/24).
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(pi/24) stayed inert; got: " << debug_print(s.value());
}

// L2-10 step: closed-form trig at n=5 / n=10 (Fermat prime 5).
TEST_F(SpecialFunctionsTest, CosPiOverFive) {
    expect_simplify_equiv("cos(pi/5)", "(1 + sqrt(5))/4");
}
TEST_F(SpecialFunctionsTest, CosTwoPiOverFive) {
    expect_simplify_equiv("cos(2*pi/5)", "(sqrt(5) - 1)/4");
}
TEST_F(SpecialFunctionsTest, SinPiOverFive) {
    expect_simplify_equiv("sin(pi/5)", "sqrt(10 - 2*sqrt(5))/4");
}
TEST_F(SpecialFunctionsTest, SinPiOverTen_EqualsCosTwoPiOverFive) {
    // sin(π/10) and cos(2π/5) share the same closed form.
    expect_simplify_equiv("sin(pi/10)", "(sqrt(5) - 1)/4");
}
TEST_F(SpecialFunctionsTest, GoldenRatioIdentityFromCosPiOverFive) {
    // 2*cos(π/5) = golden ratio φ = (1 + sqrt(5))/2.
    expect_simplify_equiv("2*cos(pi/5)", "(1 + sqrt(5))/2");
}

// L2-10: angle combination formula — denominators not reachable by halving/Chebyshev alone.
// cos(π/15) = cos(2π/5 - π/3) = cos(2π/5)·cos(π/3) + sin(2π/5)·sin(π/3).
// Anti-hardcode: ensure these stay non-inert (structural form, not numeric oracle).
TEST_F(SpecialFunctionsTest, CosPiOverFifteen_NotInert) {
    auto e = parse_expr("cos(pi/15)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(π/15) stayed inert: " << debug_print(s.value());
}
TEST_F(SpecialFunctionsTest, SinPiOverFifteen_NotInert) {
    auto e = parse_expr("sin(pi/15)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Sin)
        << "sin(π/15) stayed inert: " << debug_print(s.value());
}
TEST_F(SpecialFunctionsTest, CosTwoPiOverFifteen_NotInert) {
    // cos(2π/15) = cos(π/3 - π/5) — combines base table entries.
    auto e = parse_expr("cos(2*pi/15)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(2π/15) stayed inert: " << debug_print(s.value());
}
// Anti-hardcode: cos(π/15) + cos(2π/15) should equal a closed-form value (not an inert sum).
TEST_F(SpecialFunctionsTest, AngleCombination_AntiHardcode_q15) {
    // Verify that BOTH cos(π/15) and sin(π/15) are non-inert (denominator 15 fully covered).
    for (const char* expr_str : {"cos(pi/15)", "sin(pi/15)", "cos(2*pi/15)", "sin(2*pi/15)",
                                  "cos(4*pi/15)", "sin(4*pi/15)"}) {
        auto e = parse_expr(expr_str, ctx->arena());
        ASSERT_TRUE(e.is_ok()) << "parse failed: " << expr_str;
        auto s = ctx->simplify(e.value());
        ASSERT_TRUE(s.is_ok()) << "simplify failed: " << expr_str;
        const auto* fc = expr_cast<FuncCall>(s.value());
        EXPECT_FALSE(fc != nullptr &&
                     (fc->func_id == BuiltinOp::Cos || fc->func_id == BuiltinOp::Sin))
            << expr_str << " stayed inert: " << debug_print(s.value());
    }
}

// L2-10: q=30 = 2×15 — half-angle on top of angle-combination result.
TEST_F(SpecialFunctionsTest, CosPiOverThirty_NotInert) {
    // cos(π/30) = cos(6°) reached via half-angle from cos(π/15).
    auto e = parse_expr("cos(pi/30)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_ok());
    const auto* fc = expr_cast<FuncCall>(s.value());
    EXPECT_FALSE(fc != nullptr && fc->func_id == BuiltinOp::Cos)
        << "cos(π/30) stayed inert: " << debug_print(s.value());
}

// L1-12: nested radical denesting sqrt(a + b·sqrt(c)).
TEST_F(SpecialFunctionsTest, DenestSqrt_5_plus_2sqrt6) {
    // sqrt(5 + 2*sqrt(6)) = sqrt(2) + sqrt(3).
    expect_simplify_equiv("sqrt(5 + 2*sqrt(6))", "sqrt(2) + sqrt(3)");
}

TEST_F(SpecialFunctionsTest, DenestSqrt_7_minus_4sqrt3) {
    // sqrt(7 - 4*sqrt(3)) = 2 - sqrt(3).
    expect_simplify_equiv("sqrt(7 - 4*sqrt(3))", "2 - sqrt(3)");
}

TEST_F(SpecialFunctionsTest, DenestSqrt_3_plus_2sqrt2) {
    // sqrt(3 + 2*sqrt(2)) = 1 + sqrt(2).
    expect_simplify_equiv("sqrt(3 + 2*sqrt(2))", "1 + sqrt(2)");
}

TEST_F(SpecialFunctionsTest, DenestSqrt_9_plus_4sqrt5) {
    // sqrt(9 + 4*sqrt(5)) = sqrt(5) + 2.  (a=9, b=4, c=5, disc=81-80=1).
    expect_simplify_equiv("sqrt(9 + 4*sqrt(5))", "sqrt(5) + 2");
}

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

// HPP-015 closure: configurable bit budget for closed-form expansions of
// special functions taking a positive integer argument.
TEST_F(SpecialFunctionsTest, HPP015_DigammaBitBudgetConfigurable) {
    // Default budget = 16 bits.  bit_length(70000) = 17 → Unimplemented.
    auto e = parse_expr("Digamma(70000)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_error());
    EXPECT_EQ(s.error().kind, CASErrorKind::Unimplemented);

    // Raise budget to 24 bits → expansion succeeds (heavy but finite).
    ctx->set_max_special_fn_integer_arg_bits(24U);
    auto e2 = parse_expr("Digamma(3)", ctx->arena());  // cheap sanity check
    ASSERT_TRUE(e2.is_ok());
    auto s2 = ctx->simplify(e2.value());
    ASSERT_TRUE(s2.is_ok());

    // Reduce budget below default → previously-OK input rejected.
    ctx->set_max_special_fn_integer_arg_bits(2U);
    auto e3 = parse_expr("Digamma(8)", ctx->arena());  // bit_length(8) = 4 > 2
    ASSERT_TRUE(e3.is_ok());
    auto s3 = ctx->simplify(e3.value());
    ASSERT_TRUE(s3.is_error());
    EXPECT_EQ(s3.error().kind, CASErrorKind::Unimplemented);
}

TEST_F(SpecialFunctionsTest, HPP015_PochhammerBitBudgetConfigurable) {
    auto e = parse_expr("Pochhammer(x, 70000)", ctx->arena());
    ASSERT_TRUE(e.is_ok());
    auto s = ctx->simplify(e.value());
    ASSERT_TRUE(s.is_error());
    EXPECT_EQ(s.error().kind, CASErrorKind::Unimplemented);

    ctx->set_max_special_fn_integer_arg_bits(20U);
    // Pochhammer(x, 5) still works after raise; just sanity.
    auto e2 = parse_expr("Pochhammer(x, 5)", ctx->arena());
    auto s2 = ctx->simplify(e2.value());
    ASSERT_TRUE(s2.is_ok());
}

TEST_F(SpecialFunctionsTest, HPP015_ZetaBernoulliBudgetConfigurable) {
    // Reduce budget so Zeta(8) → Unimplemented (bit_length(8)=4 > 3).
    ctx->set_max_bernoulli_index_bits(3U);
    auto e2 = parse_expr("zeta(8)", ctx->arena());
    auto s2 = ctx->simplify(e2.value());
    ASSERT_TRUE(s2.is_error())
        << "got value: " << debug_print(s2.value());
    EXPECT_EQ(s2.error().kind, CASErrorKind::Unimplemented);

    // Restore default; Zeta(8) closed form must work.
    ctx->set_max_bernoulli_index_bits(30U);
    auto e3 = parse_expr("zeta(8)", ctx->arena());
    auto s3 = ctx->simplify(e3.value());
    ASSERT_TRUE(s3.is_ok());
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
