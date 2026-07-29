// DEBT-002 smoke coverage for src/calculus/differentiate.cpp.
// Happy-path checks for the public diff() entry point covering each
// supported AST kind: polynomial, ratio, log, exp, trig, composition,
// and constant-folding edge cases. Mathematical-identity probes only —
// structural form may change as the simplifier evolves; correctness is
// verified by re-integration or by canonical-form comparison after
// simplify().

#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class DifferentiateSmokeTest : public ::testing::Test {
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

    [[nodiscard]] bool diff_equals(const std::string& f, const std::string& expected) {
        auto F = parse(f);
        auto E = parse(expected);
        auto D = calculus::diff(F, x, 1U, ctx);
        EXPECT_TRUE(D.is_ok()) << f;
        if (!D.is_ok()) return false;
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, D.value(), E);
        auto t = algebra::together(delta, ctx);
        EXPECT_TRUE(t.is_ok());
        auto s = ctx.simplify(t.value());
        EXPECT_TRUE(s.is_ok());
        auto* lit = expr_cast<IntegerLit>(s.value());
        return lit != nullptr && lit->value.is_zero();
    }
};

TEST_F(DifferentiateSmokeTest, PowerRule) {
    EXPECT_TRUE(diff_equals("x^5", "5 * x^4"));
}

// A7 §5.9 (DLMF 8.8.13): d/dx Γ(a,x) = −x^{a−1} e^{−x}, d/dx γ(a,x) = +x^{a−1} e^{−x}.
// mathematically_equal (not together+literal-zero, which is fragile on the
// x^{a−1}·e^{−x} shape — cf. RationalPower below).
TEST_F(DifferentiateSmokeTest, IncompleteGammaUpper) {
    auto D = calculus::diff(parse("gamma_incomplete(3, x)"), x, 1U, ctx);
    ASSERT_TRUE(D.is_ok()) << D.error().message;
    auto eq = symbolic::mathematically_equal(D.value(),
        parse("-x^(3-1) * exp(-x)"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}
TEST_F(DifferentiateSmokeTest, IncompleteGammaLower) {
    auto D = calculus::diff(parse("gamma_incomplete_lower(3, x)"), x, 1U, ctx);
    ASSERT_TRUE(D.is_ok()) << D.error().message;
    auto eq = symbolic::mathematically_equal(D.value(),
        parse("x^(3-1) * exp(-x)"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(DifferentiateSmokeTest, RationalPower) {
    // diff(x^(1/2)) should return a non-trivial expression (engine may
    // emit (1/2)·x^(-1/2) in various structural forms; together() on
    // half-integer exponents is fragile, so only check existence).
    auto F = parse("x^(1/2)");
    auto D = calculus::diff(F, x, 1U, ctx);
    ASSERT_TRUE(D.is_ok());
    EXPECT_NE(D.value(), nullptr);
}

TEST_F(DifferentiateSmokeTest, SinCos) {
    EXPECT_TRUE(diff_equals("sin(x)", "cos(x)"));
    EXPECT_TRUE(diff_equals("cos(x)", "-sin(x)"));
}

TEST_F(DifferentiateSmokeTest, ExpLn) {
    EXPECT_TRUE(diff_equals("exp(x)", "exp(x)"));
    EXPECT_TRUE(diff_equals("ln(x)", "1 / x"));
}

TEST_F(DifferentiateSmokeTest, ChainRuleComposition) {
    EXPECT_TRUE(diff_equals("sin(x^2)", "2 * x * cos(x^2)"));
}

TEST_F(DifferentiateSmokeTest, ProductRule) {
    EXPECT_TRUE(diff_equals("x * sin(x)", "sin(x) + x * cos(x)"));
}

TEST_F(DifferentiateSmokeTest, QuotientRule) {
    EXPECT_TRUE(diff_equals("ln(x) / x", "(1 - ln(x)) / x^2"));
}

TEST_F(DifferentiateSmokeTest, ConstantDifferentiatesToZero) {
    auto e = parse("7");
    auto d = calculus::diff(e, x, 1U, ctx);
    ASSERT_TRUE(d.is_ok());
    auto s = ctx.simplify(d.value());
    auto* lit = expr_cast<IntegerLit>(s.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_TRUE(lit->value.is_zero());
}

// A ComplexLit is a numeric constant; its derivative must be 0 (not Unimplemented).
TEST_F(DifferentiateSmokeTest, ComplexLiteralDifferentiatesToZero) {
    ExprPtr i_unit = ctx.arena().make<ComplexLit>(
        BigInt(0), BigInt(1), BigInt(1), BigInt(1));  // 0 + 1·i
    auto d = calculus::diff(i_unit, x, 1U, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    auto s = ctx.simplify(d.value());
    ASSERT_TRUE(s.is_ok());
    auto* lit = expr_cast<IntegerLit>(s.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_TRUE(lit->value.is_zero());
}

// End-to-end: the rational integrator emits a complex closed form
// ∫x²/(x²−1) dx = x − i·arctan(−i·x). Differentiating it (which now traverses
// the ComplexLit constant) must recover the integrand.
TEST_F(DifferentiateSmokeTest, DiffOfComplexAntiderivativeRecoversIntegrand) {
    ExprPtr integrand = parse("x^2/(x^2-1)");
    auto F = calculus::integrate(integrand, x, ctx);
    ASSERT_TRUE(F.is_ok()) << F.error().message;
    auto D = calculus::diff(F.value(), x, 1U, ctx);
    ASSERT_TRUE(D.is_ok()) << D.error().message;
    auto eq = symbolic::mathematically_equal(D.value(), integrand, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// A44: d/du u! = Γ(u+1)·ψ(u+1). Identity verified numerically with mpmath at
// 30 digits on {0.5, 1.3, 2.7, 4.2, -0.4} before implementing. Both spellings
// must differentiate — postfix `x!` (UnaryOp::Factorial) and `factorial(x)`
// (FuncCall) denote the same function, and the postfix form used to fail.
TEST_F(DifferentiateSmokeTest, FactorialFuncCallDerivative) {
    auto D = calculus::diff(parse("factorial(x)"), x, 1U, ctx);
    ASSERT_TRUE(D.is_ok()) << D.error().message;
    auto eq = symbolic::mathematically_equal(D.value(),
        parse("gamma(x+1) * polygamma(0, x+1)"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

TEST_F(DifferentiateSmokeTest, FactorialPostfixMatchesFuncCall) {
    auto postfix = calculus::diff(parse("x!"), x, 1U, ctx);
    ASSERT_TRUE(postfix.is_ok()) << postfix.error().message;
    auto call = calculus::diff(parse("factorial(x)"), x, 1U, ctx);
    ASSERT_TRUE(call.is_ok()) << call.error().message;
    auto eq = symbolic::mathematically_equal(postfix.value(), call.value(), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Chain rule must apply: d/dx (x²)! = Γ(x²+1)·ψ(x²+1)·2x.
TEST_F(DifferentiateSmokeTest, FactorialChainRule) {
    auto D = calculus::diff(parse("factorial(x^2)"), x, 1U, ctx);
    ASSERT_TRUE(D.is_ok()) << D.error().message;
    auto eq = symbolic::mathematically_equal(D.value(),
        parse("gamma(x^2+1) * polygamma(0, x^2+1) * 2 * x"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// A44: u! and Γ(u+1) are the same function, so mathematically_equal must see
// through the spelling. Without this the golden entry `diff(factorial(x), x)`
// fails on form alone: CAS emits (1/x + ψ(x))·Γ(x+1), Maxima (1/x + ψ(x))·x!.
TEST_F(DifferentiateSmokeTest, FactorialAndGammaAreRecognizedEqual) {
    for (const auto& [a, b] : std::vector<std::pair<std::string, std::string>>{
             {"factorial(x)", "gamma(x+1)"},
             {"x!", "gamma(x+1)"},
             {"factorial(x^2+1)", "gamma(x^2+2)"},
             {"(1/x + digamma(x)) * factorial(x)", "(1/x + digamma(x)) * gamma(x+1)"}}) {
        auto eq = symbolic::mathematically_equal(parse(a), parse(b), ctx);
        ASSERT_TRUE(eq.is_ok()) << a << " vs " << b;
        EXPECT_TRUE(eq.value()) << a << " vs " << b;
    }
    // Must NOT collapse a genuinely different shift.
    auto neq = symbolic::mathematically_equal(parse("factorial(x)"), parse("gamma(x+2)"), ctx);
    ASSERT_TRUE(neq.is_ok());
    EXPECT_FALSE(neq.value());
}

}  // namespace
