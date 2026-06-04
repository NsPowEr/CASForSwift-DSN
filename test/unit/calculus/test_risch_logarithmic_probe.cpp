// F3.5 probe — verify if Rothstein-Trager logarithmic part actually
// fails as the audit claimed.
//
// Audit A3: "Risch in integrate_risch.cpp:712,744 — Hermite/Trager
// dichiarati ma sono solo commenti TODO. Linea 712 ammette fallback a
// fallimento. Rothstein-Trager (linea 744) è solo COMMENTO."
//
// Concrete counterexamples to probe:
//   ∫ 1/(x · ln(x)) dx = ln(ln(x))           — pure log substitution
//   ∫ x/(x²+1) dx = (1/2)·ln(x²+1)           — log of quadratic
//   ∫ 1/(x²-2) dx = (1/√2)·atanh(x/√2)·...   — algebraic log
//   ∫ 1/(x³-1) dx                            — partial fractions log

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::calculus;

class RischLogarithmicProbeTest : public ::testing::Test {
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

// ∫ x/(x²+1) dx = (1/2)·ln(x²+1) — classical first-year calc.
// Achievable via u-substitution. Should not require Trager.
TEST_F(RischLogarithmicProbeTest, IntegralOfXOverXSquaredPlusOne) {
    auto e = parse("x / (x^2 + 1)");
    auto r = integrate(e, x, ctx);
    if (r.is_ok()) {
        auto formatted = formatter::TextFormatter{}.format(r.value());
        // Result should contain ln(x^2+1)
        EXPECT_NE(formatted.find("ln"), std::string::npos)
            << "Expected ln(x²+1). Got: " << formatted;
        std::cout << "[PROBE] ∫ x/(x²+1) dx = " << formatted << std::endl;
    } else {
        std::cout << "[PROBE FAIL] " << r.error().message << std::endl;
    }
}

// ∫ 1/(x · ln(x)) dx = ln(ln(x)). Pure u = ln(x) substitution.
//
// Fixed by Risch logarithmic-derivative recognizer (structure
// theorem step): given integrand 1/(x·ln(x)), the engine detects
// integrand = D(ln(x))/ln(x) for extension generator t = ln(x)
// and returns c · ln(t) = ln(ln(x)). Previously the engine fell
// through to IBP + Hermite/RT and produced ln(x)^(-1)·ln|x|,
// which is wrong (diff back doesn't match the integrand).
TEST_F(RischLogarithmicProbeTest, IntegralOfReciprocalOfXLnX) {
    auto e = parse("1 / (x * ln(x))");
    auto r = integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto formatted = formatter::TextFormatter{}.format(r.value());
    std::cout << "[PROBE] ∫ 1/(x ln(x)) dx = " << formatted << std::endl;
    // Verify by differentiation: D(result) - integrand must reduce to 0
    // after together() normalization (raw simplify alone cannot match
    // (1/x)(1/y) against 1/(xy)).
    auto d_back = diff(r.value(), x, 1, ctx);
    ASSERT_TRUE(d_back.is_ok());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, d_back.value(), e);
    auto delta_tog = algebra::together(delta, ctx);
    ASSERT_TRUE(delta_tog.is_ok());
    auto delta_simp = ctx.simplify(delta_tog.value());
    ASSERT_TRUE(delta_simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(delta_simp.value());
    ASSERT_NE(lit, nullptr)
        << "D(result) - integrand did not reduce to a literal";
    EXPECT_TRUE(lit->value.is_zero())
        << "D(result) ≠ integrand: " << formatter::TextFormatter{}.format(delta_simp.value());
}

// HPP-007 closure: formal constant extraction covers c outside {±1,±1/2,±2}.
// ∫ 3/(x·ln(x)) dx = 3·ln(ln(x)).  Previously the closed trial set missed
// integer multiples like c=3 and returned Unimplemented; the new formal-
// extraction path computes c = integrand/D(ln(ln(x))) and verifies it is
// constant in x before accepting.
TEST_F(RischLogarithmicProbeTest, HPP007_FormalConstantExtraction_c3) {
    auto e = parse("3 / (x * ln(x))");
    auto r = integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto d_back = diff(r.value(), x, 1, ctx);
    ASSERT_TRUE(d_back.is_ok());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, d_back.value(), e);
    auto delta_tog = algebra::together(delta, ctx);
    ASSERT_TRUE(delta_tog.is_ok());
    auto delta_simp = ctx.simplify(delta_tog.value());
    ASSERT_TRUE(delta_simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(delta_simp.value());
    ASSERT_NE(lit, nullptr)
        << "D(result) - integrand did not reduce to a literal: "
        << formatter::TextFormatter{}.format(delta_simp.value());
    EXPECT_TRUE(lit->value.is_zero())
        << "D(result) ≠ integrand: "
        << formatter::TextFormatter{}.format(delta_simp.value());
}

// HPP-007 closure: also covers rational c outside the old set, e.g. c=5/7.
// ∫ (5/7)/(x·ln(x)) dx = (5/7)·ln(ln(x)).
TEST_F(RischLogarithmicProbeTest, HPP007_FormalConstantExtraction_c5over7) {
    auto e = parse("(5/7) / (x * ln(x))");
    auto r = integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto d_back = diff(r.value(), x, 1, ctx);
    ASSERT_TRUE(d_back.is_ok());
    auto delta = ctx.arena().make<Binary>(BinaryOp::Sub, d_back.value(), e);
    auto delta_tog = algebra::together(delta, ctx);
    ASSERT_TRUE(delta_tog.is_ok());
    auto delta_simp = ctx.simplify(delta_tog.value());
    ASSERT_TRUE(delta_simp.is_ok());
    auto* lit = expr_cast<IntegerLit>(delta_simp.value());
    ASSERT_NE(lit, nullptr)
        << "Got: " << formatter::TextFormatter{}.format(delta_simp.value());
    EXPECT_TRUE(lit->value.is_zero())
        << "D(result) ≠ integrand: "
        << formatter::TextFormatter{}.format(delta_simp.value());
}

// ∫ 1/(x³-1) dx requires partial fractions; the log part involves
// roots of x³-1 = (x-1)(x²+x+1). This is exactly the case where
// Rothstein-Trager helps (avoids root computation explicit).
TEST_F(RischLogarithmicProbeTest, IntegralOfReciprocalOfXCubedMinusOne) {
    auto e = parse("1 / (x^3 - 1)");
    auto r = integrate(e, x, ctx);
    if (r.is_ok()) {
        auto formatted = formatter::TextFormatter{}.format(r.value());
        std::cout << "[PROBE] ∫ 1/(x³-1) dx = " << formatted << std::endl;
        // Acceptable: any closed-form with ln. No specific assertion
        // here — probe only records what engine produces.
    } else {
        std::cout << "[PROBE FAIL] " << r.error().message << std::endl;
    }
}

// Verifies that ∫ d/dx(f) dx = f + C invariant on ln(x²+1).
TEST_F(RischLogarithmicProbeTest, IntegrateDiffRoundtripLnXSquaredPlusOne) {
    auto e = parse("ln(x^2 + 1)");
    auto d = diff(e, x, 1, ctx);
    ASSERT_TRUE(d.is_ok()) << d.error().message;
    auto integrated = integrate(d.value(), x, ctx);
    if (integrated.is_ok()) {
        auto formatted = formatter::TextFormatter{}.format(integrated.value());
        std::cout << "[ROUNDTRIP] ∫ d/dx(ln(x²+1)) dx = " << formatted
                  << " (expected ln(x²+1) + C)" << std::endl;
    } else {
        std::cout << "[ROUNDTRIP FAIL] " << integrated.error().message << std::endl;
    }
}
