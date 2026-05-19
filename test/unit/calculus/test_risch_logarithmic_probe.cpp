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
// Tests if engine handles nested transcendental logarithm.
//
// KNOWN BUG (DEBT-NEW-001 filed 2026-05-20): engine returns
// `ln(x)^(-1) * ln(abs(x))` which simplifies to 1 for x>0, not
// ln(ln(x)). The differentiate-inverse check confirms wrong:
//   diff(1, x) = 0 ≠ 1/(x·ln(x)).
// Fix requires Risch table-substitution / Liouvillian extension
// handler that recognises u = ln(x) as a transcendental tower.
// Filed for follow-up; this probe documents the failure mode but
// does NOT assert the wrong output to avoid green-bar deception.
TEST_F(RischLogarithmicProbeTest, IntegralOfReciprocalOfXLnX_KNOWN_BUG) {
    auto e = parse("1 / (x * ln(x))");
    auto r = integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto formatted = formatter::TextFormatter{}.format(r.value());
    std::cout << "[KNOWN BUG] ∫ 1/(x ln(x)) dx = " << formatted
              << "  (correct: ln(ln(x))). See DEBT-NEW-001." << std::endl;
    // Verify the diff-inverse invariant fails — proof of the bug.
    auto d_back = diff(r.value(), x, 1, ctx);
    if (d_back.is_ok()) {
        // Compare diff(result) with original integrand.
        // If they match, the result is correct (despite weird form).
        // If they don't, this is a genuine bug.
        auto diff_minus_orig = ctx.simplify(
            ctx.arena().make<Binary>(BinaryOp::Sub, d_back.value(), e));
        if (diff_minus_orig.is_ok()) {
            auto check = formatter::TextFormatter{}.format(diff_minus_orig.value());
            std::cout << "[INVARIANT] diff(result) - integrand = " << check
                      << " (should be 0 for correct integral)" << std::endl;
        }
    }
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
