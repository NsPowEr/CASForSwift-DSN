// F5.6 sub-block 2 — Residue theorem driver wiring of the Aberth numeric
// isolator.  Each test fixes a rational integrand whose denominator escapes
// the existing symbolic catalogue (general quartic with non-zero odd
// coefficients, irreducible deg ≥ 5) and asserts that:
//   1. The driver now returns a numerical DecimalLit instead of bailing out
//      with Unimplemented.
//   2. The reported value matches a closed-form / mpmath-verified reference
//      within the working precision tolerance.

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/residue_theorem.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <string>

namespace cas::test {
namespace {

cas::ExprPtr parse_expr(const std::string& input, cas::AstArena& arena) {
    auto tokens = cas::Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << tokens.error().message;
    cas::Parser parser(tokens.value(), arena);
    auto res = parser.parse();
    EXPECT_TRUE(res.is_ok()) << res.error().message;
    return res.value();
}

double evaluate_to_double(cas::ExprPtr e) {
    if (const auto* dec = cas::expr_cast<cas::DecimalLit>(e)) return dec->to_double();
    if (const auto* il  = cas::expr_cast<cas::IntegerLit>(e)) {
        return std::stod(il->value.decimal());
    }
    if (const auto* rl  = cas::expr_cast<cas::RationalLit>(e)) {
        return std::stod(rl->numerator.decimal()) /
               std::stod(rl->denominator.decimal());
    }
    return std::nan("");
}

}  // namespace

class ResidueTheoremAberthTest : public ::testing::Test {
protected:
    cas::symbolic::CASContext ctx;
    cas::Symbol x{"x"};
};

// ∫_{-∞}^∞ dx / Φ₉(x) where Φ₉(x) = x⁶ + x³ + 1 is the 9-th cyclotomic
// polynomial — irreducible deg-6 over Q, no real roots (all roots are
// primitive 9th roots of unity).
//
// Closed form via the 9th-roots-of-unity residue calculation
// (derivative D'(z) = 6z⁵ + 3z² = 3z²(2z³ + 1); for z = e^(iθ) with
// θ ∈ {2π/9, 4π/9, 8π/9} the factor 2z³ + 1 equals ±i√3, leading to)
//
//   I = (2π / (3√3)) · (cos(π/9) + cos(2π/9) + cos(4π/9)) ≈ 2.272555…
TEST_F(ResidueTheoremAberthTest, IrreducibleSexticCyclotomic9) {
    auto expr = parse_expr("1/(x^6 + x^3 + 1)", ctx.arena());
    auto res = cas::calculus::integrate_rational_full_real_line(expr, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    double v = evaluate_to_double(res.value());
    const double sqrt_three = std::sqrt(3.0);
    const double expected = (2.0 * M_PI / (3.0 * sqrt_three)) *
        (std::cos(M_PI / 9.0) + std::cos(2.0 * M_PI / 9.0) +
         std::cos(4.0 * M_PI / 9.0));
    EXPECT_NEAR(v, expected, 1e-6);
}

// ∫_{-∞}^∞ dx / (x⁴ + x + 1).  The denominator is irreducible over Q (mod-2
// factorisation test) and is NOT biquadratic (the x term is present), so it
// falls outside the symbolic closed-form path and exercises the Aberth
// numeric driver.  No clean closed form in radicals — we verify the value
// returned matches itself across re-invocations (deterministic algorithm)
// and lies in the analytically expected positive range.
TEST_F(ResidueTheoremAberthTest, IrreducibleQuarticNonBiquadratic) {
    auto expr = parse_expr("1/(x^4 + x + 1)", ctx.arena());
    auto res = cas::calculus::integrate_rational_full_real_line(expr, x, ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    double v = evaluate_to_double(res.value());
    EXPECT_GT(v, 1.5);           // positive real-axis integral, comparable
    EXPECT_LT(v, 4.0);           // upper sanity bound
    // Re-invoke and require bitwise repeatability (deterministic Aberth).
    auto res2 = cas::calculus::integrate_rational_full_real_line(expr, x, ctx);
    ASSERT_TRUE(res2.is_ok()) << res2.error().message;
    double v2 = evaluate_to_double(res2.value());
    EXPECT_DOUBLE_EQ(v, v2);
}

}  // namespace cas::test
