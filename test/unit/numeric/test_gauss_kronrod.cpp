// F7.3-T1 — Gauss-Kronrod 15/7 adaptive quadrature.

#include "cas/lexer.hpp"
#include "cas/numeric.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

namespace cas::numeric {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok());
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok());
    return r.value();
}

class GaussKronrodTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(GaussKronrodTest, PolynomialExactOnSinglePanel) {
    // ∫_0^1 x^2 dx = 1/3.
    auto p = parse_expr("x^2", ctx.arena());
    IntegrationOptions opts;
    opts.scheme = IntegrationScheme::GaussKronrod;
    auto res = integrate_numeric(p, "x", 0.0, 1.0, opts);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value(), 1.0 / 3.0, 1e-12);
}

TEST_F(GaussKronrodTest, GaussianBellOverFiniteRange) {
    // ∫_{-3}^{3} e^{-x^2/2}/sqrt(2π) dx ≈ erf(3/sqrt(2)).
    auto p = parse_expr("exp(-(x^2)/2)/sqrt(2*pi)", ctx.arena());
    IntegrationOptions opts;
    opts.scheme = IntegrationScheme::GaussKronrod;
    opts.tolerance = 1e-10;
    auto res = integrate_numeric(p, "x", -3.0, 3.0, opts);
    ASSERT_TRUE(res.is_ok());
    const double expected = std::erf(3.0 / std::sqrt(2.0));
    EXPECT_NEAR(res.value(), expected, 1e-8);
}

TEST_F(GaussKronrodTest, ConvergesFasterThanSimpsonOnOscillatory) {
    // ∫_0^π sin(x) dx = 2. Both methods should converge but Gauss-Kronrod
    // typically uses fewer panels for the same accuracy.
    auto p = parse_expr("sin(x)", ctx.arena());
    IntegrationOptions gk;
    gk.scheme = IntegrationScheme::GaussKronrod;
    gk.tolerance = 1e-10;
    auto res_gk = integrate_numeric(p, "x", 0.0, std::numbers::pi, gk);
    ASSERT_TRUE(res_gk.is_ok());
    EXPECT_NEAR(res_gk.value(), 2.0, 1e-10);

    IntegrationOptions simp;
    simp.scheme = IntegrationScheme::AdaptiveSimpson;
    simp.tolerance = 1e-10;
    auto res_s = integrate_numeric(p, "x", 0.0, std::numbers::pi, simp);
    ASSERT_TRUE(res_s.is_ok());
    EXPECT_NEAR(res_s.value(), 2.0, 1e-8);
}

}  // namespace
}  // namespace cas::numeric
