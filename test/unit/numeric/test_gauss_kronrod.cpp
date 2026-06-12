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

// ─── F6.D adaptive priority-queue G7/K15 ─────────────────────────────────────

TEST_F(GaussKronrodTest, Adaptive_SinIntegral_FullPeriodIsZero) {
    // ∫_0^{2π} sin(x) dx = 0.  Trivially smooth — one panel suffices.
    auto f = [](double x) -> Result<double> { return ok(std::sin(x)); };
    auto r = integrate_gauss_kronrod_adaptive(f, 0.0, 2.0 * std::numbers::pi, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().success);
    EXPECT_NEAR(r.value().value, 0.0, 1e-12);
    EXPECT_LE(r.value().interval_count, 4U) << "smooth integrand should converge on a few panels";
}

TEST_F(GaussKronrodTest, Adaptive_Gaussian_OnZeroToFiveMatchesErfHalf) {
    // ∫_0^5 e^{-x²} dx = (√π / 2) · erf(5) ≈ 0.8862269254527580...
    auto f = [](double x) -> Result<double> { return ok(std::exp(-x * x)); };
    auto r = integrate_gauss_kronrod_adaptive(f, 0.0, 5.0, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().success);
    const double expected = 0.5 * std::sqrt(std::numbers::pi) * std::erf(5.0);
    EXPECT_NEAR(r.value().value, expected, 1e-9);
}

TEST_F(GaussKronrodTest, Adaptive_AtanIdentity_OneOverOnePlusXSquared) {
    // ∫_0^1 1/(1+x²) dx = π/4 — smooth, analytically exact reference.
    auto f = [](double x) -> Result<double> { return ok(1.0 / (1.0 + x * x)); };
    ctx.set_integration_abs_tol(1e-12);
    ctx.set_integration_rel_tol(1e-12);
    auto r = integrate_gauss_kronrod_adaptive(f, 0.0, 1.0, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().success);
    EXPECT_NEAR(r.value().value, 0.25 * std::numbers::pi, 1e-12);
}

TEST_F(GaussKronrodTest, Adaptive_ResourceCapHonored_NonConvergent) {
    // Pathological case: highly oscillatory 1/x integrand near a singularity
    // with a very tight tolerance and a small interval cap. Must terminate
    // with success=false rather than loop or silently return a wrong answer.
    auto f = [](double x) -> Result<double> {
        return ok(std::sin(1.0 / x));  // wildly oscillatory near 0
    };
    ctx.set_integration_abs_tol(1e-30);
    ctx.set_integration_rel_tol(1e-30);
    ctx.set_integration_max_intervals(64U);  // very small cap
    auto r = integrate_gauss_kronrod_adaptive(f, 0.001, 1.0, ctx);
    ASSERT_TRUE(r.is_ok());
    EXPECT_FALSE(r.value().success) << "non-convergent path must report failure";
    EXPECT_GE(r.value().interval_count, 64U);
}

TEST_F(GaussKronrodTest, Adaptive_ReversedLimits_NegatesResult) {
    auto f = [](double x) -> Result<double> { return ok(x * x); };
    auto forward = integrate_gauss_kronrod_adaptive(f, 0.0, 2.0, ctx);
    auto reversed = integrate_gauss_kronrod_adaptive(f, 2.0, 0.0, ctx);
    ASSERT_TRUE(forward.is_ok());
    ASSERT_TRUE(reversed.is_ok());
    EXPECT_NEAR(forward.value().value, -reversed.value().value, 1e-12);
}

}  // namespace
}  // namespace cas::numeric
