// F7.2-T1 — Normal distribution unit tests.

#include "cas/ast.hpp"
#include "cas/numeric.hpp"
#include "cas/statistics.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

namespace cas::statistics {
namespace {

class NormalTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
};

TEST_F(NormalTest, PdfAtZeroStandardNormal) {
    // φ(0; 0, 1) = 1 / sqrt(2π).
    auto x = ctx.arena().make<IntegerLit>(BigInt(0));
    auto mu = ctx.arena().make<IntegerLit>(BigInt(0));
    auto sigma = ctx.arena().make<IntegerLit>(BigInt(1));
    auto pdf = normal_pdf(x, mu, sigma, ctx);
    auto v = numeric::eval(pdf);
    ASSERT_TRUE(v.is_ok());
    const double expected = 1.0 / std::sqrt(2.0 * std::numbers::pi);
    EXPECT_NEAR(v.value(), expected, 1e-12);
}

TEST_F(NormalTest, CdfAtZeroStandardNormal) {
    // Φ(0; 0, 1) = 1/2.
    auto x = ctx.arena().make<IntegerLit>(BigInt(0));
    auto mu = ctx.arena().make<IntegerLit>(BigInt(0));
    auto sigma = ctx.arena().make<IntegerLit>(BigInt(1));
    auto cdf = normal_cdf(x, mu, sigma, ctx);
    auto v = numeric::eval(cdf);
    ASSERT_TRUE(v.is_ok());
    EXPECT_NEAR(v.value(), 0.5, 1e-12);
}

TEST_F(NormalTest, QuantileRoundTripsCdf) {
    // Φ⁻¹(0.975; 0, 1) ≈ 1.959964 (canonical 97.5% z-score).
    auto q = normal_quantile(0.975, 0.0, 1.0);
    ASSERT_TRUE(q.is_ok());
    EXPECT_NEAR(q.value(), 1.959963984540054, 1e-8);
    auto q2 = normal_quantile(0.5, 0.0, 1.0);
    ASSERT_TRUE(q2.is_ok());
    EXPECT_NEAR(q2.value(), 0.0, 1e-8);
}

TEST_F(NormalTest, QuantileBadInput) {
    EXPECT_TRUE(normal_quantile(0.0, 0.0, 1.0).is_error());
    EXPECT_TRUE(normal_quantile(1.0, 0.0, 1.0).is_error());
    EXPECT_TRUE(normal_quantile(0.5, 0.0, 0.0).is_error());
    EXPECT_TRUE(normal_quantile(0.5, 0.0, -1.0).is_error());
}

}  // namespace
}  // namespace cas::statistics
