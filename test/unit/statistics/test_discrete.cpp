// F7.2-T3 — Binomial + Poisson distribution tests.

#include "cas/statistics.hpp"

#include <gtest/gtest.h>

#include <cmath>

namespace cas::statistics {
namespace {

TEST(BinomialTest, KnownPmfValues) {
    // Bin(10, 0.3): P(X = 3) = C(10,3) · 0.3³ · 0.7⁷ = 0.266827932.
    auto v = binomial_pmf(3, 10, 0.3);
    ASSERT_TRUE(v.is_ok());
    EXPECT_NEAR(v.value(), 0.266827932, 1e-9);

    // Boundary: P(X = 0) when p = 0 is 1.0; all other k give 0.
    EXPECT_NEAR(binomial_pmf(0, 5, 0.0).value(), 1.0, 1e-12);
    EXPECT_NEAR(binomial_pmf(3, 5, 0.0).value(), 0.0, 1e-12);
    EXPECT_NEAR(binomial_pmf(5, 5, 1.0).value(), 1.0, 1e-12);
}

TEST(BinomialTest, CdfMonotone) {
    double prev = 0.0;
    for (long long k = 0; k <= 20; ++k) {
        auto c = binomial_cdf(k, 20, 0.4);
        ASSERT_TRUE(c.is_ok());
        EXPECT_GE(c.value(), prev);
        prev = c.value();
    }
    EXPECT_NEAR(prev, 1.0, 1e-10);
}

TEST(BinomialTest, RejectsBadParameters) {
    EXPECT_TRUE(binomial_pmf(2, -1, 0.5).is_error());
    EXPECT_TRUE(binomial_pmf(2, 10, -0.1).is_error());
    EXPECT_TRUE(binomial_pmf(2, 10, 1.5).is_error());
    EXPECT_NEAR(binomial_pmf(-1, 10, 0.5).value(), 0.0, 1e-12);
    EXPECT_NEAR(binomial_pmf(11, 10, 0.5).value(), 0.0, 1e-12);
}

TEST(PoissonTest, KnownPmfValues) {
    // Pois(2): P(X = 0) = e^{-2}, P(X = 1) = 2·e^{-2}, P(X = 2) = 2·e^{-2}.
    const double e_neg2 = std::exp(-2.0);
    EXPECT_NEAR(poisson_pmf(0, 2.0).value(), e_neg2, 1e-12);
    EXPECT_NEAR(poisson_pmf(1, 2.0).value(), 2.0 * e_neg2, 1e-12);
    EXPECT_NEAR(poisson_pmf(2, 2.0).value(), 2.0 * e_neg2, 1e-12);
}

TEST(PoissonTest, CdfApproachesOne) {
    auto c = poisson_cdf(30, 2.0);
    ASSERT_TRUE(c.is_ok());
    EXPECT_NEAR(c.value(), 1.0, 1e-10);
}

TEST(PoissonTest, ZeroLambdaIsDeterministic) {
    EXPECT_NEAR(poisson_pmf(0, 0.0).value(), 1.0, 1e-12);
    EXPECT_NEAR(poisson_pmf(5, 0.0).value(), 0.0, 1e-12);
    EXPECT_NEAR(poisson_cdf(0, 0.0).value(), 1.0, 1e-12);
}

}  // namespace
}  // namespace cas::statistics
