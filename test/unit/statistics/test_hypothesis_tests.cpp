// F7.2-B2 — Hypothesis tests + multivariate OLS.

#include <gtest/gtest.h>

#include "cas/statistics.hpp"

#include <cmath>

using namespace cas;
using namespace cas::statistics;

namespace {

TEST(ZTestOneSample, RejectsWhenMeanFarFromMu0) {
    std::vector<double> sample{10.0, 10.5, 11.0, 9.8, 10.2, 10.7, 10.3, 9.9};
    auto r = z_test_one_sample(sample, 0.0, 1.0);
    ASSERT_TRUE(r.is_ok());
    EXPECT_GT(r.value().statistic, 5.0);
    EXPECT_LT(r.value().p_value, 0.001);
    EXPECT_TRUE(r.value().reject_at_alpha_005);
}

TEST(ZTestOneSample, DoesNotRejectWhenNullTrue) {
    std::vector<double> sample{0.05, -0.1, 0.02, -0.03, 0.08, -0.04, 0.01};
    auto r = z_test_one_sample(sample, 0.0, 1.0);
    ASSERT_TRUE(r.is_ok());
    EXPECT_FALSE(r.value().reject_at_alpha_005);
}

TEST(TTestOneSample, RejectsLargeDeparture) {
    std::vector<double> sample{10.0, 10.5, 11.0, 9.8, 10.2, 10.7, 10.3, 9.9};
    auto r = t_test_one_sample(sample, 0.0);
    ASSERT_TRUE(r.is_ok());
    EXPECT_GT(std::abs(r.value().statistic), 30.0);
    EXPECT_LT(r.value().p_value, 0.001);
    EXPECT_TRUE(r.value().reject_at_alpha_005);
}

TEST(TTestTwoSample, DetectsMeanDifference) {
    std::vector<double> s1{10.0, 10.5, 11.0, 9.8, 10.2};
    std::vector<double> s2{0.0, 0.5, 1.0, -0.2, 0.3};
    auto r = t_test_two_sample(s1, s2);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().reject_at_alpha_005);
}

TEST(ChiSquaredGoodnessOfFit, FairDieAccepted) {
    // 60 rolls of a fair die, ~10 each — expect not to reject uniformity.
    std::vector<double> observed{12, 8, 11, 9, 10, 10};
    std::vector<double> expected{10, 10, 10, 10, 10, 10};
    auto r = chi_squared_goodness_of_fit(observed, expected);
    ASSERT_TRUE(r.is_ok());
    EXPECT_FALSE(r.value().reject_at_alpha_005);
}

TEST(ChiSquaredGoodnessOfFit, BiasedDieRejected) {
    std::vector<double> observed{30, 5, 5, 5, 5, 10};
    std::vector<double> expected{10, 10, 10, 10, 10, 10};
    auto r = chi_squared_goodness_of_fit(observed, expected);
    ASSERT_TRUE(r.is_ok());
    EXPECT_TRUE(r.value().reject_at_alpha_005);
}

TEST(FTestVariance, EqualVariancesAccepted) {
    std::vector<double> s1{1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> s2{2.0, 3.0, 4.0, 5.0, 6.0};
    auto r = f_test_variance(s1, s2);
    ASSERT_TRUE(r.is_ok());
    EXPECT_FALSE(r.value().reject_at_alpha_005);
}

TEST(MultivariateOLS, RecoversLinearCombinationExactly) {
    // y = 1 + 2·x1 + 3·x2 + noise=0
    std::vector<std::vector<double>> X{
        {1.0, 1.0},
        {2.0, 1.0},
        {1.0, 2.0},
        {3.0, 2.0},
        {2.0, 3.0},
        {4.0, 4.0},
        {1.0, 5.0},
    };
    std::vector<double> y;
    for (const auto& row : X) {
        y.push_back(1.0 + 2.0 * row[0] + 3.0 * row[1]);
    }
    auto r = multivariate_linear_regression(X, y);
    ASSERT_TRUE(r.is_ok());
    ASSERT_EQ(r.value().coefficients.size(), 3U);
    EXPECT_NEAR(r.value().coefficients[0], 1.0, 1e-9);
    EXPECT_NEAR(r.value().coefficients[1], 2.0, 1e-9);
    EXPECT_NEAR(r.value().coefficients[2], 3.0, 1e-9);
    EXPECT_NEAR(r.value().r_squared, 1.0, 1e-12);
    EXPECT_NEAR(r.value().residual_sum_of_squares, 0.0, 1e-12);
}

TEST(MultivariateOLS, DimensionMismatchRejected) {
    std::vector<std::vector<double>> X{{1.0, 2.0}};
    std::vector<double> y{1.0, 2.0};
    auto r = multivariate_linear_regression(X, y);
    EXPECT_TRUE(r.is_error());
}

}  // namespace
