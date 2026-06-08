// F7.2-T2 — Linear regression OLS unit tests.

#include "cas/statistics.hpp"

#include <gtest/gtest.h>

namespace cas::statistics {
namespace {

TEST(LinearRegressionTest, PerfectLineFitsExactly) {
    // y = 2x + 1 on 4 samples → slope 2, intercept 1, R²=1.
    std::vector<double> x = {0, 1, 2, 3};
    std::vector<double> y = {1, 3, 5, 7};
    auto res = linear_regression(x, y);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value().slope, 2.0, 1e-12);
    EXPECT_NEAR(res.value().intercept, 1.0, 1e-12);
    EXPECT_NEAR(res.value().r_squared, 1.0, 1e-12);
}

TEST(LinearRegressionTest, FitsNoisyLine) {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2.1, 4.0, 5.9, 8.1, 10.0};
    auto res = linear_regression(x, y);
    ASSERT_TRUE(res.is_ok());
    EXPECT_NEAR(res.value().slope, 1.98, 0.05);
    EXPECT_NEAR(res.value().intercept, 0.1, 0.1);
    EXPECT_GT(res.value().r_squared, 0.99);
}

TEST(LinearRegressionTest, RejectsMismatchedLength) {
    auto res = linear_regression({1, 2}, {1, 2, 3});
    EXPECT_TRUE(res.is_error());
}

TEST(LinearRegressionTest, RejectsZeroVariance) {
    auto res = linear_regression({1, 1, 1}, {2, 3, 4});
    EXPECT_TRUE(res.is_error());
}

}  // namespace
}  // namespace cas::statistics
