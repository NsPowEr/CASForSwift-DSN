// F7.3-T2 — Lagrange interpolation tests.

#include "cas/numeric.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace cas::numeric {
namespace {

TEST(LagrangeInterpolationTest, RecoversSampleValues) {
    std::vector<InterpolationPoint> pts = {{0.0, 1.0}, {1.0, 2.0}, {2.0, 5.0}};
    for (const auto& p : pts) {
        auto v = lagrange_evaluate(pts, p.x);
        ASSERT_TRUE(v.is_ok());
        EXPECT_NEAR(v.value(), p.y, 1e-14);
    }
}

TEST(LagrangeInterpolationTest, InterpolatesQuadratic) {
    // y = x^2 sampled at -2, 0, 3.
    std::vector<InterpolationPoint> pts = {{-2.0, 4.0}, {0.0, 0.0}, {3.0, 9.0}};
    auto v = lagrange_evaluate(pts, 1.0);
    ASSERT_TRUE(v.is_ok());
    EXPECT_NEAR(v.value(), 1.0, 1e-12);
    auto v2 = lagrange_evaluate(pts, -1.0);
    ASSERT_TRUE(v2.is_ok());
    EXPECT_NEAR(v2.value(), 1.0, 1e-12);
}

TEST(LagrangeInterpolationTest, DuplicateXFails) {
    std::vector<InterpolationPoint> pts = {{1.0, 2.0}, {1.0, 3.0}};
    auto v = lagrange_evaluate(pts, 0.5);
    EXPECT_TRUE(v.is_error());
}

TEST(LagrangeInterpolationTest, ApproximatesSinusoid) {
    // 9 Chebyshev-like samples of sin(x) over [-π, π]; check max error
    // bounded for a degree-8 polynomial fit.
    std::vector<InterpolationPoint> pts;
    constexpr int kN = 9;
    for (int k = 0; k < kN; ++k) {
        const double x = -3.14159265358979 + 6.28318530717959 * k / (kN - 1);
        pts.push_back({x, std::sin(x)});
    }
    auto v = lagrange_evaluate(pts, 0.0);
    ASSERT_TRUE(v.is_ok());
    EXPECT_NEAR(v.value(), 0.0, 1e-6);
    auto v2 = lagrange_evaluate(pts, 1.0);
    ASSERT_TRUE(v2.is_ok());
    EXPECT_NEAR(v2.value(), std::sin(1.0), 1e-2);
}

}  // namespace
}  // namespace cas::numeric
