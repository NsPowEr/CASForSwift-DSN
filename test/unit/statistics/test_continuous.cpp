// F7.2-T4 — Chi-squared, Student-t, F-distribution tests.

#include "cas/statistics.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

namespace cas::statistics {
namespace {

// ── Chi-squared ──────────────────────────────────────────────────────────────

TEST(ChiSquaredTest, PdfAtOneDfThree) {
    // χ²_3 pdf at x=1: x^(1/2) · e^{-1/2} / (2^{3/2} · Γ(3/2)).
    // Γ(3/2) = √π/2 → denominator = 2^{3/2} · √π/2 = 2^{1/2}·√π.
    auto v = chi_squared_pdf(1.0, 3.0);
    ASSERT_TRUE(v.is_ok());
    const double expected = std::exp(-0.5) / std::sqrt(2.0 * std::numbers::pi);
    EXPECT_NEAR(v.value(), expected, 1e-12);
}

TEST(ChiSquaredTest, CdfKnownQuantile) {
    // For χ²_5, the 95% quantile is x ≈ 11.0705 (canonical statistical table).
    auto c = chi_squared_cdf(11.0705, 5.0);
    ASSERT_TRUE(c.is_ok());
    EXPECT_NEAR(c.value(), 0.95, 1e-4);
}

TEST(ChiSquaredTest, CdfMonotone) {
    double prev = 0.0;
    for (double x = 0.5; x < 20.0; x += 0.5) {
        auto c = chi_squared_cdf(x, 4.0);
        ASSERT_TRUE(c.is_ok());
        EXPECT_GE(c.value(), prev - 1e-12);
        prev = c.value();
    }
    EXPECT_NEAR(prev, 1.0, 1e-2);
}

// ── Student-t ────────────────────────────────────────────────────────────────

TEST(StudentTTest, PdfSymmetry) {
    for (double x : {0.3, 1.0, 2.0, 4.0}) {
        auto lhs = student_t_pdf(x, 5.0);
        auto rhs = student_t_pdf(-x, 5.0);
        ASSERT_TRUE(lhs.is_ok());
        ASSERT_TRUE(rhs.is_ok());
        EXPECT_NEAR(lhs.value(), rhs.value(), 1e-12);
    }
}

TEST(StudentTTest, CdfMedianIsHalf) {
    auto c = student_t_cdf(0.0, 4.0);
    ASSERT_TRUE(c.is_ok());
    EXPECT_NEAR(c.value(), 0.5, 1e-12);
}

TEST(StudentTTest, CdfKnownQuantile95Pct) {
    // For t with 10 df, the 95% one-sided critical value is t = 1.8125.
    auto c = student_t_cdf(1.8125, 10.0);
    ASSERT_TRUE(c.is_ok());
    EXPECT_NEAR(c.value(), 0.95, 1e-4);
}

// ── F-distribution ───────────────────────────────────────────────────────────

TEST(FDistributionTest, PdfZeroOutsideSupport) {
    EXPECT_NEAR(f_pdf(-1.0, 3.0, 4.0).value(), 0.0, 1e-12);
    EXPECT_NEAR(f_pdf(0.0, 3.0, 4.0).value(), 0.0, 1e-12);
}

TEST(FDistributionTest, CdfKnownPercentile) {
    // F(3, 30) 95% quantile ≈ 2.9223 (statistical table).
    auto c = f_cdf(2.9223, 3.0, 30.0);
    ASSERT_TRUE(c.is_ok());
    EXPECT_NEAR(c.value(), 0.95, 1e-4);
}

TEST(FDistributionTest, RejectsBadDfs) {
    EXPECT_TRUE(f_pdf(1.0, -1.0, 5.0).is_error());
    EXPECT_TRUE(f_pdf(1.0, 3.0, 0.0).is_error());
    EXPECT_TRUE(chi_squared_pdf(1.0, 0.0).is_error());
    EXPECT_TRUE(student_t_pdf(1.0, -2.0).is_error());
}

}  // namespace
}  // namespace cas::statistics
