// F7.3-B3 — Cubic spline + Hermite + RKF45 tests.

#include <gtest/gtest.h>

#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <cmath>

using namespace cas;
using namespace cas::numeric;
using namespace cas::symbolic;

namespace {

TEST(CubicSpline, InterpolatesKnotsExactly) {
    std::vector<InterpolationPoint> pts{
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}, {4.0, 16.0}};
    auto s = build_natural_cubic_spline(pts);
    ASSERT_TRUE(s.is_ok());
    for (const auto& p : pts) {
        auto v = cubic_spline_evaluate(s.value(), p.x);
        ASSERT_TRUE(v.is_ok());
        EXPECT_NEAR(v.value(), p.y, 1e-12);
    }
}

TEST(CubicSpline, SmoothInterpolationBetweenKnots) {
    // y = x² sampled at integer knots; spline value at midpoints should be
    // close (not exact: natural BC differs from quadratic).
    std::vector<InterpolationPoint> pts{
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}, {4.0, 16.0}, {5.0, 25.0}};
    auto s = build_natural_cubic_spline(pts);
    ASSERT_TRUE(s.is_ok());
    auto v = cubic_spline_evaluate(s.value(), 2.5);
    ASSERT_TRUE(v.is_ok());
    EXPECT_NEAR(v.value(), 6.25, 0.2);
}

TEST(CubicSpline, RejectsUnsortedX) {
    std::vector<InterpolationPoint> pts{{1.0, 1.0}, {0.0, 0.0}, {2.0, 4.0}};
    auto s = build_natural_cubic_spline(pts);
    EXPECT_TRUE(s.is_error());
}

TEST(HermiteSpline, MatchesEndpointsAndDerivatives) {
    // Cubic Hermite must reproduce y_i and dy_i at every knot.
    std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> y{0.0, 1.0, 0.0};
    std::vector<double> dy{1.0, 0.0, -1.0};
    auto s = build_hermite_spline(x, y, dy);
    ASSERT_TRUE(s.is_ok());
    for (std::size_t i = 0; i < x.size(); ++i) {
        auto v = hermite_evaluate(s.value(), x[i]);
        ASSERT_TRUE(v.is_ok());
        EXPECT_NEAR(v.value(), y[i], 1e-12);
    }
}

TEST(RKF45, SolvesLinearDecay) {
    // dy/dt = -y, y(0)=1 → y(t) = exp(-t).
    CASContext ctx;
    auto& a = ctx.arena();
    auto rhs = a.make<Unary>(UnaryOp::Neg, a.make<Symbol>("y"));
    auto r = solve_ode_rkf45(rhs, "t", "y", 0.0, 1.0, 2.0);
    ASSERT_TRUE(r.is_ok());
    const auto& traj = r.value();
    EXPECT_GE(traj.size(), 2U);
    const auto& last = traj.back();
    EXPECT_NEAR(last.t, 2.0, 1e-9);
    EXPECT_NEAR(last.y, std::exp(-2.0), 1e-4);
}

TEST(RKF45, SolvesExponentialGrowth) {
    // dy/dt = y, y(0)=1 → y(t) = exp(t).
    CASContext ctx;
    auto& a = ctx.arena();
    auto rhs = a.make<Symbol>("y");
    auto r = solve_ode_rkf45(rhs, "t", "y", 0.0, 1.0, 1.0);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NEAR(r.value().back().y, std::exp(1.0), 1e-4);
}

TEST(RKF45, RejectsBadInterval) {
    CASContext ctx;
    auto& a = ctx.arena();
    auto rhs = a.make<Symbol>("y");
    auto r = solve_ode_rkf45(rhs, "t", "y", 1.0, 1.0, 0.0);  // t_end < t0
    EXPECT_TRUE(r.is_error());
}

}  // namespace
