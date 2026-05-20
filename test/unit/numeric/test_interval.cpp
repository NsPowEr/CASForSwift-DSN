// CAS-L3-13 — Interval arithmetic tests.

#include <gtest/gtest.h>

#include "cas/interval.hpp"

using namespace cas::numeric;
using namespace cas;

namespace {

[[nodiscard]] BigFloat bf(double v) { return BigFloat::from_double(v); }

TEST(IntervalTest, ConstructionAndAccessors) {
    Interval i(bf(1.0), bf(2.0));
    EXPECT_FALSE(i.lo() < bf(1.0));
    EXPECT_FALSE(bf(1.0) < i.lo());
    EXPECT_FALSE(i.hi() < bf(2.0));
    EXPECT_FALSE(bf(2.0) < i.hi());
}

TEST(IntervalTest, NormalizesSwappedEndpoints) {
    Interval i(bf(5.0), bf(1.0));
    EXPECT_FALSE(bf(1.0) < i.lo());  // lo == 1
    EXPECT_FALSE(i.hi() < bf(5.0));  // hi == 5
}

TEST(IntervalTest, SingletonContainsItself) {
    Interval i(bf(3.14));
    EXPECT_TRUE(i.contains(bf(3.14)));
    EXPECT_FALSE(i.contains(bf(3.15)));
}

TEST(IntervalTest, AdditionMonotone) {
    Interval a(bf(1.0), bf(2.0));
    Interval b(bf(3.0), bf(4.0));
    Interval s = a + b;
    // [1+3, 2+4] = [4, 6]
    EXPECT_TRUE(s.contains(bf(5.0)));
    EXPECT_FALSE(s.contains(bf(3.9)));
    EXPECT_FALSE(s.contains(bf(6.1)));
}

TEST(IntervalTest, MultiplicationSignAware) {
    // [-2, 1] · [-3, 4] = [min(-2·-3,-2·4,1·-3,1·4), max(...)]
    //                   = [min(6,-8,-3,4), max(...)] = [-8, 6]
    Interval a(bf(-2.0), bf(1.0));
    Interval b(bf(-3.0), bf(4.0));
    Interval p = a * b;
    EXPECT_TRUE(p.contains(bf(0.0)));
    EXPECT_TRUE(p.contains(bf(-7.0)));
    EXPECT_TRUE(p.contains(bf(5.0)));
    EXPECT_FALSE(p.contains(bf(-9.0)));
    EXPECT_FALSE(p.contains(bf(7.0)));
}

TEST(IntervalTest, NegationFlipsBounds) {
    Interval a(bf(-3.0), bf(2.0));
    Interval n = -a;
    EXPECT_TRUE(n.contains(bf(-2.0)));
    EXPECT_TRUE(n.contains(bf(3.0)));
    EXPECT_FALSE(n.contains(bf(-3.1)));
}

TEST(IntervalTest, DivisionAvoidsZero) {
    Interval a(bf(1.0), bf(2.0));
    Interval b(bf(2.0), bf(4.0));
    Interval q = a / b;
    // [1/4, 2/2] = [0.25, 1.0]
    EXPECT_TRUE(q.contains(bf(0.5)));
    EXPECT_FALSE(q.contains(bf(0.1)));
    EXPECT_FALSE(q.contains(bf(1.1)));
}

TEST(IntervalTest, DivisionByZeroWidens) {
    Interval a(bf(1.0), bf(2.0));
    Interval b(bf(-1.0), bf(1.0));  // contains zero
    Interval q = a / b;
    // Wide bound — must contain anything large.
    EXPECT_TRUE(q.contains(bf(1e10)));
    EXPECT_TRUE(q.contains(bf(-1e10)));
}

TEST(IntervalTest, SqrtOnPositiveInterval) {
    Interval i(bf(4.0), bf(9.0));
    Interval s = Interval::sqrt(i);
    EXPECT_TRUE(s.contains(bf(2.5)));
    EXPECT_FALSE(s.contains(bf(1.9)));
    EXPECT_FALSE(s.contains(bf(3.1)));
}

TEST(IntervalTest, ExpMonotone) {
    Interval i(bf(0.0), bf(1.0));
    Interval e = Interval::exp(i);
    // [exp(0), exp(1)] = [1, e ≈ 2.718]
    EXPECT_TRUE(e.contains(bf(1.0)));
    EXPECT_TRUE(e.contains(bf(2.5)));
    EXPECT_FALSE(e.contains(bf(0.9)));
    EXPECT_FALSE(e.contains(bf(3.0)));
}

TEST(IntervalTest, AntiHardcodeBoundPropagation) {
    // f(x) = x² + 1 on x ∈ [-1, 2]
    // f([-1, 2]) ⊆ Interval calculation
    Interval x(bf(-1.0), bf(2.0));
    Interval x_sq = x * x;
    Interval one(bf(1.0));
    Interval f = x_sq + one;
    // Actual range of x² on [-1, 2] = [0, 4] (min at 0, max at 2)
    // But interval arithmetic x·x = [min(-1·-1, -1·2, 2·-1, 2·2), max(...)]
    //                              = [min(1,-2,-2,4), max(1,-2,-2,4)]
    //                              = [-2, 4]
    // f = x_sq + 1 = [-1, 5]. Real f range = [1, 5].
    // Anti-hardcode: f must CONTAIN [1, 5] (over-approximation acceptable).
    EXPECT_TRUE(f.contains(bf(1.0)));
    EXPECT_TRUE(f.contains(bf(5.0)));
    EXPECT_TRUE(f.contains(bf(3.0)));
}

TEST(IntervalTest, SinCriticalPointDetection) {
    // [0, π] contains π/2 (max sin = 1). Expected: hi == 1, lo = min(sin(0), sin(π)) = 0.
    Interval x(bf(0.0), BigFloat::pi());
    Interval s = Interval::sin(x);
    EXPECT_TRUE(s.contains(bf(1.0))) << "should contain critical max 1";
    EXPECT_TRUE(s.contains(bf(0.5)));
}

TEST(IntervalTest, SinMonotonicNoExpansion) {
    // [0, π/4]: monotone, no critical point inside. Result [sin(0), sin(π/4)] = [0, ~0.707].
    Interval x(bf(0.0), BigFloat::pi() / bf(4.0));
    Interval s = Interval::sin(x);
    EXPECT_TRUE(s.contains(bf(0.3)));
    EXPECT_FALSE(s.contains(bf(0.8)))  // 0.707 < 0.8
        << "no critical point → tight bound";
}

TEST(IntervalTest, IsPositiveNegativeZero) {
    Interval pos(bf(1.0), bf(3.0));
    Interval neg(bf(-3.0), bf(-1.0));
    Interval mixed(bf(-1.0), bf(1.0));
    EXPECT_TRUE(pos.is_positive());
    EXPECT_FALSE(pos.is_negative());
    EXPECT_FALSE(pos.contains_zero());
    EXPECT_TRUE(neg.is_negative());
    EXPECT_TRUE(mixed.contains_zero());
    EXPECT_FALSE(mixed.is_positive());
    EXPECT_FALSE(mixed.is_negative());
}

}  // namespace
