// A6 / CAS-L3-18 — Set-resolvent machinery tests (exact, closed-form oracles).

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "../../../src/algebra/galois_setpoly_internal.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"

using namespace cas;
using namespace cas::algebra;
using namespace cas::algebra::galois_setpoly;

namespace {

[[nodiscard]] RatPoly rp(const std::vector<std::int64_t>& low_to_high) {
    std::vector<Rational> c;
    c.reserve(low_to_high.size());
    for (const auto v : low_to_high) c.emplace_back(BigInt(v));
    RatPoly p(std::move(c));
    p.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return p;
}

TEST(GaloisSetpolyTest, TwoSetResolventQuadratic) {
    // f = (x−1)(x−2) = x² − 3x + 2 → R₂ = y − (1+2) = y − 3.
    auto r = two_set_resolvent(rp({2, -3, 1}));
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().coefficients(), rp({-3, 1}).coefficients());
}

TEST(GaloisSetpolyTest, TwoSetResolventCubicSplit) {
    // f = (x−1)(x−2)(x−3) → pair sums {3,4,5}:
    // R₂ = (y−3)(y−4)(y−5) = y³ − 12y² + 47y − 60.
    auto r = two_set_resolvent(rp({-6, 11, -6, 1}));
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().coefficients(), rp({-60, 47, -12, 1}).coefficients());
}

TEST(GaloisSetpolyTest, TwoSetResolventSexticStructural) {
    // f = x⁶ − x − 1 (irreducible, Galois group S₆ by Osada's theorem on
    // trinomials xⁿ−x−1). R₂ has degree C(6,2)=15, is monic, and for this
    // f is squarefree (distinct pair sums).
    auto r = two_set_resolvent(rp({-1, -1, 0, 0, 0, 0, 1}));
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().degree(), 15U);
    EXPECT_TRUE(r.value().leading_coeff() == Rational(BigInt(1)));
    auto sf = is_squarefree_q(r.value());
    ASSERT_TRUE(sf.is_ok());
    EXPECT_TRUE(sf.value());
}

TEST(GaloisSetpolyTest, ThreeSetResolventQuarticSplit) {
    // f = (x−1)(x−2)(x−3)(x−4): triple sums {6,7,8,9} →
    // R₃ = (y−6)(y−7)(y−8)(y−9) = y⁴ − 30y³ + 335y² − 1650y + 3024.
    auto r = three_set_resolvent(rp({24, -50, 35, -10, 1}));
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().coefficients(),
              rp({3024, -1650, 335, -30, 1}).coefficients());
}

TEST(GaloisSetpolyTest, TwoSetResolventQuarticWithCollision) {
    // f = (x−1)(x−2)(x−3)(x−4): pair sums {3,4,5,5,6,7} — the collision
    // 1+4 = 2+3 = 5 must appear as an exact double root:
    // R₂ = (y−3)(y−4)(y−5)²(y−6)(y−7).
    auto r = two_set_resolvent(rp({24, -50, 35, -10, 1}));
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(r.value().degree(), 6U);
    // Expected expansion (Maxima-verified):
    // y⁶ − 30y⁵ + 370y⁴ − 2400y³ + 8629y² − 16290y + 12600.
    EXPECT_EQ(r.value().coefficients(),
              rp({12600, -16290, 8629, -2400, 370, -30, 1}).coefficients());
    auto sf = is_squarefree_q(r.value());
    ASSERT_TRUE(sf.is_ok());
    EXPECT_FALSE(sf.value());
}

TEST(GaloisSetpolyTest, SquarefreeDetectsMultipleRoot) {
    // (y−1)²(y+2) = y³ − 3y + 2 → NOT squarefree.
    auto sf = is_squarefree_q(rp({2, -3, 0, 1}));
    ASSERT_TRUE(sf.is_ok());
    EXPECT_FALSE(sf.value());
    auto sf2 = is_squarefree_q(rp({-6, 11, -6, 1}));
    ASSERT_TRUE(sf2.is_ok());
    EXPECT_TRUE(sf2.value());
}

TEST(GaloisSetpolyTest, TschirnhausCollisionThenInjective) {
    // f = x² − 2 (roots ±√2). c=0: β = α² = 2 twice → degenerate (non-
    // squarefree). c=1: β = 2 ± √2 distinct → g squarefree of degree 2
    // with the same splitting field: g = (y−2)² − 2 = y² − 4y + 2.
    const RatPoly f = rp({-2, 0, 1});
    auto g0 = tschirnhaus_quadratic(f, BigInt(0));
    ASSERT_TRUE(g0.is_ok()) << g0.error().message;
    auto sf0 = is_squarefree_q(g0.value());
    ASSERT_TRUE(sf0.is_ok());
    EXPECT_FALSE(sf0.value());

    auto g1 = tschirnhaus_quadratic(f, BigInt(1));
    ASSERT_TRUE(g1.is_ok()) << g1.error().message;
    EXPECT_EQ(g1.value().coefficients(), rp({2, -4, 1}).coefficients());
    auto sf1 = is_squarefree_q(g1.value());
    ASSERT_TRUE(sf1.is_ok());
    EXPECT_TRUE(sf1.value());
}

}  // namespace
