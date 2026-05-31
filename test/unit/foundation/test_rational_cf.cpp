// F1.2-NEW — Rational continued-fraction tests.
//
// Property: the convergent reconstructed from [a0; a1, ..., ak] equals
// the original rational (exact round-trip for non-truncated expansions).
// Edge cases: integer rationals, negative rationals, n_max truncation.

#include <gtest/gtest.h>

#include "cas/rational.hpp"
#include "cas/bigint.hpp"

#include <cstddef>
#include <limits>
#include <vector>

using namespace cas;

namespace {

// -------------------------------------------------------------------------
// Helper: reconstruct Rational from continued-fraction coefficients.
// convergent([a0; a1, ..., ak]) = a0 + 1/(a1 + 1/(a2 + ...))
// Computed bottom-up to avoid deep recursion.
// -------------------------------------------------------------------------
[[nodiscard]] Rational reconstruct_cf(const std::vector<BigInt>& cf) {
    if (cf.empty()) return Rational{};
    Rational result(cf.back());
    for (std::size_t i = cf.size() - 1; i > 0; --i) {
        // result = cf[i-1] + 1/result
        // = (cf[i-1] * result.denom + result.numer) / result.denom
        // Use Rational arithmetic: 1/result = denom/numer, then add cf[i-1].
        Rational inv(result.denominator(), result.numerator());  // 1/result
        result = Rational(cf[i - 1]) + inv;
    }
    return result;
}

// -------------------------------------------------------------------------
// Round-trip test: CF expansion → reconstruct → must equal original.
// -------------------------------------------------------------------------
void check_roundtrip(const Rational& r) {
    auto cf = r.to_continued_fraction();
    ASSERT_FALSE(cf.empty()) << "CF of non-zero rational must be non-empty";
    Rational reconstructed = reconstruct_cf(cf);
    EXPECT_EQ(reconstructed, r)
        << "CF round-trip failed for rational "
        << r.numerator().decimal() << "/" << r.denominator().decimal();
}

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

TEST(RationalCFTest, IntegerIsMonoElement) {
    // An integer has CF [n] (single element).
    Rational r(BigInt(7));
    auto cf = r.to_continued_fraction();
    ASSERT_EQ(cf.size(), 1U);
    EXPECT_EQ(cf[0], BigInt(7));
}

TEST(RationalCFTest, OneHalf) {
    // 1/2 → [0; 2]
    Rational r(BigInt(1), BigInt(2));
    auto cf = r.to_continued_fraction();
    ASSERT_EQ(cf.size(), 2U);
    EXPECT_EQ(cf[0], BigInt(0));
    EXPECT_EQ(cf[1], BigInt(2));
    check_roundtrip(r);
}

TEST(RationalCFTest, ThreeSeventh) {
    // 3/7 → [0; 2, 3]  (since 7/3 = 2 rem 1, 3/1 = 3)
    Rational r(BigInt(3), BigInt(7));
    check_roundtrip(r);
}

TEST(RationalCFTest, NegativeRational) {
    // -7/3: floor(-7/3) = -3, remainder = -7 - (-3)*3 = 2.
    // So CF = [-3; 2] → -3 + 1/2 = -5/2? No: check arithmetic.
    // Actually floor(-7/3) = -3 (since -7/3 ≈ -2.33 → floor = -3).
    // r = -3 + 1/(3/2) = -3 + 2/3 = -7/3. ✓
    Rational r(BigInt(-7), BigInt(3));
    check_roundtrip(r);
    auto cf = r.to_continued_fraction();
    EXPECT_TRUE(cf[0].is_negative()) << "a0 of negative rational must be negative";
    // All subsequent partial quotients must be non-negative.
    for (std::size_t i = 1; i < cf.size(); ++i) {
        EXPECT_FALSE(cf[i].is_negative())
            << "Partial quotient a" << i << " must be non-negative";
    }
}

TEST(RationalCFTest, GoldenRatioApproximant) {
    // 144/89 is a Fibonacci ratio.
    // Actual CF (computed via Euclidean algorithm):
    //   144 = 1·89 + 55  → a0=1
    //   89  = 1·55 + 34  → a1=1
    //   55  = 1·34 + 21  → a2=1
    //   34  = 1·21 + 13  → a3=1
    //   21  = 1·13 +  8  → a4=1
    //   13  = 1· 8 +  5  → a5=1
    //    8  = 1· 5 +  3  → a6=1
    //    5  = 1· 3 +  2  → a7=1
    //    3  = 1· 2 +  1  → a8=1
    //    2  = 2· 1 +  0  → a9=2  (last quotient is 2, not 1)
    // CF = [1; 1, 1, 1, 1, 1, 1, 1, 1, 2]
    Rational r(BigInt(144), BigInt(89));
    check_roundtrip(r);
    auto cf = r.to_continued_fraction();
    ASSERT_EQ(cf.size(), 10U);
    // All quotients except the last must be 1.
    for (std::size_t i = 0; i + 1 < cf.size(); ++i) {
        EXPECT_EQ(cf[i], BigInt(1)) << "cf[" << i << "] should be 1";
    }
    EXPECT_EQ(cf.back(), BigInt(2)) << "last quotient of 144/89 CF is 2";
}

TEST(RationalCFTest, TruncationRoundtripProperty) {
    // For a rational with finite CF, truncating to fewer terms gives
    // a convergent, not the original. But full expansion must round-trip.
    // Test 22/7 (famous pi approximation: CF = [3; 7]).
    Rational r(BigInt(22), BigInt(7));
    auto cf = r.to_continued_fraction();
    ASSERT_EQ(cf.size(), 2U);
    EXPECT_EQ(cf[0], BigInt(3));
    EXPECT_EQ(cf[1], BigInt(7));
    check_roundtrip(r);
}

TEST(RationalCFTest, LargeNumeratorDenominator) {
    // 355/113 → known good CF = [3; 7, 16].
    Rational r(BigInt(355), BigInt(113));
    auto cf = r.to_continued_fraction();
    ASSERT_EQ(cf.size(), 3U);
    EXPECT_EQ(cf[0], BigInt(3));
    EXPECT_EQ(cf[1], BigInt(7));
    EXPECT_EQ(cf[2], BigInt(16));
    check_roundtrip(r);
}

TEST(RationalCFTest, NMaxTruncation) {
    // 355/113 has CF = [3; 7, 16]. Asking for n_max=2 gives [3; 7].
    Rational r(BigInt(355), BigInt(113));
    auto cf2 = r.to_continued_fraction(2U);
    ASSERT_EQ(cf2.size(), 2U);
    EXPECT_EQ(cf2[0], BigInt(3));
    EXPECT_EQ(cf2[1], BigInt(7));

    // Convergent 3 + 1/7 = 22/7 ≠ 355/113 (intentional truncation).
    Rational conv = reconstruct_cf(cf2);
    EXPECT_EQ(conv, Rational(BigInt(22), BigInt(7)));
}

// 20 random-ish rationals round-trip test.
TEST(RationalCFTest, TwentyRandomRoundTrips) {
    struct Case { long long num; long long den; };
    const Case cases[] = {
        {1, 1}, {1, 3}, {2, 5}, {5, 8}, {8, 13}, {13, 21},
        {100, 37}, {-5, 7}, {-11, 4}, {17, 1}, {0, 1},
        {1000, 999}, {999, 1000}, {-100, 101}, {101, -100},
        {7919, 1000}, {1, 7919}, {6, 4}, {12, 8}, {355, 113},
    };
    for (const auto& c : cases) {
        if (c.num == 0) {
            // CF of 0 is [0].
            Rational r{};
            auto cf = r.to_continued_fraction();
            ASSERT_GE(cf.size(), 1U);
            EXPECT_EQ(cf[0], BigInt(0));
            continue;
        }
        auto r_res = Rational::make(BigInt(c.num), BigInt(c.den));
        ASSERT_TRUE(r_res.is_ok());
        check_roundtrip(r_res.value());
    }
}

}  // namespace
