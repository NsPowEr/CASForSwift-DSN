#include "cas/error.hpp"
#include "cas/numtheory.hpp"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace cas::numtheory {
namespace {

[[nodiscard]] Integer parse_integer_or_fail(const char* decimal) {
    auto parsed = Integer::parse(decimal);
    EXPECT_TRUE(parsed.is_ok()) << parsed.error().message;
    return parsed.is_ok() ? std::move(parsed.value()) : Integer(0);
}

void expect_factorization_equals(
    const IntegerFactorization& factorization,
    const Integer& expected_sign,
    const std::vector<std::pair<std::string, unsigned int>>& expected_factors) {
    ASSERT_EQ(factorization.sign, expected_sign);
    ASSERT_EQ(factorization.prime_factors.size(), expected_factors.size());
    for (std::size_t index = 0; index < expected_factors.size(); ++index) {
        EXPECT_EQ(factorization.prime_factors[index].first.decimal(), expected_factors[index].first);
        EXPECT_EQ(factorization.prime_factors[index].second, expected_factors[index].second);
    }
}

TEST(NumberTheoryModularTest, ComputesPowerModExactly) {
    auto result = power_mod(Integer(4), Integer(13), Integer(497));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "445");
}

TEST(NumberTheoryModularTest, ComputesExtendedGcdCoefficients) {
    auto result = extended_gcd(Integer(240), Integer(46));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& [g, x, y] = result.value();
    EXPECT_EQ(g.decimal(), "2");
    EXPECT_EQ((Integer(240) * x + Integer(46) * y).decimal(), "2");
}

TEST(NumberTheoryModularTest, ComputesModularInverse) {
    auto result = modular_inverse(Integer(3), Integer(11));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "4");
}

TEST(NumberTheoryModularTest, RejectsMissingModularInverse) {
    auto result = modular_inverse(Integer(6), Integer(9));
    ASSERT_TRUE(result.is_error());
    EXPECT_EQ(result.error().kind, CASErrorKind::Undefined);
}

TEST(NumberTheoryModularTest, SolvesChineseRemainderTheorem) {
    auto result = chinese_remainder_theorem(
        std::vector<Integer>{Integer(2), Integer(3), Integer(2)},
        std::vector<Integer>{Integer(3), Integer(5), Integer(7)});
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "23");
}

TEST(NumberTheoryPrimalityTest, DetectsPrimeAndCompositeNumbers) {
    auto prime = is_prime(Integer(997));
    auto composite = is_prime(Integer(1001));

    ASSERT_TRUE(prime.is_ok()) << prime.error().message;
    ASSERT_TRUE(composite.is_ok()) << composite.error().message;
    EXPECT_TRUE(prime.value());
    EXPECT_FALSE(composite.value());
}

TEST(NumberTheoryPrimalityTest, MillerRabinHandlesReferenceCases) {
    auto prime = is_prime_miller_rabin(Integer(104729));
    auto composite = is_prime_miller_rabin(Integer(10403));

    ASSERT_TRUE(prime.is_ok()) << prime.error().message;
    ASSERT_TRUE(composite.is_ok()) << composite.error().message;
    EXPECT_TRUE(prime.value());
    EXPECT_FALSE(composite.value());
}

TEST(NumberTheoryPrimalityTest, FindsNextAndNthPrime) {
    auto next = next_prime(Integer(14));
    auto next_after_composite_band = next_prime(Integer(90));
    auto nth = nth_prime(Integer(10));
    auto hundredth = nth_prime(Integer(100));

    ASSERT_TRUE(next.is_ok()) << next.error().message;
    ASSERT_TRUE(next_after_composite_band.is_ok()) << next_after_composite_band.error().message;
    ASSERT_TRUE(nth.is_ok()) << nth.error().message;
    ASSERT_TRUE(hundredth.is_ok()) << hundredth.error().message;
    EXPECT_EQ(next.value().decimal(), "17");
    EXPECT_EQ(next_after_composite_band.value().decimal(), "97");
    EXPECT_EQ(nth.value().decimal(), "29");
    EXPECT_EQ(hundredth.value().decimal(), "541");
}

TEST(NumberTheoryPrimalityTest, RejectsNonPositiveNthPrimeIndex) {
    auto zero = nth_prime(Integer(0));
    ASSERT_TRUE(zero.is_error());
    EXPECT_EQ(zero.error().kind, CASErrorKind::InvalidArgument);
}

TEST(NumberTheoryArithmeticTest, ComputesBinomialWithoutOverflow) {
    auto result = binomial(Integer(30), Integer(15));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "155117520");
}

TEST(NumberTheoryArithmeticTest, SolvesLinearDiophantineEquation) {
    auto result = solve_linear_diophantine(Integer(15), Integer(21), Integer(12));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    const auto& [x, y] = result.value();
    EXPECT_EQ((Integer(15) * x + Integer(21) * y).decimal(), "12");
}

TEST(NumberTheoryFactorizationTest, FactorizesNegativeIntegers) {
    auto result = factor_integer(Integer(-12));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_factorization_equals(
        result.value(),
        Integer(-1),
        {{"2", 2U}, {"3", 1U}});
}

TEST(NumberTheoryFactorizationTest, FactorizesPollardRhoReferenceCase) {
    auto result = factor_integer(Integer(10403));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    expect_factorization_equals(
        result.value(),
        Integer(1),
        {{"101", 1U}, {"103", 1U}});
}

TEST(NumberTheoryFactorizationTest, HandlesOneAndZeroEdgeCases) {
    auto one = factor_integer(Integer(1));
    ASSERT_TRUE(one.is_ok()) << one.error().message;
    EXPECT_EQ(one.value().sign, Integer(1));
    EXPECT_TRUE(one.value().prime_factors.empty());

    auto zero = factor_integer(Integer(0));
    ASSERT_TRUE(zero.is_error());
    EXPECT_EQ(zero.error().kind, CASErrorKind::InvalidArgument);
}

TEST(NumberTheoryArithmeticTest, ComputesEulerPhiFromFactorization) {
    auto result = euler_phi(Integer(36));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_EQ(result.value().decimal(), "12");
}

TEST(NumberTheoryArithmeticTest, ComputesMoebiusMuFromFactorization) {
    auto square_free = moebius_mu(Integer(10));
    auto repeated_factor = moebius_mu(Integer(12));

    ASSERT_TRUE(square_free.is_ok()) << square_free.error().message;
    ASSERT_TRUE(repeated_factor.is_ok()) << repeated_factor.error().message;
    EXPECT_EQ(square_free.value(), 1);
    EXPECT_EQ(repeated_factor.value(), 0);
}

TEST(NumberTheoryArithmeticTest, EnumeratesPositiveDivisorsInOrder) {
    auto result = divisors(Integer(12));
    ASSERT_TRUE(result.is_ok()) << result.error().message;

    ASSERT_EQ(result.value().size(), 6U);
    EXPECT_EQ(result.value()[0].decimal(), "1");
    EXPECT_EQ(result.value()[1].decimal(), "2");
    EXPECT_EQ(result.value()[2].decimal(), "3");
    EXPECT_EQ(result.value()[3].decimal(), "4");
    EXPECT_EQ(result.value()[4].decimal(), "6");
    EXPECT_EQ(result.value()[5].decimal(), "12");
}

TEST(NumberTheoryPrimalityTest, RejectsCompositeNumbersBeyondSixtyFourBits) {
    auto result = is_prime(parse_integer_or_fail("18446744073709551617"));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_FALSE(result.value());
}

// ── P0 Boundary Tests — Miller-Rabin 12-base determinism for n < 2^64 ────────
//
// These tests guard the gap [3.8×10^18, 2^64] where the previous 9-base set
// had no determinism guarantee.  The canonical witness is:
//   3825123056546413051 — strong pseudoprime to all of {2,3,5,7,11,13,17,19,23}
//   but composite (= 149491 × 747451 × 34233211, confirmed by factorization).
//   Base 37 is the first of the extended set that detects it as composite.
// Reference: Sorenson-Webster, Mathematics of Computation 84 (2015), Table 2.
//
// Also tested: 2^64 − 59 = 18446744073709551557, a known prime near the u64 bound
// (confirmed by independent primality provers).

TEST(NumberTheoryPrimalityBoundaryTest, RejectsPseudoprimeToFirst9Bases) {
    // 3825123056546413051 is a strong pseudoprime to all bases {2,3,5,7,11,13,17,19,23}
    // but is composite (= 149491 * 747451 * 34233211).
    // With the corrected 12-base set it must be rejected.
    auto result = is_prime(parse_integer_or_fail("3825123056546413051"));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_FALSE(result.value())
        << "3825123056546413051 is composite but passed 9-base Miller-Rabin; "
           "the 12-base fix (adding 29,31,37) must detect it via base 37";
}

TEST(NumberTheoryPrimalityBoundaryTest, AcceptsLargePrimeNearU64Max) {
    // 2^64 - 59 = 18446744073709551557 is prime.
    // This exercises the final tier of the deterministic branch for n near 2^64.
    auto result = is_prime(parse_integer_or_fail("18446744073709551557"));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_TRUE(result.value())
        << "18446744073709551557 (= 2^64 - 59) is a known prime";
}

TEST(NumberTheoryPrimalityBoundaryTest, RejectsCompositeNearU64MaxGap) {
    // 18446744073709551615 = 2^64 - 1 = 3 * 5 * 17 * 257 * 65537 * 6700417
    // Must be rejected by is_prime.
    auto result = is_prime(parse_integer_or_fail("18446744073709551615"));
    ASSERT_TRUE(result.is_ok()) << result.error().message;
    EXPECT_FALSE(result.value())
        << "18446744073709551615 (= 2^64 - 1) is composite";
}

}  // namespace
}  // namespace cas::numtheory
