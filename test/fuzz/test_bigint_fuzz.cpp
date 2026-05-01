#include "cas/bigint.hpp"
#include "../helpers/property_test.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <string>

namespace cas {
namespace {

void expect_bigint_matches_int64(const BigInt& actual, std::int64_t expected) {
    const bool expected_negative = expected < 0;
    const std::uint64_t expected_magnitude = expected_negative
        ? static_cast<std::uint64_t>(-(expected + 1)) + 1U
        : static_cast<std::uint64_t>(expected);

    EXPECT_EQ(actual.decimal(), std::to_string(expected));
    EXPECT_EQ(actual.is_negative(), expected_negative && expected_magnitude != 0U);
}

TEST(BigIntFuzzTest, DeterministicArithmeticMatchesInt64OnSafeRange) {
    test::run_seeded_cases(0xB16B1AULL, 512U, [](test::DeterministicRng& rng, std::size_t) {
        const std::int64_t lhs_value = static_cast<std::int64_t>(rng.next_int(-1000000, 1000000));
        const std::int64_t rhs_value = static_cast<std::int64_t>(rng.next_int(-1000000, 1000000));

        const BigInt lhs(lhs_value);
        const BigInt rhs(rhs_value);

        expect_bigint_matches_int64(lhs + rhs, lhs_value + rhs_value);
        expect_bigint_matches_int64(lhs - rhs, lhs_value - rhs_value);
        expect_bigint_matches_int64(lhs * rhs, lhs_value * rhs_value);

        if (rhs_value != 0) {
            const BigInt quotient = lhs / rhs;
            const BigInt remainder = lhs % rhs;

            EXPECT_EQ(lhs, rhs * quotient + remainder);

            const BigInt rhs_abs = rhs.abs();
            if (!remainder.is_zero()) {
                EXPECT_LT(remainder.abs(), rhs_abs);
            }
        }
    });
}

TEST(BigIntFuzzTest, DeterministicGcdMatchesInt64Reference) {
    test::run_seeded_cases(0x6CD123ULL, 256U, [](test::DeterministicRng& rng, std::size_t) {
        const std::int64_t lhs_value = static_cast<std::int64_t>(rng.next_int(-500000, 500000));
        const std::int64_t rhs_value = static_cast<std::int64_t>(rng.next_int(-500000, 500000));

        const BigInt lhs(lhs_value);
        const BigInt rhs(rhs_value);

        std::int64_t a = lhs_value < 0 ? -lhs_value : lhs_value;
        std::int64_t b = rhs_value < 0 ? -rhs_value : rhs_value;
        while (b != 0) {
            const std::int64_t next = a % b;
            a = b;
            b = next;
        }

        EXPECT_EQ(gcd(lhs, rhs).decimal(), std::to_string(a));
    });
}

}  // namespace
}  // namespace cas
