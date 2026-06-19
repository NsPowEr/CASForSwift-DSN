#include "cas/numtheory.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace cas::numtheory {

TEST(MillerRabinTest, SmallPrimes) {
    const std::vector<int64_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (int64_t p : primes) {
        auto result = is_prime(Integer(p));
        ASSERT_TRUE(result.is_ok());
        EXPECT_TRUE(result.value()) << p << " should be prime";
    }
}

TEST(MillerRabinTest, SmallComposites) {
    const std::vector<int64_t> composites = {4, 6, 8, 9, 10, 12, 14, 15, 21, 25, 27, 33, 35, 49};
    for (int64_t c : composites) {
        auto result = is_prime(Integer(c));
        ASSERT_TRUE(result.is_ok());
        EXPECT_FALSE(result.value()) << c << " should be composite";
    }
}

TEST(MillerRabinTest, EdgeCases) {
    EXPECT_FALSE(is_prime(Integer(0)).value());
    EXPECT_FALSE(is_prime(Integer(1)).value());
    EXPECT_FALSE(is_prime(Integer(-1)).value());
    EXPECT_FALSE(is_prime(Integer(-7)).value());
}

TEST(MillerRabinTest, DeterministicBases) {
    // 2047 = 23 * 89 (Primo pseudoprimo forte in base 2)
    EXPECT_FALSE(is_prime(Integer(2047)).value());
    
    // 1373653 = 829 * 1657 (Primo pseudoprimo forte in basi 2, 3)
    EXPECT_FALSE(is_prime(Integer(1373653)).value());
    
    // Un numero primo grande (che rientra in u64)
    // 2^31 - 1 (Mersenne prime)
    EXPECT_TRUE(is_prime(Integer(2147483647)).value());
}

TEST(MillerRabinTest, LargePrime) {
    // 10^12 + 39 è primo (oltre 32 bit, sotto 64 bit)
    EXPECT_TRUE(is_prime(Integer(1000000000039LL)).value());

    // 10^12 + 1 = 73 * 137 * 99990001 (composto noto: (10^4)^3+1 = 10001·99990001, 10001=73·137)
    EXPECT_FALSE(is_prime(Integer(1000000000001LL)).value());
}

} // namespace cas::numtheory
