// F7.0-A3.6 — BigInt limb-allocation budget tests.

#include <gtest/gtest.h>

#include "cas/bigint.hpp"

using namespace cas;

namespace {

TEST(BigIntLimbBudget, DefaultUnlimited) {
    BigInt::clear_budget_exhausted();
    EXPECT_EQ(BigInt::max_limbs(), 0U);
    EXPECT_FALSE(BigInt::budget_exhausted());
}

TEST(BigIntLimbBudget, SetGet) {
    BigInt::set_max_limbs(1024U);
    EXPECT_EQ(BigInt::max_limbs(), 1024U);
    BigInt::set_max_limbs(0U);  // restore
}

TEST(BigIntLimbBudget, SmallMultiplicationUnaffected) {
    BigInt::set_max_limbs(0U);
    BigInt a(123456789);
    BigInt b(987654321);
    auto c = a * b;
    EXPECT_FALSE(c.is_zero());
    EXPECT_FALSE(BigInt::budget_exhausted());
}

TEST(BigIntLimbBudget, LargeMultiplicationBlocked) {
    // Construct two large BigInts via repeated squaring of 2 to get
    // ~512 limbs each. Then multiply with budget 16 limbs → must yield 0
    // and set the exhausted flag.
    BigInt::set_max_limbs(0U);
    BigInt::clear_budget_exhausted();

    BigInt two(2);
    BigInt big = two;
    for (int i = 0; i < 14; ++i) {  // 2^(2^14) ≈ 2^16384 ≈ 512 limbs of 32 bits
        big = big * big;
    }
    ASSERT_GT(big.limb_count(), 16U);

    BigInt::set_max_limbs(16U);
    auto product = big * big;
    EXPECT_TRUE(product.is_zero());
    EXPECT_TRUE(BigInt::budget_exhausted());

    BigInt::set_max_limbs(0U);  // restore for subsequent tests
}

TEST(BigIntLimbBudget, ClearResetsExhaustedFlag) {
    BigInt::set_max_limbs(1U);
    BigInt::clear_budget_exhausted();
    BigInt huge_a = BigInt(1000000) * BigInt(1000000);
    BigInt huge_b = BigInt(1000000) * BigInt(1000000);
    auto p = huge_a * huge_b;
    EXPECT_TRUE(p.is_zero());
    EXPECT_TRUE(BigInt::budget_exhausted());
    BigInt::clear_budget_exhausted();
    EXPECT_FALSE(BigInt::budget_exhausted());
    BigInt::set_max_limbs(0U);
}

}  // namespace
