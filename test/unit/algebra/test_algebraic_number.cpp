#include "cas/algebraic_number.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace cas::test {

using namespace cas::algebra;

TEST(AlgebraicNumberTest, BasicArithmetic) {
    // Q(sqrt(2)) -> min_poly = x^2 - 2 -> {-2, 0, 1}
    std::vector<Rational> min_poly = {Rational(BigInt(-2)), Rational(BigInt(0)), Rational(BigInt(1))};
    
    // a = 1 + sqrt(2) -> {1, 1}
    AlgebraicNumber a({Rational(BigInt(1)), Rational(BigInt(1))}, min_poly);
    // b = 2 - sqrt(2) -> {2, -1}
    AlgebraicNumber b({Rational(BigInt(2)), Rational(BigInt(-1))}, min_poly);
    
    // a + b = 3 -> {3, 0}
    auto sum = a + b;
    std::vector<Rational> expected_sum = {Rational(BigInt(3))};
    EXPECT_EQ(sum.value(), expected_sum);
    
    // a * b = (1 + sqrt(2))(2 - sqrt(2)) = 2 - sqrt(2) + 2*sqrt(2) - 2 = sqrt(2) -> {0, 1}
    auto prod = a * b;
    std::vector<Rational> expected_prod = {Rational(BigInt(0)), Rational(BigInt(1))};
    EXPECT_EQ(prod.value(), expected_prod);
}

TEST(AlgebraicNumberTest, Inversion) {
    // Q(sqrt(2))
    std::vector<Rational> min_poly = {Rational(BigInt(-2)), Rational(BigInt(0)), Rational(BigInt(1))};
    
    // a = 1 + sqrt(2)
    AlgebraicNumber a({Rational(BigInt(1)), Rational(BigInt(1))}, min_poly);
    
    // 1 / (1 + sqrt(2)) = (sqrt(2) - 1) / (2 - 1) = sqrt(2) - 1 -> {-1, 1}
    auto inv_res = a.inverse();
    ASSERT_TRUE(inv_res.is_ok());
    
    std::vector<Rational> expected_inv = {Rational(BigInt(-1)), Rational(BigInt(1))};
    EXPECT_EQ(inv_res.value().value(), expected_inv);
    
    // Verify a * inv = 1
    auto identity = a * inv_res.value();
    std::vector<Rational> expected_identity = {Rational(BigInt(1))};
    EXPECT_EQ(identity.value(), expected_identity);
}

} // namespace cas::test
