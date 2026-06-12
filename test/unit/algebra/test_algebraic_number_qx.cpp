#include "cas/algebraic_number_qx.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace cas::test {

using namespace cas::algebra;

TEST(AlgebraicNumberQxTest, BasicArithmetic) {
    // Q(x)[alpha], where alpha^2 = x.
    // r(x) = x -> num = {0, 1}, den = {1}
    std::vector<Rational> r_num = {Rational(BigInt(0)), Rational(BigInt(1))};
    std::vector<Rational> r_den = {Rational(BigInt(1))};

    // a = 1 + alpha -> a = 1, b = 1
    AlgebraicNumberQx a(
        {Rational(1)}, {Rational(1)},
        {Rational(1)}, {Rational(1)},
        r_num, r_den
    );

    // b = 1 - alpha -> a = 1, b = -1
    AlgebraicNumberQx b(
        {Rational(1)}, {Rational(1)},
        {Rational(-1)}, {Rational(1)},
        r_num, r_den
    );

    // a + b = 2
    auto sum = a + b;
    EXPECT_EQ(sum.a_num(), std::vector<Rational>({Rational(2)}));
    EXPECT_EQ(sum.a_den(), std::vector<Rational>({Rational(1)}));
    EXPECT_EQ(sum.b_num(), std::vector<Rational>()); // zero

    // a * b = 1 - alpha^2 = 1 - x
    auto prod = a * b;
    std::vector<Rational> expected_prod_num = {Rational(1), Rational(-1)};
    std::vector<Rational> expected_prod_den = {Rational(1)};
    EXPECT_EQ(prod.a_num(), expected_prod_num);
    EXPECT_EQ(prod.a_den(), expected_prod_den);
    EXPECT_EQ(prod.b_num(), std::vector<Rational>()); // zero
}

TEST(AlgebraicNumberQxTest, Inversion) {
    // Q(x)[alpha], alpha^2 = x
    std::vector<Rational> r_num = {Rational(0), Rational(1)};
    std::vector<Rational> r_den = {Rational(1)};

    // a = 1 + alpha
    AlgebraicNumberQx a(
        {Rational(1)}, {Rational(1)},
        {Rational(1)}, {Rational(1)},
        r_num, r_den
    );

    // 1 / (1 + alpha) = (1 - alpha) / (1 - x)
    // inv_a = 1 / (1 - x) -> num = {1}, den = {1, -1}
    // inv_b = -1 / (1 - x) -> num = {-1}, den = {1, -1}
    auto inv_res = a.inverse();
    ASSERT_TRUE(inv_res.is_ok());

    auto inv = inv_res.value();
    std::vector<Rational> expected_inv_a_num = {Rational(-1)};
    std::vector<Rational> expected_inv_a_den = {Rational(-1), Rational(1)};
    std::vector<Rational> expected_inv_b_num = {Rational(1)};
    std::vector<Rational> expected_inv_b_den = {Rational(-1), Rational(1)};

    EXPECT_EQ(inv.a_num(), expected_inv_a_num);
    EXPECT_EQ(inv.a_den(), expected_inv_a_den);
    EXPECT_EQ(inv.b_num(), expected_inv_b_num);
    EXPECT_EQ(inv.b_den(), expected_inv_b_den);

    // Verify a * inv = 1
    auto identity = a * inv;
    EXPECT_EQ(identity.a_num(), std::vector<Rational>({Rational(1)}));
    EXPECT_EQ(identity.a_den(), std::vector<Rational>({Rational(1)}));
    EXPECT_EQ(identity.b_num(), std::vector<Rational>());
}

} // namespace cas::test
