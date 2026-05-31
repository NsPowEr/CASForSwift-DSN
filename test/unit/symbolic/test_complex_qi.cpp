// F1.6-NEW — ComplexRational Q[i] arithmetic tests.
//
// Verifies exact arithmetic over Q[i]:
//   (a+bi)(a-bi) ≡ a²+b²  (norm identity)
//   (1+i)² ≡ 2i
//   1/i ≡ -i
//   abs²(3+4i) ≡ 25 (exact, as Rational)
//   associativity / distributivity spot-check

#include <gtest/gtest.h>

#include "cas/complex_rational.hpp"
#include "cas/bigint.hpp"
#include "cas/rational.hpp"

using namespace cas;

namespace {

// -------------------------------------------------------------------------
// Convenience constructors.
// -------------------------------------------------------------------------

[[nodiscard]] Rational Q(long long n) {
    return Rational(BigInt(n));
}

[[nodiscard]] Rational Q(long long n, long long d) {
    return Rational(BigInt(n), BigInt(d));
}

[[nodiscard]] ComplexRational C(long long re, long long im) {
    return ComplexRational(Q(re), Q(im));
}

// -------------------------------------------------------------------------
// Tests
// -------------------------------------------------------------------------

// F1.6-NEW: (a+bi)(a-bi) ≡ a²+b² for 10 inputs.
TEST(ComplexQiTest, NormIdentityTenInputs) {
    struct Case { long long a, b; };
    const Case cases[] = {
        {1,1}, {3,4}, {5,12}, {-2,3}, {7,-1},
        {0,1}, {1,0}, {-1,-1}, {100,100}, {6,8},
    };
    for (const auto& c : cases) {
        ComplexRational z(Q(c.a), Q(c.b));
        ComplexRational z_conj = z.conjugate();
        ComplexRational product = z * z_conj;

        Rational expected(BigInt(c.a) * BigInt(c.a) + BigInt(c.b) * BigInt(c.b));

        EXPECT_EQ(product.real(), expected)
            << "Re[(a+bi)(a-bi)] should equal a²+b² for a=" << c.a << " b=" << c.b;
        EXPECT_EQ(product.imag(), Rational{})
            << "Im[(a+bi)(a-bi)] must be 0 for a=" << c.a << " b=" << c.b;
    }
}

// F1.6-NEW: (1+i)² ≡ 2i.
TEST(ComplexQiTest, OneISquaredIsTwoI) {
    ComplexRational z = C(1, 1);
    ComplexRational sq = z * z;
    EXPECT_EQ(sq.real(), Q(0));
    EXPECT_EQ(sq.imag(), Q(2));
}

// F1.6-NEW: div(1, i) ≡ -i.
TEST(ComplexQiTest, OneOverIIsMinusI) {
    ComplexRational one = ComplexRational::one();
    ComplexRational i   = ComplexRational::imag_unit();
    auto result = one.divide(i);
    ASSERT_TRUE(result.is_ok());
    EXPECT_EQ(result.value().real(), Q(0));
    EXPECT_EQ(result.value().imag(), Q(-1));
}

// F1.6-NEW: abs²(3+4i) ≡ 25 (exact Rational, no floating-point).
TEST(ComplexQiTest, NormSqThreeFourIsTwentyFive) {
    ComplexRational z = C(3, 4);
    Rational ns = z.norm_sq();
    EXPECT_EQ(ns, Q(25));
}

// F1.6-NEW: div(a+bi, c+di) * (c+di) ≡ a+bi (left-inverse).
TEST(ComplexQiTest, DivisionLeftInverse) {
    ComplexRational z  = C(3, 4);
    ComplexRational w  = C(1, 2);
    auto q_res = z.divide(w);
    ASSERT_TRUE(q_res.is_ok());
    ComplexRational product = q_res.value() * w;
    EXPECT_EQ(product.real(), z.real());
    EXPECT_EQ(product.imag(), z.imag());
}

// F1.6-NEW: rational-coefficient division — (1/2 + i/3) / (1 + i).
TEST(ComplexQiTest, RationalCoefficientDivision) {
    ComplexRational z(Q(1, 2), Q(1, 3));
    ComplexRational w = C(1, 1);
    auto q_res = z.divide(w);
    ASSERT_TRUE(q_res.is_ok());
    // Product should recover z.
    ComplexRational back = q_res.value() * w;
    EXPECT_EQ(back.real(), z.real());
    EXPECT_EQ(back.imag(), z.imag());
}

// F1.6-NEW: divide by zero returns error.
TEST(ComplexQiTest, DivideByZeroIsError) {
    ComplexRational z = C(3, 4);
    ComplexRational zero = ComplexRational::zero();
    auto result = z.divide(zero);
    EXPECT_FALSE(result.is_ok());
}

// F1.6-NEW: is_unit detects ±1 and ±i.
TEST(ComplexQiTest, IsUnitCorrect) {
    EXPECT_TRUE(ComplexRational::one().is_unit());
    EXPECT_TRUE(ComplexRational::imag_unit().is_unit());
    EXPECT_TRUE((-ComplexRational::one()).is_unit());
    EXPECT_TRUE((-ComplexRational::imag_unit()).is_unit());
    EXPECT_FALSE(C(2, 0).is_unit());
    EXPECT_FALSE(C(1, 1).is_unit());
}

// F1.6-NEW: distributivity (a+b)*c ≡ a*c + b*c spot-check.
TEST(ComplexQiTest, DistributivitySpotCheck) {
    ComplexRational a = C(3, -1);
    ComplexRational b = C(2, 5);
    ComplexRational c = C(-1, 4);
    ComplexRational lhs = (a + b) * c;
    ComplexRational rhs = a * c + b * c;
    EXPECT_EQ(lhs.real(), rhs.real());
    EXPECT_EQ(lhs.imag(), rhs.imag());
}

}  // namespace
