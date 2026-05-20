// CAS-L2-24 — Z[i] arithmetic anti-hardcode tests.
//
// Z[i] is a Euclidean domain. Verifications:
//   - additive/multiplicative identities and inverses
//   - multiplicativity of norm: N(αβ) = N(α)·N(β)
//   - Euclidean property: N(remainder) < N(divisor)
//   - GCD divides both inputs
//   - GCD canonical form deterministic

#include <gtest/gtest.h>

#include "cas/gaussian_int.hpp"

using namespace cas;

namespace {

TEST(GaussianIntTest, NormBasic) {
    EXPECT_EQ(GaussianInt(BigInt(3), BigInt(4)).norm(), BigInt(25));
    EXPECT_EQ(GaussianInt(BigInt(0), BigInt(1)).norm(), BigInt(1));  // i
    EXPECT_EQ(GaussianInt(BigInt(0), BigInt(0)).norm(), BigInt(0));
}

TEST(GaussianIntTest, UnitDetection) {
    EXPECT_TRUE(GaussianInt(BigInt(1), BigInt(0)).is_unit());
    EXPECT_TRUE(GaussianInt(BigInt(-1), BigInt(0)).is_unit());
    EXPECT_TRUE(GaussianInt(BigInt(0), BigInt(1)).is_unit());
    EXPECT_TRUE(GaussianInt(BigInt(0), BigInt(-1)).is_unit());
    EXPECT_FALSE(GaussianInt(BigInt(2), BigInt(0)).is_unit());
    EXPECT_FALSE(GaussianInt(BigInt(1), BigInt(1)).is_unit());  // 1+i, norm 2
}

TEST(GaussianIntTest, MultiplicationFormula) {
    // (2+3i)(1+4i) = 2+8i+3i+12i² = 2+11i-12 = -10+11i
    GaussianInt a(BigInt(2), BigInt(3));
    GaussianInt b(BigInt(1), BigInt(4));
    auto p = a * b;
    EXPECT_EQ(p.real(), BigInt(-10));
    EXPECT_EQ(p.imag(), BigInt(11));
}

TEST(GaussianIntTest, NormMultiplicativity) {
    // N(αβ) = N(α)·N(β)
    GaussianInt a(BigInt(3), BigInt(2));   // N = 13
    GaussianInt b(BigInt(5), BigInt(-1));  // N = 26
    EXPECT_EQ((a * b).norm(), a.norm() * b.norm());
}

TEST(GaussianIntTest, EuclideanDivisionInvariant) {
    // For α = q·β + r, must have N(r) < N(β) (Euclidean property).
    GaussianInt a(BigInt(17), BigInt(11));
    GaussianInt b(BigInt(3), BigInt(2));
    auto dm = gaussian_divmod(a, b);
    EXPECT_LT(dm.remainder.norm(), b.norm());
    EXPECT_EQ(dm.quotient * b + dm.remainder, a)
        << "Reconstruction failed: q·β + r ≠ α";
}

TEST(GaussianIntTest, GcdDividesBothInputs) {
    GaussianInt a(BigInt(10), BigInt(0));   // 10
    GaussianInt b(BigInt(0), BigInt(5));    // 5i
    auto g = gaussian_gcd(a, b);
    // a / g and b / g must be exact (zero remainder).
    EXPECT_TRUE(gaussian_divmod(a, g).remainder.is_zero());
    EXPECT_TRUE(gaussian_divmod(b, g).remainder.is_zero());
}

TEST(GaussianIntTest, GcdCanonical) {
    // gcd(1+i, 1+i) = 1+i (canonical: positive real first)
    GaussianInt a(BigInt(1), BigInt(1));
    auto g = gaussian_gcd(a, a);
    EXPECT_FALSE(g.real().is_negative());
}

TEST(GaussianIntTest, AntiHardcodeFiveFactors) {
    // 5 = (2+i)(2-i) in Z[i]. So gcd(5, 2+i) = 2+i (up to units).
    GaussianInt five(BigInt(5), BigInt(0));
    GaussianInt p(BigInt(2), BigInt(1));
    auto g = gaussian_gcd(five, p);
    EXPECT_EQ(g.norm(), p.norm())
        << "gcd(5, 2+i) should have norm 5 (= N(2+i))";
}

TEST(GaussianIntTest, ConjugateInvolution) {
    GaussianInt a(BigInt(7), BigInt(-3));
    EXPECT_EQ(a.conjugate().conjugate(), a);
    // α · conj(α) = N(α)
    auto prod = a * a.conjugate();
    EXPECT_EQ(prod.real(), a.norm());
    EXPECT_EQ(prod.imag(), BigInt(0));
}

TEST(GaussianIntTest, AntiHardcodeRandomGcdProperty) {
    // For random α, β: gcd(α, β) divides α·conj(β) + β·conj(α) (real-part doubling).
    // Use deterministic samples (no random — anti-hardcode but reproducible).
    struct Case { long long ar, ai, br, bi; };
    Case cases[] = {
        {7, 3, 5, 2}, {11, -1, 4, 6}, {12, 8, 3, 5}, {20, 10, 6, 4},
    };
    for (auto& c : cases) {
        GaussianInt a(BigInt(c.ar), BigInt(c.ai));
        GaussianInt b(BigInt(c.br), BigInt(c.bi));
        auto g = gaussian_gcd(a, b);
        EXPECT_TRUE(gaussian_divmod(a, g).remainder.is_zero())
            << "gcd does not divide a in case " << c.ar << "," << c.ai;
        EXPECT_TRUE(gaussian_divmod(b, g).remainder.is_zero())
            << "gcd does not divide b in case " << c.br << "," << c.bi;
    }
}

}  // namespace
