// F1.6-NEW — Gaussian prime factorization tests.
//
// Verifies:
//   factor_gaussian(5)  = {2+i, 2-i} (split, p ≡ 1 mod 4)
//   factor_gaussian(13) = {3+2i, 3-2i}
//   factor_gaussian(7)  = {7} (inert, p ≡ 3 mod 4)
//   factor_gaussian(2)  = {1+i}^2 (ramified)
//   product of factors ≡ original (round-trip identity)

#include <gtest/gtest.h>

#include "cas/gaussian_int.hpp"
#include "cas/bigint.hpp"

#include <numeric>

using namespace cas;

namespace {

// -------------------------------------------------------------------------
// Helper: multiply a list of Gaussian prime factors (with exponents)
// by the unit and check they equal n.
// -------------------------------------------------------------------------
void check_product_equals(const GaussianFactorization& f, const BigInt& n) {
    GaussianInt product = f.unit;
    for (const auto& gf : f.factors) {
        GaussianInt pwr(BigInt(1), BigInt(0));
        for (unsigned int e = 0; e < gf.exponent; ++e) {
            pwr = pwr * gf.prime;
        }
        product = product * pwr;
    }
    // product should have zero imaginary part and real part equal to n.
    EXPECT_EQ(product.imag(), BigInt(0))
        << "Product has non-zero imaginary part";
    EXPECT_EQ(product.real().abs(), n)
        << "Product real part does not match input n = " << n.decimal();
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(5) = {2+i, 2-i} (p = 5 ≡ 1 mod 4, split).
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, Five) {
    auto res = factor_gaussian(BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();

    ASSERT_EQ(f.factors.size(), 2U) << "5 should split into two distinct Gaussian primes";

    // Both primes should have norm 5.
    for (const auto& gf : f.factors) {
        EXPECT_EQ(gf.prime.norm(), BigInt(5))
            << "Each Gaussian prime of 5 must have norm 5";
        EXPECT_EQ(gf.exponent, 1U);
    }

    // Product round-trip.
    check_product_equals(f, BigInt(5));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(13) = {3+2i, 3-2i} (p = 13 ≡ 1 mod 4).
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, Thirteen) {
    auto res = factor_gaussian(BigInt(13));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();

    ASSERT_EQ(f.factors.size(), 2U) << "13 should split into two Gaussian primes";
    for (const auto& gf : f.factors) {
        EXPECT_EQ(gf.prime.norm(), BigInt(13));
        EXPECT_EQ(gf.exponent, 1U);
    }
    check_product_equals(f, BigInt(13));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(7) = {7} (inert, p = 7 ≡ 3 mod 4).
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, SevenInert) {
    auto res = factor_gaussian(BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();

    ASSERT_EQ(f.factors.size(), 1U) << "7 is inert: should remain as a single Gaussian prime";
    EXPECT_EQ(f.factors[0].prime.norm(), BigInt(49)) // N(7) = 49
        << "Inert prime 7 in Z[i] has norm 49";
    EXPECT_EQ(f.factors[0].prime.real(), BigInt(7));
    EXPECT_EQ(f.factors[0].prime.imag(), BigInt(0));
    check_product_equals(f, BigInt(7));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(2) — ramified, (1+i)^2 · unit.
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, TwoRamified) {
    auto res = factor_gaussian(BigInt(2));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();

    // Should have exactly one Gaussian prime entry: 1+i with exponent 2.
    ASSERT_EQ(f.factors.size(), 1U) << "2 should have one Gaussian prime factor (1+i)";
    EXPECT_EQ(f.factors[0].prime.real(), BigInt(1));
    EXPECT_EQ(f.factors[0].prime.imag(), BigInt(1));
    EXPECT_EQ(f.factors[0].prime.norm(), BigInt(2));
    EXPECT_EQ(f.factors[0].exponent, 2U);
    check_product_equals(f, BigInt(2));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(1) — trivial (unit, no prime factors).
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, OneIsTrivial) {
    auto res = factor_gaussian(BigInt(1));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();
    EXPECT_TRUE(f.factors.empty()) << "1 has no prime factors";
    check_product_equals(f, BigInt(1));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(25) = {2+i}^2 · {2-i}^2.
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, TwentyFive) {
    auto res = factor_gaussian(BigInt(25));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    const auto& f = res.value();

    ASSERT_EQ(f.factors.size(), 2U) << "25 = 5^2 splits into two Gaussian prime power entries";
    for (const auto& gf : f.factors) {
        EXPECT_EQ(gf.prime.norm(), BigInt(5));
        EXPECT_EQ(gf.exponent, 2U);
    }
    check_product_equals(f, BigInt(25));
}

// -------------------------------------------------------------------------
// Test: factor_gaussian(n ≤ 0) returns error.
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, NonPositiveReturnsError) {
    EXPECT_FALSE(factor_gaussian(BigInt(0)).is_ok());
    EXPECT_FALSE(factor_gaussian(BigInt(-5)).is_ok());
}

// -------------------------------------------------------------------------
// Test: product round-trip for several composites.
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, RoundTripComposites) {
    const BigInt composites[] = {
        BigInt(4), BigInt(9), BigInt(10), BigInt(20),
        BigInt(50), BigInt(100), BigInt(65)  // 65 = 5·13
    };
    for (const auto& n : composites) {
        auto res = factor_gaussian(n);
        ASSERT_TRUE(res.is_ok()) << "factor_gaussian(" << n.decimal() << ") failed: " << res.error().message;
        check_product_equals(res.value(), n);
    }
}

// -------------------------------------------------------------------------
// Test: HPP-018 — Hermite-Serret invariant: norm(alpha) == p for each
// Gaussian prime arising from a split prime p ≡ 1 mod 4.
// Explicitly verifies a²+b² = p for the canonical Gaussian prime.
// -------------------------------------------------------------------------
TEST(GaussianFactorTest, HermiteSerretNormInvariant) {
    // Primes p ≡ 1 mod 4: 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97.
    const BigInt split_primes[] = {
        BigInt(5), BigInt(13), BigInt(17), BigInt(29), BigInt(37),
        BigInt(41), BigInt(53), BigInt(61), BigInt(73), BigInt(89), BigInt(97)
    };
    for (const auto& p : split_primes) {
        auto res = factor_gaussian(p);
        ASSERT_TRUE(res.is_ok())
            << "factor_gaussian(" << p.decimal() << ") failed";
        const auto& f = res.value();
        ASSERT_EQ(f.factors.size(), 2U)
            << "Split prime " << p.decimal() << " must yield two Gaussian primes";
        for (const auto& gf : f.factors) {
            // Core invariant: a² + b² = p (Hermite-Serret norm property).
            BigInt re = gf.prime.real();
            BigInt im = gf.prime.imag();
            BigInt norm_check = re * re + im * im;
            EXPECT_EQ(norm_check, p)
                << "Hermite-Serret invariant FAILED for p=" << p.decimal()
                << ": got prime=" << re.decimal() << "+" << im.decimal() << "i"
                << " with norm=" << norm_check.decimal();
        }
    }
}

}  // namespace
