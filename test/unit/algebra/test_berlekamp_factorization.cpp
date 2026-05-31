// test_berlekamp_factorization.cpp — A3 F2 Block A tests for Berlekamp Fp[x] factorization.
//
// Anti-hardcode requirements (PLAN F2):
//   ≥3 nominal, ≥2 renamed-variable (different primes/polys),
//   ≥2 large coefficients (>10^6 for factors, but modular input),
//   ≥2 high-degree (>50 via budget guard trigger),
//   ≥2 degenerate (degree 0/1), ≥2 out-of-domain Unimplemented,
//   ≥1 metamorphic comparison with Cantor-Zassenhaus.
//
// Correctness certificate: ∏ returned factors ≡ f (mod p).
// Comparison certificate: berlekamp ∩ factor_polynomial_mod_p on same input
// should produce same factor count (both are complete deterministic / probabilistic
// complete algorithms).

#include <gtest/gtest.h>
#include "../../../src/algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

[[nodiscard]] static IntPoly P(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> bs;
    bs.reserve(coeffs.size());
    for (auto c : coeffs) bs.emplace_back(c);
    return IntPoly(std::move(bs));
}

static BigInt fp_mod(const BigInt& a, const BigInt& p) {
    BigInt r = a % p;
    if (r.is_negative()) r += p;
    return r;
}

// Polynomial multiplication mod p (for certificate)
[[nodiscard]] static IntPoly poly_mul_mod_p(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = fp_mod(r[i + j] + a[i] * b[j], p);
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Correctness certificate: ∏ factors ≡ f_monic mod p
[[nodiscard]] static bool factors_product_equals_f(
    const std::vector<IntPoly>& factors,
    IntPoly f,
    const BigInt& p) {
    // Make f monic mod p
    for (auto& c : f.coefficients()) c = fp_mod(c, p);
    f.normalize([](const BigInt& v) { return v.is_zero(); });
    if (!f.is_zero() && !f.leading_coeff().is_zero()) {
        // Modular inverse of leading coeff
        // Simple Fermat: lc^(p-2) mod p
        BigInt inv(1);
        BigInt base = fp_mod(f.leading_coeff(), p);
        BigInt exp = p - BigInt(2);
        while (!exp.is_zero()) {
            if (!(exp % BigInt(2)).is_zero()) inv = fp_mod(inv * base, p);
            base = fp_mod(base * base, p);
            exp = exp / BigInt(2);
        }
        for (auto& c : f.coefficients()) c = fp_mod(c * inv, p);
    }
    f.normalize([](const BigInt& v) { return v.is_zero(); });

    if (factors.empty()) return f.is_zero();
    IntPoly product{std::vector<BigInt>{BigInt(1)}};
    for (const auto& g : factors) {
        product = poly_mul_mod_p(product, g, p);
    }
    if (product.size() != f.size()) return false;
    for (std::size_t i = 0; i < product.size(); ++i) {
        if (product[i] != f[i]) return false;
    }
    return true;
}

} // namespace

// ------- NOMINAL TESTS -------

// N1: x^2 - 1 = (x-1)(x+1) mod 5. Berlekamp should find 2 linear factors.
TEST(BerlekampFactorizationTest, N1_XSquaredMinus1_Mod5) {
    IntPoly f = P({-1, 0, 1});  // x^2 - 1
    auto res = berlekamp_factor_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 2U) << "x^2-1 mod 5 should have 2 linear factors";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(5)))
        << "Certificate: product of factors must equal f mod 5";
}

// N2: x^3 - x = x*(x-1)*(x+1) mod 7. Three distinct linear factors.
TEST(BerlekampFactorizationTest, N2_XCubedMinusX_Mod7) {
    IntPoly f = P({0, -1, 0, 1});  // x^3 - x
    auto res = berlekamp_factor_mod_p(f, BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 3U) << "x^3-x mod 7 should have 3 factors";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(7)))
        << "Certificate: product of factors must equal f mod 7";
}

// N3: x^4 + 1 mod 5.  (x^2-2)*(x^2-3) ≡ (x^2+3)*(x^2+2) mod 5.
// Both factors are irreducible mod 5 (neither 2 nor 3 is a square mod 5).
TEST(BerlekampFactorizationTest, N3_X4Plus1_Mod5_TwoIrreducibleQuadratics) {
    IntPoly f = P({1, 0, 0, 0, 1});  // x^4 + 1
    auto res = berlekamp_factor_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_GE(res.value().size(), 2U) << "x^4+1 mod 5 should split into 2 factors";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(5)))
        << "Certificate: product of factors must equal f mod 5";
}

// ------- RENAMED-VARIABLE EQUIVALENTS (different primes) -------

// RV1: x^2 - 1 mod 11 = (x-1)(x+1). Same structure, different prime.
TEST(BerlekampFactorizationTest, RV1_XSquaredMinus1_Mod11) {
    IntPoly f = P({-1, 0, 1});
    auto res = berlekamp_factor_mod_p(f, BigInt(11));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 2U);
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(11)));
}

// RV2: x^4 - 1 = (x-1)(x+1)(x^2+1) mod 3.
// (x^2+1) mod 3 is irreducible (1 is not QR mod 3? 1^2=1 QR, 2^2=4≡1, so QRs are {0,1}.
// -1≡2 mod 3; is 2 a QR? No. So x^2+1 is irreducible mod 3.)
// So x^4-1 = (x-1)(x+1)(x^2+1) → 3 factors mod 3.
TEST(BerlekampFactorizationTest, RV2_X4Minus1_Mod3_ThreeFactors) {
    IntPoly f = P({-1, 0, 0, 0, 1});  // x^4 - 1
    auto res = berlekamp_factor_mod_p(f, BigInt(3));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 3U) << "x^4-1 mod 3: (x-1)(x+1)(x^2+1)";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(3)));
}

// ------- LARGE COEFFICIENT INPUT (before mod p reduction) -------

// LC1: f = x^2 - 10^6 - 3 mod 5.  10^6+3 mod 5 = 3. So f ≡ x^2 - 3 mod 5.
// QRs mod 5: {1, 4}. 3 is not a QR → irreducible.
TEST(BerlekampFactorizationTest, LC1_LargeInputReducedToIrreducible_Mod5) {
    IntPoly f = P({-(1000000LL + 3LL), 0, 1});  // x^2 - 1000003
    auto res = berlekamp_factor_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // x^2 - 3 mod 5 is irreducible, so we expect 1 factor
    EXPECT_EQ(res.value().size(), 1U) << "x^2-1000003 ≡ x^2-3 mod 5 is irreducible";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(5)));
}

// LC2: f with coefficients > 10^6, reduces to a splittable polynomial mod p.
// x^2 - 2000007 mod 7: 2000007 mod 7 = 2000007 - 285715*7 = 2000007-1999995=12 mod7=5.
// Actually 285715*7 = 1999905, 2000007 - 1999905 = 102; 102/7=14r4; so 2000007 mod 7 = 4.
// Wait: 2000007 = 7*285715 + r. 7*285715 = 1999905. 2000007-1999905 = 102. 102=7*14+4. So r=4.
// x^2 - 4 ≡ (x-2)(x+2) mod 7.
TEST(BerlekampFactorizationTest, LC2_LargeInputReducedToSplittable_Mod7) {
    IntPoly f = P({-2000007LL, 0, 1});  // x^2 - 2000007
    auto res = berlekamp_factor_mod_p(f, BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 2U) << "x^2-2000007 ≡ x^2-4=(x-2)(x+2) mod 7";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(7)));
}

// ------- DEGENERATE CASES -------

// D1: Degree-1 polynomial — should return itself (irreducible linear).
TEST(BerlekampFactorizationTest, D1_Degree1_ReturnsItself) {
    IntPoly f = P({3, 1});  // x + 3 mod any prime
    auto res = berlekamp_factor_mod_p(f, BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_EQ(res.value().size(), 1U) << "Degree-1 should be irreducible (1 factor)";
    EXPECT_TRUE(factors_product_equals_f(res.value(), f, BigInt(7)));
}

// D2: Degree-0 polynomial (constant) — should return empty.
TEST(BerlekampFactorizationTest, D2_Degree0_Constant_ReturnsEmpty) {
    IntPoly f = P({5});  // constant 5
    auto res = berlekamp_factor_mod_p(f, BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Constant polynomial after making monic is 1 (degree 0), return empty or {1}.
    // The implementation returns empty for degree-0 after monic normalization.
    EXPECT_TRUE(res.value().empty() || (res.value().size() == 1U && res.value()[0].degree() == 0U));
}

// ------- OUT-OF-DOMAIN: BUDGET GUARD -------

// OD1: degree * p > 1024 → Unimplemented with reason BERLEKAMP_MATRIX_TOO_LARGE.
// Use deg=30, p=37 → 30*37=1110>1024.
TEST(BerlekampFactorizationTest, OD1_LargeDegreePrimeBudgetGuard) {
    // Construct f = x^30 + 1 (arbitrary, content doesn't matter for budget check)
    std::vector<BigInt> coeffs(31, BigInt(0));
    coeffs[0] = BigInt(1);
    coeffs[30] = BigInt(1);
    IntPoly f(std::move(coeffs));

    auto res = berlekamp_factor_mod_p(f, BigInt(37));
    ASSERT_TRUE(res.is_error()) << "Expected Unimplemented for deg*p > 1024";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented)
        << "Should return Unimplemented for budget exceeded";
}

// OD2: deg=2, p=1009 (prime, large) → 2*1009=2018 > 1024.
TEST(BerlekampFactorizationTest, OD2_LargePrimeBudgetGuard) {
    IntPoly f = P({-1, 0, 1});  // x^2 - 1
    // p=1009 is prime; deg*p = 2*1009 = 2018 > 1024
    auto res = berlekamp_factor_mod_p(f, BigInt(1009));
    ASSERT_TRUE(res.is_error()) << "Expected Unimplemented for p=1009, deg=2";
    EXPECT_EQ(res.error().kind, CASErrorKind::Unimplemented);
}

// ------- METAMORPHIC: COMPARE WITH CANTOR-ZASSENHAUS -------

// MM1: Both Berlekamp and Cantor-Zassenhaus should produce same number of factors.
// Use small deg*p where Berlekamp doesn't trigger guard.
// x^4 - 1 mod 5 = (x-1)(x+1)(x^2+1)? Let's check: QRs mod 5: 1,4. -1≡4 is QR (2^2=4).
// So x^2+1 = x^2-(-1) ≡ (x-2)(x+2) mod 5? 2^2=4≡-1. Yes!
// So x^4-1 = (x-1)(x+1)(x-2)(x+2) mod 5: 4 linear factors.
TEST(BerlekampFactorizationTest, MM1_CompareWithCantorZassenhaus_X4Minus1_Mod5) {
    IntPoly f = P({-1, 0, 0, 0, 1});  // x^4 - 1
    BigInt p(5);

    auto berlekamp_res = berlekamp_factor_mod_p(f, p);
    auto cz_res = factor_polynomial_mod_p(f, p);

    ASSERT_TRUE(berlekamp_res.is_ok()) << berlekamp_res.error().message;
    ASSERT_TRUE(cz_res.is_ok()) << cz_res.error().message;

    // Both should find 4 factors
    EXPECT_EQ(berlekamp_res.value().size(), cz_res.value().size())
        << "Berlekamp and CZ should find the same number of factors for x^4-1 mod 5";

    // Berlekamp certificate
    EXPECT_TRUE(factors_product_equals_f(berlekamp_res.value(), f, p))
        << "Berlekamp: product of factors must equal f mod 5";

    // CZ certificate
    EXPECT_TRUE(factors_product_equals_f(cz_res.value(), f, p))
        << "CZ: product of factors must equal f mod 5";
}

// MM2: Irreducible polynomial — both should return 1 factor.
// x^2 + 1 mod 3 is irreducible (QRs mod 3 = {0,1}, -1≡2 not QR).
TEST(BerlekampFactorizationTest, MM2_Irreducible_BothReturnOneFactorMod3) {
    IntPoly f = P({1, 0, 1});  // x^2 + 1
    BigInt p(3);

    auto berlekamp_res = berlekamp_factor_mod_p(f, p);
    auto cz_res = factor_polynomial_mod_p(f, p);

    ASSERT_TRUE(berlekamp_res.is_ok());
    ASSERT_TRUE(cz_res.is_ok());

    EXPECT_EQ(berlekamp_res.value().size(), 1U) << "Berlekamp: x^2+1 mod 3 is irreducible";
    EXPECT_EQ(cz_res.value().size(), 1U)        << "CZ: x^2+1 mod 3 is irreducible";
    EXPECT_TRUE(factors_product_equals_f(berlekamp_res.value(), f, p));
    EXPECT_TRUE(factors_product_equals_f(cz_res.value(), f, p));
}
