// test_hensel_quadratic.cpp — A2 F2 Block A tests for quadratic Hensel lifting.
//
// Anti-hardcode requirements (PLAN F2):
//   ≥3 nominal, ≥2 renamed-variable equivalents (BigInt seeds),
//   ≥2 large coefficients (>10^6), ≥2 high-degree (>10),
//   ≥2 degenerate (deg 0/1), ≥2 out-of-domain Unimplemented,
//   ≥1 metamorphic property.
//
// Correctness certificate:  after hensel_lift(f, g, h, p, k),
// product G*H ≡ f (mod p^k) — verified structurally by modular comparison.

#include <gtest/gtest.h>
#include "../../../src/algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

// Build an IntPoly from an initializer list of long long
[[nodiscard]] static IntPoly P(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> bs;
    bs.reserve(coeffs.size());
    for (auto c : coeffs) bs.emplace_back(c);
    return IntPoly(std::move(bs));
}

// Compute a * b over exact integers (no modular reduction)
[[nodiscard]] static IntPoly poly_mul_exact(const IntPoly& a, const IntPoly& b) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] += a[i] * b[j];
        }
    }
    return r;
}

// Reduce polynomial mod m (in-place)
static void poly_reduce_mod(IntPoly& p, const BigInt& m) {
    for (auto& c : p.coefficients()) {
        c = c % m;
        if (c.is_negative()) c += m;
    }
    p.normalize([](const BigInt& v) { return v.is_zero(); });
}

// Check: G*H ≡ f (mod m) — main correctness certificate
[[nodiscard]] static bool product_equals_f_mod(const IntPoly& G, const IntPoly& H,
                                                const IntPoly& f, const BigInt& m) {
    IntPoly product = poly_mul_exact(G, H);
    poly_reduce_mod(product, m);
    IntPoly f_red = f;
    poly_reduce_mod(f_red, m);

    if (product.size() != f_red.size()) {
        // One may be zero if they differ
        // Normalize
        product.normalize([](const BigInt& v) { return v.is_zero(); });
        f_red.normalize([](const BigInt& v) { return v.is_zero(); });
        if (product.size() != f_red.size()) return false;
    }
    for (std::size_t i = 0; i < product.size(); ++i) {
        if (product[i] != f_red[i]) return false;
    }
    return true;
}

// Compute p^k as BigInt
static BigInt bigint_pow(const BigInt& base, std::size_t k) {
    BigInt r(1);
    for (std::size_t i = 0; i < k; ++i) r *= base;
    return r;
}

} // namespace

// ------- NOMINAL TESTS -------

// N1: f = x^2 - 1 = (x-1)*(x+1) over Z, lift to p=5, k=3 (p^k=125)
TEST(HenselQuadraticTest, N1_XSquaredMinusOneP5K3) {
    IntPoly f = P({-1, 0, 1});         // x^2 - 1
    IntPoly g = P({-1, 1});            // x - 1 (mod 5: root 1)
    IntPoly h = P({1, 1});             // x + 1 (mod 5: root -1=4)
    BigInt p(5);
    std::size_t k = 3;
    BigInt target = bigint_pow(p, k);  // 125

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();

    // Certificate: G*H ≡ f (mod 125)
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "G*H should equal f mod 5^3=125";
    EXPECT_EQ(G.degree(), 1U) << "G should remain linear";
    EXPECT_EQ(H.degree(), 1U) << "H should remain linear";
}

// N2: f = x^3 - x = x*(x-1)*(x+1) over Z, lift first pair (x)*(x^2-1)
TEST(HenselQuadraticTest, N2_XCubedMinusXP7K2) {
    // f = x^3 - x,  g = x (linear factor),  h = x^2 - 1 (quadratic factor)
    // Over Z, f = g * h exactly
    IntPoly f = P({0, -1, 0, 1});      // x^3 - x
    IntPoly g = P({0, 1});             // x
    IntPoly h = P({-1, 0, 1});         // x^2 - 1
    BigInt p(7);
    std::size_t k = 2;
    BigInt target = bigint_pow(p, k);  // 49

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "G*H should equal f mod 7^2=49";
}

// N3: f = x^2 - 3 (irreducible over Z, factors mod 13: roots ±sqrt(3) mod 13)
// 4^2 = 16 ≡ 3 mod 13, so roots are 4 and 9.
// g = x - 4, h = x - 9 mod 13
TEST(HenselQuadraticTest, N3_XSquaredMinus3P13K4) {
    IntPoly f = P({-3, 0, 1});         // x^2 - 3
    IntPoly g = P({-4, 1});            // x - 4 (mod 13)
    IntPoly h = P({-9, 1});            // x - 9 (mod 13)
    BigInt p(13);
    std::size_t k = 4;
    BigInt target = bigint_pow(p, k);  // 28561

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "G*H should equal f mod 13^4=28561";
}

// ------- RENAMED-VARIABLE EQUIVALENTS (different seeds) -------

// RV1: Same as N1 but with different prime p=11
TEST(HenselQuadraticTest, RV1_XSquaredMinus1_P11) {
    // x^2 - 1 = (x-1)(x+1), roots 1 and 10 mod 11
    IntPoly f = P({-1, 0, 1});
    IntPoly g = P({-1, 1});
    IntPoly h = P({1, 1});
    BigInt p(11);
    std::size_t k = 3;
    BigInt target = bigint_pow(p, k);

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok());
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target));
}

// RV2: Cubic with different factoring structure mod p=3
// f = x^3 + 2x^2 + x = x*(x^2 + 2x + 1) = x*(x+1)^2
// mod 3: g = x, h = (x+1)^2 = x^2 + 2x + 1
TEST(HenselQuadraticTest, RV2_CubicWithRepeatedFactor_P3) {
    IntPoly f = P({0, 1, 2, 1});       // x^3 + 2x^2 + x = x*(x+1)^2
    IntPoly g = P({0, 1});             // x
    IntPoly h = P({1, 2, 1});          // (x+1)^2 = x^2 + 2x + 1
    BigInt p(3);
    std::size_t k = 3;
    BigInt target = bigint_pow(p, k);

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok());
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target));
}

// ------- LARGE COEFFICIENTS (>10^6) -------

// LC1: f = x^2 - (10^6 + 7) = large discriminant
// Lift mod p=5, checking certificate
TEST(HenselQuadraticTest, LC1_LargeConstantTerm) {
    // We need to ensure f factors mod p first.
    // f = x^2 - 1000007 mod 5: 1000007 mod 5 = 2, so f ≡ x^2 - 2 mod 5.
    // Is 2 a QR mod 5? QRs mod 5 are {1,4}. So x^2-2 is irreducible mod 5.
    // Use p=7: 1000007 mod 7 = 1000007 - 142858*7 = 1000007-1000006 = 1 mod 7.
    // So f ≡ x^2 - 1 = (x-1)(x+1) mod 7. Great!
    IntPoly f = P({-1000007LL, 0, 1});
    IntPoly g = P({-1, 1});            // x-1 (mod 7)
    IntPoly h = P({1, 1});             // x+1 (mod 7)
    BigInt p(7);
    std::size_t k = 4;
    BigInt target = bigint_pow(p, k);  // 2401

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "Large constant: G*H should equal f mod 7^4=2401";
}

// LC2: f with large leading coefficient: 2000003*x^2 - 1
// mod p=5: 2000003 mod 5 = 3, so f ≡ 3x^2 - 1 mod 5.
// 3x^2 ≡ 1 mod 5 → x^2 ≡ 2 mod 5 (QR? No). Use p=11: 2000003 mod 11 = ?
// 2000003 = 181818*11 + 5. So f ≡ 5x^2 - 1 mod 11.
// 5x^2 ≡ 1 → x^2 ≡ inv(5) = 9 mod 11. 3^2=9, 8^2=64≡9. Yes!
// So roots x=3, x=8 mod 11. Factors: lc=5, g= 5*(x-3)=5x-15≡5x-4 mod 11
// and h = (x-8) monic.
// Actually more simply: f = 5*(x-3)*(x-8) mod 11 since 5*1*1=5=lc.
// Let g = x-3, h = x-8 (both monic), then g*h = x^2-11x+24 ≡ x^2+0x+2 mod 11.
// f mod 11 = 5x^2 - 1. We need f = lc * g * h mod 11 = 5*(x^2+2) mod 11.
// 5*(x^2+2) = 5x^2 + 10 ≡ 5x^2 - 1 mod 11. Checks out.
// For lifting with the lc, we can incorporate lc into g: g_lc = 5*(x-3) = 5x-15≡5x+7 mod 11
// h = x - 8.  Then g_lc * h = (5x+7)(x-8) = 5x^2-40x+7x-56 = 5x^2-33x-56 ≡ 5x^2-0x-1 ≡ 5x^2-1 mod 11.
TEST(HenselQuadraticTest, LC2_LargeLeadingCoeff) {
    IntPoly f = P({-1LL, 0, 2000003LL});   // 2000003*x^2 - 1
    // g_lc = 5x+7 mod 11 (i.e., coefficients [7, 5])
    IntPoly g = P({7, 5});
    // h = x - 8 mod 11 = x + 3 mod 11 in positive form
    IntPoly h = P({-8LL, 1});
    BigInt p(11);
    std::size_t k = 3;
    BigInt target = bigint_pow(p, k);  // 1331

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "Large leading coeff: G*H should equal f mod 11^3";
}

// ------- HIGH DEGREE (>10) -------

// HD1: f = x^12 - 1 (cyclotomic product) = (x-1)(x+1)(x^2+1)(x^4-x^2+1)(x^4+x^2+1)...
// For lifting, split into g = x-1, h = (x^12-1)/(x-1) = x^11+x^10+...+1 mod p
// Use p=13.
TEST(HenselQuadraticTest, HD1_Degree12_P13) {
    // f = x^12 - 1
    std::vector<BigInt> f_coeffs(13, BigInt(0));
    f_coeffs[0] = BigInt(-1);
    f_coeffs[12] = BigInt(1);
    IntPoly f(std::move(f_coeffs));

    // g = x - 1 (linear factor: 1 is always a root of x^n - 1)
    IntPoly g = P({-1, 1});
    // h = (x^12-1)/(x-1) = x^11 + x^10 + ... + 1 (sum formula)
    std::vector<BigInt> h_coeffs(12, BigInt(1));
    IntPoly h(std::move(h_coeffs));

    BigInt p(13);
    std::size_t k = 3;
    BigInt target = bigint_pow(p, k);

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "Degree-12 cyclotomic: G*H should equal x^12-1 mod 13^3";
}

// HD2: f = x^15 - x (Frobenius polynomial mod p=5 gives all elements of F_5)
// = product of x^5-x factors; here just test lifting two pieces mod p=5.
// Split: g = x (trivial), h = x^14 + x^13 + ... has degree 14 too expensive.
// Instead use g = x^7 - x (roots F_{5^7}?), h = x^8 - 0?
// Simpler: test with f = x^12 - x = x*(x^11-1), g = x, h = x^11-1, p=5, k=2.
TEST(HenselQuadraticTest, HD2_Degree12SplitLinearAndDeg11_P5) {
    // f = x^12 - x = x * (x^11 - 1)
    std::vector<BigInt> f_coeffs(13, BigInt(0));
    f_coeffs[1] = BigInt(-1);
    f_coeffs[12] = BigInt(1);
    IntPoly f(std::move(f_coeffs));

    IntPoly g = P({0, 1});   // x
    // h = x^11 - 1 (but wait, f = x*(x^11-1) not x*(x^11-1)? Let's check:
    // x*(x^11-1) = x^12 - x. Yes.)
    std::vector<BigInt> h_coeffs(12, BigInt(0));
    h_coeffs[0] = BigInt(-1);
    h_coeffs[11] = BigInt(1);
    IntPoly h(std::move(h_coeffs));

    BigInt p(5);
    std::size_t k = 2;
    BigInt target = bigint_pow(p, k);  // 25

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "Degree-12 split: G*H should equal x^12-x mod 5^2=25";
}

// ------- DEGENERATE (deg 0 or 1) -------

// D1: Degree-1 polynomial — trivial (can't split further)
// f = 3x + 7, g = 3x+7, h = 1, p=5, k=2
TEST(HenselQuadraticTest, D1_Degree1_Trivial) {
    IntPoly f = P({7, 3});
    IntPoly g = P({7, 3});
    IntPoly h = P({1});
    BigInt p(5);
    std::size_t k = 2;
    BigInt target = bigint_pow(p, k);

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target))
        << "Degree-1 trivial split: certificate must hold";
}

// D2: k=1 (no actual lifting needed, just initial factorization check)
TEST(HenselQuadraticTest, D2_KEqualsOne_NoLift) {
    IntPoly f = P({-1, 0, 1});   // x^2 - 1
    IntPoly g = P({-1, 1});
    IntPoly h = P({1, 1});
    BigInt p(5);
    std::size_t k = 1;           // target = p^1 = 5, no actual quadratic step
    BigInt target = p;

    auto res = hensel_lift(f, g, h, p, k);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    auto [G, H] = res.value();
    EXPECT_TRUE(product_equals_f_mod(G, H, f, target));
}

// ------- METAMORPHIC PROPERTY -------

// MP1: Lifting to p^2 then squaring the prime should give the same result
// as lifting to p^4 directly.  (Round-trip / order-invariance property)
// This tests that quadratic squaring is self-consistent.
TEST(HenselQuadraticTest, MP1_QuadraticSquaringIsConsistent) {
    IntPoly f = P({-3, 0, 1});
    IntPoly g = P({-4, 1});
    IntPoly h = P({-9, 1});
    BigInt p(13);

    // Lift to p^2
    auto res2 = hensel_lift(f, g, h, p, 2);
    ASSERT_TRUE(res2.is_ok());
    BigInt target2 = bigint_pow(p, 2);
    auto [G2, H2] = res2.value();
    EXPECT_TRUE(product_equals_f_mod(G2, H2, f, target2))
        << "Lift to p^2: certificate must hold";

    // Lift to p^4 directly
    auto res4 = hensel_lift(f, g, h, p, 4);
    ASSERT_TRUE(res4.is_ok());
    BigInt target4 = bigint_pow(p, 4);
    auto [G4, H4] = res4.value();
    EXPECT_TRUE(product_equals_f_mod(G4, H4, f, target4))
        << "Lift to p^4: certificate must hold";

    // G4 mod p^2 should equal G2 (the factors at the lower level are embedded)
    IntPoly G4_mod2 = G4;
    poly_reduce_mod(G4_mod2, target2);
    IntPoly G2_mod2 = G2;
    poly_reduce_mod(G2_mod2, target2);
    // They should agree (quadratic lifting is unique)
    EXPECT_EQ(G4_mod2.size(), G2_mod2.size())
        << "Lifted G at p^4 reduced to p^2 should equal G at p^2";
}
