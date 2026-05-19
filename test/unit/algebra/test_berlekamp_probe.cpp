// F3.2 probe — verify if Berlekamp/Cantor-Zassenhaus factorization
// is a real implementation or a hardcoded FACADE.
//
// Audit C7 claimed: "factorization_integers.cpp — only trial divisors
// + hardcoded Kronecker quadratic, Berlekamp Q-matrix zero traccia,
// distinct-degree NON esiste, equal-degree via random polinomi mai
// testato".
//
// Code inspection of polynomial_modular.cpp:111 reveals real
// distinct_degree_factorization() and equal_degree_factorization()
// implementations. Probes verify they correctly factor mod p.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "../../../src/algebra/polynomial_internal.hpp"

using namespace cas;
using namespace cas::algebra;

namespace {

[[nodiscard]] IntPoly poly(std::initializer_list<long long> coeffs) {
    std::vector<BigInt> bs;
    bs.reserve(coeffs.size());
    for (auto c : coeffs) bs.emplace_back(c);
    return IntPoly(std::move(bs));
}

[[nodiscard]] BigInt poly_eval_mod(const IntPoly& f, const BigInt& x, const BigInt& p) {
    BigInt result(0);
    BigInt power(1);
    for (std::size_t i = 0; i < f.size(); ++i) {
        result = (result + f[i] * power) % p;
        power = (power * x) % p;
    }
    if (result.is_negative()) result = result + p;
    return result;
}

}  // namespace

// f = x² - 1 = (x-1)(x+1) mod 5. Berlekamp should find 2 factors.
TEST(BerlekampProbeTest, XSquaredMinusOneMod5) {
    IntPoly f = poly({-1, 0, 1});
    auto res = factor_polynomial_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_GE(res.value().size(), 2U)
        << "Expected ≥2 factors of x²-1 mod 5";
}

// f = x³ - x = x·(x-1)·(x+1) mod 5. Three distinct linear factors.
TEST(BerlekampProbeTest, XCubedMinusXMod5) {
    IntPoly f = poly({0, -1, 0, 1});
    auto res = factor_polynomial_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_GE(res.value().size(), 3U)
        << "Expected 3 distinct linear factors of x³-x mod 5";
    // Each factor must evaluate to 0 at the expected root.
    // Roots: 0, 1, -1 ≡ 4 mod 5.
    std::vector<int> roots = {0, 1, 4};
    for (int r : roots) {
        bool found = false;
        for (const auto& g : res.value()) {
            if (poly_eval_mod(g, BigInt(r), BigInt(5)) == BigInt(0)) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "Root " << r << " not found in factorisation";
    }
}

// f = x⁴ + 1 mod 5 has 4 roots: 2, 3, ... actually let me check.
// x^4 + 1 ≡ 0 mod 5: x^4 ≡ -1 ≡ 4 mod 5.
// 1^4=1, 2^4=16≡1, 3^4=81≡1, 4^4=256≡1. None give 4.
// So x^4+1 is irreducible mod 5? Or factors into quadratics.
// (2,3 roots of x²+1? 2²=4≡-1 mod 5 yes! 3²=9≡4≡-1 mod 5 yes!)
// So x²+1 = (x-2)(x-3) mod 5 = (x+3)(x+2). Hence x^4+1 = (x²+1)·(x²-1)? No that's x^4-1.
// Skip this test — algebra non trivial.

// Verifies the distinct-degree path: f = x^p - x mod p has all linear
// factors (Fermat's little theorem). For p=7: x^7 - x = x·(x-1)·...·(x-6).
TEST(BerlekampProbeTest, FermatLittleX7MinusXMod7) {
    // x^7 - x  mod 7
    IntPoly f = poly({0, -1, 0, 0, 0, 0, 0, 1});
    auto res = factor_polynomial_mod_p(f, BigInt(7));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Should give 7 linear factors (or product structure equivalent).
    EXPECT_GE(res.value().size(), 7U)
        << "x^7 - x ≡ ∏(x-i) for i=0..6 mod 7: expected 7 linear factors";
}

// Verifies equal-degree path: f = (x² - 2) · (x² - 3) mod 5 has two
// irreducible quadratic factors mod 5 (neither 2 nor 3 is square
// mod 5: 2 is QR? 1²=1, 2²=4, 3²=4, 4²=1 → squares {1,4}; 2,3 not.
// So both x²-2 and x²-3 irreducible mod 5).
TEST(BerlekampProbeTest, TwoQuadraticFactorsMod5) {
    // (x²-2)·(x²-3) = x^4 - 5x² + 6 ≡ x^4 + 0x² + 1 mod 5
    IntPoly f = poly({1, 0, 0, 0, 1});  // x^4 + 1 mod 5 ?
    // Actually x^4 - 5x² + 6 mod 5 = x^4 + 0·x² + 1 = x^4 + 1.
    // x^4 + 1 mod 5 should factor as (x²-2)(x²-3) ≡ (x²+3)(x²+2) mod 5.
    auto res = factor_polynomial_mod_p(f, BigInt(5));
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    EXPECT_GE(res.value().size(), 2U)
        << "x^4+1 mod 5 should split into 2 quadratic factors";
}
